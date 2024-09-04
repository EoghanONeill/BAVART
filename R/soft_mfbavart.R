


#' @title Mixed Frequency Soft Bayesian Additive Vector Autoregression Trees (Huber et al. 2023)
#'
#' @description Mixed Frequency Soft Bayesian Additive Vector Autoregression Trees (Huber et al. 2023).
#' @import dbarts
#' @import MASS
#' @import stochvol
#' @import mfbvar
#' @import abind
#' @import collapse
#' @import SoftBart
#' @param data list of data frames for each variable.
#' @param ite option for constructing the loading matrix. Can be "m", "lvl", or "grw".
#' @param p number of lags
#' @param fhorz number of forecast horizons
#' @return The following objects are returned:
#' \item{Y}{Array of in-sample posterior draws. First dimension is number of MCMC iterations, second dimension is number of training time periods, third dimension is number of variables (length of input data list).}
#' \item{fcst}{Array of out-of-sample posterior draws (including noise form error term). First dimension is number of MCMC iterations, second dimension is number of horecast horizons, third dimension is number of variables.}
#' \item{Yq}{Array of out-of-sample quarterly aggregates. First dimension is number of MCMC iterations, second dimension is number of horecast horizons, third dimension is number of variables.}
#' @examples
#' # -----------------------------------------------------------------------------
#' # load dataset
#' library(alfred) # load dataset from FRED (Fed St. Louis)
#'
#' variables <- c("CPIAUCSL","UNRATE","GDPC1")
#' out <- lapply(variables, get_alfred_series,observation_start = "1980-01-01",observation_end = "2020-12-01",realtime_start = "2020-12-01",realtime_end = "2020-12-01")
#' alfred_to_ts <- function(x, freq){
#'   ts(x[,3],start=c(1980,1),frequency=freq)
#' }
#' mf_list <- mapply(alfred_to_ts, x = out, freq = c(12, 12, 4))
#' names(mf_list) <- variables
#' log_diff <- function(x) {
#'   freq <- frequency(x)
#'   100 * freq * diff(log(x))}
#' mf_list[c("CPIAUCSL", "GDPC1")] <- lapply(mf_list[c("CPIAUCSL", "GDPC1")], log_diff)
#' mf_list <- mapply(window, x = mf_list, start = list(c(1980, 4), c(1980, 4), c(1980, 2)))
#' data <- mf_list
#'
#' est_obj <- mfbavart(data,itr="grw",fhorz=2,prior.sig=c(200,0.75))
#'
#' Yq <- est_obj$Yq
#' GDPC1_post <- t(apply(Yq,c(2,3),quantile,probs=c(0.16,0.5,0.84),na.rm=TRUE)[,,"GDPC1"])
#' ts.plot(GDPC1_post)
#' @export
soft_mfbavart <- function(data,itr,p=5,fhorz=0,cons=FALSE,exact=FALSE,sv=FALSE,var.thrsh=10,max.count.var=10,
                     cgm.level=0.95,cgm.exp=2,sd.mu=2,num.trees=250,#prior.sig,
                     nburn=1000,nsave=1000,thinfac=1,
                     quiet=FALSE,
                     SB_group = NULL,
                     SB_alpha = 1,
                     SB_beta = 2,
                     SB_gamma = 0.95,
                     SB_k = 2,
                     SB_sigma_hat = NULL,
                     SB_shape = 1,
                     SB_width = 0.1,
                     # SB_num_tree = 20,
                     SB_alpha_scale = NULL,
                     SB_alpha_shape_1 = 0.5,
                     SB_alpha_shape_2 = 1,
                     SB_tau_rate = 10,
                     SB_num_tree_prob = NULL,
                     SB_temperature = 1,
                     SB_weights = NULL,
                     SB_normalize_Y = TRUE){


  # construct design matrices
  Ylist <- list_to_matrix(data)
  Yraw <- Ylist$Yraw
  Y_fq <- Ylist$freq

  # standardize data
  Ymu <- apply(Yraw,2,mean,na.rm=TRUE)
  Ysd <- apply(Yraw,2,sd,na.rm=TRUE)
  Yraw <- apply(Yraw,2,function(x){(x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)})

  # MCMC setup
  nthin <- round(thinfac * nsave)
  ntot <- nburn + nthin
  thin.set <- floor(seq(nburn+1,ntot,length.out=nsave))
  in.thin <- 0

  # selection and frequencies -----------------------------------------------------------------------------
  # selection and frequencies
  M_h <- sum(Y_fq=="h")
  M_l <- sum(Y_fq=="l")
  M <- M_h+M_l

  sl_h <- which(Y_fq=="h")
  sl_l <- which(Y_fq=="l")

  Y_obs <- 1-is.na(Yraw)

  # create design matrices Y/X
  Y <- Yact <- Yraw # leave original input matrix unchanged
  Y <- fill_na(Y) # fill missing values
  if(cons){
    X <- cbind(mlag(Y,p),1)[(p+1):nrow(Y),]
  }else{
    X <- cbind(mlag(Y,p))[(p+1):nrow(Y),]
  }
  Yinit <- Y[1:p,]
  Y <- Y[(p+1):nrow(Y),]
  Yact <- Yact[(p+1):nrow(Yact),]
  Y_obs <- Y_obs[(p+1):nrow(Y_obs),]

  # for naming output
  dates <- rownames(Yraw)[(p+1):nrow(Yraw)]
  dates_fc <- as.character(seq(as.Date(dates[length(dates)]),by="month",length.out=fhorz+1)[2:(fhorz+1)])

  # dimensions
  Tnobs <- nrow(Y)
  K <- ncol(X)

  # define last balanced observation
  if(sum(Y_obs[,1:M_h]==0)){
    T_b <- as.numeric(which(apply(Y_obs[,1:M_h]==0,1,sum)>0))[1]-1
  }else{
    T_b <- Tnobs
  }

  # get loadings/lambda for sampling of states with c++ function
  lvl <- c(1,1,1)/3 # aggregation scheme for levels
  grw <- c(1,2,3,2,1)/9 # aggregation scheme for log-growth rates
  Lambda <- build_Lambda2(itr,p)
  Loadings <- matrix(0,M_l,p)
  for(mm in seq_len(M)){
    if(mm>M_h){
      if(itr[mm-M_h]=="grw"){
        Loadings[mm-M_h,] <- grw
      }else if(itr[mm-M_h]=="lvl"){
        Loadings[mm-M_h,] <- lvl
      }
    }
  }

  # initialize draws -----------------------------------------------------------------------------
  # initialize draws
  A_draw <- A_OLS <- solve(crossprod(X))%*%crossprod(X,Y)
  Sig_OLS <- crossprod(Y-X%*%A_OLS)/Tnobs

  Sig_t <- array(0,dim=c(Tnobs,M,M))
  for(tt in 1:Tnobs){
    Sig_t[tt,,] <- Sig_OLS
  }

  # covariance related objects
  eta <- matrix(NA,Tnobs,M)
  H <- matrix(-3,Tnobs,M)
  A0_draw <- diag(M)

  # stochastic volatility using stochvol package
  sv_draw <- list()
  sv_latent <- list()
  for (mm in seq_len(M)){
    sv_draw[[mm]] <- list(mu = 0, phi = 0.99, sigma = 0.01, nu = Inf, rho = 0, beta = NA, latent0 = 0)
    sv_latent[[mm]] <- rep(0,Tnobs)
  }

  # construct priors for SV
  sv_priors <- list()
  if(sv){
    for(mm in 1:M){
      sv_priors[[mm]] <- specify_priors(
        mu = sv_normal(mean = 0, sd = 10),
        phi = sv_beta(shape1 = 5, shape2 = 1.5),
        sigma2 = sv_gamma(shape = 0.5, rate = 10),
        nu = sv_infinity(),
        rho = sv_constant(0)
      )
    }
  }else{
    for(mm in 1:M){
      sv_priors[[mm]] <- specify_priors(
        mu = sv_constant(0),
        phi = sv_constant(1-1e-12),
        sigma2 = sv_constant(1e-12),
        nu = sv_infinity(),
        rho = sv_constant(0)
      )
    }
  }
  # priors for coefficients
  A_prior <- matrix(0,K,M)
  theta_A <- matrix(1,K,M)
  theta_A0 <- matrix(1,M,M)

  # # BART initialization
  # control <- dbartsControl(verbose = FALSE, keepTrainingFits = TRUE, useQuantiles = FALSE,
  #                          keepTrees = FALSE, n.samples = ntot,
  #                          n.cuts = 100L, n.burn = nburn, n.trees = num.trees, n.chains = 1,
  #                          n.threads = 1, n.thin = 1L, printEvery = 1,
  #                          printCutoffs = 0L, rngKind = "default", rngNormalKind = "default",
  #                          updateState = FALSE)

  opts <- Opts(update_sigma = TRUE, num_print = nburn + nsave + 1)


  hyperslist <- list()


  for(i in seq_len(M)){

    hyperslist[[i]] <- Hypers(as.matrix(X), as.vector(Y[,i]),
                              num_tree = num.trees, #sigma_hat = 1,
                              group = SB_group,
                              alpha = SB_alpha,
                              beta = SB_beta,
                              gamma = SB_gamma,
                              k = SB_k,
                              sigma_hat = NULL, #sighat,
                              shape = SB_shape,
                              width = SB_width,
                              # num_tree = 20,
                              alpha_scale = SB_alpha_scale,
                              alpha_shape_1 = SB_alpha_shape_1,
                              alpha_shape_2 = SB_alpha_shape_2,
                              tau_rate = SB_tau_rate,
                              num_tree_prob = SB_num_tree_prob,
                              temperature = SB_temperature,
                              weights = SB_weights,
                              normalize_Y = SB_normalize_Y
    )

  }





  sampler.list <- list()
  svdraw.list <- list()
  for (jj in seq_len(M)){
    # sampler.list[[jj]] <- dbarts(Y[,jj]~X, control = control,tree.prior = cgm(cgm.exp, cgm.level), node.prior = normal(sd.mu), n.samples = nsave, weights=rep(1,Tnobs),
    #                              sigma=sqrt(Sig_OLS[jj,jj])#, resid.prior = chisq(prior.sig[[1]], prior.sig[[2]])
    # )
    sampler.list[[jj]] <- MakeForest(hyperslist[[jj]], opts, warn = FALSE)


  }
  # sampler.run <- list()
  sigma.mat <- matrix(NA, M, 1)
  count.mat <- matrix(0, M*p, M)

  # Initialize HS prior on covariances -----------------------------------------------------------------------------
  # Initialize HS prior on covariances (if any)
  lambda.A0 <- 1
  nu.A0 <- 1
  tau.A0 <- 1
  zeta.A0 <- 1
  prior.cov <- rep(1, M*(M-1)/2)

  # storage objects
  Y_store <- array(NA,dim=c(nthin,Tnobs,M)) # filtered data
  fcst_store <- array(NA,dim=c(nthin,fhorz,M))

  # Start Gibbs Sampler -----------------------------------------------------------------------------
  # start Gibbs sampler
  # show progress
  if(!quiet){
    pb <- txtProgressBar(min = 0, max = ntot, style = 3)
    start  <- Sys.time()
  }


  weights.new.mat <- matrix(1/nrow(X), nrow = nrow(X), ncol = M)

  ##### begin loop over MCMC iterations ###########

  for(irep in 1:ntot){
    # 1) sample model coefficients (either linear VAR or BART)
    count.var <- 0
    var.check <- TRUE

    # stability check
    while(var.check){
      count.var <- count.var + 1
      X.ginv <- MASS::ginv(X)
      for (mm in seq_len(M)){
        if (mm > 1){
          Z_mm <- eta[,1:(mm-1), drop=F]
          A0_mm <- A0_draw[mm,1:(mm-1)]
          # sampler.list[[mm]]$setResponse(Y[,mm] - Z_mm%*%A0_mm)

          if(irep ==1){
            # sampler.list[[mm]]$setResponse(y = ymat[,mm] )
            tempy <- Y[,mm] #- (ymat[,-mm]-preds.train[,-mm])%*%(solve(Sigma_mat[-mm,-mm]) %*% Sigma_mat[-mm,mm]  )
          }else{
            tempy <- Y[,mm] - Z_mm%*%A0_mm

            # sampler.list[[mm]]$setResponse(y = ymat[,mm] - (ymat[,-mm]-preds.train[,-mm])%*%(solve(Sigma_mat[-mm,-mm]) %*% Sigma_mat[-mm,mm]  ) )
            # sampler.list[[jj]]$set_sigma(sqrt(Sigma_mat[mm,mm] - Sigma_mat[mm,-mm]%*%(solve(Sigma_mat[-mm,-mm]) %*% Sigma_mat[-mm,mm] ))  )
          }
        }else{
          tempy <- Y[,mm] #- (ymat[,-mm]-preds.train[,-mm])%*%(solve(Sigma_mat[-mm,-mm]) %*% Sigma_mat[-mm,mm]  )
        }


        # rep_mm <- sampler.list[[mm]]$run(0L, 1L) # construct BART sample using dbarts (V. Dorie)
        # mutrain <- sampler.list[[mm]]$do_gibbs(as.matrix(X), as.matrix(tempy), as.matrix(X), 1)

        mutrain <- t(sampler.list[[mm]]$do_gibbs_weighted(as.matrix(X),
                                             as.vector(tempy),
                                             weights.new.mat[,mm],
                                             as.matrix(X),
                                             1))

        # print("tempy = ")
        # print(tempy)
        #
        # print("mutrain = ")
        # print(mutrain)
        #
        # print("nrow(X) = ")
        # print(nrow(X))
        #
        # print("ncol(X) = ")
        # print(ncol(X))
        #
        # print("sampler.list[[mm]]$get_sigma()= ")
        # print(sampler.list[[mm]]$get_sigma())

        # sampler.run[[mm]] <- rep_mm
        sigma.mat[mm,] <- sampler.list[[mm]]$get_sigma() # rep_mm$sigma
        if (any(is.na(mutrain))){

          print("(tempmodel@tree.prior@splitProbabilities) = ")
          print((tempmodel@tree.prior@splitProbabilities))

          print("(s_y_list[[mm]]) = ")
          print((s_y_list[[mm]]))

          print("Y[,mm] - Z_mm%*%A0_mm = ")
          print(Y[,mm] - Z_mm%*%A0_mm)

          print("var_count_y_list[[mm]] = ")
          print(var_count_y_list[[mm]])

          print(" count.mat[,mm] = ")
          print( count.mat[,mm])

          # print("rep_mm$varcount = ")
          # print(rep_mm$varcount)


          stop("NA value in tree prediction sample")
        }
        if (any(is.na(sigma.mat[mm,] ))){
          stop("NA value in sigma sample")
        }
        eta[,mm] <- Y[,mm] - mutrain#rep_mm$train

        # print("mutrain = ")
        # print(mutrain)
        #
        # print("nrow(X.ginv) = ")
        # print(nrow(X.ginv))
        #
        # print("ncol(X.ginv) = ")
        # print(ncol(X.ginv))


        A_draw[,mm] <- X.ginv%*%mutrain#rep_mm$train
        count.mat[,mm] <- sampler.list[[mm]]$get_counts() # rep_mm$varcount






        if (mm > 1){
          norm_mm <- as.numeric(exp(-.5*sv_latent[[mm]]) * 1/sigma.mat[mm,])
          u_mm <- eta[,1:(mm-1),drop=F]*norm_mm
          eta_mm <- eta[,mm]*norm_mm
          if (mm == 2) v0.inv <- 1/theta_A0[mm,1] else v0.inv <- diag(1/theta_A0[mm,1:(mm-1)])
          V.cov <- solve(crossprod(u_mm) + v0.inv)
          mu.cov <- V.cov %*% crossprod(u_mm, eta_mm)
          mu.cov.draw <- mu.cov + t(chol(V.cov)) %*% rnorm(ncol(V.cov))
          A0_draw[mm,1:(mm-1)] <- mu.cov.draw
        }
      }

      shock_norm <- eta %*% t(solve(A0_draw))
      if(sv){
        for (mm in seq_len(M)){
          svdraw_mm <- svsample_general_cpp(shock_norm[,mm]/sigma.mat[mm], startpara = sv_draw[[mm]], startlatent = sv_latent[[mm]], priorspec = sv_priors[[mm]])
          sv_draw[[mm]][c("mu", "phi", "sigma")] <- as.list(svdraw_mm$para[, c("mu", "phi", "sigma")])
          sv_latent[[mm]] <- svdraw_mm$latent

          normalizer <- as.numeric(exp(-.5*svdraw_mm$latent))
          weights.new.mat[,mm] <- as.numeric(exp(-svdraw_mm$latent))

          # dat <- dbartsData(formula = Y[,mm]~X,weights=weights.new)
          # sampler.list[[mm]]$setData(dat)
          H[,mm] <- log(sigma.mat[mm]^2) + svdraw_mm$latent
        }
      }else{
        H[,mm] <- log(sigma.mat[mm]^2)
      }

      for(tt in seq_len(Tnobs)){
        S_tmp <- exp(H[tt,])
        S.t <- t(A0_draw)%*%crossprod(diag(S_tmp),(A0_draw))

        # if(class(try(solve(S.t),silent=T))!="matrix"){
        #
        #   print("(tempmodel@tree.prior@splitProbabilities) = ")
        #   print((tempmodel@tree.prior@splitProbabilities))
        #
        #   print("(s_y_list[[mm]]) = ")
        #   print((s_y_list[[mm]]))
        #
        #   print("Y[,mm] - Z_mm%*%A0_mm = ")
        #   print(Y[,mm] - Z_mm%*%A0_mm)
        #
        #   print("var_count_y_list[[mm]] = ")
        #   print(var_count_y_list[[mm]])
        #
        #   print(" count.mat[,mm] = ")
        #   print( count.mat[,mm])
        #
        #   print("rep_mm$varcount = ")
        #   print(rep_mm$varcount)
        #
        #   stop("S.t not invertible")
        # }

        Sig_t[tt,,] <- S.t
      }

      # 2) updating shrinkage priors
      # hierarchical prior values
      hs_draw <- get.hs(bdraw=A0_draw[lower.tri(A0_draw)],lambda.hs=lambda.A0,nu.hs=nu.A0,tau.hs=tau.A0,zeta.hs=zeta.A0)
      lambda.A0 <- hs_draw$lambda
      nu.A0 <- hs_draw$nu
      tau.A0 <- hs_draw$tau
      zeta.A0 <- hs_draw$zeta
      prior.cov <- hs_draw$psi
      theta_A0[lower.tri(theta_A0)] <- prior.cov
      theta_A0[theta_A0>10] <- 10
      theta_A0[theta_A0<1e-8] <- 1e-8

      # 3) updating latent states
      if(cons){
        Phi <- cbind(t(A_draw[c(1:(K-1),K),]))
      }else{
        Phi <- cbind(0,t(A_draw))
      }
      Sigma_chol <- array(NA,dim=c(M,M,Tnobs))
      for(tt in 1:Tnobs){
        Sigma_chol[,,tt] <- t(chol(Sig_t[tt,,]))
      }

      # adaptive simulation smoother as implemented in mfbvar (by S. Ankargren)
      beta2 <- mfbvar:::rsimsm_adaptive_sv(y_=Yact,Phi=Phi,Sigma_chol=Sigma_chol,Lambda=Lambda,Z1=Yinit,n_q_=M_l,T_b=T_b)

      if(M_l==1){
        var.check <- !(sd(beta2[,(M_h+1):M]) < var.thrsh)
      }else{
        var.check <- !any(apply(beta2[,(M_h+1):M],2,sd) < var.thrsh)
      }
      if (count.var == max.count.var){
        message("Reached maximum amount of replications.")
        count.var <- FALSE
        var.check <- FALSE
      }
    }

    # create new data matrices
    YY <- rbind(Yinit,beta2[,1:M])
    if(cons){
      X <- cbind(mlag(YY,p),1)[(p+1):nrow(YY),]
    }else{
      X <- mlag(YY,p)[(p+1):nrow(YY),]
    }
    Y <- YY[(p+1):nrow(YY),]
    # for(mm in 1:M){
    #   sampler.list[[mm]]$setPredictor(X) # set predictors in BART to new latent X matrix
    # }

    if(irep %in% thin.set){
      in.thin <- in.thin+1
      if(exact){
        for(mm in seq_len(M)){
          if(mm>M_h){
            # rep_mm <- sampler.run[[mm]]
            Y_store[in.thin,,mm] <- (mutrain + sigma.mat[mm,]*rnorm(Tnobs))*Ysd[mm] + Ymu[mm]
          }else{
            Y_store[in.thin,,mm] <- (beta2[,mm]*Ysd[mm]) + Ymu[mm]
          }
        }
      }else{
        Y_store[in.thin,,] <- (beta2*t(matrix(Ysd,M,Tnobs)))+t(matrix(Ymu,M,Tnobs))
      }


    }

    Yfc <- matrix(NA,fhorz,M)
    if(fhorz>0){
      if (cons){
        X.hat <- c(Y[Tnobs,],X[Tnobs,1:(M*(p-1))],1)
      }else{
        X.hat <- c(Y[Tnobs,],X[Tnobs,1:(M*(p-1))])
      }

      Sig_T <- Sig_t[Tnobs,,] # use final observation for Sigma
      tree.pred <- matrix(0, M)
      for (hh in seq_len(fhorz)){
        for (j in seq_len(M)){
          if(! is.matrix(X.hat)){
            tree.pred[j] <- sampler.list[[j]]$do_predict(as.matrix(t(X.hat)))
          }else{
            tree.pred[j] <- sampler.list[[j]]$do_predict(X.hat)
          }
        }
        Y.tp1 <- as.numeric(tree.pred) + t(chol(Sig_T))%*%rnorm(M)

        if (cons){
          X.hat <- c(Y.tp1, X.hat[1:(M*(p-1))],1)
        }else{
          X.hat <- c(Y.tp1, X.hat[1:(M*(p-1))])
        }
        Yfc[hh,] <- Y.tp1
      }
      fcst_store[in.thin,,] <- (Yfc*t(matrix(Ysd,M,fhorz)))+t(matrix(Ymu,M,fhorz))
    }
    if(!quiet) setTxtProgressBar(pb, irep)
  }

  dimnames(Y_store) <- list(paste0("mcmc",1:nthin),dates,colnames(Y))
  if(fhorz>0){
    dimnames(fcst_store) <- list(paste0("mcmc",1:nthin),dates_fc,colnames(Y))
    Yf_store <- abind(Y_store,fcst_store,along=2) # bind in-sample and out of sample values
  }else{
    Yf_store <- Y_store
  }

  # get quarterly aggregates
  Yq_store <- Yf_store*NA
  for(irep in 1:nthin){
    for(ii in 1:M){
      if(ii <= M_h){
        Yq_store[irep,,ii] <- Yf_store[irep,,ii]
      }else{
        if(itr[ii-M_h]=="lvl"){
          for(tt in p:(Tnobs+fhorz)){
            Yq_store[irep,tt,ii] <- sum(c(Yf_store[irep,tt,ii],Yf_store[irep,tt-1,ii],Yf_store[irep,tt-2,ii])*Loadings[ii-M_h,])
          }
        }else if(itr[ii-M_h]=="grw"){
          for(tt in p:(Tnobs+fhorz)){
            Yq_store[irep,tt,ii] <- sum(c(Yf_store[irep,tt,ii],Yf_store[irep,tt-1,ii],Yf_store[irep,tt-2,ii],Yf_store[irep,tt-3,ii],Yf_store[irep,tt-4,ii])*Loadings[ii-M_h,])
          }
        }
      }
    }
  }
  Yq_store <- Yq_store[,(p+1):(Tnobs+fhorz),]
  return_obj <- list("Y"=Y_store,"fcst"=fcst_store,"Yq"=Yq_store)
  return(return_obj)
}
