# "Nowcasting in a Pandemic using Non-Parametric Mixed Frequency VARs"
# Huber, F., Koop, G., Onorante, L., Pfarrhofer, M., and J. Schreiner,
# Journal of Econometrics (forthcoming)
#
# This file creates a function mfbavart(...) to estimate the MF-BAVART model.
# In addition to the baseline model in the paper, the code also includes an option
# to introduce stochastic volatility in the error terms.


#' @title Mixed Frequency Bayesian Additive Vector Autoregression Trees (Huber et al. 2023)
#'
#' @description Mixed Frequency Bayesian Additive Vector Autoregression Trees (Huber et al. 2023).
#' @import dbarts
#' @import MASS
#' @import stochvol
#' @import mfbvar
#' @import abind
#' @import collapse
#' @param data list of data frames for each variable.
#' @param ite option for constructing the loading matrix. Can be "m", "lvl", or "grw".
#' @param p number of lags
#' @param fhorz number of forecast horizons
#' @useDynLib BAVART, .registration = TRUE
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
mfbavart <- function(data,itr,p=5,fhorz=0,cons=FALSE,exact=FALSE,sv=FALSE,var.thrsh=10,max.count.var=10,
                     cgm.level=0.95,cgm.exp=2,sd.mu=2,num.trees=250,#prior.sig,
                     nburn=1000,nsave=1000,thinfac=1,
                     quiet=FALSE,
                     sparse = FALSE,
                     alpha_a_y = 0.5,
                     alpha_b_y = 1,
                     alpha_split_prior = TRUE){


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

  # BART initialization
  control <- dbartsControl(verbose = FALSE, keepTrainingFits = TRUE, useQuantiles = FALSE,
                           keepTrees = FALSE, n.samples = ntot,
                           n.cuts = 100L, n.burn = nburn, n.trees = num.trees, n.chains = 1,
                           n.threads = 1, n.thin = 1L, printEvery = 1,
                           printCutoffs = 0L, rngKind = "default", rngNormalKind = "default",
                           updateState = FALSE)

  if(sparse){

    # s_y_mat <- matrix(NA, nrow = ncol(x.train), ncol = M)
    p_y_vec <- rep(NA, M)
    rho_y_vec <- rep(NA, M)
    alpha_s_y_vec <- rep(NA, M)
    alpha_scale_y_vec <- rep(NA, M)
    # var_count_y_mat <- matrix(NA, nrow = ncol(x.train), ncol = M)

    s_y_list <- list() # matrix(NA, nrow = ncol(x.train), ncol = M)
    var_count_y_list <- list() # matrix(NA, nrow = ncol(x.train), ncol = M)

    for (jj in 1:M){

      p_y <- ncol(X)

      s_y <- rep(1 / p_y, p_y) # probability vector to be used during the growing process for DART feature weighting
      rho_y <- p_y # For DART

      if(alpha_split_prior){
        alpha_s_y <- p_y
      }else{
        alpha_s_y <- 1
      }
      alpha_scale_y <- p_y

      var_count_y <- rep(0, p_y)


      # s_y_mat[, jj] <- s_y
      p_y_vec[jj] <- p_y
      rho_y_vec[jj] <- rho_y
      alpha_s_y_vec[jj] <- alpha_s_y
      alpha_scale_y_vec[jj] <- alpha_scale_y
      # var_count_y_mat[,jj] <- var_count_y

      s_y_list[[jj]] <- s_y
      var_count_y_list[[jj]] <- var_count_y

    }

    alpha_s_y_store_arr <- matrix(NA, nrow = M, ncol =  nthin )
    # var_count_y_store_arr <- array(NA, dim = c( p_y, M, irep ))
    # s_prob_y_store_arr <- array(NA, dim = c( p_y, M, irep ))

    var_count_y_store_list <- list()
    s_prob_y_store_list <- list()

  }



  sampler.list <- list()
  svdraw.list <- list()
  for (jj in seq_len(M)){
    sampler.list[[jj]] <- dbarts(Y[,jj]~X, control = control,tree.prior = cgm(cgm.exp, cgm.level), node.prior = normal(sd.mu), n.samples = nsave, weights=rep(1,Tnobs),
                                 sigma=sqrt(Sig_OLS[jj,jj])#, resid.prior = chisq(prior.sig[[1]], prior.sig[[2]])
                                 )


    if(sparse){
      tempmodel <- sampler.list[[jj]]$model
      # print("length(tempmodel@tree.prior@splitProbabilities) = ")
      # print(length(tempmodel@tree.prior@splitProbabilities))
      #
      # print("ncol(sampler.list[[jj]]$data@x) = ")
      # print(ncol(sampler.list[[jj]]$data@x))
      #
      # print("length(s_y_list[[mm]]) = ")
      # print(length(s_y_list[[mm]]))
      #
      # print("(tempmodel@tree.prior@splitProbabilities) = ")
      # print((tempmodel@tree.prior@splitProbabilities))
      #
      # print("(s_y_list[[mm]]) = ")
      # print((s_y_list[[mm]]))


      tempmodel@tree.prior@splitProbabilities <- s_y_list[[jj]]
      sampler.list[[jj]]$setModel(newModel = tempmodel)
    }
  }
  sampler.run <- list()
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
          sampler.list[[mm]]$setResponse(Y[,mm] - Z_mm%*%A0_mm)
        }

        if(sparse){
          tempmodel <- sampler.list[[mm]]$model

          if(length(tempmodel@tree.prior@splitProbabilities)!= length(s_y_list[[mm]])){
            print("length(tempmodel@tree.prior@splitProbabilities) = ")
            print(length(tempmodel@tree.prior@splitProbabilities))

            print("length(s_y_list[[mm]]) = ")
            print(length(s_y_list[[mm]]))
            stop("error length of split sprobs")
          }
          tempmodel@tree.prior@splitProbabilities <- s_y_list[[mm]]
          sampler.list[[mm]]$setModel(newModel = tempmodel)
        }
        rep_mm <- sampler.list[[mm]]$run(0L, 1L) # construct BART sample using dbarts (V. Dorie)

        sampler.run[[mm]] <- rep_mm
        sigma.mat[mm,] <- rep_mm$sigma
        if (any(is.na(rep_mm$train))){

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

          print("rep_mm$varcount = ")
          print(rep_mm$varcount)


          stop("NA value in tree prediction sample")
        }
        if (any(is.na(rep_mm$sigma))){
          stop("NA value in sigma sample")
        }
        eta[,mm] <- Y[,mm] - rep_mm$train
        A_draw[,mm] <- X.ginv%*%rep_mm$train
        count.mat[,mm] <- rep_mm$varcount


        if(sparse){
          tempcounts <- fcount(sampler.list[[mm]]$getTrees()$var)
          tempcounts <- tempcounts[tempcounts$x != -1, ]
          var_count_y <- rep(0, p_y_vec[mm])
          var_count_y[tempcounts$x] <- tempcounts$N
          # var_count_y_mat[,mm] <- var_count_y

          var_count_y_list[[mm]] <- var_count_y
        }

        if (sparse & (irep > floor(nburn * 0.5))) {
          # s_update_z <- update_s(var_count_z, p_z, alpha_s_z)
          # s_z <- s_update_z[[1]]

          s_update_y <- update_s(var_count_y_list[[mm]], p_y_vec[mm], alpha_s_y_vec[mm])
          s_y_list[[mm]] <- s_update_y[[1]]

          if(alpha_split_prior){
            # alpha_s_z <- update_alpha(s_z, alpha_scale_z, alpha_a_z, alpha_b_z, p_z, s_update_z[[2]])
            alpha_s_y_vec[mm] <- update_alpha(s_y_list[[mm]], alpha_scale_y_vec[mm], alpha_a_y, alpha_b_y, p_y_vec[mm], s_update_y[[2]])
          }
        }

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
          weights.new <- as.numeric(exp(-svdraw_mm$latent))

          dat <- dbartsData(formula = Y[,mm]~X,weights=weights.new)
          sampler.list[[mm]]$setData(dat)
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
    for(mm in 1:M){
      sampler.list[[mm]]$setPredictor(X) # set predictors in BART to new latent X matrix
    }

    if(irep %in% thin.set){
      in.thin <- in.thin+1
      if(exact){
        for(mm in seq_len(M)){
          if(mm>M_h){
            rep_mm <- sampler.run[[mm]]
            Y_store[in.thin,,mm] <- (rep_mm$train + rep_mm$sigma*rnorm(Tnobs))*Ysd[mm] + Ymu[mm]
          }else{
            Y_store[in.thin,,mm] <- (beta2[,mm]*Ysd[mm]) + Ymu[mm]
          }
        }
      }else{
        Y_store[in.thin,,] <- (beta2*t(matrix(Ysd,M,Tnobs)))+t(matrix(Ymu,M,Tnobs))
      }

      if(sparse){
        alpha_s_y_store_arr[,in.thin] <- alpha_s_y_vec

        var_count_y_store_list[[in.thin]] <- var_count_y_list
        s_prob_y_store_list[[in.thin]] <- s_y_list
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
        for (j in seq_len(M)) tree.pred[j] <- sampler.list[[j]]$predict(X.hat)
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

  if(sparse){
    return_obj$alpha_s_y_store_arr <- alpha_s_y_store_arr
    return_obj$var_count_y_store_list <- var_count_y_store_list
    return_obj$s_prob_y_store_list <- s_prob_y_store_list
  }
  return(return_obj)
}
