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
#' @import GIGrvg
#' @import quantreg
#' @import statmod
#' @param data list of data frames for each variable.
#' @param ite option for constructing the loading matrix. Can be "m", "lvl", or "grw".
#' @param p number of lags
#' @param fhorz number of forecast horizons
#' @return The following objects are returned:
#' \item{Y}{Array of in-sample posterior draws. First dimension is number of MCMC iterations, second dimension is number of training time periods, third dimension is number of variables (length of input data list).}
#' \item{fcst}{Array of out-of-sample posterior draws (including noise form error term). First dimension is number of MCMC iterations, second dimension is number of horecast horizons, third dimension is number of variables.}
#' \item{Yq}{Array of out-of-sample quarterly aggregates. First dimension is number of MCMC iterations, second dimension is number of horecast horizons, third dimension is number of variables.}
#' @examples
#' run <- 1
#'
#' library(dbarts)
#' library(GIGrvg)
#' library(stochvol)
#' library(quantreg)
#' library(statmod)
#'
#' load("data/data_raw.rda")
#'
#' silent <- FALSE
#' sl.norm <- TRUE
#' length.hold.out <- 130
#' sl.cN <- c("DE","FR","IT","UK","US","AT","DK","ES","FI","NL","SE")
#'
#' fm.grid <- c(FALSE,TRUE)
#' horz.grid <- 1:4
#' weights.sample <- c("parm","non-parm","mix")
#'
#' h.out <- seq(1, length.hold.out)
#' mod.grid <- c("AR-abg","AR-abg-m")
#' pool.grid <- c(FALSE,TRUE)
#'
#' # set up grids
#' grid0 <- expand.grid("h"=h.out,
#'                      "fhorz"=horz.grid,
#'                      "fm"=FALSE,
#'                      "wghs"="parm",
#'                      "mod"="AR-abg-fq",
#'                      "pool"=FALSE,
#'                      "ald"=c("fixed"),
#'                      "bart"=c("v1"),
#'                      stringsAsFactors = FALSE)
#' grid1 <- expand.grid("h"=h.out,
#'                      "fhorz"=horz.grid,
#'                      "fm"=fm.grid,
#'                      "wghs"=weights.sample,
#'                      "mod"=mod.grid,
#'                      "pool"=pool.grid,
#'                      "ald"="regular",
#'                      "bart"=c("v1"),
#'                      stringsAsFactors = FALSE)
#'
#' grid <- rbind(grid0,grid1)
#' rownames(grid) <- 1:nrow(grid)
#'
#' # select corresponding specification
#' grid.slct <- grid[run,]
#' h <- as.numeric(grid.slct[["h"]])
#' fhorz <- as.numeric(grid.slct[["fhorz"]])
#' fm <- as.logical(grid.slct[["fm"]])
#' wghs <- grid.slct[["wghs"]]
#' mod <- grid.slct[["mod"]]
#' pool <- grid.slct[["pool"]]
#'
#' # prior setup for BART and ALD
#' ald <- grid.slct[["ald"]]
#' bart.v <- grid.slct[["bart"]]
#'
#' # storage of results
#' foldername <- paste0("results")
#' dir.create(foldername, showWarnings = FALSE)
#' filen <- paste0(formatC(h,flag="0",width=3),"_h",fhorz,
#'                 "_fm",fm,"_w",wghs,"_mod",mod,"_pool",pool,
#'                 "_ald",ald,"_bart",bart.v)
#' filename1 <- paste0(foldername,"/",filen,".rda")
#'
#' # ---------------------------------------------------------------------------------------------------
#' # model selection
#' P <- 1 # lags
#' sv <- TRUE # stochastic volatility
#' set.p <- c(0.05,0.1,0.16,0.25,0.4,0.5,0.6,0.75,0.84,0.9,0.95) # quantile function argument # seq(0.05,0.95,by=0.05)
#' R <- length(set.p) # Number of latent factors = one factor per quantile (cross-section information)
#'
#' # prior settings
#' fac.var <- 1
#' beta.pool.tau <- 10
#' B_h <- 1 # prior on state innovations in H
#'
#' # settings for ALD scale
#' if(ald=="regular"){
#'   ald.scale <- TRUE
#'   n0 <- rep(1,R)
#'   s0 <- rep(1,R)
#' }else if(ald=="fixed"){
#'   ald.scale <- TRUE
#'   n0 <- rep(0,R)
#'   s0 <- rep(0,R)
#' }
#'
#' # set hyperparameters for BART
#' cgm.level <- 0.95 # default value
#' cgm.exp <- 2 # default value
#'
#' if(bart.v=="v1"){
#'   sd.mu <- 2 # default
#'   num.trees <- 250 # default
#'
#'   beta.setup <- "hs"
#' }
#'
#' # ---------------------------------------------------------------------------------------------------
#' # compute design matrix
#' design_ls <- design.matrix(data_raw=data_raw,sl.cN=sl.cN,P=P,fhorz=fhorz,mod=mod,
#'                            sl.norm=sl.norm,
#'                            length.hold.out=length.hold.out,h=h)
#' list2env(design_ls,envir=globalenv())
#' mod <- "AR-abg"
#' wghs <- "non-parm"
#' # construct prior settings object in list
#' prior_obj <- list("nsave"=5000,"nburn"=5000,"silent"=silent,
#'                   "mod"=mod,"fm"=fm,"wghs"=wghs,"pool"=pool,
#'                   "set.p"=set.p,"R"=length(set.p),"sv"=sv,
#'                   fac.var=fac.var,"n0"=n0,"s0"=s0,"ald.scale"=ald.scale,"B_h"=B_h,
#'                   "cgm.level"=cgm.level,"cgm.exp"=cgm.exp,"sd.mu"=sd.mu,"num.trees"=num.trees,
#'                   "beta.setup"=beta.setup,"beta.pool.tau"=beta.pool.tau
#' )
#'
#' if(!file.exists(filename1)){
#'   ret.list <- qfbart(Y=Y,X=X,sl.X=sl.X,X.out=X.out,train.start=train.start,
#'                      Ymu=Ymu,Ysd=Ysd,prior_obj=prior_obj)
#'   save(file=filename1, ret.list)
#' }
#' @export
qfbart <- function(Y,X,sl.X,X.out,train.start,
                   Ymu,Ysd,prior_obj,
                   sparse = FALSE,
                   alpha_a_y = 0.5,
                   alpha_b_y = 1,
                   alpha_split_prior = TRUE){
  list2env(prior_obj,envir=globalenv())
  ntot <- nsave+nburn

  Tnum <- NROW(Y)
  K <- ncol(sl.X)
  mcmc_output <- TRUE

  # frequentist variant of the model
  if(mod == "AR-abg-fq"){
    pe.quants <- array(NA,dim=c(4,length(set.p),ncol(Y)))
    quant.mean <- array(NA,dim=c(4,length(set.p),Tnum,ncol(Y)))
    dimnames(pe.quants) <- list(c("low","med","high","mean"),set.p,colnames(Y))
    dimnames(quant.mean) <- list(c("low","med","high","mean"),set.p,NULL,colnames(Y))
    for(i in 1:N){
      Y.i <- Y[,i]
      X.i <- X[,sl.X[i,],drop=F]
      X.out.i <- X.out[sl.X[i,]]
      sim <- rq(Y.i~X.i-1,tau=set.p)

      for(mm in 1:4){
        pe.quants[mm,,i] <- (apply(as.matrix(sim$coefficients)*X.out.i,2,sum) * Ysd[i]) + Ymu[i]
        quant.mean[mm,,,i] <- (apply(sim$fitted.values,1,sort) * Ysd[i]) + Ymu[i]
      }
    }

    ret.list <- list("quants"=quant.mean,            # insample fitted quantiles
                     "predq"=pe.quants               # quantile estimates of quantile forecast
    )
  }

  # bayesian implementation and nested models
  if(mod != "AR-abg-fq"){
    # prior on the weights
    M <- length(set.p)
    b <- rep(1, M)
    a <- rep(1, M)

    if (as.character(wghs) == "non-parm"){
      est.weights <- FALSE
      omega.start <- 1
    }else if (as.character(wghs) == "parm"){
      est.weights <- FALSE
      omega.start <- 0
    }else if(as.character(wghs) == "mix"){
      est.weights <- FALSE
      omega.start <- 0
    }
    omega <- matrix(omega.start, M, N)

    # factor model specification
    if (R>0 & fm){
      f <- matrix(rnorm(R*Tnum, 0, 0.01), Tnum, R)
      Lambda <- matrix(0, N*M, R)
    }else{
      f <- matrix(0, Tnum, R)
      Lambda <- matrix(0, N*M, R)
    }

    list.tree.eq <- list()
    list.sampler.run  <- list()

    # pooling
    sl.dom <- matrix(FALSE,K,N)
    for(i in 1:N){
      sl.dom[,i] <- grepl(paste(c(substr(colnames(Y)[i],1,3),"cons"),collapse="|"),colnames(X)[sl.X[i,]])
    }
    K.dom <- unique(apply(sl.dom,2,sum))

    # initialize BART
    control <- dbartsControl(verbose = FALSE, keepTrainingFits = TRUE, useQuantiles = FALSE,
                             keepTrees = FALSE, n.samples = 1,
                             n.cuts = 100L, n.burn = 0, n.trees = num.trees, n.chains = 1,
                             n.threads = 1, n.thin = 1L, printEvery = 1,
                             printCutoffs = 0L, rngKind = "default", rngNormalKind = "default",
                             updateState = FALSE)

    if(sparse){

      list_p_y_vec <- list()
      list_rho_y_vec <- list()
      list_alpha_s_y_vec <- list()
      list_alpha_scale_y_vec <- list()

      list_s_y_list <- list()
      list_var_count_y_list<- list()


      for (ii in seq_len(N)){

        # s_y_mat <- matrix(NA, nrow = ncol(x.train), ncol = M)
        p_y_vec <- rep(NA, M)
        rho_y_vec <- rep(NA, M)
        alpha_s_y_vec <- rep(NA, M)
        alpha_scale_y_vec <- rep(NA, M)
        # var_count_y_mat <- matrix(NA, nrow = ncol(x.train), ncol = M)

        s_y_list <- list() # matrix(NA, nrow = ncol(x.train), ncol = M)
        var_count_y_list <- list() # matrix(NA, nrow = ncol(x.train), ncol = M)


        for (jj in seq_len(M)){

          p_y <- ncol(X[,sl.X[ii,]])

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

        list_p_y_vec[[ii]] <- p_y_vec
        list_rho_y_vec[[ii]] <- rho_y_vec
        list_alpha_s_y_vec[[ii]] <- alpha_s_y_vec
        list_alpha_scale_y_vec[[ii]] <- alpha_scale_y_vec

        list_s_y_list[[ii]] <- s_y_list
        list_var_count_y_list[[ii]] <- var_count_y_list
      }
      alpha_s_y_store_list <- list()
      # var_count_y_store_arr <- array(NA, dim = c( p_y, M, irep ))
      # s_prob_y_store_arr <- array(NA, dim = c( p_y, M, irep ))

      var_count_y_store_list <- list()
      s_prob_y_store_list <- list()

    }


    for (ii in seq_len(N)){
      sampler.list <- list()
      sampler.run <- list()

      # s_y_list <- list_s_y_list[[ii]]

      # p_y_vec  <- list_p_y_vec[[ii]]
      # rho_y_vec  <- list_rho_y_vec[[ii]]
      # alpha_s_y_vec  <- list_alpha_s_y_vec[[ii]]
      # alpha_scale_y_vec  <- list_alpha_scale_y_vec[[ii]]
      #
      # list_s_y_list[[ii]] <- s_y_list
      # list_var_count_y_list[[ii]] <- var_count_y_list

      for (jj in seq_len(M)){
        prior.sig = c(10000^50, 0.5)
        sigma.init <- 1

        sampler.list[[jj]] <- dbarts(Y[,ii]~X[,sl.X[ii,]],
                                     control = control,
                                     tree.prior = cgm(cgm.exp, cgm.level),
                                     node.prior = normal(sd.mu), n.samples = nsave, weights=rep(1,Tnum),
                                     sigma=sigma.init,
                                     resid.prior = chisq(prior.sig[[1]], prior.sig[[2]]))

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
          # print("Line 332. ii = ")
          # print(ii)
          #
          # print("jj = ")
          # print(jj)
          #
          # print("ncol(X[,sl.X[ii,]]) = ")
          # print(ncol(X[,sl.X[ii,]]))
          #
          # print("tempmodel = ")
          # print(tempmodel)
          #
          # print("(tempmodel@tree.prior@splitProbabilities) = ")
          # print((tempmodel@tree.prior@splitProbabilities))
          #
          # print("list_s_y_list[[ii]][[jj]] = ")
          # print(list_s_y_list[[ii]][[jj]])
          #
          # print("tempmodel = ")
          # print(tempmodel)
          #
          # print("ncol(sampler.list[[jj]]$data@x)")
          # print(ncol(sampler.list[[jj]]$data@x))
          #
          # print("X[,sl.X[ii,]] = ")
          # print(X[,sl.X[ii,]])


          p_y <- ncol(sampler.list[[jj]]$data@x)

          s_y <- rep(1 / p_y, p_y) # probability vector to be used during the growing process for DART feature weighting
          rho_y <- p_y # For DART

          if(alpha_split_prior){
            alpha_s_y <- p_y
          }else{
            alpha_s_y <- 1
          }
          alpha_scale_y <- p_y

          var_count_y <- rep(0, p_y)


          list_p_y_vec[[ii]][jj] <- p_y
          list_rho_y_vec[[ii]][jj] <- rho_y
          list_alpha_s_y_vec[[ii]][jj] <- alpha_s_y
          list_alpha_scale_y_vec[[ii]][jj] <- alpha_scale_y

          list_s_y_list[[ii]][[jj]] <- s_y
          list_var_count_y_list[[ii]][[jj]] <- var_count_y


          tempmodel@tree.prior@splitProbabilities <- list_s_y_list[[ii]][[jj]]
          sampler.list[[jj]]$setModel(newModel = tempmodel)

          # print("line 359")
        }
      }
      list.tree.eq[[ii]] <- sampler.list
      list.sampler.run[[ii]] <- sampler.run
    }

    # initialize stochvol and set priors
    if (sv){
      sv.draw.idio <- list()
      sv.latent.idio <- list()

      sv.draw.fac <- list()
      sv.latent.fac <- list()
      for (mm in seq_len(R)){
        sv.draw.fac[[mm]] <- list(mu = 0, phi = 0.95, sigma = 0.01, nu = Inf, rho = 0, beta = NA, latent0 = 0)
        sv.latent.fac[[mm]] <- rep(0,Tnum)
      }

      sv_priors <- list()
      for(mm in 1:R){
        sv_priors[[mm]] <- specify_priors(
          mu = sv_constant(0), # -4, 1e-3
          phi = sv_beta(shape1 = 5, shape2 = 1.5),
          sigma2 = sv_gamma(shape = 0.5, rate = 1/(2*B_h)),
          nu = sv_infinity(),
          rho = sv_constant(0)
        )
      }
    }
    H <- matrix(0, Tnum, R)

    beta.mat <- array(NA, c(K, M, N))
    Vprior.array <- array(1, c(K, M, N))
    for (ii in seq_len(N)){
      for (jj in seq_len(M)){
        beta.mat[,jj,ii] <- lm(Y[,ii]~X[,sl.X[ii,]]-1)$coefficients
      }
    }
    beta.prior <- beta.mat*0 # stays zero if not updated in mode pool=TRUE

    # compute theta and tau^2  (this is the same across equations)
    theta <- (1-2*set.p)/(set.p * (1 - set.p))
    tau2 <- 2/(set.p * (1-set.p))
    sigma <- matrix(1, M, N)

    # Starting values for auxiliary variable in Gaussian representation of ALD
    v <- array(1, c(Tnum, M, N))

    # storage
    omega.store <- array(NA, c(nsave, M,N))
    beta.store <- array(NA, c(nsave, K, M, N))
    quant.store <- array(NA, c(nsave, Tnum, M, N))
    f.store <- array(NA, c(nsave, Tnum, R))
    Lambda.store <- array(NA, c(nsave, M, R*N))
    pred.store <- array(NA, c(nsave,M,N))
    count.store <- array(NA, c(nsave, K-1, N, M))

    y.quantiles <- array(NA, c(Tnum, M, N))
    y.quantiles_parts <- array(NA, c(Tnum, M, N, 2))
    eta  <- array(NA, c(Tnum, M, N))

    #Starting values for the factors and loadings
    sigma2.draw <- rep(1, M*N)

    # starting values for the HS prior
    bigpi <- 1
    nu <- 1
    xi <- 1
    zeta <- 1

    if(beta.setup=="flat"){
      V.prior.mat <- array(10^10, c(K, M,N))
    }else if(beta.setup=="tight"){
      V.prior.mat <- array(1, c(K, M,N))
    }else{
      V.prior.mat <- array(1, c(K, M,N))
    }
    ind.restr <- matrix(seq(N*M), N, M)

    # starting sampling algorithm ------------------------------------------------------------------------------
    # starting sampling algorithm
    pb <- txtProgressBar(min = 0, max = ntot, style = 3) #start progress bar
    for (irep in seq_len(ntot)){
      # This block samples the parameters of the equations associated with the different quantiles
      var.mat <- array(NA, c(Tnum, M, N))
      count.mat <- array(NA, c(K-1,M,N))

      for (i in seq_len(N)){


        # p_y_vec  <- list_p_y_vec[[i]]
        # rho_y_vec  <- list_rho_y_vec[[i]]
        # alpha_s_y_vec  <- list_alpha_s_y_vec[[i]]
        # alpha_scale_y_vec  <- list_alpha_scale_y_vec[[i]]
        #
        # s_y_list <- list_s_y_list[[i]]
        # var_count_y_list <- list_var_count_y_list[[i]]



        X.i <- X[, sl.X[i,],drop=F]
        Lambda.i <- Lambda[seq(i, N*M, by=N),]
        for (q in seq_len(M)){
          # These are quantities necessary to render the regression and BART part conditionally Gaussian and homoscedastic
          v.i <- v[,q,i]
          var.i <- tau2[q]*sigma[q,i]*v.i
          var.mat[,q,i] <- var.i

          if(as.character(wghs) == "mix"){
            ytildeq <- (Y[,i] - (X.i%*%beta.mat[, q,i]) - theta[q]*v.i - (f%*%t(Lambda.i))[,q])
            var.qt <- var.i

            # We need to re-define the response and the scaling parameter for BART to work
            list.tree.eq[[i]][[q]]$setResponse(ytildeq)
            list.tree.eq[[i]][[q]]$setWeights(1/var.qt)


            if(sparse){
              tempmodel <- list.tree.eq[[i]][[q]]$model
              # print("Line 466. i = ")
              # print(i)
              #
              # print("q = ")
              # print(q)
              #
              # print("(tempmodel@tree.prior@splitProbabilities) = ")
              # print((tempmodel@tree.prior@splitProbabilities))
              #
              # print("list_s_y_list[[ii]][[jj]] = ")
              # print(list_s_y_list[[ii]][[jj]])

              tempmodel@tree.prior@splitProbabilities <- list_s_y_list[[i]][[q]]
              list.tree.eq[[i]][[q]]$setModel(newModel = tempmodel)
            }


            rep_mm <- list.tree.eq[[i]][[q]]$run(0L, 1L)
            list.sampler.run[[i]][[q]] <- rep_mm

            count.mat[,q,i] <- t(rep_mm$varcount)/num.trees

            # print("Line 547 count.mat[,q,i] = ")
            # print(count.mat[,q,i])


            if(sparse){
              tempcounts <- fcount(list.tree.eq[[i]][[q]]$getTrees()$var)
              tempcounts <- tempcounts[tempcounts$x != -1, ]
              var_count_y <- rep(0, list_p_y_vec[[i]][q])
              var_count_y[tempcounts$x] <- tempcounts$N
              # var_count_y_mat[,mm] <- var_count_y

              list_var_count_y_list[[i]][[q]] <- var_count_y
            }

            if (sparse & (irep > floor(nburn * 0.5))) {
              # s_update_z <- update_s(var_count_z, p_z, alpha_s_z)
              # s_z <- s_update_z[[1]]

              s_update_y <- update_s(list_var_count_y_list[[i]][[q]], list_p_y_vec[[i]][q], list_alpha_s_y_vec[[i]][q])
              list_s_y_list[[i]][[q]] <- s_update_y[[1]]

              if(alpha_split_prior){
                # alpha_s_z <- update_alpha(s_z, alpha_scale_z, alpha_a_z, alpha_b_z, p_z, s_update_z[[2]])
                list_alpha_s_y_vec[[i]][q] <- update_alpha(list_s_y_list[[i]][[q]], list_alpha_scale_y_vec[[i]][q], alpha_a_y, alpha_b_y, list_p_y_vec[[i]][q], s_update_y[[2]])
              }
            }

            # Now, conditional on the tree sample the regression coefficients
            yhatq <- (Y[,i] - rep_mm$train - theta[q]*v.i  - (f%*%t(Lambda.i))[,q])/(sqrt(var.i))
            Xhatq <- X.i/(sqrt(var.i))

            if (K > 1) V.prior.inv <- diag(1/V.prior.mat[,q,i]) else V.prior.inv <- 1/V.prior.mat[,q,i]
            V.q <- solve(crossprod(Xhatq) + V.prior.inv)
            m.q <- V.q %*% (V.prior.inv %*% beta.prior[,q,i] + crossprod(Xhatq,yhatq))
            m.draw <- m.q + t(chol(V.q))%*%rnorm(K)
            beta.mat[,q,i] <- m.draw

            # Draw the parameters used for the Gauss approximation part
            d.q.2 <- as.numeric(sqrt(theta[q]^2 + 2*tau2[q]) / abs(Y[,i] - rep_mm$train - (X.i%*%m.draw) - (f%*%t(Lambda.i))[,q]))
            g.q.2 <- (theta[q]^2 + 2*tau2[q]) / (sigma[q,i] * tau2[q])
            v.i <- 1/rinvgauss(Tnum,mean=d.q.2,dispersion=1/g.q.2)
            v[,q,i] <- v.i

            # Draw sigma
            if(ald.scale){
              n.tilda <- (n0 + 3*Tnum)/2
              s.tilda <- (s0 + 2*sum(v.i) + sum((Y[,i] - rep_mm$train - (X.i%*%m.draw) - theta[q]*v.i - (f%*%t(Lambda.i))[,q])^2/(tau2[q] * v.i)))/2
              sigma[q,i] <- 1/rgamma(1, n.tilda, s.tilda)
            }else{
              sigma[q,i] <- 1
            }

            y.quantiles[ ,q,i] <- rep_mm$train + (X.i%*%m.draw) + (f%*%t(Lambda.i))[,q]
            y.quantiles_parts[,q,i,1] <- rep_mm$train
            y.quantiles_parts[,q,i,2] <- (X.i%*%m.draw)

            if (R > 0){
              # Construct shocks for estimating the latent factor
              eta[, q,i] <- Y[,i] - rep_mm$train - (X.i%*%m.draw)  - v.i*theta[q]
            }
          }else{
            if (omega[q,i]==1){
              ytildeq <- (Y[,i] - (X.i%*%beta.mat[, q,i]) * (1-omega[q,i]) - theta[q]*v.i - (f%*%t(Lambda.i))[,q])/omega[q,i]
              var.qt <- var.i/omega[q,i]^2

              # We need to re-define the response and the scaling parameter for BART to work
              list.tree.eq[[i]][[q]]$setResponse(ytildeq)
              list.tree.eq[[i]][[q]]$setWeights(1/var.qt)

              if(sparse){
                tempmodel <- list.tree.eq[[i]][[q]]$model

                tempmodel@tree.prior@splitProbabilities <- list_s_y_list[[i]][[q]]
                list.tree.eq[[i]][[q]]$setModel(newModel = tempmodel)
              }

              rep_mm <- list.tree.eq[[i]][[q]]$run(0L, 1L)
              list.sampler.run[[i]][[q]] <- rep_mm

              count.mat[,q,i] <- t(rep_mm$varcount)/num.trees
              m.draw <- rep(0,K)

              # print("Line 547 count.mat[,q,i] = ")
              # print(count.mat[,q,i])



              if(sparse){
                tempcounts <- fcount(list.tree.eq[[i]][[q]]$getTrees()$var)
                tempcounts <- tempcounts[tempcounts$x != -1, ]
                var_count_y <- rep(0, list_p_y_vec[[i]][q])
                var_count_y[tempcounts$x] <- tempcounts$N
                # var_count_y_mat[,mm] <- var_count_y

                list_var_count_y_list[[i]][[q]] <- var_count_y
              }

              if (sparse & (irep > floor(nburn * 0.5))) {
                # s_update_z <- update_s(var_count_z, p_z, alpha_s_z)
                # s_z <- s_update_z[[1]]

                s_update_y <- update_s(list_var_count_y_list[[i]][[q]], list_p_y_vec[[i]][q], list_alpha_s_y_vec[[i]][q])
                list_s_y_list[[i]][[q]] <- s_update_y[[1]]

                if(alpha_split_prior){
                  # alpha_s_z <- update_alpha(s_z, alpha_scale_z, alpha_a_z, alpha_b_z, p_z, s_update_z[[2]])
                  list_alpha_s_y_vec[[i]][q] <- update_alpha(list_s_y_list[[i]][[q]], list_alpha_scale_y_vec[[i]][q], alpha_a_y, alpha_b_y, list_p_y_vec[[i]][q], s_update_y[[2]])
                }
              }


            }else{
              # print("no trees")

              rep_mm <- list()
              rep_mm$train <- rep(0, Tnum)
            }

            if(omega[q,i]==0){
              # Now, conditional on the tree sample the regression coefficients
              yhatq <- (Y[,i] - rep_mm$train * omega[q,i] - theta[q]*v.i  - (f%*%t(Lambda.i))[,q])/(sqrt(var.i))
              Xhatq <- X.i * (1-omega[q,i])/(sqrt(var.i))

              if (K > 1) V.prior.inv <- diag(1/V.prior.mat[,q,i]) else V.prior.inv <- 1/V.prior.mat[,q,i]
              V.q <- solve(crossprod(Xhatq) + V.prior.inv)
              m.q <- V.q %*% (V.prior.inv %*% beta.prior[,q,i] + crossprod(Xhatq,yhatq))
              m.draw <- m.q + t(chol(V.q))%*%rnorm(K)
              beta.mat[,q,i] <- m.draw
            }

            # Draw the parameters used for the Gauss approximation part
            d.q.2 <- as.numeric(sqrt(theta[q]^2 + 2*tau2[q]) / abs(Y[,i] - rep_mm$train*omega[q,i] - (X.i%*%m.draw)*(1-omega[q,i]) - (f%*%t(Lambda.i))[,q]))
            g.q.2 <- (theta[q]^2 + 2*tau2[q]) / (sigma[q,i] * tau2[q])
            v.i <- 1/rinvgauss(Tnum,mean=d.q.2,dispersion=1/g.q.2)
            v[,q,i] <- v.i

            # Draw sigma
            if(ald.scale){
              n.tilda <- (n0 + 3*Tnum)/2
              s.tilda <- (s0 + 2*sum(v.i) + sum((Y[,i] - rep_mm$train*omega[q,i] - (X.i%*%m.draw)*(1-omega[q,i]) - theta[q]*v.i - (f%*%t(Lambda.i))[,q])^2/(tau2[q] * v.i)))/2
              sigma[q,i] <- 1/rgamma(1, n.tilda, s.tilda)
            }else{
              sigma[q,i] <- 1
            }

            y.quantiles[ ,q,i] <- rep_mm$train*omega[q,i] + (X.i%*%m.draw)*(1-omega[q,i]) + (f%*%t(Lambda.i))[,q]
            y.quantiles_parts[,q,i,1] <- rep_mm$train
            y.quantiles_parts[,q,i,2] <- (X.i%*%m.draw)

            if (R > 0){
              # Construct shocks for estimating the latent factor
              eta[, q,i] <- Y[,i] - rep_mm$train*omega[q,i] - (X.i%*%m.draw)*(1-omega[q,i])  - v.i*theta[q]
            }
          }
        }
      }


      if(sparse & (irep>nburn) ){
        alpha_s_y_store_list[[irep-nburn]] <- list_alpha_s_y_vec

        var_count_y_store_list[[irep-nburn]] <- list_var_count_y_list
        s_prob_y_store_list[[irep-nburn]] <- list_s_y_list
      }

      # Step II: Sample shrinkage parameters for the HS
      if(!(beta.setup %in% c("tight","flat"))){
        if(pool){
          for(q in 1:M){
            V.q <- matrix(V.prior.mat[,q,][sl.dom],K.dom,N)
            beta.q <- matrix(beta.mat[,q,][sl.dom],K.dom,N)

            V.q.inv <- diag(1/(apply(1/V.q,1,sum) + 1/beta.pool.tau))
            beta.q.mu <- V.q.inv %*% apply(beta.q/V.q,1,sum)
            beta.pool <- beta.q.mu + t(chol(V.q.inv)) %*% rnorm(K.dom,0,1)
            beta.prior[,q,][sl.dom] <- beta.pool
          }
        }

        hs_draw <- get.hs(bdraw=(as.vector(beta.mat)-as.vector(beta.prior)),lambda.hs=bigpi,nu.hs=nu,tau.hs=xi,zeta.hs=zeta)
        bigpi <- hs_draw$lambda
        nu <- hs_draw$nu
        xi <- hs_draw$tau
        zeta <- hs_draw$zeta
        V.prior.mat <- array(hs_draw$psi, c(K, M, N))
        V.prior.mat[V.prior.mat<1e-8] <- 1e-8
        V.prior.mat[V.prior.mat>10] <- 10
      }

      # Step III: Sample the common factor in included in X
      # Before we do this, create matrices eta.mat and Var.mat, we need to do this such that eta.mat = (eta_1p1, eta_2p1, .., eta_Np1, eta_1p2, eta_2p2, ..)
      eta.mat <- NULL
      Var.mat <- NULL
      for (jj in seq_len(M)){
        eta.mat <- cbind(eta.mat, eta[,jj,])
        Var.mat <- cbind(Var.mat, var.mat[,jj,])
      }

      if (R>0 & fm){
        for (t in seq_len(Tnum)){
          normalizer <- 1/sqrt(Var.mat[t, ])
          Lt <- Lambda * normalizer
          yt <- eta.mat[t, ] * normalizer

          if (R > 1) Q <- diag(1/exp(H[t,])) else Q <- 1/exp(H[t])
          fac.Sigma <-  solve(crossprod(Lt)+Q)
          fac.mean <- fac.Sigma%*%crossprod(Lt,yt)

          if (R > 1 ){
            f[t,] <- fac.mean+t(chol(fac.Sigma))%*%rnorm(R)
          } else {
            f[t] <-  fac.mean+t(chol(fac.Sigma))%*%rnorm(R)
          }
        }
        f <- apply(f,2,function(x){(x-mean(x))/sd(x)}) # normalize factor draw

        if (sv){
          # Step IV: Sample the factor volatilities
          sv.para <- matrix(NA, R, 3) # first dim = mu, second = rho, third = sigma
          for (jj in 1:R){
            svdraw_mm <- svsample_general_cpp(f[,jj], startpara = sv.draw.fac[[jj]], startlatent = sv.latent.fac[[jj]], priorspec = sv_priors[[jj]])
            sv.draw.fac[[jj]][c("mu", "phi", "sigma")] <- as.list(svdraw_mm$para[, c("mu", "phi", "sigma")])
            sv.latent.fac[[jj]] <- svdraw_mm$latent
            H[,jj] <- svdraw_mm$latent
            sv.para[jj,] <- svdraw_mm$para[, c("mu", "phi", "sigma")]
          }
        }else{
          H[] <- 0
        }

        for (j in seq_len(N*M)){
          sl.fac <- which(ind.restr==j,arr.ind = TRUE)[[2]]
          normalizer <- 1/sqrt(Var.mat[, j])
          yj <- eta.mat[,j]*normalizer
          fj <- f[,sl.fac]*normalizer
          prior.v <- fac.var

          V.f <- solve(crossprod(fj)+prior.v)
          m.lambda.f <- V.f %*% crossprod(fj, yj)
          lambda.draw <- m.lambda.f + t(chol(V.f))%*%rnorm(1)
          Lambda[j,sl.fac] <- lambda.draw
        }

        # identify sign of the factor
        for(rr in 1:R){
          L.tmp <- Lambda[,rr]
          L.sign <- sign(mean(L.tmp[L.tmp!=0]))
          Lambda[,rr] <- Lambda[,rr]*L.sign
          f[,rr] <- f[,rr]*L.sign
        }
      }

      pred.mat <- array(0, c(M,N)); dimnames(pred.mat) <- list(set.p,colnames(Y))
      if (irep > nburn){
        pred.mat <- array(NA, c(M,N)); dimnames(pred.mat) <- list(set.p,colnames(Y))

        # Shock DFM and predict log-volas
        if (fm){
          if (sv){
            HT1 <- sv.para[,1] + sv.para[,2]*(H[Tnum,] - sv.para[,1]) + sv.para[,3]*rnorm(R)
          }else{
            HT1 <- H[Tnum,]
          }
          HT1[exp(HT1/2)>20] <- 2*log(20) # offsetting for pandemic
          f.shock <- Lambda%*%rnorm(R, 0, exp(HT1/2))
        }else{
          f.shock <- Lambda%*%rep(0,R)
        }

        # use samples from holdout
        X.p <- X.out
        for (i in seq_len(N)){
          f.i <- f.shock[seq(i, N*M, by=N)]

          for (j in seq_len(M)){
            pred.tree <- list.tree.eq[[i]][[j]]$predict(X.p[sl.X[i,-ncol(sl.X)]])
            pred.reg <- X.p[sl.X[i,]]%*%beta.mat[,j,i]

            if(as.character(wghs) == "mix"){
              pred.t <- pred.tree + pred.reg + f.i[j]
            }else{
              pred.t <- pred.tree*omega[j,i] + pred.reg*(1-omega[j,i]) + f.i[j]
            }

            pred.mat[j,i] <- pred.t
          }
        }

        # rescaling
        for(i in 1:N){
          pred.mat[,i] <- (pred.mat[,i] * Ysd[i]) + Ymu[i]
          y.quantiles[,,i] <- (y.quantiles[,,i] * Ysd[i]) + Ymu[i]
        }

        # storage
        pred.store[irep-nburn,,] <- pred.mat

        if(as.character(wghs) == "mix"){
          var_quantiles <- apply(y.quantiles_parts,c(2,3,4),sd)^2
          omega.store[irep-nburn,,] <- var_quantiles[,,1] / apply(var_quantiles,c(1,2),sum)
        }else{
          omega.store[irep-nburn,,] <- omega
        }

        beta.store[irep-nburn,,,] <- beta.mat
        quant.store[irep-nburn,,,] <- y.quantiles
        f.store[irep-nburn,,] <- f
        Lambda.store[irep-nburn,,] <- Lambda
        count.store[irep-nburn,,,] <- count.mat
      }

      # progress bar
      if(!silent) setTxtProgressBar(pb, irep)
    }

    pe.quants <- apply(pred.store, c(2,3), quantmean.fun)
    dimnames(pe.quants) <- list(c("low", "med", "high", "mean"),
                                set.p,
                                colnames(Y))

    quant.mean <- apply(quant.store, c(2,3,4), quantmean.fun) # insample estimate of the quantiles
    dimnames(quant.mean) <- list(c("low","med","high","mean"),
                                 as.character(rownames(Y)),
                                 set.p,
                                 colnames(Y))

    weights.quant <- apply(omega.store,c(2,3),quantmean.fun)
    beta.quant <- apply(beta.store, c(2,3,4), quantmean.fun)
    f.quant <- apply(f.store, c(2,3), quantmean.fun)
    L.quant <- apply(Lambda.store, c(2,3), quantmean.fun)
    count.quant <- apply(count.store, c(2,3,4), quantmean.fun)

    dimnames(weights.quant) <- list(c("low","med","high","mean"),set.p,colnames(Y))
    dimnames(beta.quant) <- list(c("low","med","high","mean"),
                                 NULL,
                                 set.p,
                                 colnames(Y))
    dimnames(f.quant) <- list(c("low","med","high","mean"),
                              as.character(rownames(Y)),
                              set.p)

    ret.list <- list("weights"=weights.quant,        # weight on the non-parametric part
                     "beta"=beta.quant,              # "linear" part of the coefficients
                     "quants"=quant.mean,            # insample fitted quantiles
                     "predq"=pe.quants,              # quantile estimates of quantile forecast
                     "factors"=f.quant,              # latent factors wrt. cross-section
                     "loadings"=L.quant,             # loadings on the factors
                     "count"=count.quant             # splitting rule count
                     )

    if(sparse){
      ret.list$alpha_s_y_store_list <- alpha_s_y_store_list
      ret.list$var_count_y_store_list <- var_count_y_store_list
      ret.list$s_prob_y_store_list <- s_prob_y_store_list
    }

  }
  return(ret.list)
}
