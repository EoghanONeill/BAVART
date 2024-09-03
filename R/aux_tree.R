
update_s <- function(var_count, p, alpha_s) {
  # s_ = rdirichlet(1, as.vector((alpha_s / p ) + var_count))

  # // Get shape vector
  # shape_up = alpha_s / p
  # print("var_count = ")
  # print(var_count)
  # print("alpha_s = ")
  # print(alpha_s)
  # print("p = ")
  # print(p)
  #
  # print("alpha_s / p = ")
  # print(alpha_s / p)

  shape_up = as.vector((alpha_s / p ) + var_count)

  # print("shape_up = ")
  # print(shape_up)

  # // Sample unnormalized s on the log scale
  templogs = rep(NA, p)
  for(i in 1:p) {
    templogs[i] = SoftBart:::rlgam(shape = shape_up[i])
  }

  if(any(templogs== -Inf)){
    print("alpha_s = ")
    print(alpha_s)
    print("var_count = ")
    print(var_count)
    print("templogs = ")
    print(templogs)
    stop('templogs == -Inf')
  }

  # // Normalize s on the log scale, then exponentiate
  # templogs = templogs - log_sum_exp(hypers.logs);
  max_log = max(templogs)
  templogs2 = templogs - (max_log + log(sum(exp( templogs  -  max_log ))))


  s_ = exp(templogs2)

  # if(any(s_==0)){
  #   print("templogs2 = ")
  #   print(templogs2)
  #   print("templogs = ")
  #   print(templogs)
  #   print("alpha_s = ")
  #   print(alpha_s)
  #   print("var_count = ")
  #   print(var_count)
  #   print("s_ = ")
  #   print(s_)
  #   stop('s_ == 0')
  # }

  ret_list <- list()
  ret_list[[1]] <- s_
  ret_list[[2]] <- mean(templogs2)


  return(ret_list)
}



update_alpha <- function(s, alpha_scale, alpha_a, alpha_b, p, mean_log_s) {

  # create inputs for likelihood

  # log_s <- log(s)
  # mean_log_s <- mean(log_s)
  # p <- length(s)
  # alpha_scale   # denoted by lambda_a in JRSSB paper

  rho_grid <- (1:1000)/1001

  alpha_grid <- alpha_scale * rho_grid / (1 - rho_grid )

  logliks <- alpha_grid * mean_log_s +
    lgamma(alpha_grid) -
    p*lgamma(alpha_grid/p) + # (alpha_a - 1)*log(rho_grid) + (alpha_b-1)*log(1- rho_grid)
    dbeta(x = rho_grid, shape1 = alpha_a, shape2 = alpha_b, ncp = 0, log = TRUE)


  # logliks <- log(ddirichlet( t(matrix(s, p, 1000))  , t(matrix( rep(alpha_grid/p,p) , p , 1000)  ) ) ) +
  #   (alpha_a - 1)*log(rho_grid) + (alpha_b-1)*log(1- rho_grid)
  # # dbeta(x = rho_grid, shape1 = alpha_a, shape2 = alpha_b, ncp = 0, log = TRUE)

  # logliks <- rep(NA, 1000)
  # for(i in 1:1000){
  #   logliks[i] <- log(ddirichlet(s  , rep(alpha_grid[i]/p,p) ) ) +
  #     (alpha_a - 1)*log(rho_grid[i]) + (alpha_b-1)*log(1- rho_grid[i])
  # }

  max_ll <- max(logliks)
  logsumexps <- max_ll + log(sum(exp( logliks  -  max_ll )))

  # print("logsumexps = ")
  # print(logsumexps)

  logliks <- exp(logliks - logsumexps)

  if(any(is.na(logliks))){
    print("logliks = ")
    print(logliks)

    print("logsumexps = ")
    print(logsumexps)

    print("mean_log_s = ")
    print(mean_log_s)

    print("lgamma(alpha_grid) = ")
    print(lgamma(alpha_grid))

    print("p*lgamma(alpha_grid/p) = ")
    print(p*lgamma(alpha_grid/p))

    print("(alpha_a - 1)*log(rho_grid) + (alpha_b-1)*log(1- rho_grid) = ")
    print((alpha_a - 1)*log(rho_grid) + (alpha_b-1)*log(1- rho_grid))

    print("max_ll = ")
    print(max_ll)

    # print("s = ")
    # print(s)


  }

  # print("logliks = ")
  # print(logliks)

  rho_ind <- sample.int(1000,size = 1, prob = logliks)


  return(alpha_grid[rho_ind])
}

