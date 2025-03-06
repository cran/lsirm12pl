
check.datatype <- function(data, missing.val = 99, ...){
  if(exists("missing.val")){ # If missing.var is imputed
    missing.temp = missing.val
  }else{missing.temp = 99}

  check.v <- c(0, 1, missing.temp)
  check.r <- length(setdiff(unique(c(data.matrix(data))), check.v)) == 0 # If the unique values are exist except check.v
  return(check.r) # TRUE: Binary, FALSE: Continuous

  # if(!exists('missing_data')){
  #   missing.temp = c()
  # }else{
  #   missing.temp = missing.var
  # }
  # check.v <- c(0,1,missing.temp)
  # check.r <- sum(check.v %in% unique(c(data.matrix(data)))) == length(check.v)  # If data have (0,1,90), cannot catch continuous
  # return(check.r)
}

# Distance
distance <- function(x, y){
  mat <- matrix(NA, nrow=nrow(x), ncol=nrow(y))

  for(i in 1:nrow(x)){
    for(j in 1:nrow(y)){
      mat[i,j] <- sqrt((x[i,1]-y[j,1])^2 + (x[i,2]-y[j,2])^2)
    }
  }

  return(mat)
}

# Clustering
clustering <- function(data, center){
  mat <- distance(data, center)

  distance <- c()
  group <- c()

  for(i in 1:nrow(data)){
    distance[i] <- min(mat[i,]) #min값
    group[i] <- which.min(mat[i,]) #min ind
  }

  res <- data.frame(cbind(distance, group))
  return(res)
}

# Ranking
ranking <- function(data){
  dist <- data$distance
  rank <- rank(dist)
  alpha <- (length(dist) - rank + 1) / length(dist)
  # dist 가까울수록 alpha값 큼
  data <- data.frame(cbind(data, rank, alpha))
}

# Credibility
cred <- function(data, center, parent){
  mat1 <- distance(data, center)
  mat2 <- distance(parent, center)

  mat <- rbind(mat1, mat2)
  rank <- c()
  for(i in 1:nrow(center)){
    rank <- cbind(rank, rank(mat[,i]))
  }

  rank <- rank[1:nrow(data),]
  alpha <- (nrow(data)+nrow(parent)-rank+1)/(nrow(data)+nrow(parent))
  # 각 center까지의 거리(열)에 대해서 rank 계산 -> 거리 가까울수록 rank 작고 alpha 큼
  return(alpha)
}


# Silhouette
silhouette <- function(data, center){
  mat <- distance(data, center)

  group <- c()
  group_near <- c()
  for(i in 1:nrow(data)){
    dist_tmp <- mat[i,]
    group[i] <- which.min(dist_tmp)
    min_second <- sort(dist_tmp)[2]
    group_near[i] <- which(dist_tmp==min_second)
  }

  dat <- data.frame(cbind(data, group, group_near))

  res <- c()
  for(i in 1:nrow(data)){
    inner <- dat[which(dat$group==dat$group[i]),]
    outer <- dat[which(dat$group==dat$group_near[i]),]

    dist_in <- c()
    for(j in 1:nrow(inner)){
      dist_in[j] <- ((data[i,1] - inner[j,1])^2 + (data[i,2] - inner[j,2])^2)
    }
    dist_out <- c()
    for(j in 1:nrow(outer)){
      dist_out[j] <- ((data[i,1] - outer[j,1])^2 + (data[i,1] - outer[j,1])^2)
    }

    a <- sum(dist_in)/(length(dist_in)-1)
    b <- mean(dist_out)
    res[i] <- (b - a)/max(a, b)
  }

  res <- mean(res)

  return(res)
}

bic <- function(df, parentnum, logllh){
  n <- dim(df)[1]
  k <- parentnum+3
  lnL <- logllh
  return(log(n)*k-2*lnL)
}

MCMCestThomas <- function(X, xlim, ylim, NStep = 10000, DiscardStep = 1000, Jump=10,
                          AreaW=1, sd_alpha=0.1, sd_omega=0.015, ...)
{
  # Update the alpha, omega, CC
  W = owin(xrange=c(0,1), yrange=c(0,1))
  # pBX <- pBetaX(X, NStep, Jump, xlim, ylim, W, AreaW, sd_alpha, sd_omega)

  ####################################################################################
  # pBetaX
  ####################################################################################
  # Initial value
  lambda <- dim(X)[1]/AreaW # p/|S|: overall intensity of items over the interaction map domain
  kappa <- runif(1, 3, 7) # initial value of kappa
  alpha <- lambda/kappa # alpha = p/(|S|*kappa)
  omega <- sqrt(AreaW)/kappa

  # initial value of cluster center CC
  CC <- matrix()
  while(dim(CC)[1]<=1){
    # rpoispp(lambda, win): generate poisson point pattern
    # lambda: Intensity of the Poisson process, win: sindow in which to simulate tha pattern
    CC <- rpoispp(kappa, win = W)
    CC <- cbind(CC$x, CC$y) # location
  }

  # integrate the Gaussian kerneal with center at ci with variance omega
  integral <- Kumulavsech(CC, omega, xlim, ylim)

  # calculate the log likelihood
  logP <- logpXCbeta(X, CC, alpha, omega, AreaW, integral)
  #
  pbx <- c(kappa, alpha, omega, logP, integral)
  accept <- 0

  oldk = kappa


  # Update function for the alpha and omega
  StepBeta <- function(kappa, alpha, omega, sd_alpha, sd_omega, X, CC, logP, integral, AreaW){
    # sd_alpha <- 0.1; sd_omega <- 0.015
    Newalpha <- rnorm(1,alpha, sd_alpha)
    Newomega <- rnorm(1,omega, sd_omega)
    Newkappa <- dim(X)[1]/AreaW/Newalpha

    if( Newalpha < dim(X)[1]/AreaW/10 | Newalpha > dim(X)[1]/AreaW/2 | Newomega < 0){
      logratio=-Inf }

    else{
      #while((Newalpha <- rnorm(1, alpha, sd_alpha)) < dim(X)[1]/AreaW/9){}
      #while((Newomega <- rnorm(1, omega, sd_omega)) < 0){}
      #Newkappa <- dim(X)[1]/AreaW/Newalpha;

      int2 <- Kumulavsech(CC, Newomega, xlim, ylim)
      logP1 <- logpXCbeta(X, CC, Newalpha, Newomega, AreaW, int2)

      logratio <-  (logP1 - logP + kappa*AreaW - Newkappa*AreaW +
                      dim(CC)[1]*log((alpha/Newalpha)) +
                      log((1 - pnorm(0, alpha, sd_alpha))/(1-pnorm(0, Newalpha, sd_alpha)
                                                           *(1 - pnorm(0, omega, sd_omega))/(1-pnorm(0, Newomega, sd_omega)))))
    }

    u <- runif(1)
    if(log(u) < logratio){
      Vystup <- c(Newkappa, Newalpha, Newomega, logP1, int2)
    }
    else{
      Vystup <- c(kappa, alpha, omega, logP, integral)
    }

    Vystup
  }

  StepMovePoint <- function(kappa, alpha, omega, X, CC, logP, integral, AreaW){

    # Discard + Propose a point
    if(runif(1) < 1/3){
      Discard <- ceiling(runif(1, min = 0, max = dim(CC)[1])) # select the point to discard
      int2 <- integral - Kumulavsech(t(as.matrix(CC[Discard,])), omega, xlim, ylim)
      CC1 <- CC[-Discard,]
      NewCenter <- t(NewPoint(xlim, ylim))
      CC1 <- rbind(CC1, NewCenter)
      int2 <- int2 + Kumulavsech(as.matrix(NewCenter), omega, xlim, ylim)
      logP2 <- logpXCbeta(X, CC1, alpha, omega, AreaW, int2)
      if(log(runif(1)) < (logP2 - logP)){
        Vystup <- list(CC1, logP2, int2)
      }
      else{
        Vystup <- list(CC, logP, integral)
      }
    }
    else{
      # Only propose the parent points
      if(runif(1) < 1/2 || dim(CC)[1] < 3){
        NewCenter <- t(NewPoint(xlim, ylim))
        CC1 <- rbind(CC, NewCenter)
        int2 <- integral + Kumulavsech(as.matrix(NewCenter), omega, xlim, ylim)
        logP2 <- logpXCbeta(X, CC1, alpha, omega, AreaW, int2)
        if(log(runif(1)) < (logP2 - logP + log(kappa*(AreaW)) - log(dim(CC1)[1]))){
          Vystup <- list(CC1, logP2, int2)
        }
        else{
          Vystup <- list(CC, logP, integral)
        }
      }
      else{
        # Only discard the parent points
        Discard <- ceiling(runif(1, min = 0, max = dim(CC)[1]))
        int2 <- integral - Kumulavsech(t(as.matrix(CC[Discard,])), omega, xlim, ylim)
        CC1 <- CC[-Discard,]
        logP2 <- logpXCbeta(X, CC1, alpha, omega, AreaW, int2)
        if(log(runif(1)) < (logP2 - logP - log(kappa*(AreaW)) + log(dim(CC1)[1]))){
          Vystup <- list(CC1, logP2, int2)
        }
        else{
          Vystup <- list(CC, logP, integral)
        }
      }
    }

    Vystup
  }

  for(step in 1:NStep){
    oldk <- kappa

    # Update function for the alpha and omega
    S <- StepBeta(kappa, alpha, omega, sd_alpha, sd_omega, X, CC, logP, integral, AreaW)

    kappa <- S[1]
    alpha <- S[2]
    omega <- S[3]
    logP <- S[4]
    integral <- S[5]

    if(oldk!=kappa){accept=accept+1}

    # Fitting Thomas process with Metropolis-Hastings algorithm and birth-death MCMC
    T <- StepMovePoint(kappa, alpha, omega, X, CC, logP, integral, AreaW)
    CC <- T[[1]]
    logP <- T[[2]]
    integral <- T[[3]]

    if((step %% Jump)==0){
      pbx <- rbind(pbx, S)
    }
  }
  pBX <- list(pbx=pbx, CC=CC, accept=accept/NStep)
  ####################################################################################


  pbx <- pBX$pbx
  Postkappa <- pbx[(round(DiscardStep/Jump)+1):(dim(pbx)[1]),1]
  Postalpha <- pbx[(round(DiscardStep/Jump)+1):(dim(pbx)[1]),2]
  Postomega <- pbx[(round(DiscardStep/Jump)+1):(dim(pbx)[1]),3]
  PostlogP <- pbx[(round(DiscardStep/Jump)+1):(dim(pbx)[1]),5]

  accept <- pBX$accept

  Kappahat <- median(Postkappa)
  Alphahat <- median(Postalpha)
  Omegahat <- median(Postomega)
  KappaCI <- quantile(Postkappa, probs = c(0.025,0.975))
  AlphaCI <- quantile(Postalpha, probs = c(0.025,0.975))
  OmegaCI <- quantile(Postomega, probs = c(0.025,0.975))

  CC <- pBX$CC

  res <- list(Postkappa = Postkappa,
              Postalpha=Postalpha,
              Postomega=Postomega,
              Kappahat=Kappahat,
              Alphahat=Alphahat,
              Omegahat=Omegahat,
              KappaCI=KappaCI,
              AlphaCI=AlphaCI,
              OmegaCI=OmegaCI,
              Accept=accept,
              pBX=pbx[(round(DiscardStep/Jump)+1):(dim(pbx)[1]),],
              CC=CC
  )

  res
}

