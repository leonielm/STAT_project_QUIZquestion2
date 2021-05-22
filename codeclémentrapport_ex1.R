#Exercice 1

#set function : 
{#CTRL ENTER HERE ----------------------------
  # Generate Bivariate Gaussian Observations
  gen_normal <- function(n, muX=0, muY=0, sX=1, sY=1, rho=0){
    x <- rnorm(n, muX, sX)
    z <- rnorm(n, muY, sY * sqrt(1-rho^2))
    y <- sY/sX * rho * (x-muX) + z
    df <- data.frame(x=x, y=y)
    return(df)
  }
  
  # Generate Discretized Bivariate Gaussian Observations
  gen_discrete <- function(n, muX=0, muY=0, sX=1, sY=1, rho=0){
    x <- rnorm(n, muX, sX)
    z <- rnorm(n, muY, sY * sqrt(1-rho^2))
    y <- sY/sX * rho * (x-muX) + z
    df <- data.frame(x=trunc(x*5)/5, y=trunc(y*5)/5)
    return(df)
  }
  
  
  # Generate data with outliers
  gen_outliers <- function(n, muX=0, muY=0, sX=1, sY=1, rho=0, out=0.05, dev=1){
    df <- gen_normal(n, muX, muY, sX, sY, rho)
    
    if(exists(".Random.seed")){
      INITIALSEED <- .Random.seed
      lockBinding("INITIALSEED", environment())
      on.exit(.Random.seed <<- INITIALSEED)
    }
    outliers=sample(n,ceiling(out*n))
    values_from <- sample(n, ceiling(out*n))
    
    df$y[outliers] <- df$y[values_from] + (df$y[values_from]-mean(df$y)) * dev
    return(df)
  }
  
  # Generate Bivariate Data with Nonlinear Relationship
  gen_nonlinear <- function(n, muX, muY, sX, sY, angle){
    rotmat <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)),2,2)
    x <- rnorm(n, muX, sX)
    y <- x**2 + rnorm(n, muY, sY)
    df <- cbind(x,y)
    df <- df %*% rotmat
    df <- as.data.frame(df)
    names(df) <- c("x", "y")
    return(df)
  }
  
  boot.CI <- function(x,y, B, method){
    n <- length(x)
    muX.est <- mean(x)
    muY.est <- mean(y)
    sX.est  <- sd(x)
    sY.est  <- sd(y)
    rho.est <- cor(x, y)
    
    rXY.boot <- rep(0, B)
    set.seed(SIMULATIONSEED)
    for(i in 1:B){
      data.boot <-   gen_normal(n=n, 
                                muX=muX.est, 
                                muY=muY.est, 
                                sX=sX.est, 
                                sY=sY.est, 
                                rho=rho.est)
      rXY.boot[i] <- cor(data.boot$x, data.boot$y, method=method)
    }
    
    # Get the .95 bootstrap percentile confidence interval
    boot.CI <- quantile(rXY.boot, probs=c(0.025, 0.975))
    boot.CI
  }
  np.boot.CI <- function(x,y, B, method){
    n <- length(x)
    
    rXY.boot <- rep(0, B)
    set.seed(SIMULATIONSEED)
    for(i in 1:B){
      s <- sample(n, replace = T)
      rXY.boot[i] <- cor(x[s], y[s], method=method)
    }
    
    # Get the .95 bootstrap percentile confidence interval
    np.boot.CI <- quantile(rXY.boot, probs=c(0.025, 0.975))
    np.boot.CI
  }
  
}

#set valeurs :
{#Insérez vos valeurs de l'énoncé -------------------------------
  SAMPLESEED = 41696 #ici
  SIMULATIONSEED = 68239 #ici
  muX = -2.4 #ici
  muY = -1.3 #ici
  sX = 0.6 #ici
  sY = 0.7 #ici
  rho = 0.44 #ici
  angle = 0.45 #ici
  dev = 3 #ici
  out = 0.07 #ici
}

#set other function :
{#CTRL ENTER HERE ----------------------------
  askN = function(){
    print("Entrer le n :")
    n = scan(what = double(), nmax = 1)
    return(n)
  }
  
  askMTD = function() {
    print("entrer le chiffre pour votre methode (1 --> pearson, 2 --> spearman)")
    MTDc = scan(what = double(), nmax = 1)
    if (MTDc == 1) MTD = "pearson" else MTD = "spearman"
    return(MTD)
  }
  
  askB = function(){
    print("entrer B : ")
    B = scan(what = double(), nmax = 1)
    return(B)
  }
  
  askBoot = function() {
    print("Parametric Bootstrap --> 1 / nonParametric Bootstrap --> 2")
    bootsN = scan(what = double(), nmax = 1)
    return(bootsN)
  }
}




Question1 = function(seed = SAMPLESEED, muX1 = muX, muY1 = muY, sX1 = sX, sY1 = sY, rho1 = rho){
  print("question 1 :")
  
  n = askN()
  set.seed(seed)
  df1 = gen_normal(n, muX1, muY1, sX1, sY1, rho1)
  
  print("val generated")
  
  MTD = askMTD()
  rhoF1 = cor.test(df1[,1],df1[,2], method = MTD)[[4]]
  
  
  if(MTD == "spearman") {
    B = askB()
    isBoot = askBoot()
    
    if (isBoot == 1) L1 = boot.CI(df1[,1], df1[,2], B, MTD)[[1]] else L1 = np.boot.CI(df1[,1], df1[,2], B, MTD)[[1]]
  } else L1 = cor.test(df1[,1],df1[,2], method = MTD)[[9]]
  
  
  print(c("Rho1 = ", round(rhoF1, 5)))
  print(c("L1 = ", round(L1[1], 5)))
  
  
  
}
Question1()

Question2 = function(seed = SAMPLESEED, muX2 = muX, muY2 = muY, sX2 = sX, sY2 = sY, rho2 = rho, out2 = out,dev2 = dev ){
  print("question 2 :")
  
  n = askN()
  set.seed(seed)
  df1 = gen_outliers(n, muX2, muY2, sX2, sY2, rho2, out2, dev2)
  
  print("val generated")
  
  MTD = askMTD()
  rhoF1 = cor.test(df1[,1],df1[,2], method = MTD)[[4]]
  
  
  if(MTD == "spearman") {
    B = askB()
    isBoot = askBoot()
    
    if (isBoot == 1) L1 = boot.CI(df1[,1], df1[,2], B, MTD)[[1]] else L1 = np.boot.CI(df1[,1], df1[,2], B, MTD)[[1]]
    
  } else L1 = cor.test(df1[,1],df1[,2], method = MTD)[[9]]
  
  
  
  print(c("Rho2 = ", round(rhoF1, 5)))
  print(c("L2 = ", round(L1[1], 5)))
  
  
  
}
Question2()

Question3 = function(seed = SAMPLESEED, muX3 = muX, muY3 = muY, sX3 = sX, sY3 = sY, rho3 = rho){
  print("question 3 :")
  
  n = askN()
  set.seed(seed)
  df1 = gen_discrete(n, muX3, muY3, sX3, sY3, rho3)
  
  print("val generated")
  
  MTD = askMTD()
  rhoF1 = cor.test(df1[,1],df1[,2], method = MTD)[[4]]
  
  
  if(MTD == "spearman") {
    B = askB()
    isBoot = askBoot()
    
    if (isBoot == 1) L1 = boot.CI(df1[,1], df1[,2], B, MTD)[[1]] else L1 = np.boot.CI(df1[,1], df1[,2], B, MTD)[[1]]
  } else L1 = cor.test(df1[,1],df1[,2], method = MTD)[[9]]
  
  
  print(c("Rho3 = ", round(rhoF1, 5)))
  print(c("L3 = ", round(L1[1], 5)))
  
  
  
}
Question3()

Question4 = function(seed = SAMPLESEED, muX4 = muX, muY4 = muY, sX4 = sX, sY4 = sY, angle4 = angle){
  print("question 4 :")
  
  n = askN()
  set.seed(seed)
  df1 = gen_nonlinear(n, muX4, muY4, sX4, sY4, angle4)
  
  print("val generated")
  
  MTD = askMTD()
  rhoF1 = cor.test(df1[,1],df1[,2], method = MTD)[[4]]
  
  
  if(MTD == "spearman") {
    B = askB()
    isBoot = askBoot()
    
    if (isBoot == 1) L1 = boot.CI(df1[,1], df1[,2], B, MTD)[[1]] else L1 = np.boot.CI(df1[,1], df1[,2], B, MTD)[[1]]
  } else {
    B = askB()
    isBoot = askBoot()
    
    if (isBoot == 1) L1 = boot.CI(df1[,1], df1[,2], B, MTD)[[1]] else L1 = np.boot.CI(df1[,1], df1[,2], B, MTD)[[1]]
  }
  
  
  print(c("Rho4 = ", round(rhoF1, 5)))
  print(c("L4 = ", round(L1[1], 5)))
  
  
  
}
Question4()
