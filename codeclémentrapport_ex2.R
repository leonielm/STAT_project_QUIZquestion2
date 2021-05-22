#set functions
{#CTRL ENTER HERE ----------------------------
  
  cor_nonlinear <- function(muX, muY, sX, sY, angle){
    # muX and muY have no effect
    sQ <- sqrt(cos(angle)**2 * (2*sX**4 + 4*sX**2*muX**2 +  sY**2) + sin(angle)**2 * sX**2 + 2 * cos(angle)*sin(angle)*(2*sX**2*muX))
    sR <- sqrt(sin(angle)**2 * (2*sX**4 + 4*sX**2*muX**2 + sY**2) + cos(angle)**2 * sX**2 - 2 * cos(angle)*sin(angle)*(2*sX**2*muX))
    cov <- sin(angle)*cos(angle) * (sX**2 - (2*sX**4 + 4*sX**2*muX**2 + sY**2)) +
      (cos(angle)**2 - sin(angle)**2) * (2*sX**2*muX)
    cor <- round(cov/(sQ*sR), 4)
    return(cor)
  }
  
  
  CI.sims <- function(n, setting, method, type, B, M, 
                      muX, muY, sX, sY, rho=NULL, out=NULL, dev=NULL, angle=NULL, level=0.8){
    if(!setting %in% c("normal", "discrete", "outliers", "nonlinear")) stop("Wrong argument `setting`")
    if(!method %in% c("pearson", "spearman")) stop("Wrong argument `method`")
    if(!type %in% c("boot", "npboot")) stop("Wrong argument `type`")
    
    if(setting=="normal"){
      if(is.null(rho)) stop("Set a value for rho.")
      gen_data <- function() gen_normal(n, muX, muY, sX, sY, rho)
    } else if(setting=="outliers"){
      if(is.null(rho)) stop("Set a value for rho.")
      if(is.null(out)) stop("Set a value for out.")
      if(is.null(dev)) stop("Set a value for dev.")
      if(out>1) stop("The value for out must be between 0 and 1.")
      gen_data <- function() gen_outliers(n, muX, muY, sX, sY, rho, out=out, dev=dev)
    } else if(setting=="discrete"){
      gen_data <- function() gen_discrete(n, muX, muY, sX, sY, rho)
    } else if(setting=="nonlinear"){
      if(is.null(angle)) stop("Set a value for angle.")
      gen_data <- function() gen_nonlinear(n, muX, muY, sX, sY, angle)
    }
    probs <- c((1-level)/2, level + (1-level )/2)
    
    if(type=="boot"){
      CI <- t(sapply(1:M, function(ni){
        df <- gen_data()
        muX <- mean(df$x)
        muY <- mean(df$y)
        sX <- sd(df$x)
        sY <- sd(df$y)
        rho <- cor(df$x,df$y)
        quantile(sapply(1:B, function(na){
          cor(gen_normal(n, muX, muY, sX, sY, rho), method=method)[2]
        }), probs=probs)}
      ))
    }else if (type=="npboot"){
      CI <- t(sapply(1:M, function(ni){
        df <- gen_data()
        quantile(sapply(1:B, function(na){
          s <- sample(n, replace=T)
          cor(df[s,], method=method)[2]
        }), probs=probs)}
      ))
    } else {
      stop("type must be one of ('boot', 'npboot'")
    }
    CI
  }
  
  
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
}

#ctrl enter chaque ligne (ligne 115 à 264) 
#Mettez les valeurs qui vous sont demandées dans la console (en bas)

findRho1_2 = function(muX, muY, sX, sY, angle1, angle2){
  print("Looking for Rho1 and Rho2 :")
  
  rho1 = cor_nonlinear(muX, muY, sX, sY, angle1)
  rho2 = cor_nonlinear(muX, muY, sX, sY, angle2)
  
  print(c("Rho1 = ", rho1))
  print(c("Rho1 = ", rho2))
  
  return(rho2)
} 

findBeta = function(seed, n, muX, muY, sX, sY, angle, rho2, M){
  
  res = 0
  set.seed(seed)
  for (i in 1:M){
    dfi = gen_nonlinear(n, muX, muY, sX, sY, angle)
    stock = cor.test(dfi[,1],dfi[,2], method = "pearson")
    res[i] = stock[[4]]
  }
  
  return(round(abs(rho2 - mean(res)), 5))
  
}


doPage2 = function(){
  
  
  print("enter muX")
  muX = scan(what = double(), nmax = 1)
  print("Enter muY")
  muY = scan(what = double(), nmax = 1)
  print("Enter sX")
  sX = scan(what = double(), nmax = 1)
  print("Enter sY")
  sY = scan(what = double(), nmax = 1)
  print("Enter the 1rst n of question 1")
  n1 = scan(what = double(), nmax = 1)
  print("Enter 1rst angle")
  angle1 = scan(what = double(), nmax = 1)
  print("Enter 2nd angle")
  angle2 = scan(what = double(), nmax = 1)
  
  rho2 = findRho1_2(muX, muY, sX, sY, angle1, angle2)
  
  
  print("Enter 2nd n (for beta)")
  n2 = scan(what = double(), nmax = 1)
  print("Enter SIMULATIONSEED")
  SIMULATIONSEED = scan(what = double(), nmax = 1)
  print("Enter M")
  M = scan(what = double(), nmax = 1)
  
  
  
  Beta = findBeta(SIMULATIONSEED, n2, muX, muY, sX, sY, angle2, rho2, M)
  print(c("Beta = ", Beta ))
  
  
  print("Now second part : ")
  print("Enter muX")
  muX = scan(what = double(), nmax = 1)
  print("Enter muY")
  muY = scan(what = double(), nmax = 1)
  print("Enter sX")
  sX = scan(what = double(), nmax = 1)
  print("Enter sY")
  sY = scan(what = double(), nmax = 1)
  print("Enter rho")
  rho = scan(what = double(), nmax = 1)
  print("Enter n for 1rst setting")
  n1 = scan(what = double(), nmax = 1)
  print("Enter n for 2nd setting")
  n2 = scan(what = double(), nmax = 1)
  print("Enter n for 3rd setting")
  n3 = scan(what = double(), nmax = 1)
  print("Enter n for 4th setting")
  n4 = scan(what = double(), nmax = 1)
  print("Enter angle for 4th setting")
  angle4 = scan(what = double(), nmax = 1)
  print("Enter out for 2nd setting")
  out = scan(what = double(), nmax = 1)
  print("Enter dev for 2nd setting")
  dev = scan(what = double(), nmax = 1)
  print("Enter SAMPLESEED")
  SAMPLESEED = scan(what = double(), nmax = 1)
  print("Enter B")
  B = scan(what = double(), nmax = 1)
  print("Enter M")
  M = scan(what = double(), nmax = 1)
  population_correlation <- rho
  
  
  print("enter ''1'' for boot or ''2'' for npboot in C1")
  cMtd = scan(what = double(), nmax = 1)
  if (cMtd == 1) Mtd = "boot" else Mtd =  "npboot"
  
  set.seed(SAMPLESEED)
  CI.sim1 <- CI.sims(n1, "normal", "spearman", Mtd, B, M, muX, muY, sX, sY, rho)
  
  
  coverage1 <- mean((CI.sim1[,1]< population_correlation) * (CI.sim1[,2]> population_correlation))
  print(c("C1 = ", coverage1))
  
  
  
  print("enter ''1'' for boot or ''2'' for npboot in C2")
  cMtd = scan(what = double(), nmax = 1)
  if (cMtd == 1) Mtd = "boot" else Mtd =  "npboot"
  
  set.seed(SAMPLESEED)
  CI.sim2 <- CI.sims(n2, "outliers", "spearman", Mtd, B, M, muX, muY, sX, sY, rho, out, dev)
  
  
  coverage2 <- mean((CI.sim2[,1]< population_correlation) * (CI.sim2[,2]> population_correlation))
  print(c("C2 = ", coverage2))
  
  
  
  print("enter ''1'' for boot or ''2'' for npboot in C3")
  cMtd = scan(what = double(), nmax = 1)
  if (cMtd == 1) Mtd = "boot" else Mtd =  "npboot"
  
  set.seed(SAMPLESEED)
  CI.sim3 <- CI.sims(n3, "discrete", "spearman", Mtd, B, M, muX, muY, sX, sY, rho)
  
  
  coverage3 <- mean((CI.sim3[,1]< population_correlation) * (CI.sim3[,2]> population_correlation))
  print(c("C3 = ", coverage3))
  
  
  
  
  print("enter ''1'' for boot or ''2'' for npboot in C4")
  cMtd = scan(what = double(), nmax = 1)
  if (cMtd == 1) Mtd = "boot" else Mtd =  "npboot"
  
  set.seed(SAMPLESEED)
  CI.sim4 <- CI.sims(n4, "nonlinear", "spearman", Mtd, B, M, muX, muY, sX, sY, angle = angle4)
  population_correlation = cor_nonlinear(muX, muY, sX, sY, angle4)
  
  
  coverage4 <- mean((CI.sim4[,1]< population_correlation) * (CI.sim4[,2]> population_correlation))
  print(c("C4 = ", coverage4))
  
  
}
doPage2()


