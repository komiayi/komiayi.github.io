rm(list=ls())
library(rpart)  
library(rpart.plot)
library(tibble)
library(tidyverse)
library(ipred)
library(randomForest)

###################################################
#####  Figure 1

plot.sigmabar <- function(sigma,B1){
  rho <- seq(0, 1, length = 100)
  plot(1, xlim=c(0.0001,1), ylim = c(0.5,4.1), type="n", xlab = expression(rho), 
       ylab = expression(bar(sigma)**2))
  abline(h=sigma**2, col='red')
  for(i in B1){
    b <- 1+i*rho
    sigmabar <-rho*sigma**2+sigma**2*(1-rho)/b
    lines(rho, sigmabar,type ='l', lwd = 3, col = i, lty = i)
  }
  legend("bottomright", legend=c("B=25", "B=50", "B=75", "B=100"),
         col=c(25, 50, 75, 100), cex=0.6, lty = c(25, 50, 75, 100), lwd = 3)
}

sigma <- 2
B1 <- c(25,50,75,100)
plot.sigmabar(sigma, B1)
######################################################




# erreur quadrarique moyenne
mse <- function(actual,pred) {mean((actual-pred)^2)}

# création de grille 
grille <- expand.grid(
  minsplit = c(3, 9,24),
  maxdepth = c(1, 10, 15)
)

# fonction pour déterminer le cp otimal
optimal.cp <- function(x) {
  min    <- which.min(x$cptable[, "xerror"])
  cp <- x$cptable[min, "CP"] 
  return(cp)
}

# fonction pour détrminer l'erruer minimale

optimal.error <- function(x) {
  min    <- which.min(x$cptable[, "xerror"])
  xerror <- x$cptable[min, "xerror"] 
  return(xerror)
}

# fonction qui calcule les paramètres optimaux de toutes combinaisons de la grille
#+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20

parameters<- function(grille, train, test){
  results <- data.frame(NULL)
  for(i in seq_len(nrow(grille))) {
    controls <- list(cp = 0, minsplit = grille$minsplit[i], maxdepth = grille$maxdepth[i])
    tree.fit <- rpart(y ~ X1+X2+X3+X4+X5,
                      data = train, method="anova", control = controls) 
    
    predictions <- test %>%
      mutate( minsplit = factor(paste("Observation=", grille$minsplit[i]), ordered = TRUE),
              maxdepth = factor(paste("Profondeur =", grille$maxdepth[i]), ordered = TRUE))
    
    predictions$pred <- predict(tree.fit, newdata= test, type = "vector")
    grille$cp[i] <- optimal.cp(tree.fit)
    grille$error[i] <- optimal.error(tree.fit)
    grille$MSE[i] <- mse(test$y,predictions$pred)
    results <- rbind(results, predictions)
  }
  return(list(results = results, grille = grille))
}

# validation croisé de m

m.cv <- function(data, m, B, fold){
  
  dat <- data
  dat$group <- sample.int(fold, nrow(dat), replace = TRUE)
  mse_m <- NULL
  
  for (u in 1:fold){
    mse_m.i <- c()
    for(i in 1:length(m)){
      rf.m<- randomForest(y ~ X1+X2+X3+X4+X5,
                          data = dat[dat$group==u,],      
                          mtry=i, ntree = B)
      pred <- predict(rf.m, newdata = dat[dat$group!=u,])
      mse_m.i <- c(mse_m.i, mse(dat[dat$group!=u,]$y,pred))
    }
    mse_m <- rbind(mse_m,mse_m.i)
  }
  dat <- NULL
  ind.m <- as.vector(apply(mse_m,2,mean))
  return(m[which.min(ind.m)])
}

# fonction aui calcule MSE moyen selon la méthode choisie 
mse.comparaison <- function(Data, method="arbre", m, B=100, fold=3, J){
  MSE <- c()
  ### Arbre de régression
  if(method=="arbre"){
    for(j in 1:J){
      #data <- Data[,j]
      data <- Data
      train <- data$train
      test <- data$test
      
      Finals <- parameters(grille, train, test)
      grilles <- Finals$grille %>%
        arrange(cp)
      controls_1 <- grilles[which.min(grilles$error),]
      
      controls_2 <- list(minsplit = controls_1[[1]], 
                         maxdepth = controls_1[[2]],
                         cp = controls_1[[3]])
      tree.optim <- rpart(y ~ X1+X2+X3+X4+X5,
                          data = train, method ="anova", control = controls_2) 
      MSE[j] <- mse(test$y, predict(tree.optim, newdata = test))
    }
    return(list(MSE = MSE))
  }else if(method =="rf"){ # Random Forest
    m.opt <- c()
    for(j in 1:J){
      #data <- Data[,j]
      data <- Data
      train <- data$train
      test <- data$test
      
      M <- m.cv(rbind(train,test), m, B, fold)
      rf.optim <- randomForest(y ~ X1+X2+X3+X4+X5,
                               data = train,      
                               mtry = M, ntree = B)
      m.opt <- c(m.opt, M)
      MSE[j] <- mse(test$y, predict(rf.optim, newdata = test))
    }
    return(list(MSE = MSE, m.opt= m.opt))#, m.opt= m.opt
    
  }else if(method =="bagging"){ # Bagging
    for(j in 1:J){
      #data <- Data[,j]
      data <- Data
      train <- data$train
      test <- data$test
      
      bg.optim <- bagging(y ~ X1+X2+X3+X4+X5,
                          data = train, nbagg = B)
      MSE[j] <- mse(test$y, predict(bg.optim, newdata = test))
    }
    return(list(MSE = MSE))
  }
}

##### stabilité arbre de régression 
# paramètre "stab = oui" , remplace z par z'

stability.arbre <- function(Data, Z, p){
  #X1+X2+X3+X4+X5

  data <- Data           # décommenter pour calcul avec donnée réel
  train.1 <- data$train 
  test.1 <- data$test
  Finals.1 <- parameters(grille, train.1, test.1[1:nrow(train.1),])
  grilles.1 <- Finals.1$grille %>%
    arrange(cp)
  controls_1 <- grilles.1[which.min(grilles.1$error),]
  
  controls_2 <- list(minsplit = controls_1[[1]], 
                     maxdepth = controls_1[[2]],
                     cp = controls_1[[3]])
  ### arbre optimal 
  
  tree.optim.1 <- rpart(y ~ X1+X2+X3+X4+X5,
                        data = train.1[1:nrow(train.1),], method ="anova",
                        control = controls_2) 
  Predict.1 <- predict(tree.optim.1, newdata = test.1[1:nrow(train.1),])
  l.1 <- (test.1$y[1:nrow(train.1)]-Predict.1)^2
  l.2 <- c()
  for(i in 1:nrow(data$train)){
    
    train.2 <- data$train
    test.2 <- data$test
    Finals.2 <- NULL; grilles.2 <- NULL; tree.optim.2 <- NULL; Predict.2 <- NULL
    controls_1.2 <- NULL; controls_2.2 <- NULL
    
    for(ind in 1:(p+1)){
      train.2[i,ind] <- Z[ind]
    }
    
    Finals.2 <- parameters(grille, train.2, test.1[1:nrow(train.1),])
    grilles.2 <- Finals.2$grille %>%
      arrange(cp)
    controls_1.2 <- grilles.2[which.min(grilles.2$error),]
    
    controls_2.2 <- list(minsplit = controls_1.2[[1]], 
                       maxdepth = controls_1.2[[2]],
                       cp = controls_1[[3]])
    tree.optim.2 <- rpart(y ~ X1+X2+X3+X4+X5,
                          data = train.2, method ="anova",
                          control = controls_2.2) 
    Predict.2 <- predict(tree.optim.2, newdata = test.1[1:nrow(train.1),])
    l.2 <- c(l.2, ((test.1$y[1:nrow(train.1)]-Predict.2)^2))
  }
  return(sum(l.2-l.1))
}

##### stabilité bagging 

stability.bagging <- function(Data, Z, p, B){
  #X1+X2+X3+X4+X5

  data <- Data
  train.1 <- data$train 
  test.1 <- data$test
  bg.optim.1 <- bagging(y ~ X1+X2+X3+X4+X5,
                      data = train.1, nbagg = B)
  Predict.1 <- predict(bg.optim.1, newdata = test.1[1:nrow(train.1),])
  l.1 <- (test.1$y[1:nrow(train.1)]-Predict.1)^2
  l.2 <- c()
  for(i in 1:nrow(data$train)){
    bg.optim.2 <- NULL; Predict.2 <- NULL
    train.2 <- data$train
    for(ind in 1:(p+1)){
      train.2[i,ind] <- Z[ind]   #remplace zi par z'
    }
    bg.optim.2 <- bagging(y ~ X1+X2+X3+X4+X5,
                          data = train.2, nbagg = B)
    Predict.2 <- predict(bg.optim.2, newdata = test.1[1:nrow(train.1),])
    l.2 <- c(l.2, ((test.1$y[1:nrow(train.1)]-Predict.2)^2))
  }
  return( sum(l.2-l.1))
}


#####   stabilité de random Forest 

stability.rf <- function(Data, Z, p, m, B, fold=3){

  data <- Data
  train.1 <- data$train 
  test.1 <- data$test
  M <- m.cv(rbind(train.1,test.1[1:nrow(train.1),]), m,B, fold)
  rf.optim.1 <- randomForest(y ~ X1+X2+X3+X4+X5,
                           data = train.1,      
                            mtry = M, ntree = B)
  Predict.1 <- predict(rf.optim.1, newdata = test.1[1:nrow(train.1),])
  l.1 <- (test.1$y[1:nrow(train.1)]-Predict.1)^2
  l.2 <- c()
  for(i in 1:nrow(data$train)){
    rf.optim.2 <- NULL; Predict.2 <- NULL
    
    train.2 <- data$train
    test.2 <- data$test
    for(ind in 1:(p+1)){
      train.2[i,ind] <- Z[ind]   #remplace zi par z'
    }
    M <- m.cv(rbind(train.2,test.1[1:nrow(train.1),]), m, B, fold)
    rf.optim.2 <- randomForest(y ~ X1+X2+X3+X4+X5,
                               data = train.2,      
                               mtry = M, ntree = B)
    Predict.2 <- predict(rf.optim.2, newdata = test.1[1:nrow(train.1),])
    l.2 <- c(l.2, ((test.1$y[1:nrow(train.1)]-Predict.2)^2))
  }
  return( sum(l.2-l.1))
}

####################################################
## Génération de la base ####

# data Friedman #1 

n <- 1200
p <- 5
m <- 1:p
B <- 100


dataset <- function(n,p){
  
  data1 <- tibble::tibble(
    x <- data.frame(replicate(p,runif(n))),
    y <-  10*sin(x[,1]*x[,2])+20*(x[,3]-0.5)^2+10*x[,4]+5*x[,5] + rnorm(n)
  )
  
  data <- data1 %>%
    rename( "y"= p+1 ) %>%
    mutate(id = row_number())
  
  train <- data %>%
    sample_frac(.1666667)  # entrainement 
  
  test <- anti_join(data, train, by = 'id')  # test
  return(list(train = train, test = test))
}

# réplication de 1000 jeu de données 

set.seed(123)
Data <- replicate(1000,dataset(n,p), simplify = "array")

########################################################
#####  Calcul et comparaison des MSE (Table 2) ########
J <- dim(Data)[2]

Arbre <- mse.comparaison(Data, method = "arbre", J=J)

Bagging <- mse.comparaison(Data, method = "bagging", J=J)

Random.Forest <- mse.comparaison(Data, method = "rf", m, B, fold = 3,J=J)

A1 <- c(mean(Arbre$MSE),var(Arbre$MSE))
b1 <- c(mean(Bagging$MSE),var(Bagging$MSE))
R1 <- c(mean(Random.Forest$MSE),var(Random.Forest$MSE))

tableau <- rbind(A1,b1,R1)
rownames(tableau) <- c("Arbre", "Bagging", "Random Forest")
colnames(tableau) <- c("MSE", "Var(MSE)")

# m optimal obtenu par CV (Table 1) 
table.m <- table(Random.Forest$m.opt)


##########################  Stabilité ########################

## calcule de z'
set.seed(524)
X <- replicate(p,runif(1))
Y <-  10*sin(X[1]*X[2])+20*(X[3]-0.5)^2+10*X[4]+5*X[5] + rnorm(1)
Z <- c(X,Y)

## stabilité 

set.seed(123)
mse.stab <- stability.arbre(Data, Z, p)    # arbre 
mse.stab.bg <- stability.bagging(Data, Z, p, B=B)  # bagging

set.seed(524)
mse.stab.rf <- stability.rf(Data, Z, p=p, m =1:p,B=B, fold=3) # random forest

### Table 4  

stabilite.table <- rbind(mse.stab, mse.stab.bg, mse.stab.rf)
rownames(stabilite.table) <- c("Arbre", "Bagging", "Random Forest")
stabilite.table
#############################################################################
################################## real data ###############################
############################################################################

#X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20+X21+X22+X23+X24+X25+X26+X27


stability.rf <- function(Data, Z, p, m, B, fold=3){
  
  data <- Data
  train.1 <- data$train 
  test.1 <- data$test
  M <- m.cv(rbind(train.1,test.1[1:nrow(train.1),]), m,B, fold)
  rf.optim.1 <- randomForest(y ~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20+X21+X22+X23+X24+X25+X26+X27,
                             data = train.1,      
                             mtry = M, ntree = B)
  Predict.1 <- predict(rf.optim.1, newdata = test.1[1:nrow(train.1),])
  l.1 <- (test.1$y[1:nrow(train.1)]-Predict.1)^2
  l.2 <- c()
  for(i in 1:nrow(data$train)){
    rf.optim.2 <- NULL; Predict.2 <- NULL
    
    train.2 <- data$train
    test.2 <- data$test
    for(ind in 1:(p+1)){
      train.2[i,ind] <- Z[ind]   #remplace zi par z'
    }
    M <- m.cv(rbind(train.2,test.1[1:nrow(train.1),]), m, B, fold)
    rf.optim.2 <- randomForest(y ~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20+X21+X22+X23+X24+X25+X26+X27,
                               data = train.2,      
                               mtry = M, ntree = B)
    Predict.2 <- predict(rf.optim.2, newdata = test.1[1:nrow(train.1),])
    l.2 <- c(l.2, ((test.1$y[1:nrow(train.1)]-Predict.2)^2))
  }
  return( sum(l.1-l.2))
}

library(readxl)

Residential_Building_Data_Set <- read_excel("Residential-Building-Data-Set.xlsx",
                                            sheet = "Data", skip =1)

#View(Residential_Building_Data_Set)
#colnames(Residential_Building_Data_Set)

Data <- Residential_Building_Data_Set %>%
  select(c(5:31,109))

p <- ncol(Data)-1
n <- nrow(Data)

colnames(Data) <- c(paste("X",1:p, sep = ""),"y")

## création de données d'entrainement et de test 

dataset <- function(data){
  data <- data%>%
    mutate(id = row_number())
  train <- data %>%
    sample_frac(.35)  # données d'entrainement
  test <- anti_join(data, train, by = 'id')
  return(list(train = train,test=test))
}

datase <- dataset(Data)  # base à utiliser

m <- 1:p
B <- 50
fold <- 3

###########  Calcul des MSE (Table 2) ###################

J <- 1

set.seed(123)
Arbre_d <- mse.comparaison(datase, method = "arbre", J=J)

set.seed(123)
Random.Forest_d <- mse.comparaison(datase, method = "rf", m, B=B, fold = 3, J=J)

set.seed(123)
Bagging_d <- mse.comparaison(datase, method = "bagging", B = B,J=J)

tableau <- rbind(Arbre_d$MSE, Bagging_d$MSE, Random.Forest_d$MSE)
rownames(tableau) <- c("Arbre", "Bagging", "Random Forest")
colnames(tableau) <- c("MSE")
tableau

##########################  Stabilité ########################

## calcule de z'
set.seed(524)
X <- replicate(p,runif(1))
Y <- sum(X)*100 + rnorm(1)
Z <- c(X,Y)

## stabilité 

set.seed(524)
mse.stab <- stability.arbre(datase, Z, p)    # arbre 

set.seed(524)
mse.stab.bg <- stability.bagging(datase, Z, p, B=B)  # bagging

set.seed(524)
mse.stab.rf <- stability.rf(datase, Z, p=p, m =1:p,B=B, fold=3) # random forest

### Table 4  

stabilite.table <- rbind(mse.stab, mse.stab.bg, mse.stab.rf)
rownames(stabilite.table) <- c("Arbre", "Bagging", "Random Forest")
stabilite.table



