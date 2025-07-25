
#' MEGB algorithm
#'
#' MEGB is an adaptation of the gradient boosting regression method to longitudinal data similar to MERF developed by Hajjem et. al. (2014) <doi:10.1080/00949655.2012.741599> which was implemented by Capitaine et. al. (2020) <doi:10.1177/0962280220946080>.
#' The algorithm estimates the parameters of a semi-parametric mixed-effects model: \deqn{Y_i(t)=f(X_i(t))+Z_i(t)\beta_i+\epsilon_i}
#' with \eqn{Y_i(t)} the output at time \eqn{t} for the \eqn{i}th individual; \eqn{X_i(t)} the input predictors (fixed effects) at time \eqn{t} for the \eqn{i}th individual;
#' \eqn{Z_i(t)} are the random effects at time \eqn{t} for the \eqn{i}th individual;
#'  \eqn{\epsilon_i} is the residual error.
#'
#' @param X [matrix]: A \code{N}x\code{p} matrix containing the \code{p} predictors of the fixed effects, column codes for a predictor.
#' @param Y [vector]: A vector containing the output trajectories.
#' @param id [vector]: Is the vector of the identifiers for the different trajectories.
#' @param Z [matrix]: A \code{N}x\code{q} matrix containing the \code{q} predictor of the random effects.
#' @param iter [numeric]: Maximal number of iterations of the algorithm. The default is set to \code{iter=100}
#' @param mtry [numeric]: Number of variables randomly sampled as candidates at each split. The default value is \code{p/3}.
#' @param ntree [numeric]: Number of trees to grow. This should not be set to too small a number, to ensure that every input row gets predicted at least a few times. The default value is \code{ntree=500}.
#' @param time [vector]: Is the vector of the measurement times associated with the trajectories in \code{Y},\code{Z} and \code{X}.
#' @param shrinkage [numeric]: a shrinkage parameter applied to each tree in the expansion. Also known as the learning rate or step-size reduction. The default value is set to 0.05.
#' @param interaction.depth [numeric]: The maximum depth of variable interactions: 1 builds an additive model, 2 builds a model with up to two-way interactions, etc. The default value is set to 1.
#' @param n.minobsinnode [numeric]: minimum number of observations (not total weights) in the terminal nodes of the trees. The default value is set to 5.
#' @param cv.folds [numeric]: Number of cross-validation folds to perform. If cv.folds>1 then gbm, in addition to the usual fit, will perform a cross-validation and calculate an estimate of generalization error returned in cv_error. The default value is set to 3.
#' @param delta [numeric]: The algorithm stops when the difference in log likelihood between two iterations is smaller than \code{delta}. The default value is set to O.O01
#' @param verbose [boolean]: If TRUE, MEGB will print out number of iterations to achieve convergence. Default is TRUE.
#'
#' @import gbm
#' @import stats
#' @return A fitted MEGB model which is a list of the following elements: \itemize{
#' \item \code{forest:} GBMFit obtained at the last iteration.
#' \item \code{random_effects :} Predictions of random effects for different trajectories.
#' \item \code{id_btilde:} Identifiers of individuals associated with the predictions \code{random_effects}.
#' \item \code{var_random_effects: } Estimation of the variance covariance matrix of random effects.
#' \item \code{sigma: } Estimation of the residual variance parameter.
#' \item \code{time: } The vector of the measurement times associated with the trajectories in \code{Y},\code{Z} and \code{X}.
#' \item \code{LL:} Log-likelihood of the different iterations.
#' \item \code{id: } Vector of the identifiers for the different trajectories.
#' \item \code{OOB: } OOB error of the fitted random forest at each iteration.
#' }
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' data <-simLong(n = 20,p = 6,rel_p = 6,time_points = 10,rho_W = 0.6, rho_Z=0.6,
#'               random_sd_intercept = sqrt(0.5),
#'               random_sd_slope = sqrt(3),
#'               noise_sd = 0.5,linear=TRUE)  # Generate the data composed by n=20 individuals.
#' # Train a MEGB model on the generated data. Should take ~ 7 seconds
#' megb <-   MEGB(X=as.matrix(data[,-1:-5]),Y=as.matrix(data$Y),
#' Z=as.matrix(data[,4:5]),id=data$id,time=data$time,ntree=500,cv.folds=3,verbose=TRUE)
#' megb$forest # is the fitted gradient boosting (GBMFit) (obtained at the last iteration).
#' megb$random_effects # are the predicted random effects for each individual.
#' plot(megb$LL,type="o",col=2) # evolution of the log-likelihood.
#' megb$OOB # OOB error at each iteration.
#'
#'
MEGB = function (X, Y, id, Z, iter = 100, ntree = 500, time, shrinkage=0.05, 
                 interaction.depth=1, n.minobsinnode=5, cv.folds=3, delta = 0.001,verbose=TRUE)
{
  
  
  q <- dim(Z)[2]
  nind <- length(unique(id))
  btilde <- matrix(0, nind, q)
  sigmahat <- 1
  Btilde <- diag(rep(1, q))
  epsilonhat <- rep(0, length(Y))
  id_btilde <- unique(id)
  Tiime <- sort(unique(time))
  sigma2 <- 1
  Vrai <- NULL
  inc <- 1
  OOB <- NULL
  
  
  for (i in 1:iter) {
    ystar <- rep(NA, length(Y))
    for (k in 1:nind) {
      indiv <- which(id == unique(id)[k])
      ystar[indiv] <- Y[indiv] - Z[indiv, , drop = FALSE] %*%
        btilde[k, ]
    }
    set.seed(i)
    forest <- gbm3::gbm(ystar ~ ., data = data.frame(X,ystar),distribution = "gaussian", 
                        n.trees = ntree, shrinkage = shrinkage,             
                        interaction.depth = interaction.depth, bag.fraction = 1, train.fraction = 1,  
                        n.minobsinnode = n.minobsinnode, cv.folds = cv.folds, keep.data = TRUE, 
                        verbose = FALSE)
    
    fhat <- as.numeric(predict(forest,n.trees=ntree,newdata=data.frame(X,ystar)))
    OOB[i] <- forest$train.error[ntree]
    for (k in 1:nind) {
      indiv <- which(id == unique(id)[k])
      V <- Z[indiv, , drop = FALSE] %*% Btilde %*%
        t(Z[indiv, , drop = FALSE]) + diag(as.numeric(sigmahat),
                                           length(indiv), length(indiv))
      btilde[k, ] <- Btilde %*% t(Z[indiv, , drop = FALSE]) %*%
        solve(V) %*% (Y[indiv] - fhat[indiv])
      epsilonhat[indiv] <- Y[indiv] - fhat[indiv] -
        Z[indiv, , drop = FALSE] %*% btilde[k, ]
    }
    sigm <- sigmahat
    sigmahat <- sig(sigma = sigmahat, id = id, Z = Z,
                    epsilon = epsilonhat, Btilde = Btilde)
    Btilde <- bay(bhat = btilde, Bhat = Btilde, Z = Z,
                  id = id, sigmahat = sigm)
    Vrai <- c(Vrai, logV(Y, fhat, Z, time, id, Btilde,
                         0, sigmahat))
    if (i > 1)
      inc <- abs((Vrai[i - 1] - Vrai[i])/Vrai[i -
                                                1])
    if (inc < delta) {
      if(verbose){print(paste0("stopped after ", i, " iterations."))}
      sortie <- list(forest = forest, random_effects = btilde,
                     var_random_effects = Btilde, sigma = sigmahat,
                     id_btilde = unique(id), LL = Vrai,
                     id = id, time = time, OOB = OOB)
      class(sortie) <- "MEGB"
      return(sortie)
    }
  }
  
}


#' Predict with longitudinal trees and random forests.
#'
#' @param object : a \code{longituRF} output of (S)MERF; (S)REEMforest; (S)MERT or (S)REEMtree function.
#' @param X [matrix]: matrix of the fixed effects for the new observations to be predicted.
#' @param Z [matrix]: matrix of the random effects for the new observations to be predicted.
#' @param id [vector]: vector of the identifiers of the new observations to be predicted.
#' @param time [vector]: vector of the time measurements of the new observations to be predicted.
#' @param ntree [numeric]: Number of trees to be used in prediction not less than number of trees used in the model object MEGB. The default value is \code{ntree=500}.
#' @param ... : low levels arguments.
#'
#' @import stats
#' @import gbm
#'
#' @return vector of the predicted output for the new observations.
#'
#' @export
#'
#' @examples 
#' set.seed(1)
#' data <-simLong(n = 20,p = 6,rel_p = 6,time_points = 10,rho_W = 0.6, rho_Z=0.6,
#'               random_sd_intercept = sqrt(0.5),
#'               random_sd_slope = sqrt(3),
#'               noise_sd = 0.5,linear=TRUE)  # Generate the data composed by n=20 individuals.
#' # Train a MEGB model on the generated data. Should take ~ 7 seconds
#' megb <-   MEGB(X=as.matrix(data[,-1:-5]),Y=as.matrix(data$Y),
#' Z=as.matrix(data[,4:5]),id=data$id,time=data$time,ntree=500,cv.folds=3,verbose=TRUE)
#' # Then we predict on the learning sample :
#' pred.MEGB <- predict(megb, X=as.matrix(data[,-1:-5]), Z=as.matrix(data[,4:5]),
#' id=data$id, time=data$time,ntree=500)
#' # Let's have a look at the predictions
#' # the predictions are in red while the real output trajectories are in blue:
#' par(mfrow=c(4,5),mar=c(2,2,2,2))
#' for (i in unique(data$id)){
#'   w <- which(data$id==i)
#'   plot(data$time[w],data$Y[w],type="l",col="blue")
#'   lines(data$time[w],pred.MEGB[w], col="red")
#' }
#'
predict.MEGB <- function(object,X,Z,id,time,ntree,...){
  dfX = data.frame(X)
  colnames(dfX) = colnames(X)
  n <- length(unique(id))
  id_btilde <- object$id_btilde
  f <- as.numeric(gbm3:::predict.GBMFit(object$forest,newdata=dfX,n.trees=ntree))
  Time <- object$time
  id_btilde <- unique(id)
  Ypred <- rep(0,length(id))
  id.app=object$id
  for (i in 1:length(unique(id))){
    w <- which(id==unique(id)[i])
    k <- which(id_btilde==unique(id)[i])
    Ypred[w] <- f[w] + Z[w,, drop=FALSE]%*%object$random_effects[k,]
  }
  return(Ypred)
  
}


#' Title
#'
#'
#' @import stats
#'
#' @keywords internal
sig <- function(sigma,id,Z, epsilon, Btilde){ #### fonction d'actualisation du param?tre de la variance des erreurs
  nind <- length(unique(id))
  Nombre <- length(id)
  sigm <- 0
  for (j in 1:nind){
    w <- which(id==unique(id)[j])
    V <- Z[w,, drop=FALSE]%*%Btilde%*%t(Z[w,, drop=FALSE])+diag(as.numeric(sigma),length(w),length(w))
    sigm <- sigm + t(epsilon[w])%*%epsilon[w] + sigma*(length(w)-sigma*(sum(diag(solve(V)))))
  }
  sigm <- sigm/Nombre
  return(sigm)
}

#' Title
#'
#' @import stats
#'
#' @keywords internal
bay <- function(bhat,Bhat,Z,id, sigmahat){ #### actualisation des param?tres de B
  nind <- length(unique(id))
  q <- dim(Z)[2]
  Nombre <- length(id)
  D <- 0
  for (j in 1:nind){
    w <- which(id==unique(id)[j])
    V <- Z[w,, drop=FALSE]%*%Bhat%*%t(Z[w,, drop=FALSE])+diag(as.numeric(sigmahat),length(w),length(w))
    D <- D+ (bhat[j,]%*%t(bhat[j,]))+ (Bhat- Bhat%*%t(Z[w,, drop=FALSE])%*%solve(V)%*%Z[w,, drop=FALSE]%*%Bhat)
  }
  D <- D/nind
  return(D)
}

#' Title
#'
#' @import stats
#'
#' @keywords internal
logV <- function(Y,f,Z,time,id,B,gamma,sigma){ # Maximization of variance of Y at M stage of EM
  Vraisem <- 0
  for (i in 1:length(unique(id))){
    w <- which(id==unique(id)[i])
    V <- Z[w,,drop=FALSE]%*%B%*%t(Z[w,,drop=FALSE])+diag(as.numeric(sigma),length(w),length(w))
    Vraisem <- Vraisem + log(det(V))+ t(Y[w]-f[w])%*%solve(V)%*%(Y[w]-f[w])
  }
  return(Vraisem)
}



######################## simulate longitudinal data#############################
################################################################################

## Function: simLong from MEGB package
## Simulates linear/nonlinear and low/high-dimensional longitudinal data.
## Parameters:
### - n: number of subjects
### - p: number of predictors
### - rel_p: number of relevant predictors
### - time_points: number of time points
### - rho: AR(1) correlation structure
### - random_sd_intercept: SD of random intercepts
### - random_sd_slope: SD of random slopes
### - noise_sd: SD of random noise
### - linear: boolean to determine if it is linear or nonlinear temporal relevant predictors.
### Requires the `MASS` package.



#### Prediction for longituRF: MERF & REEMForest

predict.longituRF <- function(object, X,Z,id,time,...){
  n <- length(unique(id))
  id_btilde <- object$id_btilde
  f <- predict(object$forest,X)
  Time <- object$time
  id_btilde <- unique(id)
  Ypred <- rep(0,length(id))
  id.app=object$id
  if (object$sto=="none"){
    for (i in 1:length(unique(id))){
      w <- which(id==unique(id)[i])
      k <- which(id_btilde==unique(id)[i])
      Ypred[w] <- f[w] + Z[w,, drop=FALSE]%*%object$random_effects[k,]
    }
    return(Ypred)
  }
  
  if (object$sto=="exp"){
    for (i in 1:length(unique(id))){
      w <- which(id==unique(id)[i])
      k <- which(id_btilde==unique(id)[i])
      om <- which(id.app==unique(id)[i])
      Ypred[w] <- f[w] + Z[w,, drop=FALSE]%*%object$random_effects[k,] + predict.exp(object$omega[om],Time[om],time[w], object$alpha)
    }
    return(Ypred)
  }
  
  if (object$sto=="fbm"){
    for (i in 1:length(unique(id))){
      w <- which(id==unique(id)[i])
      k <- which(id_btilde==unique(id)[i])
      om <- which(id.app==unique(id)[i])
      Ypred[w] <- f[w] + Z[w,, drop=FALSE]%*%object$random_effects[k,] + predict.fbm(object$omega[om],Time[om],time[w], object$Hurst)
    }
    return(Ypred)
  }
  
  for (i in 1:length(unique(id))){
    w <- which(id==unique(id)[i])
    k <- which(id_btilde==unique(id)[i])
    om <- which(id.app==unique(id)[i])
    Ypred[w] <- f[w] + Z[w,, drop=FALSE]%*%object$random_effects[k,] + predict.sto(object$omega[om],Time[om],time[w], object$sto)
  }
  return(Ypred)
}



###### Prediction function for glmmLasso


predict.glmmLasso = function(object,test_data){
  n <- length(unique(test_data$id))
  id = test_data$id
  random_effects = matrix(object$ranef,nrow=n,ncol=2,byrow=T)
  Z = as.matrix(test_data[,4:5])
  
  f <- as.numeric(as.matrix(test_data[,c(-1:-3,-5)])%*%object$coef)
  Time <- test_data$time
  id_btilde <- unique(id)
  Ypred <- rep(0,length(id))
  id.app=test_data$id
  for (i in 1:length(unique(id))){
    w <- which(id==unique(id)[i])
    k <- which(id_btilde==unique(id)[i])
    Ypred[w] <- f[w] + Z[w,, drop=FALSE]%*%random_effects[k,]
  }
  return(Ypred)
  
}

##################### Modelling

library(lme4)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(caret)
library(randomForestSRC)
library(randomForest)
library(gbm)
library(gpboost)
library(glmmLasso)
#library(LongituRF)
library(MEGB) ### install.packages("MEGB") or load the function above


# Define the model fit functions
fit_lm <- function(data, y_col="y", predictors=predictors) {
  formula <- as.formula(paste(y_col, "~", paste(predictors, collapse = "+")))
  lm(formula, data = data)
}

fit_lmer <- function(data, y_col="y", predictors=predictors, id_col="id") {
  formula <- as.formula(paste(y_col, "~", paste(predictors[-1:-2], collapse = "+"), "+ (1 |", id_col, ")"))
  lmer(formula, data = data)
}


fit_rf <- function(data, y_col="y", predictors=predictors) {
  formula <- as.formula(paste(y_col, "~", paste(predictors, collapse = "+")))
  rfsrc(formula, data = data,ntree = 500)
}

fit_rf2 <- function(data, y_col="y", predictors=predictors) {
  formula <- as.formula(paste(y_col, "~", paste(predictors, collapse = "+")))
  randomForest(formula, data = data,ntree = 500)
}


fit_gbm <- function(data, y_col="y", predictors=predictors) {
  formula <- as.formula(paste(y_col, "~", paste(predictors, collapse = "+")))
  gbm3::gbm(formula, data = data,distribution = "gaussian",
            n.trees = 500, shrinkage = 0.05,
            interaction.depth = 1, bag.fraction = 1, train.fraction = 1,
            n.minobsinnode = 2, cv.folds = 3, keep.data = TRUE,
            verbose = FALSE)
}

fit_merf <- function(data, y_col="y", predictors=predictors) {
  LongituRF::MERF(X=as.matrix(data[,-1:-5]),Y=as.matrix(data[,3]),
                  Z=as.matrix(data[,4:5]),id=data$id,
                  time=data$time,ntree=500,sto="none")
}

fit_megb <- function(data, y_col="y", predictors=predictors) {
  MEGB(X=as.matrix(data[,-1:-5]),Y=as.matrix(data[,3]),
       Z=as.matrix(data[,4:5]),id=data$id,
       time=data$time,ntree=500)
}


fit_reemf <- function(data, y_col="y", predictors=predictors) {
  LongituRF::REEMforest(X=as.matrix(data[,-1:-5]),Y=as.matrix(data[,3]),
                        Z=as.matrix(data[,4:5]),id=data$id,
                        time=data$time,ntree=500,sto="none",
                        mtry=floor(ncol(data[,-1:-5])/3))
}


fit_gpboost <- function(data, y_col="y", predictors=predictors, id_col="id") {
  gp_model <- GPModel(
    group_data   = data[,1],
    group_rand_coef_data = as.matrix(data[,4:5]),
    likelihood           = "gaussian",
    ind_effect_group_rand_coef = c(1,1)
  )
  
  gpboost_model = gpboost(
    data         = as.matrix(data[, -1:-5]),
    label        = data$y,
    gp_model     = gp_model,
    nrounds      = 200,
    learning_rate= 0.05,
    max_depth    = 3,
    objective    = "regression_l2",
    verbose = -1
  )
  
  gpboost_model
  
}

fit_glmmlasso <- function(data, y_col="y", predictors=predictors, id_col="id") {
  data$id = factor(data$id)
  fixed_effects <- paste(predictors[-1:-2], collapse = "+")
  fixed_formula <- as.formula(paste(y_col, "~", fixed_effects))
  glmmLasso(fix = fixed_formula, rnd = list(id=~1+RandomSlope), data = data, 
            lambda = 1,family=gaussian(link="identity"),switch.NR = TRUE,
            control=list(center = FALSE, standardize = FALSE),
            final.re=TRUE)
  
}





################################################################################
##################### MEGB: paper###############################################
################################################################################


lmerImp = function(lmerModel){
  mest = fixef(lmerModel)[-1]
  vest = sqrt(diag(vcov(lmerModel)))[-1]
  pest = 2*(1-pnorm(abs(mest/vest)))
  
  dest = data.frame(mest,vest,pest)
  dest[order(dest$pest, decreasing = FALSE), ]
}


cross_validate_blocked <- function(data, id_col, time_col, y_col, k = 5, repeats = 1,
                                   model_fit_fn, predictors, rel_p) {
  # Arguments:
  # data: A data frame containing the longitudinal data.
  # id_col: The name of the column identifying subjects.
  # time_col: The name of the time column (if applicable).
  # y_col: The name of the outcome column.
  # k: Number of folds for k-fold cross-validation (used for "k-Fold" and "Stratified").
  # repeats: Number of repetitions for repeated cross-validation.
  # model_fit_fn: A function to fit the model, accepting training data and returning predictions.
  # predictors: A vector of predictor column names.
  
  library(dplyr)
  library(caret)
  library(purrr)
  
  #data = data; id_col = "id";time_col = "time";y_col = "y";repeats = 1;k = 10;model_fit_fn = fit_gpboost;
  #predictors = predictors; rel_p = 2
  
  
  # Relevant Fixed Effect Predictors for variable importance
  
  fixedPred = colnames(data[,-1:-5])[1:rel_p]
  total_irrelevant <- ifelse((ncol((data[,-1:-5])) - rel_p)==0,1e+4,ncol((data[,-1:-5])) - rel_p) # to avoid 0/0
  
  #set.seed(123)
  
  # Initialize results storage
  all_results <- list()
  
  results <- replicate(repeats, {
    data <- arrange(data, !!sym(time_col))
    block_size <- ceiling(nrow(data) / k)
    blocks <- split(data, (seq_len(nrow(data)) - 1) %/% block_size)
    
    map(seq_along(blocks), function(i) {
      # Track time for each fold
      fold_start_time <- Sys.time()
      #i=1
      train_data <- bind_rows(blocks[-i])
      test_data <- blocks[[i]]
      
      # Fit model and predict
      model <- model_fit_fn(data = train_data, predictors = predictors)
      if (class(model)[1] == "lmerMod") {
        
        preds <- predict(model, test_data, allow.new.levels = TRUE)
        
        vimp_internal = lmerImp(model)
        vimpVar = sum((rownames(vimp_internal)[1:rel_p] %in% fixedPred))/rel_p
        fprVar = sum(!(rownames(vimp_internal)[1:rel_p] %in% fixedPred))/total_irrelevant
        
        
      } else {
        if (class(model)[1] == "rfsrc") {
          preds <- predict(model, test_data)$predicted
          
          vimp_internal = sort(vimp(model)$importance, decreasing = TRUE)
          vimpVar = sum((names(vimp_internal[1:rel_p]) %in% fixedPred))/rel_p
          fprVar = sum(!(names(vimp_internal[1:rel_p]) %in% fixedPred))/total_irrelevant
          
        } else {
          if(class(model)[1] == "longituRF") {
            preds = predict(model, X = as.matrix(test_data[,-1:-5]),
                            Z = as.matrix(test_data[,4:5]), id = test_data$id,
                            time = test_data$time)
            
            vimp_data = model$forest$importance
            vimp_internal = vimp_data[order(vimp_data[,2], decreasing = TRUE), ]
            vimpVar = sum((rownames(vimp_internal)[1:rel_p] %in% fixedPred))/rel_p
            fprVar = sum(!(rownames(vimp_internal)[1:rel_p] %in% fixedPred))/total_irrelevant
            
          } else {
            if(class(model)[1] == "GBMFit") {
              preds <- predict(model, test_data, n.trees = 500)
              
              vimp_internal = summary(model,plot_it=FALSE)
              vimpVar = sum((rownames(vimp_internal)[1:rel_p] %in% fixedPred))/rel_p
              fprVar = sum(!(rownames(vimp_internal)[1:rel_p] %in% fixedPred))/total_irrelevant
              
            } else {
              if(class(model)[1] == "longituGBM") {
                preds <- predict(model, dfX = test_data[,-1:-5], X = as.matrix(test_data[,-1:-5]),
                                 Z = as.matrix(test_data[,4:5]), id = test_data$id,
                                 time = test_data$time, n.trees = 500)
                
                vimp_internal = summary(model$forest,plot_it=FALSE)
                vimpVar = sum((rownames(vimp_internal)[1:rel_p] %in% fixedPred))/rel_p
                fprVar = sum(!(rownames(vimp_internal)[1:rel_p] %in% fixedPred))/total_irrelevant
                
              }else{ 
                
                if(class(model)[1] == "gpb.Booster") {
                  preds <- predict(model,data= as.matrix(test_data[, -1:-5]),group_data_pred = test_data[,1], 
                                   group_rand_coef_data_pred = as.matrix(test_data[,4:5]))$response_mean
                  
                  vimp_internal = gpboost::gpb.importance(model)$Feature
                  vimpVar = sum(vimp_internal[1:rel_p] %in% fixedPred)/rel_p
                  fprVar = sum(!(vimp_internal[1:rel_p] %in% fixedPred))/total_irrelevant
                  
                }else{
                  
                  if(class(model)[1] == "glmmLasso") {
                    preds <- predict(model,test_data)
                    
                    vimp_internal = sort(abs(model$coefficients[-1]),decreasing=TRUE)
                    vimpVar = sum((names(vimp_internal[1:rel_p]) %in% fixedPred))/rel_p
                    fprVar = sum(!(names(vimp_internal[1:rel_p]) %in% fixedPred))/total_irrelevant
                    
                  }else {
                    preds <- predict(model, test_data)
                    
                  }
                }
              }
            }
          }
        }
      }
      mse <- mean((test_data[[y_col]] - preds)^2)
      
      # Track fold computation time
      fold_end_time <- Sys.time()
      fold_time <- as.numeric(fold_end_time - fold_start_time, units = "secs")
      
      # Return both MSE and fold computation time
      return(list(mse = mse, fold_time = fold_time,vimpVar=vimpVar,fprVar=fprVar))
    })
  }, simplify = FALSE)
  
  # Combine results across repetitions
  results_combined <- map(results, ~bind_rows(.x))
  
  # Format and return
  list(
    mse = map(results_combined, ~.x$mse),
    fold_times = map(results_combined, ~.x$fold_time),
    vimpVar = map(results_combined, ~.x$vimpVar),
    fprVar = map(results_combined, ~.x$fprVar)
  )
}




# Data generation for different p
generate_data <- function(p) {
  
  set.seed(12345)
  simLong(n = 20, p = p, rel_p = 2, time_points = 10, rho = 0.6,
          random_sd_intercept = sqrt(.5), random_sd_slope = sqrt(3),
          noise_sd = 0.5, linear = FALSE)
}

# Initialize results for p = 6, 100, 2000
p_values <- c(6,170,2000)
#p_values <- c(6)
results <- list()

for (p in p_values) {
  data <- generate_data(p)
  names(data)[names(data) == "Y"] <- "y"
  predictors <- c("time", "RandomSlope", colnames(data[, -1:-5]))
  
  # Cross-validation for different models
  if(p<=170){
    
    model_fits <- list(
      lmer = fit_lmer,
      rf = fit_rf,
      gbm = fit_gbm,
      #merf = fit_merf,
      #reemf = fit_reemf,
      megb = fit_megb,
      gpboost = fit_gpboost,
      glmmlasso = fit_glmmlasso
    )
  } else{
    
    model_fits <- list(
      rf = fit_rf,
      gbm = fit_gbm,
      #merf = fit_merf,
      #reemf = fit_reemf,
      megb = fit_megb,
      gpboost = fit_gpboost,
      glmmlasso = fit_glmmlasso)
    
  }
  for (model_name in names(model_fits)) {
    model_fit_fn <- model_fits[[model_name]]
    
    res <- cross_validate_blocked(
      data = data,
      id_col = "id",
      time_col = "time",
      y_col = "y",
      repeats = 10,
      k = 10,
      model_fit_fn = model_fit_fn,
      predictors = predictors, rel_p = 2
    )
    
    model_results <- data.frame(
      p = p,
      model = model_name,
      fold = seq_along(unlist(res$mse)),
      mse = unlist(res$mse),
      time = unlist(res$fold_times),
      vimpVar = 100*unlist(res$vimpVar),
      fprVar = 100*unlist(res$fprVar)
    )
    
    results[[paste(p, model_name, sep = "_")]] <- model_results
    cat("Processed p =", p, "for model =", model_name, "...\n")
  }
}

# Combine results

long_results = do.call(rbind, results)

## save and recall results for future use

#library(here)
#here() starts at C:/Users/USER/Documents
#saveRDS(long_results, file = here("long_results.rds"))
#saveRDS(long_results, file = here("long_results_nonlinear.rds"))
#saveRDS(long_results, file = here("long_results_nonlinear2.rds"))
#write.csv(long_results,"long_results_megb_paper.csv")

##library(here)

#long_results = readRDS(here("long_results.rds"))
#long_results = readRDS(here("long_results_nonlinear.rds"))


#long_results1 = readRDS(here("long_results_nonlinear.rds"))
#long_results2 = readRDS(here("long_results_nonlinear2.rds"))

#long_results_m = long_results1

#long_results_m$model[long_results_m$model == "glmmlasso"] <- long_results2$model[long_results2$model == "glmmlasso"]

long_results = long_results1

long_results$Dimension = factor(ifelse(long_results$p==6,1,
                                       ifelse(long_results$p==170,2,3)),
                                labels=c("Low-Dimensional (p = 6)",
                                         "Medium-Dimensional (p = 170)",
                                         "High-Dimensional (p = 2000)"))

# Visualization
library(ggplot2)

# MSE plot
long_results <- long_results %>%
  group_by(Dimension,model) %>%
  mutate(median_mse = mean(mse)) %>%
  ungroup()

# Reorder the 'model' factor based on median MSE within each p
long_results <- long_results %>%
  mutate(model = reorder(model, median_mse, decreasing=TRUE))

# Visualization with sorted MSE


ggplot(long_results, aes(x = model, y = mse, fill = model)) +
  geom_boxplot() +
  facet_wrap(~Dimension, scales = "free") +
  labs(
    title = "",
    x = "Models",
    y = "Mean Squared Error (MSE)",
    fill = "Model"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 9, face = "bold")
  )




# Computational time plot
#long_results <- do.call(rbind, results)

long_results$Dimension = factor(ifelse(long_results$p==6,1,
                                       ifelse(long_results$p==170,2,3)),
                                labels=c("Low-Dimensional (p = 6)",
                                         "Medium-Dimensional (p = 170)",
                                         "High-Dimensional (p = 2000)"))
# Visualization
library(ggplot2)

# Time plot
long_results <- long_results %>%
  group_by(p, model) %>%
  mutate(median_time = mean(time)) %>%
  ungroup()

# Reorder the 'model' factor based on median Time within each p
long_results <- long_results %>%
  mutate(model = reorder(model, median_time, decreasing=TRUE))

# Visualization with sorted MSE

ggplot(long_results, aes(x = model, y = time, fill = model)) +
  geom_boxplot() +
  facet_wrap(~Dimension, scales = "free") +
  labs(
    title = "",
    x = "Models",
    y = "Computational Time (seconds)",
    fill = "Model"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 9, face = "bold")
  )


# Variable Importance plot
#long_results <- do.call(rbind, results)
long_results$Dimension = factor(ifelse(long_results$p==6,1,
                                       ifelse(long_results$p==170,2,3)),
                                labels=c("Low-Dimensional (p = 6)",
                                         "Medium-Dimensional (p = 170)",
                                         "High-Dimensional (p = 2000)"))
# Visualization
library(ggplot2)

# Time plot
long_results <- long_results %>%
  group_by(p, model) %>%
  mutate(median_vimp = median(vimpVar)) %>%
  ungroup()

# Reorder the 'model' factor based on median Time within each p
long_results <- long_results %>%
  mutate(model = reorder(model, median_vimp, decreasing=TRUE))

# Visualization with sorted MSE

ggplot(long_results, aes(x = model, y = vimpVar, fill = model)) +
  geom_boxplot() +
  facet_wrap(~Dimension, scales = "free") +
  labs(
    title = "CV Variable Importance for Different Models and p",
    x = "Models",
    y = "Variable Importance (%)",
    fill = "Model"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 9, face = "bold")
  )


aggregate(long_results$mse,list(long_results$model,long_results$p),
          function(x){paste0(round(mean(x),2),"(",round(sd(x),3),")")})

aggregate(long_results$time,list(long_results$model,long_results$p),
          function(x){paste0(round(mean(x),2),"(",round(sd(x),3),")")})

aggregate(long_results$vimpVar,list(long_results$model,long_results$p),
          function(x){paste0(round(mean(x),2),"(",round(sd(x),3),")")})

aggregate(long_results$fprVar,list(long_results$model,long_results$p),
          function(x){paste0(round(mean(x),2),"(",round(sd(x),3),")")})







################################################################################
############################# Real dataset analysis ############################
################################################################################

library(ggplot2)
library(dplyr)


#fdata2 <- read.csv("C:/Users/USER/Documents/maternal_fetal_rna_dataset.csv")

fdata2 <- read.csv("./maternal_fetal_rna_dataset.csv")

fdata2 <- fdata2 %>% rename(sample = id)

fdata2 <- fdata2 %>% mutate(sample=as.factor(sample))
fdata2 <- fdata2 %>% mutate(time = factor(fdata2$time, labels = c(1:4), 
                    levels = c("Trimester 1",
                               "Trimester 2","Trimester 3","Post Partum")))


# Create a spaghetti plot with GAM smoothing
ggplot(fdata2, aes(x = as.numeric(time), y = y, group = sample)) +
  #geom_line(alpha = 0.6, size = 1) +  # Spaghetti plot for individual samples
  geom_point(size = 2, alpha = 0.8) +  # Add points for each observation
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4), se = FALSE, linetype = "solid", size = 1, color = "blue") +  # GAM smoothing
  facet_wrap(~ sample, ncol = 4, scales = "free_y", labeller = "label_both") +  # Facet for each gene/sample
  theme_bw() +
  labs(title = "",
       x = "Pregnancy Stage",
       y = "Fetal RNA Score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_x_discrete(limits=c("Trimester 1",
                            "Trimester 2","Trimester 3","Post Partum"))


# Plot the Fetal RNA Score across time points

fdata2 <- fdata2 %>% mutate(time = factor(time, labels = c("Trimester 1",
                                                 "Trimester 2","Trimester 3","Post Partum")))
set.seed(123) # for reproducibility in geom_jitter
ggplot(fdata2, aes(x = time, y = y)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(aes(col = sample),width = 0.2, size = 2, alpha = 0.7) +
  #geom_text(aes(label = sample), vjust = -1, size = 4, color = "black") +
  labs(title = "",
       x = "Pregnancy Stage",
       y = "Fetal RNA Score") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


################################################################################
######################## Model Analysis ########################################
################################################################################


fdata2 <- fdata2 %>% rename(id = sample)

data.select2 <-  fdata2


######## KFold Cross-validation function


cross_validate_kFold <- function(data, id_col, time_col, y_col, k = 5, repeats = 1,
                                 model_fit_fn, predictors, rel_p) {
  # Arguments:
  # data: A data frame containing the longitudinal data.
  # id_col: The name of the column identifying subjects.
  # time_col: The name of the time column (if applicable).
  # y_col: The name of the outcome column.
  # k: Number of folds for k-fold cross-validation (used for "k-Fold" and "Stratified").
  # repeats: Number of repetitions for repeated cross-validation.
  # model_fit_fn: A function to fit the model, accepting training data and returning predictions.
  # predictors: A vector of predictor column names.
  
  library(dplyr)
  library(caret)
  library(purrr)
  
  data$time = as.numeric(data$time)
  
  # Prepare subject-level data
  subjects <- unique(data[[id_col]])
  
  # Relevant Fixed Effect Predictors for variable importance
  
  fixedPred = rel_p
  rel_p_n = length(fixedPred)
  
  #set.seed(123)
  
  # Initialize results storage
  all_results <- list()
  
  set.seed(123)
  results <- replicate(repeats, {
    folds <- createFolds(subjects, k = k, list = TRUE)
    map(folds, function(fold) {
      #fold=folds$Fold1
      fold_start_time <- Sys.time()
      train_data <- filter(data, !(!!sym(id_col) %in% fold))
      test_data <- filter(data, !!sym(id_col) %in% fold)
      
      # Fit model and predict
      
      #predictors_preg = colnames(data.select2[,-1:-4])
      #model <- fit_gbm(data = data.select2, predictors = predictors_preg)
      
      model <- model_fit_fn(data = train_data, predictors = predictors)
      
      if (class(model)[1] == "lmerMod") {
        preds <- predict(model, test_data, allow.new.levels = TRUE)
        
        vimp_internal = lmerImp(model)
        #vimpVar = sum((rownames(vimp_internal)[1:rel_p_n] %in% fixedPred))/rel_p_n
        vimpVar = rownames(vimp_internal)[1:rel_p_n]
      } else {
        if (class(model)[1] == "rfsrc") {
          preds <- predict(model, test_data)$predicted
          
          vimp_internal = sort(vimp(model)$importance, decreasing = TRUE)
          #vimpVar = sum((names(vimp_internal[1:rel_p_n]) %in% fixedPred))/rel_p_n
          vimpVar = names(vimp_internal[1:rel_p_n])
        } else {
          if(class(model)[1]=="longituRF"){
            preds = predict(model,X=as.matrix(test_data[,-1:-3]),
                            Z=as.matrix(test_data[,3]),id=test_data$id,
                            time=test_data$time)
            vimp_data = model$forest$importance
            vimp_internal = vimp_data[order(vimp_data[,2], decreasing = TRUE), ]
            #vimpVar = sum((rownames(vimp_internal)[1:rel_p_n] %in% fixedPred))/rel_p_n
            vimpVar = rownames(vimp_internal)[1:rel_p_n]
            
          }else{
            if(class(model)[1]=="GBMFit"){
              preds <- predict(model, test_data,n.trees = 500)
              
              vimp_internal = summary(model,plot_it=FALSE)
              #vimpVar = sum((rownames(vimp_internal)[1:rel_p_n] %in% fixedPred))/rel_p_n
              vimpVar = rownames(vimp_internal)[1:rel_p_n]
              
            }else{
              if(class(model)[1]=="longituGBM"){
                preds <- predict(model, dfX=test_data[,-1:-3],X=as.matrix(test_data[,-1:-3]),
                                 Z=as.matrix(test_data[,3]),id=test_data$id,
                                 time=test_data$time,n.trees = 500)
                vimp_internal = summary(model$forest,plot_it=FALSE)
                #vimpVar = sum((rownames(vimp_internal)[1:rel_p_n] %in% fixedPred))/rel_p_n
                vimpVar = rownames(vimp_internal)[1:rel_p_n]
              }else{ 
                
                if(class(model)[1] == "gpb.Booster") {
                  preds <- predict(model,data= as.matrix(test_data[, -1:-3]),group_data_pred = test_data[,1]
                  )$response_mean
                  
                  vimp_internal = gpboost::gpb.importance(model)$Feature
                  #vimpVar = sum(vimp_internal[1:rel_p_n] %in% fixedPred)/rel_p_n
                  vimpVar = vimp_internal[1:rel_p_n]
                  
                }else{
                  
                  if(class(model)[1] == "glmmLasso") {
                    preds <- predict(model,test_data)
                    
                    vimp_internal = sort(abs(model$coefficients[-1]),decreasing=TRUE)
                    #vimpVar = sum((names(vimp_internal[1:rel_p_n]) %in% fixedPred))/rel_p_n
                    vimpVar = names(vimp_internal[1:rel_p_n])
                    
                  }
                  else{
                    preds <- predict(model, test_data)
                    vimp_internal = summary(model$forest,plot_it=FALSE)
                    #vimpVar = sum((rownames(vimp_internal)[1:rel_p_n] %in% fixedPred))/rel_p_n
                    vimpVar = rownames(vimp_internal)[1:rel_p_n]
                    
                  }
                }
              }
            }
          }
        }
      }
      mse <- mean((test_data[[y_col]] - preds)^2)
      
      # Track fold computation time
      fold_end_time <- Sys.time()
      fold_time <- as.numeric(fold_end_time - fold_start_time, units = "secs")
      
      # Return both MSE and fold computation time
      return(list(mse = mse, fold_time = fold_time,vimpVar=vimpVar))
    })
  }, simplify = FALSE)
  
  
  # Combine results across repetitions
  results_combined <- map(results, ~bind_rows(.x))
  
  # Format and return
  list(
    mse = map(results_combined, ~.x$mse),
    fold_times = map(results_combined, ~.x$fold_time),
    vimpVar = map(results_combined, ~.x$vimpVar)
  )
  
}





######################################### model fit functions

fit_lmer2 <- function(data, y_col="y", predictors=predictors, id_col="id") {
  formula <- as.formula(paste(y_col, "~", paste(predictors, collapse = "+"), "+ (1 |", id_col, ")"))
  lmer(formula, data = data)
}



fit_gbm <- function(data, y_col="y", predictors=predictors) {
  formula <- as.formula(paste(y_col, "~", paste(predictors, collapse = "+")))
  gbm3::gbm(formula, data = data,distribution = "gaussian",
            n.trees = 500, shrinkage = 0.05,
            interaction.depth = 1, bag.fraction = 1, train.fraction = 1,
            n.minobsinnode = 2, cv.folds = 3, keep.data = TRUE,
            verbose = FALSE)
}

fit_merf2 <- function(data, y_col="y", predictors=predictors) {
  LongituRF::MERF(X=as.matrix(data[,predictors]),Y=as.matrix(data[,"y"]),
                  Z=as.matrix(data[,"Z"]),id=data$id,
                  time=data$time,ntree=500,sto="none")
}


fit_megb3 <- function(data, y_col="y", predictors=predictors) {
  MEGB(X=as.matrix(data[,predictors]),Y=as.matrix(data[,"y"]),
       Z=as.matrix(data[,"Z"]),id=data$id,
       time=data$time,ntree=500)
}

fit_reemf2 <- function(data, y_col="y", predictors=predictors) {
  LongituRF::REEMforest(X=as.matrix(data[,predictors]),Y=as.matrix(data[,"y"]),
                        Z=as.matrix(data[,"Z"]),id=data$id,
                        time=data$time,ntree=500,sto="none",
                        mtry=floor(ncol(data[,predictors])/3))
}


fit_gpboost2 <- function(data, y_col="y", predictors=predictors, id_col="id") {
  gp_model <- GPModel(
    group_data   = data[,1],
    likelihood           = "gaussian"
  )
  
  gpboost_model = gpboost(
    data         = as.matrix(data[, -1:-3]),
    label        = data$y,
    gp_model     = gp_model,
    nrounds      = 200,
    learning_rate= 0.05,
    max_depth    = 3,
    objective    = "regression_l2",
    verbose = -1
  )
  
  gpboost_model
  
}

fit_glmmlasso2 <- function(data, y_col="y", predictors=predictors, id_col="id") {
  data$id = factor(data$id)
  fixed_effects <- paste(predictors, collapse = "+")
  fixed_formula <- as.formula(paste(y_col, "~", fixed_effects))
  glmmLasso(fix = fixed_formula, rnd = list(id=~1), data = data, 
            lambda = 0.01,family=gaussian(link="identity"),switch.NR = TRUE,
            control=list(center = TRUE, standardize = FALSE),
            final.re=FALSE)
  
}


p_dim <- ncol(data.select2)-4
results <- list()
for (p in p_dim) {
  data <-data.select2
  names(data)[names(data) == "Y"] <- "y"
  
  predictors <- colnames(data.select2[,-1:-4])
  
  # Cross-validation for different models
  
  model_fits <- list(
    rf = fit_rf,
    gbm = fit_gbm,
    merf = fit_merf2,
    reemf = fit_reemf2,
    megb = fit_megb3,
    gpboost = fit_gpboost2
  )
  
  
  for (model_name in names(model_fits)) {
    model_fit_fn <- model_fits[[model_name]]
    
    res <- cross_validate_kFold(
      data = data,
      id_col = "id",
      time_col = "time",
      y_col = "y",
      repeats = 10,
      k = 10,
      model_fit_fn = model_fit_fn,
      predictors = predictors, rel_p = top_fetal_genes2
    )
    
    model_results <- data.frame(
      p = p,
      model = model_name,
      fold = seq_along(unlist(res$mse)),
      mse = unlist(res$mse),
      time = unlist(res$fold_times),
      vimpVar = unlist(res$vimpVar)
    )
    
    results[[paste(p, model_name, sep = "_")]] <- model_results
    cat("Processed p =", p, "for model =", model_name, "...\n")
  }
}

# Combine results

long_results = do.call(rbind, results)

#saveRDS(long_results, file = here("long_results_reallife.rds"))

# Visualization
library(ggplot2)

# MSE plot
long_results <- long_results %>%
  group_by(p,model) %>%
  mutate(median_mse = mean(mse)) %>%
  ungroup()

# Reorder the 'model' factor based on median MSE within each p
long_results <- long_results %>%
  mutate(model = reorder(model, median_mse, decreasing=TRUE))

# Visualization with sorted MSE


p1 = ggplot(long_results, aes(x = model, y = mse, fill = model)) +
  geom_boxplot() +
  facet_wrap(~p, scales = "free",labeller="label_both") +
  labs(
    title = "",
    x = "Models",
    y = "Mean Squared Error (MSE)",
    fill = "Model"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 9, face = "bold")
  )



# Computational time plot
long_results <- do.call(rbind, results)

# Visualization
library(ggplot2)

# Time plot
long_results <- long_results %>%
  group_by(p, model) %>%
  mutate(median_time = mean(time)) %>%
  ungroup()

# Reorder the 'model' factor based on median Time within each p
long_results <- long_results %>%
  mutate(model = reorder(model, median_time, decreasing=TRUE))

# Visualization with sorted MSE

p2=ggplot(long_results, aes(x = model, y = time, fill = model)) +
  geom_boxplot() +
  facet_wrap(~p, scales = "free",labeller="label_both") +
  labs(
    title = "",
    x = "Models",
    y = "Computational Time (seconds)",
    fill = "Model"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 9, face = "bold")
  )


cowplot::plot_grid(p1,p2)


# Variable Importance plot
long_results <- do.call(rbind, results)



aggregate(long_results$mse,list(long_results$model,long_results$p),
          function(x){paste0(round(median(x),2),"(",round(sd(x),3),")")})

aggregate(long_results$time,list(long_results$model,long_results$p),
          function(x){paste0(round(median(x),2),"(",round(sd(x),3),")")})


get_top_vimp_vars_all_models <- function(data, top_n = 9) {
  if (!"model" %in% names(data) || !"vimpVar" %in% names(data)) {
    stop("Data must contain 'model' and 'vimpVar' columns.")
  }
  
  # Get all unique models
  models <- unique(data$model)
  
  # Initialize empty list to store top vars per model
  top_vars_list <- lapply(models, function(mod) {
    # Count and sort variable importance for the model
    var_counts <- table(data$vimpVar[data$model == mod])
    sorted_vars <- sort(var_counts, decreasing = TRUE)[1:min(top_n, length(var_counts))]
    
    # Create a data frame with model, variable, and frequency
    df <- data.frame(
      model = mod,
      variable = names(sorted_vars),
      freq = as.integer(sorted_vars),
      rank = seq_along(sorted_vars)
    )
    return(df)
  })
  
  # Combine all into a single data frame
  top_vars_all <- do.call(rbind, top_vars_list)
  
  return(top_vars_all)
}

top_vars_df <- get_top_vimp_vars_all_models(long_results, top_n = 9)
print(top_vars_df)

# Reorder variables within each model for better plotting
top_vars_df <- top_vars_df %>%
  group_by(model) %>%
  mutate(variable = reorder(variable, -freq))

# Create faceted bar plot
ggplot(top_vars_df, aes(x = variable, y = freq, fill = model)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  facet_wrap(~ model, scales = "free_y") +
  labs(
    title = "Top Variable Importance Features by Model",
    x = "Variable",
    y = "Frequency"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )



# Step 1: Get top variables
top_vars_df <- get_top_vimp_vars_all_models(long_results, top_n = 9)

# Step 2: Identify how many models each variable appears in
var_model_counts <- top_vars_df %>%
  distinct(model, variable) %>%
  count(variable, name = "model_count")

# Step 3: Mark shared vs unique variables
top_vars_df <- top_vars_df %>%
  left_join(var_model_counts, by = "variable") %>%
  mutate(shared_status = ifelse(model_count > 1, "Shared", "Unique"))

# Step 4: Reorder variables within each model for clean plotting
top_vars_df <- top_vars_df %>%
  group_by(model) %>%
  mutate(variable = reorder(variable, -freq))

# Step 5: Plot with fill color based on shared_status
ggplot(top_vars_df, aes(x = variable, y = freq, fill = shared_status)) +
  geom_col(show.legend = TRUE) +
  coord_flip() +
  facet_wrap(~ model, scales = "free_y") +
  labs(
    title = "",
    x = "Top Gene Transcripts Selected",
    y = "Frequency",
    fill = "Gene Transcript Status"
  ) +
  scale_fill_manual(values = c("Shared" = "#D95F02", "Unique" = "#1B9E77")) +  # Custom colors
  theme_bw() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )





###### Transcripts curve

#plot_data = read.csv("C:/Users/USER/Documents/maternal_fetal_rna_dataset.csv")

plot_data = read.csv("./maternal_fetal_rna_dataset.csv")

top_megb = c("X7933084","X8142120", "X8019842", "X7940996","X7940216","X8042391", 
             "X7893518","X8149109","X8128123")


splot_data = plot_data[,c("id","y","time",top_megb)]

# Required libraries
library(ggplot2)
library(dplyr)
library(tidyr)


# Vector of selected variable names
genes <-c("X7933084","X8142120", "X8019842", "X7940996","X7940216","X8042391", 
          "X7893518","X8149109","X8128123")

# Convert to long format
long_data <- splot_data %>%
  pivot_longer(cols = all_of(genes), 
               names_to = "gene", 
               values_to = "expression")

# Plot with mean trajectory, 95% CI, raw data points, and faceting
ggplot(long_data, aes(x = time, y = expression)) +
  stat_summary(fun = mean, geom = "line", aes(group = 1), color = "blue", size = 1) +
  stat_summary(fun.data = mean_cl_normal, geom = "ribbon", 
               aes(group = 1), fill = "blue", alpha = 0.2) +
  geom_point(aes(color = time), position = position_jitter(width = 0.1), size = 1.5) +
  facet_wrap(~ gene, scales = "free_y") +
  theme_bw(base_size = 14) +
  labs(title = "",
       x = "Pregnancy Stage", y = "Expression Level") +
  theme(legend.position = "none",axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_x_discrete(limits=c("Trimester 1",
                            "Trimester 2","Trimester 3","Post Partum"))


