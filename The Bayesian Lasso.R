#install.packages(c('glmnet','MASS','lars','statmod'))
library(glmnet)
library(lars)
library(MASS)
library(statmod)
#----------------upload data----------------
data(diabetes)
attach(diabetes)
n <- nrow(diabetes)
X <- diabetes$x
y <- diabetes$y
X <- scale(X)/ sqrt(n - 1)#standarization
y <- as.numeric(scale(y,center=TRUE,scale = FALSE))#centered
dim(X) # 442 10
length(y) #442
n <- nrow(X) #the number of observations
p <- ncol(X)# the number of predictors
lambda_grid <- exp(seq(log(1000), log(0.001), length.out = 50))
#--------------- OLS Method-----------------
ols <- lm(y~X)# if you write data=diabetes, then it corresponding to you are doing diabetes$y, so the centered option doesn't work
beta_ols <- coef(ols)[-1]
#---------------Ridge Regression------------
#lambda <- exp(seq(log(1000),log(0.001),length.out=50))
#fit_ridge <- glmnet(x=X,y=y,alpha=0,lambda=lambda,intercept=FALSE,standardize = FALSE)
#beta_ridge <- fit_ridge$beta#10 50
cv_ridge <- glmnet(x=X,y=y,alpha=0,lambda=lambda_grid,intercept=FALSE,standardize=FALSE)

beta_ridge <- as.matrix(coef(cv_ridge))[-1,]#return the path of beta
#beta_ridge_best <- coef(cv_ridge,s='lambda.min')[-1,]

#---------------Lasso Regression------------
#fit_lasso <- glmnet(x+X,y=y,alpha=1,lambda=lambda,intercept = FALSE,standardize = FALSE)
#beta_lasso <- fit_lasso$beta#10 50
cv_lasso <- cv.glmnet(x=X,y=y,alpha=1,intercept=FALSE,standardize=FALSE)
beta_lasso <- as.matrix(coef(cv_lasso$glmnet.fit))[-1,]
beta_lasso_best <- coef(cv_lasso,s='lambda.min')[-1,]
cv_lasso$lambda.min
#---------------Bayesian Approach-----------
lasso_gibbs <- function(X,y,lambda,n_iter,burn){
  n <- nrow(X)
  p <- ncol(X)
  #setup
  XtX <- t(X)%*%X
  Xty <- t(X)%*%y
  beta <- rep(0,p)
  sigma2 <- var(y)
  tau2 <- rep(1,p)
  
  #saving
  
  beta_save <- matrix(NA,nrow=n_iter-burn,ncol=p)
  sigma2_save <- numeric(n_iter-burn)
  tau2_save <- matrix(NA,nrow=n_iter-burn,ncol=p)
  eps <- 1e-8
  
  for(iter in 1:n_iter){
    #update beta,follow a multivariate normal distribution(A^(-1)X^Ty,sigam2A^(-1)),where A=X^TX+D^(-1)
    D_inv <- diag(1/tau2,p,p)
    A <- XtX+D_inv
    L <- chol(A)
    # 计算均值：A * beta = Xty -> beta = A^-1 * Xty
    mean_beta <- backsolve(L, forwardsolve(t(L), Xty))
    # 抽取多元正态分布：beta = mean + sqrt(sigma2) * L^-1 * Z
    beta <- mean_beta + sqrt(sigma2) * backsolve(L, rnorm(p))
    #A_inv <- solve(A)
    #mean_beta <- A_inv%*%Xty
    #var_beta <- sigma2*A_inv
    #beta <- as.numeric(MASS::mvrnorm(1,mu=mean_beta,Sigma = var_beta))
    #update sigma2: IG(,)
    resid <- y-X%*%beta
    shape <- (n-1)/2+p/2
    rate <- sum(resid^2)/2+sum(beta^2/tau2)/2
    sigma2 <- 1/rgamma(1,shape=shape,rate=rate)
    
    #update tau2
    for(i in 1:p){
      
      mu_prime <- sqrt((lambda^2*sigma2)/(beta[i]^2+eps))
      inv_tau2_i <- statmod::rinvgauss(1,mean=mu_prime,shape=lambda^2)
      tau2[i] <- 1/inv_tau2_i
      
      
    }
    
    #saving
    if(iter>burn){
      beta_save[iter-burn,] <- beta
      sigma2_save[iter-burn] <- sigma2
      tau2_save[iter-burn,] <- tau2
      
    }
   
    
    
  }
  return(list(beta_samples=beta_save,
              sigma2_samples =sigma2_save,
              tau2_samples = tau2_save,
              beta_median=apply(beta_save, 2, median),
              beta_mean=colMeans(beta_save)
  ))
  
}
##define the lambda 



fit_ls <- lm(y ~ X - 1)#ols
beta_ls <- coef(fit_ls)
sigma2_ls <- sum(residuals(fit_ls)^2) / (n - p)
    
lambda_init <-p * sqrt(sigma2_ls) / sum(abs(beta_ls))

n <- 10
lambda_path <- numeric(n)

lambda_path[1] <- lambda_init

gibbs_results <- vector('list',n)


tol <- 1e-4


for(k in 1:n){
  cat("EM iteration:", k, " lambda =", lambda_path[k], "\n")
  fit_k <- lasso_gibbs(X,y,lambda = lambda_path[k],n_iter=11000,burn=1000)
  tau2_mean <- colMeans(fit_k$tau2_samples)
  lambda_new <- sqrt(2 * p / sum(tau2_mean))
  lambda_path[k + 1] <- lambda_new
  gibbs_results[[k]] <- list(
    lambda = lambda_path[k],
    fit    = fit_k
  )
  #if (abs(lambda_new - lambda_path[k]) < tol) {
   #cat("Converged at iteration", k, "\n")
   #break
  #}
}
#obtain the optimal lambda of bayesian
score_bayes <- numeric(length(gibbs_results))
for(k in 1:length(gibbs_results)){
  #posterior mean log-likelihood
  beta_mean_k <- colMeans(gibbs_results[[k]]$fit$beta_samples)
  mu_k <- as.numeric(X%*%beta_mean_k)
  s2_k <- mean(gibbs_results[[k]]$fit$sigma2_samples)
  score_bayes[k] <- sum(dnorm(y,mean=mu_k,sd=sqrt(s2_k),log=TRUE))
}

best_ind_bayes <- which.max(score_bayes)

bayes_best <- gibbs_results[[best_ind_bayes]]$fit$beta_median


beta_bayes <- matrix(NA,nrow=p,ncol=length(lambda_grid))
sigma2_bayes <- numeric(length(lambda_grid))

for(k in 1:length(lambda_grid)){
  cat('lambda=',lambda_grid[k],'\n')
 fit_b <- lasso_gibbs(X,y,lambda_grid[k],n_iter=11000,burn=1000)
  beta_bayes[,k] <- fit_b$beta_median

  sigma2_bayes[k] <- mean(fit_b$sigma2_samples)
}

#---------------------------------------

relative_l1 <- function(B)
{
  l1 <- colSums(abs(B))
  
  l1/max(l1)
}

x_lasso <- relative_l1(beta_lasso)

x_ridge <- relative_l1(beta_ridge)

x_bayes<- relative_l1(beta_bayes)
dim(beta_bayes)

#lasso vertical line
inx_lasso_line <- which.min(abs(cv_lasso$glmnet.fit$lambda-cv_lasso$lambda.min))
vline_lasso <- x_lasso[inx_lasso_line]
#ridge vertical line
#inx_ridge_line <- which.min(abs(cv_ridge$glmnet.fit$lambda-cv_ridge$lambda.min))

#vline_ridge <- x_ridge[inx_ridge_line]
#bayesian vertical line
vline_bayes <- sum(abs(bayes_best))/max(colSums(abs(beta_bayes)))

#-----------------Plot-------------------

par(mfrow = c(1, 3), mar = c(5, 4, 3, 2), oma = c(0, 1, 0, 0))


# panel (1):lasso
matplot(x_lasso,t(beta_lasso),
     type='l',lty=1,lwd=0.5,
     xlab='relative',
     ylab='Standarized Coefficients',
     main='(a):lasso'
     #col=1:p
          )




abline(v=vline_lasso,lty=2,lwd=0.5)

for (j in 1:p) {
  text(x = max(x_lasso), y = beta_lasso[j, ncol(beta_lasso)],
       labels = j, pos = 4, cex = 0.8)
}
#panel(2):bayesian

matplot(x_bayes,t(beta_bayes),type='l',lty=1,lwd=0.5,
        xlab='relative',
        ylab='Standarized Coefficients',
        main='(b):bayesian')
abline(v=vline_bayes,lty=2,lwd=0.5)

for (j in 1:p) {
  text(x = max(x_bayes), y = beta_bayes[j, ncol(beta_bayes)],
       labels = j, pos = 4, cex = 0.8)
}

#panel(3):ridge
#beta_ridge_scaled <- beta_ridge / apply(abs(beta_ridge), 1, max)
matplot(x_ridge,t(beta_ridge),type='l',lty=1,lwd=0.5,
        xlab='relative',
        ylab='Standarized Coefficients',
        main='(c):ridge')

#trace plot of bayesian approach
#par(mfrow=c(2,5))
#for(i in 1:10){plot(fit_b$beta_samples[,i],type = 'l')}
#plot similar figure 2






p <- 10

ypos <- 1:p
bayes_median <- gibbs_results[[best_ind_bayes]]$fit$beta_median
target_l1 <- sum(abs(bayes_median))
l1_norms <- colSums(abs(beta_lasso))
j <- which.min(abs(l1_norms-target_l1))

lasso_l1 <- beta_lasso[,j]
ci_low <- apply(gibbs_results[[best_ind_bayes]]$fit$beta_samples,2,quantile,probs=0.025)
ci_up <- apply(gibbs_results[[best_ind_bayes]]$fit$beta_samples,2,quantile,probs=0.975)

plot(NA,NA,xlim=range(c(ci_low,ci_up,beta_ols,beta_lasso,lasso_l1,bayes_median)),
     ylim=c(p+0.5,0.5),
     xlab='',
     ylab='',
     yaxt='n',
     bty='l'
     
     
     )
axis(2,at=ypos,labels=ypos,las=1)
abline(v=0,lty=3)
#bayesian credible intervals
arrows(ci_low,ypos,ci_up,ypos,angle=90,code=3,length=0.01)


points(bayes_median,ypos,pch=1,cex=1)
points(bayes_median,ypos,pch=3,cex=1)
points(beta_ols,ypos,pch=4,cex=1.8)
points(beta_lasso_best,ypos+0.18,pch=2,cex=1.2)
points(lasso_l1,ypos-0.18,pch=6,cex=1.2)
#obtain the 95% credible intervals of lambda
lambda_hat <- 0.237
fit_hat <- lasso_gibbs(X,y,lambda=lambda_hat,n_iter=11000,burn=1000)
tau2_samples <- fit_hat$tau2_samples
  
log_lr <- function(lambda, lambda_hat, tau2_samples){
  M <- nrow(tau2_samples)
  p <- ncol(tau2_samples)
  
  log_w <- numeric(M)
  
  for(m in 1:M){
    tau2_m <- tau2_samples[m, ]
    
    log_w[m] <- p * log(lambda^2 / lambda_hat^2) -
      (lambda^2 - lambda_hat^2) * sum(tau2_m) / 2
  }
  
  # log mean exp trick
  max_log_w <- max(log_w)
  log_mean_w <- max_log_w + log(mean(exp(log_w - max_log_w)))
  
  return(log_mean_w)
}

LR_lambda <- function(lambda){
  log_ratio <- log_lr(lambda, lambda_hat, tau2_samples)
  D <- -2 * log_ratio
  return(D)
}
lambda_grid <- seq(0.05, 0.6, length.out = 100)

LR_vals <- sapply(lambda_grid, LR_lambda)

cutoff <- qchisq(0.95, df = 1)

lambda_ci <- lambda_grid[LR_vals <= cutoff]

range(lambda_ci)
  
  
  
