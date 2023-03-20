getwd()
setwd('E:/Learn/MallowsMA_code_review/')
library(quadprog)
library(MASS)

rsquare = seq(.1,.9,.1)
alpha = 1
n = 400
m = 1000
Risk = matrix(0,6,length(rsquare))
for (i in 1:length(rsquare)) {
  theta = sqrt(1/(1-rsquare[i])-1)*sqrt(2*alpha)*(1:m)^(-alpha-.5)
  for (j in 1:100) { # repeat 100
    print(j)
    ## DGP
    X = mvrnorm(n, rep(1,m), diag(m))
    epsilon = X[,m]
    X = cbind(rep(1,n),X[,-m])
    mu = X %*% theta
    Y = mu + epsilon
    ## Model
    M = round(3*n^(1/3))
    Y_fit = matrix(NA,n,M)
    e_fit = matrix(NA,n,M)
    sigmas = numeric(M)
    for (k in 1:M) {
      newX = X[,1:k]
      beta = solve(t(newX)%*%newX)%*%t(newX)%*%Y
      Y_fit[,k] = newX %*% beta
      e_fit[,k] = Y - Y_fit[,k]
      sigmas[k] = (t(e_fit[,k]) %*% e_fit[,k])/(n-k)
    }
    sighat <- sigmas[M]
    
    a1 <- t(e_fit) %*% e_fit
    if (qr(a1)$rank<ncol(e_fit)) a1 <- a1 + diag(M)*1e-10
    a2 <- matrix(c(-sighat*(1:M)),M,1)  
    a3 <- t(rbind(matrix(1,nrow=1,ncol=M),diag(M),-diag(M)))
    a4 <- rbind(1,matrix(0,nrow=M,ncol=1),matrix(-1,nrow=M,ncol=1))
    
    w0 <- matrix(1,nrow=M,ncol=1)/M
    QP <- solve.QP(a1,a2,a3,a4,1)
    w <- QP$solution
    w <- as.matrix(w)
    w <- w*(w>0)
    w <- w/sum(w0)
    
    Yfore = matrix(NA,n,7)
    Yfore[,1] = Y_fit %*% w
    
    sigmaMLE = sigmas*(n-1:M)/n
    ## AIC
    AIC = n*log(sigmaMLE) + 2 * (1:M)
    i_aic = which(AIC==min(AIC),arr.ind = T)
    i_aic = i_aic[1]
    Yfore[,2] = Y_fit[,i_aic]
    
    ## sAIC
    w_aic = exp(-AIC/2)/sum(exp(-AIC/2))
    Yfore[,3] = Y_fit %*% w_aic
    
    ## BIC
    BIC = n*log(sigmaMLE) + log(n) * (1:M)
    i_bic = which(BIC==min(BIC),arr.ind = T)
    i_bic = i_bic[1]
    Yfore[,4] = Y_fit[,i_bic]
    
    ## sAIC
    w_bic = exp(-BIC/2)/sum(exp(-BIC/2))
    Yfore[,5] = Y_fit %*% w_bic
    
    ## SA
    Yfore[,6] = Y_fit %*% w0
    
    ## full model
    i_fm = which(sigmaMLE==min(sigmaMLE),arr.ind = T)
    i_fm = i_fm[1]
    Yfore[,7] = Y_fit[,i_fm]
    
    # compute risk
    risk = numeric(7)
    for (s in 1:7) {
      risk[s] = sum((Yfore[,s]-mu)^2)
    }
    for (s in 1:6) {
      Risk[s,i] = Risk[s,i] + risk[s]/risk[7]
    }
  }
}
Risk = Risk / 100
Risk
