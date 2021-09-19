
# Function Definitions #


# Stationary Subspace Analysis of Nonstationary Processes 

# Sundararajan and Pourahmadi (2018), Journal of Time Series Analysis, 39 (3), 338-355 


######################################################################################


# Functions and Notations #

# Input 'x' : a T x p multivariate (zero-mean second-order nonstationary) time series

# T is the series length and p is the dimension

# m is the no. of DFT covariance lags to be used in the method; See Eqn (9).

# d.s is the dimension of the stationary subspace (should be less than input dimension p)


# Usage #

# stysol() returns the transpose of B_1 : (p x d.s) matrix that gives the stationary subspace Y_t^s = B_1 X_t


#########################################################################################

# Load the necessary library files
# Use install.packages("name-of-the-package") to install

library(expm)
library(rstiefel)
library(pracma)


########################################################################################

stysol=function(x,m,d.s){
  
  n=length(x[,1])
  p=length(x[1,])
  
  b.mats = list()
  nlist.bmats = 1
  obj.fn.vals = NULL
  
  
  objfn=function(B1){
    
    B1.n <- t(B1)
    
    n=length(x[,1])
    p=length(x[1,])
    
    y = B1.n%*%t(x)
    y=t(y)
    
    dft=mvfft(y)/sqrt(2*pi*n)
    dft=rbind(dft,dft)
    
    cr=array(0,c(d.s,d.s,m))
    for(r in 1:m)
    {
      sum1=matrix(0,d.s,d.s)
      
      for(k in 1:(n)){
        sum1=sum1+ t(t(dft[k,]))%*%Conj(dft[k+r,])/n  
      }
      
      cr[,,r]=sum1  
      
    } #end loop over lags
    
    sum2=0
    
    for(r in 1:m) sum2=sum2+sum( Re(cr[,,r])^2)+sum(Im(cr[,,r])^2 )
    
    return(sum2)
    
  } #end function objfn  
  
  # Gradient matrix
  grad_objfn=function(B1) {
    
    B1.m <- t(B1)
    
    n=length(x[,1])
    p=length(x[1,])
    
    dft=mvfft(x)/sqrt(2*pi*n)
    dft=rbind(dft,dft)
    
    cr=array(0,c(p,p,m))
    
    for(r in 1:m)
    {
      sum1=matrix(0,p,p)
      for(k in 1:(n)){
        sum1=sum1+ t(t(dft[k,]))%*%Conj(dft[k+r,])/n  
      }
      
      cr[,,r]=sum1  
    } # end loop over lags
    
    G1=matrix(0,d.s,p)
    G2=matrix(0,d.s,p)
    
    for(r in 1:m){   #sum2=sum2+cr[,,r]  
      
      A1 = Re(cr[,,r])
      A2 = Im(cr[,,r])
      
      G1 = G1 + B1.m%*%t(A1)%*%t(B1.m)%*%B1.m%*%A1 + B1.m%*%A1%*%t(B1.m)%*%B1.m%*%t(A1) +
        B1.m%*%A1%*%t(B1.m)%*%B1.m%*%t(A1) + B1.m%*%t(A1)%*%t(B1.m)%*%B1.m%*%A1
      
      G2 = G2 + B1.m%*%t(A2)%*%t(B1.m)%*%B1.m%*%A2 + B1.m%*%A2%*%t(B1.m)%*%B1.m%*%t(A2) +
        B1.m%*%A2%*%t(B1.m)%*%B1.m%*%t(A2) + B1.m%*%t(A2)%*%t(B1.m)%*%B1.m%*%A2
      
    } # end loop over r
    
    return(t(G1)+t(G2))  
    
  } # function to evaluate gradient
  
  ############################################################################
  
  ### (Optional) ###
  # Statrting B_0 for the gradient descent algorithm
  mat.start = matrix(0,d.s,d.s)
  sum.mat.tmp = matrix(0,p,p)
  cr.n = array(0,c(p,p,m))
  
  dft=mvfft(x)/sqrt(2*pi*n)
  dft=rbind(dft,dft)
  
  for(r in 1:m)
  {
    sum1=matrix(0,p,p)
    for(k in 1:(n)){
      sum1=sum1+ t(t(dft[k,]))%*%Conj(dft[k+r,])/n  
    }
    
    cr.n[,,r]=sum1  
  } # end loop over lags
  
  for(r in 1:m){
    eigen.r = eigen( Re(cr.n[,,r]) ) 
    eigen.i = eigen( Im( cr.n[,,r]) ) 
    
    sum.mat.tmp = sum.mat.tmp + eigen.r$vectors%*%diag(c(abs(eigen.r$values)))%*%t(eigen.r$vectors)
    + eigen.i$vectors%*%diag(c(abs(eigen.i$values)))%*%t(eigen.i$vectors)  
  }
  mat.start = eigen(sum.mat.tmp)$vectors[,1:d.s] 
  if(d.s==1)   mat.start = matrix(mat.start)
  
  ## End optional part for statrting B_0 for the gradient descent algorithm
  
  #################################################################################
  
  for(nbmat in 1:50){
    
    stiefel.m.rand = rustiefel(p , d.s) 
    
    if(nbmat==1) stiefel.m.rand = mat.start # (Optional)
    
    b.mats[[nlist.bmats]] = optStiefel(objfn, grad_objfn, Vinit = stiefel.m.rand,
                                       method="curvilinear",
                                       searchParams=list(rho1=0.1, rho2=0.9, tau=1),tol=1e-4,
                                       maxLineSearchIters=20)
    
    obj.fn.vals = c(obj.fn.vals,objfn(b.mats[[nlist.bmats]]))
    
    nlist.bmats = nlist.bmats + 1
    
  } # end loop over nbmats
  
  
  return(b.mats[[which.min(obj.fn.vals)]])
  
} # end function stysol


#########################################################################################
#########################################################################################
##################################### END ###############################################
