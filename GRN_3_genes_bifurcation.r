# Bifurcation in three-gene dynamics using ODE
# iv) Determine how many steady states it can have in different parameters
# v) Explore the bifurcation of the statess
# Author: Baihan Lin
# Date:   July 2016

# Modified from:
# Simulate three-gene dynamics using ODE
# i) Calculate trajectories
# ii) Plot gene time profiles
# iii) Simulate different parameters input
# iv) Determine how many steady states it can have in different parameters
# v) Explore the bifurcation of the statess
# Author: Baihan Lin
# Date:   July 2016

# Modified from: 
# Model three-gene dynamics using ODE
# i) Calculate trajectories
# ii) Plot gene time profiles
# iii) how many steady states it can have
# iv) write a code to scan the parameters for different steady state.
# Author: Baihan Lin
# Date:   July 2016

# Modified from: 
# Model three-species food web dynamics using ODE
# i) Calculate trajectories
# ii) Plot species time profiles
# Author: Joseph X. Zhou
# Date:   12/11/2014
# Modify History:
# 1) Date: 06/04/2014


rm(list=ls())
require(deSolve) # load the ode package

n = 3
trial = 1000

x0 = mat.or.vec(trial,n)   # startvalues (genes)
x0[1,] <- c(x1=0,x2=0,x3=0)      
x0[2,] <- c(x1=10,x2=0,x3=0)  
x0[3,] <- c(x1=0,x2=10,x3=0) 
x0[4,] <- c(x1=0,x2=0,x3=10)     
x0[5,] <- c(x1=10,x2=10,x3=0)     
x0[6,] <- c(x1=10,x2=0,x3=10)     
x0[7,] <- c(x1=0,x2=10,x3=10)     
x0[8,] <- c(x1=10,x2=10,x3=10)  

set.seed(1234)
randset = 1000
for (i in 9:trial){
  x1t <- runif(1, 0, randset)
  x2t <- runif(1, 0, randset)
  x3t <- runif(1, 0, randset)
  x0[i,] <- c(x1t,x2t,x3t) 
}

ph = 20
times <- seq(0,ph,0.1)  # time steps for output
parms <- c()          # parameter (if necesarry)

# Direct parameter declaration:
n = 3
X = mat.or.vec(n,1)   # gene number
A = mat.or.vec(n,1)   # gene synthesis rate
K = mat.or.vec(n,n)   # gene-gene interaction (repress or activate) 
M = mat.or.vec(n,1)   # gene degradation rate

# Here we made a few assumptions:
# 1) The synthesis of gene 1 and gene 2 are the same.
# 2) In interactions, all activations are the same.
# 3) In interactions, all repressions are the same.
# 4) Degradation rate are the same for all 3 genes.

# init values

Ptrial = 10
P0 = mat.or.vec(Ptrial,5)   # startvalues (genes)
P0[1,] <- c(P1=2,P2=15,P3=1,P4=1,P5=100)      

set.seed(1457)
Prandset = 1000
for (i in 2:Ptrial){
  P1t <- runif(1, 0, Prandset)
  P2t <- runif(1, 0, Prandset)
  P3t <- runif(1, 0, Prandset)
  P4t <- runif(1, 0, Prandset)
  P5t <- runif(1, 0, Prandset)
  P0[i,] <- c(P1t,P2t,P3t,P4t,P5t) 
}

for (j in 1:Ptrial){
  A[1] = A[2] = P0[j,1]
  A[3] = P0[j,2]
  M = M + P0[j,3]
  K[2,1]=K[3,1]=K[3,2]=P0[j,4]
  K[1,3]=P0[j,5]
  
  # ODE function
  func <- function(t,xx,p)
  {
  x1 <- xx[1]   
  x2 <- xx[2]   
  x3 <- xx[3]   
  
  dx1 <-  A[1]/(1 + K[1,3]*x3) - M[1]*x1
  dx2 <-  A[2]*K[2,1]*x1/(1 + K[2,1]*x1) - M[2]*x2
  dx3 <-  A[3]*K[3,1]*x1*K[3,2]*x2 /((1 + K[3,1]*x1)*(1 + K[3,2]*x2)) - M[3]*x3
  
    list(c(dx1,dx2,dx3))    # give the change rates to the solver
  }
  
  # solve ODE
  
  x1r = mat.or.vec(trial,length(times))
  x2r = mat.or.vec(trial,length(times))
  x3r = mat.or.vec(trial,length(times))
  
  #res0 <- lsoda(x00,times, func, parms)   # solve it
  #res0 <- as.data.frame(res0)             # make a data frame
  
  for (i in 1:trial){
    x0r <- x0[i,]
    names(x0r) = c("x1", "x2", "x3")
    res <- lsoda(x0r,times, func, parms)   # solve it
    res <- as.data.frame(res)             # make a data frame
    x1r[i,] <- res$x1                        
    x2r[i,] <- res$x2                         
    x3r[i,] <- res$x3                   
  }
  
  
  # plot gene time profiles                                                       
  #graphics.off()
  #windows(xpos=1,ypos=-50,width=n,height=4)
  
  maxY = max( c(max(x1r[1,]),max(x2r[1,]),max(x3r[1,])))
  #maxY = max( c(max(x1r),max(x2r),max(x3r)))
  
  png(filename = paste("Bif_trial", trial,"_P", P0[j,1], "_", P0[j,2],
                       "_", P0[j,3], "_", P0[j,4], "_", P0[j,5], 
                       "_rand0_", randset, "_Prand0_", Prandset,".png", sep=""), 
      width = 480, height = 480, 
      units = "px", pointsize = 12, bg = "white")
  
  plot(times,x1r[1,],type="l",lwd=2,col=rgb(0,1,0), ylim=c(0, maxY) )
  for (i in 1:trial){
    lines(times,x1r[i,],lwd=2,col=rgb(0,1,0))
    lines(times,x2r[i,],lwd=2,col=rgb(1,0,0))
    lines(times,x3r[i,],lwd=2,col=rgb(0,0,1))
  }
  
  dev.off()

}

