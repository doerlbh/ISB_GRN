# Simulate nine-gene dynamics using ODE
# i) Calculate trajectories
# ii) Plot gene time profiles
# iii) Simulate different parameters input
# iv) Determine how many steady states it can have in different parameters
# v) Explore the bifurcation of the statess
# Author: Baihan Lin
# Date:   July 2016

# Modified from
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
for (i in 9:trial){
  x1t <- runif(1, 0, 1000)
  x2t <- runif(1, 0, 1000)
  x3t <- runif(1, 0, 1000)
  x0[i,] <- c(x1t,x2t,x3t) 
}

ph = 20
times <- seq(0,ph,0.1)  # time steps for output
parms <- c()          # parameter (if necesarry)

# Direct parameter declaration:
n = 9
X = mat.or.vec(n,1)   # gene number
A = mat.or.vec(n,1)   # gene synthesis rate
K = mat.or.vec(n,n)   # gene-gene interaction (repress or activate) 
M = mat.or.vec(n,1)   # gene degradation rate

# init values
A[1] = 
A[2] = 
A[3] = 
M = M + 1
K[2,1]=K[3,1]=K[3,2]=1
K[1,3]=100

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
#time = mat.or.vec(trial,length(times))
#x1r = mat.or.vec(trial,length(times))
#x2r = mat.or.vec(trial,length(times))
#x3r = mat.or.vec(trial,length(times))

#for (i in 1:trial){
#  res <- lsoda(x0[i,],times, func, parms)   # solve it
#  res <- as.data.frame(res)             # make a data frame
#  time[i,] <- res$time                      # time
#  x1r[i,] <- res$x1                        
#  x2r[i,] <- res$x2                         
#  x3r[i,] <- res$x3                   
#}

#res = data.frame(trial,length(times),length(times),length(times),length(times))
#res[1,,,,] = res0

#res = list()
#res$x00 =res0

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
graphics.off()

#maxY = max( c(max(X01),max(X02),max(X03)))
#maxY = max( c(max(X01),max(X02),max(X03),max(X11),max(X12),max(X13),
#              max(X21),max(X22),max(X23),max(X31),max(X32),max(X33),
#              max(X41),max(X42),max(X43),max(X51),max(X52),max(X53),
#              max(X61),max(X62),max(X63),max(X71),max(X72),max(X73)))

maxY = max( c(max(x1r[1,]),max(x2r[1,]),max(x3r[1,])))
#maxY = max( c(max(x1r),max(x2r),max(x3r)))

plot(times,x1r[1,],type="l",lwd=2,col=rgb(0,1,0), ylim=c(0, maxY) )
for (i in 1:trial){
  lines(times,x1r[i,],lwd=2,col=rgb(0,1,0))
  lines(times,x2r[i,],lwd=2,col=rgb(1,0,0))
  lines(times,x3r[i,],lwd=2,col=rgb(0,0,1))
}
