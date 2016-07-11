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

x0 <- c(x1=0,x2=0,x3=0)       # startvalues    (genes)
x0_1 = x0
ph = 120
times <- seq(0,ph,0.1)  # time steps for output
parms <- c()          # parameter (if necesarry)

# Direct parameter declaration:
n = 3
X = mat.or.vec(n,1)   # gene number
A = mat.or.vec(n,1)   # gene synthesis rate
K = mat.or.vec(n,n)   # gene-gene interaction (repress or activate) 
M = mat.or.vec(n,1)   # gene degradation rate

# init values
A[1] = A[2] = 2
A[3] = 15
M = M + 1
K[2,1]=K[3,1]=K[3,2]=1
K[1,3]=100

# ODE function
func <- function(t,xx,p)
{
x1 <- xx[1]   # actual prey value
x2 <- xx[2]   # actual pred value
x3 <- xx[3]   # actual prey value


dx1 <-  A[1]/(1 + K[1,3]*x3) - M[1]*x1
dx2 <-  A[2]*K[2,1]*x1/(1 + K[2,1]*x1) - M[2]*x2
dx3 <-  A[3]*K[3,1]*x1*K[3,2]*x2 /((1 + K[3,1]*x1)*(1 + K[3,2]*x2)) - M[3]*x3

  list(c(dx1,dx2,dx3))    # give the change rates to the solver
}

# solve ODE
res <- lsoda(x0,times, func, parms)   # solve it
res <- as.data.frame(res)             # make a data frame

time <- res$time                      # time
X1 <- res$x1                         # prey
X2 <- res$x2                         # predator
X3 <- res$x3                         # prey

# plot gene time profiles                                                       
graphics.off()
windows(xpos=1,ypos=-50,width=n,height=4)

maxY = max( c(max(X1),max(X2),max(X3) ) )
         
plot(time,X1,type="l",lwd=2,col=rgb(0,1,0), ylim=c(0, maxY) )
lines(time,X2,lwd=2,col=rgb(1,0,0))
lines(time,X3,lwd=2,col=rgb(0,0,1))
#legend(30, 50, c("X1", "X2", "X3"))

windows(xpos=50,ypos=50,width=4,height=4)
plot(X1,X2,type="l",lwd=1,col=rgb(0,0,1),xlim=c(0,max(X1)),ylim=c(0,max(X2)))
#points(x0_1[1],x0_1[2])
points(x0[1],x0[2],col=rgb(0,1,0))
n = length(X1)
points(X1[n],X2[n],col=rgb(1,0,0))

windows(xpos=450,ypos=50,width=4,height=4)
plot(X1,X3,type="l",lwd=1,col=rgb(0,0,1),xlim=c(0,max(X1)),ylim=c(0,max(X3)))
#points(x0_1[1],x0_1[2])
points(x0[1],x0[3],col=rgb(0,1,0))
n = length(X1)
points(X1[n],X3[n],col=rgb(1,0,0))

windows(xpos=800,ypos=50,width=4,height=4)
plot(X2,X3,type="l",lwd=1,col=rgb(0,0,1),xlim=c(0,max(X2)),ylim=c(0,max(X3)))
#points(x0_1[1],x0_1[2])
points(x0[2],x0[3],col=rgb(0,1,0))
n = length(X1)
points(X2[n],X3[n],col=rgb(1,0,0))

