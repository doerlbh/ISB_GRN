# Model three-gene dynamics using ODE
# i) Calculate trajectories
# ii) Plot species time profiles
# Author: Baihan Lin
# Date:   July 2016

# Modified from: 
# Model three-species food web dynamics using ODE
# i) Calculate trajectories
# ii) Plot species time profiles
# Author: Joseph X. Zhou
# Date:   12/11/2014
# Modify History:
# 1) Date: 06/04/2014,


rm(list=ls())
require(deSolve) # load the ode package

x0 <- c(x1=0.5,x2=1.1,x3=1.0)       # startvalues    (prey / predator )
x0_1 = x0
ph = 120
times <- seq(0,ph,0.1)  # time steps for output
parms <- c()          # parameter (if necesarry)

# I like more the direct parameter declaration:

n = 3
X = mat.or.vec(n,1)   # species number
A = mat.or.vec(n,1)   # species growth rate
B = mat.or.vec(n,n)   # prey death rate by predation 
C = mat.or.vec(n,1)   # growth saturation 
K = mat.or.vec(n,1)   # death saturation 
eps = mat.or.vec(n,n) # predator growth rate by predation
M = mat.or.vec(n,1)   # predator death rate
# init values
contl = 1
A <- A + 2*contl    # prey growth rate
B <- B + 2*contl     # prey death rate by predation
C <- C + 1.5
K <- K + 0.5
eps <- eps + 0.75 # predator growth rate by predation
M <- M + 0.5

#A[1] = 2.2
#A[2] = 1.9
#B[1,3] = 2
#B[2,3] = 2.5

func <- function(t,xx,p)
{
x1 <- xx[1]   # actual prey value
x2 <- xx[2]   # actual pred value
x3 <- xx[3]   # actual prey value

xs1 =  x1 + x2
dx1 <-  A[1]*x1*(C[1]-x1) - B[1,3]*x1*x3/(xs1 +K[1])  # calc prey growth
dx2 <-  A[2]*x2*(C[2]-x2) - B[2,3]*x2*x3/(xs1 +K[2]) 
dx3 <-  eps[1,3]* B[1,3]*x1*x3/(xs1 +K[1]) +  eps[2,3]*B[2,3]*x2*x3/(xs1 +K[2]) -   M[3]*x3  # calc pred growth

  list(c(dx1,dx2,dx3))    # give the change rates to the solver
}

res <- lsoda(x0,times, func, parms)   # solve it
res <- as.data.frame(res)             # make a data frame

time <- res$time                      # time
X1 <- res$x1                         # prey
X2 <- res$x2                         # predator
X3 <- res$x3                         # prey

                                                       
# from oscillation back to stable
#nn = length(X[1])
#contl = 0.5
#A <- 2#*contl    # prey growth rate
#B <- 2*contl    # prey death rate by predation
#x0 <- c(x=prey[nn],y=pred[nn])         # startvalues    (prey /predator )
#times <- seq(ph,100,0.1)  # time steps for output
#res <- lsoda(x0,times, func, parms)   # solve it
#res <- as.data.frame(res)             # make a data frame
#time1 <- res$time                      # time
#prey1 <- res$x                         # prey
#pred1 <- res$y                         # predator
#
#time = c(time,time1)
#prey = c(prey,prey1)
#pred = c(pred,pred1)

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
