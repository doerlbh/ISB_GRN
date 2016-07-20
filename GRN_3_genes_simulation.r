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
X = mat.or.vec(n,)   # gene number

x00 <- c(x1=0,x2=0,x3=0)       # startvalues    (genes)
x01 <- c(x1=10,x2=0,x3=0)       # startvalues    (genes)
x02 <- c(x1=0,x2=10,x3=0)       # startvalues    (genes)
x03 <- c(x1=0,x2=0,x3=10)       # startvalues    (genes)
x04 <- c(x1=10,x2=10,x3=0)       # startvalues    (genes)
x05 <- c(x1=10,x2=0,x3=10)       # startvalues    (genes)
x06 <- c(x1=0,x2=10,x3=10)       # startvalues    (genes)
x07 <- c(x1=10,x2=10,x3=10)       # startvalues    (genes)

#x0_1 = x0
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

# init values
A[1] = A[2] = 2
A[3] = 15
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
res0 <- lsoda(x00,times, func, parms)   # solve it
res0 <- as.data.frame(res0)             # make a data frame
time0 <- res0$time                      # time
X01 <- res0$x1                        
X02 <- res0$x2                         
X03 <- res0$x3                   

res1 <- lsoda(x01,times, func, parms)   # solve it
res1 <- as.data.frame(res1)             # make a data frame
time1 <- res1$time                      # time
X11 <- res1$x1                        
X12 <- res1$x2                
X13 <- res1$x3                        

res2 <- lsoda(x02,times, func, parms)   # solve it
res2 <- as.data.frame(res2)             # make a data frame
time2 <- res2$time                      # time
X21 <- res2$x1                       
X22 <- res2$x2                  
X23 <- res2$x3                        

res3 <- lsoda(x03,times, func, parms)   # solve it
res3 <- as.data.frame(res3)             # make a data frame
time3 <- res3$time                      # time
X31 <- res3$x1                    
X32 <- res3$x2                        
X33 <- res3$x3                       

res4 <- lsoda(x04,times, func, parms)   # solve it
res4 <- as.data.frame(res4)             # make a data frame
time4 <- res4$time                      # time
X41 <- res4$x1                    
X42 <- res4$x2                        
X43 <- res4$x3                       

res5 <- lsoda(x05,times, func, parms)   # solve it
res5 <- as.data.frame(res5)             # make a data frame
time5 <- res5$time                      # time
X51 <- res5$x1                    
X52 <- res5$x2                        
X53 <- res5$x3                       

res6 <- lsoda(x06, times, func, parms)   # solve it
res6 <- as.data.frame(res6)             # make a data frame
time6 <- res6$time                      # time
X61 <- res6$x1                    
X62 <- res6$x2                        
X63 <- res6$x3                       

res7 <- lsoda(x07,times, func, parms)   # solve it
res7 <- as.data.frame(res7)             # make a data frame
time7 <- res7$time                      # time
X71 <- res7$x1                    
X72 <- res7$x2                        
X73 <- res7$x3                       


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

