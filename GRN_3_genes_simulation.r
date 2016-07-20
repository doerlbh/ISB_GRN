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
trial = 8

x0 = mat.or.vec(trial,n)   # startvalues (genes)
x0[1,] <- c(x1=0,x2=0,x3=0)      
x0[2,] <- c(x1=10,x2=0,x3=0)  
x0[3,] <- c(x1=0,x2=10,x3=0) 
x0[4,] <- c(x1=0,x2=0,x3=10)     
x0[5,] <- c(x1=10,x2=10,x3=0)     
x0[6,] <- c(x1=10,x2=0,x3=10)     
x0[7,] <- c(x1=0,x2=10,x3=10)     
x0[8,] <- c(x1=10,x2=10,x3=10)     

x00 <- c(x1=0,x2=0,x3=0)     
x01 <- c(x1=10,x2=0,x3=0)     
x02 <- c(x1=0,x2=10,x3=0)     
x03 <- c(x1=0,x2=0,x3=10)     
x04 <- c(x1=10,x2=10,x3=0)     
x05 <- c(x1=10,x2=0,x3=10)     
x06 <- c(x1=0,x2=10,x3=10)     
x07 <- c(x1=10,x2=10,x3=10)     

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

res = data.frame(trial,length(times),length(times),length(times),length(times))
res[1,,,,] = res0

res = list()
res$x00 =res0

x1r = mat.or.vec(trial,length(times))
x2r = mat.or.vec(trial,length(times))
x3r = mat.or.vec(trial,length(times))

res0 <- lsoda(x00,times, func, parms)   # solve it
res0 <- as.data.frame(res0)             # make a data frame

for (i in 1:trial){
  x0r <- x0[i,]
  names(x0r) = c("x1", "x2", "x3")
  
  tempres = res$x0r
  x1r[i,] <- tempres$x1                        
  x2r[i,] <- tempres$x2                         
  x3r[i,] <- tempres$x3                   
}

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

maxY = max( c(max(X01),max(X02),max(X03)))
#maxY = max( c(max(X01),max(X02),max(X03),max(X11),max(X12),max(X13),
#              max(X21),max(X22),max(X23),max(X31),max(X32),max(X33),
#              max(X41),max(X42),max(X43),max(X51),max(X52),max(X53),
#              max(X61),max(X62),max(X63),max(X71),max(X72),max(X73)))
         
plot(times,X01,type="l",lwd=2,col=rgb(0,1,0), ylim=c(0, maxY) )
lines(times,X02,lwd=2,col=rgb(1,0,0))
lines(times,X03,lwd=2,col=rgb(0,0,1))

lines(times,X11,lwd=2,col=rgb(0,1,0))
lines(times,X12,lwd=2,col=rgb(1,0,0))
lines(times,X13,lwd=2,col=rgb(0,0,1))

lines(times,X21,lwd=2,col=rgb(0,1,0))
lines(times,X22,lwd=2,col=rgb(1,0,0))
lines(times,X23,lwd=2,col=rgb(0,0,1))

lines(times,X31,lwd=2,col=rgb(0,1,0))
lines(times,X32,lwd=2,col=rgb(1,0,0))
lines(times,X33,lwd=2,col=rgb(0,0,1))

lines(times,X41,lwd=2,col=rgb(0,1,0))
lines(times,X42,lwd=2,col=rgb(1,0,0))
lines( timesX43,lwd=2,col=rgb(0,0,1))

lines( timesX51,lwd=2,col=rgb(0,1,0))
lines( timesX52,lwd=2,col=rgb(1,0,0))
lines( timesX53,lwd=2,col=rgb(0,0,1))

lines( timesX61,lwd=2,col=rgb(0,1,0))
lines( timesX62,lwd=2,col=rgb(1,0,0))
lines( timesX63,lwd=2,col=rgb(0,0,1))

lines( timesX71,lwd=2,col=rgb(0,1,0))
lines( timesX72,lwd=2,col=rgb(1,0,0))
lines( timesX73,lwd=2,col=rgb(0,0,1))

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

