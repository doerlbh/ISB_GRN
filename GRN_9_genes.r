# Creating n-gene networks
# i) Create a nine-gene network
# ii) Determine how many steady states it can have in different parameters
# iii) only keep those who have more than one steady states
# Author: Baihan Lin
# Date:   July 2016

# Wild Bifurcation without huge assumptions in three-gene dynamics using ODE
# iv) Determine how many steady states it can have in different parameters
# v) Explore the bifurcation of the statess
# Author: Baihan Lin
# Date:   July 2016

# Modified from:
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

set.seed(5784)

num = 1      # network number 
n = 3        # gene
trial = 100  # trial
trshd = 0.01 # threshold for different states

# Starting values
randset = 1000  # range of 
Prandset = 1000 # range of parameters

x0 = mat.or.vec(trial,n)   # startvalues (genes)
x0 = randset*matrix(round(runif(trial*n),randset), trial, n) 

# Time series
ph = 20
times <- seq(0,ph,0.1)  # time steps for output
parms <- c()          # parameter (if necesarry)

# Direct parameter declaration:
X = mat.or.vec(n,1)   # gene number
A = mat.or.vec(n,n)   # gene sythesis rate
K = mat.or.vec(n,n)   # gene-gene interaction (repress or activate) 
M = mat.or.vec(n,1)   # gene degradation rate

# init parameters
P0 = mat.or.vec(n,n)  # startvalues (genes)
N = mat.or.vec(n,n)   # interaction types

for (p in 1:num) {
  P0 = Prandset*matrix(round(runif(n*n),Prandset), n, n) 
  A = diag(P0);
  
  for (node in 1:n) {
    N[node,] = sample(c(1, 0,-1),n,TRUE)
    N[node,node] = 0
  }
  
  
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
  
  png(filename = paste("./data/N", n, "-T", trial,"-P",p, ".png", sep=""), 
      width = 480, height = 480, 
      units = "px", pointsize = 12, bg = "white")
  
  plot(times,x1r[1,],main=paste("N",n,"-T",trial,"-Trajectory",sep=""),
       type="l",xlab="t",ylab="x",lwd=2,col=rgb(0,1,0), ylim=c(0, maxY) )
  legend("topleft", legend=c("x1", "x2","x3"),col=c("green","red", "blue"), lty=1:1)
  states = 1;
  for (i in 1:trial){
    lines(times,x1r[i,],lwd=2,col=rgb(0,1,0))
    lines(times,x2r[i,],lwd=2,col=rgb(1,0,0))
    lines(times,x3r[i,],lwd=2,col=rgb(0,0,1))
    add = (abs(x1r[i,length(x1r)/trial] - x1r[1,length(x1r)/trial])<trshd)
    add = add*(abs(x2r[i,length(x1r)/trial] - x2r[1,length(x1r)/trial])<trshd)
    add = add*(abs(x3r[i,length(x1r)/trial] - x3r[1,length(x1r)/trial])<trshd)
    states = states + !add
  }
  mtext(paste(n,"-gene ",states, "-state network #", p, sep=""))
  dev.off()
}
}

