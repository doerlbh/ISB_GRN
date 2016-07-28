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
# Date:   Dec 2014

rm(list=ls())
require(deSolve) # load the ode package

set.seed(34)

num = 3      # network number 
n = 3        # gene
trial = 100  # trial
trshd = 0.001 # threshold for different states

# Starting values
randset = 1000  # range of 
Prandset = 100 # range of parameters

x0 = mat.or.vec(trial,n)   # startvalues (genes)
x0 = randset*matrix(round(runif(trial*n),randset), trial, n) 

# Time series
ph = 50
times <- seq(0,ph,0.1)  # time steps for output
parms <- c()          # parameter (if necesarry)
temp1 = mat.or.vec(n,1)
temp2 = mat.or.vec(n,1)

# Direct parameter declaration:
x = mat.or.vec(n,1)   # gene number
A = mat.or.vec(n,n)   # gene sythesis rate
M = mat.or.vec(n,1)   # gene degradation rate
P = mat.or.vec(n,n)   # parameters, i.e. gene-gene interaction (repress or activate) 
N = mat.or.vec(n,n)   # interaction types

for (p in 1:num) {
  P = Prandset*matrix(round(runif(n*n),Prandset), n, n) 
  M = 10*matrix(round(runif(n*1), Prandset),n,1)
  A = diag(P)
  
  for (node in 1:n) {
    N[node,] = sample(c(1, 0,-1),n,TRUE)
    N[node,node] = 0
    #print(N)
  }
  
  
  # ODE function
  func <- function(t,xx,p)
  {
    dx = mat.or.vec(n,1)
    #print(xx)
    x = xx  
    for (node in 1:n) {
      temp1 = (N>0)*P[node,]*x
      temp1[temp1 == 0] = 1
      temp2 = (N<0)*P[node,]*x
      temp2[temp2 == 0] = 1
      dx[node] = A[node]*prod(temp1)/(prod(1+temp1)*prod(temp2)) - M[node]*x[node]
    }
    
    list(dx)    # give the change rates to the solver
  }
  
  # solve ODE
  
  for (i in 1:trial){
    xr = mat.or.vec(n,length(times))
    x0r <- x0[i,]
    res <- lsoda(x0r,times, func, parms)   # solve it
    res <- as.data.frame(res)             # make a data frame
    for (node in 1:n) {
    xr[node,] <- res[,node+1]  
    }
  }
  
  # plot gene time profiles                                                       
  #graphics.off()
  #windows(xpos=1,ypos=-50,width=n,height=4)
  
  maxY = 1.5*max( c(max(x1r[,length(x1r)/trial]),max(x2r[,length(x2r)/trial]),
                max(x3r[,length(x1r)/trial])))
  #maxY = max( c(max(x1r),max(x2r),max(x3r)))
  
  png(filename = paste("./data/N", n, "-T", trial,"-P",p, ".png", sep=""), 
      width = 480, height = 480, 
      units = "px", pointsize = 12, bg = "white")
  
  plot(times,x1r[1,],ylim=c(0, maxY), main=paste("N",n,"-T",trial,"-Trajectory",sep=""),
       type="l",xlab="t",ylab="x",lwd=2,col=rgb(0,1,0))
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


