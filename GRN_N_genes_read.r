# Perturbing n-gene networks
# i) read a n-gene network
# ii) generate time frame
# iii) only keep 
# Author: Baihan Lin
# Date:   August 2016

# Modified from:
# Creating n-gene networks
# i) Create a nine-gene network
# ii) Determine how many steady states it can have in different parameters
# iii) only keep those who have more than one steady states
# iv) output the network states
# Author: Baihan Lin
# Date:   August 2016

# Creating 9-gene networks
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

set.seed(1234)
Ipath="./data_20160807/"
Opath="./data_20160807/"

pert = 20    # time steps of perturbation
num = 2      # num of network
trail = 10   # num of change of parameters
Amprg = -1000:1000  # Amplication range
chg = 1      # perturb how many parameters 

n = 9        # gene

# for (i in 1:num) {
#SS = as.matrix(read.table(paste("N9-X", i,"-SS.txt",sep=""), header=F))
#PP = as.matrix(read.table(paste("N9-X", i,"-Para.txt",sep=""), header=F))

c = 1;
while (c < trail) { 
  Amp = sample(Amprg,1,TRUE);    # Amplication
  
  SS = as.matrix(read.table(paste("N9-X1-SS.txt",sep=""), header=F));
  PP = as.matrix(read.table(paste("N9-X1-Para.txt",sep=""), header=F));
  PPc = PP;
  n = length(PP[1,]);
  
  for (i in 1:chg) {
    seq = sample(1:n*(n+1),1,TRUE);
    PPc[sample(1:n*(n+1),1,TRUE)] = Amp*PPc[sample(1:n*(n+1),1,TRUE)]
  }
  
  # normal parameters
  P = PP[1:n, 1:n];
  A = 100*diag(P);
  M = PP[n+1, 1:n]
  N = PP[n+2, 1:n]
  
  # perturbation parameters
  Pc = PPc[1:n, 1:n];
  Ac = 100*diag(Pc);
  Mc = PPc[n+1, 1:n]
  Nc = PPc[n+2, 1:n]
  
  # ODE function
  func <- function(t,xx,p)
  {
    dx = mat.or.vec(n,1)
    #print(length(xx))
    x = xx  
    for (node in 1:n) {
      temp1 = (N>0)*P[node,]*x
      temp1[temp1 == 0] = 1
      temp2 = (N<0)*P[node,]*x
      temp2[temp2 == 0] = 1
      dx[node] = A[node]*prod(temp1)/(prod(1+temp1)*prod(1+temp2)) - M[node]*x[node]
    }
    
    list(dx)    # give the change rates to the solver
  }
  
  # perturbation solve ODE
  parms = c()          # parameter (if necesarry)
  timesp = seq(0,pert,0.1) 
  
  x0 = 
  
  x0r <- x0[t,]
  res <- lsoda(x0r,times, func, parms)   # solve it
  res <- as.data.frame(res)             # make a data frame
  for (node in 1:n) {
    xr[node,] <- res[,node+1]  
  }
  
  
  
  xr = array(rep(1, trial*n*length(times)), dim=c(trial,n,length(times)))
  
  
  # solve ODE
  ph = 100
  parms = c()          # parameter (if necesarry)
  times = seq(0,ph,0.1) 
  xr = array(rep(1, trial*n*length(times)), dim=c(trial,n,length(times)))
  
  t = 1;
  while (t <= trial) {
    #print("t")
    #print(t);
    x0r <- x0[t,]
    res <- lsoda(x0r,times, func, parms)   # solve it
    res <- as.data.frame(res)             # make a data frame
    for (node in 1:n) {
      xr[t,node,] <- res[,node+1]  
    }
    #print(dim(res))
    oldt = length(times);
    pht = 2*ph;
    timest = seq(0,pht,0.1);
    xrt = array(rep(1, trial*n*length(timest)), dim=c(trial,n,length(timest)));
    xrt[,,1:oldt] = xr[,,1:oldt];
    for (trialt in 1:t){
      for (tt in oldt:length(timest)) {
        for (node in 1:n) {
          xrt[trialt,node,tt]=xr[trialt,node,oldt];
        }
      }
    }
    if (mean(abs(xrt[t,,length(times)]-xrt[t,,length(times)-1])) > sstrshd) {
      ph = pht;
      xr = xrt; 
      times = timest;
      t = t - 1;
    }
    
    t = t + 1;
    #print(ph)
  }
  
  # plot gene time profiles                                                       
  #graphics.off()
  #windows(xpos=1,ypos=-50,width=n,height=4)
  
  Nss = cbind(Nss, xr[1,,length(xr)/(n*trial)]);
  for (i in 1:trial){
    for (j in 1:round(length(Nss)/node)) {
      if (mean(abs(xr[i,,length(xr)/(n*trial)] - Nss[,j]))>trshd) {
        ssc = ssc + 1;
        Nss = cbind(Nss, xr[i,,length(xr)/(n*trial)]);        
      }
    }
  }
  
  if (ssc == 1) {
    p = p - 1;
  } else {
    maxY = max(xr)
    #maxY = max( c(max(x1r),max(x2r),max(x3r)))
    
    png(filename = paste(Opath,"N", n, "-T", trial,"-P",p, ".png", sep=""), 
        width = 480, height = 480, 
        units = "px", pointsize = 12, bg = "white")
    
    plot(times,xr[1,1,],ylim=c(0, maxY), main=paste("N",n,"-T",trial,"-Trajectory",sep=""),
         type="l",xlab="t",ylab="x",lwd=2,col=rgb(0,0,1/n))
    #legend("topleft", lty=1:1)
    
    for (i in 1:trial){
      for (node in 1:n) {
        lines(times,xr[i,node,],lwd=2,col=rgb(0,0,node/n));
      }
    }
    
    mtext(paste(n,"_gene_",ssc, "_state_network_#", p, sep=""))
    
    dev.off()
    
    
  }
  c = c + 1;
}
#}



