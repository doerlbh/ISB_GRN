# Creating n-gene networks
# i) Create a n-gene network
# ii) Determine how many steady states it can have in different parameters
# iii) only keep those who have more than one steady states
# iv) output the network states
# Author: Baihan Lin
# Date:   August 2016

# Modified from:
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

set.seed(1)
Opath="./data_20160817/"

num = 5         # network number 
n = 4           # gene
trial = 50      # trial
trshd = 1       # threshold for different states
sstrshd = 1     # threshold for equilibrium of steady states
tend = 3000     # threshold for equilibrium vs. non-equlibrium
zthrs = 1    # threshold for non-trivial states

Nstate = mat.or.vec(num,1)  # store how many states each network can have

# Starting values
randset = 20  # range of initial states
Prandset = 200 # range of parameters

x0 = mat.or.vec(trial,n)   # startvalues (genes)
x0 = randset*matrix(round(runif(trial*n),randset), trial, n) 

temp1 = mat.or.vec(n,1)
temp2 = mat.or.vec(n,1)

# Direct parameter declaration:
x = mat.or.vec(n,1)   # gene number
A = mat.or.vec(n,n)   # gene sythesis rate
M = mat.or.vec(n,1)   # gene degradation rate
P = mat.or.vec(n,n)   # parameters, i.e. gene-gene interaction (repress or activate) 
N = mat.or.vec(n,n)   # interaction types

p = 1;
countcall = 1;
while (p <= num) {
  Nss = NULL      # store all the steady states 
  ssc = 1;
  
  P = Prandset*matrix(round(runif(n*n),Prandset), n, n) 
  M = matrix(round(runif(n*1), Prandset),n,1)
  A = 100*diag(P)
  
  for (node in 1:n) {
    N[node,] = sample(c(1, 0,-1),n,TRUE)
    N[node,node] = 0
    #print(N)
  }
  
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
    if (pht < tend) {
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
      if (mean(abs(xrt[t,,length(times)]-xrt[t,,length(times)-3])) < sstrshd) {
        if (mean(abs(xrt[t,,length(times)]-xrt[t,,length(times)-5])) < sstrshd) { 
          if(mean(abs(xrt[t,,length(times)]-xrt[t,,length(times)-7])) < sstrshd) { 
            #if(mean(abs(xrt[t,,length(times)]-xrt[t,,length(times)-9])) < sstrshd) { 
            equ = 1;
          } else {
            equ = 0;
          }
        } else {
          equ = 0;
        }
        #} else {
        # equ = 0;
        #}
      } else {
        equ = 0;
      }
      
      if (equ == 0) {
        ph = pht;
        xr = xrt; 
        times = timest;
        t = t - 1;
      }
      
      t = t + 1;
      #print(ph)
      end = 0;
    } else {
      end = 1;
    }
  }
  
  if (end == 0) {
    
    # plot gene time profiles                                                       
    #graphics.off()
    #windows(xpos=1,ypos=-50,width=n,height=4)
    
    Nss = cbind(Nss, xr[1,,length(xr)/(n*trial)]);
    for (i in 1:trial){
      for (j in 1:round(length(Nss)/node)) {
        if (mean(abs(xr[i,,length(xr)/(n*trial)] - Nss[,j]))>trshd) {
          if (mean(abs(xr[i,,length(xr)/(n*trial)]))>zthrs) {
            ssc = ssc + 1;
            Nss = cbind(Nss, xr[i,,length(xr)/(n*trial)]);        
          }
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
      
      sink(paste(Opath,"N",n,"-X",p,"-S",ssc,".txt",sep=""))
      
      cat(sprintf("Network %d with %d nodes has %d states:\n", p, n, ssc))
      cat("=============================\n")
      cat("Network parameters:\n")
      cat("\n P is \n")
      
      for(node in 1:n) {
        cat(P[node,]);
      }
      
      cat("\n M is \n")
      cat(M);
      
      cat("\n A is \n")
      cat(A);
      
      cat("\n N is \n")
      cat(N)
      
      cat("\n Network steady states: \n")    
      
      for(j in 1:ssc) {
        cat(Nss[,j]);
        cat("\n")
      }
      
      cat("\n=============================\n")
      sink()
      
      # For better reading
      sink(paste(Opath,"N",n,"-X",p,"-Para.txt",sep=""))
      for(node in 1:n) {
        cat(P[node,]);
        cat("\n")
      }
      cat(M);
      cat("\n")
      cat(N)
      sink()
      
      sink(paste(Opath,"N",n,"-X",p,"-SS.txt",sep=""))
      for(j in 1:ssc) {
        cat(Nss[,j]);
        cat("\n")
      }
      sink()
      
      Nstate[p] = ssc;
    }
    
    p = p + 1;
    cat("round ", countcall, "\n");
    countcall = countcall + 1;
  }
}

sink(paste(Opath,"SScount.txt",sep=""))
cat("ss = ", Nstate, "\n");
sink()
