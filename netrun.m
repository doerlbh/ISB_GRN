% Baihan Lin, August 2016

function [y,nph,ss] = netrun(n,P,N,A,M,ph,eq,y0,tend,sstrshd,zthrs)
% run networks

if ph > tend
    ss = 0;
else
    if ~eq
        ph = 2*ph;
    end
    
    tspan = 0:0.1:ph;
    
    [t,y] = ode45(@(t,y) netode(t,y,P,N,A,M,n), tspan, y0);
    y = nontrivial(y,zthrs);
    c = size(y,1);
    neq = (abs(mean(y(c-2,:)-y(c,:))) < sstrshd);
    neq = neq*(abs(mean(y(c-10,:)-y(c,:))) < sstrshd);
    neq = neq*(abs(mean(y(c-20,:)-y(c,:))) < sstrshd);
    nph = ph;
    ss = 1;
    if ~neq
        [y,nph,ss] = netrun(n,P,N,A,M,nph,1,y0,tend,sstrshd,zthrs);
    end
end
end