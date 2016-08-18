% Baihan Lin, August 2016

function [y,nph,ss] = netrun(P,N,A,M,ph,eq,y0,tend)
% run networks

if ph > tend
    ss = 0;
else
    if ~eq
        ph = 2*ph;
    end
    
    tspan = 0:0.1:ph;
    
    [t,y] = ode45(@(t,y) netode(t,y,P,N,A,M), tspan, y0);
    
    c = size(y,2);
    neq = (abs(mean(y(:,c)-y(:,c-2))) < sstrshd);
    neq = neq*(abs(mean(y(:,c)-y(:,c-4))) < sstrshd);
    neq = neq*(abs(mean(y(:,c)-y(:,c-6))) < sstrshd);
    nph = ph;
    
    if ~neq
        [y,t,nph] = netrun(P,N,A,M,nph,1,y0);
    end
    ss = 1;
end
end