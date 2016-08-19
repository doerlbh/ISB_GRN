% Baihan Lin, August 2016

function y = pertrun(pet,px,n,P,N,A,M,ph,eq,y0,sstrshd,zthrs,pow,Pt,Nt,At,Mt)
% perturb networks

if ~px
    
    petspan = 0:0.1:ph;
    
    [t,yst] = ode45(@(t,y) netodepownorm(t,y,Pt,Nt,At,Mt,n,pow), petspan, y0);
    %     [t,yst] = ode45(@(t,y) netodepow(t,y,Pt,Nt,At,Mt,n), petspan, y0);
    %     [t,yst] = ode45(@(t,y) netode(t,y,Pt,Nt,At,Mt,n), petspan, y0);
    
else
    yst = y0;
end

if ph > tend
    ss = 0;
else
    if ~eq
        ph = 2*ph;
    end
    
    tspan = 0:0.1:ph;
    
    [t,y] = ode45(@(t,y) netodepownorm(t,y,P,N,A,M,n,pow), tspan, yst);
    %     [t,y] = ode45(@(t,y) netodepow(t,y,P,N,A,M,n), tspan, yst);
    %     [t,y] = ode45(@(t,y) netode(t,y,P,N,A,M,n), tspan, yst);
    y = nontrivial(y,zthrs);
    c = size(y,1);
    neq = (abs(mean(y(c-2,:)-y(c,:))) < sstrshd);
    neq = neq*(abs(mean(y(c-10,:)-y(c,:))) < sstrshd);
    neq = neq*(abs(mean(y(c-20,:)-y(c,:))) < sstrshd);
    nph = ph;
    ss = 1;
    if ~neq
y = pertrun(pet,1,n,P,N,A,M,nph,1,yst,sstrshd,zthrs,pow,Pt,Nt,At,Mt)
    end
end
end