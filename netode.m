% Baihan Lin, August 2016

function dydt = netode(t,y, P, N, A, M)
% ode function used to generate networks

    dydt = zeros(n,1); 
    parfor node = 1:n
      temp1 = (N>0).*P(node,:).*y;
      temp1(temp1 == 0) = 1;
      temp2 = (N<0).*P(node,:).*y;
      temp2(temp2 == 0) = 1;
      dydt(node) = A(node).*prod(temp1)/(prod(1+temp1).*prod(1+temp2)) - M(node)*y(node);
    end

end