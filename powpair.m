% Baihan Lin, August 2016

function temp1 = powpair(y, (N>0).*P(node,:))
% ode function used to generate networks

    dydt = zeros(n,1); 
    parfor node = 1:n
      temp1 = powpair(y, (N>0).*P(node,:));
      temp1(temp1 == 0) = 1;
      temp2 = powpair(y, (N<0).*P(node,:));
      temp2(temp2 == 0) = 1;
      dydt(node) = A(node).*prod(temp1)/(prod(1+temp1).*prod(1+temp2)) - M(node).*y(node);
    end

end