% Baihan Lin, August 2016

function dydt = netodepownorm(t, y, P, N, A, M, n, pow)
% ode function used to generate networks

% P = frameTo(P, 2, 2);

dydt = zeros(n,1);
parfor node = 1:n
    temp1 = ((N>0).*y).^pow;
    temp1k = ((N>0).*P(node,:)).^pow;
    temp1c = temp1 + temp1k;
    temp1 = nontrivial(temp1,1);
    temp1(temp1 == 0) = 1;
    temp2 = ((N<0).*y).^pow;
    temp2k = ((N<0).*P(node,:)).^pow;
    temp2c = temp2 + temp2k;
    temp2c = nontrivial(temp2c,1);
    temp2c(temp1 == 0) = 1;
    dydt(node) = A(node).*prod(temp1)/(prod(temp1c).*prod(temp2c)) - M(node).*y(node);
end

end