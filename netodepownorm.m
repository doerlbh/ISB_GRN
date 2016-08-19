% Baihan Lin, August 2016

function dydt = netodepownorm(t, y, P, N, A, M, n, pow)
% ode function used to generate networks

P = frameTo(P, 2, 2);

dydt = zeros(n,1);
parfor node = 1:n
    temp1 = ((N>0).*y).^pow;
    temp1k = ((N>0).*P(node,:)).^pow;
    temp2 = ((N<0).*y).^pow;
    temp2k = ((N<0).*P(node,:)).^pow;
    dydt(node) = A(node).*prod(temp1)/(prod(temp1k+temp1).*prod(temp2k+temp2)) - M(node).*y(node);
end

end