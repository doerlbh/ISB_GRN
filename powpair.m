% Baihan Lin, August 2016

function x = powpair(y, p)
% power function by pair

n = length(y);
x = ones(1,n);
parfor c = 1:n
    x(c) = y(c)^p(c);
end
    
end