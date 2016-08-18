% Baihan Lin, August 2016

function y = nontrivial(y,zthrs)
% only keep nontrivial number

for r = 1:size(y, 1)
    y(r,:) = y(r,:).*(y(r,:)>zthrs);
end

end