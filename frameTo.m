% Baihan Lin, August 2016

function y = frameTo(y,thrs, big)
% only keep big enough number

for r = 1:size(y, 1)
    y(r,:) = y(r,:).*(y(r,:)>thrs);
end

y(y == 0) = big;

end