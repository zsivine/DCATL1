function y = unitnorm(y)
m = size(y,2);
for i = 1:m
    y(:,i) = y(:,i)./norm( y(:,i));
end