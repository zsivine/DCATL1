function y = frac_norm(x,p)
% compute fraction norm p for vector x
x = abs(x);
x = x.^p;
y = sum(x);
y = y^(1/p);