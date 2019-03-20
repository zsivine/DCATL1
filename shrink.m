function z = shrink(x,r)
z = max(0,x - r) - max(0,-x - r);
end