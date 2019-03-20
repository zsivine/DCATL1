% l0 norm function
function z = norm0(x,option)
if ~exist('option','var'), 
    option = []; 
end
eps = 0;
if isfield(option,'eps'),        eps = option.eps;             end

x = abs(x);
x = x > eps;
z = sum(x);
