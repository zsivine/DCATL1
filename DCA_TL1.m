%--------------------------------------------------------------------------
% Difference of Convex Functions Methods for 1-D Compressed Sensing Problem. 
%
% Solves
%           min  P_a(x)
%           s.t. Ax = b
%
% Reference: "Minimization of Transformed L_1 Penalty: Theory, 
%             Difference of Convex Function Algorithm, 
%             and Robust Application in Compressed Sensing." 
%             Shuai Zhang, Jack Xin 
%             IEEE transactions on Information Theory, 2014
% Available at: 
%             http://arxiv.org/abs/1411.5735
%             ftp://ftp.math.ucla.edu/pub/camreport/cam14-86.pdf
% 
% Author: Shuai Zhang  
% Date: Feb 07. 2015
%--------------------------------------------------------------------------

% min_x  .5||Ax-b||^2 + gam(sum  (a+1)|x_i| / (a+|x_i|) )
% here g(x) =  .5||Ax-b||^2 + C||x||^2 + gam* (a+1)/a * ||x||_1  
%      h(x) = gam*{(a+1)/a *||x||_1 - sum (a+1)|x_i|/(a+|x_i|)} + C||x||^2
%      F(x) = g(x) - h(x)
%Input: dictionary A, data b, parameters option and initial value x0
% recommended value 
% x0: initial value, better to chosed as zeros vector, see reference
% pm.C       1.0e-9    
% pm.del     1.0e-4
% pm.maxoit  1000
% pm.gam     1.0e-5
% pm.tol     1.0e-5
% pm.a       1
%Output: computed coefficients x


function x = DCA_TL1(A,b,pm,x0)


if ~exist('pm','var'), 
    pm = []; 
end
C = 1.0e-9;    
% del = 1.0e-4;
maxoit = 20;
gam = 1.0e-6;
tol = 1.0e-5;
del = gam*10;
a = 1;
if isfield(pm,'C'),        C = pm.C;             end
if isfield(pm,'del'),      del = pm.del;         end
if isfield(pm,'maxoit'),   maxoit = pm.maxoit;   end
if isfield(pm,'gam'),      gam = pm.gam;         end
if isfield(pm,'tol'),      tol = pm.tol;         end
if isfield(pm,'a'),        a = pm.a;             end

[M,N] = size(A);
x = x0;
p_init = zeros(N,1);
%precompute inverse to be used in inner problem
L = chol( speye(M) + 1/(2*C+del)*(A*A'), 'lower' ); % cholesky factorization
L = sparse(L);
U = sparse(L');
Atb = A'*b;

for it = 1:maxoit
    if norm(x) < eps
        f = -Atb;
    else
        V = gam*(a+1)/a * sign(x) - gam*sign(x)*(a+1)./( a+abs(x) ) ...
            + gam*(a+1)*x./( a+abs(x) ).^2 + 2*C*x;
        f = -Atb - V;
    end
    [x_new,p_init] = DCA_TL1_sub(A,L,U,f,pm,x,p_init);
    %tolerance based stopping condition
    err = max(abs(x-x_new))/max(max(abs(x_new)),1);
%     err = sqrt(sum((x-x_new).^2))/max(sqrt(sum(x.^2)),1)
    x = x_new;
    if err < tol     
        disp(['tolerance met after ' num2str(it) ' iterations with a = ' num2str(a)]);
        break;
    end 
end

if it == maxoit, fprintf('Unconstrained DCATL1 fails ...... \n'); end

end