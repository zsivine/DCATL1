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

function [x,p] = DCA_TL1_sub(A,L,U,f,pm,x_init,p_init)

if ~exist('pm','var'), 
    pm = []; 
end
C = 1.0e-9;    
gam = 1.0e-6;
del = gam*10;
a = 1;
% del = 1.0e-4;
if isfield(pm,'C'),        C = pm.C;             end
if isfield(pm,'del'),      del = pm.del;         end
if isfield(pm,'gam'),      gam = pm.gam;         end
if isfield(pm,'a'),        a = pm.a;             end
if isfield(pm,'maxit'),    maxit = pm.maxit;     end


%initialize variables for inner problem
p = p_init;
aa = (a + 1)/a;
z = x_init;
N = size(A,2);
%iterate to solve L1 regularized program
maxIter = min(round(2.5*N), maxit);
for i = 1:maxIter
    %x update
    q = -f + del*(z - p);   
    x = q/(2*C+del) - (A'*(U\(L\(A*q))))/(2*C+del)^2;
    
    %z update
    z = shrink(x + p,aa*gam/del);
    
    %p update
    p = p + x - z;
        
end

end