clc, close all, clear;
M = 64; N = 1024;
F = 2;
K = 10;

% % linear system
% A = zeros(M,N);
% r = rand(M,1);
% l = 1:N;
% for k = 1:M
%     A(k,:) = sqrt(1/M)*cos(2*pi*r(k)*(l-1)/F);
% end
% 
% supp = randsample_separated(N,K,2*F);
% x = zeros(N,1);
% xs = randn(K,1);
% x(supp) = xs;
% b = A*x;

r = 0.0;
sigma = r*ones(N);
sigma = sigma + (1-r)*eye(N);
mu = zeros(1,N);
A = mvnrnd(mu,sigma,M);
A = unitnorm(A);
xz = randn(K,1);
xzs = [xz; zeros(N - size(xz,1),1) ];
ind = randperm(N);
x = xzs(ind);
b = A*x;

% % normalizing ... 
% d = 1./sqrt(sum(A.^2,2));
% A = sparse(1:M,1:M,d)*A;
% y = d.*y;


% load data.mat
% N = size(A,2);
% parameters  
pm.gam = 1e-8;
pm.C = 0;
pm.ep = 1e-6;
pm.maxit = 5000;
pm.maxoit = 20;
pm.otol = 1e-8;
pm.tol = 1e-5;
pm.del = pm.gam*10;
pm.a = 1;


%opts1.x0 = randn(N,1);
% opts1.so = round(K);
% opts1.maxit = 5000;
% opts1.tol = 1e-06;
% opts1.gamma = 0.5; opts1.eps0 = 1;
% lam1 = 1e-6; 


% tic;
% x1 = Cons_DCATL1(A,b,[],zeros(N,1));
% time1 = toc;
% tic;
% x2 = Cons_DCATL1_Linprog(A,b,[],zeros(N,1));
% time2 = toc;
% tic;
% x3 = Cons_DCATL1_Linprog_v2(A,b,[],zeros(N,1));
% time3 = toc;


tic;
x1 = DCA_TL1(A,b,pm,zeros(N,1));
time1 = toc;
tic;


time1 

relerr1 = norm(x1-x)/norm(x)
subplot(2,1,1);
plot(x);
subplot(2,1,2);
plot(x1);





