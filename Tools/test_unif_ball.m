function [succ,mu] = test_unif_ball(M,N,k)

% test M*N matrix, column vectors are generated uniformly in unit M-d sphere
% matrix A uniformly over the sphere in R^M
% k = sparsity 

test_num = 100;
thresh = 1e-3;
succ = zeros(4,1);
mu = zeros(test_num,1);

% coefficient for l1-l2
pm.gam = 1e-6;
pm.C = 1e-9;
pm.maxoit = 20;
pm.tol = 1e-5;
pm.del = pm.gam*10;
pm.a = 1;
% coefficient for IRucLq_v
q = 0.5;
lam2 = 1e-6;

% coefficient for yall1 
opyall.tol = 1e-5;
opyall.rho = 1e-5;
opyall.nonorth = 1;

for i = 1:test_num
    
    % intialize 
    A = unif_ball(M,N);
    mu(i) = mucohere_matrix(A); 
    fprintf('The mutual coherence of matrix A is %d \n', mu(i));
    
    % separate two nonzeros elements with certain distance ex. N/3
    d = round(N/(2*k)); 
    supp = randsample_separated(N,k,d);
    x = zeros(N,1);
    x(supp) = sign(randn(k,1));
    y = A*x;
    
    % solve by different methods
    tic;
    z1 = DCA_USC(A,y,pm,zeros(N,1)); 
    time1 = toc;
    tic;
    [z2,~] = IRucLq_v(A,y,lam2,q);
    time2 = toc;
    tic;
    z3 = DCAunl1_l2(A,y,pm,zeros(N,1));
    time3 = toc;
    tic;
    z4 = yall1(A,y,opyall);
    time4 = toc;



    % display error 
    err1 = norm(z1-x)/norm(x);
    if err1 < thresh
        succ(1) = succ(1) + 1;
        fprintf('relative error for DCATL1: %4.3e, iteration time: %4.3e \n',err1,time1);
    end

    err2 = norm(z2-x)/norm(x);
    if err2 < thresh
        succ(2) = succ(2) + 1;
        fprintf('relative error for IRucLq_v: %4.3e, iteration time: %4.3e \n',err2,time2);
    end

    err3 = norm(z3-x)/norm(x);
    if err3 < thresh
        succ(3) = succ(3) + 1;
        fprintf('relative error for DCA l1-l2: %4.3e, iteration time: %4.3e \n',err3,time3);
    end

    err4 = norm(z4-x)/norm(x);
    if err4 < thresh
        succ(4) = succ(4) + 1;
        fprintf('relative error for yall1: %4.3e, iteration time: %4.3e \n',err4,time4);
    end
    
%     % sparsity 
%     op_spar.eps = 5e-2;
%     s1 = norm0(z1,op_spar);
%     fprintf('sparsity difference: %d \n',s1-k );
%     s2 = norm0(z2,op_spar);
%     fprintf('sparsity difference: %d \n',abs(s2-k));
%     s3 = norm0(z3,op_spar);
%     fprintf('sparsity difference: %d \n',abs(s3-k));
    
    
%     subplot(2,1,1);
%     plot(x,'o'); 
%     subplot(2,1,2);
%     plot(z1,'*');
%     pause;
    
end

succ = succ./test_num;
mu = sum(mu(:))/test_num;


