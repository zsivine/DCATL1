function z = DCATL1_RL_sparsity(l)
% Over-sampled DCT matrix 
% s: sparsity
% F: RL
% l: interger, such that separation length = l*F 

N = 1024 ; % # of columns
M = 64   ; % # of rows
thresh = 1e-3; % standard for success 
testnum = 10;

F = 10;

k_seq = 5:3:20;
len = length(k_seq);


succ = zeros(1,len);


% coefficient for l1-l2
pm.gam = 1e-6;
pm.C = 1e-9;
pm.ep = 1e-6;
pm.maxit = 5000;
pm.maxoit = 20;
pm.otol = 1e-8;
pm.tol = 1e-5;
pm.del = pm.gam*10;

pm.a = 1;


for ii = 1:len
    ii
    k = k_seq(ii);
    
    for j = 1:testnum
         
        %% generate matrix A and right hand y         
        % cosine matrix
        A = coherentdic(M,N,F);
        supp = randsample_separated(N,k,l*F);
        x = zeros(N,1);
        xx = randn(k,1);
        x(supp) = xx;
        y = A*x;

        % normalizing ... 
        d = 1./sqrt(sum(A.^2,2));
        A = sparse(1:M,1:M,d)*A;
        y = d.*y;            
        
        
        
        %% computing ...
        tic;
        z1 = DCA_USC(A,y,pm,zeros(N,1));
        time1 = toc;

        err1 = norm(z1-x)/norm(x);
        if err1 < thresh
            succ(1,ii) = succ(1,ii) + 1;
            fprintf('relative error for DCATL1: %4.3e, iteration time: %4.3e \n',err1,time1);
        end

    end
end
z = succ

