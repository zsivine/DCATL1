function test_3Nmatrix(N)

% test 3*N matrix, column vectors are generated uniformly in unit sphere
% matrix A uniformly over the sphere in R^3
% write the vectors in the angle variables, then sample the angles by uniform distribution
test_num = 50;
thresh = 1e-3;
succ = zeros(3,1);

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


for i = 1:test_num
    
    % intialize 
    theta = pi*(rand(1,N) -1/2);
    fa = 2*pi*rand(1,N);
    A = [cos(theta).*cos(fa);cos(theta).*sin(fa);sin(theta)];
    fprintf('The mutual coherence of matrix A is %d \n', mucohere_matrix(A));
    
    % separate two nonzeros elements with certain distance ex. N/3
    d = round(N/3); 
    supp = randsample_separated(N,2,d);
    x = zeros(N,1);
    x(supp) = sign(randn(2,1));
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

%     tic;
%     pm.p = 4;
%     z4 = l1_lP(A,y,pm,zeros(N,1));
%     time4 = toc;
%     tic;
%     pm.p = 6;
%     z5 = l1_lP(A,y,pm,zeros(N,1));
%     time5 = toc;
%     tic;
%     pm.p = 8;
%     z6 = l1_lP(A,y,pm,zeros(N,1));
%     time6 = toc;

    % display error 
    err1 = norm(z1-x)/norm(x)
    if err1 < thresh
        succ(1) = succ(1) + 1;
        fprintf('relative error for DCA_USC: %4.3e, iteration time: %4.3e \n',err1,time1);
    end

    err2 = norm(z2-x)/norm(x);
    if err2 < thresh
        succ(2) = succ(2) + 1;
        fprintf('relative error for IRucLq_v: %4.3e, iteration time: %4.3e \n',err2,time2);
    end

    err3 = norm(z3-x)/norm(x);
    if err3 < thresh
        succ(3) = succ(3) + 1;
        fprintf('relative error for l1 - l2: %4.3e, iteration time: %4.3e \n',err3,time3);
    end

%     err4 = norm(z4-x)/norm(x);
%     if err4 < thresh
%         succ(4) = succ(4) + 1;
%         fprintf('relative error for l1 - l4: %4.3e, iteration time: %4.3e \n',err4,time4);
%     end
%     err5 = norm(z5-x)/norm(x);
%     if err5 < thresh
%         succ(5) = succ(5) + 1;
%         fprintf('relative error for l1 - l6: %4.3e, iteration time: %4.3e \n',err5,time5);
%     end
%     err6 = norm(z6-x)/norm(x);
%     if err6 < thresh
%         succ(6) = succ(6) + 1;
%         fprintf('relative error for l1 - l8: %4.3e, iteration time: %4.3e \n',err6,time6);
%     end    
    
    subplot(2,1,1);
    plot(x,'o'); 
    subplot(2,1,2);
    plot(z1,'*');
    pause;
    
end
succ = succ./test_num;
succ








