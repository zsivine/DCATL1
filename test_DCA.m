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
%             Mathematical Programming
% Available at: 
%             http://arxiv.org/abs/1411.5735
%             ftp://ftp.math.ucla.edu/pub/camreport/cam14-86.pdf
% 
% Author: Shuai Zhang  
% Date: Feb 07. 2015
%--------------------------------------------------------------------------

function test_DCA(fr,kind,M,N,k_seq)
% kind: 1. Cosine matrix  
%       2. Gaussian matrix 

thresh = 1e-3; % standard for success 
testnum = 10;

if kind == 1, 
    F = fr;
elseif kind == 2
    r = fr;
end

% k_seq = 5:2:36;
len = length(k_seq);

% % % %  
% how many methods you want to use? 
succ = zeros(1,len);
% % % % 

% coefficient for l1-l2
pm.gam = 1e-8;
pm.C = 0;
pm.ep = 1e-6;
pm.maxit = 5000;
pm.maxoit = 20;
pm.otol = 1e-8;
pm.tol = 1e-5;
pm.del = pm.gam*10;

pm.a = 1;


% coefficient for IRucLq_v
q = 0.5;
lam2 = 1e-8;
% opts2.s0 = round(M/2);
% opts2.gamma = 0.5;

% coefficient for yall1 
opyall.tol = 1e-5;
opyall.rho = 1e-5;
opyall.nonorth = 1;

% coefficient for CEL0
pm_cel0.lamb = 1e-5;


for ii = 1:len
    [kind, fr, ii]
    k = k_seq(ii);
    
    for j = 1:testnum
         
        %% generate matrix A and right hand y 
        if kind == 1
            % cosine matrix
            A = coherentdic(M,N,F);
            supp = randsample_separated(N,k,2*F);
            x = zeros(N,1);
            xx = randn(k,1);
            x(supp) = xx;
            y = A*x;
            
            % normalizing ... 
            d = 1./sqrt(sum(A.^2,2));
            A = sparse(1:M,1:M,d)*A;
            y = d.*y;
            
        elseif kind == 2        
            % Gaussian matrix
            sigma = r*ones(N);
            sigma = sigma + (1-r)*eye(N);

            mu = zeros(1,N);
            A = mvnrnd(mu,sigma,M);
            A = unitnorm(A);

            xz = randn(k,1);
            xzs = [xz; zeros(N - size(xz,1),1) ];
            ind = randperm(N);
            x = xzs(ind);
            y = A*x;
            
            % normalizing ... 
            d = 1./sqrt(sum(A.^2,2));
            A = sparse(1:M,1:M,d)*A;
            y = d.*y;
            
        end
        
        
        
        %% computing ...
        tic;
        z1 = DCA_TL1(A,y,pm,zeros(N,1));
        time1 = toc;
%         tic;
%         [z2,~] = IRucLq_v(A,y,lam2,q);
%         time2 = toc;
%         tic;
%         z3 = DCAunl1_l2(A,y,pm,zeros(N,1));
%         time3 = toc;
% %         tic;
% %         z4 = yall1(A,y,opyall);
% %         time4 = toc;
%         tic; 
%         [z5,~] = macro_cel0(A,y,zeros(N,1),pm_cel0);
%         time5 = toc;
        
        

        err1 = norm(z1-x)/norm(x);
        if err1 < thresh
            succ(1,ii) = succ(1,ii) + 1;
            fprintf('relative error for DCATL1: %4.3e, iteration time: %4.3e \n',err1,time1);
        end

%         err2 = norm(z2-x)/norm(x);
%         if err2 < thresh
%             succ(2,ii) = succ(2,ii) + 1;
%             fprintf('relative error for IRucLq_v: %4.3e, iteration time: %4.3e \n',err2,time2);
%         end
% 
%         err3 = norm(z3-x)/norm(x);
%         if err3 < thresh
%             succ(3,ii) = succ(3,ii) + 1;
%             fprintf('relative error for l1 - l2: %4.3e, iteration time: %4.3e \n',err3,time3);
%         end

%         err4 = norm(z4-x)/norm(x);
%         if err4 < thresh
%             succ(4,ii) = succ(4,ii) + 1;
%             fprintf('relative error for yall1: %4.3e, iteration time: %4.3e \n',err4,time4);
%         end
                
%         err5 = norm(z5-x)/norm(x);
%         if err5 < thresh
%             succ(4,ii) = succ(4,ii) + 1;
%             fprintf('relative error for CEL0: %4.3e, iteration time: %4.3e \n',err5,time5);
%         end
        
%         tic;
%         pm.p = 4;
%         z4 = l1_lP(A,y,pm,zeros(N,1));
%         time4 = toc;
%         tic;
%         pm.p = 6;
%         z5 = l1_lP(A,y,pm,zeros(N,1));
%         time5 = toc;
%         tic;
%         pm.p = 8;
%         z6 = l1_lP(A,y,pm,zeros(N,1));
%         time6 = toc;
%         err4 = norm(z4-x)/norm(x);
%         if err4 < thresh
%             succ(4,i) = succ(4,i) + 1;
%             fprintf('relative error for l1 - l4: %4.3e, iteration time: %4.3e \n',err4,time4);
%         end
%         err5 = norm(z5-x)/norm(x);
%         if err5 < thresh
%             succ(5,i) = succ(5,i) + 1;
%             fprintf('relative error for l1 - l6: %4.3e, iteration time: %4.3e \n',err5,time5);
%         end
%         err6 = norm(z6-x)/norm(x);
%         if err6 < thresh
%             succ(6,i) = succ(6,i) + 1;
%             fprintf('relative error for l1 - l8: %4.3e, iteration time: %4.3e \n',err6,time6);
%         end

        fprintf('\n');
    end
end

succ = succ./testnum;
succ

%% graph
% figure
% x = 1:len;
% 
% plot(x,succ(1,:),'-ro',  x,succ(2,:),'-*b',  x,succ(3,:),'-g',...
%      x,succ(4,:),'-.k', 'linewidth',3);
% axis([0 , len+1, 0, 1]);
% 
% hleg = legend('DCATL1','IRucLq_v','DCA l1-l2','CEL0');
% 
% set(hleg,'Location','SouthWest')
% set(hleg,'Interpreter','none')
% set(gca,'XTick',x);
% 
% if kind == 1
%     k_label = { '5';   '7';  '9'; '11'; '13'; '15';
%                 '17'; '19'; '21'; '23'; '25'; '27';
%                 '29'; '31'; '33'; '35'};
%     title('Over-sampled matrix with M = 100, N = 1500, F = ');
% elseif kind == 2
%     k_label = { '5';   '7';  '9'; '11'; '13'; '15';
%                 '17'; '19'; '21'; '23'; '25'};  
%     title('Gaussian Matrix with M = 64, N = 1024, r = ');
% end
% 
% set(gca,'XTickLabel',k_label);
% 
% 
% xlabel('sparsity k');






