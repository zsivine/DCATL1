%% numerical test 1: Gaussian matrix 
% M = 64, N = 1024, r = {0:0.2:1}
% k = 5:2:25
M = 64;
N = 1024;
k = 5:2:25;
for p = [0, 0.2, 0.4, 0.6 , 0.8]
    test_DCA(p,2,M,N,k);
end

%% numerical test 2: Over-sampled matrix 
M = 100; 
N = 1500;
k = 25:2:35;
% F = 2:2:20

% for p = [2,4,6,10,12,16,20]
% for p = [6,10,12,16,20]
for p = [20]
    test_DCA(p,1,M,N,k);
end


% succnum = zeros(5,6);
% for i = 1:5
%     succnum(i,:) = DCATL1_RL_sparsity(i);
% end









