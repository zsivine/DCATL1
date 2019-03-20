%% test different a for DCATL1 

M = 64; 
N = 256;
testnum = 100;
k_seq = 8:2:32;
thresh = 1e-3;

pm.gam = 1e-6;
pm.C = 1e-9;
pm.ep = 1e-6;
pm.maxit = 5000;
pm.maxoit = 20;
pm.otol = 1e-8;
pm.tol = 1e-5;
pm.del = pm.gam*10;

a = [0.1 0.3 1 2 10];
z = zeros(N,length(a));
succ = zeros(length(a),length(k_seq));

for j = 1:length(k_seq)
    j
    k = k_seq(j);
    for i = 1:testnum
        
        A = randn(M,N);
        x = zeros(N,1);
        x(randperm(N,k)) = randn(k,1);
        y = A*x;
        
        for l = 1:length(a)
            pm.a = a(l);
            z(:,l) = DCA_USC(A,y,pm,zeros(N,1));
            err = norm(z(:,l)-x)/norm(x);
            if err < thresh
                succ(l,j) = succ(l,j) + 1;
                fprintf('relative error : %4.3e with a = %4.3e \n',err,pm.a);
            end            
        end
        
        
        
        
    end
end

succ = succ./testnum;

% graph
figure
x = 1:length(k_seq);

plot(x,succ(1,:),'-or',  x,succ(2,:),'-*b',  x,succ(3,:),'-pg',...
     x,succ(4,:),'-dm',  x,succ(5,:),'-.k','linewidth',3);
axis([0 , length(k_seq)+1, 0, 1]);

hleg = legend('a=0.1','a=0.3','a=1','a=2','a=10');

set(hleg,'Location','SouthWest')
set(hleg,'Interpreter','none')
set(gca,'XTick',x);
k_label = { '8';   '10';  '12'; '14'; '16'; '18';
            '20'; '22'; '24'; '26'; '28'; '30'; '32'};
        
set(gca,'XTickLabel',k_label);

xlabel('sparsity k');
ylabel('success rate');







