%% level lines for l0 , l1 and transformed l1 \rho_a
lamda = 0.5:0.1:1;
N = length(lamda);
M = 50;
xx = zeros(N+1,2*M+2);
yy = xx;

%% l1 
l = 2/M;
x1 = -1:l:1;
x2 = 1:-l:-1;
y1 = 1-abs(x1);
y2 = abs(x2) - 1;
xx(1,:) = [x1 x2]; yy(1,:) = [y1 y2]; 
figure; 
plot(xx(1,:),yy(1,:));

%% transformed l1 with different parameters
for i = 1:N
    L = lamda(i)/(2-lamda(i));
    l = 2*L/M;
    x1 = -L:l:L;
    x2 = L:-l:-L;
    lam = lamda(i) - 2.*abs(x1)./(1+abs(x1));
    y1 = lam./(2-lam);
    lam = lamda(i) - 2.*abs(x2)./(1+abs(x2));
    y2 = lam./(lam-2);
    xx(i+1,:) = [x1 x2]; yy(i+1,:) = [y1 y2]; 
end

plot(xx(2:end,:)', yy(2:end,:)')

axis([-1 1 -1 1]);



