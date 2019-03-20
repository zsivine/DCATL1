function  A = unif_ball(m,n) 
% generate m*n matirx, where column vectors are distributed uniformly on a
% m-dimensional unit ball. 

%% method 1: fail 
% % generate spherical coordinates 
% theta = rand(m-1,n,'double');
% theta(1:end-1,:) = theta(1:end-1,:)*pi;
% theta(end,:) = theta(end,:)*2*pi;
% 
% z = ones(1,n); 
% A = zeros(m,n);
% for i = 1:m-1 
%     A(i,:) = cos(theta(i,:)).*z;
%     z = z.*sin(theta(i,:));
% end
% A(end,:) = z;




%% method 2: 
% generate n Gaussian random variables x. 
% Then the distribution of unitized vectors is uniform over the surface $S^(n-1)$
x = randn(m,n);
A = bsxfun(@rdivide,x,sqrt(sum(x.^2,1)));




