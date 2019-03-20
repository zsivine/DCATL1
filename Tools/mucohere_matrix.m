function z = mucohere_matrix(A)
% computer mutual coherence of matrix A
% normlize columns of A
[m,n] = size(A);
for j = 1:n
    t = norm(A(:,j),2);
    A(:,j) = A(:,j)./t;
end

mc = A'*A;
mc(eye(n)~=0) = 0;
mc = abs(mc(:));
z = max(mc);


