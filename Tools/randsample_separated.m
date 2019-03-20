% generate N_dimension vector with sparsity K, the distance of each nonzero
% element is greater than L
function supp = randsample_separated(N,K,L)
supp = randsample(N-L*(K-1),K);
supp = sort(supp);
supp = supp + (0:K-1)'*L;
end