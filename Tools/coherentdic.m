function A = coherentdic(M,N,F)
% Generate an M x N oversampled partial Fourier matrix A,
% F is the refinement factor, see Fannjiang's paper:
% 'coherence-pattern guided compressive sensing with unresolved grids'
A = zeros(M,N);
r = rand(M,1);
l = 1:N;
for k = 1:M
    A(k,:) = sqrt(2/M)*cos(2*pi*r(k)*(l-1)/F);
end
end



