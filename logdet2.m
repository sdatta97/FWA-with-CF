function c = logdet2(X)
%log2 det(I + X) for X Hermitian PSD (up to roundoff), via Cholesky:
%the MIMO capacity term of the compute_link_rates_* functions
X = (X + X')/2;
c = 2*sum(log2(real(diag(chol(eye(size(X,1)) + X)))));
end
