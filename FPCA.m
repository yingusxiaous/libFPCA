function V = FPCA(X, d, dampen)
% Fourier PCA by recursive partitioning of the space of the MN^{-1} matrix 
% and picking best basis for the subspace. This is for the under-determined 
% case. The difference between this and the advanced version is that we are
% not partitioning on the singular vectors here

% This function carries inside it a Bernoulli generative model for the
% data. For release, should strip it out and make it general.

% This only works for n source variables and n signal variables.

% P: matrix of basis vectors
% A: Original mixing matrix we're generating samples from (we're just
%       pretending we know it here).
% n: dimensionsality of the signal.
% m: number of samples

% What are the major differences here from recursiveFPCA?
% 1. Columns of A are no longer orthogonal.
% 2. We have to use SVD instead of EIG.

n = size(X,2);

% Pick our random vectors.
u = randn( [n, 1]);
u = dampen * u / norm(u);

M = buildFourthDerivative(X, u);

[W, ~, ~] = svd(real(M)); % This implies that W should be real
assert( all(all(W == real(W)))); % We will use this as the first map which 
%will hit all fourth moment tensors first.

% Changes the basis back to rankMap
V = W(:, 1:d) * FPCAhelper(W(:, 1:d), {eye(d)},{eye(d)}, X, d, dampen);
end