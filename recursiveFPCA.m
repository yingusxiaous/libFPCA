function V = recursiveFPCA(P, X, dampen)
% Fourier PCA by recursive partitioning of the space and picking best
% splits of the subspace. This is dramatically better than naive FPCA.

% This function carries inside it a Bernoulli generative model for the
% data. For release, should strip it out and make it general.

% This only works for n source variables and n signal variables.

% P: matrix of basis vectors
% A: Original mixing matrix we're generating samples from (we're just
%       pretending we know it).
% n: dimensionsality of the signal.
% m: number of samples

d = size(P, 2);
m = size(X,1);

if d == 1
    V = P;
    return;
end

Y = X * P;

% Now run our usual Fourier PCA routine.
u =  randn( [d,1]);
u = dampen * u / norm(u);

weights = exp( 1i * Y * u);
reweighted = bsxfun(@times, Y, weights);

zeroth = sum( weights )/m; % zeroth order
first = sum( reweighted ).'/m;
second = reweighted.' * Y/m;

psi = second * zeroth - first * first.';

[U,D] = eig(real(psi)); % U are the eigenvectors.

% Sort the eigenvalues and eigenvectors in ascending order.
[eigenvalues,perm] = sort(diag(D));
U = U(:,perm);

diff = eigenvalues(2:end) - eigenvalues(1:end-1);
[~,index] = max(diff);

smallBasis = U(:,1:index);
bigBasis = U(:,index+1:end);

smallProblem = recursiveFPCA( P * smallBasis, X, dampen);
bigProblem = recursiveFPCA( P * bigBasis, X, dampen);

V = [ smallProblem bigProblem ];
end