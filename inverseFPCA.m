function V = inverseFPCA(X, dampen)

% Inverse Fourier PCA algorithm where we compute two second derivative
% matrices and then invert one of them. This is a necessary intermediate
% step because I'm not sure of the conditioning of the matrices.

% This function carries inside it a Bernoulli generative model for the
% data. For release, should strip it out and make it general.

% This only works for n source variables and n signal variables.

% A: Original mixing matrix we're generating samples from (we're just
%       pretending we know it
% n: dimensionsality of the signal.
% m: number of samples

[m, n] = size(X);

% Hessian 1
u =  randn( [n,1]);
u = dampen * u / norm(u);

weights = exp( 1i * X * u);
reweighted = bsxfun(@times, X, weights);

zeroth = sum( weights )/m; % zeroth order
first = sum( reweighted ).'/m;
second = reweighted.' * X/m;

psi = second * zeroth - first * first.';

% Hessian 2
v =  randn( [n,1]);
v = dampen * v / norm(v);

weights = exp( 1i * X * v);
reweighted = bsxfun(@times, X, weights);

zeroth = sum( weights )/m; % zeroth order
first = sum( reweighted ).'/m;
second = reweighted.' * X/m;

phi = second * zeroth - first * first.';

M = psi / phi;

[V, ~] = eig(M);
end