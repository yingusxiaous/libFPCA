function V = naiveFPCA(X, dampen)
% Naive Fourier PCA algorithm where we try to find the entire basis all at
% once. This is dramatically worse (in terms of sample complexity) than the
% recursive FPCA function.

% This function carries inside it a Bernoulli generative model for the
% data. For release, should strip it out and make it general.

% This only works for n source variables and n signal variables.

% A: Original mixing matrix we're generating samples from (we're just
%       pretending we know it
% n: dimensionsality of the signal.
% m: number of samples

[m, n] = size(X);

% Now run our usual Fourier PCA routine.
u =  randn( [n,1]);
u = dampen * u / norm(u);

weights = exp( 1i * X * u);
reweighted = bsxfun(@times, X, weights);

zeroth = sum( weights )/m; % zeroth order
first = sum( reweighted ).'/m;
second = reweighted.' * X/m;

psi = second * zeroth - first * first.';

[V, ~] = eig(psi);
end