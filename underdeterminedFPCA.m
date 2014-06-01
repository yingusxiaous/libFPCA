function V = underdeterminedFPCA(X, d, dampen)
% This does the under-determined fourth order Fourier PCA calculation in 
% our paper. The problem with this code is that hitting all the cols of A
% simultaneously seems to be a low-probability event. My guess is that the
% SVD to compute W is a low probability event. So we have to do a lot of
% engineering around this.

% X is the dataset: m rows of data, each of dimension n.

% d is the number of latent variables, n is the number of observed dims.

% dampen is the norm of the random vectors chosen (perhaps these should
% even be different?!?

% V is the returned basis, we have not corrected them for complex phases.

% Pick our random vectors.
n = size(X, 2);

u =  randn( [n,1]);
u = dampen * u / norm(u);

v = randn([n, 1]);
v = dampen * v / norm(v);

P = buildFourthDerivative(X,u);
Q = buildFourthDerivative(X,v);

% SVD step
[W,~,~] = svd(P);
W = W(:,1:d);

P = W.' * P * W;
Q = W.' * Q * W;

% Inverse and 
M = P / Q;
[V,~] = eig(M);
V = W * V;
end