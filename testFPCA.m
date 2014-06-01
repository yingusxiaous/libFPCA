% testFPCA
% This should be the FPCA problem from the article.

% m is the number of samples
n = 3;
d = 4;
m = 1000000;
norm = 3;
numIterations = 100;

% Matrix with unit norm columns.
A = randn(3, 4);
A = bsxfun(@rdivide, A, sum(A.^2).^0.5);
X = generateData(m, A);

data = zeros( [numIterations,1]);
% Compute and evaluate basis using under-determined algorithm.
tic;
for j=1:numIterations
    j
    V = FPCA(X, d, norm);
    W = bsxfun(@rdivide, V, sum(V.* conj(V)).^0.5);
    [score, ~] = basisEvaluation(kr(A,A),W);
    data(j) = score;
    toc;
end
