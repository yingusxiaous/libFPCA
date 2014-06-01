% testUnderdetermined

% This test tries to determined what the best choice for the dampening
% parameter is. It seems pretty consistent even if I changed n and d.

% m is the number of samples
n = 4;
d = 6;
m = 1000000;

% Matrix with unit norm columns.

norms = 0.1:0.1:5;
numIterations = 1;
data = zeros( [ length(norms), numIterations]);

A = randn([n, d]);
A = bsxfun(@rdivide, A, sum(A.^2).^0.5);
X = generateData(m, A);

tic;
for i=1:length(norms)
    
    dampen = norms(i);
    
    for j=1:numIterations
        V = underdeterminedFPCA(X, d, dampen);
        [score, ~] = basisEvaluation(kr(A,A),V);
        data(i,j) = score;
        toc;
    end
end
