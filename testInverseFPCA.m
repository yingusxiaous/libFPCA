n = 10;
m = 100000;

A = randn(n+1);
[A,~] = qr(A);
B = A(:,1:end-1);

X = generateData(m, A);

% Over-determined if we stick it through over-determined
V = overdeterminedFPCA(X, n, 1.5);
basisEvaluation(B,V)

% SVD on A first
[C, ~, ~] = svd(B);
C = C(:, 1:n);

W = inverseFPCA(X * C, 1.5);
basisEvaluation(B, C * W)