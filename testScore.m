% Test the scoring functions and bijection code at the end.
n = 10;
m = 40000;
dampen = 1.5;

numIterations = 100;

A = randn(n);
[B,~] = qr(A); % random unitary matrix

X = generateData(m, eye(n));

% Test for recursiveFPCA
V = recursiveFPCA( eye(n), X * B.', dampen);
[recursiveScore, ~] = basisEvaluation(B, V);

% Test for naiveFPCA
baseline = 0;
trueScore = 0;

for j=1:numIterations
    V = naiveFPCA( X, dampen);
    [score, ~] = basisEvaluation(eye(n), V);
    baseline = baseline + score;
    
    V = naiveFPCA(X * B.', dampen);
    [score,~] = basisEvaluation(B, V);
    trueScore = trueScore + score;
end
[baseline/numIterations trueScore/numIterations recursiveScore]