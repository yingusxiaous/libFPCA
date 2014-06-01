% testGetRankOnes

A = randn(10);
A = bsxfun(@rdivide, A, sum(A.^2).^0.5);

B = kr(A,A);

C = getRankOnes(B);

[score, D] = basisEvaluation(A, C);

% Yeah, this works.