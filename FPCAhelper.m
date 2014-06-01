function V = FPCAhelper(rankMap, basisChanges, projections, X, d, dampen)
% This is the helper function for FPCA. The point is that we've thrown away
% the excess dimensions at this point, and can now just run the FPCAhelper
% function recursively.

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

% Check for termination.
if size(projections{end},1) <= 1
   V = 1;
   return;
end

% Now do the recursive case.
[m, n] = size(X);

% Pick our random vectors.
u = randn( [n, 1]);
u = dampen * u / norm(u);

v = randn([n, 1]);
v = dampen * v / norm(v);

% Much more complicated tensor to build-up now.
S = rankMap.' * buildFourthDerivative(X,u) * rankMap;
T = rankMap.' * buildFourthDerivative(X,v) * rankMap;

M = S / T;

for j=1:length(projections)
    M = basisChanges{j} \ M * basisChanges{j};
    M = projections{j} * M * projections{j}.';
end

% Now we can try to find the top eigenvalues etc.
[U, D] = eig(real(M));
[eigenvalues, perm] = sort(diag(D));
U = U(:,perm); % Make sure to sort everything in this order too.

diff = eigenvalues(1:end-1) - eigenvalues(2:end);
[~,index] = max(diff);

basisChanges{end+1} = U;
smallProjections = projections;
bigProjections = projections;

identity = eye(length(eigenvalues)); % Are there going to be rank issues here?
smallProjections{end+1} = identity(1:index,:);
bigProjections{end+1} = identity(index+1:end,:);

smallProblem = FPCAhelper( rankMap, basisChanges, smallProjections, X, d, dampen);
bigProblem = FPCAhelper( rankMap, basisChanges, bigProjections, X, d, dampen);

% Reconstituting the solutions now. We're only unwinding one level of maps.
smallV = smallProjections{end}.' * smallProblem;
bigV = bigProjections{end}.' * bigProblem;

V = [ U\smallV U\bigV ];
end