function [score, B] = basisEvaluation(A, V, phaseCorrection)
% This function constructs the best bijection between vectors A and V
% according to the objective function max <A_i, V_j>^2. This works for
% complex vectors too. Note the inner product squared.

% A is the real basis that we're comparing against. 
% V is the (approximate) basis that we've computed.

% Score is the sum of <A_i,V_j>^2 for the optimal overlap
% B is the corrected up to sign and phase and permutation V

[~,m] = size(V);

weights = -abs(A' * V) .^ 2; % This has to be a Hermitian inner product.

% Hungarian algorithm to compute the max weight matching
[assignment, score] = munkres(weights);
score = - score / m;

B = V(:, assignment);

% Phase vs sign correction
if( nargin == 3 && phaseCorrection )
    % This is currently not implemented because we do realProjection and
    % getRankOnes outside the scope of this function.
    assert(false);
else
    % Sign correction only
    innerProducts = diag(real(A' * B) < 0);
    B(:,innerProducts) = -B(:,innerProducts);
end
end