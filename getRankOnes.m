function V = getRankOnes(U)
% This function takes a complex matrix U of size n^2 * d, packs the columsn
% of U into a n*n matrices and extracts the best rank 1 vv^T approximation
% to the columns of U.

[n, d] = size(U);
n = sqrt(n);

V = zeros([n, d]);

for j=1:d    
    u = real(realProjection(U(:, j)));
    
    % Symmetrise u.
    A = reshape(u, n, n);
    A = A + A.';
    
    [P, ~] = eig(A);
    
    V(:,j) = P(:,1);    
end
end