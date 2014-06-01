function M = buildFourthDerivative(X, u)
% A whole bunch of fancy Matlab tricks to build the fourth derivative
% tensor of the logarithm of the Fourier transform evaluated at u.

[zeroth, moments] = buildFourierDerivatives(X, u, 4);

n = size(X,2);

% And now for the deadly formula...

% (4) term
M = reshape(zeroth^3 * moments{4}, n, n, n, n);

% The rest is some clever Matlab coding -- we look at everything as a
% tensor, and just permute the arguments to achieve the right symmetries.
% It seems pretty damn fast.

% (3,1) term
tensor = reshape(-1 * zeroth^2 * kr(moments{3},moments{1}), n, n, n, n);
M = M + tensor ...
    + permute(tensor, [1 2 4 3]) ...
    + permute(tensor, [1 4 2 3]) ...
    + permute(tensor, [4 1 2 3]);

% (2,2) term
tensor = reshape(-1 * zeroth^2 * kr(moments{2}, moments{2}), n, n, n, n);
M = M + tensor ...
    + permute(tensor, [1 3 2 4]) ...
    + permute(tensor, [1 3 4 2]);

% (2,1,1) term.
tensor = reshape(2 * zeroth * kr(moments{2}, moments{1}, moments{1}), n, n, n, n);
M = M + tensor ...
    + permute(tensor, [1 3 2 4]) ...
    + permute(tensor, [1 3 4 2]) ...
    + permute(tensor, [3 1 2 4]) ...
    + permute(tensor, [3 1 4 2]) ...
    + permute(tensor, [3 4 1 2]);

% (1,1,1,1) term
tensor = reshape( -6 * kr(moments{1}, moments{1}, moments{1}, moments{1}), n, n, n, n);
M = M + tensor;

% Pack it into an n^2 * n^2 matrix and we're done, finally!
M = reshape(M, n^2, n^2).'; 
clear('moments');
end