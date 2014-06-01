function [zeroth, M] = buildFourierDerivatives(X, u, k, chunkSize)
% This computes terms of the form E{ x_i ...x_j exp(i u^T x)} and stores
% them as vectors.

% X is the dataset (rows of data).
% u is the direction (vector) we're interested in.
% k is maximum number of Fourier derivatives we'd like to compute
% chunkSize is the batch size of our averaging operation. It's very
% expensive to build all the outer product tensors simultaneously, so we do
% it in batches instead.

% M is a cell array of Fourier moment tensors of order 1 to k.

[m, n] = size(X);

if nargin < 4
    % Really, it's a much better idea to let the program decide how big to
    % make the chunks. If this is done incorrectly, the code slows down
    % orders of magnitudes due to paging from disk.
    [~, systemView] = memory();
    chunkSize = floor(systemView.PhysicalMemory.Available/ (8 * 1.2 * 2 * n^k));
    
    % Magic constants all round: 8 is for the size of the doubles, n^k is
    % for the size of the largest tensor, 2 is because we have to make all
    % previous tensors as well, 1.2 is for overhead. The point is that we
    % want to make everything fit in memory, once is goes to hdd, this
    % becomes much slower.
end

M = cell(k,1);

% Generates the value of the Fourier transform
weights = exp( 1i * X * u);
reweighted = bsxfun(@times, X, weights).';
Y = X.';

zeroth = sum(weights)/m;

M{1} = sum(reweighted,2)/m;

% Now compute all the higher derivatives. We have to do this in chunks,
% otherwise this uses up all our physical memory and goes to disk.
for j=2:k
    M{j} = zeros([n^j,1]);
end

numChunks = ceil(m/chunkSize);

for i=1:numChunks
    start = 1 + (i-1)*chunkSize;
    finish = min(i * chunkSize, m);
    Z=Y(:,start:finish);
    reweightedZ = reweighted(:,start:finish);
    
    for j=2:k
        % Using the Khatri-Rao product seems to be the fastest way to do this. 
        reweightedZ = kr(reweightedZ, Z);
        M{j} = M{j} + sum(reweightedZ, 2)/m;
    end
end
end