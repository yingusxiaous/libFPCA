% testMain

% This is the main function to test the under-determined ICA code. See the
% readme for why the code is written this way.

% Parameters
n = 10; % Dimensionality of measurements.
d = 20; % Hidden latent signals.
m = 75000; % Total number of samples.
dampen = 1.6; % Coefficient used to change norm of random vectors.

% Mixing matrix with unit norm columns.
A = randn([n, d]);
A = bsxfun(@rdivide, A, sum(A.^2).^0.5);

X = generateData(m, A);

numIterations = 200;

% Now the metrics we keep to evaluate the algorithm:
data = zeros( [ numIterations, 1]); % This stores average overlap^2.

Ahat = zeros( size(A)); % This'll be our "averaged" basis. 

storedA = zeros( [n, numIterations * d]); 
% Keep track of all the cols of A found and then we'll cluster them at the
% end.

% Main loop to run the algorithm. This can be quite slow.
tic;
for j=1:numIterations  
    j
    
    V = underdeterminedFPCA(X, d, dampen);
    U = getRankOnes(V);
    
    % Store out data to evaluate the algorithm:
    [data(j), W] = basisEvaluation(Ahat, U);
    
    Ahat = Ahat + W;
    
    startIdx = (j-1) * d + 1;
    endIdx = j * d;
    
    storedA(:,startIdx:endIdx) = U;
    toc;
end

% Now compute these evaluation metrics:
mean(data)

Ahat = bsxfun( @rdivide, Ahat, sum( Ahat.^2).^0.5);
[score, final] = basisEvaluation(A, Ahat);
score

% Auto-pruning
final = [];

leftOver = storedA;

while size(final,2) < d
    size(final, 2)
    [basis, leftOver] = cleanUp(leftOver, d - size(final, 2));
    
    final = [ final basis];    
end
[b,c] = basisEvaluation(A, final);
median(diag(A' * c))
