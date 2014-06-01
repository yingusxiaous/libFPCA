function data = generateData(m, A)
% Generates Bernoulli data for an ICA model
n = size(A,2);
samples = 2 * ((rand([m,n]) < 0.5)) - 1;
data = samples * A.'; % Because all the values are transposed.