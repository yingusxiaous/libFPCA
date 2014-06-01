% testFourthDerivative

% Test buildFourthDerivative function using the 4th cumulant (we know what
% it's meant to be a priori).
n = 4;
d = 4;

m = 1000000;

% Matrix with unit norm columns.
A = randn([n, d]);
A = bsxfun(@rdivide, A, sum(A.^2).^0.5);
A = eye(n);

X = generateData(m, A);

u = zeros([n,1 ]);

N = buildFourthDerivative(X, u);

Y = generateData(m, A);

data = [];

for i=100000:100000:1000000
    M = buildFourthDerivative( Y(1:i,:), u);
    data(end+1) = sqrt(sum( sum( (M-round(N)).^2))); 
end