% testFourierDerivatives

% Test the convergence of the Fourier derivatives.

n = 4;
d = 4;

m = 1000000;

% Matrix with unit norm columns.
A = randn([n, d]);
A = bsxfun(@rdivide, A, sum(A.^2).^0.5);
A = eye(n);

X = generateData(m, A);

u = zeros([n,1 ]);

[~, N] = buildFourierDerivatives(X, u, 4);

Y = generateData(m, A);

data = [];
data2 = [];

indices = 100000:100000:1000000;

for i= indices
    Z = Y(1:i,:);
    
    [~, M] = buildFourierDerivatives( Z, u, 2);
    data(end+1) = sqrt(sum( (M{2}-round(N{2})).^2));

    P = Z.' * Z / i;
    P = reshape(P, n^2, 1);
    data2(end+1) = sqrt(sum( (P-round(N{2})).^2));
end

clf;
hold on;
plot(data, 'ro');
plot(data2, 'b');