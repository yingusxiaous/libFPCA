% testRealProjection

% This test whether the realProjection function actually gives us the
% optimal answer. Honestly, the iterative approach is probably just as
% good, but the SVD is much faster.

n = 20;

u = randn([n, 1]) + 1i * randn([n,1]);

v = realProjection(u);

range = 0:0.01:2 * pi;
data = zeros(size(range));

for j=1:length(range)
    theta = range(j);
    
    data(j) = norm( real( u * exp(1i * theta)));
end

plot(range, data);

[ sqrt(v' * v) sqrt(real(v).' * real(v)) max(data)]