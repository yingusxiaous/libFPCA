function v = realProjection( u, rescale)
% Input is a complex vector u. The output is (e^i \theta) * v where \theta
% is the optimum angle (in the interval [0, 2\pi]) that yields the largest
% projection to the reals R^n. We do this essentially by SVD

% The option rescale tells us whether we should rescale v to the norm of u
% or not

A = [ real(u) -imag(u)];

[~, ~, V] = svd(A); % Top right singular vector 
assert( all(size(V) == [2 2]));

phase1 = V(1,1) - 1i * V(1,2);
phase2 = V(1,1) + 1i * V(1,2);

% Why do I need to do this? Seems unnecessary, but oh well.
if norm( real(u * phase1)) > norm( real(u * phase2))
    v = phase1 * u;
else
    v = phase2 * u;
end

if nargin == 2 && rescale == true
    v = v / norm(v) * norm(u)
end

end
