% testDampening

% This computes the optimal dampening coefficient to be used in all
% subsequent codes.
n = 10;
m = 40000;
repetitions = 10; 

norms = 0.1:0.1:5.0;
data = zeros([length(norms),2]);

A = randn(n);
[A,~] = qr(A);

for k=1:repetitions
    k
    j=1;
    X = generateData(m, A);
    
    for dampen=norms
        V = naiveFPCA( X, dampen);
        [scoring,~] = basisEvaluation(A, V);
        data(j,1) = data(j,1) + scoring;
        
        
        V = recursiveFPCA( eye(n), X, dampen);
%        V = inverseFPCA(X, dampen);
        [scoring,~] = basisEvaluation(A, V);
        data(j,2) = data(j,2) + scoring;
        
        j = j + 1;
    end
end

data = data/repetitions;

clf;
hold on;
plot(norms,data(:,1));
plot(norms, data(:,2),'r');
xlabel('Norm of random vector');
ylabel('Average overlap with basis vectors');
title('Choice of norm of random vector using 40000 samples in 10 dimensions');
legend( 'Naive', 'Recursive');
