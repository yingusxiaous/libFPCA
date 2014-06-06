
libFPCA
Fourier PCA package
(Ying Xiao)

Introduction:

This library contains a variety of Matlab implementations for the
Fourier PCA algorithm. The main paper describing the approach can be
found at http://arxiv.org/abs/1306.5825

Code Structure:

The code is roughly divided into three parts:

1) Algorithms
I implemented a number of versions of the FPCA algorithms, in
ascending order of sophistication.

NaiveFPCA: Fully determined. This is essentially the 2nd order version
from the paper. This finds complex eigenvalues.

InverseFPCA: Fully determined. This takes the product of two matrices
MN^-1 and then goes after the eigenvalues.

RecursiveFPCA: Fully determined. This recursively goes after the best
eigenvalue gaps of the real part of the derivative and then partitions
the problem into subproblems over a smaller subspace.

underdeterminedFPCA: Underdetermined. This is the analogue of
InverseFPCA for the underdetermined case. We use this as a major
building block for the tests below. In many ways, this is our
canonical underdetermined algorithm.

FPCA: Underdetermined. This is a recursive partitioning version, much
like RecursiveFPCA except for the underdetermined case.

2) Evaluation metrics 
There are a number of functions here that we use to evaluate our
algorithm. The major one is basisEvaluation which constructs the best
possible bijection, in terms of inner product squared (whether it's a
scalar or Hermitian product is inconsequential now) between two bases
and then rearranges the second basis to be similar to the first
one. We use a max weight bipartite matching code here.

As auxiliary functions, we also have realProjections for individual
complex vectors which uses a SVD to compute the best possible rotation
in the complex plane (multiplication by exp(i \theta) to move as much
of the norm onto the real line as possible. Similary, getRankOnes does
this for a "squared" vector.

3) Test
The test scripts show you how to run the algorithms -- they are
essentially the unit tests for this package. Of particular interest is
testMain which runs the underdeterminedFPCA algorithm, note that
getRankOnes and evaluation metrics have to be run separately from the
FPCA algorithm.

Usage:

Run testMain and see what happens. Really, this contains the most
interesting part of the package.

Thoughts:

The implementation of the underdetermined code is highly nontrivial,
and its performance is actually a little disappointing. The best
performance I seem to be able to extract out of it is to iterate the
FPCA underdetermined algorithm many times, and then to run k-means
clustering on the output columns to obtain the most accurate estimates
of the columns.
