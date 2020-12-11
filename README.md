# PMedian: Parallel Selection and Median Finding

A common statistical problem is that of finding the median element in
a set of data. This parallel MPI code implements an efficient,
randomized high-level parallel algorithms for finding the median given
a set of elements distributed across a parallel machine. In fact, our
algorithm solves the general selection problem that requires the
determination of the element of rank k, for an arbitrarily given
integer k. We use efficient techniques for distributing and coalescing
data as well as efficient combinations of task and data
parallelism. The algorithms have been coded in the message passing
standard MPI, and our experimental results from the IBM SP-2
illustrate the scalability and efficiency of our algorithm and improve
upon all the related experimental results known to the authors.

References:

D.A. Bader and J. Já Já. "Practical Parallel Algorithms for Dynamic Data Redistribution, Median Finding, and Selection," presented at the 10th International Parallel Processing Symposium (IPPS 96) Conference, Honolulu, HI, pp. 292-301, April 15-19, 1996.

D.A. Bader, "An Improved Randomized Selection Algorithm With an
Experimental Study," presented at the Second Workshop on Algorithm
Engineering and Experiments (ALENEX00), (sponsored by DIMACS, ACM
SIGACT, and SIAM), San Francisco, CA, January 7-8, 2000.
