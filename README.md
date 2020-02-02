# pcg

Implementation of a Preconditioned Conjugate Gradient. Actual code:
	- Is parallelized
	- Expands and powers patterns
	- Solves small row CG with Lapack or implemented algorithm

TODO:
	- Optimize parallel CSC SpMV to avoid atomic operations
	- Implement Kaporin algorithm
	- Check for errors in powering and extending


To compile use make.
To run use: 
	- ./cg.out [][][][][][][][]


