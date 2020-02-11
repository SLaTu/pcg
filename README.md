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
- ./cg.out [RHS][Matrix][Repetitions][IMax][Tol][CG Mode][Pattern Mode] (Optional) [Add Percentage][Power]
- Matrix: .mtx file located in a folder in $../Inputs/Matrices/
  - RHS: 1 to input RHS file located in $../Inputs/RHS/
  - Repetitions: Number of times the CG is performed. The higher the better results.
  - IMax: Maximum number of iterations.
  - Tol: Tolerance value (e.g.: 1E-08)
  - CG Mode (1: Basic CG; 2: D⁻1 Preconditioning; 3: PCG)
  - Pattern Mode (1: Lower Triangle Pattern; 2: Lower Triangle; 3: Lower Triangle Pattern + Power + Extending; 4: Lower Triangle Pattern + Extending)
  - Add Percentage: Percentage of elements that will be added to the lower triangle pattern
  - Power: Powers the initial pattern to the Power value 
