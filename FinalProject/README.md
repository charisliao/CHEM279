# Welcome to the Final Project Repository for CHEM279! 

This README file will serve as a guide to the repository and the project.

## Repository Structure
The repository is organized as follows:
- `Bindir`: Contains the executable files for the CNDO/2 code and the input files for the CNDO/2 code. 
- `Include`: Contains the header files for the CNDO/2 code.
- `Libdir`: Contains the library files for the CNDO/2 code.
- `Source`: Contains the source files for the CNDO/2 code.
- `Test`: Contains the test files for the CNDO/2 code.
    - Some text files that has molecule names in the front are similar to the test files in hw4; however, functions are changed to optimized ones to ensure the correctness of the code.     
- `Makefile`: Contains the makefile for the CNDO/2 code. 


## Project Description
The goal of my project is to optimize CNDO/2 code so that it is efficient enough to treat very large molecules by eliminating overhead and optimizing the construction of `h` and `f` so that the compuation time is dominated by the linear algebra of diagonalizing `f`.     

## Project Plan
1. [Code Profiling](#code-profiling) 
2. [Memory Management](#memory-management) 
3. [Linear Algebra Optimization](#linear-algebra-optimization) 
4. [Compiler Optimization](#compiler-optimization) 
5. [Test and Validate](#test-and-validate)     

### Code Profiling 
Using profiling tools liek `gprf` or `perf` to identify the bottlenecks in my current CNDO/2 code. This could help pinpoint areas that I need to optimize. 

### Memory Management
There are two main approaches for improving memory management: 
- Minimize Memory
    - Optimize memory allocation and deallocation. Minimize unnecessary copies and reallocations, especially in loops where matrices are constructed.      
- Preallocate Memory   
    -  For matrices that grow during the computation, preallocate memory to their final size before performing operations. This reduces the need for dynamic memory resizing.

### Linear Algebra Optimization
Leverage optimized linear algebra libraries (e.g., BLAS, LAPACK) for matrix operations. Armadillo itself is designed for efficiency, but external libraries may provide additional speed.

### Compiler Optimization
Compilers can apply various optimizations, and compiling with different optimization levels (-O1, -O2, -O3) can influence runtime.

### Test and Validate
Demonstrate the effectiveness of your optimizations by testing the code on molecules of increasing size, ranging up to 1000 atoms.



