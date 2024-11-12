# README

## Project Title: Parallelization of a 3D Poisson Problem

### Overview
This project involves parallelizing a 3D Poisson equation solver using a Preconditioned Conjugate Gradient (PCG) approach. The work focuses on optimizing the solver's performance and convergence for various mesh sizes by employing preconditioning techniques like Symmetric Successive Over-Relaxation (SSOR) and Additive Schwarz. 

### Key Components
1. **Preliminary Profiling & Optimization**: Analyzed computational bottlenecks, identifying preconditioning steps, matrix-vector multiplications, and dot products as key areas for optimization.
2. **Parallelization Strategy**:
   - **Unpreconditioned Code**: Optimized core functions (e.g., matrix-vector multiplication) using OpenMP, ensuring scalability for larger meshes.
![image](https://github.com/user-attachments/assets/3d08dd70-ce1a-4137-8928-c3e4f32a9f84)
   - **Additive Schwarz Preconditioner**: Implemented Schwarz domain decomposition, efficiently managing boundaries and data copying to optimize multi-threaded execution.
![image](https://github.com/user-attachments/assets/8b548108-4aba-43c7-9212-8381beb571c9)
3. **Scalability Analysis**: Conducted strong and weak scaling tests, evaluating performance across thread counts and problem sizes to identify efficiency thresholds.
4. **Advanced Analysis on Convergence Properties**: Explored convergence behaviors with different RHS configurations, highlighting the impact of frequency content on solver efficiency.
![image](https://github.com/user-attachments/assets/0b20141f-7d46-40e5-ab3b-6ae887e4d67d)
![image](https://github.com/user-attachments/assets/72ec4213-24e7-422a-a982-8e39287abee2)
### Results
The solver demonstrated significant speed-ups with parallelization, particularly in SSOR preconditioning. Scalability limitations were noted at higher thread counts, mainly due to synchronization and memory overhead. Convergence analyses showed that high-frequency modes are reduced quickly, while low-frequency components require more iterations, underlining the need for tailored preconditioning strategies.

### Conclusion
This project offers an efficient 3D Poisson solver optimized for parallel execution, showcasing the benefits of combining advanced preconditioning with careful profiling and optimization. The code's adaptability to various mesh sizes and thread counts makes it suitable for large-scale scientific applications.

### Instructions
- **Dependencies**: OpenMP, C compiler with parallel support.
- **Execution**: Adjust thread count and mesh size via configuration file or command line arguments.
- **Usage**: Run `make` to build the project, then execute `./poisson_solver` with desired options.
