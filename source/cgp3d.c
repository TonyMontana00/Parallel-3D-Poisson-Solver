#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <omp.h>


#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "pcg.h"
#include "params.h"
// #include <cstring>  // For memset()
#include "timing.h"


//#define CLOCK CLOCK_REALTIME


void mul_poisson3d(int N, void* data, 
                   double* restrict Ax, 
                   double* restrict x)
{
    #define X(i,j,k) (x[((k)*n+(j))*n+(i)])
    #define AX(i,j,k) (Ax[((k)*n+(j))*n+(i)])



    int n = *(int*) data;
    int inv_h2 = (n-1)*(n-1);
    #ifdef PARALLEL
         #pragma omp parallel for shared(Ax)
    #endif
    for (int k = 0; k < n; ++k) {
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                double xx = X(i,j,k);
                double xn = (i > 0)   ? X(i-1,j,k) : 0;
                double xs = (i < n-1) ? X(i+1,j,k) : 0;
                double xe = (j > 0)   ? X(i,j-1,k) : 0;
                double xw = (j < n-1) ? X(i,j+1,k) : 0;
                double xu = (k > 0)   ? X(i,j,k-1) : 0;
                double xd = (k < n-1) ? X(i,j,k+1) : 0;
                AX(i,j,k) = (6*xx - xn - xs - xe - xw - xu - xd)*inv_h2;
            }
        }
    }

    #undef AX
    #undef X
}


// void mul_poisson3d(int N, void* data, 
//                    double* restrict Ax, 
//                    double* restrict x)
// {
    
//     #define X(i,j,k) (x[((k)*n+(j))*n+(i)])
//     #define AX(i,j,k) (Ax[((k)*n+(j))*n+(i)])

//     int n = *(int*) data;
//     int inv_h2 = (n-1)*(n-1);

//     #ifdef PARALLEL
//         #pragma omp parallel shared(Ax, x, n, inv_h2)
//         {
//             // int thread_id = omp_get_thread_num();
//             // int num_threads = omp_get_num_threads();

//             // #pragma omp single
//             // {
//             //     printf("Parallel execution with %d threads\n", num_threads);
//             //     fflush(stdout);
//             // }

//             //#pragma omp for collapse(2)
//             for (int k = 0; k < n; ++k) {
//                 for (int j = 0; j < n; ++j) {
//                     for (int i = 0; i < n; ++i) {
//                         // Computation
//                         double xx = X(i,j,k);
//                         double xn = (i > 0)   ? X(i-1,j,k) : 0;
//                         double xs = (i < n-1) ? X(i+1,j,k) : 0;
//                         double xe = (j > 0)   ? X(i,j-1,k) : 0;
//                         double xw = (j < n-1) ? X(i,j+1,k) : 0;
//                         double xu = (k > 0)   ? X(i,j,k-1) : 0;
//                         double xd = (k < n-1) ? X(i,j,k+1) : 0;
//                         AX(i,j,k) = (6*xx - xn - xs - xe - xw - xu - xd)*inv_h2;
//                     }
//                 }
//             }
//         }
//     #else
//         // Sequential version without OpenMP
//         for (int k = 0; k < n; ++k) {
//             for (int j = 0; j < n; ++j) {
//                 for (int i = 0; i < n; ++i) {
//                     double xx = X(i,j,k);
//                     double xn = (i > 0)   ? X(i-1,j,k) : 0;
//                     double xs = (i < n-1) ? X(i+1,j,k) : 0;
//                     double xe = (j > 0)   ? X(i,j-1,k) : 0;
//                     double xw = (j < n-1) ? X(i,j+1,k) : 0;
//                     double xu = (k > 0)   ? X(i,j,k-1) : 0;
//                     double xd = (k < n-1) ? X(i,j,k+1) : 0;
//                     AX(i,j,k) = (6*xx - xn - xs - xe - xw - xu - xd)*inv_h2;
//                 }
//             }
//         }
//     #endif

//     #undef AX
//     #undef X
// }



/*@T
 * \section{Preconditioners for the Laplacian}
 *
 * \subsection{The identity preconditioner}
 *
 * The simplest possible preconditioner is the identity ([[pc_identity]]):
 *@c*/

#ifdef PARALLEL
void parallel_copy(int n, double *Ax, double *x) 
{ 
#pragma omp parallel 
    { 
        int n_threads = omp_get_num_threads(); 
        int id = omp_get_thread_num(); 
        int share = (n + n_threads - 1) / n_threads; 
        if (id * share < n && (id + 1) * share >= n) 
        { 
            int rem = n - id * share; 
            memcpy(&(Ax[id * share]), &(x[id * share]), rem * sizeof(double)); 
        } 
        else if (id * share <= n) 
        { 
            memcpy(&(Ax[id * share]), &(x[id * share]), share * 
sizeof(double)); 
        } 
    } 
}
#endif

void pc_identity(int n, void *data, double *Ax, double *x) 
{ 
#ifdef PARALLEL 
    parallel_copy(n, Ax, x); 
#else 
    memcpy(Ax, x, n * sizeof(double)); 
#endif 
}


typedef struct pc_ssor_p3d_t {
    int n;          /* Number of points in one direction on mesh */
    double omega;   /* SSOR relaxation parameter */
} pc_ssor_p3d_t;


#define BS 20

void ssor_forward_sweep(int n, int i1, int i2, int j1, int j2, int k1, int k2,
                        double* restrict Ax, double w)
{
    #define AX(i,j,k) (Ax[((k)*n+(j))*n+(i)])
    
#ifdef PARALLEL
    // Parallel version using OpenMP tasks with dependencies
    const int num_blocks_i = (i2 - i1 + BS - 1) / BS;
    const int num_blocks_j = (j2 - j1 + BS - 1) / BS;
    const int num_blocks_k = (k2 - k1 + BS - 1) / BS;

    // Allocate and initialize dependency matrix
    char *dep_matrix = malloc(num_blocks_i * num_blocks_j * num_blocks_k * sizeof(char));
    memset(dep_matrix, 0, num_blocks_i * num_blocks_j * num_blocks_k * sizeof(char));

    // Dummy variable for non-existent dependencies
    char dummy = 0;

    // Parallel region
    #pragma omp parallel
    {
        #pragma omp single nowait
        {
            for (int k = 0; k < num_blocks_k; ++k)
            {
                for (int j = 0; j < num_blocks_j; ++j)
                {
                    for (int i = 0; i < num_blocks_i; ++i)
                    {
                        // Compute dependency pointers
                        char *dep_in_i = (i > 0) ? &dep_matrix[((k) * num_blocks_j + (j)) * num_blocks_i + (i - 1)] : &dummy;
                        char *dep_in_j = (j > 0) ? &dep_matrix[((k) * num_blocks_j + (j - 1)) * num_blocks_i + (i)] : &dummy;
                        char *dep_in_k = (k > 0) ? &dep_matrix[((k - 1) * num_blocks_j + (j)) * num_blocks_i + (i)] : &dummy;
                        char *dep_out = &dep_matrix[((k) * num_blocks_j + (j)) * num_blocks_i + (i)];

                        // Create task with dependencies
                        #pragma omp task depend(in: dep_in_i[0], dep_in_j[0], dep_in_k[0]) depend(out: dep_out[0])
                        {
                            // Compute block boundaries
                            int ibegin = i1 + i * BS;
                            int iend = (ibegin + BS < i2) ? ibegin + BS : i2;

                            int jbegin = j1 + j * BS;
                            int jend = (jbegin + BS < j2) ? jbegin + BS : j2;

                            int kbegin = k1 + k * BS;
                            int kend = (kbegin + BS < k2) ? kbegin + BS : k2;

                            // Perform the sweep on the block
                            for (int kk = kbegin; kk < kend; ++kk) {
                                for (int jj = jbegin; jj < jend; ++jj) {
                                    for (int ii = ibegin; ii < iend; ++ii) {
                                        double xx = AX(ii,jj,kk);
                                        double xn = (ii > 0)   ? AX(ii-1,jj,kk) : 0;
                                        double xe = (jj > 0)   ? AX(ii,jj-1,kk) : 0;
                                        double xu = (kk > 0)   ? AX(ii,jj,kk-1) : 0;
                                        AX(ii,jj,kk) = (xx + xn + xe + xu) / 6 * w;
                                    }
                                }
                            }

                            // Mark the block as completed
                            dep_out[0] = 1;
                        } // End of task
                    }
                }
            }
        } // End of single
    } // End of parallel region

    // Free dependency matrix
    free(dep_matrix);

#else
    // Sequential version
    for (int k = k1; k < k2; ++k) {
        for (int j = j1; j < j2; ++j) {
            for (int i = i1; i < i2; ++i) {
                double xx = AX(i,j,k);
                double xn = (i > 0)   ? AX(i-1,j,k) : 0;
                double xe = (j > 0)   ? AX(i,j-1,k) : 0;
                double xu = (k > 0)   ? AX(i,j,k-1) : 0;
                AX(i,j,k) = (xx + xn + xe + xu) / 6 * w;
            }
        }
    }
#endif

    #undef AX
}


// void ssor_backward_sweep(int n, int i1, int i2, int j1, int j2, int k1, int k2, 
//                          double* restrict Ax, double w)
// {
//     #define AX(i,j,k) (Ax[((k)*n+(j))*n+(i)])
//     for (int k = k2-1; k >= k1; --k) {
//         for (int j = j2-1; j >= j1; --j) {
//             for (int i = i2-1; i >= i1; --i) {
//                 double xx = AX(i,j,k);
//                 double xs = (i < n-1) ? AX(i+1,j,k) : 0;
//                 double xw = (j < n-1) ? AX(i,j+1,k) : 0;
//                 double xd = (k < n-1) ? AX(i,j,k+1) : 0;
//                 AX(i,j,k) = (xx+xs+xw+xd)/6*w;
//             }
//         }
//     }
//     #undef AX
// }

void ssor_backward_sweep(int n, int i1, int i2, int j1, int j2, int k1, int k2,
                         double* restrict Ax, double w)
{
    #define AX(i,j,k) (Ax[((k)*n+(j))*n+(i)])

#ifdef PARALLEL
    // Parallel version using OpenMP tasks with dependencies
     // Block size, adjust as needed
    int num_blocks_i = (i2 - i1 + BS - 1) / BS;
    int num_blocks_j = (j2 - j1 + BS - 1) / BS;
    int num_blocks_k = (k2 - k1 + BS - 1) / BS;

    int total_blocks = num_blocks_i * num_blocks_j * num_blocks_k;

    char *dep_matrix = malloc(total_blocks * sizeof(char));
    memset(dep_matrix, 0, total_blocks * sizeof(char));

    char dummy = 0; // Dummy variable for dependencies at boundaries

    #pragma omp parallel
    {
        #pragma omp single nowait
        {
            for (int kk = num_blocks_k - 1; kk >= 0; --kk)
            {
                for (int jj = num_blocks_j - 1; jj >= 0; --jj)
                {
                    for (int ii = num_blocks_i - 1; ii >= 0; --ii)
                    {
                        // Compute dependency pointers
                        char *dep_in_i = (ii < num_blocks_i - 1) ? &dep_matrix[(kk * num_blocks_j + jj) * num_blocks_i + (ii + 1)] : &dummy;
                        char *dep_in_j = (jj < num_blocks_j - 1) ? &dep_matrix[(kk * num_blocks_j + (jj + 1)) * num_blocks_i + ii] : &dummy;
                        char *dep_in_k = (kk < num_blocks_k - 1) ? &dep_matrix[((kk + 1) * num_blocks_j + jj) * num_blocks_i + ii] : &dummy;
                        char *dep_out = &dep_matrix[(kk * num_blocks_j + jj) * num_blocks_i + ii];

                        // Create the task with dependencies
                        #pragma omp task depend(in: dep_in_i[0], dep_in_j[0], dep_in_k[0]) depend(out: dep_out[0])
                        {
                            int ibegin = i1 + ii * BS;
                            int iend = (ibegin + BS < i2) ? ibegin + BS : i2;

                            int jbegin = j1 + jj * BS;
                            int jend = (jbegin + BS < j2) ? jbegin + BS : j2;

                            int kbegin = k1 + kk * BS;
                            int kend = (kbegin + BS < k2) ? kbegin + BS : k2;

                            // Perform the sweep on the block
                            for (int k = kend - 1; k >= kbegin; --k)
                            {
                                for (int j = jend - 1; j >= jbegin; --j)
                                {
                                    for (int i = iend - 1; i >= ibegin; --i)
                                    {
                                        double xx = AX(i,j,k);
                                        double xs = (i < n - 1) ? AX(i + 1,j,k) : 0;
                                        double xw = (j < n - 1) ? AX(i,j + 1,k) : 0;
                                        double xd = (k < n - 1) ? AX(i,j,k + 1) : 0;
                                        AX(i,j,k) = (xx + xs + xw + xd) / 6 * w;
                                    }
                                }
                            }

                            // Mark the block as completed
                            dep_out[0] = 1;
                        } // End of task
                    }
                }
            }
        }
    }

    free(dep_matrix);

#else
    // Sequential version
    for (int k = k2 - 1; k >= k1; --k) {
        for (int j = j2 - 1; j >= j1; --j) {
            for (int i = i2 - 1; i >= i1; --i) {
                double xx = AX(i,j,k);
                double xs = (i < n - 1) ? AX(i + 1,j,k) : 0;
                double xw = (j < n - 1) ? AX(i,j + 1,k) : 0;
                double xd = (k < n - 1) ? AX(i,j,k + 1) : 0;
                AX(i,j,k) = (xx + xs + xw + xd) / 6 * w;
            }
        }
    }
#endif

    #undef AX
}


void ssor_diag_sweep(int n, int i1, int i2, int j1, int j2, int k1, int k2, 
                     double* restrict Ax, double w)
{
    #define AX(i,j,k) (Ax[((k)*n+(j))*n+(i)])
    #ifdef PARALLEL
    #pragma omp parallel for shared(Ax)
    #endif
    for (int k = k1; k < k2; ++k)
        for (int j = j1; j < j2; ++j)
            for (int i = i1; i < i2; ++i)
                AX(i,j,k) *= (6*(2-w)/w);
    #undef AX
}

/* @T
 *
 * Finally, the [[pc_ssor_poisson3d]] function actually applies the
 * preconditioner.
 *@c*/
void pc_ssor_poisson3d(int N, void* data, 
                       double* restrict Ax, 
                       double* restrict x)
{
    #define FORWARD_SWEEP_TIME 10
    #define DIAG_SWEEP_TIME    11
    #define BACKWARD_SWEEP_TIME 12
    #define COPY_TIME 13

    // double forward_time = 0;
    // double diag_time = 0;
    // double backward_time = 0;

    pc_ssor_p3d_t* ssor_data = (pc_ssor_p3d_t*) data;
    int n = ssor_data->n;
    double w = ssor_data->omega;

    tic(COPY_TIME);
    // #ifdef PARALLEL
    // parallel_copy(n, Ax, x);
    // #else
    memcpy(Ax, x, N*sizeof(double));
    // #endif
    total_copy_time_ssor += toc(COPY_TIME);
    
    // ssor_forward_sweep (n, 0,n, 0,n, 0,n, Ax, w);
    // ssor_diag_sweep    (n, 0,n, 0,n, 0,n, Ax, w);
    // ssor_backward_sweep(n, 0,n, 0,n, 0,n, Ax, w);

    // Time the forward sweep
    tic(FORWARD_SWEEP_TIME);
    ssor_forward_sweep(n, 0, n, 0, n, 0, n, Ax, w);
    total_forward_time += toc(FORWARD_SWEEP_TIME);

    // Time the diagonal sweep
    tic(DIAG_SWEEP_TIME);
    ssor_diag_sweep(n, 0, n, 0, n, 0, n, Ax, w);
    total_diag_time += toc(DIAG_SWEEP_TIME);

    // Time the backward sweep
    tic(BACKWARD_SWEEP_TIME);
    ssor_backward_sweep(n, 0, n, 0, n, 0, n, Ax, w);
    total_backward_time += toc(BACKWARD_SWEEP_TIME);

}

typedef struct pc_schwarz_p3d_t {
    int n;           /* Number of mesh points on a side */
    int overlap;     /* Number of points through the overlap region */
    double omega;    /* SSOR relaxation parameter */
    double* scratch; /* Scratch space used by the preconditioner */
} pc_schwarz_p3d_t;


/*@T
 *
 * In order to compute independently on overlapping subdomains, we first
 * get the local data from the vector to which we're applying the
 * preconditioner; then we do an inexact solve on the local piece of
 * the data; and then we write back the updates from the solve.
 * The data motion is implemented in [[schwarz_get]] and [[schwarz_add]].
 *@c*/
// void schwarz_get(int n, int i1, int i2, int j1, int j2, int k1, int k2, 
//                  double* restrict x_local, 
//                  double* restrict x)
// {
//     #define X(i,j,k) (x[((k)*n+(j))*n+(i)])
//     #define XL(i,j,k) (x_local[((k)*n+(j))*n+(i)])

//     for (int k = k1; k < k2; ++k)
//         for (int j = j1; j < j2; ++j)
//             for (int i = i1; i < i2; ++i)
//                 XL(i,j,k) = X(i,j,k);

//     if (k1 > 0)
//         for (int j = j1; j < j2; ++j)
//             for (int i = i1; i < i2; ++i)
//                 XL(i,j,k1-1) = 0;
//     if (j1 > 0)
//         for (int k = k1; k < k2; ++k)
//             for (int i = i1; i < i2; ++i)
//                 XL(i,j1-1,k) = 0;
//     if (i1 > 0)
//         for (int k = k1; k < k2; ++k)
//             for (int j = j1; j < j2; ++j)
//                 XL(i1-1,j,k) = 0;

//     if (k2 < n)
//         for (int j = j1; j < j2; ++j)
//             for (int i = i1; i < i2; ++i)
//                 XL(i,j,k2) = 0;
//     if (j2 < n)
//         for (int k = k1; k < k2; ++k)
//             for (int i = i1; i < i2; ++i)
//                 XL(i,j2,k) = 0;
//     if (i2 < n)
//         for (int k = k1; k < k2; ++k)
//             for (int j = j1; j < j2; ++j)
//                 XL(i2,j,k) = 0;

//     #undef XL
//     #undef X
// }

void schwarz_get(int n, int i1, int i2, int j1, int j2, int k1, int k2,
                 double* restrict x_local, double* restrict x)
{
    #define X(i,j,k) (x[((k)*n+(j))*n+(i)])
    #define XL(i,j,k) (x_local[((k)*n+(j))*n+(i)])

    // Copy the data from x to x_local
    #ifdef PARALLEL
        // Parallel version: Use OpenMP to parallelize the nested loops
        #pragma omp parallel for collapse(3)
    #endif
    for (int k = k1; k < k2; ++k)
        for (int j = j1; j < j2; ++j)
            for (int i = i1; i < i2; ++i)
                XL(i,j,k) = X(i,j,k);

    // Set boundary slices to zero
    if (k1 > 0) {
        #ifdef PARALLEL
            #pragma omp parallel for collapse(2)
        #endif
        for (int j = j1; j < j2; ++j)
            for (int i = i1; i < i2; ++i)
                XL(i,j,k1-1) = 0;
    }
    if (j1 > 0) {
        #ifdef PARALLEL
            #pragma omp parallel for collapse(2)
        #endif
        for (int k = k1; k < k2; ++k)
            for (int i = i1; i < i2; ++i)
                XL(i,j1-1,k) = 0;
    }
    if (i1 > 0) {
        #ifdef PARALLEL
            #pragma omp parallel for collapse(2)
        #endif
        for (int k = k1; k < k2; ++k)
            for (int j = j1; j < j2; ++j)
                XL(i1-1,j,k) = 0;
    }

    if (k2 < n) {
        #ifdef PARALLEL
            #pragma omp parallel for collapse(2)
        #endif
        for (int j = j1; j < j2; ++j)
            for (int i = i1; i < i2; ++i)
                XL(i,j,k2) = 0;
    }
    if (j2 < n) {
        #ifdef PARALLEL
            #pragma omp parallel for collapse(2)
        #endif
        for (int k = k1; k < k2; ++k)
            for (int i = i1; i < i2; ++i)
                XL(i,j2,k) = 0;
    }
    if (i2 < n) {
        #ifdef PARALLEL
            #pragma omp parallel for collapse(2)
        #endif
        for (int k = k1; k < k2; ++k)
            for (int j = j1; j < j2; ++j)
                XL(i2,j,k) = 0;
    }

    #undef XL
    #undef X
}


// void schwarz_add(int n, int i1, int i2, int j1, int j2, int k1, int k2, 
//                  double* restrict Ax_local, 
//                  double* restrict Ax)
// {
//     #define AX(i,j,k) (Ax[((k)*n+(j))*n+(i)])
//     #define AXL(i,j,k) (Ax_local[((k)*n+(j))*n+(i)])
//     for (int k = k1; k < k2; ++k)
//         for (int j = j1; j < j2; ++j)
//             for (int i = i1; i < i2; ++i)
//                 AX(i,j,k) += AXL(i,j,k);
//     #undef AXL
//     #undef AX
// }

void schwarz_add(int n, int i1, int i2, int j1, int j2, int k1, int k2,
                 double* restrict Ax_local, double* restrict Ax)
{
    #define AX(i,j,k)  (Ax[((k)*n + (j))*n + (i)])
    #define AXL(i,j,k) (Ax_local[((k)*n + (j))*n + (i)])

    #ifdef PARALLEL
        // Parallel version: Use OpenMP to parallelize the nested loops
        #pragma omp parallel for collapse(3)
    #endif
    for (int k = k1; k < k2; ++k)
        for (int j = j1; j < j2; ++j)
            for (int i = i1; i < i2; ++i)
                AX(i,j,k) += AXL(i,j,k);

    #undef AXL
    #undef AX
}


/*@T
 *
 * The [[pc_schwarz_poisson3d]] function applies a preconditioner by
 * combining independent SSOR updates for the bottom half plus an overlap
 * region (a slap $n/2+o/2$ nodes thick), then updating the top half plus
 * an overlap region (another $n/2+o/2$ node slab).  The same idea could
 * be applied to more regions, or to better approximate solvers.
 *@c*/

#ifdef PARALLEL
void parallel_memset(double *start_p, double value, int n) 
{ 
#pragma omp parallel 
    { 
        int n_threads = omp_get_num_threads(); 
        int id = omp_get_thread_num(); 
        int share = (n + n_threads - 1) / n_threads; 
        if (id * share < n && (id + 1) * share >= n) 
        { 
            int rem = n - id * share; 
            memset(&(start_p[id * share]), value, rem * sizeof(double)); 
        } 
        else if (id * share <= n) 
        { 
            memset(&(start_p[id * share]), value, share * sizeof(double)); 
        } 
    } 
} 
#endif






void pc_schwarz_poisson3d(int N, void* data, 
                          double* restrict Ax, 
                          double* restrict x)
{
    pc_schwarz_p3d_t* ssor_data = (pc_schwarz_p3d_t*) data;
    double* scratch = ssor_data->scratch;
    int n = ssor_data->n;
    int o = ssor_data->overlap/2;
    double w = ssor_data->omega;

    #define COPY_TIME_AS 14
    #define SCHWARZ_GET_TIME 15
    #define SSOR_TIME    16
    #define SCHWARZ_ADD_TIME 17
    #define SSOR_FORWARD_TIME 18
    #define SSOR_DIAG_TIME 19
    #define SSOR_BACKWARD_TIME 20

    int n1 = n/2+o;
    int n2 = n/2-o;

    tic(COPY_TIME_AS);
    #ifdef PARALLEL
    //#parallel_memset(&(Ax[n1 * n * n]), 0, (N - n1 * n * n));
    parallel_memset(Ax, 0, N);
    #else
    memset(Ax, 0, N*sizeof(double));
    #endif
    total_copy_time_as += toc(COPY_TIME_AS);

   

    tic(SCHWARZ_GET_TIME);  
    schwarz_get        (n, 0,n, 0,n, 0,n1, scratch, x);
    total_get_time += toc(SCHWARZ_GET_TIME);

    tic(SSOR_TIME);
    
    tic(SSOR_FORWARD_TIME);
    ssor_forward_sweep (n, 0,n, 0,n, 0,n1, scratch, w);
    partial_ssor_forward += toc(SSOR_FORWARD_TIME);
    
    tic(SSOR_DIAG_TIME);
    ssor_diag_sweep    (n, 0,n, 0,n, 0,n1, scratch, w);
    partial_ssor_diag += toc(SSOR_DIAG_TIME);

    tic(SSOR_BACKWARD_TIME);
    ssor_backward_sweep(n, 0,n, 0,n, 0,n1, scratch, w);
    partial_ssor_backward += toc(SSOR_BACKWARD_TIME);
    
    total_ssor_time += toc(SSOR_TIME);

    tic(SCHWARZ_ADD_TIME);
    schwarz_add        (n, 0,n, 0,n, 0,n1, scratch,Ax);
    total_add_time += toc(SCHWARZ_ADD_TIME);

    tic(SCHWARZ_GET_TIME);
    schwarz_get        (n, 0,n, 0,n, n2,n, scratch, x);
    total_get_time += toc(SCHWARZ_GET_TIME);

    tic(SSOR_TIME);
    tic(SSOR_FORWARD_TIME);
    ssor_forward_sweep (n, 0,n, 0,n, n2,n, scratch, w);
    partial_ssor_forward += toc(SSOR_FORWARD_TIME);

    tic(SSOR_DIAG_TIME);
    ssor_diag_sweep    (n, 0,n, 0,n, n2,n, scratch, w);
    partial_ssor_diag += toc(SSOR_DIAG_TIME);

    tic(SSOR_BACKWARD_TIME);
    ssor_backward_sweep(n, 0,n, 0,n, n2,n, scratch, w);
    partial_ssor_backward += toc(SSOR_BACKWARD_TIME);

    total_ssor_time += toc(SSOR_TIME);

    tic(SCHWARZ_ADD_TIME);
    schwarz_add        (n, 0,n, 0,n, n2,n, scratch,Ax);
    total_add_time += toc(SCHWARZ_ADD_TIME);
}

/*@T
 * \section{Forcing functions}
 *
 * The convergence of CG depends not only on the operator and the
 * preconditioner, but also on the right hand side.  If the error
 * is very high frequency, the convergence will appear relatively
 * fast.  Without a good preconditioner, it takes more iterations
 * to correct a smooth error.  In order to illustrate these behaviors,
 * we provide two right-hand sides: a vector with one nonzero 
 * (computed via [[setup_rhs0]]) and a vector corresponding to a smooth
 * product of quadratics in each coordinate direction ([[setup_rhs1]]).
 *@c*/
void setup_rhs0(int n, double* b)
{
    int N = n*n*n;
    memset(b, 0, N*sizeof(double));
    b[0] = 1;
}

void setup_rhs1(int n, double* b)
{
    int N = n*n*n;
    memset(b, 0, N*sizeof(double));
    for (int i = 0; i < n; ++i) {
        double x = 1.0*(i+1)/(n+1);
        for (int j = 0; j < n; ++j) {
            double y = 1.0*(i+1)/(n+1);
            for (int k = 0; k < n; ++k) {
                double z = 1.0*(i+1)/(n+1);
                b[(k*n+j)*n+i] = x*(1-x) * y*(1-y) * z*(1-z);
            }
        }
    }
}


void setup_rhs2(int n, double* b)
{
    int N = n * n * n;
    for (int i = 0; i < N; ++i) {
        b[i] = 1.0;  // Constant function over the domain
    }
}


void setup_rhs3(int n, double* b)
{
    int idx = 0;
    for (int k = 0; k < n; ++k) {
        int kz = k % 2;
        for (int j = 0; j < n; ++j) {
            int jy = j % 2;
            for (int i = 0; i < n; ++i) {
                int ix = i % 2;
                int parity = (ix + jy + kz) % 2;
                b[idx++] = parity == 0 ? 1.0 : -1.0;  // Checkerboard pattern
            }
        }
    }
}


// Low-frequency sinusoidal RHS setup
void setup_rhs_low_freq(int n, double* b)
{
    int N = n * n * n;
    memset(b, 0, N * sizeof(double));

    const double pi = 3.14159265358979323846;
    const int m = 1.0;

    for (int i = 0; i < N; ++i)
    {
        // General low-frequency sinusoidal pattern
        b[i] = sin(2 * pi * m * i / N);
    }
}

// High-frequency sinusoidal RHS setup
void setup_rhs_high_freq(int n, double* b)
{
    int N = n * n * n;
    memset(b, 0, N * sizeof(double));

    const double pi = 3.14159265358979323846;
    const int m = 10;  // Frequency multiplier for high frequency
    const int A = 1.0;

    for (int i = 0; i < N; ++i)
    {
        // General high-frequency sinusoidal pattern
        b[i] = A *sin(2 * pi * m * i / N);
    }
}



/*@T
 * \section{The [[main]] event}
 *
 * The main driver is pretty simple: read the problem and solver parameters
 * using [[get_params]] and then run the preconditioned solve.
 *
 * I originally thought I'd write your own option routine, but I changed my
 * mind!
 *@c*/
int main(int argc, char** argv)
{
    solve_param_t params;
    if (get_params(argc, argv, &params))
        return -1;

    int n = params.n;
    int N = n*n*n;

    double* b = malloc(N*sizeof(double));
    double* x = malloc(N*sizeof(double));
    double* r = malloc(N*sizeof(double));
    memset(b, 0, N*sizeof(double));
    memset(x, 0, N*sizeof(double));
    memset(r, 0, N*sizeof(double));

    /* Set up right hand side */
    //setup_rhs_low_freq(n, b);

    //setup_rhs_high_freq(n, b);



    switch (params.rhs_type) {
        case 0:
            setup_rhs0(n, b);
            break;
        case 1:
            setup_rhs1(n, b);
            break;
        case 2:
            setup_rhs2(n, b);
            break;
        case 3:
            setup_rhs3(n, b);
            break;
        case 4:
            setup_rhs_low_freq(n, b);
            break;
        case 5:
            setup_rhs_high_freq(n, b);
            break;
        default:
            fprintf(stderr, "Invalid RHS type in params. Choose between 0 and 5.\n");
            free(b);
            free(x);
        
            return 1;
    }



    //  setup_rhs0(n,b);

    /* Solve via PCG */
    int maxit   = params.maxit;
    double rtol = params.rtol;

    if (params.ptype == PC_SCHWARZ) {
        double* scratch = malloc(N*sizeof(double));
        pc_schwarz_p3d_t pcdata = {n, params.overlap, params.omega, scratch};
        pcg(N, pc_schwarz_poisson3d, &pcdata, mul_poisson3d, &n, x, b, maxit, rtol);
        free(scratch);
        printf("  Copy time: %g seconds\n", total_copy_time_as);
        printf("  Get time: %g seconds\n", total_get_time);
        printf("  SSOR time: %g seconds\n", total_ssor_time);
        printf("  Forward SSOR time: %g seconds\n", partial_ssor_forward);
        printf("  Diagonal SSOR time: %g seconds\n", partial_ssor_diag);
        printf("  Backward SSOR time: %g seconds\n", partial_ssor_backward);
    
        printf("  Add time: %g seconds\n", total_add_time);
    } else if (params.ptype == PC_SSOR) {
        pc_ssor_p3d_t ssor_data = {n, params.omega};
        pcg(N, pc_ssor_poisson3d, &ssor_data, mul_poisson3d, &n, x, b,  maxit, rtol);
        printf("  Copy time: %g seconds\n", total_copy_time_ssor);
        printf("  Forward sweep time:   %g seconds\n", total_forward_time);
        printf("  Diagonal sweep time:  %g seconds\n", total_diag_time);
        printf("  Backward sweep time:  %g seconds\n", total_backward_time);
    } else {
        pcg(N, pc_identity, NULL, mul_poisson3d, &n, x, b, maxit, rtol);
    }

    /* Check answer */
    mul_poisson3d(N, &n, r, x);
    double rnorm2 = 0;
    for (int i = 0; i < n; ++i) r[i] = b[i]-r[i];
    for (int i = 0; i < n; ++i) rnorm2 += r[i]*r[i];
    printf("rnorm = %g\n", sqrt(rnorm2));

    free(r);
    free(x);
    free(b);
}
