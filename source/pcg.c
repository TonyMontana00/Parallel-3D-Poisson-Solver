#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "timing.h"
#include "pcg.h"


double dot(int n, const double* x, const double* y)
{
    double result = 0;
    #ifdef PARALLEL
    #pragma omp parallel for reduction(+: result)
    #endif
    for (int i = 0; i < n; ++i)
        result += x[i]*y[i];
    return result;
}

/*@T
 * \section{Preconditioned CG}
 *
 * The PCG routine multiplies by $A$ and $M^{-1}$ through the
 * [[Mfun]] and [[Afun]] function pointers (taking a vector length $n$,
 * an opaque data object, an output buffer, and an input vector as arguments).
 * We also pass [[Mdata]] and [[Adata]] as a way of getting context into
 * the function calls%
 * \footnote{This could admittedly be more convenient in C}.
 * In addition, we take storage for the solution (set to an initial guess
 * on input) and the right hand side, as well as the maximum number of
 * iterations allowed and a relative error tolerance.
 *
 * The relative error tolerance is actually slightly subtle; we terminate
 * the iteration when
 * \[
 *   \frac{\|r^{(k)}\|_{M^{-1}}}
 *        {\|r^{(0)}\|_{M^{-1}}} < \mathrm{tol},
 * \]
 * where $\|\cdot\|_{M^{-1}}$ refers to the norm induced by the $M^{-1}$
 * inner product, i.e. $\|z\|_{M^{-1}}^2 = z^T M^{-1} z$.  This may or
 * may not be the norm anyone actually cares about... but it surely is
 * cheap to compute.
 *@c*/
double pcg(int n,
           mul_fun_t Mfun, void* Mdata,
           mul_fun_t Afun, void* Adata,
           double* restrict x,
           const double* restrict b,
           int maxit,
           double rtol)
{
    void parallel_copy(int n, double *Ax, double *x);

    double* r = malloc(n*sizeof(double));
    double* z = malloc(n*sizeof(double));
    double* q = malloc(n*sizeof(double));
    double* p = malloc(n*sizeof(double));

    double rho0     = 0;
    double rho      = 0;
    double rho_prev = 0;
    double rtol2 = rtol*rtol;
    int is_converged = 0;
    int step;

    double total_time = 0;
    double A_fun_time = 0;
    double res_time = 0;
    double M_fun_time = 0;
    double dot_time =   0;
    double search_dir_time = 0;
    double A_fun_2_time = 0;
    double dot_2_time = 0;
    double sol_time = 0;
    double res_time_2 = 0;

    double tot_MatVec_time = 0;
    double tot_dot_time = 0; 
    double tot_res_time = 0;
    double tot_prec_time = 0;
    double tot_search_dir_time = 0;
    double tot_sol_time = 0;




    #define TOTAL_TIME 0
    #define A_fun_TIME 1
    #define RES_TIME   2
    #define M_fun_TIME 3
    #define DOT_TIME   4
    #define SEARCH_DIR_TIME  5
    #define A_fun_2_TIME 6
    #define DOT_2_TIME 7
    #define SOL_TIME 8
    #define RES_TIME_2 9




    tic(TOTAL_TIME);

    /* Form residual */

    // A*x
    tic(A_fun_TIME);
    Afun(n, Adata, r, x);
    A_fun_time += toc(A_fun_TIME);

   tic(RES_TIME);
   
   #ifdef PARALLEL
    #pragma omp parallel for
    #endif
   for (int i = 0; i < n; ++i)
    {
        r[i] = b[i]-r[i];     
    }
    res_time += toc(RES_TIME);

    for (step = 0; step < maxit && !is_converged; ++step) {
        
        // M^-1*r
        tic(M_fun_TIME);
        Mfun(n, Mdata, z, r);
        M_fun_time += toc(M_fun_TIME);

        rho_prev = rho;

        // z^T*r
        tic(DOT_TIME);
        rho = dot(n, r, z);
        dot_time += toc(DOT_TIME);


        if (step == 0) {
            rho0 = rho;
            #ifdef PARALLEL
            parallel_copy(n, p, z);
            #else
            memcpy(p, z, n*sizeof(double));
            #endif
        } else {
            double beta = rho/rho_prev;
            
            tic(SEARCH_DIR_TIME);
            // p = z + beta*p
            #ifdef PARALLEL
            #pragma omp parallel for
            #endif
            for (int i = 0; i < n; ++i)
            {   
                p[i] = z[i] + beta*p[i];  
            }
            search_dir_time += toc(SEARCH_DIR_TIME);
        }
        
        // A*p
        tic(A_fun_2_TIME);
        Afun(n, Adata, q, p);
        A_fun_2_time += toc(A_fun_2_TIME);

        // p^T*A*p
        tic(DOT_2_TIME);
        double alpha = rho/dot(n, p, q);
        dot_2_time += toc(DOT_2_TIME);

        
        tic(SOL_TIME);
        // x = x + alpha*p
        #ifdef PARALLEL
        #pragma omp parallel for
        #endif
        for (int i = 0; i < n; ++i) 
        {
            x[i] += alpha*p[i];
        }
        sol_time += toc(SOL_TIME);

        tic(RES_TIME_2);
        // r = r - alpha*q
        #ifdef PARALLEL
        #pragma omp parallel for
        #endif
        for (int i = 0; i < n; ++i) 
        {     
            r[i] -= alpha*q[i];
        }
        res_time_2 += toc(RES_TIME_2);
        
        is_converged = (rho/rho0 < rtol2);
    }

    total_time = toc(TOTAL_TIME);

    tot_MatVec_time = A_fun_time + A_fun_2_time;
    tot_dot_time = dot_time + dot_2_time;
    tot_prec_time = M_fun_time;
    tot_res_time = res_time + res_time_2;
    tot_search_dir_time = search_dir_time;
    tot_sol_time = sol_time;
    

    printf("%d steps, residual reduction %g (%s tol %g); time %g\n",
        step, sqrt(rho/rho0), is_converged ? "<=" : ">", rtol, total_time);

    //printf("Total time: %g\n", total_time);
    // printf("A_fun time: %g\n", A_fun_time);
    // printf("Residual time: %g\n", res_time);
    // printf("M_fun time: %g\n", M_fun_time);
    // printf("Dot time: %g\n", dot_time);
    // printf("Search dir time: %g\n", search_dir_time);
    // printf("A_fun_2 time: %g\n", A_fun_2_time);
    // printf("Dot_2 time: %g\n", dot_2_time);
    // printf("Sol time: %g\n", sol_time);
    // printf("Residual time 2: %g\n", res_time_2);

    printf("Total MatVec time: %g\n", tot_MatVec_time);
    printf("Total dot time: %g\n", tot_dot_time);
    printf("Total prec time: %g\n", tot_prec_time);
    printf("Total residual time: %g\n", tot_res_time);
    printf("Total search dir time: %g\n", tot_search_dir_time);
    printf("Total sol time: %g\n", tot_sol_time);




    free(p);
    free(q);
    free(z);
    free(r);

    return rho/rho0;
}
