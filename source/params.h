#ifndef PARAMS_H
#define PARAMS_H

/*@T
 * \section{Solver parameters}
 * 
 * The [[solve_param_t]] structure holds the parameters that
 * describe the simulation. These parameters are filled in
 * by the [[get_params]] function. Details of the parameters
 * are described elsewhere in the code.
 *@c*/
enum {                /* Types of preconditioners available: */
    PC_ID = 1,        /* 1. Identity                         */
    PC_SSOR = 2,      /* 2. SSOR (Symmetric Successive Over-Relaxation) */
    PC_SCHWARZ = 3    /* 3. Additive Schwarz (overlapping) */
};

enum {                /* Types of RHS configurations available: */
    RHS_POINT_SOURCE = 0,     /* 0. Point source at origin */
    RHS_VARYING_FUNCTION = 1, /* 1. Function varying with spatial coordinates */
    RHS_CONSTANT = 2,         /* 2. Constant function */
    RHS_CHECKERBOARD = 3,     /* 3. Checkerboard pattern */
    RHS_LOW_FREQ = 4,         /* 4. Low-frequency sinusoidal */
    RHS_HIGH_FREQ = 5         /* 5. High-frequency sinusoidal */
};

typedef struct solve_param_t {
    int    n;            /* Mesh size */
    int    maxit;        /* Maximum PCG iterations */
    double rtol;         /* Relative residual convergence tolerance */
    int    ptype;        /* Preconditioner type */
    double omega;        /* SSOR relaxation parameter */
    int    overlap;      /* Overlap size */
    int    rhs_type;     /* RHS configuration type */
} solve_param_t;

/* Function to retrieve parameters from command line arguments */
int get_params(int argc, char** argv, solve_param_t* params);

#endif /* PARAMS_H */
