#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "params.h"

static void print_usage()
{
    fprintf(stderr, 
            "cgp3d\n"
            "\t-h: print this message\n"
            "\t-n: mesh size (100)\n"
            "\t-M: maximum iteration count (200)\n"
            "\t-r: relative residual tolerance (1e-6)\n"
            "\t-p: preconditioner: as, ssor, or id (id)\n"
            "\t-w: SSOR relaxation parameter (1.9)\n"
            "\t-o: Schwarz method overlap (10)\n"
            "\t-t: RHS type: point, func, const, check, low, high (point)\n");
}

static void default_params(solve_param_t* params)
{
    params->n       = 100;
    params->maxit   = 200;
    params->rtol    = 1e-6;
    params->ptype   = PC_ID;
    params->omega   = 1.9;
    params->overlap = 10;
    params->rhs_type = RHS_VARYING_FUNCTION;  // Default RHS type
}

int get_params(int argc, char** argv, solve_param_t* params)
{
    extern char* optarg;
    const char* optstring = "hn:M:r:p:w:o:t:";
    int c;

    #define get_int_arg(c, field) \
        case c: params->field = atoi(optarg); break
    #define get_flt_arg(c, field) \
        case c: params->field = (float) atof(optarg); break

    default_params(params);
    while ((c = getopt(argc, argv, optstring)) != -1) {
        switch (c) {
        case 'h': 
            print_usage(); 
            return -1;
        case 'p':
            if (strcmp(optarg, "id") == 0) {
                params->ptype = PC_ID;
            } else if (strcmp(optarg, "ssor") == 0) {
                params->ptype = PC_SSOR;
            } else if (strcmp(optarg, "as") == 0) {
                params->ptype = PC_SCHWARZ;
            } else {
                fprintf(stderr, "Unknown preconditioner type: %s\n", optarg);
                return -1;
            }
            break;
        case 't':  // New option for RHS type
            if (strcmp(optarg, "point") == 0) {
                params->rhs_type = RHS_POINT_SOURCE;
            } else if (strcmp(optarg, "func") == 0) {
                params->rhs_type = RHS_VARYING_FUNCTION;
            } else if (strcmp(optarg, "const") == 0) {
                params->rhs_type = RHS_CONSTANT;
            } else if (strcmp(optarg, "check") == 0) {
                params->rhs_type = RHS_CHECKERBOARD;
            } else if (strcmp(optarg, "low") == 0) {
                params->rhs_type = RHS_LOW_FREQ;
            } else if (strcmp(optarg, "high") == 0) {
                params->rhs_type = RHS_HIGH_FREQ;
            } else {
                fprintf(stderr, "Unknown RHS type: %s\n", optarg);
                return -1;
            }
            break;
        get_int_arg('n', n);
        get_int_arg('M', maxit);
        get_flt_arg('r', rtol);
        get_flt_arg('w', omega);
        get_int_arg('o', overlap);
        default:
            fprintf(stderr, "Unknown option\n");
            return -1;
        }
    }
    return 0;
}
