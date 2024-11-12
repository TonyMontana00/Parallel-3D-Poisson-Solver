#include "timing.h"
#include <time.h>
#include <stdio.h>

/* If the clock macro is defined, we use the POSIX clock_gettime,
 * and CLOCK defines which timer should be used.  Otherwise, we use
 * the system clock() command.
 */


#ifdef CLOCK
static struct timespec watches[NWATCHES];
#else
static clock_t watches[NWATCHES];
#endif


void tic(int watch)
{
#ifdef CLOCK
    clock_gettime(CLOCK, watches+watch);
#else
    watches[watch] = clock();
#endif
}


double toc(int watch)
{
    double elapsed;
#ifdef CLOCK
    struct timespec now;
    clock_gettime(CLOCK, &now);
    elapsed = now.tv_nsec - (double) watches[watch].tv_nsec;
    elapsed *= 1.0E-9L;
    elapsed += now.tv_sec - (double) watches[watch].tv_sec;
#else
    clock_t now = clock();
    elapsed = (double) (now-watches[watch])/CLOCKS_PER_SEC;
#endif    
    return elapsed;
}

// Global timing variables SSOR
double total_copy_time_ssor = 0.0;
double total_forward_time = 0.0;
double total_diag_time = 0.0;
double total_backward_time = 0.0;


// Global timing variables AS
double total_copy_time_as = 0.0;
double total_get_time = 0.0;
double total_ssor_time = 0.0;
double total_add_time = 0.0;

double partial_ssor_forward = 0.0;
double partial_ssor_diag = 0.0;
double partial_ssor_backward = 0.0;


// Function to reset SSOR timings
void reset_ssor_timing() {
    total_forward_time = 0.0;
    total_diag_time = 0.0;
    total_backward_time = 0.0;
    total_copy_time_ssor = 0.0;
}

// Function to print SSOR timings
// void print_ssor_timing() {
//     printf("SSOR Preconditioning Timings:\n");
//     printf("  Forward sweep time:   %g seconds\n", total_forward_time);
//     printf("  Diagonal sweep time:  %g seconds\n", total_diag_time);
//     printf("  Backward sweep time:  %g seconds\n", total_backward_time);
//     printf("  Copy time:            %g seconds\n", total_copy_time);
// }
