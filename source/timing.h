#ifndef TIMING_H
#define TIMING_H

#define NWATCHES 8

void   tic(int watch);
double toc(int watch);

// Global timing variables

// SSOR
extern double total_forward_time;
extern double total_diag_time;
extern double total_backward_time;
extern double total_copy_time_ssor;

// AS
extern double total_copy_time_as;
extern double total_get_time;
extern double total_ssor_time;
extern double total_add_time;

extern double partial_ssor_forward;
extern double partial_ssor_diag;
extern double partial_ssor_backward;

// Functions to reset and print timings
void reset_ssor_timing();
//void print_ssor_timing();


#endif /* TIMING_H */




