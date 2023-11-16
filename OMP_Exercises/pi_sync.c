// this is the solutons for the pi.c

// Compile with OpenMP 
// gcc -fopenmp pi_sync.c -o pi_sync
//
// execute the program after compile
// ./pi_sync

#include <stdio.h>
#include <omp.h>

#define OMP_NUM_THREADS 2

static long num_steps = 100000000;
double step;

// SPMD algorithm strategy
// each thread has it's own copy of code
// unintended sharing of dtat causes race conditions
// to control the race conditions use synchronization to protect data conflicts
// change how data is accessed to minize the need for sync

// barrier sync
// mutual exclusion sync: one at a time

// constructs
// critical, atomic, barrier (these are simple)
// ordered, flush (low-level), locks

// this is an ugly solution
#define PAD 8 // assume 64byte L1 cache line size

/*
Exercise 

In previous ex. I used an array to create space for each thread to store its partial sum.

If array elements happen to share a cache line, this leads to false sharing

modify pi program to avoid false sharing using sync

*/

void main()
{
    int i, nthreads;
    double x, pi, sum[OMP_NUM_THREADS][PAD];
    double start_time, run_time;

    step = 1.0 / (double)num_steps;

    // here we promote the scalar sum from the pi.c to an array
    // it is of the same size as the number of threads we are going to use
    // Although we do not really have a truely shared memory
    // there re going to be problems when we try to write the sum into the cache
    // if some other thread is writing the values, the other thread has to wait, degrading performance
    // so for now we add a padding to sum. this padding will save avoid the problem of writing into one memory location

    omp_set_num_threads(OMP_NUM_THREADS);

    #pragma omp parallel
    {
        int i, id, nthrds;
        id = omp_get_thread_num();
        nthrds = omp_get_num_threads();

        if(id == 0) nthreads = nthrds;
        for(i = id; i < num_steps; i = i + nthrds)
        {
            sum[id][0] = 0.0;
            x = (i + 0.5) * step;
            sum[id][0] += 4.0/(1.0 + x*x);
        }
    }
    pi = 0.0;
    for(i = 0; i < nthreads; i++)pi +=sum[i][0] * step;
    printf("\n pi with %i num threads is %lf in %lf seconds\n ", nthreads, pi, run_time);
}