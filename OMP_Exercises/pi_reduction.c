// this is the solutons for the pi.c

// Compile with OpenMP 
// gcc -fopenmp pi_reduction.c -o pi_reduction
//
// execute the program after compile
// ./pi_sol

#include <stdio.h>
#include <omp.h>

#define OMP_NUM_THREADS 2

static long num_steps = 100000000;
double step;

void main() {
    int i;
    double x, pi, sum = 0.0;
    double start_time, run_time;

    step = 1.0 / (double)num_steps;

    int max_threads_to_loop = 128;

    for (int i_num_threads = 1; i_num_threads <= max_threads_to_loop; i_num_threads ++)
    {
        int max_threads = omp_get_max_threads();

        start_time = 0;
        start_time = omp_get_wtime();

        sum = 0;
        omp_set_num_threads(i_num_threads);
        // omp_set_num_threads(OMP_NUM_THREADS);

        //#pragma omp parallel: This pragma creates a parallel region, 
        // which means that the code block within this pragma will be 
        // executed by multiple threads simultaneously. The number of 
        // threads created is determined by the environment or the number 
        // set explicitly using omp_set_num_threads().

        // private(i, x): This clause declares private variables i and x. 
        // Private variables in OpenMP are specific to each thread, 
        // meaning each thread will have its own copy of these variables. 
        // In this context, each thread will have its independent i and x
        // variables for loop iteration and computation.

        // reduction(+:sum): The reduction clause is used in parallel 
        // computing to perform a reduction operation on a variable 
        // across all threads and store the result in a single variable.
        // In this case, sum is a reduction variable, and +: specifies 
        // that the reduction operation is summation.

        // possible reduction operators
        // + (initial 0)
        // * (initial 1)
        // - (initial 0)
        // min (largest positive num)
        // max (most negative num)
        // C/C++ only &, |, ^, &&, ||

        #pragma omp parallel private(i, x) reduction(+:sum)
        {
            int id = omp_get_thread_num();
            int num_threads = omp_get_num_threads();

            for (i = id + 1; i <= num_steps; i += num_threads) {
                x = (i - 0.5) * step;
                sum = sum + 4.0 / (1.0 + x * x);
            }
        }

        pi = step * sum;
        run_time = omp_get_wtime() - start_time;
        printf("\n pi with %i num threads is %lf in %lf seconds\n ", i_num_threads, pi, run_time);
    }
}

// collect data for multiple runs

// Run1
//  pi with 1 steps is 3.141593 in 1.952806 seconds
 
//  pi with 2 steps is 3.141593 in 0.834938 seconds
 
//  pi with 3 steps is 3.141593 in 0.493269 seconds
 
//  pi with 4 steps is 3.141593 in 0.372953 seconds
 
//  pi with 5 steps is 3.141593 in 0.294046 seconds
 
//  pi with 6 steps is 3.141593 in 0.240796 seconds
 
//  pi with 7 steps is 3.141593 in 0.215361 seconds
 
//  pi with 8 steps is 3.141593 in 0.200445 seconds

// Run2
// pi with 1 steps is 3.141593 in 1.856094 seconds
 
// pi with 2 steps is 3.141593 in 0.808347 seconds

// pi with 3 steps is 3.141593 in 0.452218 seconds

// pi with 4 steps is 3.141593 in 0.316549 seconds

// pi with 5 steps is 3.141593 in 0.253433 seconds

// pi with 6 steps is 3.141593 in 0.224982 seconds

// pi with 7 steps is 3.141593 in 0.214976 seconds

// pi with 8 steps is 3.141593 in 0.224022 seconds

// Run3
// pi with 1 steps is 3.141593 in 1.843879 seconds

// pi with 2 steps is 3.141593 in 0.819182 seconds

// pi with 3 steps is 3.141593 in 0.457144 seconds

// pi with 4 steps is 3.141593 in 0.355077 seconds

// pi with 5 steps is 3.141593 in 0.261044 seconds

// pi with 6 steps is 3.141593 in 0.225336 seconds

// pi with 7 steps is 3.141593 in 0.221622 seconds

// pi with 8 steps is 3.141593 in 0.219872 seconds

//Run4
// pi with 1 steps is 3.141593 in 1.859892 seconds

// pi with 2 steps is 3.141593 in 0.846166 seconds

// pi with 3 steps is 3.141593 in 0.472267 seconds

// pi with 4 steps is 3.141593 in 0.365122 seconds

// pi with 5 steps is 3.141593 in 0.265104 seconds

// pi with 6 steps is 3.141593 in 0.239265 seconds

// pi with 7 steps is 3.141593 in 0.262894 seconds

// pi with 8 steps is 3.141593 in 0.283924 seconds

