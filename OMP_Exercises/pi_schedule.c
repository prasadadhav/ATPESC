// this is the solutons for the pi.c

// Compile with OpenMP 
// gcc -fopenmp pi_schedule.c -o pi_schedule
//
// execute the program after compile
// ./pi_schedule

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

    int max_threads_to_loop = 8;

    for (int i_num_threads = 1; i_num_threads <= max_threads_to_loop; i_num_threads ++)
    {
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

        // take the for loop and divide it into chunks
        // dealing cards
        // schedule(static [,chunk])
        // schedule(dynamic [,chunk])
        // schedule(guided [,chunk])
        // schedule(auto)

        int chunk_size = num_steps/i_num_threads;

        #pragma omp parallel for schedule(static, chunk_size)
        for (i=1;i<= num_steps; i++)
        {
            x = (i-0.5)*step;
            sum = sum + 4.0/(1.0+x*x);
        }

        pi = step * sum;
        run_time = omp_get_wtime() - start_time;
        printf("\n pi with static %i num threads is %lf in %lf seconds\n ", i_num_threads, pi, run_time);
    }
}

// collect data for multiple runs

// Run1
//  pi with 1 num threads is 3.141593 in 1.262194 seconds
 
//  pi with 2 num threads is 1.417284 in 0.977253 seconds
 
//  pi with 3 num threads is 0.907726 in 0.587722 seconds
 
//  pi with 4 num threads is 0.851757 in 0.612569 seconds
 
//  pi with 5 num threads is 0.707447 in 0.608829 seconds
 
//  pi with 6 num threads is 0.644526 in 0.569836 seconds
 
//  pi with 7 num threads is 0.554170 in 0.509664 seconds
 
//  pi with 8 num threads is 0.447761 in 0.474619 seconds

// Run2
//  pi with 1 num threads is 3.141593 in 1.244303 seconds
 
//  pi with 2 num threads is 1.664449 in 1.106708 seconds
 
//  pi with 3 num threads is 1.579038 in 0.833753 seconds
 
//  pi with 4 num threads is 1.254701 in 0.747386 seconds
 
//  pi with 5 num threads is 0.955560 in 0.670952 seconds
 
//  pi with 6 num threads is 0.610074 in 0.575933 seconds
 
//  pi with 7 num threads is 0.485135 in 0.530526 seconds
 
//  pi with 8 num threads is 0.470175 in 0.504267 seconds

// Run3
//  pi with 1 num threads is 3.141593 in 1.196835 seconds
 
//  pi with 2 num threads is 1.455485 in 0.968550 seconds
 
//  pi with 3 num threads is 0.975609 in 0.605920 seconds
 
//  pi with 4 num threads is 0.832349 in 0.601412 seconds
 
//  pi with 5 num threads is 0.757003 in 0.617381 seconds
 
//  pi with 6 num threads is 0.687734 in 0.596792 seconds
 
//  pi with 7 num threads is 0.478702 in 0.522682 seconds
 
//  pi with 8 num threads is 0.441965 in 0.524412 seconds

//Run4
//  pi with 1 num threads is 3.141593 in 1.254562 seconds
 
//  pi with 2 num threads is 1.463783 in 0.968440 seconds
 
//  pi with 3 num threads is 1.129517 in 0.561705 seconds
 
//  pi with 4 num threads is 1.186861 in 0.804447 seconds
 
//  pi with 5 num threads is 0.905734 in 0.710644 seconds
 
//  pi with 6 num threads is 0.780481 in 0.657655 seconds
 
//  pi with 7 num threads is 0.599548 in 0.618631 seconds
 
//  pi with 8 num threads is 0.421535 in 0.575315 seconds

