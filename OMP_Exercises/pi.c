/*

This program will numerically compute the integral of

                  4/(1+x*x) 
				  
from 0 to 1.  The value of this integral is pi -- which 
is great since it gives us an easy way to check the answer.

The is the original sequential program.  It uses the timer
from the OpenMP runtime library

History: Written by Tim Mattson, 11/99.

*/

// Compile with OpenMP 
// gcc -fopenmp pi.c -o pi
//
// execute the program after compile
// ./pi

#include <stdio.h>
#include <omp.h>

#define OMP_NUM_THREADS 2

static long num_steps = 100000000;
double step;
void main ()
{
	int i;
	double x, pi, sum = 0.0;
	double start_time, run_time;

	step = 1.0/(double) num_steps;
	
	start_time = omp_get_wtime();

	for (i=1;i<= num_steps; i++){
		x = (i-0.5)*step;
		sum = sum + 4.0/(1.0+x*x);
	}

	pi = step * sum;
	run_time = omp_get_wtime() - start_time;
	printf("\n pi with %ld steps is %lf in %lf seconds\n ",num_steps,pi,run_time);

	// start parallelisations
	double x_omp, pi_omp, sum_omp[OMP_NUM_THREADS];
	double start_time_omp, run_time_omp;
	int nthreads;
	
	// set the number of threads
    omp_set_num_threads(OMP_NUM_THREADS);

	// #pragma omp parallel
	// {
	// 	// returns the rank of the thread
	// 	int id, nthrds;

	// 	// check number of threads is 
	// 	// same as what is asked for
	// 	nthrds = omp_get_num_threads();
	// 	id = omp_get_thread_num();

	// 	start_time_omp = omp_get_wtime();
	// 	if (id == 0) nthreads = nthrds;
	// 	for (i=id; sum_omp[id] = 0.0; i<= num_steps; i = i + nthrds)
	// 	{
	// 		sum_omp[id] = 0.0;
	// 		x_omp = (i+0.5)*step;
	// 		sum_omp[id] = sum_omp[id] + 4.0/(1.0+x_omp*x_omp);
	// 	}		
	// }
	// pi_omp = 0.0;
	// for(i = 0; i < nthreads; i++)
	// {
	// 	pi_omp += sum_omp[i] * step;
	// }
	// run_time_omp = omp_get_wtime() - start_time_omp;

	// printf("\n pi with %ld steps is %lf in %lf seconds\n ",num_steps,pi_omp,run_time_omp);

	// -----------------------------------------------------
	// another way of doing it is dividing the number of iterations
	// by the numer of threads

	int chunk_size = num_steps/OMP_NUM_THREADS;

	#pragma omp parallel
	{
		// return rank of the thread
		int id, nthrds;

		nthrds = omp_get_num_threads();
		id = omp_get_thread_num();

		

		start_time_omp = omp_get_wtime();

		for (i)
	}
}	  





