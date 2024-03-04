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

// static long num_steps = 100000000;
// static long num_steps = 250000000;
static long num_steps = 500000000;
// static long num_steps = 1000000000;


double step;
int pi (int big_loop_itr)
{
	int i;
	double x, pi, sum = 0.0;
	double start_time, run_time;

	step = 1.0/(double) num_steps;

        	 
	start_time = omp_get_wtime();

	for (i=1;i<= num_steps; i++){
		x = (i+0.5)*step;
		sum = sum + 4.0/(1.0+x*x);
	}

	pi = step * sum;
	run_time = omp_get_wtime() - start_time;
	printf("Serial Execution, %lf, %lf, 1, %d\n", pi,run_time, big_loop_itr);
	// printf("\n Serial execution pi with %ld steps is %lf in %lf seconds\n ",num_steps,pi,run_time);
	
	// printf("\n----------------------------------------------------------\n");

	// ----------------------------------OMP Simple SMPD----------------------------------
	// Here we start the OMP parallelization
	// use #pragma omp parallel only
	//
	//  int omp_get_num_threads();  number of threads in team
	//  int omp_get_thread_num();	thread ID or Rank
	// 
	// get OMP stuff
	// we want to know how many threads were made available

	// These are the ones that we define 
	#define OMP_NUM_THREADS 128

	// declare
	int i_omp_simple, j_omp_simple;
	double pi_omp_simple, full_sum_omp_simple = 0.0;
	double start_time_omp_simple, run_time_omp_simple;
	double sum_omp_simple[OMP_NUM_THREADS];

	// loop over the threads
	for (int i_thrd = 1; i_thrd <= OMP_NUM_THREADS; i_thrd++) {
		omp_set_num_threads(i_thrd);
		pi_omp_simple = 0; 

		omp_set_num_threads(i_thrd);
      	full_sum_omp_simple=0.0;

		// we use the same step as serial execution
		start_time_omp_simple = omp_get_wtime();

		// we start the parallelisation
		#pragma omp parallel
		{
			int i_omp_simple_inner;
	  		int id_omp_simple = omp_get_thread_num();
	  		int numthreads_omp_simple = omp_get_num_threads();
	  		double x_omp_simple;

			sum_omp_simple[id_omp_simple] = 0.0;

        	// if (id_omp_simple == 0) 
            // printf(" num_threads = %d", numthreads_omp_simple);

			for (i_omp_simple_inner = id_omp_simple; i_omp_simple_inner < num_steps; i_omp_simple_inner += numthreads_omp_simple){
		  		x_omp_simple = (i_omp_simple_inner + 0.5) * step;
		  		sum_omp_simple[id_omp_simple] = sum_omp_simple[id_omp_simple] + 4.0/(1.0 + x_omp_simple*x_omp_simple);
	  		}
		}

		// sum up all the outputs from threads
		// calculate pi
		for (int id_itr = 0; id_itr <= i_thrd; id_itr++)
		{
			// pi_omp_simple = pi_omp_simple + step * sum_omp_simple[id_itr];
			full_sum_omp_simple +=  sum_omp_simple[id_itr];
		}

		pi_omp_simple = step * full_sum_omp_simple;


	// end timer
	// print timer
	run_time_omp_simple = omp_get_wtime() - start_time_omp_simple;
	// printf("\n Simple Parallel execution pi with %ld steps is %lf in %lf seconds\n ", num_steps, pi_omp_simple, run_time_omp_simple);
		printf("Simple Parallel execution, %lf, %lf, %d, %d\n", pi_omp_simple, run_time_omp_simple, i_thrd, big_loop_itr);

	}

	// the performance of this is very poor due to false sharing
	// printf("\n----------------------------------------------------------\n");
	// ---------------------------------- OMP Padded SMPD----------------------------------
	// When an int is promoted to an array, it  is shared on L! cache
	// This is false sharing
	// We will pad the array so that it exceeds the cache 
	// and is distributed
	// Padding the arrays to avoid false sharing
	//  
	#define PADDING 8

	// declare
	int i_omp_sim_pad, j_omp_sim_pad;
	double pi_omp_sim_pad, full_sum_omp_sim_pad = 0.0;
	double start_time_omp_sim_pad, run_time_omp_sim_pad;
	double sum_omp_sim_pad[OMP_NUM_THREADS][PADDING];

	// loop over the threads
	for (int i_thrd = 1; i_thrd <= OMP_NUM_THREADS; i_thrd++) {
		omp_set_num_threads(i_thrd);
		pi_omp_sim_pad = 0; 

		omp_set_num_threads(i_thrd);
      	full_sum_omp_sim_pad=0.0;

		// we use the same step as serial execution
		start_time_omp_sim_pad = omp_get_wtime();

		// we start the parallelisation
		#pragma omp parallel
		{
			int i_omp_sim_pad_inner;
	  		int id_omp_sim_pad = omp_get_thread_num();
	  		int numthreads_omp_sim_pad = omp_get_num_threads();
	  		double x_omp_sim_pad;

			sum_omp_sim_pad[id_omp_sim_pad][0] = 0.0;

        	// if (id_omp_sim_pad == 0) 
            // printf(" num_threads = %d", numthreads_omp_sim_pad);

			for (i_omp_sim_pad_inner = id_omp_sim_pad; i_omp_sim_pad_inner < num_steps; i_omp_sim_pad_inner += numthreads_omp_sim_pad){
		  		x_omp_sim_pad = (i_omp_sim_pad_inner + 0.5) * step;
		  		sum_omp_sim_pad[id_omp_sim_pad][0] = sum_omp_sim_pad[id_omp_sim_pad][0] + 4.0/(1.0 + x_omp_sim_pad*x_omp_sim_pad);
	  		}
		}

		// sum up all the outputs from threads
		// calculate pi
		for (int id_itr = 0; id_itr <= i_thrd; id_itr++)
		{
			// pi_omp_sim_pad = pi_omp_sim_pad + step * sum_omp_sim_pad[id_itr];
			full_sum_omp_sim_pad +=  sum_omp_sim_pad[id_itr][0];
		}

		pi_omp_sim_pad = step * full_sum_omp_sim_pad;


	// end timer
	// print timer
	run_time_omp_sim_pad = omp_get_wtime() - start_time_omp_sim_pad;
	// printf("\n Padded Parallel execution pi with %ld steps is %lf in %lf seconds\n ", num_steps, pi_omp_sim_pad, run_time_omp_sim_pad);
	printf("Padded Parallel Execution, %lf, %lf, %d, %d\n", pi_omp_sim_pad, run_time_omp_sim_pad, i_thrd, big_loop_itr);
	}

	// printf("\n----------------------------------------------------------\n");
	// ---------------------------------- OMP Critical ----------------------------------
	// Here we synchronize everything
	// use #pragma omp parallel only
	//
	//  int omp_get_num_threads();  number of threads in team
	//  int omp_get_thread_num();	thread ID or Rank
	// 
	// get OMP stuff
	// we want to know how many threads were made available

	// declare
	int i_omp_sim_crit;
	double pi_omp_sim_crit;
	double start_time_omp_sim_crit, run_time_omp_sim_crit;

	// start to loop over threads
	for (int i_thrd = 1; i_thrd <= OMP_NUM_THREADS; i_thrd++) {
		omp_set_num_threads(i_thrd);

		// reset vars
		start_time_omp_sim_crit = 0.0;
		run_time_omp_sim_crit = 0.0;

		// if I put pi = 0 here somehow the simple parallel exe pat, thrds = 5, returns pi = -nan
		// also this leads to the final pi values from each loop to wildy increased
		// as opposed to increasing by only 3.14
		// pi_omp_sim_crit = 0.0;

		// start time
		start_time_omp_sim_crit = omp_get_wtime();

		// OMP
		#pragma omp parallel
		{
			int i_omp_sim_crit_inner;
	  		int id_omp_sim_crit = omp_get_thread_num();
	  		int numthreads_omp_sim_crit = omp_get_num_threads();

			// declare local variables
			double sum_omp_sim_crit = 0.0, x_omp_sim_crit = 0.0;

			// #pragma omp single
			// if (id_omp_sim_crit == 0) 
            // printf(" num_threads = %d", numthreads_omp_sim_crit);

			for (i_omp_sim_crit_inner = 1; i_omp_sim_crit_inner <= num_steps; i_omp_sim_crit_inner++){
				x_omp_sim_crit = (i_omp_sim_crit_inner + 0.5)*step;
				sum_omp_sim_crit = sum_omp_sim_crit + 4.0/(1.0 + x_omp_sim_crit*x_omp_sim_crit);
			}

			// synochornize one at a time
			#pragma omp critical
			// setting pi = 0 does not work outside the critical zone
			// also same bug, when declaring pi = 0 here leads to simple parallel
			// thrd 5 to retunr pi -nan
			pi_omp_sim_crit = 0.0;
			pi_omp_sim_crit += step * sum_omp_sim_crit;
		}	
		run_time_omp_sim_crit = omp_get_wtime() - start_time_omp_sim_crit;

	// end timer
	// print timer
	run_time_omp_sim_crit = omp_get_wtime() - start_time_omp_sim_crit;
	// printf("\n Critical Parallel execution pi with %ld steps is %lf in %lf seconds\n ", num_steps, pi_omp_sim_crit, run_time_omp_sim_crit);
	printf("Critical Parallel Execution, %lf, %lf, %d, %d\n", pi_omp_sim_crit, run_time_omp_sim_crit, i_thrd, big_loop_itr);
	}

	// here somehow pi is not reset when starting a new for loop
	// for more number of threads
	// the bug goes away on 2nd execution

	// printf("\n----------------------------------------------------------\n");
	// ---------------------------------- OMP for loop ----------------------------------
	// last time we put the critical region inside OMP pragma
	// what if we put a critical zone inside the for loop computing sum
	// this is done with #pragma omp parallel for
	int i_omp_for;
	double pi_omp_for, sum_omp_for = 0.0;
	double start_time_omp_for, run_time_omp_for;

	// loop over the threads
	for (int i_thrd = 1; i_thrd <= OMP_NUM_THREADS; i_thrd++) {
		omp_set_num_threads(i_thrd);
		// start time
		start_time_omp_for = omp_get_wtime();

		// OMP
		#pragma omp parallel private(sum_omp_for)
		{
	  		int id_omp_for = omp_get_thread_num();
	  		int numthreads_omp_for = omp_get_num_threads();

			double x_omp_for;

			// #pragma omp single
			// if (id_omp_for == 0) 
            // printf(" num_threads = %d", numthreads_omp_for);

			#pragma omp for
			for (i_omp_for = 1; i_omp_for <= num_steps; i_omp_for++){
				x_omp_for = (i_omp_for + 0.5)*step;
				sum_omp_for = sum_omp_for + 4.0/(1.0 + x_omp_for * x_omp_for);
			}	
		}

		pi_omp_for = step * sum_omp_for;
		run_time_omp_for = omp_get_wtime() - start_time_omp_for;
		// printf("\n FOR parallel execution pi with %ld steps is %lf in %lf seconds\n ",num_steps,pi_omp_for, run_time_omp_for);
	printf("FOR parallel Execution, %lf, %lf, %d, %d\n", pi_omp_for, run_time_omp_for, i_thrd, big_loop_itr);
	}

	// printf("\n----------------------------------------------------------\n");

	// ---------------------------------- OMP for loop with reduction----------------------------------
	// In the previous problem there is a problem with sum
	// We are combining values into a single accumulation variable (ave) ... 
	// there is a true dependence between loop iterations that can’t be trivially removed.
	// This is a very common situation ... it is called a “reduction”
	// 
	// A local copy of each list variable is made and initialized depending on the “op” (e.g. 0 for “+”).
	// Updates occur on the local copy
	// Local copies are reduced into a single value and 
	// combined with the original global value.

	int i_omp_for_red;
	double pi_omp_for_red, sum_omp_for_red = 0.0;
	double start_time_omp_for_red, run_time_omp_for_red;

	// loop over the threads
	for (int i_thrd = 1; i_thrd <= OMP_NUM_THREADS; i_thrd++) {
		omp_set_num_threads(i_thrd);

		// declare local variables
		double sum_omp_for_red = 0.0, x_omp_for_red = 0.0;

		// start time
		start_time_omp_for_red = omp_get_wtime();

		// OMP
		#pragma omp parallel
		{
	  		int id_omp_for_red = omp_get_thread_num();
	  		int numthreads_omp_for_red = omp_get_num_threads();

			double x_omp_for_red;

			// #pragma omp single
			// if (id_omp_for_red == 0) 
            // printf(" num_threads = %d", numthreads_omp_for_red);

			#pragma omp for reduction(+:sum_omp_for_red)
			for (i_omp_for_red = 1; i_omp_for_red <= num_steps; i_omp_for_red++){
				x_omp_for_red = (i_omp_for_red + 0.5)*step;
				sum_omp_for_red = sum_omp_for_red + 4.0/(1.0 + x_omp_for_red * x_omp_for_red);
			}	
		}

		pi_omp_for_red = step * sum_omp_for_red;
		run_time_omp_for_red = omp_get_wtime() - start_time_omp_for_red;
		// printf("\n Reduced FOR parallel execution pi with %ld steps is %lf in %lf seconds\n ",num_steps,pi_omp_for_red, run_time_omp_for_red);
		printf("Reduced FOR parallel Execution, %lf, %lf, %d, %d\n", pi_omp_for_red, run_time_omp_for_red,  i_thrd, big_loop_itr);
	}

	// printf("\n----------------------------------------------------------\n");



	return 0;
}	  



int main()
{
	int num_max_iterations = 10;
	for (int big_loop_itr = 0; big_loop_itr <= num_max_iterations; big_loop_itr++)
	{
		printf("Label, PI value, Time [s], OMP Threads, Iteration No.\n");
		// printf("\n----------------------------------------------------------\n");
		// printf("\n----------------------------------------------------------\n");
		// printf("\n The iteration %d results are", big_loop_itr);
		pi(big_loop_itr);
		// printf("\n----------------------------------------------------------\n");
		// printf("\n----------------------------------------------------------\n");
	}

	return 0;
}
