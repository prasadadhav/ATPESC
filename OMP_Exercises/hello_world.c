#include<omp.h>
#include<stdio.h>

void printValues(int ID, double A) 
{
    printf("ID: %d\n", ID);
    printf("A: %lf\n", A);
}

// Compile with OpenMP 
// gcc -fopenmp hello_world.c -o hello_world
//
// execute the program after compile
// ./hello_world

int main()
{
    // printf("hello \n");
    // printf("world \n");
    // printf("from OpenMP \n");

    double A[1000];

    // set the number of threads
    omp_set_num_threads(4);

    //need to check the number of threads allocated
    // this is for large scale systems

    // each block executes the same copy of code  within  the structured block

    #pragma omp parallel
    {
        // returns the rank of the thread
        int ID = omp_get_thread_num();
        printValues(ID,A[ID]);
    }
}