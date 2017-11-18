// module purge
// module load gcc-4.6.2
// module load mvapich2-1.9a2/gnu-4.6.2



#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include  <time.h>
#include <stdlib.h>
// #include "MyMPI.h"


// n is the number of elements and id is the randk of each and p is the number of process
//for example we have 13 blocks and we want to divide between 5 process so 
//for process 0 we should start from 0*13/5=0 for process 1 we will get 1*13/5 =2 , process 2 will be 2*13/5 = 5 and ...
#define BLOCK_LOW(id,p,n)  ((unsigned long long)(id)*(unsigned long long)(n)/(unsigned long long)(p)) 
#define BLOCK_HIGH(id,p,n)  (BLOCK_LOW((id)+1,p,n)-1) 
#define BLOCK_SIZE(id,p,n)  (BLOCK_LOW((id)+1,p,n) - BLOCK_LOW(id,p,n))
#define BLOCK_OWNER(index,p,n) ((((unsigned long long p)*(unsigned long long index)+1)-1 ) / (unsigned long long n) )
#define MIN(a,b)  ((a)<(b)?(a):(b))


int main (int argc, char *argv[])
{
	double theTime=0,timeBegin=0,timeEnd=0,theGlobalTime=0;


	int p , id;
	unsigned long long low_value=0, high_value=0, size =0, prime = 3, n = 10000000000, index = 0;
	char *marked;
   
   	//Basically with MPI_Init we are forking our program so before this line every thing is the same for all the process. 
   	//MPI_init and MPI_Finalize are the part that we fork multiple process and exit them after we finished our job. These two parts are mandatory.
   	MPI_Init (&argc, &argv);
   	//now we are going to tell for each process what the identification of each process is with using MPI_Comm_rank,size 
   	//we are letting for all the process to know how many process participating in a given MPI communicator
	MPI_Comm_rank (MPI_COMM_WORLD, &id);//which one am I? 
	//MPI_COMM_WORLD contains all processes launched as a part of MPI programm
	MPI_Comm_size (MPI_COMM_WORLD, &p);//how many processes participate in a given MPI communicator
	MPI_Barrier(MPI_COMM_WORLD);
	// elapsed_time = -MPI_Wtime();

	double elapsed_time = -MPI_Wtime();


	   
	// if (argc != 2) 
	// {
	//     if (!id) 
	//     printf ("Command line: %s <m>\n", argv[0]);
	//     MPI_Finalize(); 
	// 	exit (1);
	// }

	//atoi is changing the string to value so we are passing the argv[1] value to n. if we run our program ./test 123 then n = 123
	if(argc>=2) n = atoll(argv[1]);

	low_value = 3 + 2 * BLOCK_LOW(id,p,n/2-1);  
	// high_value = 3 + 2 * BLOCK_HIGH(id,p,n/2-1);
	size = BLOCK_SIZE(id,p,n/2-1);

	unsigned long long proc0_size = (n/2-1)/p;
	if ((3 + 2 * proc0_size) < (int) sqrt((double) n)) 
	{
	    if (!id) 
	    printf ("Too many processes\n");
	    MPI_Finalize();
	    exit (1);
	}


	marked = (char *) malloc (size);//this is the array which will show if we have prime marked[2]=0 or not prime marked[4]=1 number
	if (marked == NULL) 
	{
		printf ("Cannot allocate enough memory\n");
		MPI_Finalize();
		exit (1);
	}

	unsigned long long i;
	for (i = 0; i < size; i++) marked[i] = 0;
	//we initializing marked array to 0 first and then if we find a non prime number we will change it to 1
		if (!id) index = 0;
		 // prime = 2;
		 unsigned long long first;
			do {
				
	   			if (prime * prime > low_value)     
	      		first = (prime * prime - low_value)/2;
	   			else {
	      			if (!(low_value % prime)) first = 0;
	      			else 
	      				{
	      					first = prime - (low_value % prime);
	      					if((low_value + first)%2 == 0) first += prime;
							first /= 2;
						}
	   				}
	   	
	   			for (i = first; i < size; i += prime) marked[i] = 1;//we are changing marked to 1 for non prime number
	   			if (!id) 
	   			{
	      		while (marked[++index]);//findng the next prime
	      		prime = 2 * index + 3;
	   			}
	   			MPI_Bcast (&prime,  1, MPI_INT, 0, MPI_COMM_WORLD);
				} while (prime * prime <= n); // we will continue to looking until we have this while condition


	unsigned long long count = 0;
	unsigned long long global_count = 0;
    for (i = 0; i < size; i++)
    if (!marked[i]) count++;//here we are going to count the prime numbers if marked is zero that's a prime number
    

    MPI_Reduce (&count, &global_count, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
   //  elapsed_time += MPI_Wtime();

   // if (!id)
   //  {
   //  	printf ("%d primes are less than or equal to %d\n",global_count, n);
   //    	printf ("Total elapsed time: %10.6f\n", elapsed_time);
   //  }
	elapsed_time += MPI_Wtime();
	
	
	if (id==0)
	{
		global_count++;
		printf ("\n%llu primes are less than or equal to %llu\n and the time is %d\n",global_count, n);
		printf ("Total elapsed time: %10.6f\n", elapsed_time);
		fflush(stdout);
	}

   free(marked);
   MPI_Finalize ();
   return 0;
}
