// module purge
// module load gcc-4.6.2
// module load mvapich2-1.9a2/gnu-4.6.2

// ssh siman003@bolt.cs.ucr.edu
// ssh siman003@tardis.cs.ucr.edu
// S44802855s}

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
	unsigned long long  low_value=0, high_value=0, size =0, prime = 3, n = 10000000000, index = 0, primeSize=0;
	char *marked;
	char *primeCollection;
   
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

	timeBegin = MPI_Wtime();


	   
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
	high_value  = 3 + 2 * BLOCK_HIGH (id,p,n/2-1); 
	// high_value = 3 + 2 * BLOCK_HIGH(id,p,n/2-1);
	size = BLOCK_SIZE(id,p,n/2-1);
	primeSize = sqrt(n);

	unsigned long long proc0_size = (n/2-1)/p;
	if ((3 + 2 * proc0_size) < (int) sqrt((double) n)) 
	{
	    if (!id) 
	    printf ("Too many processes\n");
	    MPI_Finalize();
	    exit (1);
	}


	marked = (char *) malloc (size);//this is the array which will show if we have prime marked[2]=0 or not prime marked[4]=1 number
	primeCollection = (char *) malloc (primeSize);
	if (marked == NULL) 
	{
		printf ("Cannot allocate enough memory\n");
		MPI_Finalize();
		exit (1);
	}

	unsigned long long i;
	unsigned long long cacheSize=1<<19,last=0;	
	for (i = 0; i < size; i++) marked[i] = 0;
	for(i = 0; i < primeSize; i++) primeCollection[i] = 0;
	//we initializing marked array to 0 first and then if we find a non prime number we will change it to 1
		do{
		 prime = 3;
		 index = 0;
		 
			do {
				unsigned long long first, primeFirst;
				
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
	   			primeFirst = (prime * prime - 3)/2;


	   			for (i = first + last; i < MIN(last+cacheSize,size); i += prime) marked[i] = 1;//we are changing marked to 1 for non prime number
	   			for (i = primeFirst; i < primeSize; i += prime) primeCollection[i] = 1;

		
		 		while (primeCollection[++index]);
		 		prime = 2 * index + 3;

				} while (prime * prime <= high_value);
		
				low_value += cacheSize * 2;
		   		last += cacheSize;
				} while (last < size); // we will continue to looking until we have this while condition


	unsigned long long count = 0;
	unsigned long long global_count = 0;
    for (i = 0; i < size; i++)
    if (!marked[i]) count++;//here we are going to count the prime numbers if marked is zero that's a prime number
    

    MPI_Reduce (&count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
   //  elapsed_time += MPI_Wtime();

   // if (!id)
   //  {
   //  	printf ("%d primes are less than or equal to %d\n",global_count, n);
   //    	printf ("Total elapsed time: %10.6f\n", elapsed_time);
   //  }
	timeEnd = MPI_Wtime();
	theTime = timeEnd - timeBegin;
	if(id == 0)
	{
		printf("\nElapsed Time For Version 3: %10.6f seconds\n",theTime);
		fflush(stdout);
	}
	MPI_Reduce (&theTime, &theGlobalTime, 1, MPI_DOUBLE, MPI_MAX,0, MPI_COMM_WORLD);
	if (id==0)
	{
		global_count++;
		printf ("\n%llu primes are less than or equal to %llu\n",global_count, n);
		fflush(stdout);
		printf("\nElapsed Time For Version 1: %10.6f seconds\n",theTime);
		fflush(stdout);
	}
	free(marked);
	free(primeCollection);

   MPI_Finalize ();
   return 0;
}
