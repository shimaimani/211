#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define BLOCK_LOW(id,p,n) 	  ((id)*(n)/(p))

#define BLOCK_HIGH(id,p,n)        ( BLOCK_LOW((id)+1,p,n)-1 ) 

#define BLOCK_SIZE(id,p,n)        (BLOCK_LOW( (id)+1, p, n) -  BLOCK_LOW( (id), p, n  ) )

#define BLOCK_OWNER(index,p,n)    ( ( ((p)*(index)+1)-1 ) / (n) )
#define MIN(a,b)  		  ((a)<(b)?(a):(b))

int main(int argc,char* argv[])
{
	int  coreNum, coreID, rc; 

	double theTime=0,timeBegin=0,timeEnd=0,theGlobalTime=0;
	unsigned long long N=10000000000,lowVal=0,highVal=0,prime=3,index=0,size=0;
	char* primeMarkArray;



	rc = MPI_Init(&argc,&argv);

	if (rc != MPI_SUCCESS) 
	{
		printf ("Error starting MPI program. Terminating...\n");
		MPI_Abort(MPI_COMM_WORLD, rc);
	}
	MPI_Comm_size(MPI_COMM_WORLD,&coreNum);
	MPI_Comm_rank(MPI_COMM_WORLD,&coreID);
	MPI_Barrier(MPI_COMM_WORLD);
	timeBegin = MPI_Wtime();



	if(argc>=2) N=atoll(argv[1]);


	lowVal  = 3 + 2*BLOCK_LOW (coreID,coreNum,N/2-1);

	size = BLOCK_SIZE(coreID,coreNum,N/2-1);

	unsigned long long P0_size=(N/2-1)/coreNum;
	if(2*P0_size+3<sqrt((double)N))
	{
		if(coreID==0) printf("Too many processors.\n");
		MPI_Finalize();
		exit(1);
	}


	primeMarkArray=new char[size];
	//Clear array
	for(unsigned long long i=0;i<size;i++) primeMarkArray[i]=0;
	//Commence search
	
	if (coreID==0) index = 0;
	//prime=3;
	do {
		unsigned long long first;
		if (prime * prime > lowVal)
		 first = (prime * prime - lowVal)/2;
		else {
		 if (lowVal % prime == 0) first = 0;
		 else 
			{
				first = (prime - (lowVal % prime));
				if((lowVal+first)%2==0) first+=prime;
				first/=2;
			}
		}
		for (unsigned long long i = first; i < size; i += prime) primeMarkArray[i] = 1;//prime; //For debugging
		if (coreID==0) {
		 while (primeMarkArray[++index]);
		 prime = 2*index + 3;
		}
		MPI_Bcast (&prime,  1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
	} while (prime * prime <= N);

	//End of search
	unsigned long long count = 0,gCount=0;
	for (unsigned long long i = 0; i < size; i++)
		if (!primeMarkArray[i]) 
		{
			count++;
			//printf("\nCore #%d [%d]: %d",coreID,i,2*i+lowVal);
			//fflush(stdout);
		}
//Debugging portion
//		else
//		{
//			printf("\nCore #%d [%d]: %d, divisible by %d",coreID,i,2*i+lowVal, primeMarkArray[i]);
//			fflush(stdout);
//		}
	MPI_Reduce (&count, &gCount, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM,0, MPI_COMM_WORLD);
	timeEnd=MPI_Wtime();
	theTime=timeEnd-timeBegin;
	MPI_Reduce (&theTime, &theGlobalTime, 1, MPI_DOUBLE, MPI_MAX,0, MPI_COMM_WORLD);
	if (coreID==0)
	{
		gCount++; //For 2
		printf ("\n%llu primes are less than or equal to %llu\n",gCount, N);
		fflush(stdout);
		printf("\nElapsed Time For Version 1: %10.6f seconds\n",theTime);
		fflush(stdout);
	}




	delete primeMarkArray;
	MPI_Finalize();
}
