//IMortant change A and B they are not ranodm any more



/*
#PBS -l nodes=1:nogpu:ppn=16,walltime=00:30:00

JOB_PATH=Your program path

cd$pwd
./project

module load gcc-4.7.2
g++ -o project project.c -I/opt/lapack/include /opt/lapack/lib/liblapacke.a /opt/lapack/lib/liblapack.a /opt/lapack/lib/librefblas.a -lgfortran -lm


*/

// #include <iostream.h>
// #include "lapacke.h"
// #include "blas.h"
// using namespace std;

// #include <iostream>
// #include "lapacke.h"
// #include "blas.h"
// using namespace std;
// // #include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdint.h> 

#include <stdio.h>


#include <string.h>
// #include <time.h>
// #include <stdlib.h>
// Mhoi91aGt
// module load gcc-4.7.2
 // -I/usr/local/opt/lapack/include
// gcc -o run hw2.cpp -I/usr/local/opt/lapack/include /opt/lapack/lib/liblapacke.a /opt/lapack/lib/liblapack.a /opt/lapack/lib/librefblas.a -lgfortran -lm
// JOB_PATH=
// cd $JOB_PATH
// ./run
// qsub jobfile_1
// cat jobfile_1.ohear
// """

//  1000, 2000, 3000, 4000, 5000

#define BILLION 4096L

/////////////////////////////


void cache1(double *a, double *bb, double *c, int n, int end, int b, int ib)

{  

    int B =32;
    int i,j,k;
    int i1,j1,k1;
// printf("line 1here %i \n",end);
// fflush(stdout);
    for(i = end; i < n; i+=B)
    {

        for(j = end; j < n; j+=B)
        {
            for(k = ib; k< end ; k+= B)

            {
                for(i1 = i; i1< i+B; i1++)
                    {
                      
                        for(j1 = j; j1 < j+B; j1++)
                        {
                            register double r=c[i1*n+j1];
                            for(k1 = k; k1 < k+B; k1++)
                            {
                                r += a[i1*n + k1]*bb[k1*n + j1]; 

                                // printf("end = %i , n = % i, a[%i*n+%i]=  %f   b[%i*n + %i]= %f\n",end , n,i1,k1,a[i1*n+k1],k1,j1, bb[k1*n + j1]);
                            }   
                            c[i1*n+j1]=r;
                           
                        }

                    }

            }
        }
    }
}


// //////////////////////////

void printB(double* B){
    for(int k = 0; k < 16; k++){
        printf("A[%i] = %lf\n",k, B[k]);
    }
}

int diff(double *a1, double *a2, int n)
{
    double difference = fabs(a1[0] - a2[0]);
    int i;
    for(i = 0;i < n*n; i++)
        {
            if(fabs(a1[i] - a2[i]) > difference)
                {
                    difference = fabs(a1[i] - a2[i]);
                }
        }
    
    return difference;
}


void randMatrixGen(double *M,int n, int m)
{
    int i,j;
    srand(time(NULL));
    for(i = 0; i < n * m; i++)
        M[i] = 2.0 * rand()/RAND_MAX - 1.0;


}

void randMatrixGenB(double *M,int n)
{
    int i;
    for(i = 0; i < n; i++)
        M[i] = 2.0 * rand()/RAND_MAX - 1.0;

    for(i = 0; i < n; i++)
        M[i] = i;

}

void randMatrixGenpvt(int *M,int n)
{

    int i;
   
    for(i = 0; i < n; i++)
        M[i] = i;
    

}

void mydgetfr(int *pvt, double *A, int n)
{
    
    int i, t;
    int maxind; 
    double max;
    for(i = 0; i < n; i++)
    {
        maxind = i;
        max = fabs(A[i*n+i]);
        for(t = i + 1 ; t < n; t++)
        {
            if(fabs(A[t*n+i]) > max)
            {
                maxind = t;
                max = fabs(A[t*n+i]);
             
            }
        }
        if(max == 0)
        {
           
            printf("LUfactoration failed: coefficient matrix is singular\n");
            // return ;
        }
        else
        {
            if(maxind != i)
            {
                int temps = pvt[i];
                pvt[i] = pvt[maxind];
                pvt[maxind] = temps;

                int j;
                double tempv;
                //i think it should be double tempv[];
                for(j = 0; j< n; j++)
                {
                    tempv = A[i*n+j];
                    A[i*n+j] = A[maxind*n+j];
                    A[maxind*n+j] = tempv;
                }
                
            }
        }

        
        int j, k;
        for(j = i+1; j < n; j++)
        {
           A[j*n+i] = A[j*n+i]/A[i*n+i];


            for(k = i+1; k < n; k++)
            {
                A[j*n+k] = A[j*n+k] - A[j*n+i] * A[i*n+k];
              
            }
        }
    }

}
// """
// %forward substitution. Your task: transform this part to mydtrsm(). 
////////////////////'''''''''''
void mydtsrsm(int *pvt, double *A, double *B, int n)
{ 

    double y[n];

    y[0] = B[pvt[0]];
    int i,j;
   

    for(i = 1; i < n; i++)
    {    

        double result = 0;   
        for(j = 0; j < i; j++)
            result += y[j] * A[i*n+j];

        y[i] =B[pvt[i]] - result;
    
    }
    double x[n];  
    x[n-1]= y[n-1]/A[(n-1)*n+n-1]; 

    for(i = n-2;i > -1; i--)
    {
        double result = 0;
        for(j = i+1; j < n; j++)
        result += x[j] * A[i*n+j];

        x[i] = (y[i]- result)/A[i*n+i];
    }
  
}


//////////////////////////''''''''''''


int main()
{   
    int b = 32;



    struct timespec begin, end;
    uint64_t t0;
    

    int n;
    n = BILLION;
     
   
    double *A = (double *)calloc(n*n,sizeof(double));


    randMatrixGen(A, n, n);
    double *B = (double *)calloc(n,sizeof(double));
    // printB(A);
 
    int sizeofC = n;
    double *c0 = (double *)calloc(sizeofC*sizeofC,sizeof(double));


    int *pvt = (int *)calloc(n,sizeof(int));

    randMatrixGenB(B, n);
    randMatrixGenpvt(pvt,n);
    int i,j;

    int ib=0;
   
    
    int t;
    int maxind;
  
    if(b > (n * n))
    {
        printf("Your Block Size is too big for this matrix");

    } 
  

    double max;



  clock_gettime(CLOCK_MONOTONIC,&begin);



     for(ib = 0; ib < n;ib = ib+b)
    {

        int end = fmin(ib+b,n);
        for(i = ib; i < end; i++)
        {
            maxind = i;
            max = fabs(A[i*n+i]);
            for(t = i + 1 ; t < n; t++)
            {
                if(fabs(A[t*n+i]) > max)
                {
                    maxind = t;
                    max = fabs(A[t*n+i]);
             
                }
            }
            if(max == 0)
            {
                printf("LUfactoration failed: coefficient matrix is singular\n");
            // return ;
            }
            else
            {
                if(maxind != i)
                {
                    int temps = pvt[i];
                    pvt[i] = pvt[maxind];
                    pvt[maxind] = temps;

               
                    double tempv;
                //i think it should be double tempv[];
                    for(j = i; j< end; j++)
                    {
                        tempv = A[i*n+j];
                        A[i*n+j] = A[maxind*n+j];
                        A[maxind*n+j] = tempv;


                    }
                
                }
            }
        }
        int  k;
        for(i = ib;i < end; i++)
        {
            for(j = i+1; j < end; j++)
            {
               A[j*n+i] = A[j*n+i]/A[i*n+i];

               for(k = i+1; k < end; k++)
                {
                    A[j*n+k] = A[j*n+k] - A[j*n+i] * A[i*n+k];
                   
                }
            }

        }


        for(i = ib;i < end; i++)
        {
            for(j = end; j < n; j++)

            {

               A[j*n+i] = A[j*n+i]/A[i*n+i];


               for(k = i+1; k < end; k++)
                {
                    A[j*n+k] = A[j*n+k] - A[j*n+i] * A[i*n+k];

                }
            }

        }


                    
                
        for(i = ib; i < end; i++)
        {
            if(pvt[i] > i)
            {
                // printf("pvt[%i] = %i\n", i, pvt[i]);
                double tempv;
                for(j = end; j< n; j++)
                {
                        tempv = A[i*n+j];
                        A[i*n+j] = A[pvt[i]*n+j];
                        // printf("A[04] = %lf\n", A[3]);
                        A[pvt[i]*n+j] = tempv;
                      
                }
            }
        }   

        int s;
        // ib=0;
        for(j = ib+1; j < end; j++)
        {
            for(k = end; k < n ; k++)
            {
                double result = 0;
                for(s = ib; s < j; s++)
                {
                    result += A[j*n+s] * A[s*n+k];
                }

                A[j*n+k] = A[j*n+k] - result;
            }
        }
        

    cache1(A,A,c0,n, end, b, ib);
     for(i = end; i < n; i++)
    {
        for(j = end; j < n; j++)
        {
            A[i*n+j] -=c0[i*n+j];
           
        }
    }

  

}


mydtsrsm(pvt,A,B,n);
// end2 = clock();
clock_gettime(CLOCK_MONOTONIC,&end);
t0 = 1000000000L *(end.tv_sec - begin.tv_sec) + (end.tv_nsec - begin.tv_nsec);
printf ("time_spent1 runtime: %llu nanoseconds\n", (long long unsigned int) t0);



return 0;

}


