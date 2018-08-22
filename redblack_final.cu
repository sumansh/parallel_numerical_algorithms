


// Red Black Gauss Seidel  CUDA Code



#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cuda_runtime_api.h>
#include <device_functions.h>
#include <cuda.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#define size 1000
#define m 7




//Convergence check kernel


int convergence_test(double* T, double* T_old)
{
	double sum = 0.0;
	for(int i=0;i < size; i++)
		sum = sum + (T[i]-T_old[i])*(T[i]-T_old[i]);
	sum = sum/size;
	if(sum < 10e-8)
		return 0;
	else 
		return 1;
}



//jacobi kernel code r odd b even 
__global__ void jacobi_b(float* x_now, float* x_next, float* const output, const float* const coef, const int* const jcoef)
{	  
      float sum=0.0;
	  float aii = 1.0 ;

	  
	  int id = threadIdx.x + 1000 * blockIdx.x;
	  id = 2*id ;
	  
	  for(int i=0;i< m ;i++)
		 {
		if(jcoef[id*m+i]!= (size +1))
			{
				if( jcoef[id*m+i]!= id)
					sum += coef[id*m+i] * x_now[(int)jcoef[id*m+i]];
				else
					aii = coef[id*m+i];
			}
		else break ;
		 }
	   x_next[id] = (output[id]-sum)/aii;	  
}


__global__ void jacobi_r(float* x_now, float* x_next, float* const output, const float* const coef, const int* const jcoef)
{	  
      float sum=0.0;
	  float aii = 1.0 ;

	  
	  int id = threadIdx.x + 1000 * blockIdx.x;
	  id = 2*id +1;
	  
	  for(int i=0;i< m ;i++)
		 {
		if(jcoef[id*m+i]!= (size +1))
			{
				if( jcoef[id*m+i]!= id)
					sum += coef[id*m+i] * x_now[(int)jcoef[id*m+i]];
				else
					aii = coef[id*m+i];
			}
		else break ;
		 }
	   x_next[id] = (output[id]-sum)/aii;	  
}




int main(int argc, char **argv)
{	
    //initialize the variables
	
int si, sj;
	int h_ab = size;
	float *x_now, *x_next, *coef;
	int *jcoef, *ab ;
	float *output;
	float h_xnow[size] = { 0.0 };
	float h_output[size] ;
	float h_coef[size][m] = { 0.0 } ;
	int h_jcoef[size][m];
	for(int i=0; i<size;i++)
	for(int j=0; j<m;j++)
	 h_jcoef[i][j] =   size + 1  ;
	float temp = 0.0 ;
	int count = 0;
	clock_t start, end;
     double cpu_time_used;
	
	
	//entering values to matrix
	FILE *amatrix;
	FILE *bvector;
	amatrix = fopen("amatrix.dat", "r");
      if (amatrix == NULL) {
         printf("I couldn't open amatrix.dat for reading.\n");
         exit(0);
      }
	  
   
	for(int a=0; a<size;a++)
		{
			count = 0;
			for(int b=0; b<size; b++)
				{
					fscanf(amatrix, "%f", &temp);
					if(temp!= 0.0 )
					{
						h_coef[a][count] = temp;
						h_jcoef[a][count] = b;
						count++;
					}
				}							
		}
    fclose(amatrix);

	
	 bvector = fopen("bvector.dat", "r");
      if (bvector == NULL) 
	  {
         printf("I couldn't open bvector.dat for reading.\n");
         exit(0);
      }
	  	
	
   while (!feof(bvector)  )             /* Check for the end of file*/
         /* To avoid memory corruption */
    {
	for( int a=0; a<size;a++)
		{
		fscanf(bvector, "%f", &h_output[a]);
		}		
	}
	fclose(bvector);

	start = clock();	 




    //Allocating memory on GPU
	(cudaMalloc((void **)&x_now, sizeof(float) * size));
	(cudaMalloc((void **)&x_next, sizeof(float) * size));	
	(cudaMalloc((void **)&output, sizeof(float) * size));
	(cudaMalloc((void **)&coef, sizeof(float) * size* m ));
	(cudaMalloc((void **)&jcoef, sizeof(int) * size* m ));
	(cudaMalloc((void **)&ab, sizeof(int)));
	
	(cudaMemcpy( coef, h_coef, sizeof(float) * size* m , cudaMemcpyHostToDevice));
	(cudaMemcpy( jcoef, h_jcoef, sizeof(int) * size* m , cudaMemcpyHostToDevice));
	(cudaMemcpy(x_now, h_xnow, sizeof(float) * size, cudaMemcpyHostToDevice));
	(cudaMemcpy(x_next, h_xnow, sizeof(float) * size, cudaMemcpyHostToDevice));
	(cudaMemcpy(output, h_output, sizeof(float)*size, cudaMemcpyHostToDevice));

	
if(size%2!=0)
	{
		si = size/2;
		sj = si + 1; 
	}
else
	{
		si = size/2;
		sj = size/2;
	}



	//Iteration for updating x
    for (int k = 1; k<30; k++)
	{ 
	  if(k%2!=0)
			{
					jacobi_r<<< (int)ceil(si/1000.0), 1000 >>>(x_now, x_next, output, coef, jcoef);
					jacobi_b<<< (int)ceil(sj/1000.0), 1000 >>>(x_now, x_next, output, coef, jcoef);					
					cudaDeviceSynchronize();
			
			
			}
	  else
			{
					jacobi_b<<< (int)ceil(sj/1000.0), 1000 >>>(x_next, x_now, output, coef, jcoef);
					jacobi_r<<< (int)ceil(si/1000.0), 1000 >>>(x_next, x_now, output, coef, jcoef);					
					cudaDeviceSynchronize();			
			}


						
	}
    (cudaMemcpy(h_xnow, x_now, sizeof(float) * size, cudaMemcpyDeviceToHost));
	
	end = clock();

	cudaFree(x_now);
	(cudaFree(output));
	(cudaFree(coef));
	(cudaFree(jcoef));
	
	if(h_xnow[0]==10.0 && h_xnow[size-1]==5.0)
		printf("SUCCESS");
	printf("\n");
	
	printf("x[0]= %f      x[last]=%f\n",h_xnow[0], h_xnow[size-1]);
	for(int i =0; i<size; i++)
          printf("h_xnow[%d] = %f \n", i, h_xnow[i]);

cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

printf("Kernel took %f seconds to execute \n", cpu_time_used);
	return 0;
}
