


// Jacobi Iteration C Code



#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cuda_runtime_api.h>
#include <device_functions.h>
#include <cuda.h>
#include <stdio.h>
#include <math.h>
#define size 1728


//jacobi kernel code
__global__ void jacobi(float* x_now,
	                   float* const output,
	                   const float* const filter, int* d_flag)
{	  
      float sum=0.0;
      float temp=0.0;
	  int id = threadIdx.x + blockDim.x * blockIdx.x;
	  for(int i=0;i<size;i++)
		 {	if(i!=id)
			 {sum=sum+filter[id*size+i]*x_now[i];}
		 }
	  temp = (output[id]-sum)/filter[id*size+id];
	  if(abs(temp-x_now[id])>0 && (abs(temp-x_now[id]))<0.0000001)
            *d_flag=1;
	  x_now[id] = temp;

}




int main(int argc, char **argv)
{	
    //initialize the variables
	
	float *h_A = (float*)malloc(sizeof(float)*size*size);
    int flag=0; int *d_flag;
	float *x_now, *filter;
	float *output;
	float h_xnow[size] = { 0.0 };
	float h_output[size] ;
	float h_filter[size][size];
	
	
	
	//entering values to matrix
	FILE *amatrix;
	FILE *bvector;
	amatrix = fopen("amatrix.dat", "r");
      if (amatrix == NULL) {
         printf("I couldn't open amatrix.dat for reading.\n");
         exit(0);
      }
	  

    while (!feof(amatrix) )              /* Check for the end of file*/
         /* To avoid memory corruption */
    {
	for(int a=0; a<size;a++)
		{
        for(int b=0; b<size; b++)
		  {
		
		fscanf(amatrix, "%f", &h_filter[a][b]);
		
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
		 

		 

		 
	//dividing the 2D matrix into 1D matrix
	int ctr = 0;
	for (int i = 0; i<size; i++)
	{
		for (int j = 0; j<size; j++)
		{
		    h_A[ctr] = h_filter[i][j];
		    ctr++;
		}
	}



    //Allocating memory on GPU
	(cudaMalloc((void **)&x_now, sizeof(float) * size));
	(cudaMalloc((void **)&output, sizeof(float) * size));
	(cudaMalloc((void **)&filter, sizeof(float) * size*size));
	(cudaMalloc((void **)&d_flag, sizeof(int)));

	
	(cudaMemcpy(filter, h_A, sizeof(float) * size*size, cudaMemcpyHostToDevice));
	(cudaMemcpy(x_now, h_xnow, sizeof(float) * size, cudaMemcpyHostToDevice));
	(cudaMemcpy(output, h_output, sizeof(float)*size, cudaMemcpyHostToDevice));
	cudaMemcpy(d_flag, &flag, sizeof(int), cudaMemcpyHostToDevice);

	//Iteration for updating x
	int num=0;
    while(flag==0)
	{ 
		jacobi<<<ceil(size/1024.0), 1000>>>(x_now, output, filter, d_flag);
		cudaDeviceSynchronize();
		cudaMemcpy(&flag, d_flag, sizeof(int), cudaMemcpyDeviceToHost);
	}
    (cudaMemcpy(h_xnow, x_now, sizeof(float) * size, cudaMemcpyDeviceToHost));

	cudaFree(x_now);
	(cudaFree(output));
	(cudaFree(filter));

	printf("\n");

	return 0;
}
