


// 2D heat conduction Problem CUDA code



#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cuda_runtime_api.h>
#include <device_functions.h>
#include <cuda.h>
#include <stdio.h>
#define size 30




//Matrix calculation Kernel code
__global__ void matrix(float* A)
{
	  int id = threadIdx.x + blockDim.x * blockIdx.x;
	  int tid = threadIdx.x;
	  int bid = blockIdx.x;

	  //for lowermost layer
	  if((bid<size) && (tid==bid))
	  	 {A[id] = 1.0;}
	  

	  //for uppermost layer
	  else if((bid>=(size*(size-1)))&& (tid==bid))
	    {  A[id]=1.0; }
	  
	  //for right boundary
	  else if((bid%size==0) && (tid==bid))
	       { A[id]=-4.0;
	         A[id+1]=2.0;
	         A[id+size]=1.0;
	         A[id-size]=1.0;
	       }
	  
	  //for leftmost boundary
	  else if(((bid+1)%size==0) && (tid==bid))
		   {
		    A[id]=-4.0;
		    A[id-1]=2.0;
		    A[id+size]=1.0;
		    A[id-size]=1.0;
	       }
	    
	  //for innermost points
	  else if(tid==bid)
		   {
		    A[id]=-4.0;
		    A[id-1]=1.0;
		    A[id+1]=1.0;
		    A[id+size]=1.0;
		    A[id-size]=1.0;
	       }
}


//jacobi kernel code
__global__ void jacobi(float* x_now,
	                   float* output,
	                   float* filter)
{
      __shared__ float temp[size*size];
      __shared__ float sum;	
      sum=0.0;

	  int id = threadIdx.x + blockDim.x * blockIdx.x;
	  int tid = threadIdx.x;
	  int bid = blockIdx.x;
	  
	  if(tid!=bid)
	       {temp[tid] = filter[id] * x_now[tid];}
	  else
	       {temp[tid]=0.0;}
      __syncthreads();
  
	for(int i=0; i<size*size;i++)
	     {sum = sum + temp[i];}

	if(tid == bid)
	{  x_now[bid] = ( output[bid] - sum ) / filter[id]; }
}



int main(int argc, char **argv)
{	
	int n=size*size*size*size;
	float *h_A= (float*)malloc(sizeof(float)*n);
	float *A;

	float *x_now, *output;
	float h_xnow[size*size] = {0.0};
	float h_output[size*size] = {0.0};

	for(int i=0;i<size;i++)
	{
	  h_output[size*size-i-1]=1.0;
	}

    //Allocating memory on GPU
	cudaMalloc((void **)&A, sizeof(float)*n);
	cudaMalloc((void **)&x_now, sizeof(float) * size*size);
	cudaMalloc((void **)&output, sizeof(float) * size*size);


	cudaMemcpy(A, h_A, sizeof(float)*n, cudaMemcpyHostToDevice);
	cudaMemcpy(x_now, h_xnow, sizeof(float) * size*size, cudaMemcpyHostToDevice);
	cudaMemcpy(output, h_output, sizeof(float)*size*size, cudaMemcpyHostToDevice);


	matrix<<<size*size, size*size>>>(A);
	cudaDeviceSynchronize();


	//Iteration for updating x
	for (int k = 1; k<30; k++)
	{   
		jacobi<<<size*size, size*size>>>(x_now, output, A);
		cudaDeviceSynchronize();
	}



	cudaMemcpy(h_A, A, sizeof(float)*n, cudaMemcpyDeviceToHost);
	cudaMemcpy(h_xnow, x_now, sizeof(float)*size*size, cudaMemcpyDeviceToHost);
	
	cudaFree(A);
	cudaFree(output);
	cudaFree(x_now);



      
      FILE *fp;
    
      /* open the file */
      fp = fopen("results.dat", "w");
      if (fp == NULL) {
         printf("I couldn't open resultsCUDA.dat for writing.\n");
         exit(0);
      }
    
      /* write to the file */
      for (int i=0; i<size*size; ++i)
	  {
		  fprintf(fp, "%f ", h_xnow[i]);
	  }
      /* close the file */
      fclose(fp);

	

	return 0;
}
