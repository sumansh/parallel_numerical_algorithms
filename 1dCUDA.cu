


// 1D Diffusion Problem CUDA Code




#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cuda_runtime_api.h>
#include <device_functions.h>



#define size 10
#define alpha 15*10e-6
#define dt 0.001
#define beta 0.0







//convergence test
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



//jacobi kernel code
__global__ void jacobi(double* a, double* b, double* x)
{
    	double sum = 0;	

	  int id = threadIdx.x + blockDim.x * blockIdx.x;

	  
	  for(int j=0;j<size;j++)
	  {
		if(id!=j)  
		sum = sum + a[id*size+j]*x[j];
	  }

	  x[id] = (b[id]-sum)/a[id*size+id];

}


//matrix generation kernel using finite difference
__global__ void matrix(double *a, double *b, double *T, double dx)
	{		 
		int i = threadIdx.x + blockDim.x * blockIdx.x;


		if(i==0)
		{
			a[0]=1;
			b[0]= T[0];
		}
		else if(i == (size-1))
		{
			a[size*size-1]= 1;
			b[i]= T[size-1];
		}
		else
		{
			a[i*size+i] = -((2*alpha)/(dx*dx)) - (1/dt);
			a[i*size+i-1] = alpha/(dx*dx);
			a[i*size+i+1] = alpha/(dx*dx);
			b[i] = -(T[i]/dt);
		}
	}











int main()
{
	
	//dynamic allocation
	double *h_T= (double*)malloc(sizeof(double)*size);
	double *T;
	
	double *A;

	double *B;

	double *h_To= (double*)malloc(sizeof(double)*size);
	h_T[size-1]=1000.0;
	
	double dx = 0.1;

	for(int i=0;i<size;i++)
	{
	  h_T[i]=0.0;
	}
	h_T[size-1]=1000.0;
    //Allocating memory on GPU
	cudaMalloc((void **)&A, sizeof(double)*size*size);
	cudaMalloc((void **)&B, sizeof(double)*size);
	cudaMalloc((void **)&T, sizeof(double)*size);


	


	

	//Iteration for updating x
	while(convergence_test(h_T,h_To))
	{   
		cudaMemcpy(T, h_T, sizeof(double)*size, cudaMemcpyHostToDevice);
		cudaMemcpy(h_To, T, sizeof(double)*size, cudaMemcpyDeviceToHost);

		matrix<<<1, size>>>(A, B, T, dx);
		cudaDeviceSynchronize();
		for(int i=0;i<50;i++)
		{
		jacobi<<<1,size>>>(A, B, T);
		cudaDeviceSynchronize();
		}
		
		cudaMemcpy(h_T, T, sizeof(double)*size, cudaMemcpyDeviceToHost);

		double sum=0.0;
		for(int i=0;i<size;i++)
		{
			sum = sum + dx*(h_T[i]-h_To[i])*beta;
		}

		//mapping of temperature values to new points
		for(int i=1;i<size-1;i++)
		{
			h_T[i] = h_T[i] + (( 2 * i +1)/2)* (sum/size) * (h_T[i+1]-h_T[i])*(1/dx);
		}

		//updating dx 
		dx = dx + sum/size;
	}

	cudaMemcpy(h_T, T, sizeof(double)*size, cudaMemcpyDeviceToHost);
	
	cudaFree(A);
	cudaFree(B);
	cudaFree(T);




for(int i=0;i<size;i++)
	printf("T[%d]= %lf ", i, h_T[i]);


printf("increase in length = %d\n", (dx*size - 1));
return 0;

}

