


// BiCG Stab CUDA code




#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cuda_runtime_api.h>
#include <device_functions.h>
#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define n 4000



double findL2norm(double *vector)
{
	double v_norm = 0.;
    double var = 0.;
		

	for (int i = 0; i < n; i++)
		var += vector[i]*vector[i];

	v_norm = sqrt(var);		
	return v_norm;
}


__global__ void work1( double *d_x, double *d_AA, double *d_rs, double *d_b, double *d_r, double *d_p, double *d_v, double *d_s,
                      double *d_t, double *d_rho_prev, double *d_alpha, double *d_omega)
{
	int id = threadIdx.x + blockDim.x * blockIdx.x;


double rho = 0.0;
double beta;

for(int i=0;i<n;i++)
  rho += d_rs[i] * d_r[i]; 
 
if(*d_rho_prev!=0 && *d_omega!=0)
 { beta = (rho/(*d_rho_prev))*(*d_alpha/(*d_omega)); }
 else
 { beta=0.0; }
 
*d_rho_prev = rho;

d_p[id] = d_r[id] + beta * (d_p[id] - ((*d_omega) * d_v[id]));
 __syncthreads();
 
}


 
__global__ void multiply1(double *d_AA, double *d_p, double *d_v)
{
	int id = threadIdx.x + blockDim.x * blockIdx.x;


double sum=0.0;

for (int i = 0; i < n; ++i)
	sum = sum + d_AA[id*n + i]*d_p[i];
d_v[id]=sum; 
}




__global__ void work2( double *d_x, double *d_AA, double *d_rs, double *d_b, double *d_r, double *d_p, double *d_v, double *d_s,
                      double *d_t, double *d_rho_prev, double *d_alpha, double *d_omega)
{
	int id = threadIdx.x + blockDim.x * blockIdx.x;

 
d_t[id]=0.0 ;
__syncthreads();
 
//CALCULATE ALPHA

double sum=0.0;
for(int i=0;i<n;i++)
  {sum+=d_rs[i]*d_v[i];}
 
if(sum!=0)
{ *d_alpha = *d_rho_prev/sum; }
else
{ *d_alpha = 0; }
 
//CALCULATE S
 d_s[id] = d_r[id] - ((*d_alpha)*d_v[id]);

 
}
 
 
__global__ void multiply2(double *d_AA, double *d_s, double *d_t)
{
	int id = threadIdx.x + blockDim.x * blockIdx.x;


double sum=0.0;

for (int i = 0; i < n; ++i)
	sum = sum + d_AA[id*n + i]*d_s[i];
d_t[id]=sum; 
}
 
 
 
__global__ void work3( double *d_x, double *d_AA, double *d_rs, double *d_b, double *d_r, double *d_p, double *d_v, double *d_s,
                      double *d_t, double *d_rho_prev, double *d_alpha, double *d_omega)
{
int id = threadIdx.x + blockDim.x * blockIdx.x;
	



// COMPUTE OMEGA

double sum1=0.0 ;
for(int i=0;i<n;i++)
     sum1+=d_t[i]*d_s[i];

 
double sum2=0.0;
for(int i=0;i<n;i++)
     sum2+=d_t[i]*d_t[i];

 
 if (sum2!=0.0)
{ *d_omega = sum1/sum2; }
else
{*d_omega = 0;}
__syncthreads();


 
 //COMPUTE X
{d_x[id]+=(*d_alpha) * d_p[id] + (*d_omega) * d_s[id];}
 
 
 //COMPUTE R
{d_r[id] = d_s[id] - (*d_omega) * d_t[id];}

 
}








int main()
{
	double *AA = (double*)malloc(sizeof(double)*n*n);

	//double AA[n*n];
	double b[n];
	double x[n] = { 0.0 };
	double r[n], rs[n], v[n], p[n], s[n], t[n];
	double rho_prev=1.0 ; double alpha = 1.0 ; double omega = 1.0 ;
	double norm_r=0.0 ; double norm_b=0.0 ;
	double eps = 1e-8;
	
	
	//entering values to matrix
	FILE *amatrix;
	FILE *bvector;
	amatrix = fopen("amatrix.dat", "r");
      if (amatrix == NULL) {
         printf("I couldn't open amatrix.dat for reading.\n");
         exit(0);
      }

   
	for(int i=0; i<n;i++)
		{
        for(int j=0; j<n; j++)
		{
		fscanf(amatrix, "%lf", &AA[i*n+j]);
		}
		}	
    fclose(amatrix);	
	 
	 
	 bvector = fopen("bvector.dat", "r");
      if (bvector == NULL) 
	  {
         printf("I couldn't open bvector.dat for reading.\n");
         exit(0);
      }
	  	
 
    
	for( int i=0; i<n;i++)
		{
		fscanf(bvector, "%lf", &b[i]);
		}				
	fclose(bvector);
	
	
	// Initialisation of other vectors
	for (int i = 0; i < n; i++)
	{   
		r[i] = b[i];
		rs[i] = r[i];
		v[i] = 0.;
		p[i] = 0.;
		s[i] = 0.;
		t[i] = 0.;
	}


	double *d_AA, *d_b, *d_normB, *d_normR, *d_p, *d_r, *d_rho_prev, *d_rs, *d_s;
	double *d_t, *d_v, *d_x, *d_alpha, *d_omega;
	
	//Allocating memory on GPU
	cudaMalloc((void **)&d_AA, sizeof(double)*n*n);
	cudaMalloc((void **)&d_x, sizeof(double)*n);
	cudaMalloc((void **)&d_r, sizeof(double)*n);
	cudaMalloc((void **)&d_b, sizeof(double)*n);
	cudaMalloc((void **)&d_normR, sizeof(double));
	cudaMalloc((void **)&d_normB, sizeof(double));
	cudaMalloc((void **)&d_v, sizeof(double)*n);
	cudaMalloc((void **)&d_s, sizeof(double)*n);
	cudaMalloc((void **)&d_t, sizeof(double)*n);
	cudaMalloc((void **)&d_rs, sizeof(double)*n);
	cudaMalloc((void **)&d_p, sizeof(double)*n);
    cudaMalloc((void **)&d_rho_prev, sizeof(double));
	cudaMalloc((void **)&d_alpha, sizeof(double));
	cudaMalloc((void **)&d_omega, sizeof(double));
	
	cudaMemcpy(d_AA, AA, sizeof(double)*n*n, cudaMemcpyHostToDevice);
	cudaMemcpy(d_x, x, sizeof(double)*n, cudaMemcpyHostToDevice);
	cudaMemcpy(d_r, r, sizeof(double)*n, cudaMemcpyHostToDevice);
	cudaMemcpy(d_v, v, sizeof(double)*n, cudaMemcpyHostToDevice);
	cudaMemcpy(d_s, s, sizeof(double)*n, cudaMemcpyHostToDevice);
	cudaMemcpy(d_t, t, sizeof(double)*n, cudaMemcpyHostToDevice);
	cudaMemcpy(d_rs, rs, sizeof(double)*n, cudaMemcpyHostToDevice);
	cudaMemcpy(d_p, p, sizeof(double)*n, cudaMemcpyHostToDevice);
	cudaMemcpy(d_rho_prev, &rho_prev, sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_alpha, &alpha, sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_omega, &omega, sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_b, b, sizeof(double)*n, cudaMemcpyHostToDevice);
	

	norm_b = findL2norm(b);	
	norm_r= norm_b;
	

while ((norm_r) > eps*(norm_b))	
	{	
	work1<<< ceil(n/1000.0), 1000>>>(d_x, d_AA, d_rs, d_b, d_r, d_p, d_v, d_s, d_t, d_rho_prev, d_alpha, d_omega);
	cudaDeviceSynchronize();
	multiply1<<< ceil(n/1000.0), 1000>>>(d_AA, d_p, d_v);
	cudaDeviceSynchronize();
	work2<<<ceil(n/1000.0), 1000>>>(d_x, d_AA, d_rs, d_b, d_r, d_p, d_v, d_s, d_t, d_rho_prev, d_alpha, d_omega);
	cudaDeviceSynchronize();
	multiply2<<< ceil(n/1000.0), 1000>>>(d_AA, d_s, d_t);
	cudaDeviceSynchronize();
	work3<<< ceil(n/1000.0), 1000>>>(d_x, d_AA, d_rs, d_b, d_r, d_p, d_v, d_s, d_t, d_rho_prev, d_alpha, d_omega);
	cudaDeviceSynchronize();
	
	norm_r=findL2norm(r);

	}
	
	cudaMemcpy(x, d_x, sizeof(double)*n, cudaMemcpyDeviceToHost);
	
	for(int i=0;i<n;i++)
	{  printf("x[%d] = %0.2f \n", i, x[i] ); }
	
	cudaFree(d_AA);
	cudaFree(d_p);
	cudaFree(d_s);
	cudaFree(d_t);
	cudaFree(d_v);
	cudaFree(d_b);
	cudaFree(d_rs);
	cudaFree(d_r);
	cudaFree(d_x);
	cudaFree(d_rho_prev);
	cudaFree(d_alpha);
	cudaFree(d_omega);
	
	return 0;
}
	
	
	