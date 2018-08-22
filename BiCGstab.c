


// BiCG Stab C code




#include<stdio.h>
#include<stdlib.h>
#include<math.h>


double findL2norm(double *vector, int n);


int main()
{
	int n  = 5000;			// size of matrics
	float A[n][n];
	double b[n], x[n];
	float output[n];
	double r[n], rs[n], v[n], p[n], s[n], t[n];
	double rho, rho_prev, alpha, omega, beta;
	double norm_r, norm_b;
	double eps = 1e-8;
	int i, j, n_iter = 0;
	double sum, temp;

	
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
	for(i=0; i<n;i++)
		{
        for(j=0; j<n; j++)
		  {
		fscanf(amatrix, "%f", &A[i][j]);
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
	for(i=0; i<n;i++)
		{
		  fscanf(bvector, "%f", &output[i]);
		}		
	}

	fclose(bvector);
		 

	// Initialisation of other vectors
	
	for (i = 0; i < n; i++)
	{
		sum = 0.;
		for (j = 0; j < n; j++)
		{
			sum += A[i][j]*x[j];
		}
		b[i] = (double)output[i];
		r[i] = b[i] - sum;
		rs[i] = r[i];
		v[i] = 0.;
		p[i] = 0.;
		s[i] = 0.;
		t[i] = 0.;
		x[i] = 0.;
	}

	// Initialisation of other parameters

	rho = 1.;
	alpha = 1.;
	omega = 1.;
	beta = 0.;

	norm_r = findL2norm(r, n);
	norm_b = findL2norm(b, n);



	while (norm_r > eps*norm_b)
	{

		// upate rho_prev and compute rho

		rho_prev = rho;
		rho = 0.;
		
		for (i = 0; i < n; i++)
			rho += rs[i]*r[i];

		// compute beta
        if(rho_prev!=0 && omega!=0)
		    beta = (rho/rho_prev) * (alpha/omega);
	    else
			beta=0;
		// compute p

		for (i = 0; i < n; i++)
		{
			p[i] = r[i] + beta * (p[i] - omega*v[i]);
		}		

		// compute v

		for (i = 0; i < n; i++)
		{
			sum = 0.;
			for (j = 0; j < n; j++)
			{
				sum += A[i][j]*p[j];
			}
			v[i] = sum;
		}

		// compute alpha
		
		sum = 0.;
		for (i = 0; i < n; i++)
			sum += rs[i]*v[i];
        
		if(sum!=0)
		  alpha = rho/ sum;
        else
	      alpha=0;
		// compute s

		for (i = 0; i < n; i++)
		{
			s[i] = r[i] - alpha*v[i];
		}

		// compute t
		
		for (i = 0; i < n; i++)		
		{
			sum = 0.;
			for (j = 0; j < n; j++)
			{			
				sum += A[i][j]*s[j];
			}
			t[i] = sum;
		}

		// compute omega

		temp = 0.;
		for (i = 0; i < n; i++)
			temp += t[i]*s[i];		
	
		
		sum = 0.;
		for (i = 0; i < n; i++)
			sum += t[i]*t[i];
         
		 if(sum!=0)
		  omega = temp/ sum;
        else
	      omega=0;
		

		// compute x

		for (i = 0; i < n; i++)
			x[i] += alpha*p[i] + omega*s[i];

		// compute r 

		for (i = 0; i < n; i++)
			r[i] = s[i] - omega*t[i];
        
		
		norm_r = findL2norm(r, n);
		n_iter++;
	}

	for (i = 0; i < n; i++)		
		printf("x[%d] = %0.2f \n", i, x[i]);
	
	printf("No. of iterations = %d\n", n_iter);	
			
	return 0;
}


double findL2norm(double *vector, int n)
{
	double v_norm = 0.;
    double var = 0.;
		

	for (int i = 0; i < n; i++)
		var += vector[i]*vector[i];

	v_norm = sqrt(var);		
	return v_norm;
}

