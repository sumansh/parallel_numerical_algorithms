


// 2D heat conduction Problem C code




#define size 30
#include<stdio.h>
#include<stdlib.h>
#define n 900


int main()
{
	
	float B[size*size] = { 0.0 };
	float A[size*size][size*size] = { 0.0 } ;
	
	for( int i = 0;i  < size*size; i++)
	{  

	  if(i<size)
	  	{A[i][i] = 1.0;}
	  

	  else if(i>=size*(size-1))
	    {  A[i][i]=1.0; 
	       B[i]=1.0;}
	  
	  else if(i%size==0) 
	       { 
	       	 A[i][i]=-4.0;
	         A[i][i+1]=2.0;
	         A[i][i-size]=1.0;
	         A[i][i+size]=1.0;
	       }
	 
	  else if((i+1)%size==0)
		   {
		    A[i][i]=-4.0;
		    A[i][i-1]=2.0;
		    A[i][i+size]=1.0;
		    A[i][i-size]=1.0;
	       }
	    
	  else
		   {
		    A[i][i]=-4.0;
		    A[i][i-1]=1.0;
		    A[i][i+1]=1.0;
		    A[i][i+size]=1.0;
		    A[i][i-size]=1.0;
	       }
	}



	float x[n] = {0.0} ;
	float sum =0.;
	int temp = 0;

 
	for( int k =0; k<30; k++)
	{	
	for( int i =0; i<n; i++)
		{
		  sum = 0.;
		  for( int j =0; j<n; j++)
			  {
                if(i!=j)
			       sum = sum+ A[i][j]*x[j];
			  }

	       temp = x[i];	
	       x[i]=(B[i]- sum)/(A[i][i]) ;
		}	
	}


FILE *fp;
   
      /* open the file */
      fp = fopen("resultsC.dat", "w");
      if (fp == NULL) {
         printf("I couldn't open resultsC.dat for writing.\n");
         exit(0);
      }
   
      /* write to the file */
      for (int i=0; i<(size*size); ++i)
	  {
		  fprintf(fp, "%f ", x[i]);
	  }
      /* close the file */
      fclose(fp);
}
	
	
	

	
