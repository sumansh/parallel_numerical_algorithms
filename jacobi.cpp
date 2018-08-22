


// Jacobi Iteration C Code



#include <bits/stdc++.h>
#include<iostream>
using namespace std;







inline int convergence ( vector<double> &r, int size) 
		{
			 double con =0.0;
			 for( int i = 0; i< size; i++)
			 con= con + r[i]*r[i];

			 if(con < 10e-8 )
			 return 0;
			 else
			 return 1;
		}











void input ( vector < vector <double> > &a , vector <double> &b , int size)
		{
			
	//entering values to matrix
	FILE *amatrix;
	FILE *bvector;
	

	amatrix = fopen("amatrix.dat", "r");
      if (amatrix == NULL) {
         printf("I couldn't open amatrix.dat for reading.\n");
         exit(0);
      }
  //  while (!feof(amatrix) )              /* Check for the end of file*/
         /* To avoid memory corruption */
    {
	for(int i=0; i<size;i++)
		{
        for(int j=0; j<size; j++)
		  {
		
		fscanf(amatrix, "%lf", &a[i][j]);

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
 //   while (!feof(bvector)  )             /* Check for the end of file*/
         /* To avoid memory corruption */
    {
	for( int i=0; i<size;i++)
		{
		fscanf(bvector, "%lf", &b[i]);

		}		
	}
	fclose(bvector);


		
		}





int main()
{
	int size = 5;
	cout << "Enter size" << endl;
	cin >> size ;

	vector<vector<double> > a(size, vector<double>(size));
	vector<vector<double> > g(size, vector<double>(size));
	vector <double> b(size,0);
	vector <double> x(size,0.0);
	vector <double> x2(size,0.0);
	vector <double> r(size,0.0);
	vector <double> p(size,0.0);
	double  alpha;
	int count =0;
	double sum = 0;
	double temp = 0;
	input(a,b,size);

	FILE * (file);
	file=fopen("gmatrix.dat","w");

// write g matrix

	for(int i=0;i<size;i++)
	{
		for(int j=0;j<size;j++)
		{if(i!=j)
		g[i][j] = a[i][j]/a[i][i];
		else
		g[i][j] = 0.0;
		fprintf (file, "%lf ",g[i][j]);
		
		}
	fprintf (file, "\n");
	}
fclose (file);

	
	
		
			
	do
	{	
		
		for( int i=0; i< size ; i++)
		{ sum =0;
			for( int j= 0; j<size; j++)
			{
				if(i!=j)
				sum = sum + a[i][j] * x[j];
			}
		temp = x[i];
		x[i] = (b[i] - sum )/a[i][i];
		r[i] = x[i]-temp;
		cout <<  "x[i]" << x[i] << endl;
		}
	count++;
	} while ( convergence(r, size) != 0 );
for( int i=0; i<size; i++)
cout << "x[" << i << "] = " << x[i] << endl;	
cout << "No of iterations = " << count << endl;		
return 0;



}
