


// 1D Diffusion Problem C Code




#include<bits/stdc++.h>
using namespace std;








// convergence test

int convergence_test(vector <double> &T , vector <double> &T_old)
{
	double sum = 0;
	for(int i=0;i<T.size();i++)
	{
		sum = sum + (T[i]-T_old[i])*(T[i]-T_old[i]);
	}
	sum = sum /T.size();
	if(sum < 10e-8)
	return 0;
	else 
	return 1;
}








void jacobi(vector < vector <double> > &a , vector <double> &b, vector <double> &x , int size)
{
vector <double> x_new(size,0) ;
double sum =0;
while(convergence_test(x,x_new))
{
	x_new = x;
	for(int i=0;i<size;i++)
	{
		sum = 0;
		for(int j=0;j<size;j++)
		{
			if(i!=j)
			sum =sum +a[i][j]*x[j];
		}
		x[i] = (b[i]- sum)/a[i][i];
	}
	
}


	
}



int update ( vector <double> &T, double dx,const double alpha,const double dt,const double beta, const int size)
	{

		vector < vector <double> > a(size,vector<double> (size,0));
		vector <double> b(size,0);
		vector <double> To = T;
		//cout << "T[5]= " << To[5] << endl;
		
		
		//initializing a matrix
		for(int i=1; i< size-1; i++)
			{
			a[i][i] = -((2*alpha)/(dx*dx)) - (1/dt);
			a[i][i-1] = alpha/(dx*dx);
			a[i][i+1] = alpha/(dx*dx);
			b[i] = -(T[i]/dt);
			}
		a[size-1][size-1] = 1;
		a[0][0] = 1;
		b[0] = T[0];
		b[size-1] = T[size-1];

		
		//solving for new temperatures after timestep
		jacobi(a, b, T, size);


		//calculation of total elongation
		double sum =0;
		for(int i=0;i<size;i++)
		{
			sum = sum + dx*(T[i]-To[i])*beta;
		}
		cout << "sum = " << sum << endl;
		vector <double> pl(size,0);
		pl[0]=dx;
		for(int i=0; i<size;i++)
		{
			pl[i] = pl[i-1] + dx * beta * (T[i]-To[i]);
		}
		for(int i=1; i<size;i++)
		{
			pl[i] = pl[i-1] + (pl[i]-pl[i-1])/2;
		}
		vector <double> pn(size,0);
		for(int i=1; i<size;i++)
		{
			pn[i] = (i+1)/2 *(dx+sum/size);
		}		
		


		//mapping of temperature values to new points
		for(int i=1;i<size-1;i++)
		{
			T[i] = T[i] + ((T[i]-T[i-1])*(pn[i]-pl[i]))/(pl[i]-pl[i-1]);
		}

		//updating dx 
		dx = dx + sum/size;
		//cout << "dx = " << dx << endl;
	}











int main()
{

	//declaring variables
	int size = 10;
	double alpha,dx,dt,beta;
	alpha = 15* 10e-6;
	dx = 1.0/size;
	dt = 0.001;
	beta = 12.0;
	//beta = 0;



	vector < vector <double> > a(size,vector<double> (size,0));
	vector <double> T(size,0.0);
	vector <double> T_old(size,0.0);
	T[size-1] = 1000.0;


	//loop for each timestep to update temperatures
	while(convergence_test(T,T_old))
	{
		T_old = T;
		update(T,dx,alpha,dt,beta,size);
				
	}




for(int i=0;i<size;i++)
{
	cout << "T[" << i   << "]= " << T[i] << endl;
}

cout << "increase in length =  " << (dx*size - 0.1*size) << endl;
return 0;

}

