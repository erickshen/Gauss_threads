#include <iostream>
#include <cmath>
#include <fstream>
#include "Definitions.h"
#include <pthread.h>
#include <sys/time.h>
#include <time.h>
#include "sync.h"
double get_time(void)
{
	timeval t;
	gettimeofday(&t,0);
	return t.tv_sec+t.tv_usec/1000000.0;
}

struct pthrData
{
	int n;
	int p;
	int* result;
	int* Column;
        double* Row;
        double* massive;
        double* inverse;
	int i; // Номер треда
	double* time;
	int* index;
	double* max;
};
int Rationing(double* x, double* inverse, int n);
void PermuteColumn(double* x, int index, int step, double*Row, int n);
int FindMax(double* x, int step, int n);
void PermuteRow(double* x, int index, int step, double*Row, int n);
void Substitute(double* x, double* inverse, int step, int n);
int Back(double* x, double* inverse, int n);
int Rationing(double* x, double* inverse, int n);
double mach_eps(void);
static double eps;
void* GaussMethod(void* Bata)
{
	pthrData *Data;
	Data = ((pthrData*)Bata);
	int n = Data->n;
	int p = Data->p; //количество тредиков
	int* result = Data->result; 
	int* Column = Data->Column;
        double* Row = Data->Row;
        double* massive = Data->massive; //матричька
        double* inverse = Data->inverse;
	int thread_num = Data->i; // Номер треда
	double Norm = 1;
	double* time = Data->time;
	int* index = Data->index;
	double* max = Data->max;
	 eps = 10*mach_eps();
	synchronize(p);
	if(thread_num == 0)
	{
	
	Norm = MatrixNorm(massive, n);
        for(int i =0; i<= n*n-1; i++)
            massive[i] = massive[i]/Norm; //Нормировка accident...
	
	}
	synchronize(p); 
	if(thread_num == 0)
		*time = get_time(); //СТАРТ
        for(int step = 0; step <= n-1; step++)
    { synchronize(p);
		if(thread_num == 0) 
			for(int i = 0; i<p; i++)
				max[i] = 0.0;
	synchronize(p);
    	for(int i = thread_num; i <= n-1; i=i+p)
      		{ 
		if(i >= step)
		  for(int j = step; j < n; j++)
	  		{
			if(sqrt(massive[i*n + j]*massive[i*n + j]) > sqrt(max[thread_num]*max[thread_num]))
        		{
           			max[thread_num] = massive[i*n+j];
            		 	index[thread_num] = i*n+j;	
        		}	
						
      		} 
} 
	synchronize(p);
	if(thread_num == 0)
       	{
		int ind_max = n*n;
		double k = 0.0;
		
		for(int i = 0; i<p; i++)
			{
			if(sqrt(max[i]*max[i]) > sqrt(k*k))
				{ 
				   k = max[i];
				   ind_max = index[i]; 
				}
			else if(sqrt(max[i]*max[i]) == sqrt(k*k) && index[i] < ind_max)
				ind_max = index[i];
				
			}

		
        	if(sqrt(massive[ind_max]*massive[ind_max])< eps)
				{*result = -1;}
                Column[step] = ind_max;
                PermuteRow(massive, ind_max/n, step, Row, n);
                PermuteColumn(massive, ind_max%n, step, Row, n);
                PermuteRow(inverse, ind_max/n, step, Row, n);

	} 
	synchronize(p); 
	if(*result == -1)
		return NULL;
	
     		for(int l = thread_num; l < n-1; l = l+p)
		{	double a = massive[(l+1)*n + step]; 
			if(l+1 > step)
            		{for(int j = step; j <= n-1; j++)
                    	massive[(l+1)*n  + j] = massive[(l+1)*n + j] - massive[step*n + j]*a/massive[step*n+step];
			for(int j = 0; j <= n-1; j++)
                    	inverse[(l+1)*n + j] = inverse[(l+1)*n + j] - inverse[step*n + j]*a/massive[step*n+step];
                		 }
		}
	synchronize(p);
	
    }
	synchronize(p);
	
   for(int i = n-1; i >= 0; i--)
        {
for(int j = i-1-thread_num; j >= 0 && j < i; j = j-p){

            if(sqrt(massive[i*n+i]*massive[i*n+i])  > eps)
                for(int k = 0; k <= n-1; k++)
                    inverse[j*n + k] = inverse[j*n+k] - inverse[i*n+k]*massive[j*n+i]/massive[i*n+i];
            else  *result = -1;
              }
	synchronize(p);
}
            
	synchronize(p);
if(thread_num == 0) { int a = 0;
    a = Rationing(massive, inverse, n);
        if(a == -1)
           {*result = -1; return NULL;}
    for(int i = n-1; i >= 0; i--)
        PermuteRow(inverse, Column[i]%n, i, Row, n);

	        for(int i =0; i<= n*n-1; i++)
        inverse[i] = inverse[i]/Norm;
	*time = get_time() - *time; //K0NEЦ
   }
    return 0;

}
int FindMax(double* x, int step, int n)
{
    double k = 0;
    int index = 0;
    for(int i = n*step+step; i <= (n)*(n)-1; i = ((i+1)/n)*n + std::max(step, (i+1)%n))
      {  if(x[i]*x[i] > k*k)
        {
            k = x[i];
            index = i;
		
		
        }
      }
    return index;
}

void PermuteRow(double* x, int index, int step, double*Row, int n)
{
       for(int i = 0; i<= n-1; i++)
       {
           Row[i] = x[n*(index) + i];
           x[n*(index) + i] = x[n*step + i];
           x[n*step + i] = Row[i];
       }
}

 void PermuteColumn(double* x, int index, int step, double*Row, int n)
 {

     for(int j = 0; j <= n-1; j++)
            {
            Row[j] = x[index + n*j];
            x[index + n*j] = x[step + n*j];
            x[step  + n*j] = Row[j];

            }
}


int Rationing(double* x, double* inverse, int n)
{
    for(int i =0; i<=n-1; i++)
    {
        if(sqrt(x[i*n+i]*x[i*n+i]) > eps)
            for(int j = 0; j<=n-1; j++)
                inverse[i*n+j] = inverse[i*n+j]/x[i*n+i];
        else
            return -1;
    }
        return 0;
}
double mach_eps(void)
{ double eps = 1.0;
 while(1.0 + eps/2.0 > 1.0)
{ eps /= 2.0;}
return eps;
}
