#include <iostream>
#include <cmath>
#include <fstream>
#include <time.h>
#include <pthread.h>
#include "Definitions.h"
double get_time(void);
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

using namespace std;
int main(int argc, char* argv[])
{
if(argc < 5)
	{ cout << "Wrong number of arguments"; return -1;
}
    int result = 0;
  char* Filename;
	int n = atoi(argv[1]);
	int m = atoi(argv[2]); 
	int k = atoi(argv[3]); //function parameter
	int p = atoi(argv[4]); //num of threads
	if(m > n || m < 0 || n <= 0 || p >= n || p <=0)
		{cout << "Incorrect arguments" << endl; return -1;}
	if(k == 0)
		Filename = argv[5];
	



    double* massive = new double[n*n];
    double* inverse = new double[n*n];
    int* Column = new int[n];
    double* Row = new double[n];
    double* Max = new double[p];
    int* Index = new int[p];
    pthread_t* threads = new pthread_t[p];
    pthrData *Data = new pthrData[p];
    double time_taken = 0.0;
    double residual = 0;
    double elapsed = 0;
	

    if(k == 0)
    {
    int a = MatrixInputFile(n*n, massive, Filename);
        if(a == -1)
            {cout << "Incorrect arguments" << endl;
		return -1;}
    for(int i = 0; i <= n*n-1; i++)
            {

                if(i/n == i%n)
                    inverse[i] = 1;
                else inverse[i] = 0;
            }
    }


    if(k > 0 && k <5)
        for(int i = 0; i <= n*n-1; i++)
            {
                massive[i] = FunctionInput(k, n, i/n, i%n);
                if(i/n == i%n)
                    inverse[i] = 1;
                else inverse[i] = 0;
            }
	if(k < 0 || k > 4)
	{ cout << "Wrong argument" << endl; return -1;}
	cout << "The matrix is of the form: " << endl;
	MatrixOutput(m, m, n, massive);
	cout << endl;
	
	for(int i =0; i < p; i++)
	{
		Data[i].n = n;
		Data[i].p = p;
		Data[i].Row = Row;
		Data[i].Column = Column;
		Data[i].massive = massive;
		Data[i].inverse = inverse;
		Data[i].result = &result;
		Data[i].i = i;
		Data[i].time = &time_taken;
		Data[i].index = Index;
		Data[i].max = Max;
	}
	for(int i = 0; i < p; i++)
	{
		//cout << "creating thread number: " << i << endl;
		int j = pthread_create(&(threads[i]), NULL, GaussMethod, &(Data[i]));
		if(j) 
		{
			cout << "Error occured while ctreating threads. Process terminated" << endl;
			return 0;
		}
	}
	for(int i = 0; i < p; i++)
		pthread_join(threads[i], NULL);
 	cout<<"\n Time taken: " << time_taken << endl;

	if(*Data[0].result == -1)
		{  cout << "Singular\n";   printf("%s : residual = %e elapsed = %.2f k = %d n = %d m = %d p = %d\n", argv[0], residual, elapsed, k, n, m, p); return 0;}

        
    if(*Data[0].result == 0)
    {
        cout << "Inverse matrix is of the form: " << endl;

    MatrixOutput(m, m, n, inverse);
    }
   
    cout << "\n";
    if(k != 0)
        for(int i = 0; i <= n*n-1; i++)
                massive[i] = FunctionInput(k, n, i/n, i%n);

    if(k == 0)
    {
        MatrixInputFile(n*n, massive, Filename);
        }
    residual = ErrorNorm(inverse, massive, n);
    elapsed = time_taken;
    printf("%s : residual = %e elapsed = %.2f k = %d n = %d m = %d p = %d\n", argv[0], residual, elapsed, k, n, m, p);
	delete[] massive;
	delete[] inverse;
	delete[] Column;
        delete[] Row;
    return 0;


}

