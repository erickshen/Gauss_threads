#include <iostream>
#include <cmath>
void MatrixOutput(int l, int m, int n, double* x)
{
	std::cout << std::endl;
    std::cout.setf(std::ios::scientific);
    for(int i = 0; i <= l-1; i++)
    {
        for(int j = 0; j <= m-1; j++)
            std::cout << x[n*i+j] << " " ;
        std::cout << std::endl;
    }
}
