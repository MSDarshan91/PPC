#include "cp.h"
#include<math.h>
void correlate(int ny, int nx, const float* data, float* result)
{
	double* inter = new double[nx*ny];
	#pragma omp parallel for schedule(static,1)
        for (int y = 0; y < ny; ++y)
        {
                double mean = 0.0;
                double sd = 0.0;
                double temp;
                //Finding the mean
		for (int x = 0; x < nx; ++x)
                {
                        mean += data[x + y*nx];
                }
                mean = mean/nx;
                //cout<<mean<<endl;
                //Finding the Standard Deviation
                for (int x = 0; x < nx; ++x)
                {
                        temp = data[x + y*nx] - mean;
                        inter[x + y*nx] = temp;
                        sd += (temp)*(temp);
                }
                sd= sqrt(sd);
                //cout<<nx<<" "<<sd<<endl;
                //Finding zero mean and unit variance
                for (int x = 0; x < nx; ++x)
                {
                        inter[x + y*nx] = inter[x + y*nx] / sd;
                }
        }
	#pragma omp parallel for schedule(static,1)
        for (int i = 0; i < ny; i++)
        {
                for (int j = i ; j < ny; j++)
                {
                        double sum = 0.0;
                        for (int k = 0; k < nx; k++)
                                sum = sum + inter[i * nx + k] * inter[j * nx + k];
                        result[i * ny + j] = sum;
                }
        }
}
