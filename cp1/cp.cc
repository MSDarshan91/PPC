#include "cp.h"
#include<math.h>
#include<iostream>
using namespace std;
void correlate(int ny, int nx, const float* data, float* result) 
{
	double* inter = new double[ny*nx];
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
		//Finding the Standard Deviation
		for (int x = 0; x < nx; ++x) 
		{
			temp = data[x + y*nx] - mean;
			inter[x + y*nx] = temp;
			sd += (temp)*(temp);
		}
		sd= sqrt(sd);
		//Finding zero mean and unit variance
		for (int x = 0; x < nx; ++x) 
		{
			inter[x + y*nx] = inter[x + y*nx] / sd;
		}
    }
	//for(int i =0;i< nx*ny;i++)
		//cout<< inter[i]<<" ";
	//cout<<endl;
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
	//for(int i =0;i< ny*ny;i++)
	//	cout<< result[i]<<" ";	
	//cout<<endl;
}

