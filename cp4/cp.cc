#include "cp.h"
#include<math.h>
#include "vector.h"
#include<iostream>
using namespace std;
#define VSIZE 4
void correlate(int ny, int nx, const float* data, float* result)
{
	for(int i =0;i<ny*ny;i++)
		result[i] = 0;
  int vx = 0;
  int pad= nx%4;
  if(nx%4)
    vx = (nx/4)+1;
  else
    vx = nx/4;  
  double4_t inter[vx*ny];
  for(int i=0; i<ny; i++)
    for(int j=0; j<nx; j++)
      inter[(i*vx)+(j/VSIZE)][j%VSIZE]= data[i*nx+j];
  if(pad)
    for(int j=1; j<=ny; j++)
      for(int i=pad; i<VSIZE; i++)
        inter[j*vx-1][i]=0.0;

  #pragma omp parallel for schedule(static,1)
        for (int y = 0; y < ny; ++y)
        {
                double mean = 0.0;
                double sd = 0.0;
                double4_t temp=double4_0;
    double4_t sd_temp=double4_0;
    double4_t mean_temp = double4_0;
                //Finding the mean
                for (int x = 0; x < vx; ++x)
                {
                        mean_temp += inter[x + y*vx];
                }
    for(int i= 0;i<4;i++)
      mean += mean_temp[i];
                mean = mean/nx;
                //Finding the Standard Deviation
    int t = nx%4;
    if(t!=0)
    {
      for(int i=t;i<4;i++)
        inter[(vx-1) + y*vx ][i] = mean;
    }
  //  for(int i=0;i<ny*vx;i++)
       //                 cout<<inter[i][0]<<" "<<inter[i][1]<<" "<<inter[i][2]<<" "<<inter[i][3]<<" ";
                for (int x = 0; x < vx; ++x)
                {
                        temp = inter[x + y*vx] - mean;
                        inter[x + y*vx] = temp;
                        sd_temp += (temp)*(temp);
                }
    for(int i =0;i<4;i++)
      sd += sd_temp[i];
                sd= sqrt(sd);
                //cout<<nx<<" "<<sd<<endl;
                //Finding zero mean and unit variance
                for (int x = 0; x < vx; ++x)
                {
                        inter[x + y*vx] = inter[x + y*vx] / sd;
                }
        }
		int s = 10;
		#pragma omp parallel for schedule(static,1)
     	for(int jj=0;jj<ny;jj+= s)
		{
			for(int ii=jj;ii<ny;ii+= s)
			{
				int t = ((jj+s)>ny?ny:(jj+s));
				int tt = ((ii+s)>ny?ny:(ii+s));
				for(int j = jj; j< t; j++)
				{
					for(int i = ii; i < tt; i++)
					{
						double4_t sum = double4_0;
						double res =0.0;
						for(int k = 0; k<vx; k++)
						{
							
							sum = sum + inter[i * vx + k] * inter[j*vx + k];
						}
						for(int i =0;i<4;i++)
						{
							res += sum[i];
						}
						result[j * ny + i] = res;
					}
				}
			}
		}
}
