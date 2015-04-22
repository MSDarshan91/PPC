#include "cp.h"
#include<math.h>
#include "vector.h"
#include<iostream>
using namespace std;

void correlate(int ny, int nx, const float* data, float* result)
{
  int vx = 0;
  if(nx%4)
    vx = (nx/4)+1;
  else
    vx = nx/4;  
  double4_t inter[vx*ny] ;
  for(int y=0; y<ny; y++)
    for(int x=0; x<nx; x++)
      inter[(y*vx)+(x/4)][x%4]= data[y*nx+x];
  int p = nx%4;
  if(p)
    for(int j=1; j<=ny; j++)
      for(int i=p; i<4; i++)
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
        #pragma omp parallel for schedule(static,1)
        for (int i = 0; i < ny; i++)
        {
                for (int j = i ; j < ny; j++)
                {
                        double4_t sum = double4_0;
      double res = 0.0;
                        for (int k = 0; k < vx; k++)
                                sum = sum + inter[i * vx + k] * inter[j*vx + k];
      for(int i =0;i<4;i++)
        res += sum[i];
                        result[i * ny + j] = res;
                }
        }
}


