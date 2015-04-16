#include "mf.h"
#include <iostream>
#include <algorithm>
using namespace std;
void mf(int ny, int nx, int hy, int hx, const float* in, float* out)
{
        #pragma omp parallel for
        for (int y = 0; y < ny; ++y)
        {
                int t = max(0,y-hy);
                int b = min(ny-1,y+hy);
        	float temp,arr[(2*hx+1)*(2*hy+1)];
                #pragma omp parallel for
                for (int x = 0; x < nx; ++x)
                {
                        int l = max(0,x-hx);
                        int r = min(nx-1,x+hx);
                        int index=0;
                        for(int j=t;j<=b;j++)
                                for(int i=l;i<=r;i++)
                                        arr[index++] = in[i+j*nx];
                        std::nth_element(arr, arr + (index/2) ,arr+index);
                        temp = arr[index/2];
                        if(index%2!=0)
                                out[x + nx*y] = temp;
                        else
                        {
                                std::nth_element(arr, arr + (index/2) - 1 ,arr+index);
                                out[x + nx*y] = (temp + arr[(index/2)-1])/2.0;
                        }
                }
        }
}
