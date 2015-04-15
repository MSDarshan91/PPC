#include "mf.h"
#include <iostream>     
#include <algorithm>    
#include <set>
using namespace std;

void mf(int ny, int nx, int hy, int hx, const float* in, float* out) {
    int t,b,l,r,index =0;
    
    float temp,arr[(2*hx+1)*(2*hy+1)];
    for (int y = 0; y < ny; ++y) 	
    	for (int x = 0; x < nx; ++x) 
   	{
	 	t = max(0,y-hy);
		b = min(ny-1,y+hy);
 		l = max(0,x-hx);
		r = min(nx-1,x+hx);

		index=0;
 	        for(int j=t;j<=b;j++)
                   for(int i=l;i<=r;i++)
			arr[index++] = in[i+j*nx];
		std::nth_element(arr, arr + (index/2) ,arr+index);
		temp = arr[index/2];
//		printf("ind/2--%.1f ",temp);
		if(index%2!=0)
			out[x + nx*y] = temp;
		else
		{
			std::nth_element(arr, arr + (index/2) - 1 ,arr+index);
//			printf("ind/2-1--%.1f",arr[(index/2)-1]);
			out[x + nx*y] = (temp + arr[(index/2)-1])/2.0;
		}
    	}
//	for (int y = 0; y < ny; ++y) 
//	{
//		for (int x = 0; x < nx; ++x) 
//	 			printf("%.1f ",out[x + nx*y] );
//		printf("\n");
//	}
//	printf("\n");
}

