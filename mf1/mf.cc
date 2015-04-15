#include "mf.h"
#include <iostream>     
#include <algorithm>    
#include <vector>  
using namespace std;

int compare (const void * a, const void * b)
{
  return ( *(float*)a - *(float*)b );
}

void mf(int ny, int nx, int hy, int hx, const float* in, float* out) {
    int t,b,l,r,index =0;
    float arr[(2*hx+1)*(2*hy+1)];
    for (int y = 0; y < ny; ++y) 	
    	for (int x = 0; x < nx; ++x) 
   	{
	 	t = max(0,y-hy);
		b = min(ny-1,y+hy);
 		l = max(0,x-hx);
		r = min(nx-1,x+hx);
		//n = (b-t+1)*(r-l+1);
//printf("%d %d %d %d %d %d\n",x,y,t,b,l,r);
		index=0;
 	        for(int j=t;j<=b;j++)
                   for(int i=l;i<=r;i++)
			arr[index++] = in[i+j*nx];
		 //std::nth_element(arr, arr + (index/2) ,arr+index);
		//qsort (arr, index , sizeof(float), compare);
		//out[x + nx*y] = arr[ index/2];
		std::sort(arr,arr+index); 
		if(index%2!=0)
			out[x + nx*y] = arr[index/2];
		else
			out[x + nx*y] = (arr[index/2]+arr[(index/2)-1])/2;
    	}
//	for (int y = 0; y < ny; ++y) 
//   			printf("%.1f ",out[x + nx*y] );
//		printf("\n");
//	}
//	printf("\n");
}

