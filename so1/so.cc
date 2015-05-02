#include "so.h"
#include <algorithm>
#include <iostream>
#include <omp.h>
using namespace std;
void merge(data_t A[], data_t B[], int m, int n) 
 {
  	int i=0, j=0, k=0;
  	int size = m+n;
  	data_t *C = (data_t *)malloc(size*sizeof(data_t));
  	while (i < m && j < n)
	{
              if (A[i] <= B[j]) 
				C[k++] = A[i++];
              else
				C[k++] = B[j++];
              
  	}
  	if (i < m) 
		for (int p = i; p < m; p++,k++) 
			C[k] = A[p];
  	else 
		for (int p = j; p < n; p++,k++) 
			C[k] = B[p];
  	for( i=0; i<size; i++ ) 
		A[i] = C[i];
  	free(C);
  }
  
  /* Merges N sorted sub-sections of array a into final, fully sorted array a */
 void arraymerge(data_t *a, int size, int *index, int N)
  {
  	int i;
  	while ( N>1 ) 
	{
  	    for( i=0; i<N; i++ ) 
			index[i]=i*size/N; 
		index[N]=size;
		#pragma omp parallel for num_threads(i)
  	    for( i=0; i<N; i+=2 ) 
		{
			merge(a+index[i],a+index[i+1],index[i+1]-index[i],index[i+2]-index[i+1]);
  	    }
  	    N /= 2;
  	}
  }

void psort(int n, data_t* data) {
    // FIXME: make this more efficient with parallelism
	int threads = omp_get_max_threads();
	int *index = (int *)malloc((threads+1)*sizeof(int));
	for(int i=0; i<threads; i++) 
	{
		index[i]=i*n/threads;
		//cout<<index[i]<<" "<<endl;
	}
	index[threads]=n;
	#pragma omp parallel for 
	for(int i=0; i<threads; i++)
		 std::sort(data+index[i], data +index[i+1]);
	if( threads>1 ) 
		arraymerge(data,n,index,threads);
    
}
