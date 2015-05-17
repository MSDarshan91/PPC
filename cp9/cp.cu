#include "cp.h"
#include <math.h>
#include <cuda_runtime.h>
#include <chrono>
#include <iostream>
#include <unistd.h>
using namespace std;
using namespace std::chrono;
#define BLOCK_SIZE 16

#define CHECK_CUDA_ERROR(call) do { \
        cudaError_t result_ = (call); \
        if (result_ != cudaSuccess) { \
            fprintf(stderr, #call " failed: %s\n", \
                    cudaGetErrorString(result_)); \
            exit(1); \
        } \
    } while(0)
__global__ void var(float *input,float *output, int N, float mean)
	{
 
      int idx=threadIdx.x+(blockDim.x*blockIdx.x);
      if (idx < N) output[idx] = (input[idx]-mean)*(input[idx]-mean);
    }
__global__ void norm(float *input, int N,float mean,float sd)
	{
 
      int idx=threadIdx.x+(blockDim.x*blockIdx.x);
      if (idx < N) input[idx] =  (input[idx]-mean)/sd;
    }


__global__ void matrixMultiply(float * A, float * C, int ny, int nx) 
{
	if(blockIdx.x < blockIdx.y)
		return;
  __shared__ float ds_M[BLOCK_SIZE][BLOCK_SIZE];
  __shared__ float ds_N[BLOCK_SIZE][BLOCK_SIZE];
  int bx= blockIdx.x; 
  int by= blockIdx.y;
  int tx= threadIdx.x; 
  int ty= threadIdx.y;
  int Row= by * BLOCK_SIZE + ty; 
  int Col= bx * BLOCK_SIZE + tx; 
  float Pvalue= 0;
  int n_tiles = nx/BLOCK_SIZE + (nx%BLOCK_SIZE == 0 ? 0:1);
  for (int m= 0; m < n_tiles; ++m) 
  {
    if (Row < ny && m*BLOCK_SIZE+tx < nx)
      ds_M[ty][tx] = A[Row*nx + m*BLOCK_SIZE+tx];
    else
      ds_M[ty][tx] = 0;
 
    if (Col < ny && m*BLOCK_SIZE+ty < nx)
      ds_N[ty][tx] = A[Col*nx + m*BLOCK_SIZE+ty];
    else
      ds_N[ty][tx] = 0;
 
    __syncthreads();
	if(Col<Row)
      continue;
    for (int k = 0; k < BLOCK_SIZE; ++k)
      Pvalue += ds_M[ty][k] * ds_N[k][tx];
    __syncthreads();
  }
  if (Row < ny && Col < ny)
    C[Row*ny+Col] = Pvalue;
}
 
void correlate(int ny, int nx, const float* data, float* result) 
{
	float* inter = new float[ny*nx];
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for (int y = 0; y < ny; ++y) 
	{
		float *a_d_input;
		float *a_h;
		float *a_d_output;	
		float mean = 0.0;
		float sd = 0.0;
		size_t size = nx * sizeof(float);
		//Finding the mean
		for (int x = 0; x < nx; ++x) 
		{
			inter[x + y*nx] = data[x + y*nx];
			mean += inter[x + y*nx];
			//cout<<data[x+y*nx];
		}
		mean = mean/nx;
		//cout<<"Mean=="<<" "<<mean<<endl;
		//Finding the Standard Deviation
		a_h =  new float[nx];
		CHECK_CUDA_ERROR(cudaMalloc((void **) &a_d_input, size)); 
		CHECK_CUDA_ERROR(cudaMalloc((void **) &a_d_output, size)); 
		CHECK_CUDA_ERROR(cudaMemcpy(a_d_input, &inter[y*nx], size, cudaMemcpyHostToDevice));
		int block_size = 400;
		int n_blocks = nx/block_size + (nx%block_size == 0 ? 0:1);
		var<<< n_blocks, block_size >>> (a_d_input,a_d_output,nx,mean);
		CHECK_CUDA_ERROR(cudaGetLastError());
		CHECK_CUDA_ERROR(cudaMemcpy(a_h, a_d_output, size, cudaMemcpyDeviceToHost));
		//cout<<"Inter"<<endl;
		for (int x = 0; x < nx; ++x) 
		{
			sd += a_h[x];
		}
		//cout<<endl;
		sd= sqrt(sd);
		
		//cudaMalloc((void **) &a_d, size); 
		CHECK_CUDA_ERROR(cudaMemcpy(a_d_output, &inter[y*nx], size, cudaMemcpyHostToDevice));
		norm<<< n_blocks, block_size >>> (a_d_output,nx,mean,sd);
		CHECK_CUDA_ERROR(cudaGetLastError());
		CHECK_CUDA_ERROR(cudaMemcpy(&inter[y*nx], a_d_output, size, cudaMemcpyDeviceToHost));
		//cout<<"SD=="<<" "<<sd<<endl;
		//Finding zero mean and unit variance
		
   	}
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
	//cout << duration<<endl;
//cout<<"Done"<<endl;
	float *d_A;
	float *d_C;
	size_t size = nx*ny*sizeof(float);
	//float* inter_t = new float[ny*nx];
	//cout<<"Inter"<<endl;
	//for(int i =0;i< nx*ny;i++)
		//cout<< inter[i]<<" ";
	//cout<<endl;
	CHECK_CUDA_ERROR(cudaMalloc((void**) &d_A, size));
//	CHECK_CUDA_ERROR(cudaMalloc((void**) &d_B, size));
	CHECK_CUDA_ERROR(cudaMalloc((void**) &d_C, ny*ny*sizeof(float)));
	CHECK_CUDA_ERROR(cudaMemcpy(d_A, inter, size,cudaMemcpyHostToDevice));
	dim3 threads(BLOCK_SIZE,BLOCK_SIZE);
	int n_blocks = ny/BLOCK_SIZE + (ny%BLOCK_SIZE == 0 ? 0:1);
   	dim3 grid(n_blocks,n_blocks);
	high_resolution_clock::time_point t3 = high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::microseconds>( t3 - t2 ).count();
	
	matrixMultiply<<< grid, threads >>>(d_A, d_C,ny, nx);
	CHECK_CUDA_ERROR(cudaGetLastError());
	CHECK_CUDA_ERROR(cudaMemcpy(result, d_C, ny*ny*sizeof(float), cudaMemcpyDeviceToHost));
	high_resolution_clock::time_point t4 = high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::microseconds>( t4 - t3 ).count();
	cout << duration<<endl;
	//cout<<endl;	
	//for(int i =0;i< ny*ny;i++)
	//	cout<< result[i]<<" ";	
	//cout<<endl;
}
