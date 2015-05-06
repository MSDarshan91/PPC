#include "cp.h"
#include <math.h>
#include <cuda_runtime.h>
#include<iostream>
using namespace std;
#define BLOCK_SIZE 10
#define TILE_SIZE  2
#define CHECK_CUDA_ERROR(call) do { \
        cudaError_t result_ = (call); \
        if (result_ != cudaSuccess) { \
            fprintf(stderr, #call " failed: %s\n", \
                    cudaGetErrorString(result_)); \
            exit(1); \
        } \
    } while(0)
__global__ void var(double *input,double *output, int N, double mean)
	{
 
      int idx=threadIdx.x+(blockDim.x*blockIdx.x);
      if (idx < N) output[idx] = (input[idx]-mean)*(input[idx]-mean);
    }
__global__ void norm(double *input, int N,double mean,double sd)
	{
 
      int idx=threadIdx.x+(blockDim.x*blockIdx.x);
      if (idx < N) input[idx] =  (input[idx]-mean)/sd;
    }

__global__ void matrixMul( float* C, double* A, int ny,int nx)
{
   int tx = threadIdx.x + (blockDim.x * blockIdx.x);
   int ty = threadIdx.y + (blockDim.y * blockIdx.y);
  
   if(tx>= ny || ty>=ny)
    return;
   double value = 0;
   for (int i = 0; i < nx; ++i)
   {
      double elementA = A[ty * nx + i];
      double elementB = A[tx * nx + i];
      value += elementA * elementB;
   }
    C[ty * ny + tx] = value;
}
 
void correlate(int ny, int nx, const float* data, float* result) 
{
	double* inter = new double[ny*nx];
	for (int y = 0; y < ny; ++y) 
	{
		double *a_d_input;
		double *a_h;
		double *a_d_output;	
		double mean = 0.0;
		double sd = 0.0;
		size_t size = nx * sizeof(double);
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
		a_h =  new double[nx];
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
	//cout<<"Done"<<endl;
	double *d_A;
	float *d_C;
	size_t size = nx*ny*sizeof(double);
//	double* inter_t = new double[ny*nx];
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
	//cout<<n_blocks<<" "<< nx<<endl;
   	dim3 grid(n_blocks,n_blocks);
 	matrixMul<<< grid, threads >>>(d_C, d_A, ny,nx);
	CHECK_CUDA_ERROR(cudaGetLastError());
	CHECK_CUDA_ERROR(cudaMemcpy(result, d_C, ny*ny*sizeof(float), cudaMemcpyDeviceToHost));
	//cout<<endl;	
	//for(int i =0;i< ny*ny;i++)
	//	cout<< result[i]<<" ";	
	//cout<<endl;
}
