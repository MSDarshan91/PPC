#include "cp.h"
#include <math.h>
#include <cuda_runtime.h>
#include<iostream>
using namespace std;
#define BLOCK_SIZE 8
#define CHECK_CUDA_ERROR(call) do { \
        cudaError_t result_ = (call); \
        if (result_ != cudaSuccess) { \
            fprintf(stderr, #call " failed: %s\n", \
                    cudaGetErrorString(result_)); \
            exit(1); \
        } \
    } while(0)
__global__ void normalize(float *input,float *output, int ny, int nx)
	{
		float mean = 0.0;
		float sd = 0.0;
		int tx = threadIdx.x + (blockDim.x * blockIdx.x);
		int ty = threadIdx.y + (blockDim.y * blockIdx.y);
		int row = ty*nx;
		if(tx>= nx || ty>=ny)
			return;
		output[ty * nx + tx] = input[ty * nx + tx];
		for(int i=0;i<nx;i++)
			mean+=input[row+i];
		mean = mean/nx;
		float temp=0.0;
		for(int i=0;i<nx;i++)
			temp+=((input[row+i]-mean)*(input[row+i]-mean));
		sd = sqrt(temp);
		output[ty * nx + tx] = (output[ty * nx + tx] - mean)/sd;
    }

__global__ void matrixMul( float* C, float* A, int ny,int nx)
{
   int tx = threadIdx.x + (blockDim.x * blockIdx.x);
   int ty = threadIdx.y + (blockDim.y * blockIdx.y);
  
   if(tx>= ny || ty>=ny)
    return;
   float value = 0;
   for (int i = 0; i < nx; ++i)
   {
      float elementA = A[ty * nx + i];
      float elementB = A[tx * nx + i];
      value += elementA * elementB;
   }
    C[ty * ny + tx] = value;
}
 
void correlate(int ny, int nx, const float* data, float* result) 
{
	float* inter = new float[ny*nx];
	float *a_d_input;
	float *a_d_output;
	CHECK_CUDA_ERROR(cudaMalloc((void **) &a_d_input, nx*ny*sizeof(float))); 
	CHECK_CUDA_ERROR(cudaMalloc((void **) &a_d_output, nx*ny*sizeof(float))); 
	CHECK_CUDA_ERROR(cudaMemcpy(a_d_input, data, nx*ny*sizeof(float), cudaMemcpyHostToDevice));
	dim3 threads(BLOCK_SIZE,BLOCK_SIZE);
	int nx_blocks = nx/BLOCK_SIZE + (nx%BLOCK_SIZE == 0 ? 0:1);
	int ny_blocks = ny/BLOCK_SIZE + (ny%BLOCK_SIZE == 0 ? 0:1);
   	dim3 grid(nx_blocks,ny_blocks);
	normalize<<< grid, threads >>>(a_d_input, a_d_output, ny,nx);
	CHECK_CUDA_ERROR(cudaGetLastError());
	CHECK_CUDA_ERROR(cudaMemcpy(inter, a_d_output, nx*ny*sizeof(float), cudaMemcpyDeviceToHost));
	float *d_A;
	float *d_C;
	size_t size = nx*ny*sizeof(float);
//	float* inter_t = new float[ny*nx];
	//cout<<"Inter"<<endl;
	//for(int i =0;i< nx*ny;i++)
		//cout<< inter[i]<<" ";
	//cout<<endl;
	CHECK_CUDA_ERROR(cudaMalloc((void**) &d_A, size));
//	CHECK_CUDA_ERROR(cudaMalloc((void**) &d_B, size));
	CHECK_CUDA_ERROR(cudaMalloc((void**) &d_C, ny*ny*sizeof(float)));
	CHECK_CUDA_ERROR(cudaMemcpy(d_A, inter, size,cudaMemcpyHostToDevice));
	//cout<<n_blocks<<" "<< nx<<endl;
   	dim3 grid_1(ny_blocks,ny_blocks);
 	matrixMul<<< grid_1, threads >>>(d_C, d_A, ny,nx);
	CHECK_CUDA_ERROR(cudaGetLastError());
	CHECK_CUDA_ERROR(cudaMemcpy(result, d_C, ny*ny*sizeof(float), cudaMemcpyDeviceToHost));
	//cout<<endl;	
	//for(int i =0;i< ny*ny;i++)
	//	cout<< result[i]<<" ";	
	//cout<<endl;
}
