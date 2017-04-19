#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
using namespace std;
#include "commonheader.h"
#include "feature.hh"
#include "bkmodel.hh"
#include "bkmodel_gpu.hh"
#define NLAT 10

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
	   if (code != cudaSuccess) 
		      {
				        fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
						      if (abort) exit(code);
							     }
}




__global__ void testgpu(float * a, float* b)
{
	int it = threadIdx.x + 16 * threadIdx.y;
}




template<typename T>
void copying_to_gpu1d(T** a, vector<T>b)
{
	gpuErrchk(cudaMalloc((void**)a, sizeof(T)*b.size()));
	gpuErrchk(cudaMemcpy(*a,&(b[0]), sizeof(T)*b.size(),cudaMemcpyHostToDevice));
}

template<typename T>
void copying_to_gpu2d(T** a, vector<vector<T>> b)
{
	gpuErrchk(cudaMalloc((void**)a, sizeof(T)*b[0].size()*b.size()));
	auto it = b.begin();
	T* p = *a;
	while(it != b.end())
	{
		T* src = &((*it)[0]);
		size_t sz = it->size();
		gpuErrchk(cudaMemcpy(p, src, sizeof(T)*sz, cudaMemcpyHostToDevice));
		p+= sz;
		it++;
	}
}

bkmodel_gpu::bkmodel_gpu():bkmodel()
{
	copying_to_gpu1d(&d_bu,bu);
	copying_to_gpu1d(&d_bt,bt);
	copying_to_gpu1d(&d_bm,bm);
	copying_to_gpu2d(&d_btm,btm);
	copying_to_gpu1d(&d_btu,btu);
		
}

void bkmodel_gpu::test()
{
	float a[16];
/*	testgpu<<<16,16>>>(d_bm,d_bt);
	cudaMemcpy(a,d_bm, 16*sizeof(float),cudaMemcpyDeviceToHost);
	cout << a[10] <<"\t"<< a[11] <<"\t"<< a[12] <<endl;
	cout << bm[1][0] << "\t"<< bm[1][1] << "\t" << bm[1][2] << endl;
	gpuErrchk(cudaPeekAtLastError());
*/}
