#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
using namespace std;
#include "commonheader.h"
#include "feature.hh"
#include "bkmodel.hh"
#include "bkmodel_gpu.hh"
#define NLAT 10
template<typename T>
void copying_to_gpu1d(void** a, vector<T>b)
{
	cudaMalloc(a, sizeof(T)*b.size());
	cudaMemcpy(*a,&(b[0]), sizeof(T)*b.size(),cudaMemcpyHostToDevice);
}

template<typename T>
void copying_to_gpu2d(void** a, vector<vector<T>> b)
{
	cudaMalloc(a, sizeof(T)*b[0].size()*b.size());
	auto it = b.begin();
	void* p = *a;
	while(it != b.end())
	{
		T* src = &((*it)[0]);
		size_t sz = it->size();
		cudaMemcpy(p, src, sizeof(T)*sz, cudaMemcpyHostToDevice);
		p+= sz;
		it++;
	}
}

bkmodel_gpu::bkmodel_gpu():bkmodel()
{
	copying_to_gpu1d((void**)&d_bu,bu);
	copying_to_gpu1d((void**)&d_btu,btu);
	copying_to_gpu2d((void**)&d_bm,bm);
	copying_to_gpu2d((void**)&d_bs,bs);
	copying_to_gpu2d((void**)&d_bt,bt);
		
}

void bkmodel_gpu::test()
{
	testgpu<<<1,128,128>>>(d_bm,d_bt);
}
