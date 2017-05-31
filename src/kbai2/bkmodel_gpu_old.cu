#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
using namespace std;
#include "commonheader.h"
#include "feature.hh"
#include "bkmodel.hh"
#include <stdio.h>
#include "bkmodel_gpu_old.hh"

#define d_pm(i,j) d_pm[NLAT*i+j]
#define d_pu(i,j) d_pu[NLAT*i +j]
#define d_pu1(i,j) d_pu1[NLAT*i+j]
#define d_btm(i,j) d_btm[30*i+j]
#define d_bfm(i,j) d_bfm[8*i+j]
#define PUT_FUNCTION_G \
	for(int i = 0 ; i < NLAT;i++)	\
	{\
		error -= d_pm(im,i) * (d_pu(iu,i) + d_pu1(iu,i)* tt);\
	};\
	error -= MEAN + d_bu[iu] + d_bm[im]  + d_btm(im,it)*d_btu[iu] + d_bt[it] + d_bta[ita] + d_bf[ife]+d_bfm(im,ife);


#define gpuErrchk(ans,i) { gpuAssert((ans), i, __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, int i, const char *file, int line, bool abort=true)
{
	   if (code != cudaSuccess) 
		      {
				        fprintf(stderr,"GPUassert: %s %s %d %d\n", cudaGetErrorString(code), file, line,i);
						      if (abort) exit(code);
							     }
}

__device__ double atomicAdd(double* address, double val)
{
	    unsigned long long int* address_as_ull = (unsigned long long int*)address;

		    unsigned long long int old = *address_as_ull, assumed;

			    do{ assumed = old;
							old = atomicCAS(address_as_ull, assumed,__double_as_longlong(val +__longlong_as_double(assumed)));
							    } while (assumed != old);

				    return __longlong_as_double(old);
}


__device__ float ug(float a)
{
	if(a>0) return 0.0363636*pow((double)a,0.4);
	else return -0.0363636*pow((double)abs(a),0.4);

}
__global__ void kernel_sgd(
	//model
 	float*d_bu,
 	float*d_bm,
 	float*d_pm,
 	float*d_pu,
 	float*d_pu1,
 	float*d_btu,
 	float*d_bf,
 	float*d_bt,
	float*d_bta,
 	float*d_bfm,
	float*d_btm,
	//data
	int*du,
	int*dm,
	int*dt,
	int*df,
	float*dtt,
	int*rate,
	float _lr,
	int shift,
	int sz
	)
{
	int ith =( shift + threadIdx.x + blockIdx.x * blockDim.x );
	while(ith < sz)
	{
	
	int iu = du[ith];
	int im = dm[ith];
	int it = (int)(dt[ith]/75);
	int ita = dt[ith];
	int ife = df[ith];
//	printf("gpu:%d\t%d\t%d\t%d\t%d\n",iu,im,it,ita,ife);
	float tt =ug(dtt[ith]);
	int ir = rate[ith];
	float tmp;
	float error = (float)ir;
	PUT_FUNCTION_G

	tmp =_lr*(error - 0.008*d_bu[iu]); //no need to regulating thses terms
	atomicAdd(&d_bu[iu],tmp);

	tmp =_lr*(error - 0.008*d_bm[im]);
	atomicAdd(&d_bm[im],tmp);


	for(int i = 0; i < NLAT; i++)
	{
		tmp = _lr*(error * (d_pu(iu,i) + d_pu1(iu,i)*tt) - 0.015 * d_pm(im,i)) ;
		atomicAdd(&(d_pm(im,i)),tmp);

		tmp = _lr*(error * d_pm(im,i) -0.015 * d_pu(iu,i));
		atomicAdd(&(d_pu(iu,i)),tmp);
	//	tmp = 0.5*_lr*(error*tt*d_pm(im,i)-0.015*d_pu1(iu,i));//		}
// step3 	
		tmp =1.0*_lr*(error*tt*d_pm(im,i) - 0.015*d_pu1(iu,i));//		}
		atomicAdd(&(d_pu1(iu,i)),tmp);

	}
	tmp = _lr*(error*d_btu[iu] - 0.015 * d_btm(im,it));
	atomicAdd(&(d_btm(im,it)), tmp);
	tmp = _lr*(error*d_btm(im,it) - 0.015* d_btu[iu]);
	atomicAdd(&(d_btu[iu]), tmp);
	tmp = 0.1*_lr*( error - 0.008 * d_bt[it]);
	atomicAdd(&(d_bt[it]), tmp);
	tmp = 0.5*_lr*(error - 0.008 * d_bf[ife]);
	atomicAdd(&(d_bf[ife]), tmp);

	tmp = _lr*(error - 0.008 * d_bfm(im,ife));
	atomicAdd(&(d_bfm(im,ife)),tmp);
	tmp = _lr*(error - 0.008 * d_bta[ita]);
	atomicAdd(&(d_bta[ita]),tmp);
	ith += blockDim.x * gridDim.x;
	__syncthreads();
}
	
//	atomicAdd(d_mean, 0.002*error);
//

	return;
}

__global__ void kernel_rmse(
	//model
 	float*d_bu,
 	float*d_bm,
 	float*d_pm,
 	float*d_pu,
 	float*d_pu1,
 	float*d_btu,
 	float*d_bf,
 	float*d_bt,
	float*d_bta,
 	float*d_bfm,
	float*d_btm,
	//data
	int*du,
	int*dm,
	int*dt,
	int*df,
	float*dtt,
	int*rate,
	int sz,
	double *sum
	)
{
	int ith =( blockIdx.x * blockDim.x + threadIdx.x );
	
	while(ith < sz)
	{

		int iu = du[ith];
		int im = dm[ith];
		int it = (int)(dt[ith]/75);
		int ita = dt[ith];
		int ife = df[ith];
//	printf("%d\n",ith);
		float tt =ug(dtt[ith]);
		int ir = rate[ith];
		double error = (double)ir;

		PUT_FUNCTION_G
//		printf("gpu:%d\t %f\n",iu, error);

		atomicAdd(sum, (1.0*error*error));
		ith += gridDim.x* blockDim.x;
//
	}

	return;
}



template<typename T>
void copying_to_gpu1d(T** a, vector<T>&b,int err = 0)
{
	gpuErrchk(cudaMalloc((void**)a, sizeof(T)*b.size()),err);
	gpuErrchk(cudaMemcpy(*a,&(b[0]), sizeof(T)*b.size(),cudaMemcpyHostToDevice),err);
}

template<typename T>
void copying_to_gpu2d(T** a, vector<vector<T>> &b, int err=0)
{
	gpuErrchk(cudaMalloc((void**)a, sizeof(T)*b[0].size()*b.size()), err);
	auto it = b.begin();
	T* p = *a;
	while(it != b.end())
	{
		T* src = &((*it)[0]);
		size_t sz = it->size();
		gpuErrchk(cudaMemcpy(p, src, sizeof(T)*sz, cudaMemcpyHostToDevice), err);
		p+= sz;
		it++;
	}
}

template<typename T>
void copying_to_cpu1d(T** a, vector<T>&b, int err=0)
{
	cout << b.size() << endl;
	gpuErrchk(cudaMemcpy(&(b[0]),*a, sizeof(T)*b.size(),cudaMemcpyDeviceToHost), err);
}

template<typename T>
void copying_to_cpu2d(T** a, vector<vector<T>>& b, int err=0)
{
	auto it = b.begin();
	T* p = *a;
	while(it != b.end())
	{
		T* src = &((*it)[0]);
		size_t sz = it->size();
		gpuErrchk(cudaMemcpy(src, p, sizeof(T)*sz, cudaMemcpyDeviceToHost),err);
		p+= sz;
		it++;
	}
}
bkmodel_gpu::bkmodel_gpu():bkmodel()
{
	copying_to_gpu1d(&d_bu,bu);
	copying_to_gpu1d(&d_btu,btu);
	copying_to_gpu1d(&d_bt,bt);
	copying_to_gpu1d(&d_bf,bf);
	copying_to_gpu1d(&d_bm,bm);
	copying_to_gpu1d(&d_bta,bta);
	copying_to_gpu2d(&d_pm,pm);
	copying_to_gpu2d(&d_pu,pu);
	copying_to_gpu2d(&d_pu1,pu1);
	copying_to_gpu2d(&d_bfm, bfm);
	copying_to_gpu2d(&d_btm, btm);

}
void bkmodel_gpu::retrieve_gpu()
{
	copying_to_cpu1d(&d_bu,bu,1001);
	copying_to_cpu1d(&d_btu,btu,1002);
	copying_to_cpu1d(&d_bt,bt,1003);
	copying_to_cpu1d(&d_bf,bf,1004);
	copying_to_cpu1d(&d_bm,bm,1005);
	copying_to_cpu1d(&d_bta,bta,1011);
	copying_to_cpu2d(&d_pm,pm,1006);
	copying_to_cpu2d(&d_pu,pu,1007);
	copying_to_cpu2d(&d_pu1,pu1,1008);
	copying_to_cpu2d(&d_bfm, bfm,1009);
	copying_to_cpu2d(&d_btm, btm,1010);
}

void bkmodel_gpu::loaddata(feature &a0, feature &a1, feature &a2)
{
	sz0 = a0.viu.size();
	copying_to_gpu1d(&d_iu,a0.viu);
	copying_to_gpu1d(&d_im,a0.vim);
	copying_to_gpu1d(&d_it,a0.vita); //ita is from 1-2243
	copying_to_gpu1d(&d_if,a0.vif);
	copying_to_gpu1d(&d_tb,a0.vtb);
	copying_to_gpu1d(&d_rate,a0.vrate);

	sz1 = a1.viu.size();
	copying_to_gpu1d(&d_iu1,a1.viu);
	copying_to_gpu1d(&d_im1,a1.vim);
	copying_to_gpu1d(&d_it1,a1.vita);
	copying_to_gpu1d(&d_if1,a1.vif);
	copying_to_gpu1d(&d_tb1,a1.vtb);
	copying_to_gpu1d(&d_rate1,a1.vrate);

	sz2 = a2.viu.size();
	copying_to_gpu1d(&d_iu2,a2.viu);
	copying_to_gpu1d(&d_im2,a2.vim);
	copying_to_gpu1d(&d_it2,a2.vita);
	copying_to_gpu1d(&d_if2,a2.vif);
	copying_to_gpu1d(&d_tb2,a2.vtb);
	copying_to_gpu1d(&d_rate2,a2.vrate);

}
double bkmodel_gpu::compute_error()
{
	double sum1 = 0, *dev_sum1;
	cudaMalloc((void**) &dev_sum1,sizeof(double));
	cudaMemcpy(dev_sum1, &sum1, sizeof(double), cudaMemcpyHostToDevice);

	double sum2 = 0, *dev_sum2;
	cudaMalloc((void**) &dev_sum2,sizeof(double));
	cudaMemcpy(dev_sum2, &sum2, sizeof(double), cudaMemcpyHostToDevice);

	kernel_rmse<<<32,32>>>(
	//model
	d_bu,
 	d_bm,
 	d_pm,
 	d_pu,
 	d_pu1,
 	d_btu,
 	d_bf,
 	d_bt,
	d_bta,
 	d_bfm,
 	d_btm,
	//data
	d_iu1,
	d_im1,
	d_it1,
	d_if1,
	d_tb1,
	d_rate1,
	sz1,
	dev_sum1);
	gpuErrchk(cudaPeekAtLastError(),400);

	kernel_rmse<<<32,32>>>(
	//model
	d_bu,
 	d_bm,
 	d_pm,
 	d_pu,
 	d_pu1,
 	d_btu,
 	d_bf,
 	d_bt,
	d_bta,
 	d_bfm,
 	d_btm,
	//data
	d_iu2,
	d_im2,
	d_it2,
	d_if2,
	d_tb2,
	d_rate2,
	sz2,
	dev_sum2);
	gpuErrchk(cudaPeekAtLastError(),400)
	cudaMemcpy( &sum1, dev_sum1,sizeof(double),  cudaMemcpyDeviceToHost);
	cudaMemcpy( &sum2, dev_sum2,sizeof(double),  cudaMemcpyDeviceToHost);
	cout << "in sample error: " << sqrt(sum1/sz1) << endl;
	cout << "out sample error: " <<sqrt(sum2/sz2) << endl;
	double rmse = sqrt(sum2/sz2);

	return rmse;


}

void bkmodel_gpu::test(float lr)
{
	kernel_sgd<<<32,32>>>(
	//model
	d_bu,
 	d_bm,
 	d_pm,
 	d_pu,
 	d_pu1,
 	d_btu,
 	d_bf,
 	d_bt,
	d_bta,
 	d_bfm,
 	d_btm,
	//data
	d_iu,
	d_im,
	d_it,
	d_if,
	d_tb,
	d_rate,
	lr,
	0,
	sz0);
	gpuErrchk(cudaPeekAtLastError(),300);



/*	testgpu<<<16,16>>>(d_bm,d_bt);
:	cudaMemcpy(a,d_bm, 16*sizeof(float),cudaMemcpyDeviceToHost);
	cout << a[10] <<"\t"<< a[11] <<"\t"<< a[12] <<endl;
	cout << bm[1][0] << "\t"<< bm[1][1] << "\t" << bm[1][2] << endl;
	gpuErrchk(cudaPeekAtLastError());
*/}
