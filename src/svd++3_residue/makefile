CXX = nvcc -gencode arch=compute_20,code=sm_20  -std=c++11 -O3 
executable:	latentfactor10.bin lat_gpu.bin  get_nmovienuser.bin dividedata.bin
bkmodel_gpu.o: bkmodel_gpu_old.cu
	$(CXX) -c  bkmodel_gpu_old.cu
%.o:	%.cpp
	$(CXX) -c  $*.cpp
latentfactor10.bin:  latentfactor10.o feature.o utils.o bkmodel.o bkmodel1.o 	
	$(CXX) -o latentfactor10.bin latentfactor10.o utils.o bkmodel.o bkmodel1.o feature.o
lat_gpu.bin:  latentfactor10_gpu.o feature.o utils.o bkmodel.o bkmodel1.o bkmodel_gpu.o	
	$(CXX) -o lat_gpu.bin latentfactor10_gpu.o utils.o bkmodel.o bkmodel1.o feature.o bkmodel_gpu_old.o

get_nmovienuser.bin:    get_nmovienuser.o
	$(CXX) -o  get_nmovienuser.bin get_nmovienuser.o 
dividedata.bin:	divide_data_into_5_parts.o
	$(CXX) -o dividedata.bin divide_data_into_5_parts.o
rbm:	
	mpicxx -std=c++11 -o rbm  ./RBM.cpp
clean: 
	rm *.o *.bin
