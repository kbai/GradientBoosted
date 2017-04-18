CXX = g++ -std=c++11
executable:	latentfactor10.bin get_nmovienuser.bin dividedata.bin
%.o:	%.cpp
	$(CXX) -c  $*.cpp
latentfactor10.bin:  latentfactor10.o feature.o utils.o bkmodel.o bkmodel1.o	
	$(CXX) -o latentfactor10.bin latentfactor10.o utils.o bkmodel.o bkmodel1.o feature.o
get_nmovienuser.bin:    get_nmovienuser.o
	$(CXX) -o  get_nmovienuser.bin get_nmovienuser.o 
dividedata.bin:	divide_data_into_5_parts.o
	$(CXX) -o dividedata.bin divide_data_into_5_parts.o
clean: 
	rm *.o *.bin
