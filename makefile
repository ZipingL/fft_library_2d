CXX = gcc
CXXFLAGS = -lm
INC= -I ./
INC+= -I ../alg/
INC+= -I ../../c7x_common/
OBJECTS = main.o FFTLIB_fft1dBatched_i16sc_c16sc_o16sc_cn.o FFT2d_cn.o 

all: $(OBJECTS)
	$(CXX) -o $@ $^ $(CXXFLAGS)

main.o: main.c FFT2d_cn.h
	$(CXX)  $(INC) -c main.c $(CXXFLAGS)
	
	
FFTLIB_fft1dBatched_i16sc_c16sc_o16sc_cn.o: ../alg/FFTLIB_fft1dBatched_i16sc_c16sc_o16sc_cn.c ../../c7x_common/VXLIB_types.h
	$(CXX)  $(INC) -c ../alg/FFTLIB_fft1dBatched_i16sc_c16sc_o16sc_cn.c $(CXXFLAGS)


FFT2d_cn.o: FFT2d_cn.c FFT2d_cn.h
	$(CXX)  $(INC) -c FFT2d_cn.c $(CXXFLAGS)
	
clean:
	rm all *.o