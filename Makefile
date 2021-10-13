COMP=g++
NVCOMP=nvcc
SRCDIR=src
ODIR=bin
LIBS=-L/usr/local/cuda/lib64 -L/usr/lib/x86_64-linux-gnu/hdf5/serial
INC=-Iinclude -I/usr/lib/x86_64-linux-gnu/openmpi/include -I/usr/local/cuda/include
CFLAGS=$(INC) $(LIBS) -Wall -g -std=c++11 -fopenmp -lm -lhdf5 
NVFLAGS=$(INC) $(LIBS) -g -G -arch=sm_70 -lm -lhdf5

_OBJ = main.o push.o grid.o
OBJ = $(patsubst %, $(ODIR)/%, $(_OBJ))

$(ODIR)/%.o: $(SRCDIR)/%.cpp
	$(COMP) -lcuda -lcudart $(CFLAGS) -c -o $@ $^

$(ODIR)/%.o: $(SRCDIR)/%.cu
	$(NVCOMP) -Iinclude -arch=sm_70 -c -o $@ $^ 

cupush:
	$(NVCOMP) $(NVFLAGS) src/main.cpp src/grid.cpp src/push.cu -o $@

.phony: clean

clean:
	rm cupush

.phony: reset

reset:
	make clean
	make
