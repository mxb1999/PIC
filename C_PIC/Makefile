COMP=g++
NVCC=/usr/local/cuda/bin/nvcc
CFLAGS=-g -Og  -lm -std=c++11
SRCDIR=src
OBJDIR=bin
CUODIR=bin/cuda
NAME=fdtd_c

INC=-Iinclude -I/usr/include/hdf5
LIB=-L/usr/lib/x86_64-linux-gnu/hdf5/serial/ -lhdf5 -lcuda -lcudart

_OBJ = time_loops.o setup_grid.o fdtd_main.o  initialize_fields.o hdf5_write.o rng.o push.o gather_current.o
OBJ= $(patsubst %, $(OBJDIR)/%, $(_OBJ))

_CUOBJ = push.o
CUOBJ= $(patsubst %, $(CUODIR)/%, $(_CUOBJ))

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(COMP) -c -MMD -MP $(CFLAGS) $(INC) $(LIB) -o $@ $^

$(CUODIR)/%.o: $(SRCDIR)/%.cu
	$(NVCC) -c -Iinclude $^ -o $@


$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(COMP) -c -MMD -MP $(CFLAGS) $(INC) $(LIB) -o $@ $^

$(NAME): $(OBJ) $(CUOBJ)
	$(NVCC) -dlink $(CUOBJ) -o $(OBJDIR)/all_cu.o
	h5c++ $(CFLAGS) $(INC) $(LIB) $^ $(OBJDIR)/all_cu.o -o $@

.phony: clean

clean:
	rm -f $(OBJDIR)/*.o
	rm -f $(NAME)

.phony: run

run:
	make
	./$(NAME)

.phony: reset

reset:
	make clean
	make