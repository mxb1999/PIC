COMP=gcc
CFLAGS=-g -Og  -lm -std=c99
SRCDIR=src
OBJDIR=bin
NAME=fdtd_c

INC=-Iinclude -I/usr/include/hdf5
LIB=-L/usr/lib/hdf5 -lhdf5

_OBJ = time_loops.o setup_grid.o fdtd_main.o  initialize_fields.o hdf5_write.o
OBJ= $(patsubst %, $(OBJDIR)/%, $(_OBJ))


$(OBJDIR)/%.o: $(SRCDIR)/%.c 
	$(COMP) -c -MMD -MP $(CFLAGS) $(INC) $(LIB) -o $@ $^  

$(NAME): $(OBJ)
	$(COMP) $(CFLAGS) $(INC) $(LIB) $^ -o $@

.phony: clean

clean:
	rm $(OBJDIR)/*.o
	rm $(NAME)

.phony: run

run:
	make
	./$(NAME)

.phony: reset

reset:
	make clean
	make