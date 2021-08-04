COMP=gcc
CFLAGS= -g -Og -Wall -Werror -std=c99 -MMD -MP -lm -Wno-unused-variable
INC=-Iinclude
ODIR=bin
SDIR=src
IDIR=include



_COBJ = pic.o
COBJ = $(patsubst %, $(ODIR)/%,$(_COBJ))

$(ODIR)/%.o: $(SDIR)/%.c
	$(COMP) -c $(CFLAGS) $(INC) $^ -o $@

pic_exec: $(COBJ)
	$(COMP) $(CFLAGS) $(INC) $^ -o $@


clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~
	rm -f $(ODIR)/*.d *~ core $(INCDIR)/*~

.phony: reset

reset:
	make clean
	make



.phony:run

run:
	make
	./pic_exec
