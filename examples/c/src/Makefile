

IDIR=../include
ODIR=obj
LDIR=../lib

LIBS=-lm

CC=gcc
CFLAGS=-I$(IDIR)

_HEADERS = hello.h
HEADERS = $(patsubst %,$(IDIR)/%,$(_HEADERS))

_OBJS = hello.o hellofunc.o 
OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))

# build the executable from the object files
hello: $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

# compile a single .c file to an .o file
$(ODIR)/%.o: %.c $(HEADERS)
	$(CC) -c -o $@ $< $(CFLAGS)


# clean up temporary files
.PHONY: clean
clean:
	rm -f $(ODIR)/*.o *~ core $(IDIR)/*~ 
