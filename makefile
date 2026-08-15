# Makefile of FastqCheck (2022/07/29)

CC = gcc
CFLAGS = -std=c99 -D_GNU_SOURCE
LIBS = -lz -lbz2 -lm -lpthread
INCLUDE =

DEBUG = 0

ifeq ($(DEBUG), 1)
	CFLAGS += -g -O0 # enable debugging
else
	CFLAGS += -O3
endif


OBJECT = bloom_filter.o \
		 file_read.o \
		 file_type.o \
		 kqueue.o \
		 params_parse.o \
		 threadpool.o\
		 fastq_check.o

PROG = fastq_check

all: $(PROG)

$(PROG): $(OBJECT) $(HTSLIB)
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ $^ $(LIBS)

# generate object file (*.o) for each source file (*.c)
%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ -c $<


.PHONY : clean
clean:
	rm -rf $(OBJECT) $(PROG) 
