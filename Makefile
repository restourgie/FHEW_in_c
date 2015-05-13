# the compiler to use.
CC=gcc
# options I'll pass to the compiler
CFLAGS= -Wall -O3 -std=c11

all: cmd

cmd: cmd/gen cmd/enc cmd/nand cmd/dec

FFT.o: FFT.h FFT.c params.h FHEW.h
	$(CC) $(CFLAGS) -c FFT.c

LWE.o: LWE.h LWE.c FFT.h params.h distrib.h
	$(CC) $(CFLAGS) -c LWE.c

FHEW.o: FHEW.h FHEW.c FFT.h LWE.h params.h
	$(CC) $(CFLAGS) -c FHEW.c

common.o: cmd/common.c cmd/common.h 
	$(CC) $(CFLAGS) -c cmd/common.c 

cmd/gen: cmd/gen.c common.o 
	$(CC) $(CFLAGS) -o cmd/gen cmd/gen.c common.o

cmd/enc: cmd/enc.c common.o 
	$(CC) $(CFLAGS) -o cmd/enc cmd/enc.c common.o

cmd/nand: cmd/nand.c common.o 
	$(CC) $(CFLAGS) -o cmd/nand cmd/nand.c common.o

cmd/dec: cmd/dec.c common.o 
	$(CC) $(CFLAGS) -o cmd/dec cmd/dec.c common.o

clean:
	rm *.o cmd/*.o fhewTest cmd/gen cmd/enc cmd/dec cmd/nand