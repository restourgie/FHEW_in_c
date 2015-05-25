# the compiler to use.
CC=gcc
# options I'll pass to the compiler -fno-stack-protector
CFLAGS= -Wall -O3 -std=c11
LIBS= -lm


all: cmd

cmd: cmd/gen cmd/enc cmd/nand cmd/dec

FFT.o: FFT.h FFT.c params.h FHEW.h
	$(CC) $(CFLAGS) -c FFT.c $(LIBS)

LWE.o: LWE.h LWE.c FFT.h params.h distrib.h
	$(CC) $(CFLAGS) -c LWE.c $(LIBS)

FHEW.o: FHEW.h FHEW.c FFT.h LWE.h params.h
	$(CC) $(CFLAGS) -c FHEW.c $(LIBS)

common.o: cmd/common.c cmd/common.h 
	$(CC) $(CFLAGS) -c cmd/common.c $(LIBS)

distrib.o: distrib.c distrib.h 
	$(CC) $(CFLAGS) -c distrib.c $(LIBS)

params.o: params.c params.h
	$(CC) $(CFLAGS) -c params.c $(LIBS)

cmd/gen: cmd/gen.c common.o LWE.o FHEW.o distrib.o params.o FFT.o
	$(CC) $(CFLAGS) -o cmd/gen cmd/gen.c LWE.o FHEW.o common.o FFT.o distrib.o params.o $(LIBS)

cmd/enc: cmd/enc.c common.o FHEW.o FFT.o distrib.o params.o LWE.o
	$(CC) $(CFLAGS) -o cmd/enc cmd/enc.c common.o FHEW.o FFT.o distrib.o params.o LWE.o $(LIBS)

cmd/nand: cmd/nand.c common.o FHEW.o FFT.o params.o distrib.o LWE.o
	$(CC) $(CFLAGS) -o cmd/nand cmd/nand.c common.o FHEW.o FFT.o params.o distrib.o LWE.o $(LIBS)

cmd/dec: cmd/dec.c common.o LWE.o distrib.o params.o
	$(CC) $(CFLAGS) -o cmd/dec cmd/dec.c LWE.o common.o distrib.o params.o $(LIBS)

clean:
	rm *.o cmd/*.o cmd/gen cmd/enc cmd/dec cmd/nand

	
