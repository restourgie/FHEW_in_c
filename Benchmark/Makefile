CC=gcc
#-fno-stack-protector
CFLAGS= -Wall -O3 -g
LDFLAGS= -lfftw3 -lm

all: benchmark

benchmark:  benchmark.c
	$(CC) $(CFLAGS) -o benchmark benchmark.c $(LDFLAGS)

clean: 
	rm benchmark