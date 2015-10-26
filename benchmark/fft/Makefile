CC=/usr/bin/gcc
CFLAGS=-Wall -Wextra -O3 -std=c11 -mavx2 -march=native 
LIBS=-lm -lfftw3
DEPS = mul.h fftw.h negacyclic.h
CFILES = test_mul.c fftw.c negacyclic.c

test_mul: $(CFILES) $(DEPS)
	$(CC) $(CFLAGS) $(CFILES) $(LIBS) -o $@

.PHONY: clean

clean:
	-rm test_mul
