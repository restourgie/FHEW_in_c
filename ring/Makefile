CC=/usr/bin/gcc
CFLAGS=-Wall -Wextra -g
LIBS=-lm

test/test_mul: ring.c ring.h fft.c fft.h test/test_mul.c
	$(CC) $(CFLAGS) ring.c fft.c test/test_mul.c $(LIBS) -o $@

.PHONY: clean

clean:
	-rm test/test_mul
