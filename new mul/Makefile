CC=/usr/bin/gcc
CFLAGS=-Wall -Wextra -g -std=c11
LIBS=-lm

test/test_mul: mul.c mul.h test/test_mul.c
	$(CC) $(CFLAGS) mul.c test/test_mul.c $(LIBS) -o $@

.PHONY: clean

clean:
	-rm test/test_mul
