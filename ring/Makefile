CC=/usr/bin/gcc
CFLAGS=-Wall -Wextra -g

test/test_mul: ring.c ring.h test/test_mul.c
	$(CC) $(CFLAGS) ring.c test/test_mul.c -o $@

.PHONY: clean

clean:
	-rm test/test_mul
