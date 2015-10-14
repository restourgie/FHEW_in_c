CC=/usr/bin/gcc
CFLAGS=-Wall -Wextra -O3 -std=c11
LIBS=-lm
DEPS = mul.h fft/normal_fft.h fft/split_radix_fft.h fft/twisted_fft.h fft/tangent_fft.h
CFILES = mul.c test/test_mul.c fft/normal_fft.c fft/split_radix_fft.c fft/twisted_fft.c fft/tangent_fft.c

test/test_mul: $(CFILES) $(DEPS)
	$(CC) $(CFLAGS) $(CFILES) $(LIBS) -o $@

.PHONY: clean

clean:
	-rm test/test_mul
