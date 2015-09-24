CC=/usr/bin/gcc
CFLAGS=-Wall -Wextra -g -std=c11 -O0 -v -da -Q
LIBS=-lm
DEPS = mul.h fft/fft_negacyc.h fft/split_radix_fft.h fft/split_radix_non_rec.h
CFILES = mul.c test/test_mul.c fft/fft_negacyc.c fft/split_radix_fft.c fft/split_radix_non_rec.c

test/test_mul: $(CFILES) $(DEPS)
	$(CC) $(CFLAGS) $(CFILES) $(LIBS) -o $@

.PHONY: clean

clean:
	-rm test/test_mul
