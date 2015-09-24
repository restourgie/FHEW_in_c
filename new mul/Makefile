CC=/usr/bin/gcc
CFLAGS=-Wall -Wextra -g -std=c11
LIBS=-lm
DEPS = mul.h support.h fft/fft_negacyc.h fft/split_radix_fft.h fft/split_radix_fast.h fft/split_radix_non_rec.h
CFILES = mul.c test/test_mul.c fft/fft_negacyc.c fft/split_radix_fft.c fft/split_radix_non_rec.c fft/split_radix_fast.c support.c

test/test_mul: $(CFILES) $(DEPS)
	$(CC) $(CFLAGS) $(CFILES) $(LIBS) -o $@

.PHONY: clean

clean:
	-rm test/test_mul
