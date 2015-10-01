CC=/usr/bin/gcc
CFLAGS=-Wall -Wextra -g -std=c11 -mavx
LIBS=-lm
DEPS = mul.h fft/fft_negacyc.h fft/split_radix_fft.h fft/sr_precomp.h fft/sr_vector.h
CFILES = mul.c test/test_mul.c fft/fft_negacyc.c fft/split_radix_fft.c fft/sr_precomp.c fft/sr_vector.c

test/test_mul: $(CFILES) $(DEPS)
	$(CC) $(CFLAGS) $(CFILES) $(LIBS) -o $@

.PHONY: clean

clean:
	-rm test/test_mul
