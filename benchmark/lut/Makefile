CC=/usr/bin/gcc
CFLAGS=-Wall -Wextra -O3 -std=c11 -mavx2
LIBS=-lm -lfftw3
DEPS = mul.h fft/sr_precomp.h fft/fft_negacyc_lut.h fft/tangent_fft.h
CFILES = mul.c test/test_mul.c fft/sr_precomp.c fft/fft_negacyc_lut.c fft/tangent_fft.c

test/test_mul: $(CFILES) $(DEPS)
	$(CC) $(CFLAGS) $(CFILES) $(LIBS) -o $@

.PHONY: clean

clean:
	-rm test/test_mul
