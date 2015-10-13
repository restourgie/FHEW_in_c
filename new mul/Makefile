CC=/usr/bin/gcc
CFLAGS=-Wall -Wextra -O3 -std=c11 -mavx2
LIBS=-lm -lfftw3
DEPS = mul.h fft/fft_negacyc.h fft/split_radix_fft.h fft/sr_precomp.h fft/sr_vector.h fft/fftw.h fft/sr_vec_nonrec.h fft/twisted_fft.h fft/fft_negacyc_lut.h
CFILES = mul.c test/test_mul.c fft/fft_negacyc.c fft/split_radix_fft.c fft/sr_precomp.c fft/sr_vector.c fft/fftw.c fft/sr_vec_nonrec.c fft/twisted_fft.c fft/fft_negacyc_lut.c

test/test_mul: $(CFILES) $(DEPS)
	$(CC) $(CFLAGS) $(CFILES) $(LIBS) -o $@

.PHONY: clean

clean:
	-rm test/test_mul
