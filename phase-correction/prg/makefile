# Intel Compiler + Intel MKL
CC = icc
CFLAGS = -O3 -std=c99 -xHost
CPPFLAGS = -I${MKLROOT}/include/fftw
# END of Intel Compiler

fftw3-test: fftw3-test.c
	$(LINK.c) fftw3-test.c -mkl -o fftw3-test
fftw3: fftw3.c
	$(LINK.c) fftw3.c -mkl -o fftw3
dft: dft.c
	$(LINK.c) dft.c -mkl -o dft
dft_c2r: dft_c2r.c
	$(LINK.c) dft_c2r.c -mkl -o dft_c2r
dft_r2c: dft_r2c.c
	$(LINK.c) dft_r2c.c -mkl -o dft_r2c
clean:
	rm -rf *.o a.out fftw3-test fftw3 dft dft_c2r dft_r2c
