/*
 * DFT forward from real to complex  -- double precision
 */

#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>

#define N 1024 

int main(int argc, char *argv[])
{
    double in1[N];
    int i, ndata,ndata2,npow;

    fftw_complex  out1[N / 2 + 1];
    fftw_plan     p1;

	i=0;
	while ( fscanf(stdin,"%lg",&(in1[i]) )!= EOF ) 
	{
		i++;
		if (i > N) 
		{
		fprintf(stderr,"%s @%u in %s() line #=%d exceeds the limit %d.\n",
		__FILE__,__LINE__,__func__,i,N);
		exit(-1);
		}
	}
	ndata=i;
	if (ndata >0){
		npow=(int)ceil(log2((double)(ndata)));
		ndata2=(int)pow((double)2,(double)npow);
	}else
		ndata2=0;

	printf("Total number of data : %6d\n",ndata);
	printf("Total number of zero padded data : %6d\n",ndata2);
	printf("power of 2 : %d ( 2^%d =%6d )\n",npow, npow,ndata2);

// zero pad

	for(i=ndata;i<ndata2;i++) in1[i]=0.0;

// prepare for plans

    p1 = fftw_plan_dft_r2c_1d(ndata2, in1, out1, FFTW_ESTIMATE);

// real to complex DFT
    fftw_execute(p1);

#ifdef DBG
// compare original and forward+backward DFT 
	printf("##    Input            ##\n");
    for (int i = 0; i < ndata2; i++) {
          printf("%2d %15.10f\n", i, in1[i]);
    }
#endif

	printf("##    Re              Im               Amp           Phs ##\n");
    for (int i = 0; i < ndata2/2 +1; i++) {
          printf("%2d %15.10f %15.10f %15.10f %15.10f\n", i, creal(out1[i]), cimag(out1[i]),cabs(out1[i]),carg(out1[i]));
    }
    for (int i = ndata2/2+1; i < ndata2; i++) {
          printf("%2d %15.10f %15.10f %15.10f %15.10f\n", i, creal(out1[ndata2-i]), -cimag(out1[ndata2-i]),cabs(out1[ndata2-i]),-carg(out1[ndata2-i]));
	}


    fftw_destroy_plan(p1);

    return 0;
}
