/*
 * DFT backward from complex to real -- double precision
 */

#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>

#define N 1024 

int main(int argc, char *argv[])
{
    double in1[N], amp, phs;
    int i, ndata,ndata2,npow;

    fftw_complex  out1[N / 2 + 1], out2[N/2+1];
    fftw_plan     p1;

	i=0;
	while ( fscanf(stdin,"%lg %lg",&amp, &phs )!= EOF ) 
	{
		out1[i]=amp*cexp(I*phs);
#ifdef DBG
	printf("%15.10f %15.10f %15.10f %15.10f\n",amp,phs,creal(out1[i]),cimag(out1[i]));
#endif
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
	for(i=ndata;i<ndata2/2+1;i++) out1[i]=0.0+0.0*I;

	for(i=0;i<ndata2/2+1;i++) out2[i]=out1[i];
// prepare for plans

    p1 = fftw_plan_dft_c2r_1d(ndata2, out1, in1, FFTW_ESTIMATE);

//  complex to real DFT
    fftw_execute(p1);

#ifdef DBG
// compare original and forward+backward DFT 
	printf("## Input Re           Im               Amp           Phs ##\n");
    for (int i = 0; i < ndata2/2 +1; i++) {
          printf("%2d %15.10f %15.10f %15.10f %15.10f\n", i, creal(out2[i]), cimag(out2[i]),cabs(out2[i]),carg(out2[i]));
    }
#endif

	printf("##    Output         ##\n");
    for (int i = 0; i < ndata2; i++) {
          printf("%2d %15.10f\n", i, in1[i]/ndata2);
    }


    fftw_destroy_plan(p1);

    return 0;
}
