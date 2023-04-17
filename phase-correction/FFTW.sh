#!/bin/bash
#
gmt_version=0
if command -v gmt &> /dev/null; then
	version=`gmt --version`
	gmt_version=${version:0:1}
	if [ $gmt_version -eq 5 ]; then
		HEADER_OPTION="h"
		NAN_OPTION="-s"
	fi
elif command -v GMT &> /dev/null ; then
	gmt_version=4
	HEADER_OPTION="H"
	NAN_OPTION=""
fi
#
if [ $gmt_version -eq 0 ]; then
	echo "### ERROR : GMT is not found in the path ###"
	exit
elif [ $gmt_version ]; then
	echo "### GMT version $gmt_version is found ###"
elif [ $gmt_version -gt 5 ]; then
	echo " ### ERROR : GMT version 6 or later is not supported ###"
exit
fi
#
if ! command -v sac $> /dev/null; then
	echo "SAC is not found in the PATH."
	exit
fi
#
# prepare phase data
#
sac << END
cut 40000 60000
r tsunami_4km_constv_9000km.sac
taper
w alpha tsunami_4km_constv_9000km_t40000sec.cut.sac
quit
END
#
# extract one-colummns time-series file from  SAC ascii file
#
tail -134 tsunami_4km_constv_9000km_t40000sec.cut.sac | awk '{for(i=1;i<=NF;i++) print $i}' > j1
#
# discreate fourier transform (DFT) by fftw3 program and measure amp and phase
#
prg/dft_r2c < j1 | awk 'NR > 4 {print $1,$4,$5}' > j2
#
# prepare knot table for interpolation
# df(Hz)=1.0/(30(sec)*1024(points))
#
seq 0 1023 |  awk 'BEGIN{pi=4.0*atan2(1.0,1.0);df=1.0/(30.0*1024.0);print "# omega"}{printf "%16.14lf\n", $1*2.0*pi*df}' > omg_table.txt
#
# interpolate the phase dispersion
#
sample1d mode.dat_4km_yn -${HEADER_OPTION}1 -Fc -Nomg_table.txt -T1 ${NAN_OPTION} > sampled_table.txt
#
# compute the phase difference
#
awk 'BEGIN{gd=9.8231*4000.0;print 0,0.0,0.0,0.0;pi=4.0*atan2(1.0,1.0)} $1 !~/^#/ {c=(sqrt(gd)-$4)*$2/gd*9000000; while ( c < -2*pi ) c += 2*pi; while ( c > 2*pi ) c-=2*pi; print NR-1, $2,$4,c}' sampled_table.txt > phase_correction.txt
#
# apply phase coorection 
#
join  phase_correction.txt j2 | awk 'BEGIN{CONVFMT="%.7g"}{print NR-1,$5,$6-$4} END{for(i=NR;i<513;i++) print i, 0.0,0.0}' > j3
#
# prepare amp phase spectrum data
#
awk '{print $2, $3}' j3 > j4
#
# inverse DFT by fftw3 program, which converts from complex spectra to real time series
#
prg/dft_c2r < j4 | awk 'NR>4' | awk '{printf "%15.6e", $2;if (NR%5==0) printf "\n"}END{ if (NR%5) printf "\n"}'  > j5
#
# construct time series 
#
head -30 tsunami_4km_constv_9000km_t40000sec.cut.sac > j6
cat j6 j5 > tsunami_4km_constv_9000km_t40000sec.cut_correctedbyFFTW.sac
#
#
rm j1 j2 j3 j4 j5 j6 sampled_table.txt phase_correction.txt  omg_table.txt tsunami_4km_constv_9000km_t40000sec.cut.sac
#rm tsunami_4km_constv_9000km_t40000sec.cut_correctedbyFFTW.sac 
