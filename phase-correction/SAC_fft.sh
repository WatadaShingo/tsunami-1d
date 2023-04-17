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
#gmtmath -T0/1023/1 -N1 = time.txt
sac << END
cut 40000 60000
r tsunami_4km_constv_9000km.sac
taper
w alpha tsunami_4km_constv_9000km.cut.sac
fft
w alpha tsunami_4km_constv_9000km_sp_t40000sec.sac
quit
END
#
head -235 tsunami_4km_constv_9000km_sp_t40000sec.sac > tsunami_4km_constv_9000km_sp_t40000sec_amp.txt
tail -205 tsunami_4km_constv_9000km_sp_t40000sec_amp.txt | awk '{for(i=1;i<=NF;i++) print $i}' | awk '{print NR-1,$1/30.0}' > amp.txt
tail -205 tsunami_4km_constv_9000km_sp_t40000sec.sac > enshu2_ph.txt
awk '{for(i=1;i<=NF;i++) print $i}' enshu2_ph.txt | awk '{print NR-1,$1}' > phase.txt
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
# correct the phase
#
join phase_correction.txt phase.txt | awk 'BEGIN{CONVFMT="%.7g"}{if( NR-1 < 512 ) print NR-1, $5-$4;else print NR-1, 0.0} END{for(i=NR;i<1024;i++) print i, 0.0}' > phase_corrected.txt
#
# prepare conjugate phase for 
#
sort -n -r phase_corrected.txt | awk 'BEGIN{print 0, 0.000000;CONVFMT="%.7g"} NR < 1024 {print NR, -1.000000*$2}' > jj
join jj phase_corrected.txt |  awk '{printf "%15.6f", $2+$3; if (NR%5==0) printf "\n"}END{if (NR%5) printf "\n"}' > tsunami_4km_constv_9000km_sp_t40000sec_ph_corrected.txt
#
# create SAC ascii spectrum file
#
cat tsunami_4km_constv_9000km_sp_t40000sec_amp.txt tsunami_4km_constv_9000km_sp_t40000sec_ph_corrected.txt > tsunami_4km_constv_9000km_sp_t40000sec_corrected.sac
#
# Inverse FFT by SAC
#
sac << END
r alpha tsunami_4km_constv_9000km_sp_t40000sec_corrected.sac
ifft
w alpha tsunami_4km_constv_9000km_t40000sec.cut_correctedbySAC.sac
quit
END
#
#
rm tsunami_4km_constv_9000km_sp_t40000sec.sac tsunami_4km_constv_9000km_sp_t40000sec_amp.txt enshu2_ph.txt phase.txt omg_table.txt sampled_table.txt phase_correction.txt phase_corrected.txt jj tsunami_4km_constv_9000km_sp_t40000sec_ph_corrected.txt tsunami_4km_constv_9000km_sp_t40000sec_corrected.sac amp.txt  tsunami_4km_constv_9000km.cut.sac
#rm tsunami_4km_constv_9000km_t40000sec.cut_correctedbySAC.sac
