#!/bin/bash
#
# make a single text file containing 2 data files.
#
paste -d " " \
tsunami_timeaxis.txt tsunami_4km_constv_9000km.txt  tsunami_4km_yn_9000km.txt |\
awk '{printf "%06d %+10.3e %+10.3e\n",$1,$2,$3}' > junk
#
cat <(echo "time 9000km prem_tsunami") junk > enshu_synthetic.txt 
#
rm junk
#
# create two SAC binary files.
#
if ! command -v sac $> /dev/null; then
	echo "SAC is not found in the PATH."
	exit
fi
#
sac << END
readtable header 1 content xyy. enshu_synthetic.txt
ch LEVEN TRUE DELTA 30.0
w tsunami_4km_constv_9000km.sac tsunami_4km_prem_9000km.sac
quit
END
#
rm enshu_synthetic.txt

