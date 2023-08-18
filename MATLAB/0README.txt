A sample MATLAB scripts for dispersive Tsunami simulation 
See equation 1 in Watada et al. (2014).

Table of Contents

1. Files
2. Running Matlab scripts

1. Files

-rw-r--r--@ 1 shingo  staff   1917 Apr 17 12:46 0README.txt

-rw-r--r--@ 1 shingo  staff   5003 Apr 16 21:43 tsunami_4km_yn.m
-rw-r--r--@ 1 shingo  staff    785 Aug  6  2014 gauss.m
-rw-r--r--@ 1 shingo  staff    853 Nov 28  2011 ricker.m

-rw-r--r--  1 shingo  staff   6077 Aug 11  2012 mode.dat_4km_yn
-rw-r--r--  1 shingo  staff  23244 Apr 16 21:58 tsunami_4km_constv_3000km.txt
-rw-r--r--  1 shingo  staff  23244 Apr 16 21:58 tsunami_4km_constv_6000km.txt
-rw-r--r--  1 shingo  staff  23258 Apr 16 21:58 tsunami_4km_constv_9000km.txt
-rw-r--r--  1 shingo  staff  22979 Apr 16 21:58 tsunami_4km_yn_0km.txt
-rw-r--r--  1 shingo  staff  22931 Apr 16 21:58 tsunami_4km_yn_3000km.txt
-rw-r--r--  1 shingo  staff  22607 Apr 16 21:58 tsunami_4km_yn_6000km.txt
-rw-r--r--  1 shingo  staff  22373 Apr 16 21:58 tsunami_4km_yn_9000km.txt
-rw-r--r--  1 shingo  staff  11915 Apr 16 21:58 tsunami_timeaxis.txt

	# MATLAB scripts #
	tsunami.m:  a matlab function M-file for tsunami propagation simulation along 
		a given dispersion branch.
		Wave packet type, either ricker or gaussian can be selected.
	ricker.m: a matlab function M-file for a ricker wave packet
	gauss.m : a matlab function M-file for a gaussian wave packet

	# Input data file #
	mode.dat_4km_yn: theoretical tsunami mode dispersion relation for a 4km deep ocean of
		the PREM earth model including the effects of the spatio-temporal change in gravity,
		water compressibility, and elasticity of the Earth.
	mode.dat_4km_yn is an output from the mode-tsunami package that is available in GitHub.

2. Running MATLAB scripts

	# Run the tsunami propagation demonstration,
	# in the MATLAB command line prompt (>>), please type 
	>> tsunami_4km_yn ricker
	or
	>> tsunami_4km_yn gauss  
	##  NOTE for MATLAB beginners ##
	Please type without suffix ".m" in the command line.
	## END of NOTE for MATLAB beginners ##

	# Output files 
	tsunami_4km_yn_0km.txt: original tsunami waveform at 0 km
	tsunami_4km_yn_3000km.txt: dispersed waveform at 3000 km
	tsunami_4km_yn_6000km.txt: dispersed waveform at 6000 km
	tsunami_4km_yn_9000km.txt(*1): dispersed waveform at 9000 km 
	tsunami_constv_3000km.txt: non-dispersed waveform at 3000 km
	tsunami_constv_6000km.txt: non-dispersed waveform at 6000 km
	tsunami_constv_9000km.txt(*1): non-dispersed waveform at 9000 km
	tsunami_timeaxis.txt: time axis

	(*1) files were copied to the phase_correction directory

Reference
Watada, S., Kusumoto, S., and Satake, K. (2014), 
Traveltime delay and initial phase reversal of distant tsunamis coupled with the self-gravitating elastic Earth,
J. Geophys. Res. Solid Earth, 119, 4287– 4310, doi:10.1002/2013JB01084

Watada,S. (2023),
Progress and Application of the Synthesis of Trans-oceanic Tsunamis,
Progress in Earth and Planetary Science, 10, 26, doi:10.1186/s40645-023-00555-1
