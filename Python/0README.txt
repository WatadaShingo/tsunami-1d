A sample Python script for dispersive Tsunami simulation 
See equation 1 in Watada et al. (2014).

Table of Contents

1. Files
2. Running Python scripts

1. Files

-rw-r--r--. 1 watada ftp  2954  8月 15 21:48 2023 0README.txt

-rw-r--r--. 1 watada ftp  5561  8月 15 21:34 2023 tsunami_4km_yn.py

-rw-r--r--. 1 watada ftp  6077  8月 11 16:53 2012 mode.dat_4km_yn

-rw-r--r--. 1 watada ftp 46531  8月 15 21:34 2023 tsunami_4km_constv_3000km.txt
-rw-r--r--. 1 watada ftp 46499  8月 15 21:34 2023 tsunami_4km_constv_6000km.txt
-rw-r--r--. 1 watada ftp 46511  8月 15 21:34 2023 tsunami_4km_constv_9000km.txt
-rw-r--r--. 1 watada ftp 10499  8月 15 21:34 2023 tsunami_4km_yn_0km.txt
-rw-r--r--. 1 watada ftp 46519  8月 15 21:34 2023 tsunami_4km_yn_3000km.txt
-rw-r--r--. 1 watada ftp 46185  8月 15 21:34 2023 tsunami_4km_yn_6000km.txt
-rw-r--r--. 1 watada ftp 45900  8月 15 21:34 2023 tsunami_4km_yn_9000km.txt
-rw-r--r--. 1 watada ftp 16011  8月 15 21:34 2023 tsunami_timeaxis.txt

	# Python scripts #
	tsunami.py:  a python script for tsunami propagation simulation along 
		a given dispersion branch.
		Wave packet type, either ricker or gaussian can be selected.

	# Input data file #
	mode.dat_4km_yn: theoretical tsunami mode dispersion relation for a 4km deep ocean of
		the PREM earth model including the effects of the spatio-temporal change in gravity,
		water compressibility, and elasticity of the Earth.
	mode.dat_4km_yn is an output from the mode-tsunami package that is available in GitHub.

2. Running Python scripts

	# Remove all existing *.txt files. Python script stops if they exist.
	% rm *.txt

	# Run the tsunami propagation demonstration,
	% python tsunami_4km_yn.py

	# Output files 
	tsunami_4km_yn_0km.txt: original tsunami waveform at 0 km
	tsunami_4km_yn_3000km.txt: dispersed waveform at 3000 km
	tsunami_4km_yn_6000km.txt: dispersed waveform at 6000 km
	tsunami_4km_yn_9000km.txt: dispersed waveform at 9000 km 
	tsunami_constv_3000km.txt: non-dispersed waveform at 3000 km
	tsunami_constv_6000km.txt: non-dispersed waveform at 6000 km
	tsunami_constv_9000km.txt: non-dispersed waveform at 9000 km
	tsunami_timeaxis.txt: time axis

Reference
Watada, S., Kusumoto, S., and Satake, K. (2014), 
Traveltime delay and initial phase reversal of distant tsunamis coupled with the self-gravitating elastic Earth,
J. Geophys. Res. Solid Earth, 119, 4287– 4310, doi:10.1002/2013JB01084

Watada,S.,
Progress and Application of the Synthesis of Trans-oceanic Tsunamis,
Progress in Earth and Planetary Science, submitted
