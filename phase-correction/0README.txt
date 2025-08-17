This directory constains sample data files and programs 
to learn tsunami waveform calculation and phase corrections.

Table of Contents

1. Requirements
2. Files
3. Running programs
4. Compare phase-corrected files
5. Final note

1. Requirements 
	Sample programs are Bash scripts that call other utilities and programs.
	The utilities are single-function and can be replaced by other 
	tools or written in other programming languages, as you like.

	1) SAC (Seismic Analysis Code, http://ds.iris.edu/ds/nodes/dmc/software/downloads/sac/)

	SAC commands are used in convert2sac.sh, FFTW.sh, and SAC_fft.sh

	2) Interpolation of table data 

	A Generic Mapping Tools (GMT, https://www.generic-mapping-tools.org) command,
	sample1d, is used in FFTW.sh and SAC_fft.sh to interpolate tabled data.
	GMT Version 4 or 5.  Version 6 is not supported yet.

	3) Discrete Fourier Transform (DFT)

	The phase of time series data is manipulated in the frequency domain.
	The following two methods are provided. Either one works.

	3-1) SAC (Seismic Analysis Code,http://ds.iris.edu/ds/nodes/dmc/software/downloads/sac/)
	SAC commands are called to perform DFT in convert2sac.sh, FFTW.sh, and SAC_fft.sh

	3-2) Intel MKL (Math Kernel Library) FFTW functions are called to perform
	DFT in prg/dft_r2c.c and prg/dft_2r.c.

2. Files

-rw-r--r-- 1 watada watada  6839 Aug  6 22:36 0README.txt
-rwxr--r-- 1 watada watada  2422 Apr 17 11:56 FFTW.sh
-rwxr--r-- 1 watada watada  3139 Apr 17 10:28 SAC_fft.sh
-rwxr-xr-x 1 watada watada   516 Apr 17 11:43 convert2sac.sh

prg:
-rw-r--r-- 1 watada watada    145 Aug  7  2014 data_c.in
-rw-r--r-- 1 watada watada     64 Aug  6  2014 data_r.in
-rwxrwxr-x 1 watada watada 306520 Apr  6 22:29 dft_c2r
-rwxr--r-- 1 watada watdaa   1721 Jul 30  2016 dft_c2r.c
-rwxrwxr-x 1 watada watada 311848 Apr  6 22:29 dft_r2c
-rwxr--r-- 1 watada watada   1707 Jul 30  2016 dft_r2c.c
-rwxrwxr-x 1 watada watada 279552 Apr 17 12:51 fftw3-test
-rw-r--r-- 1 watada watada    671 Jun  5  2015 fftw3-test.c
-rw-r--r-- 1 watada watada    459 Apr  6 15:59 makefile

inputs: copied from either the MATLAB or the Python directory
-rw-r--r-- 1 watada watada 23258 Apr 16 21:58 tsunami_4km_constv_9000km.txt
-rw-r--r-- 1 watada watada 22373 Apr 16 21:58 tsunami_4km_yn_9000km.txt
-rw-r--r-- 1 watada watada 11915 Apr 16 21:58 tsunami_timeaxis.txt

inputs: created by running solid_rk.c 
-rw-r--r-- 1 watada watada  6077 May  2  2022 mode.dat_4km_yn

outputs: created by running convert2sac.sh
-rw-rw-rw- 1 watada watada  8824 Apr 17 11:43 tsunami_4km_constv_9000km.sac
-rw-rw-rw- 1 watada watada  8824 Apr 17 11:43 tsunami_4km_prem_9000km.sac

output: created by running FFTW.sh
-rw-rw-r-- 1 watada watada 17237 Apr 17 11:56 tsunami_4km_constv_9000km_t40000sec.cut_correctedbyFFTW.sac

output: created by running SAC_fft.sh
-rw-rw-r-- 1 watada watada 21980 Apr 17 10:57 tsunami_4km_constv_9000km_t40000sec.cut_correctedbySAC.sac

	# input data files #
	tsunami_4km_constv_9000km.txt (*1): non-dispersed tsunami waveform at 9000 km
		to which phase correction is applied 
	tsunami_4km_yn_9000km.txt (*1): dispersed tsunami waveform at 9000 km
	tsunami_timeaxis.txt (*1): time axis

	(*1) copied from either the MATLAB or the Python directory

	mode.dat_4km_yn (*2): theoretical tsunami mode dispersion relation for a 4km deep ocean of
		the PREM earth model including the effects of the spatio-temporal change in gravity,
      		water compressibility, and elasticity of the Earth.

	(*2) an output from the mode-tsunami package available in GitHub.

	# output data files #
	tsunami_4km_constv_9000km_t40000sec.cut_correctedbyFFTW.sac (*3): phase-corrected tsunami waveform
	tsunami_4km_constv_9000km_t40000sec.cut_correctedbySAC.sac (*3): phase-corrected tsunami waveform

	(*3) in SAC ascii format

	tsunami_4km_constv_9000km.sac (*4): conveted from tsunami_4km_constv_9000km.txt
	tsunami_4km_prem_9000km.sac (*4): converted from tsunami_4km_yn_9000km.txt

	(*4) in SAC binary format

3. Running programs

	Compute phase-corrected tsunami using the method proposed by Watada et al. (2014)

	# Compile C programs #
	# forward and inverse FFT
	# sample programs to test FFTW package are under prg directory
	% cd prg
	% make

	prg/dft_r2c converts real time series to complex spectrum 
	prg/dft_c2r converts complex spectrum to real time 

	# Tsunami wavefom data files were copied from the MALTAB directory
	tsunami_4km_constv_9000km.txt: 
	tsunami_4km_yn_9000km.txt
	tsunami_timeaxis.txt

	# Bash shell scripts to compute phase-corrected waveforms
	# convert *.txt files to SAC binary files
	% convert2sac.sh

	# Gaussian waveform tsunami_4km_constv_9000km.sac is windowed and phase correction is applied
	# using SAC package for forward and inverse FFT, 
	# and GMT command sample1d for interpolation of the phase velocity table
	% SAC_fft.sh

	# Another sample shell script to do forward and inverse FFT
	# by using the FFTW functions in the Intel MKL library
	% FFTW.sh


4. Compare phase-corrected files

	# Compare the two resulting phase corrected SAC ascii files
	# tsunami_4km_constv_9000km_t40000sec.cut_correctedbyFFTW.sac and
	# tsunami_4km_constv_9000km_t40000sec.cut_correctedbySAC.sac,
	# which are originally tsunami_4km_constv_9000km.sac before phase correction,
	# with PREM synthetic in the SAC binary format
	# tsunami_4km_prem_9000km.sac
	% sac
	SAC> r alpha tsunami_4km_constv*.sac
	SAC> r more tsunami_4km_prem_9000km.sac
	SAC> bd x;qdp off;p1

	Two phase-correction methods (SAC_fft.sh and FFTW.sh) provide nearly the same
	dispersed dwaveforms. tsunami_4km_prem_9000km.sac and two phase-corrected files
	are very close at the begining.  However, toward the end of the dispersed small
	amplitude waveforms, the phase shift between the two increases.
	This is because, in the later part of the dispersed packet, the waves are
	propagating slowly and their phase velocity deviates much from
	the long-wave speed. The phase correction method by Watada et al. (2014)
	assumes that the phase-speed change is much smaller than the long-wave velocity.

4. Final note
	The programs and scripts are provided for demonstration purposes only.
	To apply the phase correction to other time series data, the FFTW.sh or SAC_fft.sh scripts
	must be modified to incorporate the new sampling interval and number of data points.
	The basic preprocessing steps required for real time series data analysis
	have been omitted.

References
Watada, S., Kusumoto, S., and Satake, K. (2014), 
Traveltime delay and initial phase reversal of distant tsunamis coupled with the self-gravitating elastic Earth,
J. Geophys. Res. Solid Earth, 119, 4287– 4310, doi:10.1002/2013JB01084

Watada,S. (2023),
Progress and Application of the Synthesis of Trans-oceanic Tsunamis,
Progress in Earth and Planetary Science, 10, 26, doi:10.1186/s40645-023-00555-1
