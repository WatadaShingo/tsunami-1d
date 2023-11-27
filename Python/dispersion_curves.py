"""
Matplotlib tsunami dispersion curves for a constant depth ocean (4km)
PREM, surface gravity wave, long-wave simulations
and amplitude spectra of initial gaussian waveforms

input:
        mode.dat_4km_yn
output:
	dispersion_curves.jpg
        

% python dispersion_curves.py < wave.py.in

Shingo Watada
2023/09/06
"""
import sys
import time
import numpy as np
from numpy import sqrt,inf, pi, cos,tanh,cosh,sinh,exp,abs
import scipy as sp
from scipy.interpolate import CubicSpline
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import animation

print("python = ",sys.version)
print("numpy = ",np.version.full_version)
print("scipy = ",sp.version.full_version)
print("matplotlib = ",matplotlib.__version__)
t=time.time()

dispersion_file=input('dispersion table file name >>')
print('')

def set_dispersion(filename):
   data=np.loadtxt(filename,dtype='float64',skiprows=1)
   el_table=data[:,0]
   omg_table=data[:,1]
   period_table=data[:,2]
   vp_table=data[:,3]
   vg_table=data[:,4]
   return data[:,0],data[:,1],data[:,3]


def gauss(i,npts,freq,dt,timesh):
   '''
   gausian waveform
   i: time index t=i*dt
   npts: number of data points
   freq: central frequency of gaussian waveform in Hz
   dt: sampling interval in second
   timesh: timeshift of the center of gaussian in second
   '''

   const=2*pi*freq
   t1=np.arange(npts)*dt
   t2=const*(t1-timesh-i*dt)
   s=0.5*exp(-t2*t2)
   return s

def vp_water(depth,omg):
    '''
    compute phase velocity of water wave
    depth: ocean depth in meter
    omg: angular frequency of water waves 
    '''
    k0=0
    c1=sqrt(g0*depth)
    k1=omg/c1
    while abs(k1-k0)/k1 > 1e-6 :
       c0=c1
       k0=k1
       c1=sqrt(g0/k0*tanh(k0*depth))
       k1=omg/c1
#
    return c1


source_npts=128
total_npts=2048
center_freq1=0.0015
center_freq2=0.000375
water_depth=4000
g0=9.8231
distance=50e3 # 50 km distance unit from source
dt=30.0
t0=1920.0
time_axis=np.arange(total_npts)*dt
#time_source=np.arange(source_npts)*dt
source_waveform1=gauss(0,source_npts,center_freq1,dt,t0)
source_waveform2=gauss(0,source_npts,center_freq2,dt,t0)
y1=np.zeros(total_npts)
y2=np.zeros(total_npts)
vp=np.zeros(total_npts)
y1[0:source_npts]=source_waveform1
y2[0:source_npts]=source_waveform2
domg=1.0/(total_npts*dt)*2*pi
omg=np.arange(total_npts)*domg
# omg[total_npts/2]: nyqiest angular frequency

el_table,omg_table,vp_table=set_dispersion(dispersion_file)
nl_table=el_table.size
print(f'number of l= {nl_table} first l={el_table[0]} last l={el_table[-1]}')
#
# check the minimum and maximum frequencies are within the dispersion table
#
min_omg_i=0

if omg[1] < omg_table[0]:
    print('minmum freq in the table = %.6e' % omg_table[0])
    print('domg = %.6f' % omg[1])
    sys.exit(0)
    #
    # find i omg[i-1]<omg_table[0]<=omg[i]
    #
    min_omg_i=np.searchsorted(omg,omg_table[0])

max_omg_i=total_npts//2

if omg_table[-1] < omg[total_npts//2]:
    print('maximum freq in the table = %.6e' % omg_table[-1])
    print('nyquist freq = %.6f' % omg[total_npts//2])
    #
    # find i omg[i-1]<omg_table[-1]<=omg[i]
    #
    max_omg_i=np.searchsorted(omg,omg_table[-1])

#
# amplitude spectrum of omg[i] is set to zero for all max_omg_i <= i
#

print(f'min_omg_i={min_omg_i}, omg[min_omg_i]={omg[min_omg_i]:e}')
print(f'max_omg_i={max_omg_i}, omg[max_omg_i]={omg[max_omg_i]:e}')

v_const=np.full(total_npts,sqrt(water_depth*g0))

v_phase=np.zeros(total_npts)
v_phase[0]=sqrt(water_depth*g0)

for i in range(1,total_npts):
   v_phase[i]=vp_water(water_depth,omg[i])

vp[min_omg_i:max_omg_i]=CubicSpline(omg_table,vp_table)(omg[min_omg_i:max_omg_i])
print(f'type of vp: {type(vp)} size of vp: {vp.size}.')
vp[:min_omg_i]=1
vp[max_omg_i:]=1

#
# size of s1 is total_npts//2+1
#
s1=np.fft.rfft(y1)
s2=np.fft.rfft(y2)
print(f'type of s1[0]: {type(s1[0])} size of s1: {s1.size}.')
# size of omg, vp, phi0, phi1 is total_npts
#

# First set up the figure, the axis, and the plot element we want to animate
#fig = plt.figure(figsize=(6.0,4.0))
fig = plt.figure(figsize=(3.0,4.0))

#
# Both ax=plt.axes() ax=fig.add_subplot() are OK
#ax = plt.axes(xlim=(0, 2), ylim=(-4, 2))
#
#ax = fig.add_subplot(xlim=(0, 2), ylim=(-4, 2))
#ax = fig.add_subplot(xlim=(0, 960), ylim=(-1.7, 0.8))
# ax.plot return a line object

#ax=plt.axes(xlim=(2e-5,1e-2),ylim=(100,220))
ax=fig.add_subplot(2,1,1,xlim=(2e-5,1e-2),ylim=(150,200))
line1,=ax.semilogx(omg/(2*pi),vp,'blue',lw=1)
line2,=ax.semilogx(omg/(2*pi),v_phase,'red',lw=1)
line3,=ax.semilogx(omg/(2*pi),v_const,'black',lw=1)

ax.legend([line1, line2,line3], \
['PREM','gravity wave','linear long-wave'], \
bbox_to_anchor=(0.05,0.0),ncol=1, loc='lower left', \
bbox_transform=ax.transAxes)
ax.set_title('Tsunami dispersion relation')
ax.set_ylabel('phase velocity (m/s)')
ax.set_yticks([160,180,200],minor=False)
ax.set_yticks([150,170,190],minor=True)
ax.grid()

ax=fig.add_subplot(2,1,2,xlim=(2e-5,1e-2),ylim=(0,15))
line1,=ax.semilogx(omg[0:total_npts//2+1]/(2*pi),abs(s1),'gray',lw=2)
line2,=ax.semilogx(omg[0:total_npts//2+1]/(2*pi),abs(s2),'black',lw=2)
ax.set_xlabel('frequency (Hz)')
ax.legend([line1,line2],['T=677 sec', 'T=2667 sec'],ncol=1,\
loc='center',bbox_to_anchor=(0.85,0.5))

ax.set_title('Amplitude spectra of initial waveforms')
ax.set_ylabel('spectral amplitude')


ax.grid()
fig.tight_layout()
plt.show()
fig.savefig("dispersion_curves.jpg",dpi=360)

print( 'clock time=%.5f' sec % (time.time()-t))

sys.exit(0)
