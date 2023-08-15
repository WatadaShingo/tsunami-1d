"""
long-wave and disperive tsunami simulations
input:
	mode.dat_4km_yn
outputs:
	tsunami_4km_yn_0km.txt
	tsunami_4km_yn_3000km.txt
	tsunami_4km_yn_6000km.txt
	tsunami_4km_yn_9000km.txt
	tsunami_4km_constv_3000km.txt
	tsunami_4km_constv_6000km.txt
	tsunami_4km_constv_9000km.txt
	tsunami_timeaxis.txt
"""

import sys
import numpy as np
from numpy import sqrt, pi, exp
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt

dispersion_file='mode.dat_4km_yn'
#dispersion_file=input('dispersion table file name >>')
#print('')

def set_dispersion(filename):
   data=np.loadtxt(filename,dtype='float64',skiprows=1)
   el_table=data[:,0]
   omg_table=data[:,1]
   period_table=data[:,2]
   vp_table=data[:,3]
   vg_table=data[:,4]
   return data[:,0],data[:,1],data[:,3]

def ricker(i,npts,freq,dt,timesh):
   '''
   Ricker waveform 
   i: time index t=i*dt
   npts: number of data points
   freq: predominant frequency of ricker waveform in Hz
   dt: sampling interval in second
   timesh: timeshift of the center of ricker waveform in second
   '''

   pi2=sqrt(pi)/2
   b=sqrt(6)/(pi*freq)
   const=2*sqrt(6)/b
   t1=np.arange(npts)*dt
   t2=t1-timesh-i*dt
   u=const*t2
   s=(u*u/4 - 0.5)*pi2*exp(-u*u/4)
   return s

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

source_npts=128
total_npts=2048
center_freq=0.0015
distance=3000e3 # 9000 km distance from source
dt=30.0
t0=1920.0
time=np.arange(total_npts)*dt
source_waveform=gauss(0,source_npts,center_freq,dt,t0)
#source_waveform=-ricker(0,source_npts,center_freq,dt,t0)
y=np.zeros(total_npts)
vp=np.zeros(total_npts)
y[0:source_npts]=source_waveform
domg=1.0/(total_npts*dt)*2*pi
omg=np.arange(total_npts)*domg
# omg[total_npts/2]: nyquist angular frequency

el_table,omg_table,vp_table=set_dispersion(dispersion_file)
nl_table=el_table.size
#
# check the minimum and maximum frequencies are within the dispersion table 
#
min_omg_i=0

if omg[1] < omg_table[0]:
    print('minmum freq the table = %.6e' % omg_table[0])
    print('domg = %.6f' % omg[1])
    sys.exit(0)
    #
    # find i omg[i-1]<omg_table[0]<=omg[i]
    #
    min_omg_i=np.searchsorted(omg,omg_table[0])

max_omg_i=total_npts//2

if omg_table[-1] < omg[total_npts//2]:
    print('maximum freq the table = %.6e' % omg_table[-1])
    print('nyquist freq = %.6f' % omg[total_npts//2])
    #
    # find i omg[i-1]<omg_table[-1]<=omg[i]
    #
    max_omg_i=np.searchsorted(omg,omg_table[-1])

#
# amplitude spectrum of omg[i] is set to zero for all max_omg_i <= i
#

vp[min_omg_i:max_omg_i]=CubicSpline(omg_table,vp_table)(omg[min_omg_i:max_omg_i])
print(f'type of vp: {type(vp)} size of vp: {vp.size}.')
vp[:min_omg_i]=1
vp[max_omg_i:]=1

#
# size of s0 is total_npts//2+1
#
s0=np.fft.rfft(y)
# size of omg, vp, phi0, phi1 is total_npts
phi_disperse=distance*omg/vp
phi_const   =distance*omg/sqrt(4000*9.8231)
#
# do not change phase for i=0, total_npts//2
#
phi_disperse[0]=0
phi_const[0]     =0
x_disperse=exp(-phi_disperse*1j)
x_const   =exp(-phi_const*1j)
x_disperse[total_npts//2]=0
x_const[total_npts//2]   =0
#
# set zero amplitude above i=max_omg_i
#
x_disperse[max_omg_i:]=0
x_const[max_omg_i:]   =0
#
# 
#
s1=s0*x_disperse[0:total_npts//2+1]
y1=np.fft.irfft(s1)
s2=s1*x_disperse[0:total_npts//2+1]
y2=np.fft.irfft(s2)
s3=s2*x_disperse[0:total_npts//2+1]
y3=np.fft.irfft(s3)


# f.write does not write '\n' at the end of the file
#l_y=[format(f,'e') for f in y]
#
#with open('tsunami_4km_yn_0km.txt',mode='x') as f:
#    f.write('\n'.join(l_y))

with open('tsunami_timeaxis.txt',mode='x') as f:
    print(*time,sep='\n',file=f)

with open('tsunami_4km_yn_0km.txt',mode='x') as f:
    print(*y,sep='\n',file=f)
with open('tsunami_4km_yn_3000km.txt',mode='x') as f:
    print(*y1,sep='\n',file=f)
with open('tsunami_4km_yn_6000km.txt',mode='x') as f:
    print(*y2,sep='\n',file=f)
with open('tsunami_4km_yn_9000km.txt',mode='x') as f:
    print(*y3,sep='\n',file=f)

fig = plt.figure()
ax = plt.axes(xlim=(0, 960), ylim=(-0.2, 0.6))
plt.xticks(np.arange(0,1020,step=60))
ax.set_title('dispersive tsunami')
ax.set_xlabel('time, min')
ax.set_ylabel('amplitude')
line1,=ax.plot(time/60,y,'k',lw=2)
line2,=ax.plot(time/60,y1,'dimgray',lw=2)
line3,=ax.plot(time/60,y2,'gray',lw=2)
line4,=ax.plot(time/60,y3,'lightgray',lw=2)

ax.legend([line1,line2,line3,line4],['0 km', '3000 km','6000 km','9000 km'],loc='upper right')
plt.show()

s1=s0*x_const[0:total_npts//2+1]
y1=np.fft.irfft(s1)
s2=s1*x_const[0:total_npts//2+1]
y2=np.fft.irfft(s2)
s3=s2*x_const[0:total_npts//2+1]
y3=np.fft.irfft(s3)

with open('tsunami_4km_constv_3000km.txt',mode='x') as f:
    print(*y1,sep='\n',file=f)
with open('tsunami_4km_constv_6000km.txt',mode='x') as f:
    print(*y2,sep='\n',file=f)
with open('tsunami_4km_constv_9000km.txt',mode='x') as f:
    print(*y3,sep='\n',file=f)

fig = plt.figure()
ax = plt.axes(xlim=(0, 960), ylim=(-0.2, 0.6))
plt.xticks(np.arange(0,1020,step=60))
ax.set_title('long-wave tsunami')
ax.set_xlabel('time, min')
ax.set_ylabel('amplitude')
ax.plot(time/60,y,'k',lw=2)
ax.plot(time/60,y1,'dimgray',lw=2)
ax.plot(time/60,y2,'gray',lw=2)
ax.plot(time/60,y3,'lightgray',lw=2)
ax.legend([line1,line2,line3,line4],['0 km', '3000 km','6000 km','9000 km'],loc='upper right')
plt.show()
