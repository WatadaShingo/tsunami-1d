"""
Matplotlib animation for tsunami propagation over a constant depth ocean (4km)
PREM, surface gravity wave, long-wave simulations
animation is a series of sntapshots of tsunami waveforms
This is a conmbined version of wave_snapshots.py and wave_moving_snapshots.py

input:
        mode.dat_4km_yn
output:
        tsunami_animation_snapshots.mp4

% python wave_animation_snapshots.py < wave.py.in

Shingo Watada
2023/09/01
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
#
# conversion from distance to time is 200m/s
#
distance_unit=50e3 # 50 km distance unit from source location
dt=30.0
t0=dt*64 # 1920.0
time_unit=250 # 250 sec time unit from origin time
dx=6e3
x0=dx*64 # 384.0e3
time_axis=np.arange(total_npts)*dt
x_axis=np.arange(total_npts)*dx
#
# angular frequency domain spectra
source_waveform1=gauss(0,source_npts,center_freq1,dt,t0)
source_waveform2=gauss(0,source_npts,center_freq2,dt,t0)
domg=1.0/(total_npts*dt)*2*pi
omg=np.arange(total_npts)*domg
#
# wavenumber domain spectra
source_waveform1=gauss(0,source_npts,center_freq1/200,dx,x0)
source_waveform2=gauss(0,source_npts,center_freq2/200,dx,x0)
dk=1.0/(total_npts*dx)*2*pi
k=np.arange(total_npts)*dk
#
y1=np.zeros(total_npts)
y2=np.zeros(total_npts)
vp=np.zeros(total_npts)
y1[0:source_npts]=source_waveform1
y2[0:source_npts]=source_waveform2

# omg[total_npts/2]: nyqiest angular frequency

el_table,omg_table,vp_table=set_dispersion(dispersion_file)
k_table=omg_table/vp_table
nl_table=el_table.size
print('type of k_table[]=', type(k_table),'size=',k_table.size )
print('type of k[]=', type(k),'size=',k.size )
print(f'number of l= {nl_table} first l={el_table[0]} last l={el_table[-1]}')
#
# check the minimum and maximum frequencies are within the dispersion table
#
"""
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
    print('maximum freq in the table = %.6e' % omg_table[-1])
    print('nyquist freq = %.6f' % omg[total_npts//2])
    #
    # find i omg[i-1]<omg_table[-1]<=omg[i]
    #
    max_omg_i=np.searchsorted(omg,omg_table[-1])
"""
#
# check the minimum and maximum wavenumbers are within the dispersion table
#
min_k_i=0

if k[1] < k_table[0]:
    print('minmum wavenumber the table = %.6e' % k_table[0])
    print('k = %.6f' % k[1])
    sys.exit(0)
    #
    # find i k[i-1]<k_table[0]<=k[i]
    #
    min_k_i=np.searchsorted(k,k_table[0])

max_k_i=total_npts//2

if k_table[-1] < k[total_npts//2]:
    print('maximum wavenumber in the table = %.6e' % k_table[-1])
    print('nyquist wavenumber = %.6f' % k[total_npts//2])
    #
    # find i k[i-1]<k_table[-1]<=k[i]
    #
    max_k_i=np.searchsorted(k,k_table[-1])

#
# amplitude spectrum of k[i] is set to zero for all max_k_i <= i
#

print(f'min_k_i={min_k_i}, k[min_k_i]={k[min_k_i]:e}')
print(f'max_k_i={max_k_i}, k[max_k_i]={k[max_k_i]:e}')

v_phase=np.zeros(total_npts)
v_phase[0]=sqrt(water_depth*g0)

#for i in range(1,total_npts):
#   v_phase[i]=vp_water(water_depth,omg[i])

for i in range(1,total_npts):
   v_phase[i]=sqrt(g0/k[i]*tanh(k[i]*water_depth))


vp[min_k_i:max_k_i]=CubicSpline(k_table,vp_table)(k[min_k_i:max_k_i])
vp[:min_k_i]=1
vp[max_k_i:]=1

print(f'type of vp: {type(vp)} size of vp: {vp.size}.')

#
# size of s0 is total_npts//2+1
#
s1=np.fft.rfft(y1)
s2=np.fft.rfft(y2)
print(f'type of s1[0]: {type(s1[0])} size of s1: {s1.size}.')
# size of omg, k, vp, phi0, phi1 is total_npts
#

# First set up the figure, the axis, and the plot element we want to animate
#fig = plt.figure()
#fig,ax=plt.subplots(figsize=(6.0,4.0))
fig=plt.figure(figsize=(6.0,4.0))
gs=fig.add_gridspec(1,2,wspace=0)
axmoving,axfixed=gs.subplots(sharey=True)
#
# Both ax=plt.axes() ax=fig.add_subplot() are OK
#ax = plt.axes(xlim=(0, 2), ylim=(-4, 2))
#
#ax = fig.add_subplot(xlim=(0, 2), ylim=(-4, 2))
#ax = fig.add_subplot(xlim=(0, 960), ylim=(-1.7, 0.8))
#ax = fig.add_subplot(xlim=(0, 12288), ylim=(-1.7, 0.8))
axmoving.set(xlim=(0, 12288), ylim=(-1.7, 0.8),autoscale_on=False)
axfixed.set(xlim=(0, 12288), ylim=(-1.7, 0.8))
# ax.plot return a line object

line1, = axmoving.plot([], [],'dodgerblue', lw=1)
line2, = axmoving.plot([], [],'lightsalmon', lw=1)
line3, = axmoving.plot([], [],'gray', lw=1)
line4, = axfixed.plot([], [],'dodgerblue', lw=1)
line5, = axfixed.plot([], [],'lightsalmon', lw=1)
line6, = axfixed.plot([], [],'gray', lw=1)
line11, = axmoving.plot([], [],'b', lw=1)
line12, = axmoving.plot([], [],'r', lw=1)
line13, = axmoving.plot([], [],'k', lw=1)
line14, = axfixed.plot([], [],'b', lw=1)
line15, = axfixed.plot([], [],'r', lw=1)
line16, = axfixed.plot([], [],'k', lw=1)
line1.set_label('line1')
line2.set_label('line2')
line3.set_label('line3')
line4.set_label('line4')
line5.set_label('line5')
line6.set_label('line6')
#
# setup legend
#
first_legend=axfixed.legend( [line11, line12,line13], \
['PREM','gravity wave','linear long-wave'], \
bbox_to_anchor=(0.0,1.02),ncol=3, loc='upper center', \
framealpha=1 )

#
time_template = 'at t=%6.1f min'
#time_text = ax.text(2000.0,0.5,'',transform=ax.transData)
time_text = axmoving.text(0.1,0.85,'',transform=axmoving.transAxes)

# initialization function: plot the background of each frame
# Either one of the following return lines works interchangeablly.
# return line1,line2,line3,line4,line5,line6
# return line1,line2,line3,line4,line5,line6,
# return [line1,line2,line3,line4,line5,line6]
def init():
    line1.set_data([], [])
    line2.set_data([], [])
    line3.set_data([], [])
    line4.set_data([], [])
    line5.set_data([], [])
    line6.set_data([], [])
    line11.set_data([], [])
    line12.set_data([], [])
    line13.set_data([], [])
    line14.set_data([], [])
    line15.set_data([], [])
    line16.set_data([], [])
    return line1,line2,line3,line4,line5,line6, line11,line12,line13,line14,line15,line16

# animation function.  This is called sequentially
# Either one of the following return lines works interchangeably.
# return line1,line2,line3,line4,line5,line6
# return line1,line2,line3,line4,line5,line6,
# return [line1,line2,line3,line4,line5,line6]
def plotlines(i):
    ''' plot three lines at each frame '''
    s=s1*exp(-time_unit*i*vp*k*1j)[0:total_npts//2+1]
    s[total_npts//2]=0
    line1.set_data(x_axis/1000, np.fft.irfft(s))
    line4.set_data(x_axis/1000, np.fft.irfft(s))
    s=s2*exp(-time_unit*i*vp*k*1j)[0:total_npts//2+1]
    s[total_npts//2]=0
    line11.set_data(x_axis/1000, np.fft.irfft(s))
    line14.set_data(x_axis/1000, np.fft.irfft(s))

    s=s1*exp(-time_unit*i*v_phase*k*1j)[0:total_npts//2+1]
    s[total_npts//2]=0
    line2.set_data(x_axis/1000, np.fft.irfft(s)-0.75)
    line5.set_data(x_axis/1000, np.fft.irfft(s)-0.75)
    s=s2*exp(-time_unit*i*v_phase*k*1j)[0:total_npts//2+1]
    s[total_npts//2]=0
    line12.set_data(x_axis/1000, np.fft.irfft(s)-0.75)
    line15.set_data(x_axis/1000, np.fft.irfft(s)-0.75)

    s=s1*exp(-time_unit*i*sqrt(water_depth*g0)*k*1j)[0:total_npts//2+1]
    s[total_npts//2]=0
    line3.set_data(x_axis/1000, np.fft.irfft(s)-1.5)
    line6.set_data(x_axis/1000, np.fft.irfft(s)-1.5)
    s=s2*exp(-time_unit*i*sqrt(water_depth*g0)*k*1j)[0:total_npts//2+1]
    s[total_npts//2]=0
    line13.set_data(x_axis/1000, np.fft.irfft(s)-1.5)
    line16.set_data(x_axis/1000, np.fft.irfft(s)-1.5)
# set text in text instance 'time_text'
#
    time_text.set_text(time_template % (time_unit*i/60))
#
    axmoving.set(xlim=(200*i*time_unit/1000-1000,200*i*time_unit/1000+600))
    fig.canvas.draw()
    return [line1,line2,line3,line4,line5,line6]

fig.suptitle('tsunami simulation for 4km deep ocean,\ninitial waveform $exp(-x^2*k_a^2), k_a=2\pi/\lambda_a, \lambda_a$=133, 533 km')
fig.text(0.5,0.01,'distance, km',ha='center')


# call the animator.  blit=True means only re-draw the parts that have changed.
# interval is milli-second between each frame
anim = animation.FuncAnimation(fig, plotlines, init_func=init,\
 frames=200, interval=60, blit=True)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
anim.save('tsunami_animation_snapshots.mp4', fps=10, extra_args=['-vcodec', 'libx264'])

# display a movie in a Matplotlib figure
#plt.show()

print( 'clock time=%.5f' sec % (time.time()-t))
