function tsunami(wavetype)
% modofied to check time shift 
% modified from wavellet fft program by Braile. S. Watada 12/01/2011
% Test FFT with wavelet
% TestFFT.m  L Braile 4/11/05
npts = 128;   % Number of points in wavelet
dt = 30;    % sample interval
freq =1.5e-3;     % set ~peak frequency of wavelet
timesh = 0.0; % allow timeshift to center the wavelet
nsw = 0;     % apply cosine bell taper to ends of signal

switch wavetype
    case 'ricker' % calculate ricker wavelet
        s = ricker(npts,freq,dt,timesh,nsw);
	wave_type='Ricker';
    case 'gauss'% calculate gaussian wavelet
        s = gauss(npts,freq,dt,timesh,nsw); 
	wave_type='Gaussian';
    otherwise
        disp('wavetype is not known')
        return
end
t = 0:dt:(npts*dt - dt);

figure
plot(t,s,'-r','linewidth',1.5)
set(gca,'fontsize',16,'linewidth',2)
xlabel('Time (s)','fontsize',16)
ylabel('Amplitude','fontsize',16)
header_line=sprintf('%s Wavelet',wave_type);
title(header_line,'fontsize',16)

mode=importdata('mode.dat_4km_yn');
omg0=mode.data(:,2).';
el=mode.data(:,1).';
vp0=mode.data(:,4).';
vg0=mode.data(:,5).';
[dum size_table]=size(omg0);

%
nf = 2048; % set length of fft
S0 = fft(s,nf); % calculate FFT of wavelet, the added zeros (nf = 1024)
y= ifft(S0,nf);
tt=0:dt:dt*(nf-1);
figure
plot(tt,y);
set(gca,'fontsize',16,'linewidth',2)
xlabel('Time (sec)','fontsize',16)
ylabel('Linear amplitude','fontsize',16)
title('Tsunami at 0 km','fontsize',16)

distance=3000e3; % at 3000 km
%
% constant phase velocity problem
% 
vp=sqrt(9.8231*4000); % constant tsunami phase velocity m/s of ocean 4km
% 9.8231 m/s^2 is the gravity value at the 4km deep ocean bottom of the PREM earth model
tshift=distance/vp;
w=tshift/dt/nf;

%
% phase velocity with dispersion problem
%
% get phase velocity from input array of dispersion diagram
domg=1/(nf*dt)*2*pi;
omg=(0:nf-1)*domg;
vp=zeros(1,nf-1);
vp=interp1(omg0,vp0,omg,'spline');

for k=1:nf
 if (omg(k)>omg0(size_table))
vp(k)=1.0;
 end
end
     vp(1)=1.0;

% compute delay phase for each frequency
delaytime=distance./vp;
delaytime(1)=0.0;
phi=delaytime/dt/nf.*(0:nf-1);
for k=1:nf
 if (omg(k)>omg0(size_table))
phi(k)=0.0;
 end
end
phi(1)=0.0;
x_disperse=exp(-2*pi*1i*phi);
% set zero amplitude for omg larger than the table
for k=1:nf/2+1
 if (omg(k)>omg0(size_table))
        x_disperse(k)=0.0+1i*0.0;
 end
end
% setup wrap-around for FFT
for k=1:nf/2-1
x_disperse(nf-k+1)=conj(x_disperse(k+1));
end
x_disperse(nf/2+1)=0.0+1i*0.0;
% METHOD I
% setting of time shift operator without for loop
%
x_const=exp(-2*pi*1i*w*[0:nf/2-1 -nf/2:-1]); % tsunami without dispersion
x_const(nf/2+1)=0.0+1i*0.0;
% end of METHOD I
% METHOD II
% setting of time shift operator construction by for-loop
z=zeros(1,nf);
for k=1:nf/2
z(k)=exp(-2*pi*(k-1)*1i*w);
end
for k=1:nf/2-1
z(nf-k+1)=conj(z(k+1));
end
z(nf/2+1)=exp(2*pi*nf/2*1i*w);
z(nf/2+1)=0.0+1i*0.0;
% end of setting of time shift operator construction by for-loop
x=z;
% END of METHOD II
x=x_disperse;% this sets tsunami dispersion

S=S0.*x;
y1= real(ifft(S,nf));
SS1=S.*conj(S)/nf;
SS1=sqrt(SS1);
S=S.*x;
y2= real(ifft(S,nf));
SS2=S.*conj(S)/nf;
SS2=sqrt(SS2);
S=S.*x;
y3= real(ifft(S,nf));
SS3=S.*conj(S)/nf;
SS3=sqrt(SS3);
figure
plot(tt,y,tt,y1,tt,y2,tt,y3);
set(gca,'fontsize',16,'linewidth',2)
xlabel('Time (sec)','fontsize',16)
ylabel('Linear amplitude','fontsize',16)
title('Tsunami at 0, 3000, 6000, 9000 km with PREM dispersion','fontsize',16)
dlmwrite('tsunami_timeaxis.txt',tt.');
dlmwrite('tsunami_4km_yn_0km.txt',real(y).');
dlmwrite('tsunami_4km_yn_3000km.txt',real(y1).');
dlmwrite('tsunami_4km_yn_6000km.txt',real(y2).');
dlmwrite('tsunami_4km_yn_9000km.txt',real(y3).');

x=x_const; % this sets tsunami without dispersion
S=S0.*x;
y1= real(ifft(S,nf));
SS1=S.*conj(S)/nf;
SS1=sqrt(SS1);
S=S.*x;
y2= real(ifft(S,nf));
SS2=S.*conj(S)/nf;
SS2=sqrt(SS2);
S=S.*x;
y3= real(ifft(S,nf));
SS3=S.*conj(S)/nf;
SS3=sqrt(SS3);
figure
plot(tt,y,tt,y1,tt,y2,tt,y3);
set(gca,'fontsize',16,'linewidth',2)
xlabel('Time (sec)','fontsize',16)
ylabel('Linear amplitude','fontsize',16)
title('Tsunami at 0, 3000, 6000, 9000 km without dispersion','fontsize',16)

dlmwrite('tsunami_4km_constv_3000km.txt',real(y1).');
dlmwrite('tsunami_4km_constv_6000km.txt',real(y2).');
dlmwrite('tsunami_4km_constv_9000km.txt',real(y3).');

% do not change the spectrum (except by providing finer sampling, df)
% because the zeros don't contribute to the sum of the area in the 
% Fourier integral
fnyq = 1/(2*dt); % Nyquist frequency
df = fnyq/(nf/2); % calculate frequency sample interval
f = 0:df:fnyq; % calculate frequency variable (will be (nf/2) + 1 long
SS = S.*conj(S)/nf;  
SS = sqrt(SS); % calculate amplitude spectrum
figure
plot(f,SS(1:(nf/2)+1),f,SS1(1:nf/2+1),f,SS2(1:nf/2+1),f,SS3(1:nf/2+1),'-r','linewidth',1.5)
set(gca,'fontsize',16,'linewidth',2)
xlabel('Frequency (Hz)','fontsize',16)
ylabel('Linear amplitude','fontsize',16)
header_line=sprintf('Amplitude Spectrum of %s Wavelet',wave_type);
title(header_line,'fontsize',16)
