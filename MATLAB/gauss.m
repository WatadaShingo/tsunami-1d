function [s] = gauss(npts,freq,dt,timesh,nsw)
%  Calculates gaussian wavelet
%  s(npts) = wavelet
%  dt = sampling interval
%  freq = predominate frequency of wavelet
%  timesh = time shift to the left (s)
%  nsw = number of points on each side of wavelet to be tapered 
%        by half of a cosine bell window
if timesh<.0001; timesh=npts*dt/2; end
s=zeros(1,npts);
const=2*pi*freq;
for i = 1:npts
   tim1 = (i-1)*dt;
   tim2 = tim1 - timesh;
   u = const*tim2;
   amp = exp(-u*u);
   s(i) = amp;
end

smax = max(abs(s));
smax=smax*2;
s = s/smax;
if nsw~=0
   for i = 1:nsw
      j = nsw-i;
      fac = 0.5*(cos(pi*j/nsw) + 1);
      s(i) = s(i)*fac;
   end
   for i = 1:nsw
      j = i-1;
      fac = 0.5*(cos(pi*j/nsw) + 1);
      k = npts-nsw+i;
      s(k) = s(k)*fac;
   end
end
