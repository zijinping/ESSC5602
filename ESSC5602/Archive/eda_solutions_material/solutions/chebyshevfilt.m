% chebyshevfilt
% chebyshev IIR bandpass filter
% din - input array of data
% Dt - sampling interval
% flow - low pass frequency, Hz
% fhigh - high pass frequency, Hz
% dout - output array of data
% u - the numerator filter
% v - the denominator filter
% these filters can be used again using dout=filter(u,v,din);

function [dout, u, v] = chebyshevfilt(din, Dt, flow, fhigh)

% sampling rate
rate=1/Dt;

% ripple parameter, set to ten percent
ripple=0.1;  

% normalise frequency
fl=2.0*flow/rate;
fh=2.0*fhigh/rate;

% center frequency 
cf = 4 * tan( (fl*pi/2) ) * tan( (fh*pi/2) );

% bandwidth
bw = 2 * ( tan( (fh*pi/2) ) - tan( (fl*pi/2) ) );

% ripple parameter factor
rpf = (sqrt((1.0+1.0/(ripple*ripple))) + 1.0/ripple) ^ (0.5);
a = 0.5*(rpf-1.0/rpf);
b = 0.5*(rpf+1.0/rpf);

u=zeros(5,1);
v=zeros(5,1);
theta = 3*pi/4;
sr = a * cos(theta);
si = b * sin(theta);
es = sqrt(sr*sr+si*si);
tmp= 16.0 - 16.0*bw*sr + 8.0*cf + 4.0*es*es*bw*bw - 4.0*bw*cf*sr + cf*cf;
v(1) = 1.0;
v(2) = 4.0*(-16.0 + 8.0*bw*sr - 2.0*bw*cf*sr + cf*cf)/tmp;
v(3) = (96.0 - 16.0*cf - 8.0*es*es*bw*bw + 6.0*cf*cf)/tmp;
v(4) = (-64.0 - 32.0*bw*sr + 8.0*bw*cf*sr + 4.0*cf*cf)/tmp;
v(5) = (16.0 + 16.0*bw*sr + 8.0*cf + 4.0*es*es*bw*bw + 4.0*bw*cf*sr + cf*cf)/tmp;
tmp = 4.0*es*es*bw*bw/tmp;
u(1) = tmp;
u(2) = 0.0;
u(3) = -2.0*tmp;
u(4) = 0.0;
u(5) = tmp;

[dout]=filter(u,v,din);

return

