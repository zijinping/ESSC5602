% problem 12.05
% modifield from eda06_14 and eda12_11
% power spectral density (PSD) of the Neuse River hydrograph
% with units

D = load('neuse.txt');
Dt = D(2,1)-D(1,1);
d = D(:,2);
N=length(d);

% round off to even number of points
N=floor(N/2)*2;
d=d(1:N);

% generic time/frequency set up
t=Dt*[0:N-1]';
tmax = Dt*(N-1);
fmax=1/(2.0*Dt);
Df=fmax/(N/2);
f=Df*[0:N/2,-N/2+1:-1]'; 
Nf=N/2+1;

% compute power spectral sensity
dbar=Dt*fft(d);
s2 = (2/(N*Dt)) * abs(dbar(1:Nf)).^2;

p=2; % degrees of freedom
ff=1; % no window function, so set to unity
sd2est=std(d-mean(d))^2; % sd2true unknown, use sd2est instead, demean first
c = ( ff*sd2est ) / (2*Nf*Df);  % normalization factor

% plot data
figure(1)
clf;
subplot(2,1,1);
set(gca,'LineWidth',2);
hold on;
axis([0,tmax,0,max(d)]);
plot(t,d,'k-','LineWidth',2);
xlabel('time, days');
ylabel('discharge, cfs');

% significance levels
cv95000 = chi2inv(0.95,p);
cv99999 = chi2inv(0.99999,p);
% plot power spectral density (skip zero fequency)
% note that plot is only of low frequencies
subplot(2,1,2);
set(gca,'LineWidth',2);
hold on;
axis([0,0.1,log10(max(s2(2:Nf)))-5,log10(max(s2(2:Nf)))]);
plot(f(2:Nf),log10(s2(2:Nf)),'k-','LineWidth',2);
plot([f(1), f(Nf)], [log10(cv95000*c), log10(cv95000*c)],'k-','LineWidth',2);
plot([f(1), f(Nf)], [log10(cv99999*c), log10(cv99999*c)],'k-','LineWidth',2);
xlabel('frequency, cycles per day');
ylabel('log10 PSD, (cfs)^2 per cycle/day');

% power in the data and the PSD (which should be the same)
% note that the mean of the data zero frequency are being excluded
v1 = std(d)^2;
v2 = Df*sum(s2(2:Nf));
disp(sprintf('power in data: %f in PSD %f',v1,v2));

