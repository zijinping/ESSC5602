% problem 9.04
% modified from eda06_03
% power spectral density of the Neuse River hydrograph

% load data
D = load('neuse.txt');
Dt = D(2,1)-D(1,1);
d = D(:,2);
N=length(d);
T = N*Dt;

% remove mean
d = d - mean(d);

% round off to even number of points
N=floor(N/2)*2;
M=N;
d=d(1:N);

% set up frequencies
Nf=N/2+1;
fmax = 1/(2*Dt);
Df = fmax/(N/2);
f = Df*[0:Nf-1]';

% compute power spectral density of raw data
dtilde = Dt * fft(d);
s2 = (2/T) * abs(dtilde).^2;

% find max
[s2max, is2max] = max(s2);
disp(sprintf('largest peak at frequency %f (raw)', (is2max-1)*Df));

% plot spectral density
figure(1)
clf;

% plot vs frequency
subplot(2,1,1);
set(gca,'LineWidth',2);
hold on;
axis( [ 0, 0.05, 0, max(s2(2:Nf)) ] );
plot(f(2:Nf),s2(2:Nf),'k-','LineWidth',2);
xlabel('frequency, cycles per day');
ylabel('PSD (raw)');

% Hamming filter
W=0.54-0.46*cos(2*pi*[0:N-1]'/(N-1));

% compute power spectral density of tapered data
dtilde = Dt * fft(W.*d);
s2 = (2/T) * abs(dtilde).^2;

% find max
[s2max, is2max] = max(s2);
disp(sprintf('largest peak at frequency %f (tapered)', (is2max-1)*Df));

% plot vs frequency
subplot(2,1,2);
set(gca,'LineWidth',2);
hold on;
axis( [ 0, 0.05, 0, max(s2(2:Nf)) ] );
plot(f(2:Nf),s2(2:Nf),'k-','LineWidth',2);
xlabel('frequency, cycles per day');
ylabel('PSD (tapered)');
