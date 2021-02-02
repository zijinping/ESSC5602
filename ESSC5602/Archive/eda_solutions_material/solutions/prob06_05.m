% Problem 6.05

% load data 
D = load('brf_temp.txt');
traw = D(:,1);
draw = D(:,2);
Nraw = length(traw);

% calculate sampling interval, min, max of data
Dt=1/24;
tmin=min(traw);
tmax=max(traw);
dmin=min(draw);
dmax=max(draw);

% The dataset set contains both missing and zero data.
% Step one is to fill in the missing data to zero. The
% array, m, is the index of the the available data into
% a timeseries that is equally-spaced with the sampling, Dx

D = floor(((traw(2:Nraw)-traw(1:Nraw-1))/Dt)+0.5);
m = [1, 1+cumsum(D)']';
N = max(m);
tc=tmin+Dt*[0:N-1]';
dc=zeros(N,1);
dc(m)=draw;

% now locate cold spikes and hot spikes and set to zero
rowindex=find( (dc<-40) | (dc>38) );
dc(rowindex)=0;

% plot the data
figure(1);
clf;
set(gca,'LineWidth',2);
hold on;
axis([tmin tmax dmin dmax]);
plot(tc,dc,'k-','LineWidth',2);
xlabel('x');
ylabel('d true');

% even number of data
N=2*floor(N/2);
T = N*Dt;
dc = dc(1:N);

% frequency definition
fmax=1/(2.0*Dt);
df=fmax/(N/2);
f=df*[0:N/2,-N/2+1:-1]'; 
Nf=N/2+1;
dw=2*pi*df;
w=dw*[0:N/2,-N/2+1:-1]';
Nw=Nf;

% fft and power spectral density, pos freqs only
dtilde = Dt*fft(dc);
s2 = (2/T) * abs(dtilde).^2;
s2 = s2(2:Nf);
f = f(2:Nf);

% plot the psd on log-log scale
figure(3);
clf;
set(gca,'LineWidth',2);
hold on;
axis( [-3, 1, 0, max(log10(s2)) ] );
plot(log10(f),log10(s2),'k-','LineWidth',2);
xlabel('log10 frequency, cycles per day');
ylabel('log10 psd, K-squared per cpd');



