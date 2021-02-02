% Problem 9.05
% modified from eda09_16

% load data
D = load('brf_temp.txt');
Dt = 1/24;
traw = D(:,1);
draw = D(:,2);
Nraw = length(traw);
diy = 365.25;

% The dataset set contains both missing and zero data.
% Step one is to fill in the missing data to zero. The
% array, m, is the index of the the available data into
% a timeseries that is equally-spaced with the sampling, Dx
D = floor(((traw(2:Nraw)-traw(1:Nraw-1))/Dt)+0.5);
m = [1, 1+cumsum(D)']';
Nc = max(m);
tc=Dt*[0:Nc-1]';
dc=zeros(Nc,1);
dc(m)=draw;

% now locate cold spikes and hot spikes and set to zero
rowindex=find( (dc<-40) | (dc>38) );
dc(rowindex)=0;

% bandpass filter at frequency around 1 cycle per day
flow = 0.9;
fhigh = 1.1;
dfiltered = chebyshevfilt(dc, Dt, flow, fhigh);

figure(1);
clf;
set(gca,'LineWidth',2);
hold on;
plot( tc/diy, dfiltered, 'k-', 'LineWidth', 2 );
xlabel('time, years');
ylabel('filtered temperature, K');

