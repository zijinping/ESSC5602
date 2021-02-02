% prob09_03
% modified from eda09_03
% load the reynolds_interpolatede.txt data set
% compute and plot autocorrelation

D=load('reynolds_interpolated.txt');
t=D(:,1); % time
d=D(:,7)-mean(D(:,7)); % chlorophyll. demeaned
N=length(d);
Dt=1;

a = xcorr(d-mean(d));
tau = Dt*[-N+1:N-1]';
Nk = [1:N, N-1:-1:1]';

% correction for length of timeseries
a = a ./ Nk;

figure(1);
clf;

subplot(2,1,1);
set(gca,'LineWidth',2);
hold on;

axis( [ -30*3, 30*3, 0, 1.1*max(a) ] );
plot( tau, a, 'k-', 'LineWidth', 2 );
xlabel('lag, days');
ylabel('autocorrelation');

subplot(2,1,2);
set(gca,'LineWidth',2);
hold on;

axis( [ -365.25*6, 365.25*6, -1.1*max(a), 1.1*max(a) ] );
plot( tau, a, 'k-', 'LineWidth', 2 );
xlabel('lag, days');
ylabel('autocorrelation');

