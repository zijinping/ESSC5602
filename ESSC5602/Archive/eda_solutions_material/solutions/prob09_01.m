% prob09.01
% modified from eda09_05
% using the cross correlation to align two timeseries
% West Point temperature data

D = load('ozone_nohead.txt');
t = D(:,1);         % time in days after 1:00 AM on Aug 1, 1993
ozone = D(:,2);     % ozone in ppm
radiation = D(:,3); % solar radiation in W/m2
temp = D(:,4);      % temperature in deg C

% time series paramaters
N=length(t);
Dt=1/24;

% plot original timeseries
figure(1);
clf;

% plot solar radiation
subplot(4,1,1);
set(gca,'LineWidth',2);
hold on;
axis( [t(1), t(N), 0, max(radiation)] );
plot( t, radiation, 'k-', 'LineWidth', 2);
xlabel('time, days');
ylabel('solar, W/m2');

% plot ozone
subplot(4,1,2);
set(gca,'LineWidth',2);
hold on;
axis( [t(1), t(N), min(temp), max(temp)] );
plot( t, temp, 'k-', 'LineWidth', 2);
xlabel('time, days');
ylabel('temp, degC');

% cross correlation
c = xcorr(radiation, temp);
[cmax, icmax] = max(c);

% compute lag
% this gives the lag of the second timeseries with respect to the frst
tlag = -Dt * (icmax-N);

% plot solar radiation
subplot(4,1,3);
set(gca,'LineWidth',2);
hold on;
axis( [t(1), t(N), 0, max(radiation)] );
plot( t, radiation, 'k-', 'LineWidth', 2);
xlabel('time, days');
ylabel('solar, W/m2');

% plot lagged ozone
subplot(4,1,4);
set(gca,'LineWidth',2);
hold on;
axis( [t(1), t(N), min(temp), max(temp)] );
plot( t-tlag, temp, 'k-', 'LineWidth', 2);
xlabel('time, days');
ylabel('temp, degC');

figure(2);
clf;

days=5;
subplot(2,1,1);
set(gca,'LineWidth',2);
hold on;
axis( [t(1), t(days*24), 0, max(radiation)] );
plot( t, radiation, 'k-', 'LineWidth', 2);
xlabel('time, days');
ylabel('solar radiation, W/m2');

subplot(2,1,2);
set(gca,'LineWidth',2);
hold on;
axis( [t(1), t(days*24), min(temp), max(temp)] );
plot( t, temp, 'k-', 'LineWidth', 2);
plot( t-tlag, temp, 'k:', 'LineWidth', 2);
title(sprintf('%.2f hour lag', 24*tlag));
xlabel('time, days');
ylabel('temp, degC');

disp(sprintf('estimated lag: %.2f hours', 24*tlag ));

figure(3);
clf;
set(gca,'LineWidth',2);
hold on;

Nc = length(c);
Ncenter = (Nc+1)/2;
hours=12;
axis( [-hours, hours, 0, 1.1*max(c)] );
plot( [-hours:hours]', c(Ncenter-hours:Ncenter+hours), 'k-', 'LineWidth', 2);
xlabel('time, hours');
ylabel('cross-correlation');
