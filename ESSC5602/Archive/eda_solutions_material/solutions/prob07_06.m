% problem 7.06
% modified from prob07_09

% read in Neuse River hydrograph
D = load('neuse.txt');
t = D(:,1);
N=length(t);
Dt=t(2)-t(1);
t=Dt*[0:N-1]';
h = D(:,2);

Ntau = 3;
taus = [5, 15, 40]';

% plot
tminplot=0;
tmaxplot=1000;
figure(1);
clf;
subplot(Ntau+1,1,1);
hold on;
set(gca,'LineWidth',2);
axis( [tminplot, tmaxplot, min(h), max(h)] );
plot( t, h, 'k-', 'LineWidth', 2 );
xlabel('time, t');
ylabel('discharge, cfs');

for itau = [1:Ntau]

% smoothing length
tau=taus(itau);
c=exp(-Dt/tau);

% create component filters
u=[1-c,0]';
v=[1.0,-c];

% filter time series
q = filter(u, v, h);

% plot
subplot(Ntau+1,1,itau+1);
hold on;
set(gca,'LineWidth',2);
axis( [tminplot, tmaxplot, min(h), max(h)] );
plot( t, q, 'k-', 'LineWidth', 2 );
xlabel('time, t');
ylabel(sprintf('tau = %.1f',tau));

end
