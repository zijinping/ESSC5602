% problem 7.05
% modified from eda07_09


% create two synthetic datasets
% h1, a unit spike
% h2, random noise with zero mean and unit variance
N=101;
Dt=1;
t=Dt*[0:N-1]';
h1=zeros(N,1);
h1(floor(N/2))=1;
h2=random('Normal',0,1,N,1);

% smoothing length
tau=4*Dt;
c=exp(-Dt/tau);

% create component filters
u=[1-c,0]';
v=[1.0,-c];

% filter time series
q1 = filter(u, v, h1);
q2 = filter(u, v, h2);

% plot
figure(1);
clf;
subplot(2,1,1);
hold on;
set(gca,'LineWidth',2);
axis( [t(1), t(N), -1.1, 1.1] );
plot( t, h1, 'k.', 'LineWidth', 2 );
plot( t, q1, 'k-', 'LineWidth', 2 );
xlabel('time, t');
ylabel('h1(t) and q1(t)');
subplot(2,1,2);
hold on;
set(gca,'LineWidth',2);
axis( [t(1), t(N), -5, 5] );
plot( t, h2, 'k.', 'LineWidth', 2 );
plot( t, q2, 'k-', 'LineWidth', 2 );
xlabel('time, t');
ylabel('h2(t) and q2(t)');

disp(sprintf('sum of elements of the unit spike: %f', sum(h1)));
disp(sprintf('sum of elements of the filtered spike: %f', sum(q1)));