% problem 6.4

% i used as imaginary unit!
clear i;

% set up time
N=128;
Dt=1.0;
T=N*Dt; 
t=Dt*[0:N-1]';
ta = Dt*[0:N/2,-N/2+1:-1]'; 

% set up frequency
fmax=1/(2.0*Dt);
df=fmax/(N/2);
f=df*[0:N/2,-N/2+1:-1]'; 
Nf=N/2+1;
dw=2*pi*df;
w=dw*[0:N/2,-N/2+1:-1]';
Nw=Nf;

% method 1
sigma=10.0;
d = exp( -ta.^2 / (2*sigma) );

% plot d(f)
figure(1);
clf;
subplot(2,1,1);
set(gca,'LineWidth',2);
hold on;
axis( [0, T, 0, 1] );
plot( t, d, 'k-', 'LineWidth', 2 );
xlabel('time, s');
ylabel('d(t)');

% fft
dtilde1 = Dt * fft(d);

% plot real part of transform d_tilde(f)
figure(2);
clf;
subplot(2,1,1);
set(gca,'LineWidth',2);
hold on;
axis( [-fmax, fmax, min(real(dtilde1)), max(real(dtilde1))] );
f1 = [ f(Nf+1:N)', f(1:Nf)']';
dtilde2 = [ dtilde1(Nf+1:N)', dtilde1(1:Nf)']';
plot( f1, real(dtilde2), 'k-', 'LineWidth', 2 );
plot( f1, imag(dtilde2), 'k:', 'LineWidth', 2 );
xlabel('frequency, Hz');
ylabel('real(dt(f)) (method 1)');

% method 2
sigma=10.0;
t0 = T/2;
d = exp( -(t-t0).^2 / (2*sigma) );

% plot d(t)
figure(1);
subplot(2,1,2);
set(gca,'LineWidth',2);
hold on;
axis( [0, T, 0, 1] );
plot( t, d, 'k-', 'LineWidth', 2 );
xlabel('time, s');
ylabel('d(t)');

% fft and phase shift
dtilde3 = Dt * fft(d) .* exp(i*w*t0);

figure(2);
subplot(2,1,2);
set(gca,'LineWidth',2);
hold on;
axis( [-fmax, fmax, min(real(dtilde2)), max(real(dtilde2))] );
f1 = [ f(Nf+1:N)', f(1:Nf)']';
dtilde4 = [ dtilde3(Nf+1:N)', dtilde3(1:Nf)']';
plot( f1, real(dtilde4), 'k-', 'LineWidth', 2 );
plot( f1, imag(dtilde4), 'k:', 'LineWidth', 2 );
xlabel('frequency, Hz');
ylabel('real(dt(f)) (method 2)');


