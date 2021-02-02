% problem 6.1
% mofified from eda06_10
% derivative using fft example

D = load('neuse.txt');
Dt = D(2,1)-D(1,1);
N = length(D);
N = 2*floor(N/2); % make N even
d = D(1:N,2);

% standard set up
tmax=Dt*(N-1); 
t=Dt*[0:N-1]'; 
fmax=1/(2.0*Dt);
df=fmax/(N/2);
f=df*[0:N/2,-N/2+1:-1]'; 
Nf=N/2+1;
dw=2*pi*df;
w=dw*[0:N/2,-N/2+1:-1]';
Nw=Nf;

% derivative using finite differences
dddt1 = Dt*(d(2:N)-d(1:N-1));

% derivative using fft
dddt2=ifft(i*w.*fft(d));

% plot
figure(1);
clf;
subplot(3,1,1);
set(gca,'LineWidth',2);
hold on;
plot(t,d,'k-','LineWidth',2);
xlabel('time, t');
ylabel('d(t)');
subplot(3,1,2);
set(gca,'LineWidth',2);
hold on;
plot(t(1:N-1),dddt1,'k-','LineWidth',2);
xlabel('time, t');
ylabel('dd/dt(t) by fd');
subplot(3,1,3);
set(gca,'LineWidth',2);
hold on;
plot(t,real(dddt2),'k-','LineWidth',2);
xlabel('time, t');
ylabel('dd/dt(t) by fft');








