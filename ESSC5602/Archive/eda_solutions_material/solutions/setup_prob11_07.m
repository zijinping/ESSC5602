% create synthetic data for eda_problem_11_07
clear all

N=10001;
Dt = 1;
t = Dt*[0:N-1]';

x = zeros(N,1);
y = zeros(N,1);
f = zeros(N,1);

Np = floor(N/10);
p = random('unid',N,Np,1);
cexp = 0.5;
x(p) = random('exp',cexp,Np,1);

figure(1);
clf;
subplot(3,1,1);
set(gca,'LineWidth',2);
hold on;
axis( [0, N, 0, 2 ] );
plot(t,x,'k-','LineWidth',2);

cfmin = 0.2;
cfmax = 0.5;
cf = cfmin + (cfmax-cfmin)*[0:N-1]/(N-1)';
for i=[2:N]
    f(i) = cf(i).*(y(i-1)^2);
    dydt = x(i) - f(i);
    y(i) = y(i-1) + dydt * Dt;
    if( y(i) < 0 )
        y(i)=0;
    end
end

subplot(3,1,2);
set(gca,'LineWidth',2);
hold on;
axis( [0, N, 0, 1.1*max(f) ] );
plot(t,f,'k-','LineWidth',2);

subplot(3,1,3);
set(gca,'LineWidth',2);
hold on;
axis( [0, N, 0, 1.1*max(y) ] );
plot(t,y,'k-','LineWidth',2);

dlmwrite('datafor1107.txt', [t,f,x], 'delimiter', '\t' );


