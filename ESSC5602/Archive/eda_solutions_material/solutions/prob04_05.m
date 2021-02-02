% problem 4.05
% modifield from eda04_11
% load Black Rock Forest temperature data

% read the data
Draw=load('brf_temp.txt');
traw=Draw(:,1);
draw=Draw(:,2);
Nraw=length(draw);

% exclude bad data
n = find( (draw~=0) & (draw>-40) & (draw<38) );
t=traw(n);
d=draw(n);

N=length(d);
M=14;

% plot the data
figure(1);
clf;
subplot(3,1,1);
set(gca,'LineWidth',2);
hold on;
axis( [0, 5000, -40, 40] );
plot(t,d,'k-','LineWidth',2);
xlabel('time, days');
ylabel('obs temp, C');

% set up data kernel
Ty=365.25;
Td=1.0;
G=zeros(N,M);
G(:,1)=1;
G(:,2)=t;
G(:,3)=cos(2*pi*t/Ty);
G(:,4)=sin(2*pi*t/Ty);
G(:,5)=cos(2*pi*t/(Ty/2));
G(:,6)=sin(2*pi*t/(Ty/2));
G(:,7)=cos(2*pi*t/(Ty/3));
G(:,8)=sin(2*pi*t/(Ty/3));
G(:,9)=cos(2*pi*t/Td);
G(:,10)=sin(2*pi*t/Td);
G(:,11)=cos(2*pi*t/(Td/2));
G(:,12)=sin(2*pi*t/(Td/2));
G(:,13)=cos(2*pi*t/(Td/3));
G(:,14)=sin(2*pi*t/(Td/3));

% predict data
mest = (G'*G)\(G'*d);
dpre = G*mest;
e = d - dpre;
E = e'*e;
subplot(3,1,2);
set(gca,'LineWidth',2);
hold on;
axis( [0, 5000, -40, 40] );
plot(t,dpre,'k-','LineWidth',2);;
xlabel('time, days');
ylabel('pre temp, C')

% plot error
subplot(3,1,3);
set(gca,'LineWidth',2);
hold on;
axis( [0, 5000, -40, 40] );
plot(t,e,'k-','LineWidth',2);
xlabel('time, days');
ylabel('error, C')

sigmadA = 0.01; % prior estimate, 0.01 deg C
sigmadB = sqrt(E/(N-M)); % posterior estimate, from fit
Cm = inv(G'*G);
slope = Ty*mest(2);
errslopeA = Ty*sigmadA*sqrt(Cm(2,2)); % note conversion to years
errslopeB = Ty*sigmadB*sqrt(Cm(2,2));

slope
errslopeA
errslopeB
