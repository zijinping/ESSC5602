% Problem 11.03
% fit of two sinusoids of unknown (but approximately diurnal and annual)
% frequency to Black Rock Forest temperature data using the
% iterative least squares method
clear all;

% load Black Rock Forest temperature data
D=load('brf_temp.txt');
t=D(:,1);
dobs=D(:,2);
Ns=size(D);

% clean the data by deleting outliers and zero (= no data) values
k = find( (dobs>-35) & (dobs<40) & (dobs~=0.0) );
t = D(k,1);
dobs = D(k,2);
N = length(dobs);

% some statistics
dmean = mean(dobs);
amp = 0.5*(max(dobs)-min(dobs));

% Diurnal period, frequency, angular frequency
Td = 1;
fd = 1/Td;
wd = 2*pi*fd;

% Annual period, frequency, angular frequency
Ta = 365;
fa = 1/Ta;
wa = 2*pi*fa;

% nonlinear theory and its derivative
% d = m1*dmean + amp*m2*sin(m4 t) + amp*m3*cos(m4 t)
% dd/dm1 =  dmean
% dd/dm2 =  amp*sin(m4 t)
% dd/dm3 =  amp*cos(m4 t)
% dd/dm4 =  amp*m2*t*cos(m4 t) - amp*m3*t*sin(m4 t)

% inital guesses
M=7;
m = zeros(M,1);
m(1) =  1;          % average level is mean
% diurnal frequency
                    % amplitude quarter range of summer/winter
m(2) =  0.01*0.25;   % phase with max in midday
m(3) = -0.99*0.25;   % so only -cos() term
m(4) =  1;          % diurnal frequency 
% annual frequency
                    % amplitude full range of summer/winter
m(5) =  0.01*1;     % phase with max in summer
m(6) = -0.99*1;     % so only -cos() term
m(7) =  1;          % annual frequency 

s1 = amp*sin(m(4)*wd*t);
c1 = amp*cos(m(4)*wd*t);
s2 = amp*sin(m(7)*wa*t);
c2 = amp*cos(m(7)*wa*t);
dpre = dmean*m(1) + m(2)*s1 + m(3)*c1 + m(5)*s2 + m(6)*c2;

figure(1);
clf;
set(gca,'LineWidth',2);
hold on;
axis([0,5000,-30,50]);
plot(t,dobs,'k.','LineWidth',2);
hold on;
plot(t,dpre,'-','LineWidth',2,'Color',[0.7,0.7,0.7]);
xlabel('time (days)');
ylabel('temperature (deg C)');

figure(2);
clf;
set(gca,'LineWidth',2);
hold on;
axis([180,210,-30,50]);
plot(t,dobs,'k.','LineWidth',2);
hold on;
xlabel('time (days)');
ylabel('temperature (deg C)');

s1 = amp*sin(m(4)*wd*t);
c1 = amp*cos(m(4)*wd*t);
s2 = amp*sin(m(7)*wa*t);
c2 = amp*cos(m(7)*wa*t);
dpre = dmean*m(1) + m(2)*s1 + m(3)*c1 + m(5)*s2 + m(6)*c2;
dd = dobs - dpre;
E = dd'*dd;
for i=[1:10]
    Eprev = E;
    G = [ dmean*ones(N,1), s1, c1, m(2)*wd*t.*c1-m(3)*wd*t.*s1, ...
                           s2, c2, m(5)*wa*t.*c2-m(6)*wa*t.*s2];
    dm = (G'*G)\(G'*dd);
    m = m + dm;
    s1 = amp*sin(m(4)*wd*t);
    c1 = amp*cos(m(4)*wd*t);
    s2 = amp*sin(m(7)*wa*t);
    c2 = amp*cos(m(7)*wa*t);
    dpre = dmean*m(1) + m(2)*s1 + m(3)*c1 + m(5)*s2 + m(6)*c2;
    dd = dobs - dpre;
    E = dd'*dd;
    if( (abs(E-Eprev)/E) < 1.0e-6 )
        break
    end
end

figure(1);
plot(t,dpre,'-','LineWidth',2,'Color',[0.9,0.9,0.9]);
figure(2);
plot(t,dpre,'-','LineWidth',2,'Color',[0.9,0.9,0.9]);

fdfinal = (wd*m(4))/(2*pi);
Tdfinal = (2*pi)/(wd*m(4));
fafinal = (wa*m(7))/(2*pi);
Tafinal = (2*pi)/(wa*m(7));

% covariance
% f = f0 + df   and   T = T0 + dT    and   T0 = 1 / f0
% T = 1 / (f0 + df ) = (1/f0) (1+df/f0)^-1
%   = (1/f0) (1-df/f0) = (1/f0) - df/(f0^2)
% dT = (-1/(f0^2)) df
% var(dT) = (-1/(f0^2))^2 var(df)

G = [ dmean*ones(N,1), s1, c1, m(2)*wd*t.*c1-m(3)*wd*t.*s1, ...
                       s2, c2, m(5)*wa*t.*c2-m(6)*wa*t.*s2];
sigma2d = E/(N-M);
GTGI = inv(G'*G);
% diurnal frequency
sigma2m4 = sigma2d*GTGI(4,4);
sigma2fd = (wd^2)*sigma2m4/((2*pi)^2);
sigma2Td = sigma2fd / (fdfinal^4);
sigmaTd = sqrt(sigma2Td);
DT95d = 2*sigmaTd;

% annual frequency
sigma2m7 = sigma2d*GTGI(7,7);
sigma2fa = (wa^2)*sigma2m7/((2*pi)^2);
sigma2Ta = sigma2fa / (fafinal^4);
sigmaTa = sqrt(sigma2Ta);
DT95a = 2*sigmaTa;

fprintf('posterior sqrt(sigma2d) %f\n', sqrt(sigma2d));
fprintf('Estimated diurnal period %.3f +/- %.3f days\n', Tdfinal, DT95d );
fprintf('Estimated annual  period %.3f +/- %.3f days\n', Tafinal, DT95a );

