% Problem 4.4

% read data
D = load('LineData01.txt');
x=D(:,1);
d=D(:,2);
N=length(d);

% degree 2 fit
M=3;
G=zeros(N,M);
G(:,1)=1;
G(:,2)=x;
G(:,3)=x.^2;
mest=(G'*G)\(G'*d);

% prediction
dpre=G*mest;

% variance calculation 
sigma2d = (d-dpre)'*(d-dpre)/(N-M);
Cm = sigma2d * inv(G'*G);

% output info
disp('2nd degree fit');
disp(sprintf('variance of data: %f', sigma2d ));
disp(sprintf('m1: %f +/- %f (2 sigma)', mest(1), 2*sqrt(Cm(1,1)) ));
disp(sprintf('m2: %f +/- %f (2 sigma)', mest(2), 2*sqrt(Cm(2,2)) ));
disp(sprintf('m3: %f +/- %f (2 sigma)', mest(3), 2*sqrt(Cm(3,3)) ));
disp(' ');

% plot fit
figure(1);
clf;
set(gca,'LineWidth',2);
hold on;
plot( x, d, 'ko', 'LineWidth', 2 );
plot( x, dpre, 'k-', 'LineWidth', 2 );
xlabel('x');
ylabel('d');
title('2nd degree polynomial fit');

% degree 3 fit
M=4;
G=zeros(N,M);
G(:,1)=1;
G(:,2)=x;
G(:,3)=x.^2;
G(:,4)=x.^3;
mest=(G'*G)\(G'*d);

% prediction
dpre=G*mest;

% variance calculation 
sigma2d = (d-dpre)'*(d-dpre)/(N-M);
Cm = sigma2d * inv(G'*G);

% output info
disp('3nd degree fit');
disp(sprintf('variance of data: %f', sigma2d ));
disp(sprintf('m1: %f +/- %f (2 sigma)', mest(1), 2*sqrt(Cm(1,1)) ));
disp(sprintf('m2: %f +/- %f (2 sigma)', mest(2), 2*sqrt(Cm(2,2)) ));
disp(sprintf('m3: %f +/- %f (2 sigma)', mest(3), 2*sqrt(Cm(3,3)) ));
disp(sprintf('m4: %f +/- %f (2 sigma)', mest(4), 2*sqrt(Cm(4,4)) ));
disp(' ');

% plot fit
figure(2);
clf;
set(gca,'LineWidth',2);
hold on;
plot( x, d, 'ko', 'LineWidth', 2 );
plot( x, dpre, 'k-', 'LineWidth', 2 );
xlabel('x');
ylabel('d');
title('3rd degree polynomial fit');

% degree 4 fit
M=5;
G=zeros(N,M);
G(:,1)=1;
G(:,2)=x;
G(:,3)=x.^2;
G(:,4)=x.^3;
G(:,5)=x.^4;
mest=(G'*G)\(G'*d);

% prediction
dpre=G*mest;

% variance calculation 
sigma2d = (d-dpre)'*(d-dpre)/(N-M);
Cm = sigma2d * inv(G'*G);

% output info
disp('4th degree fit');
disp(sprintf('variance of data: %f', sigma2d ));
disp(sprintf('m1: %f +/- %f (2 sigma)', mest(1), 2*sqrt(Cm(1,1)) ));
disp(sprintf('m2: %f +/- %f (2 sigma)', mest(2), 2*sqrt(Cm(2,2)) ));
disp(sprintf('m3: %f +/- %f (2 sigma)', mest(3), 2*sqrt(Cm(3,3)) ));
disp(sprintf('m4: %f +/- %f (2 sigma)', mest(4), 2*sqrt(Cm(4,4)) ));
disp(sprintf('m5: %f +/- %f (2 sigma)', mest(5), 2*sqrt(Cm(5,5)) ));
disp(' ');

% plot fit
figure(3);
clf;
set(gca,'LineWidth',2);
hold on;
plot( x, d, 'ko', 'LineWidth', 2 );
plot( x, dpre, 'k-', 'LineWidth', 2 );
xlabel('x');
ylabel('d');
title('4th degree polynomial fit');

