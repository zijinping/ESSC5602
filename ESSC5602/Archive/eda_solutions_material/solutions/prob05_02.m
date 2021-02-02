% problem 5.02
% modified from eda05_05
% simple example of smoothness applied
% to the problem of filling in data gaps
% sparse matrix implementation

clear F;
global F;

% the model parameters are the values of a sinusoidal
% function evenly spaced along the x axis
M=101;
Dx=1.0;
Dx2 = Dx^2;
xmin=0.0;
xmax = xmin+Dx*(M-1);
x = xmin+Dx*[0:M-1]';

% the data are the vales the function measured on
% only a subset of points along the x axis
N=40;
% randomly choose N of the M points to have data
rowindex=sort(unidrnd(M,N,1));
% determine their x positions
xd=xmin+Dx*rowindex;
% evaluate a sinusoidal function at these points to
% simulate the data
freq=0.02;
d = sin(2*pi*freq*xd);

% add repeated data.  Pick the middke x, x(k) with k=N/2,
% and add ten versions of it to the end of the dataset,
% each with different noise

k=N/2;
for i = [1:10]
    N=N+1;
    xd(N) = xd(k);
    d(N) = d(k)+random('Normal',0,0.5);
    rowindex(N) = rowindex(k);
end

% plot the data
figure(1);
clf;
set(gca,'LineWidth',2);
hold on;
axis([xmin xmax -2 2]);
plot(xd,d,'ko','LineWidth',2);
xlabel('x');
ylabel('d');

% the first N rows of F are that the model parameter
% equals the observed data
L=N+M;
F=spalloc(L,M,3*L);
f=zeros(L,1);
for p = [1:N]
    F(p,rowindex(p)) = 1;
    f(p)=d(p);
end

% shi is 1/sh, where sh^2 is the variance of the smoothness
shi = 1e-1;

% implement the smoothness constraints in the interior
% of the x interval
for p = [1:M-2]
    F(p+N,p) = shi/Dx2;
    F(p+N,p+1) = -2*shi/Dx2;
    F(p+N,p+2) = shi/Dx2;
    f(p+N)=0.0;
end

% implement flatness constraint at left end of the
% x interval
F(L-1,1)=-shi*Dx;
F(L-1,2)=shi*Dx;
f(L-1)=0;

% implement flatness constraint at right end of the
% x interbal
F(L,M-1)=-shi*Dx;
F(L,M)=shi*Dx;
f(L)=0;

% the old old way: mest= (F'*F)\(F'*f);
% the new old way: mest=bicg(F'*F,F'*f,1e-10,3*L);
% the best way:
mest=bicg(@afun,F'*f,1e-10,3*L);
plot(x,mest,'k-','LineWidth',2);

e = d-mest(rowindex(:));
E=e'*e;

