
clear F;
global F;

% the model parameters are the values of a sinusoidal
% function evenly spaced along the x axis
M=1001;
Dx=1.0;
Dx2 = Dx^2;
xmin=0.0;
xmax = xmin+Dx*(M-1);
xt = xmin+Dx*[0:M-1]';
dt = 2+sin(5*pi*xt/(xmax-xmin));
dt(M/5:2*M/5)=0;

rowindex=find(dt~=0);
d=dt(rowindex);
x=xt(rowindex);
N=length(d);

% plot the data
figure(1);
clf;
subplot(2,1,1);
set(gca,'LineWidth',2);
hold on;
axis([xmin xmax -4 4]);
plot(xt,dt,'k-','LineWidth',1);
plot(x,d,'k.','LineWidth',1);
xlabel('x');
ylabel('d true');

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

subplot(2,1,2);
set(gca,'LineWidth',2);
hold on;
axis([xmin xmax -4 4]);
plot(xt,mest,'k-','LineWidth',1);
xlabel('x');
ylabel('d est');

e = d-mest(rowindex(:));
E=e'*e;

