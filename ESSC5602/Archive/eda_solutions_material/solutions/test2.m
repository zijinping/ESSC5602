
clear F;
global F;

% load data but use just a section that brackets a hole
D = load('brf_temp.txt');
xraw = D(28000:33000,1);
draw = D(28000:33000,2);

% calculate sampling interval, min, max of data
Dx=xraw(2)-xraw(1);
Dx2 = Dx^2;
xmin=min(xraw);
xmax=max(xraw);
dmin=min(draw);
dmax=max(draw);

% The dataset set contains both missing and zero data.
% Step one is to fill in the missing data to zero. The
% array, m, is the index of the the available data into
% a timeseries that is equally-spaced with the sampling, Dx
m = floor((xraw-xmin)/Dx)+1;
M = max(m);
xt=xmin+Dx*[0:M-1]';
dt=zeros(M,1);
dt(m)=draw;

% now locate zeros, cold spikes and hot spikes
rowindex=find( (dt~=0) & (dt>-40) & (dt<38) );
d=dt(rowindex);
x=xt(rowindex);
N=length(d);

% plot the data
figure(1);
clf;
subplot(2,1,1);
set(gca,'LineWidth',2);
hold on;
axis([xmin xmax dmin dmax]);
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
axis([xmin xmax dmin dmax]);
plot(xt,mest,'k-','LineWidth',1);
xlabel('x');
ylabel('d est');

e = d-mest(rowindex(:));
E=e'*e;
disp(sprintf('the error of the data %f',E));


