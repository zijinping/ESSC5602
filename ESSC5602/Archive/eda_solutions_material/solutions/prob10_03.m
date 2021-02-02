% problem 10.03
% modifield from eda10_04
% simple example of kriging

% the model parameters are the values of a sinusoidal
% function evenly spaced along the x axis
M=101;
Dx=1.0;
xmin=0;
x = xmin+Dx*[0:M-1]';
xmin=x(1);
xmax=x(M);

% the data are the vales the function measured on
% only a subset of points along the x axis
N=50;
% randomly choose N of the M points to have data
rowindex=sort(unidrnd(M,N,1));
% determine their x positions
xobs=xmin+Dx*rowindex;
% evaluate a sinusoidal function at these points to
% simulate the data
freq=0.02;
dobs = sin(2*pi*freq*xobs);

% estimate data at these points
xest=x;
dest=zeros(M,1);

% plot the observed data with circles
figure(1);
clf;

% loop over a few choices of L2

L2v = [1, 16, 900, 10000]';
NL2 = length(L2v);
for iplot = [1:NL2]

L2=L2v(iplot);

subplot(1,NL2,iplot);
set(gca,'LineWidth',2);
hold on;
axis([xmin xmax -2 2]);
plot(xobs,dobs,'ko','LineWidth',2);
xlabel('x');
ylabel('d');
title(sprintf('L=%.1f',sqrt(L2)));

% assume the autocorrelation function of the data
% is a a Normal function with variance L^2.  The normalization
% constant cancels, so just set to unity

% basics equations are Aw=v(i), dest(i)=dobs'*w(i);

% set up the kriging matrix, A
A = exp(-abs(xobs*ones(N,1)'-ones(N,1)*xobs').^2 /(2*L2));

% set up kriging vector, v, but with vector for each interplation
% point a column of the matrix, V
V = exp(-abs(xobs*ones(M,1)'-ones(N,1)* xest').^2 /(2*L2));

% solve, but add some damping, just in case A
% is near singular.  This selects the smallest w
% that solves the equation
dest=dobs'*((A+1e-6*eye(N))\V);
    
% plot estimated data with line
plot(xest,dest,'k-','LineWidth',2);

end

