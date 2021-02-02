% Problem 2.2b
% alternate version using more sophisticated
% MatLab intrinsic functions

% read data from file
D=load('brf_temp.txt');
t=D(:,1);
d=D(:,2);
N = length(t);
Dt = t(2)-t(1);

% derivative
Dti = t(2:N) - t(1:N-1);
dip1 = d(2:N);
di = d(1:N-1);
dddt = (dip1-di)./Dti;

% find good data
goodi = find( ((abs(Dti-Dt)/Dt)<0.1) & (dip1>=-40) & (dip1<=38) & (dip1~=0) & (di>=-40) & (di<=38) & (di~=0) );
goodt = t(goodi);
gooddddt = dddt(goodi);

% find largest
[largest, ilargest] = max( abs(gooddddt) );
timeoflargest = goodt(ilargest);

disp(sprintf('largest hourly change is %f K/hr occurs at time %f days',largest/24,timeoflargest));





