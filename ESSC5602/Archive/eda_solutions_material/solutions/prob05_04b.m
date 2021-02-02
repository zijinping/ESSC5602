% problem 5.04b
% modifield from eda05_09
% make simulated pressure data

KAPPA =  [0.05, 0.10, 0.20]';
SIGMAD = [0.05, 0.10, 1.00]';

for ik = [1:length(KAPPA)]
for is = [1:length(SIGMAD)]
    
% rows of grid are x
I = 41;
Dx = 1;
xmin = 0;
xmax = 40;
x = xmin + (xmax-xmin)*[0:I-1]'/(I-1);

% columns of grid are y
J = 41;
Dy = 1;
ymin = 0;
ymax = 40;
y = ymin + (ymax-ymin)*[0:J-1]'/(J-1);

% grid 10% populated with data
N = floor(I*J/10);
rowindex=unidrnd(I,N,1);
xobs = x(rowindex);
colindex=unidrnd(J,N,1);
yobs = y(colindex);

kappa = KAPPA(ik);
dtrue = 10*sin(kappa*xobs).*exp(-kappa*yobs);
sigmad = SIGMAD(is);
dobs = dtrue + random('normal',0.0,sigmad,N,1);

D=zeros(N,3);
D(:,1)=xobs;
D(:,2)=yobs;
D(:,3)=dobs;

dlmwrite( sprintf('pressure%d_%d.txt', ik, is),D,'\t');

end
end

