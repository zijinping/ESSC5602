% problem 7.2
% modified from eda07_05

clear g H;
global g H;

% read in data
D = load('neuse.txt');

% set up timeseries
t=D(:,1);
d=D(:,2);
N = length(d);
Dt = t(2)-t(1);
tmin = 0.0;
t = tmin + Dt*[0:N-1]';
tmax = t(N);

% compute pef's of these lengths
iE = [2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]';
NPEL = length(iE);
E=zeros(NPEL,1);

for ipel = [1:NPEL]
    
% sizes of various arrays
M=iE(ipel); % length of pef 
K=1;
L=N+K;

% for a prediction error filter, pef, the data is the filter, g
g = d;

% prior information that first element of the pef is (-1)
e=(10^1)*max(abs(g)); % large weight
H=zeros(K,M);
H(1,1)=e;
h=zeros(1,1);
h(1)=(-e);

% FT f needed only for biconjugate gradient method
% FTf = GT 0 + HT h = [e 0 0 ... 0]T [-e]
FTf = zeros(M,1);
FTf(1) = -e^2;

% solve for predictione error filter using generalized least squares
% constraint that pef(1)=(-1) implemented with generalized least squares
pef2=bicg(@filterfun,FTf,1e-12,3*L);

% prediction error
perror = conv(pef2,g);
E(ipel) = perror'*perror;

end

% plot errors
figure(1);
clf;
set(gca,'LineWidth',2);
plot( iE, E, 'k-', 'LineWidth', 2);
xlabel('length of filter');
ylabel('Error');








