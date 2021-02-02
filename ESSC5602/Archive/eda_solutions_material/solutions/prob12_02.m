% problem 12.2
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

% test significance of every paris
PF = zeros(NPEL,NPEL);
for ipel = [1:NPEL]
for jpel = [ipel:NPEL] 
    NA = N-iE(ipel);
    NB = N-iE(jpel);
    F = (E(ipel)/NA) / (E(jpel)/NB);
    if( F<1 )
        F=1/F;
    end
    PF(ipel,jpel) = 1 - (fcdf(F,NA,NB)-fcdf(1/F,NA,NB));
    PF(jpel,ipel) = PF(ipel,jpel);
    disp(sprintf('E(%d)=%f E(%d)=%f Fest=%f, P(F<1/Fest|F>Fest) = %f', iE(ipel), E(ipel), iE(jpel), E(jpel), F, PF(ipel,jpel)));
end
end

eda_draw(' ', PF, 'caption P(F<Fest)', ' ', PF<0.05, 'caption different' );


