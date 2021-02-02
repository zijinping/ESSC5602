% problem 12.01
% modified from eda04_11
% Black Rock Forest temperature data

% read the data
Draw=load('brf_temp.txt');
traw=Draw(:,1);
draw=Draw(:,2);
Nraw=length(draw);

% exclude bad data
n = find( (draw~=0) & (draw>-40) & (draw<38) );
t=traw(n);
d=draw(n);

diy=365.15;

Nyr=zeros(12,1);
coldmean=zeros(12,1);
coldsigma=zeros(12,1);
hotmean=zeros(12,1);
hotsigma=zeros(12,1);
yr=[0:11]';

for i=[0:11]
    n = find( (t>=i*diy) & (t<(i+1)*diy) );
    tyr = t(n);
    dyr = d(n);
    Nyr(i+1) = length(tyr);
    dsort=sort(dyr);
    coldmean(i+1)=mean(dsort(1:10));
    hotmean(i+1)=mean(dsort(Nyr(i+1)-9:Nyr(i+1)));
    coldsigma(i+1)=std(dsort(1:10));
    hotsigma(i+1)=std(dsort(Nyr(i+1)-9:Nyr(i+1)));
end

Nyr0 = 10;
dmean0  = coldmean(1);
sigmad0 = coldsigma(1);
sigmam0 = sigmad0/sqrt(Nyr0);

Nyr11 = 10;
dmean11  = coldmean(12);
sigmad11 = coldsigma(12);
sigmam11 = sigmad11/sqrt(Nyr11);

test = (dmean0-dmean11) / sqrt( sigmam0^2 + sigmam11^2 );
m = ( sigmam0^2 + sigmam11^2 )^2 / ( ((sigmam0^4)/(Nyr0-1)) + ((sigmam11^4)/(Nyr11-1)) );
P = 1 - (tcdf(abs(test),m)-tcdf(-abs(test),m));

disp('coldest:')
disp(sprintf('year  1:  N %d mean %f sigmad %f sigmam %f', Nyr0, dmean0, sigmad0, sigmam0 ));
disp(sprintf('year 12:  N %d mean %f sigmad %f sigmam %f', Nyr11, dmean11, sigmad11, sigmam11 ));
disp(sprintf('t-est %f degrees of freedom  %f', test, m ));
disp(sprintf('P(|t|> %f = %f', abs(test), P ));

Nyr0 = 10;
dmean0  = hotmean(1);
sigmad0 = hotsigma(1);
sigmam0 = sigmad0/sqrt(Nyr0);

Nyr11 = 10;
dmean11  = hotmean(12);
sigmad11 = hotsigma(12);
sigmam11 = sigmad11/sqrt(Nyr11);

test = (dmean0-dmean11) / sqrt( sigmam0^2 + sigmam11^2 );
m = ( sigmam0^2 + sigmam11^2 )^2 / ( ((sigmam0^4)/(Nyr0-1)) + ((sigmam11^4)/(Nyr11-1)) );
P = 1 - (tcdf(abs(test),m)-tcdf(-abs(test),m));

disp('hotest:')
disp(sprintf('year 1:   N %d mean %f sigmad %f sigmam %f', Nyr0, dmean0, sigmad0, sigmam0 ));
disp(sprintf('year 12:  N %d mean %f sigmad %f sigmam %f', Nyr11, dmean11, sigmad11, sigmam11 ));
disp(sprintf('t-est %f degrees of freedom  %f', test, m ));
disp(sprintf('P(|t|> %f = %f', abs(test), P ));


