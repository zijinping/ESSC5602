% problem 8.3
% modified from eda09_09
% CAC Pacific Sea Surface Temperature (SST) dataset

% load the data
D=load('cac_sst.txt');

% set up sizes of arrays
[I, J] = size(D);
NBLOCK=85;
NLON = NBLOCK-1;
NLAT = J-1;
IMAGES = floor(I/NBLOCK);
Dt = 1/12;

% arrays for SST, date, lats, lons
SST = zeros(NLON,NLAT,IMAGES);
thedate = zeros(IMAGES,1); % in month.year format
theyear = zeros(IMAGES,1); % year extracted from thedate
themonth = zeros(IMAGES,1); % month extracted from thedate
lats=zeros(NLAT,1);
lons=zeros(NLON,1);

% cut up data into SST arrays
lats = D(1,2:J)';
lons = D(2:NBLOCK,1);
for j=[1:IMAGES]
    k1 = 2+(j-1)*NBLOCK;
    k2 = k1+NLON-1;
    thedate(j)=D(k1-1,1);
    themonth(j)=floor(thedate(j));
    theyear(j)=floor(10000*(thedate(j)-themonth(j))+0.1);
    SST(:,:,j) = D(k1:k2,2:J);
end

% some definitions related to plotting
MONTHS=12;
YEARS=ceil(IMAGES/12);
YBLOCKSIZE = 6;
YBLOCKS = ceil(YEARS/YBLOCKSIZE);

% fold out into array
SSTV = zeros( IMAGES, NLAT*NLON );
for j=[1:IMAGES]
    for p = [1:NLAT]
    for q = [1:NLON]
        k = p + (q-1)*NLAT;
        SSTV(j,k) = SST( q, p, j );
    end
    end
end

% singular value decomposition
[U,S,V] = svd(SSTV,0);
s=diag(S);
Ns=length(s);
C=U*S;
F=V';

% number of high-amplitde factors
Nf=12;

% truncate to even length
IMAGES2 = 2*floor(IMAGES/2);
T = IMAGES2*Dt;

% frequenct definitions
freqmax=1/(2.0*Dt);
dfreq=freqmax/(IMAGES2/2);
freq=dfreq*[0:IMAGES2/2,-IMAGES2/2+1:-1]'; 
Nfreq=IMAGES2/2+1;
freq = freq(1:Nfreq);

% compute and plot PSD of loading timeseries
figure(1);
clf;
for f = [1:Nf]
    
    % power spectral density
    c = C(1:IMAGES2,f)';
    ctilde = Dt * fft(c);
    sc2 = (2/T) * abs(ctilde).^2;
    sc2 = sc2(1:Nfreq);
    
    % plot
    subplot(Nf/2, 2, f);
    set(gca,'LineWidth',2);
    hold on;
    axis( [ 0, 3, 0, max(sc2) ] );
    plot( freq, sc2, 'k-', 'LineWidth', 2 );
    if( (f==1) || (f==2) )
        title( 'PSD of loadings' );
    end
    if( (f==(Nf-1)) || (f==Nf) )
        xlabel( 'frequency, cycles per year' );
    end
    ylabel(sprintf('C(%d)',f));
end
