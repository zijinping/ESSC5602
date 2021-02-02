% problem 10_05
% modified from eda08_09
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

% process this image
j=9;
F=flipud(SST(:,:,j)'); % reorient so it plots as map

% truncate to even number of samples
[Nx, Ny] = size(F);
Nx = 2*floor(Nx/2);
Ny = 2*floor(Ny/2);
F = F(1:Nx,1:Ny);

% standard distance/wavenumber set up
Nkx=Nx/2+1;
Dx=1;
x=Dx*[0:Nx-1]';
xmax = Dx*(Nx-1);
kxmax=2*pi/(2*Dx);
Dkx=kxmax/(Nx/2);
kx=Dkx*[0:Nx/2,-Nx/2+1:-1]';

Nky=Ny/2+1;
Dy=1;
y=Dy*[0:Ny-1]';
ymax = Dy*(Ny-1);
kymax=2*pi/(2*Dy);
Dky=kymax/(Ny/2);
ky=Dky*[0:Ny/2,-Ny/2+1:-1]';

% unnormalized power spectral density
Ftt = fft2(F);
Fpsd = (abs(Ftt)).^2;

% put into more reasonable order
Fpsd2=zeros(Nx,Ny);
Fpsd2(Nx-Nkx+1:Nx,Ny-Nky+1:Ny) = Fpsd(1:Nkx,1:Nky);
Fpsd2(1:Nkx-2,Ny-Nky+1:Ny) = Fpsd(Nkx+1:Nx,1:Nky);
Fpsd2(1:Nkx-2,1:Nky-2) = Fpsd(Nkx+1:Nx,Nky+1:Ny);
Fpsd2(Nkx-1:Nx,1:Nky-2) = Fpsd(1:Nkx,Nky+1:Ny);

% plot
figure(1);
clf;
bw=0.9*(256-linspace(0,255,256)')/256;
colormap([bw,bw,bw]);
axis equal;
hold on;
axis ij;
set(gca,'XTick',[]); % turn off horizontal axis
set(gca,'YTick',[]); % turn off vertical axis
range=max(max(Fpsd2))-min(min(Fpsd2));
if( range==0 )
    range=1;
end
imagesc( (Fpsd2-min(min(Fpsd2)))/range );
axis tight;
xlabel('ky');
ylabel('kx');



