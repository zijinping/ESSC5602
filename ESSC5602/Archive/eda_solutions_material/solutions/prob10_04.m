% problem 10_04
% modified from eda10_08

% standard distance/wavenumber set up
Nx=32;
Nkx=Nx/2+1;
Dx=1;
x=Dx*[0:Nx-1]';
xmax = Dx*(Nx-1);
kxmax=2*pi/(2*Dx);
Dkx=kxmax/(Nx/2);
kx=Dkx*[0:Nx/2,-Nx/2+1:-1]';

Ny=32;
Nky=Ny/2+1;
Dy=1;
y=Dy*[0:Ny-1]';
ymax = Dy*(Ny-1);
kymax=2*pi/(2*Dy);
Dky=kymax/(Ny/2);
ky=Dky*[0:Ny/2,-Ny/2+1:-1]';

% cosine wave
theta = (pi/180)*30;
t1 = [cos(theta), sin(theta)]';;
kr1 = kxmax/4;
F = zeros(Nx,Ny);
for n = [1:Nx]
for m = [1:Ny]
    F(n,m) = cos( kr1*t1(1)*x(n) + kr1*t1(2)*y(m) );
end
end

% Hamming filter
Wx=0.54-0.46*cos(2*pi*[0:Nx-1]'/(Nx-1));
Wy=0.54-0.46*cos(2*pi*[0:Ny-1]'/(Ny-1));
F = F .* (Wx*Wy'); % note use of outer product

% amplitude spectral density
Ftt = fft2(F);
Fasd = ((Dx^2)/xmax) * ((Dy^2)/ymax) * sqrt((abs(Ftt).^2));

% put into more reasonable order
Fasd2=zeros(Nx,Ny);
Fasd2(Nx-Nkx+1:Nx,Ny-Nky+1:Ny) = Fasd(1:Nkx,1:Nky);
Fasd2(1:Nkx-2,Ny-Nky+1:Ny) = Fasd(Nkx+1:Nx,1:Nky);
Fasd2(1:Nkx-2,1:Nky-2) = Fasd(Nkx+1:Nx,Nky+1:Ny);
Fasd2(Nkx-1:Nx,1:Nky-2) = Fasd(1:Nkx,Nky+1:Ny);

eda_draw('  ', F, 'caption F',' ',' ',' ',  Fasd, 'caption Ftt',' ',' ',' ', Fasd2, 'caption Ftt' );


