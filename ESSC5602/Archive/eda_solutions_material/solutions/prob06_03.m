% problem 6.3
% modified from eda06_14
% make a random noise dataset

% lengths of running averages
LL = [3 5 7]';
NLL = length(LL);

% set up time
N=1024;
Dt=1.0;
T=N*Dt; 
t=Dt*[0:N-1]';

% set up frequency
fmax=1/(2.0*Dt);
df=fmax/(N/2);
f=df*[0:N/2,-N/2+1:-1]'; 
Nf=N/2+1;
dw=2*pi*df;
w=dw*[0:N/2,-N/2+1:-1]';
Nw=Nf;

% make vector d of random noise
dbar=0.0;
sigmad=1.0;
d = random('Normal',dbar, sigmad, N, 1);

% transform of d
dtilde = Dt * fft(d);

% power spectral density, + frequencies only
s2 = (2/T) * abs(dtilde(1:Nf)).^2;

% plot power spectral density
figure(1);
clf;
subplot(NLL+1,1,1);
set(gca,'LineWidth',2);
hold on;
axis( [0, fmax, 0, max(s2)] );
plot( f(1:Nf), s2, 'k-', 'LineWidth', 2 );
xlabel('frequency, Hz');
ylabel('s2');

for j = [1:NLL]
    
L = LL(j);
L2 = floor(L/2);

a=zeros(N,1);
for i=[1:N]
    l1 = i-L2;
    if( l1<1 )
        l1=1;
    end
    l2 = i+1;
    if( l2>N )
        l2=N;
    end
    a(i) = mean(d(l1:l2));
end

% transform of a
atilde = Dt * fft(a);

% power spectral density, + frequencies only
sa2 = (2/T) * abs(atilde(1:Nf)).^2;

% plot power spectral density
subplot(NLL+1,1,1+j);
set(gca,'LineWidth',2);
hold on;
axis( [0, fmax, 0, max(s2)] );
plot( f(1:Nf), sa2, 'k-', 'LineWidth', 2 );
xlabel('frequency, Hz');
ylabel(sprintf('%d pt avg',L));

end
