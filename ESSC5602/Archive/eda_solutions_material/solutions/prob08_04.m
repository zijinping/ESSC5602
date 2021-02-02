% Problem 8.04

% size of sample vectors
M = 24;

% load data
D = load('brf_temp.txt');
Dt = 1/24;
traw = D(:,1);
draw = D(:,2);
Nraw = length(traw);

% calculate sampling interval, min, max of data
Dt=1/24;
diy = 365.25;
tmin=min(traw);
tmax=max(traw);
dmin=min(draw);
dmax=max(draw);

% The dataset set contains both missing and zero data.
% Step one is to fill in the missing data to zero. The
% array, m, is the index of the the available data into
% a timeseries that is equally-spaced with the sampling, Dx

D = floor(((traw(2:Nraw)-traw(1:Nraw-1))/Dt)+0.5);
m = [1, 1+cumsum(D)']';
Nc = max(m);
tc=tmin+Dt*[0:Nc-1]';
dc=zeros(Nc,1);
dc(m)=draw;

% now locate cold spikes and hot spikes and set to zero
rowindex=find( (dc<-40) | (dc>38) );
dc(rowindex)=0;

N = floor(Nc/M);  % number of whole days in time series 
S = zeros(N,M);   % sample matrix
dayno=zeros(N,1);   % keep tracck of time of good samples
Ngood=0;          % good days of data
% loop over days of data
for is = [1:N]
    
    % grab one day
    is1 = 1+(is-1)*M;
    is2 = M+(is-1)*M;
    dday = dc(is1:is2);
    
    % discard if it contains any zeros
    dzero=find(dday==0,1);
    if( length(dzero)>0 )
        continue;
    end
    
    % copy data to S, remove its mean
    Ngood = Ngood+1;
    S(Ngood,:)=dday';
    dday=dday-mean(dday);
    dayno(Ngood)=is;
    
end

S = S(1:Ngood,:); % truncate S, dayno
dayno = dayno(1:Ngood);
N = Ngood;

% S = C F by singular value decomposition
[U, SV, V ] = svd(S);
sv = diag(SV);
C = U*SV;
F = V';

Emax=max(max(abs(S-C*F)));

% plot singular values
figure(1);
clf;
set(gca,'LineWidth',2);
axis( [1, M, 0, max(sv) ] );
plot( [1:M], sv, 'k-', 'LineWidth', 2 );

% first 6 singular values are large
p=6;

% plot p factors
figure(2);
clf;
for ip = [1:p]
    f = F(ip,:)';
    subplot(p,1,ip);
    set(gca,'LineWidth',2);
    hold on;
    axis( [1, M, min(f), max(f) ] );
    plot( [1:M]', f, 'k-', 'LineWidth', 2 );
    xlabel('time, hours');
    ylabel(sprintf('f(%d)',ip));
end

% plot loadings
figure(3);
clf;
for ip = [1:p]
    c = C(:,ip);
    subplot(p,1,ip);
    set(gca,'LineWidth',2);
    hold on;
    axis( [dayno(1)/diy, dayno(N)/diy, min(c), max(c) ] );
    plot( dayno/diy, c, 'k-', 'LineWidth', 2 );
    xlabel('time, years');
    ylabel(sprintf('c(%d)',ip));
end
    


