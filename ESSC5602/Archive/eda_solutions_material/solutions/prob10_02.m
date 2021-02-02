% prob10.2
% modified from interpolate_reynolds
% read and interpolate and output Renolds Channel Water Quality dataset
% columns are days, precip, air_temp, water_temp, salinity, turbidity,
% chloropjyl

% load data
D = load('reynolds_uninterpolated.txt');
[N, Ncols] = size(D);
Dt=1;
Di=D;

% interpolate with splines
for i = [2:Ncols]
    % insure data in first and last rows
    if( D(1,i) < 0 )
        D(1,i)=0;
    end
    if( D(N,i) < 0 )
        D(N,i)=0;
    end
    good = find( D(:,i)>=0 );
    x = D(good,1);
    y = D(good,i);
    z = spline(x,y,D(:,1));
    Di(:,i)=z;
end

t = Di(:,1);        % time, days
p = Di(:,2);        % precipitation, inches
at = Di(:,3);       % air temperature, deg C
wt = Di(:,4);       % water temperature, deg C
s = Di(:,5);        % salinity, practical salinity units
tb = Di(:,6);       % turbidity, 
c = Di(:,7);        % chlorophyl, micrograms per liter

figure(1);
clf;

subplot(6,1,1);
set(gca,'LineWidth',2);
hold on;
axis( [t(1), t(N), min(p), max(p)] );
plot( t, p, 'k-', 'LineWidth', 2 );
ylabel('precip');

subplot(6,1,2);
set(gca,'LineWidth',2);
hold on;
axis( [t(1), t(N), min(at), max(at)] );
plot( t, at, 'k-', 'LineWidth', 2 );
ylabel('T-air');

subplot(6,1,3);
set(gca,'LineWidth',2);
hold on;
axis( [t(1), t(N), min(wt), max(wt)] );
plot( t, wt, 'k-', 'LineWidth', 2 );
ylabel('T-water');

subplot(6,1,4);
set(gca,'LineWidth',2);
hold on;
axis( [t(1), t(N), min(s), max(s)] );
plot( t, s, 'k-', 'LineWidth', 2 );
xlabel('time, days');
ylabel('salinity');

subplot(6,1,5);
set(gca,'LineWidth',2);
hold on;
axis( [t(1), t(N), min(tb), max(tb)] );
plot( t, tb, 'k-', 'LineWidth', 2 );
xlabel('time, days');
ylabel('turbidity');

subplot(6,1,6);
set(gca,'LineWidth',2);
hold on;
axis( [t(1), t(N), min(c), max(c)] );
plot( t, c, 'k-', 'LineWidth', 2 );
xlabel('time, days');
ylabel('chlorophyl');

dlmwrite('reynolds_spline.txt', Di, 'delimiter','\t' );
