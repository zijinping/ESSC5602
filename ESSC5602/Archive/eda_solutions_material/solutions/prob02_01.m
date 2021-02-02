% Problem 2.1

% read data from file
D=load('brf_temp.txt');
t=D(:,1);
d=D(:,2);

% convert tine to years
daysinyear = 365.25;
ty = t/daysinyear;

% plot
figure(1);
clf;
set(gca,'LineWidth',2);
axis( [0, 14, -60, 40] );
plot(ty,d);
xlabel('time, years');
ylabel('discharge, cfs');



