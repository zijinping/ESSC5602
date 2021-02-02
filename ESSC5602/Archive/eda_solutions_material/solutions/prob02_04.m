% problem 2.4
% adaped version of eda02_03
% make a plot of a section of data centered at a mouse-click
% this version presumes that the whole dataset is already
% plotted in Figure 1, and that the data are in vectors t and d

% read the data
D=load('neuse.txt');
t=D(:,1);
d=D(:,2);
Ns=size(D);
N=Ns(1);
M=Ns(2);

% plot the data
figure(1);
clf;
set(gca,'Linewidth',2);
hold on;
plot(t,d,'k-','Linewidth',2);
title('Neuse River Hydrograph');
xlabel('time, days');
ylabel('discharge, cfs');

% plot section based on cursor click
w=100; % width of new plot in samples
figure(1);
hold on;
[tc, dc] = ginput(1);  % detect mouse click
% i=find((t>=tc),1 ); % find index i of x(i) corresponding to click
% note early versions of MatLab need this instead: i=find((t>=tc)); i=i(1);
i=find((t>=tc)); i=i(1);
figure(2);
clf;
set(gca,'LineWidth',2);
hold on;
plot( t(i-w/2:i+w/2), d(i-w/2:i+w/2), 'k-','LineWidth',2 );
plot( t(i-w/2:i+w/2), d(i-w/2:i+w/2), 'k.','LineWidth',2 );
title('Neuse River Hydrograph');
xlabel('time in days');
ylabel('discharge, cfs');
figure(1);
