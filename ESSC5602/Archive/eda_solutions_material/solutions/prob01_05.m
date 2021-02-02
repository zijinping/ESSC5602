% Problem 1.5

% read data from file
D=load('neuse.txt');
t=D(:,1);
d=D(:,2);

dbad = find(d<0);
N=length(dbad);
disp(sprintf('%d data points with negative discharge',N));

