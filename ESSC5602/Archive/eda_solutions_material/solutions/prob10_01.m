% problem 9.01
% evaluate covariance formula

N=50;
M=zeros(N,2);
M(:,2) = [0:N-1]'/(N-1);
M(:,1) = flipud(M(:,2));
sigma2d = 1;

Cm = sigma2d * M * M';

eda_draw(' ',Cm,'caption Cm');
