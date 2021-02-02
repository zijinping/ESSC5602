% Problem eda11_02

% make a random set of (x,y) data
N = 1000;
xbar = 0.6;
sigmax = 0.01;
sigma2x = sigmax^2;
x = random('Normal',xbar,sigmax,N,1);
ybar = 0.7;
sigmay = 0.02;
sigma2y = sigmay^2;
y = random('Normal',ybar,sigmay,N,1);

% transform to (u,v)
u = x.*y;
v = x./y;

% estimate covariance two different ways and compare
Cn1 = cov(u,v);
fprintf('Numerical   s2u %f s2v %f cov %f\n', Cn1(1,1), Cn1(2,2), Cn1(1,2) );
M = [ybar, xbar; 1/ybar, -xbar/(ybar^2)];
Cn2 = M * diag([sigma2x,sigma2y]') * M';
fprintf('Approximate s2u %f s2v %f cov %f\n', Cn2(1,1), Cn2(2,2), Cn2(1,2) );


