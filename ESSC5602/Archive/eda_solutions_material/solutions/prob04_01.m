% Problem 4.1

M = 40;

% set up matrix, G
r = zeros(1,M);
r(1)=1;
c = zeros(M,1);
c(1)=1;
c(2)=1;
G=toeplitz(c,r);

% propagate error
sigma2d = (1)^2;
Cm = sigma2d*inv(G'*G);
sigma2m = diag(Cm);
sigmam = sqrt(sigma2m);

figure(1);
clf;
set(gca,'LineWidth',2);
plot( [1:M], sigmam, 'k-', 'LineWidth', 2 );
xlabel('index j');
ylabel('sigma m');

