% Problem 3.2

L = 100;
Dd = 1/L;
d = Dd*[0:L-1]';

lambda1 = 5;
lambda2 = 10;

p1 = lambda1 * exp( -lambda1*d );
p2 = lambda2 * exp( -lambda2*d );
eda_draw(p1,'caption 5',' ',' ',p2,'caption 10');

area1 = Dd*sum(p1);
area2 = Dd*sum(p2);

dbar1 = Dd*sum(d.*p1);
dbar2 = Dd*sum(d.*p2);

sigma2d1 = Dd*sum( (d-dbar1).^2 .* p1);
sigma2d2 = Dd*sum( (d-dbar2).^2 .* p2);

disp(sprintf('lambda: %f %f',lambda1, lambda2));
disp(sprintf('areas: %f %f',area1, area2));
disp(sprintf('means: %f %f',dbar1, dbar2));
disp(sprintf('variances: %f %f',sigma2d1, sigma2d2));
