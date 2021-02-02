% Problem 3.1

L = 100;
Dd = 1/L;
d = Dd*[0:L-1]';
p = 2*(1-d);

area = Dd*sum(p);
area

dbar = Dd*sum(d.*p);
dbar

sigma2d = Dd*sum( (d-dbar).^2 .* p);
sigma2d