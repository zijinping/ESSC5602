% Problem 3.3

L = 100;
DE = 2/L;
E = DE*[1:L]';

p = (4/sqrt(2*pi)) * E .* exp( -0.5*(E.^4) );

eda_draw(p,' ',' ','caption p(E)');

area = DE*sum(p);

[dmax, imode] = max(p);
themode = E(imode);

P = DE * cumsum(p);
themedian = E(find(P>0.50,1));

themean = DE*sum( E .* p);

sigma2E = DE*sum( (E-themean).^2 .* p);

disp(sprintf('area: %f',area));
disp(sprintf('mode: %f',themode));
disp(sprintf('median: %f',themedian));
disp(sprintf('mean: %f',themean));
disp(sprintf('variance: %f',sigma2E));
