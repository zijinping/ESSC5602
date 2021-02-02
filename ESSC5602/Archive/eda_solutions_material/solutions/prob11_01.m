% Problem eda11_01

N = 201;
xmin = -1;
xmax = 1;
x = xmin + (xmax-xmin)*[0:N]'/(N-1);

f = exp(x);
fa = 1+x;

figure(1);
clf;
subplot(2,1,1);
set(gca,'LineWidth',2);
hold on;
axis( [xmin, xmax, 0, 4] );
plot( x, f, '-', 'LineWidth', 3, 'Color', [0.8,0.8,0.8] );
plot( x, fa, 'k-', 'LineWidth', 2 );
xlabel('x');
ylabel('f(x)');

tol=0.05;
re = abs((f-fa)./f);
k = find(re<=tol);

subplot(2,1,2);
set(gca,'LineWidth',2);
hold on;
axis( [xmin, xmax, 0, 0.1] );
plot( x, re, 'k-', 'LineWidth', 2 );
plot( x(k), re(k), 'k-', 'LineWidth', 3 );
plot( [xmin, xmax]', [tol, tol]', 'k:', 'LineWidth', 2);
xlabel('x');
ylabel('relative error');



fprintf('the relative error is <= than %f\n', tol);
fprintf('in the interval %.2f < x < %.2f\n', x(k(1)), x(k(end)) );