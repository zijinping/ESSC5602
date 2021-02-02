% problem 7.04
% modifield from eda07_06

% i used as imaginary unit!
clear i;

% define a filter
goriginal = [6, 5.5, 5.3, 5, 4, 3, 2]';
N = length(goriginal);
gN = g(N);

figure(1);
clf;

Nloop=5;

for iloop = [1:Nloop]
    
g = goriginal;
g(1) = g(1)*(1-0.13*(iloop-1));

% find roots of g
r = roots(flipud(g));

% check size of roots
badroots=0;
for j = [1:N-1]
    if( abs(r(j)) <= 1 )
        badroots=badroots+1;
    end
end

% rebuild filter from its roots
gp = zeros(1,1);
gp(1)=gN;
for j = [1:N-1]
    gp = conv(gp,[-r(j),1]');
end
% there should be no imiginary part
gp = real(gp);

% check that it worked
E1 = (g-gp)'*(g-gp);

% construct inverse filter, gi, of length Ni
Ni = 50;
gi = zeros(Ni,1);
gi(1)=1/gN;
% filter cascade, one filter per loop
for j = [1:N-1]
    % construct inverse filter of a length-2 filter
    tmp = zeros(Ni,1);
    tmp(1) = 1;
    for k = [2:Ni]
        tmp(k)=tmp(k-1)/r(j);
    end
    tmp = -tmp/r(j);
    gi = conv(gi,tmp);
    gi=gi(1:Ni);
end
% delete imaginary part (which should be zero)
gi = real(gi);

% plot results
subplot(Nloop,1,iloop);
set(gca,'LineWidth',2);
hold on;
axis( [0, Ni, min(gi), max(gi)] );
plot( [1:Ni]', gi, 'k-', 'LineWidth', 2);
xlabel('element j');
ylabel('g-inv');
title(sprintf('g(1)=%.2f %d roots not outside',g(1),badroots));
end


    
