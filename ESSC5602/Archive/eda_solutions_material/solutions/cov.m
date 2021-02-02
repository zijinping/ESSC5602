N=10;

r = zeros(1,N);
r(1)=1;
r(2)=1;
c = zeros(N,1);
c(1)=1;
G = toeplitz( r, c );

C = inv( G'*G );
