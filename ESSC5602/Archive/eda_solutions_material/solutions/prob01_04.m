% Problem 1.4

N=50;
r=zeros(1,N);
r(1)=1;
r(2)=-1;
c=zeros(N,1);
c(1)=1;
M=toeplitz(c,r);
M
