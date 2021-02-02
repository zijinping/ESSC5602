% Problem 1.3

N=4;
M=zeros(N,N);
M(1,1)=1;
M(1,2)=-1;
M(2,2)=1;
M(2,3)=-1;
M(3,3)=1;
M(3,4)=-1;
M(4,4)=1;
M

y=[1, 2, 3, 5]';
y

x = M\y;
x

error = y-M*x;
error

