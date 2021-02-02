% problem 5.6
% modified from eda05_08
% use generealized least square to fill in
% data gaps for data with the prior information
% that dd/dx=0 and dd/dy=0 (separately) in the
% interior and no prior information on the edges

% need this to communicate with afun()
clear F;
global F;

% load the data.  These data are only 'simulated'
% and were created by the script eda05_08
D = load('pressure.txt');
xobs=D(:,1);
yobs=D(:,2);
dobs=D(:,3);
N=length(dobs);
sigmad = 1.0;
rsigmad = 1/sigmad;
rsigmad2 = rsigmad^2;

% rows of grid are x
I = 41;
Dx = 1;
xmin = 0;
xmax = 40;
x = xmin + (xmax-xmin)*[0:I-1]'/(I-1);
rDx = 1/Dx;
rDx2 = rDx^2;

% columns of grid are y
J = 41;
Dy = 1;
ymin = 0;
ymax = 40;
y = ymin + (ymax-ymin)*[0:J-1]'/(J-1);
rDy = 1/Dy;
rDy2 = rDy^2;

% associate observations with points on the grid
% by scaling the range (xmin, xmax) to (1, N)
rowindex = 1+floor((I-1)*(xobs-xmin)/(xmax-xmin)+0.5);
colindex = 1+floor((J-1)*(yobs-ymin)/(ymax-ymin)+0.5);

% don't let any data fall outside the grid.  In this
% case, we just put it on the edge of the grid.  The
% alternative is to throw it away.
for i = [1:N]
    if( rowindex(i)<1 )
        rowindex(i)=1;
    elseif( rowindex(i)>I )
        rowindex(i)=I;
    end
    if( colindex(i)<1 )
        colindex(i)=1;
    elseif( colindex(i)>J )
        colindex(i)=J;
    end
end

% copy the data to the grid A for later use
Aobs = zeros(I,J);
for i = [1:N]
    Aobs(rowindex(i),colindex(i)) = dobs(i);
end

% one model parameter per grid point
M = I*J;

% Fm=f is LxM with L = N rows of data plus
% (I-1)*J rows of dd/dx plus
% I*(J-1) rows of dd/dy plus
% M rows of d=0
L = N+I*(J-1)+J*(I-1)+M;

% F will have no more than 3L non-zero elements

F = spalloc(L,M,3*L);
f = zeros(L,1);

% reminders on how to move from grid to vector & vice-versa
% k = (i-1)*J+j;
% m(k) = A(i,j);
% i = floor((k-1)/J)+1;
% j = k-(i-1)*J;
% B(i,j)=m(k);

% build F and f
% the variable q is used to keep track of the current row of F and f
q=1;

% first part of Fm=f: the data
for i = [1:N]
    k = (rowindex(i)-1)*J+colindex(i);
    F(q,k) = rsigmad;
    f(q) = rsigmad*dobs(i);
    dindex(i)=k;
    q = q+1;
end

% error of gradient prior information
sigmas = 0.2;
rsigmas = 1/sigmas;

% error of size prior information
sigmam = 1000;
rsigmam = 1/sigmam;

% second part, dd/dx is zero in the interior
for j = [1:J] % y
for i = [2:I] % x
    k = (i-1)*J+j;
    F(q,k) = -rsigmas*rDx;
    k = ((i+1)-1)*J+j;
    F(q,k) = rsigmas*rDx;
    q = q+1;
end
end

% third part, dd/dy is zero in the interior
for i = [1:I] % x
for j = [2:J] % y
    k = (i-1)*J+j;
    F(q,k) = -rsigmas*rDy;
    k = (i-1)*J+(j+1);
    F(q,k) = rsigmas*rDy;
    q = q+1;
end
end

% fourth part, m is zero everywhere
for i = [1:I] % x
for j = [2:J] % y
    k = (i-1)*J+j;
    F(q,k) = rsigmam;
    q = q+1;
end
end

% solve equation Fm=f in least squares sense
mest=bicg(@afun,F'*f,1e-10,3*L);

% reorganize estimated model parameters back into grid
Apre = zeros(I,J);
for k = [1:M]
    i = floor((k-1)/J)+1;
    j = k-(i-1)*J;
    Apre(i,j) = mest(k);
end

% plot results
eda_draw('  ',Aobs,'caption p-obs',Apre,'caption p-pre');

% compute weighted prediction error Ed
dpre=zeros(N,1);
for i = [1:N]
    k=dindex(i);
    dpre(i) = mest(k);
end
Ed = rsigmad2*(dobs-dpre)'*(dobs-dpre);
Ed




