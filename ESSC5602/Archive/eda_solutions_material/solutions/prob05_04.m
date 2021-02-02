% problem 5.4
% adapted from eda05_08
% use generealized least square to fill in
% data gaps for data that satisfy Laplace's
% equation in the in the interior of a grid,
% and the normal derivative is zero on the
% edges of the grid

% need this to communicate with afun()
clear F;
global F;

figure(1);
clf;

figure(2);
clf;

KAPPA =  [0.05, 0.10, 0.20]';
SIGMAD = [0.05, 0.10, 1.00]';

NKAPPA = length(KAPPA);
NSIGMAD = length(SIGMAD);

for ik = [1:NKAPPA]
for is = [1:NSIGMAD]

% load the data.  These data are only 'simulated'
% and were created by the script eda05_08
D = load(sprintf('pressure%d_%d.txt',ik,is));
xobs=D(:,1);
yobs=D(:,2);
dobs=D(:,3);
N=length(dobs);
sigmad = SIGMAD(is);
kappa = KAPPA(ik);
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

% Fm=f is LxM with L = N rows of data plus M rows of prior information
L = N+M;

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

% error of prior information
sigmah = 1000;
rsigmah = 1/sigmah;

% second part, laplaces equation in interior
% second derivative is a (1,-2,1)/Dx^2 in x plus a (1,-2,1)/Dy^2 in y
for i = [2:I-1]
for j = [2:J-1]
    k = (i-1)*J+j;
    F(q,k) = rsigmah*(-2*rDx2-2*rDy2);
    k = ((i-1)-1)*J+j;
    F(q,k) = rsigmah*rDx2;
    k = ((i+1)-1)*J+j;
    F(q,k) = rsigmah*rDx2;
    k = (i-1)*J+(j-1);
    F(q,k) = rsigmah*rDy2;
    k = (i-1)*J+(j+1);
    F(q,k) = rsigmah*rDy2;
    f(q) = 0.0;
    q = q+1;
end
end

% third part: top edge
for j = [2:J-1]
    i = 1;
    k = (i-1)*J+j;
    F(q,k) = -rsigmah*rDx;
    k = ((i+1)-1)*J+j;
    F(q,k) = rsigmah*rDx;
    f(q) = 0.0;
    q = q+1;
end

% fourth part: bottom edge
for j = [2:J-1]
    i = I;
    k = (i-1)*J+j;
    F(q,k) = rsigmah*rDx;
    k = ((i-1)-1)*J+j;
    F(q,k) = -rsigmah*rDx;
    f(q) = 0.0;
    q = q+1;
end

% fifth part: left edge
for i = [2:I-1]
    j = 1;
    k = (i-1)*J+j;
    F(q,k) = -rsigmah*rDy;
    k = (i-1)*J+(j+1);
    F(q,k) = rsigmah*rDy;
    f(q) = 0.0;
    q = q+1;
end;

% sixth part: right edge
for i = [2:I-1]
    j = J;
    k = (i-1)*J+j;
    F(q,k) = rsigmah*rDy;
    k = (i-1)*J+(j-1);
    F(q,k) = -rsigmah*rDy;
    f(q) = 0.0;
    q = q+1;
end; 

% seventh part: four corners
% top left
rDz = 1/sqrt((Dx^2)+(Dy^2));
i=1;
j=1;
k = (i-1)*J+j;
F(q,k) = -rsigmah*rDz;
k = ((i+1)-1)*J+(j+1);
F(q,k) = rsigmah*rDz;
f(q)=0.0;
q=q+1;
% top right
i=1;
j=J;
k = (i-1)*J+j;
F(q,k) = rsigmah*rDz;
k = ((i+1)-1)*J+(j-1);
F(q,k) = -rsigmah*rDz;
f(q)=0.0;
q=q+1;
% bottom left
i=I;
j=1;
k = (i-1)*J+j;
F(q,k) = -rsigmah*rDz;
k = ((i-1)-1)*J+(j+1);
F(q,k) = rsigmah*rDz;
f(q)=0.0;
q=q+1;
% bottom right
i=I;
j=1;
k = (i-1)*J+j;
F(q,k) = rsigmah*rDz;
k = ((i-1)-1)*J+(j-1);
F(q,k) = -rsigmah*rDz;
f(q)=0.0;

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
figure(1);
subplot(NKAPPA,NSIGMAD,ik+(is-1)*NKAPPA);
set(gca,'LineWidth',2);
axis([0, 1, 0, 1]);
hold on;
axis ij;
axis off;
bw=0.9*(256-[0:255]')/256;
colormap([bw,bw,bw]);
imagesc( [0, 1], [0, 1], Aobs);
title(sprintf('k %.2f sd %.2f', kappa, sigmad));

figure(2);
subplot(NKAPPA,NSIGMAD,ik+(is-1)*NKAPPA);
set(gca,'LineWidth',2);
axis([0, 1, 0, 1]);
hold on;
axis ij;
axis off;
bw=0.9*(256-[0:255]')/256;
colormap([bw,bw,bw]);
imagesc( [0, 1], [0, 1], Apre);
title(sprintf('k %.2f sd %.2f', kappa, sigmad));

% compute weighted prediction error Ed
dpre=zeros(N,1);
for i = [1:N]
    k=dindex(i);
    dpre(i) = mest(k);
end
Ed = rsigmad2*(dobs-dpre)'*(dobs-dpre);
disp(sprintf('kappa %f sigmad %f Ed %f',kappa,sigmad,Ed));

end
end




