% problem 11.05
% simple least squares fit of gaussian with a single tower

% grid of points in (x,y)
Nx = 21;
Ny = 21;
xmax = pi;
ymax = pi;
xmin = 0;
ymin = 0;
x = xmin + (xmax-xmin)* [0:Nx-1]'/(Nx-1);
y = ymin + (ymax-ymin)* [0:Ny-1]'/(Ny-1);

% (x,y) index arrays for easy access
kofij=zeros(Nx,Ny);  % linear index k of a particular (i,j)
iofk=zeros(Nx*Ny,1); % row index of a particular k
jofk=zeros(Nx*Ny,1); % col index of a particular k
Nxy=0;
for i=[1:Nx]
for j=[1:Ny]
    Nxy=Nxy+1;
    kofij(i,j)=Nxy;
    iofk(Nxy)=i;
    jofk(Nxy)=j;
end
end

% Exemplary function sin(x)*cos(y)
dobs = zeros(Nx,Ny);
dobsvec = zeros(Nxy,1);
xbar = 2.4;
ybar = 2.6;
amp = 5.0;
sigmax = 1.2;
sigmay = 0.8;
cx = 1/(2*sigmax*sigmax);
cy = 1/(2*sigmay*sigmay);
for i=[1:Nx]
for j=[1:Ny]
    dobs(i,j) = sin(x(i))*cos(y(j));
    dobsvec(kofij(i,j))=dobs(i,j);
end
end

% build the test network for a two 2D towers
Nc = 2;     % two towers
Xc = [pi/2 pi/2]; % Xc: the center of the towers
Yc = [pi/4 3*pi/4];
Dx = [1 1];   % Dx: the width of the towers
Dy = [1 1];
h = [1 -1];    % h: the height of the towers
slope = 5/4;  % slope: slope of the step functions

% ---Neural Net definitions---

% Nmax: number of layers;
Lmax = 4;

% Lmax: maximum number of neurons in any layer
Nmax = Nc * 4;

[N,w,b,a] = eda_net_init(Lmax,Nmax);

% daLmaxdb(i,j,k): the change in the activity of the i-th neuron in layer Lmax
% due to a change in the bias of the of the j-th neuron in the layer k
daLmaxdb = zeros(Nmax,Nmax,Lmax);

% daLmaxdw(i,j,k,l): the change in the activity of the i-th neuron in layer Lmax
% due to a change in the weight connecting the j-th neuron in the layer l
% with the k-th neuron of layer l-1
daLmaxdw = zeros(Nmax,Nmax,Nmax,Lmax);

[N,w,b] = eda_net_2dtower( Nc, slope, Xc, Yc, Dx, Dy, h, N, w, b );

w(1:2,2,2)=[.0001,.0001];
w(3:4,1,2)=[.0001,.0001];

figure(2);
clf;
eda_net_view(N,w,b);

[ Nb,Nw,layer_of_bias,neuron_of_bias,index_of_bias,layer_of_weight, ...
    r_neuron_of_weight,l_neuron_of_weight,index_of_weight ] = eda_net_map(N,w,b);

% linearized data kernel
M = Nb+Nw;
G = zeros(Nxy,M);

% predicted data
dpre = zeros(Nx,Ny);
dprevec = zeros(Nxy,1);

% initial predicted data
d0 = zeros(Nx,Ny);
d0vec = zeros(Nxy,1);

Niter=20;
for iter=[1:Niter]

% loop over (x,y)
for kx=[1:Nx]
for ky=[1:Ny]
    
    % activities
    a(1:N(1),1) = [x(kx);x(ky)];
    a = eda_net_eval( N,w,b,a);
    dpre(kx,ky) = a(1,4);
    dprevec(kofij(kx,ky))=dpre(kx,ky);
    
    [daLmaxdw,daLmaxdb] = eda_net_deriv( N,w,b,a );
    
    % biases first
    for ib=[1:Nb]
        nob=neuron_of_bias(ib);
        lob=layer_of_bias(ib);
        G(kofij(kx,ky), ib ) = daLmaxdb(1, nob, lob );
    end
    % weights second
    for iw=[1:Nw]
        rnow = r_neuron_of_weight(iw);
        lnow = l_neuron_of_weight(iw);
        low = layer_of_weight(iw);
        G(kofij(kx,ky), Nb+iw ) = daLmaxdw( 1, rnow, lnow, low );
    end
    
end % next (x,y)
end

Dd = dobsvec - dprevec;
 
if( iter==1 )
    d0 = dpre;
    d0vec = dobsvec;
    error = Dd'*Dd;
    fprintf('starting error %.2f\n', error );
end

% implement some damping to improve stability
% but relax it after a few iterations
if( iter<5 )
    epsi=1.0;
elseif (iter<10)
    epsi=0.1;
else
    epsi=0.001;
end
Dm = (G'*G+epsi*eye(M,M))\(G'*Dd);

% add in bias perturbations
for ib=[1:Nb]
    nob=neuron_of_bias(ib);
    lob=layer_of_bias(ib);
    b(nob, lob) = b(nob, lob) + Dm(ib);
end
% add in weights perturbations
for iw=[1:Nw]
    rnow = r_neuron_of_weight(iw);
    lnow = l_neuron_of_weight(iw);
    low = layer_of_weight(iw);
    w(rnow, lnow, low) = w(rnow, lnow, low) + Dm( Nb+iw );
end

end % next iteration
    
% estimate exemplary function
for kx=[1:Nx]
for ky=[1:Ny]
    
    % activities
    a(1:N(1),1) = [x(kx);x(ky)];
    a = eda_net_eval( N,w,b,a);
    dpre(kx,ky) = a(1,4);
    dprevec(kofij(kx,ky))=dpre(kx,ky);
end
end

% plot exemplary function
eda_draw( dobs, 'caption d-obs', d0, 'caption d-start',dpre, 'caption d-pre');

Dd = dobsvec - dprevec;
error = Dd'*Dd;
fprintf('ending error %.2f\n', error );

