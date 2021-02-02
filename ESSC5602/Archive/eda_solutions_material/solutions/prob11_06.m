% problem 11.6
%
% Use a neural net to approximate an parabolic function y(x)
% by superimposing 1D towers

% ---Polynomial function definitions---

% Nc: the number of centers in the approximation
Nc= 2;
% xcmax: the maximum value of the centers in the approximation
xcmax = 1;
% xcmin: the minimum value of the centers in the approximation
xcmin = 0;
% xc: the centers in the approximation
xc = xcmin + (xcmax-xcmin)* [0:Nc-1]'/(Nc-1);


% yc: the true value of the polynomial function
yc = 0.2 + 1.75*xc.^2;

% ---Step function definitions---

% Xc: the center of the step function
Xc = xc;

% Dx: the width of the step function
Dx = ones(Nc,1)*(xcmax-xcmin)/(Nc-1);

% h: the height of the step function
h = yc;

% slope: slope of the step functiona
slope = 50/4;

% ---Neural Net definitions---

% Nmax: number of layers;
Lmax = 3;

% Lmax: maximum number of neurons in any layer
Nmax = Nc * 2;

[N,w,b,a] = eda_net_init(Lmax,Nmax);

[N,w,b] = eda_net_1dtower( Nc, slope, Xc, Dx, h, N, w, b );

figure(1);
clf;
eda_net_view(N,w,b);

% ---Picking the points to plot---

% Nx: the number of points we want to plot
Nx= 100;
% xmax: the maximum value of the points we want to plot
xmax = xcmax + .1;
% xmin: the minimum value of the points we want to plot
xmin = xcmin - .1;
% x: the points we are plotting
x = xmin + (xmax-xmin)* [0:Nx-1]'/(Nx-1);
dobs = 0.2 + 1.75*x.^2;
m= zeros(Nx,1);



[ Nb,Nw,layer_of_bias,neuron_of_bias,index_of_bias,layer_of_weight, ...
    r_neuron_of_weight,l_neuron_of_weight,index_of_weight ] = eda_net_map(N,w,b);

% linearized data kernel
M = Nb+Nw;
G = zeros(Nx,M);

% predicted data
dpre = zeros(Nx,1);

% initial predicted data
d0 = zeros(Nx,1);



% ---Calculating the Activities---

Niter=50;
for iter=[1:Niter]

% loop over (x,y)
for kx=[1:Nx]
    
    % activities
    a(1:N(1),1) = [x(kx)];
    a = eda_net_eval( N,w,b,a);
    dpre(kx,1) = a(1,3);
    
    [daLmaxdw,daLmaxdb] = eda_net_deriv( N,w,b,a );
    
    % biases first
    for ib=[1:Nb]
        nob=neuron_of_bias(ib);
        lob=layer_of_bias(ib);
        G(kx, ib ) = daLmaxdb(1, nob, lob );
    end
    % weights second
    for iw=[1:Nw]
        rnow = r_neuron_of_weight(iw);
        lnow = l_neuron_of_weight(iw);
        low = layer_of_weight(iw);
        G(kx, Nb+iw ) = daLmaxdw( 1, rnow, lnow, low );
    end
    
end % next (x)

Dd = dobs - dpre;
 
if( iter==1 )
    d0 = dpre;
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
    
% loop over (x,y)
for kx=[1:Nx]
    % activities
    a(1:N(1),1) = [x(kx)];
    a = eda_net_eval( N,w,b,a);
    dpre(kx) = a(1,3);
end


figure(2);
clf;
subplot(2,1,1);
set(gca,'LineWidth',2);
hold on;
axis xy;
axis( [xmin xmax -1 2] );
plot(x,dobs,'Color',[.8,.8,.8],'LineWidth',3);
plot(x,d0,'k-','LineWidth',2); 
xlabel('x');
ylabel('d(x)');

% ---Plot true value---

dobs= zeros(Nx,1);
dobs = 0.2 + 1.75*x.^2;

subplot(2,1,2);
set(gca,'LineWidth',2);
hold on;
axis xy;
axis( [xmin xmax -1 2] );
plot(x,dobs,'Color',[.8,.8,.8],'LineWidth',3); 
plot(x,dpre,'k-','LineWidth',2);
xlabel('x');
ylabel('x(x)');

Dd = dobs - dpre;
error = Dd'*Dd;
fprintf('ending error %.2f\n', error );


