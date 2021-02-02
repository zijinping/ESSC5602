% problem 11.07
% network to match non-linear response funtion,
% continuously retrained
% modified from eda11_15
% note that this script uses synthetic data computed from
% setup_prob11_07

clear all;

epsi1 = 1;
epsi2 = 0.05;
Niter=100;

% load dataset
D = load('datafor1107.txt');
Nx = length(D);

% break out into separate variables, metric conversion
% time in days
t = D(:,1);
% discharge in cubic meters per second
dobs = D(:,2)/35.3146;
% precipitation in millimeters per day
x = D(:,3)*25.4;

% normalization
tmin = min(t);
tmax = max(t);
xnorm = 1/max(x);
dobsnorm = 1/max(dobs);
x = x * xnorm;
dobs = dobs * dobsnorm;

% build the network; length Nc=10 filter
% initialize to an exponential decay
Nc = 10;
c0 = 2;
To = 3; 
f = exp(-[0:Nc-1]'/To);
f = [f(1)/1000;f(1:end-1)];
h = c0*f/sum(f);

% empty neural net
Lmax = 4;
Nmax = Nc ;
[N,w,b,a] = eda_net_init(Lmax,Nmax);

% neural net set to linear filter type
slope1 = 0.01;
slope2 = 0.02;
[N,w,b] = eda_net_linear( Nc, slope1, slope2, 0.0, h, N, w, b );

% index arrays for eacy access to weights, biases
[ Nb,Nw,layer_of_bias,neuron_of_bias,index_of_bias,layer_of_weight, ...
    r_neuron_of_weight,l_neuron_of_weight,index_of_weight ] = eda_net_map(N,w,b);

% train neural net on first Ntrain data
Ntrain = 501;
dpre = zeros(Ntrain,1);
% linearized data kernel
M = Nb+Nw;
G = zeros(Ntrain,M);

% first stage of training, on linear filter neural net
% using first Ntrain points in dataset
for iter=[1:Niter]

% loop over times
for kx=[1:Ntrain]
    
    % activities
    for j=[1:N(1)]
        k=kx-j+1;
        if( k>0 )
            a(j,1) = x(k);
        else
            a(j,1) = 0;
        end
    end
    a = eda_net_eval( N,w,b,a);
    dpre(kx,1) = a(1,4);
    [daLmaxdw,daLmaxdb] = eda_net_deriv( N,w,b,a );
    
    % form data kernel
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

% implement some damping to improve stability
% but relax it after a few iterations
if( iter==50 )
    epsi1=epsi1 / 10;
end

% damped least squares solution
Dd = dobs(1:Ntrain) - dpre;
Dm = (G'*G+epsi1*eye(M,M))\(G'*Dd);

% update neural net
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
    
% modifty the neural net by adding new connections
for i=[1]
for j=[1:N(1)]
    if( w(i,j,2) == 0 )
        w(i,j,2) = random('Normal',0,0.0001,1,1);
    end
end
end
for i=[1]
for j=[1:N(1)]
    if( w(i,j,3) == 0 )
        w(i,j,3) = random('Normal',0,0.0001,1,1);
    end
end
end

% remap
[ Nb,Nw,layer_of_bias,neuron_of_bias,index_of_bias,layer_of_weight, ...
    r_neuron_of_weight,l_neuron_of_weight,index_of_weight ] = eda_net_map(N,w,b);

% linearized data kernel
M = Nb+Nw;
G = zeros(Ntrain,M);

% predicted data
dpre = zeros(Ntrain,1);

% re-train modified neural network
for iter=[1:Niter]

% loop over x
for kx=[1:Ntrain]
    
    % activities
    for j=[1:N(1)]
        k=kx-j+1;
        if( k>0 )
            a(j,1) = x(k);
        else
            a(j,1) = 0;
        end
    end
    
    a = eda_net_eval( N,w,b,a);
    dpre(kx,1) = a(1,4);
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
    
end % next x

% implement some damping to improve stability
% but relax it after a few iterations
if( iter==50 )
    epsi2=epsi2/10;
end

% least squares solution
Dd = dobs(1:Ntrain) - dpre;
Dm = (G'*G+epsi2*eye(M,M))\(G'*Dd);

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

% now predict entire dataset
dpre = zeros(Nx,1);
error1 = zeros(Nx,1);
for kx=[1:Nx]
    for j=[1:N(1)]
        k=kx-j+1;
        if( k>0 )
            a(j,1) = x(k);
        else
            a(j,1) = 0;
        end
    end
    a = eda_net_eval( N,w,b,a);
    dpre(kx) = a(1,4);
    error1(kx) = abs(dobs(kx)-dpre(kx));
end

% plot the error
figure(1);
clf;
hold on;
set(gca,'LineWidth',3);
axis( [tmin, tmax, 0, max(error1) ] );
plot( t, error1, '-', 'LineWidth', 3, 'Color', [0.8,0.8,0.8] );
xlabel('time (days)');
ylabel('error');

% predict entire dataset again, with Kiter
% iterations of training after each time
Kiter=1;
dpre = zeros(Nx,1);
dpre2 = zeros(Ntrain,1);
error2 = zeros(Nx,1);
for kx=[1:Nx] % loop over days
    % predict today before re-training filter
        for j=[1:N(1)]
            k=kx-j+1;
            if( k>0 )
                a(j,1) = x(k);
            else
                a(j,1) = 0;
            end
        end
        a = eda_net_eval( N,w,b,a);
        dpre(kx,1) = a(1,4);
        error2(kx) = abs(dobs(kx)-dpre(kx));
    if( (kx+1) > Ntrain ) % train neural network on moving window
                          % of last Ntrain days
        for it = [1:Kiter] % iterations
          irow=0;
          for kkx=[kx-Ntrain+1:kx] % build data kernel
            irow=irow+1;
            for j=[1:N(1)]
                k=kkx-j+1;
                if( k>0 )
                    a(j,1) = x(k);
                else
                    a(j,1) = 0;
                end
             end
             a = eda_net_eval( N,w,b,a);
             dpre2(irow,1) = a(1,4);
             [daLmaxdw,daLmaxdb] = eda_net_deriv( N,w,b,a );
             for ib=[1:Nb]
                 nob=neuron_of_bias(ib);
                 lob=layer_of_bias(ib);
                 G(irow, ib ) = daLmaxdb(1, nob, lob );
             end
             for iw=[1:Nw]
                 rnow = r_neuron_of_weight(iw);
                 lnow = l_neuron_of_weight(iw);
                 low = layer_of_weight(iw);
                 G(irow, Nb+iw ) = daLmaxdw( 1, rnow, lnow, low );
             end
          end % end build data kernel
          % least squares solution
          Dd = dobs(kx-Ntrain+1:kx) - dpre2;
          Dm = (G'*G+epsi2*eye(M,M))\(G'*Dd);
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
    end % if training
    end % next day

% plot the error
figure(1);
plot( t, error2, 'k-', 'LineWidth', 2 );

