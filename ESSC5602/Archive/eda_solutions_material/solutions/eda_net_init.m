function [N,w,b,a] = eda_net_init(Lmax,Nmax)
% initialize a network; produces only arrays of zeros
% input:
%        Lmax: number of layers
%        Nmax: maximum number of neurons on any layer
% output:
%        N: neurons on a layer array
%        w: weights
%        b: biases
%        a: activities

% N(i): column-vector with number of neurons in each layer;
N = zeros(1,Lmax);

% biases: b(1:N(i),i) is a column-vector that gives the
% biases of the N(i) neurons on the i-th layer
b = zeros( Nmax, Lmax );

% weights: w(1:N(i),1:N(i-1),i) is a matrix that gives the
% weights between the N(i) neurons in the i-th layer and the
% N(i-1) neurons on the (i-1)-th layer. The weights assigned
% to the first layer are never used and should be set to zero
w = zeros( Nmax, Nmax, Lmax );

% activity: a(1:N(i),i) is a column-vector that gives the
% activitites of the N(i) neurons on the i-th layer
a = zeros( Nmax, Lmax );


end
