function [Nb,Nw,layer_of_bias,neuron_of_bias,index_of_bias,layer_of_weight, ...
    r_neuron_of_weight,l_neuron_of_weight,index_of_weight ] = eda_net_map(N,w,b)
% builds look-up tables to allow sequential ordering
% of non-trivual biases and weights
% used to build list of model parameters used in training
% input:
%        N: neurons on layer array
%        w: weights
%        b: biases
% output
%        Nb: number of biases
%        Nw: Number of non-zero weights
%            (use near-zero weights to fool it)
%        layer_of_bias, vector(i) specifying layer of bias i
%        neuron_of_bias, vector(i) specifying neuron of bias i
%        index_of_bias, 2d array(i,j) specifying number of bias
%            of layer j, neuron i
%         layer_of_weight, vector(i) specifying layer of weight i
%         r_neuron_of_weight, vector(i) specifying right neuron of weight i
%         l_neuron_of_weight, vector(i) specifying left neuron of weight i
%         index_of_weight, array(i,j,k) specifying number of weight(i,j,k)

    Nmax = max(N);
    Lmax = length(N);

layer_of_bias = zeros(Lmax*Nmax,1);
neuron_of_bias = zeros(Lmax*Nmax,1);
index_of_bias = zeros(Nmax,Lmax);
Nb=0;  % number of nontrivial biases
for k=[2:Lmax]
    for j=[1:N(k)]
        Nb = Nb+1;
        layer_of_bias(Nb)=k;
        neuron_of_bias(Nb)=j;
        index_of_bias(j,k)=Nb;
    end
end
Nw=0;  % number of nontrivial weights
layer_of_weight = zeros(Lmax*Nmax*Nmax,1);
r_neuron_of_weight = zeros(Lmax*Nmax*Nmax,1);
l_neuron_of_weight = zeros(Lmax*Nmax*Nmax,1);
index_of_weight = zeros(Nmax,Nmax,Lmax);


for k=[2:Lmax]
    for i=[1:N(k)]
    for j=[1:N(k-1)]
        if( w(i,j,k) ~= 0 )
            Nw = Nw+1;
            layer_of_weight(Nw)=k;
            r_neuron_of_weight(Nw)=i;
            l_neuron_of_weight(Nw)=j;
            index_of_weight(j,i,k)=Nw;
        end
    end
    end
end

end

