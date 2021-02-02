function [N,w,b] = eda_net_1dtower( Nc, slope, Xc, Dx, h, N, w, b )
% eda_net_linear
% defines a network that computes a set of 1D towers
% input Nc: number of towers
%       slope: slope of all the towers
%       Xc: centers of the towers
%       Dx: width of the towers
%       h: heights of the towaers
%       N: layer array from eda_net_init()
%       w: weights array from eda_net_init()
%       b: bias array from from eda_net_init()
% output: updated N, w, b arrays

N = [1, Nc*2, 1]';

for ncvalue = 1:Nc
    w12 = 4*slope;
    w22 = 4*slope;
    w( ((ncvalue-1)*2)+1:((ncvalue-1)*2)+2, 1:N(1),2) = [w12;w22];

    w13 = h(ncvalue);
    w23 = -h(ncvalue);
    w(1:N(3), ((ncvalue-1)*2)+1:((ncvalue-1)*2)+2  ,3) = [w13;w23];

    b12 = -w12 * ( Xc(ncvalue) - Dx(ncvalue)/2);
    b22 = -w22 * ( Xc(ncvalue) + Dx(ncvalue)/2);
    b( ((ncvalue-1)*2)+1:((ncvalue-1)*2)+2 ,2) = [b12;b22];
end

end