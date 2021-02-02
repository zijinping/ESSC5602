function [N,w,b] = eda_net_2dtower( Nc, slope, Xc, Yc, Dx, Dy, h, N, w, b )
% eda_net_linear
% defines a network that superimposes several towers

% input:
%       Nc: number of towers
%       slope: slope of all the towers
%       Xc, Dc: vector of the centers of the towers in x and y
%       Dx, Dy: vector of the widths of the towers in x and y
%       h: vector of heights of the towers 
%       N: layer array from eda_net_init()
%       w: weights array from eda_net_init()
%       b: bias array from from eda_net_init()
% output: updated N, w, b arrays

N = [2, Nc*4, Nc, 1]';

for ncvalue = 1:Nc

    w12 = 4*slope;
    w22 = 4*slope;
    w32 = 4*slope;
    w42 = 4*slope;
    w( ((ncvalue-1)*4)+1:((ncvalue-1)*4)+2, 1 ,2) = [w12 ; w22];
    w( ((ncvalue-1)*4)+3:((ncvalue-1)*4)+4, 2 ,2) = [w32 ; w42];
  
    w13 = 4*slope;
    w23 = -4*slope;
    w33 = 4*slope;
    w43 = -4*slope;
    w(ncvalue, ((ncvalue-1)*4)+1:((ncvalue-1)*4)+4,3) = [w13;w23;w33;w43];
    
    w14 = h(ncvalue);
    w(1, ncvalue, 4) = w14;

    b12 = -w12 * ( Xc(ncvalue) - Dx(ncvalue)/2);
    b22 = -w22 * ( Xc(ncvalue) + Dx(ncvalue)/2);
    b32 = -w32 * ( Yc(ncvalue) - Dy(ncvalue)/2);
    b42 = -w42 * ( Yc(ncvalue) + Dy(ncvalue)/2);
    b( ((ncvalue-1)*4)+1:((ncvalue-1)*4)+4 ,2) = [b12;b22;b32;b42];
    
    b13 = -(3*(4*slope))/2;
    b( ncvalue,3) = b13;
    
end 

end