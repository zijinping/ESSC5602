function [N,w,b] = eda_net_linear( Nc, slope1, slope2, c, h, N, w, b )
% eda_net_linear
% defines a network that cpmputes the linear function
% d(x1, x2, ...) = c + h1*x1 + h2*x2 + h2*h3 ...
% (good for implementing a filter)
% input Nc: number of linear components
%       slope1: slope of the first stage (say 0.01)
%       slope2: slope of the second stage (say 0.02)
%               should be different from slope1
%       c: intercept
%       h: slopes
%       N: layer array from eda_net_init()
%       w: weights array from eda_net_init()
%       b: bias array from from eda_net_init()
% output: updated N, w, b arrays

N = [Nc, Nc, Nc, 1]';

b(1,4) = c;

for ncvalue = 1:Nc
% these two slopes are arbitrary, as long as they are small


    w12 = 4*slope1;
    w(ncvalue, ncvalue ,2) = w12;
  
    w13 = 4*slope2;
    w(ncvalue, ncvalue ,3) = w13;
    
    w2w3 = w12*w13;
    
    w14 = 16*h(ncvalue)/w2w3;
    w(1,ncvalue, 4) = w14;

    b12 = 0;
    b(ncvalue,2) = b12;
    
    b13 = -w13/2;
    b(ncvalue,3) = b13;
    
    b14 = -8*h(ncvalue)/w2w3;
    b(1,4) = b(1,4)+ b14;
end 

end