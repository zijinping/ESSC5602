function [ a ] = eda_net_eval( N,w,b,a)
% eda_eval_net evaluates a network with properites N, w and b 

    Nmax = max(N);
    Lmax = length(N);

    for i=[2:Lmax]
        z(1:N(i),i) = w(1:N(i),1:N(i-1),i) * a(1:N(i-1),i-1) + b(1:N(i),i);
        if i < Lmax
            a(1:N(i),i) = 1.0./(1.0+exp(-z(1:N(i),i)));
        else
            a(1:N(i),i) = z(1:N(i),i);
        end
    end
end

