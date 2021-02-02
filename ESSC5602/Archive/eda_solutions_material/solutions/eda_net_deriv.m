function [daLmaxdw,daLmaxdb] = eda_net_deriv( N,w,b,a )
% eda_net_deriv
% calculates the derivatives daLmax/dw and daLmax/dw
% needed to train the network
    Nmax = max(N);
    Lmax = length(N);
    %   Detailed explanation goes here
    % dada(i,j,k): the change in the activity of the i-th neuron of the layer k
    % due to a change in the activity of the j-th neuron in the layer k-1
    dada = zeros(Nmax,Nmax,Lmax);

    % daLmaxda(i,j,k): the change in the activity of the i-th neuron in layer Lmax
    % due to a change in the activity of the j-th neuron in the layer k-1
    daLmaxda = zeros(Nmax,Nmax,Lmax);

    % dadz(i,k): the change in the activity of the i-th neuron of layer k
    % due to a change in the input z of the same neuron
    dadz = zeros(Nmax,Lmax);
    
    % daLmaxdb(i,j,k): the change in the activity of the i-th neuron in layer Lmax
    % due to a change in the bias of the of the j-th neuron in the layer k
    daLmaxdb = zeros(Nmax,Nmax,Lmax);

    % daLmaxdw(i,j,k,l): the change in the activity of the i-th neuron in layer Lmax
    % due to a change in the weight connecting the j-th neuron in the layer l
    % with the k-th neuron of layer l-1
    daLmaxdw = zeros(Nmax,Nmax,Nmax,Lmax);

    % da/dz
    for i = [1:Lmax]
        if (i==Lmax)
            % sigma function is sigma(z)=z for top layer
            % d sigma / dz = 1
            dadz(1:N(i),i) = ones(1:N(i),1);
        else
            % sigma function is sigma(z)=(1+exp(-z))^(-1) for other layers
            % d sigma / dz = (-1)(-1) exp(-z)(1+exp(-z))^(-1)
            %              = exp(-z) (a^2)
            % and note exp(-z) = (1/a)-1
            % so d sigma / dz = a-a^2
            dadz(1:N(i),i) = a(1:N(i),i) - (a(1:N(i),i).^2);  
        end
    end
        
    % da/da
    for i = [Lmax:-1:2]
        if i == Lmax
            % a=z in the top of the layer, and
            % z = (sum) w a + b
            % so da(layer i)/da(layer i-1) = dz/da(layer i-1) = w
            dada(1:N(i),1:N(i-1),i) = w(1:N(i),1:N(i-1),i);
            % initialize daLmaxda to prepare for backpropagation
            daLmaxda(1:N(i),1:N(i-1),i) = dada(1:N(i),1:N(i-1),i);
        else
            % a=sigma(z) in the other layers, so use the chain rule
            % da/da = da(layer)/dz * dz/da(layer i-1)
            %    where as before dz/da(layer i-1) = w
            dada(1:N(i),1:N(i-1),i) = diag(dadz(1:N(i),i))*w(1:N(i),1:N(i-1),i);
            % backpropagate
            daLmaxda(1:N(Lmax),1:N(i-1),i) =  ...
                daLmaxda(1:N(Lmax),1:N(i),i+1)*dada(1:N(i),1:N(i-1),i);
        end            
    end
        
    % da/db
    for i = [1:Lmax]
        if( i==Lmax )
            % a=z in the top layer, and z(layer i)= w *a(layer i-1) + b
            % so da(layer i)/db = dz(layer i)/db = 1
            daLmaxdb(1:N(Lmax),1:N(Lmax),Lmax) = eye(N(Lmax),N(Lmax));
        else
            % a = sigma(z) in the other layers, so use the chain rule
            % da(layer i)/db = da(layer i)/dz * dz/db
            % where as before dzdb=1
            dzdb = 1;
            dadb = dadz(1:N(i),i) * dzdb;
            % then apply the chain rule for activities
            % da(layer Lmax)/db = da(layer Lmax)/da(layer i) * da(layer i)/db
            daLmaxdb(1:N(Lmax),1:N(i),i) = daLmaxda(1:N(Lmax),1:N(i),i+1) * diag(dadb);
        end
    end

    % da/dw
    % calculation of weight derivatives
    % since z(layer i) = w*a(layer i-1) + b
    % dz(layer i)/dw = a(layer i-1)
    for l = [2:Lmax]  % for simplicity, loop over right neuron
        if( l==Lmax )
            % in the top layer, a=z, so da/dw = dz/dw = a
            for j = [1:N(l)]
                daLmaxdw(j,j,1:N(l-1),l) = a(1:N(l-1),l-1);
            end
        else
            % in the other layers, use chain rule
            % da/dw = da/dz * dz/dw
            % and then the chain ruke again
            % da(layer Lmax)/dw = da(layer Lmax)/da(layer i) * da(layer i)/dw
            for j = [1:N(l)]
                dzdw = a(1:N(l-1),l-1);
                dadw = dadz(j,l) * dzdw;
                daLmaxdw(1:N(Lmax),j,1:N(l-1),l) = daLmaxda(1:N(Lmax),j,l+1) * dadw;
            end
        end
    end

end

