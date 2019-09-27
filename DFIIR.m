function y2 = DFIIR(x,b,negNorm)
    % John Vorhies, The University of Akron, Feb 2019
    % Computes the output of the 4D separable dual-fan IIR filter,
    % Dansereau, et. al 2007
    % Input:
    %       x:       An epipolar image from an (s,t,u,v) light field
    %                parameterization
    %       b:       The feedback filter coefficients (from DFFilterParams)
    %       negNorm: indication from DFFilterParams that the theta_c
    %                value has produced a negative normal
    % Output:
    %       y2:      filtered EPI
    
    % check for negative normals
    if negNorm(1) == true
        x = flip(x,1);
    end
    if negNorm(2) == true
        x = flip(x,2);
    end 
    % pad x with zeroes for zero-phase filtering
    s_pad = 4;
    u_pad = 100;
    x = padarray(x,[s_pad,u_pad],0,'post');
    [Ns,Nu] = size(x);
    y1 = zeros(Ns,Nu);
    y2 = y1;

    for ns = 2:Ns
        for nu = 2:Nu
            x_comp = x(ns,nu) + x(ns,nu-1) + x(ns-1,nu) + x(ns-1,nu-1);
            y_comp = b(1,2)*y1(ns,nu-1)+b(2,1)*y1(ns-1,nu)+b(2,2)*y1(ns-1,nu-1);
            y1(ns,nu) = 1/b(1,1)*(x_comp-y_comp);
        end
    end

    % Zero-phase filtering
    y1 = flip(y1,1);
    y1 = flip(y1,2);
    for ns = 2:Ns
        for nu = 2:Nu
            x_comp = y1(ns,nu) + y1(ns,nu-1) + y1(ns-1,nu) + y1(ns-1,nu-1);
            y_comp = b(1,2)*y2(ns,nu-1)+b(2,1)*y2(ns-1,nu)+b(2,2)*y2(ns-1,nu-1);
            y2(ns,nu) = 1/b(1,1)*(x_comp-y_comp);
        end
    end
    y2 = flip(y2,2);
    y2 = flip(y2,1);

    y2(end+1-s_pad:end,:) = [];
    y2(:,end+1-u_pad:end) = [];

    if negNorm(1) == true
        y2 = flip(y2,1);
    end
    if negNorm(2) == true
        y2 = flip(y2,2);
    end
end 









