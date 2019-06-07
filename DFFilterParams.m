function [Nb,b,M,h_bp,H_z,negNorm] = DFFilterParams(d,theta_c,theta_zmin,theta_zmax)
    % John Vorhies, The University of Akron, Feb 2019
    % Determines the filter parameters for the dual-fan filter bank
    % Dansereau, et. al 2007
    % Passes some parameters needed for DFFVisuals and DFIIR
    % Input: Minimum and maximum depths relative to focal length d
    % Output:
    %       Nb:      Number of subbands to divide EPI
    %       b:       Feedback filter coefficients
    %       M:       Bandpass filter length
    %       h_bp:    Nb row vectors of FIR bandpass coefficients
    %       negNorm: Boolean indicator for negative norms
    
    %Subband filter parameters
    Nb = 8; %number of subbands
    L = 2*(Nb-1);
    q = 1; %number of side lobes for sinc function
    M = 2*L*q+1;

    negNorm = zeros(1,2);

    u = linspace(0,M-1,M);
    h_lb = 1/L*sinc((u-(M-1)/2)/L);
    h_bp = zeros(Nb,M);

    for nb = 1:Nb
        if nb == 1
            h_bp(nb,:) = h_lb(1:M);
        elseif nb == Nb
            h_bp(nb,:) = h_lb(1:M).*cos(2*pi*(nb-1)*u/L);
        else
            h_bp(nb,:) = 2*h_lb(1:M).*cos(2*pi*(nb-1)*u/L);
        end
    end

    %Hyperplanar filter parameters
    zmin = d/(tan(theta_zmin)+1);
    zmax = d/(tan(theta_zmax)+1);
    N = [1 ,-tan(theta_c)]/sqrt(1+(tan(theta_c))^2);
    B = zeros(1,Nb);
    b = zeros(2,2,Nb);

    if N(1) < 0
        negNorm(1) = true;
        N(1) = abs(N(1));
    else
        negNorm(1) = false;
    end

    if N(2) < 0
        negNorm(2) = true;
        N(2) = abs(N(2));
    else
        negNorm(2) = false;
    end

    %Calculate bandwidths
    for nb = 1:Nb
        B(nb) = abs(d*(nb-1+0.5)/(2*Nb)*(1/zmin-1/zmax));
    end
    %Calculate weights
    for nb=1:Nb
        for j = 1:2
            for k = 1:2
                b(j,k,nb) = 1 + ((-1)^(j-1)*N(2) + (-1)^(k-1)*N(1))/B(nb);
            end
        end
    end
    
    %Calculate frequency-domain filter
    omega = linspace(-pi,pi,1024);
    z_s = exp(1j*omega);
    z_u = exp(1j*omega);
    H_z = zeros(length(omega),length(omega));

    for nb = 1:Nb
        for i = 1:length(omega)
            for j = 1:length(omega)
                H_z(i,j,nb) = (1+z_u(j).^-1+z_s(i).^-1+z_u(j).^-1*...
                    z_s(i).^-1)./(b(1,1,nb)+b(1,2,nb)*z_u(j).^-1+b(2,1,nb)*...
                    z_s(i).^-1+b(2,2,nb)*z_s(i).^-1*z_u(j).^-1);
            end
        end
    end
    
    %Continuous frequency response
    omega = linspace(-pi,pi,1024);
    s_su = 1j*omega;
    H_s = zeros(length(omega),length(omega));
    for nu = 1:1024
        for ns = 1:1024
            H_s(ns,nu) = (1+(N(1)*s_su(ns)+N(2)*s_su(nu))/B(4)).^(-1);
        end
    end
    figure
    mesh(omega,omega,abs(H_s))
    xlabel('u')
    ylabel('s')
    view([90 270])
    title('continuous freq response')
end