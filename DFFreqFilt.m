function x = DFFreqFilt(x,H_z,negNorm)
    % John Vorhies, The University of Akron, May 2019
    % Frequency domain implementation of the dual-fan depth filter,
    % Dansereau, et. al 2007
    % Input:
    %       x:       An epipolar image from an (s,t,u,v) light field
    %                parameterization
    %       H_z:     Frequency-domain filter response (DFFilterParams)
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
    
    [Ns,Nu] = size(x);
    
    % Frequency Filtering
    x = fftshift(fft2(x,1024,1024));
    x = x .* H_z;
    x = ifft2(fftshift(x),Ns,Nu,'symmetric');
    
    % Zero-Phase Filtering
    x = flip(x,1);
    x = flip(x,2);
    x = fftshift(fft2(x,1024,1024));
    x = x .* H_z;
    x = ifft2(fftshift(x),Ns,Nu,'symmetric');
    x = flip(x,1);
    x = flip(x,2);
    
    if negNorm(1) == true
        x = flip(x,1);
    end
    if negNorm(2) == true
        x = flip(x,2);
    end

end