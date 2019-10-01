function [theta_c,theta_zmin,theta_zmax] = findThetaC(EPI, fftsize)
    % John Vorhies, The University of Akron, Feb 2019
    % Finds the angle of the central axis, zmin and zmax
    % for each depth in an EPI.
    % Input: 
    %       EPI:        Epipolar image from an (s,t,u,v) parameterization
    %       fftsize:    Size of fft to perform for analysis
    % Output:
    %       theta_c, theta_zmin, theta_zmax:   from farthest depth to 
    %                                          closest depth.

    %Convert EPI to frequency domain and normalize
    EPI_fft = fftshift(fft2(EPI,fftsize,fftsize));
    EPI_fft = abs(EPI_fft);
    EPI_fft = log(EPI_fft+10e-6);
    EPI_fft = EPI_fft/max(EPI_fft(:));
    EPI_fft = EPI_fft';
    [Ns, Nu] = size(EPI_fft);

    x = linspace(-pi,pi,fftsize); 
    y = linspace(-pi,pi,fftsize);
    [X,Y] = meshgrid(x, y);
    sigmaS = 0.1;
    sigmaU = 0.7;
    scaleS = 1;
    scaleU = 10;
    sigmaS = scaleS*sigmaS;
    sigmaU = scaleU*sigmaU;
    theta = linspace(-pi/2,pi/2,200);
    filt_norm = zeros(1,length(theta));

    %---------------Pre-filtering and thresholding--------------
    
    % Filter omega_s axis
    theta_prefilt = 0;
    sigmaU_prefilt = sigmaU * 10;
    sigmaS_prefilt = sigmaS;
    EPI_fft_prefilt = EPI_fft;
    a = ((cos(theta_prefilt)^2) / (2*sigmaU_prefilt^2)) +...
        ((sin(theta_prefilt)^2) / (2*sigmaS_prefilt^2));
    b = -((sin(2*theta_prefilt)) / (4*sigmaU_prefilt^2)) +...
        ((sin(2*theta_prefilt)) / (4*sigmaS_prefilt^2));
    c = ((sin(theta_prefilt)^2) / (2*sigmaU_prefilt^2)) +...
        ((cos(theta_prefilt)^2) / (2*sigmaS_prefilt^2));

    gauss = exp(-(a*X.^2 + 2*b*X.*Y + c*Y.^2));
    EPI_fft_prefilt = (1-gauss) .* EPI_fft_prefilt;

    EPI_fft_prefilt = EPI_fft_prefilt/max(EPI_fft_prefilt(:));
    EPI_fft_thresh = EPI_fft_prefilt.^2;
    EPI_mean = mean(EPI_fft_thresh(:));
    EPI_std = std(EPI_fft_thresh(:));

    for nu = 1:Nu
        for ns = 1:Ns
            if EPI_fft_thresh(ns,nu) < EPI_mean - EPI_std
                EPI_fft_thresh(ns,nu) = 0;
            elseif EPI_fft_thresh(ns,nu) > EPI_mean + 2*EPI_std
                EPI_fft_thresh(ns,nu) = 1;
            end 
        end
    end

    %Apply 2D Gaussian filter at rotating angles
    j = 1;
    for k = length(theta):-1:1
        a = ((cos(theta(k))^2) / (2*sigmaS^2)) +...
            ((sin(theta(k))^2) / (2*sigmaU^2));
        b = ((sin(2*theta(k))) / (4*sigmaS^2)) -...
            ((sin(2*theta(k))) / (4*sigmaU^2));
        c = ((sin(theta(k))^2) / (2*sigmaS^2)) +...
            ((cos(theta(k))^2) / (2*sigmaU^2));

        gauss = exp(-(a*X.^2 + 2*b*X.*Y + c*Y.^2));
        filt_norm(j) = norm(gauss .* EPI_fft_thresh,'fro');
        j = j+1;
    end
    
    filt_norm = filt_norm/max(filt_norm);
    
    %Find peaks and valleys
    [~,peaks] = findpeaks(filt_norm,'MinPeakProminence',0.01);
    [~,valleys] = findpeaks(1 - filt_norm,'MinPeakProminence',0.01);
    %[~,peaks] = findpeaks(filt_norm);
    %[~,valleys] = findpeaks(1-filt_norm);
    
    %Find cutoff to determine theta_zmin, theta_zmax
    theta_zmin = zeros(1,length(peaks));
    theta_zmax = zeros(1,length(peaks));
    B = 0.9; % Cutoff bandwidth
    
    if isempty(valleys) ~= true % Check for only one global max
        for j = 1:length(peaks)
            %first iteration
            if j == 1
                for k = 1:peaks(j)
                    if filt_norm(k) > B * filt_norm(peaks(j))
                        theta_zmax(j) = theta(k);
                        break
                    end
                end
                for k = peaks(j):valleys(j)
                    if filt_norm(k) < B * filt_norm(peaks(j))
                        theta_zmin(j) = theta(k);
                        break
                    end
                    if theta_zmin(j) == 0 && k == valleys(j)
                        theta_zmin(j) = theta(valleys(j));
                    end
                end
                continue
            end
            %last iteration
            if j == length(peaks)
                for k = valleys(j-1):peaks(j)
                    if filt_norm(k) > B * filt_norm(peaks(j))
                        theta_zmax(j) = theta(k);
                        break
                    end
                end
                for k = peaks(j):length(theta)
                    if filt_norm(k) < B * filt_norm(peaks(j))
                        theta_zmin(j) = theta(k);
                        break
                    end
                    if theta_zmin(j) == 0 && k == length(theta)
                        theta_zmin(j) = theta(length(theta));
                    end
                end
                continue
            end
            %middle iterations
            for k = valleys(j-1):peaks(j)
                if filt_norm(k) > B * filt_norm(peaks(j))
                    theta_zmax(j) = theta(k);
                    break
                end
            end
            for k = peaks(j):valleys(j)
                if filt_norm(k) < B * filt_norm(peaks(j))
                    theta_zmin(j) = theta(k);
                    break
                end
                if theta_zmin(j) == 0 && k == valleys(j)
                    theta_zmin(j) = theta(valleys(j));
                end
            end

        end
    else
        for k = 1:peaks
            if filt_norm(k) > B * filt_norm(peaks)
                theta_zmax = theta(k);
                break
            end
        end
        for k = peaks:length(theta)
            if filt_norm(k) < B * filt_norm(peaks)
                theta_zmin = theta(k);
                break
            end
        end
    end
        
    theta_c = zeros(1,length(peaks));

    for k = 1:length(theta_c)
        theta_c(k) = theta(peaks(k));
    end
    
    figure
    mesh(EPI_fft_thresh)
    xlabel('$$\omega_{s}$$')
    ylabel('$$\omega_{u}$$')
    title('Thresholded EPI')
    view([0 90])
    
    figure
    plot(theta,filt_norm)
    xlabel('$$\theta$$')
    ylabel('Norm')
    title('Theta Norms')
end


