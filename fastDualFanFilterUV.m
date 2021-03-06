function filtered_image = fastDualFanFilterUV(st_uv)
    % John Vorhies, The University of Akron, Feb 2019
    % Implements the "fast" dual-fan filter bank, Dansereau, et. al 2007.
    % Only the outputs needed to filter the central light field image are 
    % computed.
    % INPUT:
    %       st_uv:  light field in (s,t,u,v) parameterization, 16-bit
    %               grayscale.
    %       d:      distance between (s,t) and (u,v) planes (focal length)
    % Output:
    %       filtered_image: depth-filtered central light field image
    %
    % ASSUMES SAME SIZE FOR S AND T DIMENSIONS

    [Nt,Ns,Nv,Nu,Nc] = size(st_uv);
    fftsize = max(size(st_uv));
    st_center = ceil(Nt/2);
    u_center = ceil(Nu/2);
    v_center = ceil(Nv/2);
    
    st_uv = single(st_uv);
    
    if Nc == 1
        EPI = squeeze(st_uv(st_center,:,v_center,:));

        [theta_c,theta_zmin,theta_zmax] = findThetaC(EPI,fftsize);
        angle = ceil(length(theta_c)/2);
        [Nb,b,M,h_bp,negNorm,N,B] = DFFilterParams(theta_c(angle),theta_zmin(angle),theta_zmax(angle));
    else
        EPI = squeeze(st_uv(st_center,:,3124,:,:));
        EPI = rgb2gray(uint16(EPI));
        
        [theta_c,theta_zmin,theta_zmax] = findThetaC(EPI,fftsize);
        angle = ceil(length(theta_c)/2);
        [Nb,b,M,h_bp,negNorm,N,B] = DFFilterParams(theta_c(angle),theta_zmin(angle),theta_zmax(angle));
    end
        
    % Visualization of frequency responses
    plotFrequencyResponse(st_uv,h_bp,b,N,B);
    
    filtered_image = zeros(Nv,Nu,Nc,'uint16');
    if Nc == 1
        filtered_image = applyFilter(st_uv,b,h_bp,negNorm,Nb,M);
    else
        for Nc = 1:3
            st_uv_channel = squeeze(st_uv(:,:,:,:,Nc));
            filtered_image(:,:,Nc) = applyFilter(st_uv_channel,b,h_bp,negNorm,Nb,M);
        end
    end
   
end












