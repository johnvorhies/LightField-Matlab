function [st_uv_wavelet_rgb,wavelet_details_rgb] = LFWaveletCompressionRGB(st_uv_rgb)
    % John Vorhies, The University of Akron, Feb 2019
    % Performs 2D discrete wavelet transform and provides analysis
    % on light fields.
    % Input:
    %       st_uv_RGB: RGB light field in the (s,t,u,v) parameterization
    % Output:
    %       st_uv_waveletGray: wavelet transformed light field in
    %                          16-bit grayscale
    
    [Nt,Ns,Nv,Nu,Nc] = size(st_uv_rgb);
    st_center = ceil(Nt/2);
    
    filter_length = 40; %db20
    Nv_wavelet_size = floor((Nv+filter_length-1)/2);
    Nu_wavelet_size = floor((Nu+filter_length-1)/2);

    st_uv_wavelet_rgb = zeros(Nt,Ns,Nv_wavelet_size,Nu_wavelet_size,3,'uint16');
    wavelet_details_rgb = zeros(Nv_wavelet_size,Nu_wavelet_size,3,3,'uint16');

    for nc = 1:Nc
        for nt = 1:Nt
            for ns = 1:Ns
                [cA,cH,cV,cD] = dwt2(squeeze(st_uv_rgb(nt,ns,:,:,nc)),'db20');
                if nt == st_center && ns == st_center
                    wavelet_details_rgb(:,:,nc,1) = cH;
                    wavelet_details_rgb(:,:,nc,2) = cV;
                    wavelet_details_rgb(:,:,nc,3) = cD;
                end
                st_uv_wavelet_rgb(nt,ns,:,:,nc) = uint16(65535*mat2gray(cA));
            end
        end
    end

end










