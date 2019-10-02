function [st_uv_wavelet,wavelet_details] = LFWaveletCompression(st_uv)
    % John Vorhies, The University of Akron, Feb 2019
    % Performs 2D discrete wavelet transform. for performance over
    % analysis.
    % Input:
    %       st_uv: 16-bit Grayscale light field in the (s,t,u,v) 
    %       parameterization.
    % Output:
    %       st_uv_wavelet:    Wavelet transformed light field in
    %       16-bit grayscale.
    %       wavelet_details:  Extra details from DWT
    
    [Nt,Ns,Nv,Nu] = size(st_uv);
    st_center = ceil(Nt/2);
    
    filter_length = 40; %db20
    Nv_wavelet_size = floor((Nv+filter_length-1)/2);
    Nu_wavelet_size = floor((Nu+filter_length-1)/2);

    st_uv_wavelet = zeros(Nt,Ns,Nv_wavelet_size,Nu_wavelet_size,'uint16');
    wavelet_details = zeros(Nv_wavelet_size,Nu_wavelet_size,3,'uint16');
    
    for nt = 1:Nt
        for ns = 1:Ns
            [cA,cH,cV,cD] = dwt2(squeeze(st_uv(nt,ns,:,:)),'db20');
            if nt == st_center && ns == st_center
                wavelet_details(:,:,1) = cH;
                wavelet_details(:,:,2) = cV;
                wavelet_details(:,:,3) = cD;
            end
            st_uv_wavelet(nt,ns,:,:) = uint16(65535*mat2gray(cA));
        end
    end
 
end










