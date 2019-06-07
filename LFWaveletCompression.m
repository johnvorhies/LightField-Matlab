function [st_uv_wavelet,wavelet_details] = LFWaveletCompression(st_uv)
    % John Vorhies, The University of Akron, Feb 2019
    % Performs 2D discrete wavelet transform. for performance over
    % analysis.
    % Input:
    %       st_uv: 8-bit Grayscale light field in the (s,t,u,v) 
    %       parameterization.
    % Output:
    %       st_uv_wavelet:    Wavelet transformed light field in
    %       8-bit grayscale.
    %       wavelet_details:  Extra details from DWT
    
    [Nt,Ns,~,~] = size(st_uv);
    st_center = ceil(Nt/2);

    st_uv_wavelet = zeros(Nt,Ns,531,531,'uint8');
    wavelet_details = zeros(531,531,3,'uint8');
    
    for nt = 1:Nt
        for ns = 1:Ns
            [cA,cH,cV,cD] = dwt2(squeeze(st_uv(nt,ns,:,:)),'db20');
            if nt == st_center && ns == st_center
                wavelet_details(:,:,1) = cH;
                wavelet_details(:,:,2) = cV;
                wavelet_details(:,:,3) = cD;
            end
            st_uv_wavelet(nt,ns,:,:) = uint8(255*mat2gray(cA));
        end
    end
 
end










