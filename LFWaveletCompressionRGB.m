function [st_uv_waveletGray,wavelet_detailsGray] = LFWaveletCompressionRGB(st_uv_RGB)
    % John Vorhies, The University of Akron, Feb 2019
    % Performs 2D discrete wavelet transform and provides analysis
    % on light fields.
    % Input:
    %       st_uv_RGB: RGB light field in the (s,t,u,v) parameterization
    % Output:
    %       st_uv_waveletGray: wavelet transformed light field in
    %                          8-bit grayscale
    
    [Nt,Ns,Nv,~,Nc] = size(st_uv_RGB);
    st_center = ceil(Nt/2);
    v_center = ceil(Nv/2);

    st_uv_wavelet = zeros(Nt,Ns,531,531,3,'uint8');
    wavelet_details = zeros(531,531,3,3,'uint8');
    v_center_wavelet = ceil(531/2);

    for nc = 1:Nc
        for nt = 1:Nt
            for ns = 1:Ns
                [cA,cH,cV,cD] = dwt2(squeeze(st_uv_RGB(nt,ns,:,:,nc)),'db20');
                if nt == st_center && ns == st_center
                    wavelet_details(:,:,nc,1) = cH;
                    wavelet_details(:,:,nc,2) = cV;
                    wavelet_details(:,:,nc,3) = cD;
                end
                st_uv_wavelet(nt,ns,:,:,nc) = uint8(255*mat2gray(cA));
            end
        end
    end
    
    st_uv_waveletGray = zeros(nt,ns,531,531,'uint8');
    for nt = 1:Nt
        for ns = 1:Ns
            st_uv_waveletGray(nt,ns,:,:,:) = rgb2gray(squeeze(st_uv_wavelet(nt,ns,:,:,:)));
        end
    end
    
    wavelet_detailsGray = zeros(531,531,3,'uint8');
    for k = 1:3
        wavelet_detailsGray(:,:,k) = rgb2gray(squeeze(wavelet_details(:,:,:,k)));
    end
    
    
    % Analyze EPI's and EPI spectra
    waveletEPI = squeeze(st_uv_wavelet(st_center,:,v_center_wavelet,:,:));
    EPI = squeeze(st_uv_RGB(st_center,:,v_center,:,:));
    waveletEPIGray = rgb2gray(waveletEPI);
    EPIGray = rgb2gray(EPI);

    fftsize = 1024;

    waveletEPIfft = fftshift(fft2(waveletEPIGray,fftsize,fftsize));
    waveletEPIfft = log(abs(waveletEPIfft)+10e-6);

    EPIfft = fftshift(fft2(EPIGray,fftsize,fftsize));
    EPIfft = log(abs(EPIfft)+10e-6);

    figure
    image(EPI)
    title('EPI')
    xlabel('U axis')
    ylabel('S axis')

    figure
    image(waveletEPI)
    title('Wavelet EPI')
    xlabel('U axis')
    ylabel('S axis')

    omega = linspace(-pi,pi,fftsize);

    figure
    mesh(omega,omega,EPIfft)
    title('EPI Spectra')
    xlabel('$$\omega_{u}$$')
    ylabel('$$\omega_{s}$$')
    view([90 270])
    xlim([-pi pi])
    ylim([-pi pi])

    figure
    mesh(omega,omega,waveletEPIfft)
    title('Wavelet EPI Spectra')
    xlabel('$$\omega_{u}$$')
    ylabel('$$\omega_{s}$$')
    view([90 270])
    xlim([-pi pi])
    ylim([-pi pi])

end










