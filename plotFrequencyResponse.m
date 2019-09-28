function plotFrequencyResponse(st_uv,h_bp,b,N,B)
    % John Vorhies, The University of Akron, Feb 2019
    % Shows the impulse and frequency responses of the FIR bandpass
    % filters, the 2D filter from the selected feedback coefficients
    % and the spectra of the input light field
    % Input parameters are outputs of DFFilterParams

    [Nt,~,Nv,~] = size(st_uv);
    st_center = ceil(Nt/2);
    v_center = ceil(Nv/2);
    fftsize = 1024;
    omega = linspace(-pi,pi,fftsize);
    
    [Nb, u] = size(h_bp);
    u = linspace(0,u-1,u);
    
    figure
    plot(u,h_bp(1,:))
    hold on
    for nb = 2:Nb
        plot(u,h_bp(nb,:))
    end
    hold off
    
    for nb = 1:Nb
        IR_legend{nb} = num2str(nb); %#ok<AGROW>
    end
    
    title('Subband Filters, Impulse Response')
    legend(IR_legend)

    h_bp_fft = zeros(Nb,fftsize);

    for nb = 1:Nb
        h_bp_fft(nb,:) = fft(h_bp(nb,:),fftsize);
    end

    figure
    plot(omega,fftshift(abs(h_bp_fft(1,:))))
    hold on
    
    for nb = 2:Nb
        plot(omega,fftshift(abs(h_bp_fft(nb,:))))
    end
    hold off
    
    for nb = 1:Nb
        FR_legend{nb} = num2str(nb);  %#ok<AGROW>
    end
    
    title('Subband Filters, Frequency Response')
    legend(FR_legend)

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

    for nb = 1:Nb
        figure
        mesh(omega,omega,abs(H_z(:,:,nb)))
        xlabel('$$\omega_{u}$$')
        ylabel('$$\omega_{s}$$')
        title_freq = strcat('Frequency Response, Band',{' '}, num2str(nb));
        title(title_freq)
        view([90 270])
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
    xlabel('$$\Omega_{u}$$')
    ylabel('$$\Omega_{s}$$')
    view([90 270])
    title('Continuous Frequency Response')

    EPI_fft = fftshift(fft2(squeeze(st_uv(st_center,:,v_center,:)),fftsize,fftsize));
    figure
    mesh(omega,omega,log(abs(EPI_fft+10e-6)))
    xlabel('$$\Omega_{u}$$')
    ylabel('$$\Omega_{s}$$')
    title('Light Field Spectra')
    view([90 270])
end