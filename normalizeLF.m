function st_uv = normalizeLF(st_uv)
    % John Vorhies, The University of Akron, Feb 2019
    % Normalizes input light field on 8-bit grayscale range
    
    [Nt,Ns,~,~] = size(st_uv);

    for ns = 1:Ns
        for nt = 1:Nt
            st_uv(nt,ns,:,:) = 255*mat2gray(squeeze(st_uv(nt,ns,:,:)));
        end
    end
end