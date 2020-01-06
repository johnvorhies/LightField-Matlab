function st_uv = normalizeLF(st_uv)
    % John Vorhies, The University of Akron, Feb 2019
    % Normalizes input light field on 16-bit grayscale or RGB range
    
    [Nt,Ns,~,~,Nc] = size(st_uv);
    
    if Nc == 1
        for ns = 1:Ns
            for nt = 1:Nt
                st_uv(nt,ns,:,:) = 65535*mat2gray(squeeze(st_uv(nt,ns,:,:)));
            end
        end
    else
        for nc = 1:Nc
            for ns = 1:Ns
                for nt = 1:Nt
                    st_uv(nt,ns,:,:,nc) = 65535*mat2gray(squeeze(...
                        st_uv(nt,ns,:,:,nc)));
                end
            end
        end
    end
        
end