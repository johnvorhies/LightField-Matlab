function filtered_image = applyFilter(st_uv,b,h_bp,negNorm,Nb,M)
    % Depth filtering operations on EPIs
    
    [Nt,Ns,Nv,Nu] = size(st_uv);
    st_center = ceil(Nt/2);
    
    % Uncomment if the filter output has a loss in dynamic range
    %st_uv = double(normalizeLF(st_uv));
    
    %-------------------------Filtering along u,v-------------------------
    % T,V subband and hyperplanar filter
    delay = floor((M-1)/2);
    delay1 = zeros(Nt,delay);
    delay2 = zeros(Ns,delay);
    st_uv_filt = zeros(Nt,Nv,Nu,Nb,'single');

    for nu = 1:Nu
        % compensate for group delay of filter
        tv = [squeeze(st_uv(:,st_center,:,nu)) delay1];
        for nb = 1:Nb
            filtered = filter(h_bp(nb,:),1,double(tv),[],2);
            filtered(:,1:delay) = [];
            st_uv_filt(:,:,nu,nb) = DFIIR(filtered,b(:,:,nb),negNorm);
        end
    end

    % sum components
    st_uv_recon1 = zeros(Nt,Nv,Nu);
    for nu = 1:Nu
        for nb = 1:Nb
            st_uv_recon1(:,:,nu) = st_uv_recon1(:,:,nu)+st_uv_filt(:,:,nu,nb);
        end
    end

    st_uv(:,st_center,:,:) = st_uv_recon1;
    
    % Uncomment if the filter output has a loss in dynamic range
    %st_uv = normalizeLF(st_uv);

    st_uv_recon1 = []; %#ok<NASGU>
    
    % S,U subband and hyperplanar filter
    st_uv_filt = zeros(Ns,Nv,Nu,Nb,'single');
    for nv = 1:Nv
        % compensate for group delay of filter
        su = [squeeze(st_uv(st_center,:,nv,:)) delay2];
        for nb = 1:Nb
            filtered = filter(h_bp(nb,:),1,double(su),[],2);
            filtered(:,1:delay) = [];
            st_uv_filt(:,nv,:,nb) = DFIIR(filtered,b(:,:,nb),negNorm);
        end
    end

    % sum components
    st_uv_recon2 = zeros(Ns,Nv,Nu);
    for nv = 1:Nv
        for nb = 1:Nb
            st_uv_recon2(:,nv,:) = st_uv_recon2(:,nv,:)+st_uv_filt(:,nv,:,nb);
        end
    end

    st_uv(st_center,:,:,:) = st_uv_recon2;

    filtered_image = uint16(65535*mat2gray(squeeze(st_uv(st_center,st_center,:,:))));
end

