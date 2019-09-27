function filtered_image = fastDualFanFilterUV(st_uv,d)
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

    [Nt,Ns,Nv,Nu] = size(st_uv);
    st_center = ceil(Nt/2);
    v_center = ceil(Nv/2);

    st_uv = double(normalizeLF(st_uv));
    
    EPI = squeeze(st_uv(st_center,:,v_center,:));
    
    [theta_c,theta_zmin,theta_zmax] = findThetaC(EPI,1024);
    
    angle = ceil(length(theta_c)/2);
    [Nb,b,M,h_bp,negNorm] = DFFilterParams(d,theta_c(2),theta_zmin(2),theta_zmax(2));

    % Visualization of frequency responses
    %DFFVisuals(st_uv,h_bp,b);

    %-------------------------Filtering along u,v-------------------------
    % T,V subband and hyperplanar filter
    delay = floor((M-1)/2);
    delay1 = zeros(Nt,delay);
    delay2 = zeros(Ns,delay);
    st_uv_filt = zeros(Nt,Nv,Nu,Nb);

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
    st_uv = normalizeLF(st_uv);

    st_uv_recon1 = []; %#ok<NASGU>
    
    % S,U subband and hyperplanar filter
    st_uv_filt = zeros(Ns,Nv,Nu,Nb);
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












