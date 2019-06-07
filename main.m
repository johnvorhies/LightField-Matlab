st_uv = load('st_uvTarotSmallRectified.mat');
st_uv = st_uv.st_uv;

% st_uv_RGB = load('st_uvTarotSmallRectifiedRGB.mat');
% st_uv_RGB = st_uv_RGB.st_uv_RGB;
% 
% [st_uv_wavelet,wavelet_details] = LFWaveletCompressionRGB(st_uv_RGB);
%%
% [st_uv_wavelet,wavelet_details] = LFWaveletCompression(st_uv);
% cH = squeeze(wavelet_details(:,:,1));
% cV = squeeze(wavelet_details(:,:,2));
% cD = squeeze(wavelet_details(:,:,3));

d = 22e-3;                  % focal length of camera
t = cputime;
filtered_image = fastDualFanFilterUV(st_uv,d);
%filtered_image = uint8(255*mat2gray(idwt2(filtered_image,cH,cV,cD,'db20')));
runtime = cputime - t;
figure
imshow(filtered_image)

%Test for commit
