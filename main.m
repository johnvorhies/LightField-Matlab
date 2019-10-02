% John Vorhies, The University of Akron, Sept 2019

% input file for depth filtering in light fields. Assumes a 4-D gray scale 
% light field input named "st_uv" or 5-D RGB light field input named 
% "st_uv_rgb" with d defined as the camera focal length. The output is the 
% depth filtered central image of the light field.

%---------------------- RGB -----------------------------
tic
filtered_image_rgb = fastDualFanFilterUV(st_uv_rgb,d);
toc

figure
imshow(filtered_image_rgb)

%%
%--------------------- Grayscale ------------------------
tic
filtered_image = fastDualFanFilterUV(st_uv,d);
toc

figure
imshow(filtered_image)

%%
%---------------------- Wavelet Compression ---------------------
% Work in progress
tic
[st_uv_wavelet,wavelet_details] = LFWaveletCompression(st_uv);
cH = squeeze(wavelet_details(:,:,1));
cV = squeeze(wavelet_details(:,:,2));
cD = squeeze(wavelet_details(:,:,3));

filtered_image = fastDualFanFilterUV(st_uv_wavelet,d);
filtered_image = uint16(65535*mat2gray(idwt2(filtered_image,cH,cV,cD,'db20')));
toc

figure
imshow(filtered_image)



