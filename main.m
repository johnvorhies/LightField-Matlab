% John Vorhies, The University of Akron, Sept 2019

% input file for depth filtering in light fields. Assumes a 4-D grayscale
% or 5-D RGB light field input. The output is the depth filtered central
% image of the light field.

tic
filtered_image = fastDualFanFilterUV(st_uv_rgb);
toc

figure
imshow(filtered_image)



