# LightField-Matlab
LightField-Matlab is a Matlab library for light field image processing. It implements the 4-D IIR dual-fan filter bank for depth filtering (Dansereau, et. al 2007), an algorithm for finding minimum and maximum depths from a light field image, and an algorithm for wavelet compression of light fields.

# Main
Calls fastDualFanFilterUV for either grayscale or RGB light field depth filtering and also calls LFWaveletCompression. It requires a .mat file of the light field to be loaded into the Matlab workspace with either the name "st_uv" or "st_uv_rgb" and the focal length of the camera "d".

# filtered_image = fastDualFanFilterUV(st_uv,d)
Main control function for depth filtering. Calls findThetaC, DFFilterParams, plotFrequencyResponse and applyFilter. This function also checks whether the input light field is grayscale or RGB. If the light field is grayscale, an EPI is extracted from the central light field image for findThetaC and applyFilter is called once. If the light field is RGB, an EPI is extracted from the central light field image and converted to grayscale for findThetaC and applyFilter is called three times, once for each color channel. For RGB, parfor is used to call applyFilter three times in parallel. This requires the Parallel Computing Toolbox. This can be changed to a for loop if this toolbox is not available to you.

# [theta_c,theta_zmin,theta_zmax] = findThetaC(EPI,fftsize)
Searches the Fourier domain of an EPI for energy content that resembles the dual-fan shape. This function returns three vectors: theta_c, theta_zmin, theta_zmax. These values are sorted from farthest depth to closest depth. At the end of the function is a section for plotting the thresholded Fourier domain that is being searched and the norm for each theta value.

# [Nb,b,M,h_bp,negNorm,N,B] = DFFilterParams(d,theta_c,theta_zmin,theta_zmax)
Creates the filter coefficients for the dual-fan filter and the coefficients for the band-pass filters in the dual-fan filter bank. This function also checks if the orientation of the dual-fan filter involves a negative normal value for use by applyFilter.

# plotFrequencyResponse(st_uv,h_bp,b,N,B)
Plots the frequency responses for the dual-fan filter bank. It also provides a single plot of the center sub-band for the continuous-time filter. This can be compared to the frequency response of the discrete-time implementation to view the frequency warping effects caused by the discretization process.

# filtered_image = applyFilter(st_uv,h_bp,b,N,B)
Implementation of the filter bank. Extracts each EPI from the central light field image, divides them into Nb sub-bands and applys the frequency-planar filter at varying bandwidths.

# y2 = DFIIR(x,b,negNorm)
Implementation of the frequency-planar filter. Takes an input x (EPI) and filter coefficients b to perform convolution via direct form I. To ensure practical-BIBO stability of the filter, the direction of iteration should be reversed for dimensions with a negative normal. Because an IIR filter is used, the filter output does not have linear phase. This can cause artifacts in the depth-filtered image. To account for this, the filter is made zero-phase by applying the convolution again with the direction of iteration reversed for both dimensions. The input signal is also padded in both directions to account for the transient response of the filter. Initial conditions for the filter are assumed to be zero.

# st_uv = normalizeLF(st_uv)
Takes each image in the light field and normalizes them for 16-bit grayscale. If the output image of the filter bank has a loss in dynamic range, this should be uncommented in applyFilter. Not recommended for RGB as it does not weight the individual channels.
