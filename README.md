# LightField-Matlab
LightField-Matlab is a Matlab library for light field image processing. It implements the 4-D IIR dual-fan filter bank for depth filtering (Dansereau, et. al 2007), and an algorithm for finding minimum and maximum depths from a light field image.

## Getting Started
Depth filtering can be accomplished from **Main**. Load the .mat file for the light field, add the file as an input argument for **fastDualFanFilterUV** and run.

All of the graphs used Latex markup. Add these 3 lines to your Matlab startup file:
```Matlab
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
```

### Reading Material

[A 4-D Dual-Fan Filter Bank for Depth Filtering in Light Fields](https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=4063539)

### Datasets

[JPEG Pleno Database](https://jpeg.org/plenodb/)

## Functions

### Main
Calls fastDualFanFilterUV for either grayscale or RGB light field depth filtering. It requires a .mat file of the light field to be loaded into the Matlab workspace as a 4-D grayscale or 5-D RGB tensor.

### filtered_image = fastDualFanFilterUV(st_uv)
Main control function for depth filtering. Calls **findThetaC**, **DFFilterParams**, **plotFrequencyResponse** and **applyFilter**. This function also checks whether the input light field is grayscale or RGB. If the light field is grayscale, an epipolar image (EPI) is extracted from the central light field image for **findThetaC** and **applyFilter** is called once. If the light field is RGB, an EPI is extracted from the central light field image and converted to grayscale for **findThetaC** and **applyFilter** is called three times, once for each color channel. The choice of *theta_c* for **DFFilterParams** from **findThetaC** is set to the middle element of the returned *theta_c array*. More information is available in the documentation for the **findThetaC** function.

#### Inputs:
* st_uv: a light field in two-plane parameterization form.

#### Outputs:
* filtered_image: The depth-filtered central image of the input light field.


### [theta_c,theta_zmin,theta_zmax] = findThetaC(EPI,fftsize)
Searches the Fourier domain of an EPI for energy content that resembles the dual-fan shape. This function returns three vectors: *theta_c*, *theta_zmin*, and *theta_zmax*. These values are sorted from farthest depth to closest depth. At the end of the function is a section for plotting the thresholded Fourier domain that is being searched and the norm for each theta value.

#### Inputs:
* EPI: An epipolar image to search for depth information
* fftsize: Size of 2-D Fast Fourier Transform (FFT) to perform on the EPI. This should be atleast the size of the longest dimension of the EPI.

#### Outputs:
* theta_c: the angle of the central axis of the depth content to the omega_u axis.
* theta_zmin: the angle of the right-most edge of the depth content to the omega_u axis.
* theta_zmax: the angle of the left-most edge of the depth content to the omega_u axis.

### [Nb,b,M,h_bp,negNorm,N,B] = DFFilterParams(theta_c,theta_zmin,theta_zmax)
Creates the filter coefficients for the dual-fan filter and the coefficients for the sub-band filters in the dual-fan filter bank. This function also checks if the orientation of the dual-fan filter involves a negative normal value for use by **applyFilter**.

#### Inputs:
* theta_c, theta_zmin, theta_zmax: a value selected from the output of *findThetaC*.

#### Outputs:
* Nb: Number of sub-bands for the filter bank. This is an internal parameter of **DFFilterParams**.
* b: The filter coefficients for convolution used in **DFIIR**.
* M: The length of the sub-band filters used in **applyFilter**.
* h_bp: The filter coefficients for the sub-band filters used in **applyFilter**.
* negNorm: Boolean value to indicate whether a negative normal was found
* N: The normal vectors used to determine the orientation of the dual-fan filter.
* B: The bandwidth used for each sub-band in the filter bank.

#### Usage:
To create the filter parameters corresponding to the farthest depth found by **findThetaC**,
```Matlab
[Nb,b,M,h_bp,negNorm,N,B] = DFFilterParams(theta_c(1),theta_zmin(1),theta_zmax(1))
```
Or to create filter parameters corresponding to the closest depth found by **findThetaC**,
```Matlab
[Nb,b,M,h_bp,negNorm,N,B] = DFFilterParams(theta_c(length(theta_c)),theta_zmin(length(theta_zmin)),theta_zmax(length(theta_zmax)))
```

### plotFrequencyResponse(st_uv,h_bp,b,N,B)
Plots the frequency responses for the dual-fan filter bank. It also provides a single plot of the center sub-band for the continuous-time filter. This can be compared to the frequency response of the discrete-time implementation to view the frequency warping effects caused by the discretization process.

#### Inputs:
* st_uv: a light field in two-plane parameterization form.
* h_bp: The filter coefficients for the sub-band filters from **DFFilterParams**.
* b: The filter coefficients for convolution used in **DFIIR** from **DFFilterParams**.
* N: The normal vectors used to determine the orientation of the dual-fan filter from **DFFilterParams**.
* B: The bandwidth used for each sub-band in the filter bank from **DFFilterParams**.

### filtered_image = applyFilter(st_uv,h_bp,negNorm,Nb,M)
Implementation of the filter bank. Extracts each EPI from the central light field image, divides them into Nb sub-bands and applys the frequency-planar filter at varying bandwidths.

#### Inputs:
* st_uv: a light field in two-plane parameterization form.
* h_bp: The filter coefficients for the sub-band filters from **DFFilterParams**.
* negNorm: Boolean value to indicate whether a negative normal was found
* Nb: Number of sub-bands for the filter bank. This is an internal parameter of **DFFilterParams**.
* M: The length of the sub-band filters from **DFFilterParams**.

#### Outputs:
* filtered_image: The depth-filtered central image of the light field.

### y2 = DFIIR(x,b,negNorm)
Implementation of the frequency-planar filter. Takes an input signal *x* (EPI) and filter coefficients *b* to perform convolution via direct form I. To ensure practical-BIBO stability of the filter, the direction of iteration should be reversed for dimensions with a negative normal. Because an IIR filter is used, the filter output does not have linear phase. This can cause artifacts in the depth-filtered image. To account for this, the filter is made zero-phase by applying the convolution again with the direction of iteration reversed for both dimensions. The input signal is also padded in both directions to account for the transient response of the filter. Initial conditions for the filter are assumed to be zero.

#### Inputs:
* x: A band-passed version of the EPI from **applyFilter** to be filtered
* b: The filter coefficients for convolution from **DFFilterParams**
* negNorm: Boolean indicating if the direction of iteration should be reversed

#### Outputs:
* y2: The filtered sub-band signal to be summed with other filtered sub-band signal components by **applyFilter**.

### st_uv = normalizeLF(st_uv)
Takes each image in the light field and normalizes them for 16-bit grayscale or RGB. If the output image of the filter bank has a loss in dynamic range, this should be uncommented in **applyFilter**.
