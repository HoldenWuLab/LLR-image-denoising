# Improved Liver Fat and R2* Quantification at 0.55T using Locally Low-Rank Denoising

This page shares the code for the paper published in Magnetic Resonance in Medicine under the title "Improved Liver Fat and R2* Quantification at 0.55T using Locally Low-Rank Denoising".

Please cite this work if you use the code:
Shih, S-F, Tasdelen, B, Yagiz, E, Zhang, Z, Zhong, X, Cui, SX, Nayak, KS and Wu, HH. Improved liver fat and R2* quantification at 0.55 T using locally low‐rank denoising. Magnetic Resonance in Medicine. https://doi.org/10.1002/mrm.30324

## Overview
Locally low-rank principal component analysis (PCA)-based denoising is one popular approach to suppress noise in multi-contrast MR images. By suppressing principal components associated with smaller coefficients, noise can be reduced while signal can be largely preserved. Here, we implemented two different denoising approaches that can accurately estimate the underlying noise variance and suppress noise.

### Method 1: Robust Locally Low-Rank (RLLR) denoising
This method was first proposed in Lugauer et al., MICCAI 2015 (doi: 10.1007/978-3-319-24571-3_80). Using samples of random matrices from a known Gaussian distribution, the noise level in the multi-echo images can be estimated. Based on Stein’s unbiased risk estimate (SURE), the singular value threshold can be objectively obtained for noise suppression. 

### Method 2: Random Matrix Theory (RMT) denoising
Previous works, such as Verrat et al., NeuroImage 2016 (doi: 10.1016/j.neuroimage.2016.08.016), have shown that Marchenko-Pastur Law in the field of random matrix theory can be used to accurately estimate noise level in low-rank image patches. Using this information, noise can be suppressed by using an optimal singular value shrinkage appraoch. 

## Requirement
(1) MATLAB (MathWorks, Natick, Massachusetts) version later than R2022a. \
(2) MATLAB Parallel Computing Toolbox (https://www.mathworks.com/products/parallel-computing.html).

## Usage
Run "denoising_experiment.m" in MATLAB. 
* This script load the sample k-space data, perform RLLR and RMT denoising and display the image results.
* Parameters: (1) dataset index (phantom or in vivo liver) and (2) patch size for denoising (e.g., [5,5,5])

## Dataset
We provide two example datasets to test the denoising performance: (1) a reference phantom containing proton-density fat fraction (PDFF)-only and R2*-only vials, and (2) an in vivo liver dataset. Both datasets were acquired using a multi-echo gradient-echo Cartesian Dixon sequence on a prototype 0.55T MRI system. Details on the acquisition parameters can be found in our paper.

## Acknowledgements
(1) We thank the Dynamic Imaging Science Center (DISC) at the University of Southern California for supporting data acquisition. 

(2) We acknowledge the use or adaptation of the following code and toolboxes. \
a. GRAPPA Reconstruction Tools by Mark Chiew (Link: https://github.com/mchiew/grappa-tools) \
b. im2row and row2im functions in ESPIRiT toolbox by Michael Lustig (Link: https://people.eecs.berkeley.edu/~mlustig/Software.html) \
c. sure_svt function (Stein's Unbiased Risk Estimate (SURE) for Singular Value Thresholding) by Emmanuel Candès (Link: https://candes.su.domains/software/SURE/data.html) \
d. optimal_shrinkage function (optimal singular value shrinkage) by Matan Gavish and David Donoho (Link: https://purl.stanford.edu/kv623gt2817) 

## Contact
Please contact Shu-Fu Shih (sshih@mednet.ucla.edu) or Holden Wu (holdenwu@mednet.ucla.edu) if you have any questions.
