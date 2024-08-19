function k = iFT_along_specific_dim(k, dim)

% Perform inverse Fourier transform along one specific dimension

k = fftshift( ifft(  ifftshift(k, dim),  [], dim ),dim);

