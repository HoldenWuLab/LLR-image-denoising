function data = ifftdim(data, dims)

for a = 1:length(dims)

    dim = dims(a);

    data = fftshift( ifft(  ifftshift(data, dim),  [], dim ),dim); 

end