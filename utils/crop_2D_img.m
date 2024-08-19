function res = crop_2D_img(img, new_size)

% Crop images to a certain size

orix = size(img,1);
oriy = size(img,2);
nx = new_size(1);
ny = new_size(2);

idx1 = floor((orix - nx)/2+1);
idx2 = idx1 + nx - 1;
idy1 = floor((oriy - ny)/2+1);
idy2 = idy1 + ny - 1;
    
res = img(idx1:idx2,idy1:idy2,:,:,:);

end