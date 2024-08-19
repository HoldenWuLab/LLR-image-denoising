function res = vol2row(im, kernelSize)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Acknowledgements: We adpated the functions "im2row.m" and "row2im.m" from 
%                the ESPIRiT toolbox provided by Michael Lustig 
%                (https://people.eecs.berkeley.edu/~mlustig/Software.html).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[sx,sy,sz,coil] = size(im);

res = zeros((sx-kernelSize(1)+1)*(sy-kernelSize(2)+1)*(sz-kernelSize(3)+1),prod(kernelSize),coil);
count = 0;

for z = 1:kernelSize(3)
for y = 1:kernelSize(2)
for x = 1:kernelSize(1)
    count = count+1;
    res(:,count,:) = reshape(im(x:sx-kernelSize(1)+x,y:sy-kernelSize(2)+y,z:sz-kernelSize(3)+z,:),...
        (sx-kernelSize(1)+1)*(sy-kernelSize(2)+1)*(sz-kernelSize(3)+1),1,coil);
end
end
end