function [res,W] = row2vol(mtx,imSize, winSize)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Acknowledgements: We adpated the functions "im2row.m" and "row2im.m" from 
%                the ESPIRiT toolbox provided by Michael Lustig 
%                (https://people.eecs.berkeley.edu/~mlustig/Software.html).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

coil = size(mtx,3);
sx = imSize(1); 
sy = imSize(2);
sz = imSize(3);

res = zeros(imSize(1),imSize(2),imSize(3),coil);
W = res;

count = 0;
for z = 1:winSize(3)
for y = 1:winSize(2)
for x = 1:winSize(1)
    count = count+1;
    res(x:sx-winSize(1)+x,y:sy-winSize(2)+y,z:sz-winSize(3)+z,:) = res(x:sx-winSize(1)+x,y:sy-winSize(2)+y,z:sz-winSize(3)+z,:) + reshape(mtx(:,count,:),(sx-winSize(1)+1),(sy-winSize(2)+1),(sz-winSize(3)+1),coil);
    W(x:sx-winSize(1)+x,y:sy-winSize(2)+y,z:sz-winSize(3)+z,:) = W(x:sx-winSize(1)+x,y:sy-winSize(2)+y,z:sz-winSize(3)+z,:)+1;
end
end
end

res = res./W;