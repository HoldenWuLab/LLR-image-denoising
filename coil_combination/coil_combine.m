function [ccImg, varargout] = coil_combine(img, dimStr, combineMethod, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function performs coil combination on multi-coil data. It supports
% 2D, 2D multislice and 3D data. It provides different coil combination
% methods. It also supports acceleration throught interpolation.
%
% INPUTS:
% (1) img: Multi-coil image data you want to do coil combination
%
% (2) dimStr: A string that explains the dimension of the 'img' data. If
% your data is (Nx,Ny,Nsl,Ncon,Nch), you should input dimStr as 'xyznc'. If
% the data dimension is (Nx,Ny,Ncon,Nsl,Nch,Nave). You should input dimStr 
% as 'xynzca'.
%           'x' -> image x dimension
%           'y' -> image y dimension
%           'z' -> image z dimension (only for 3D data)
%           'c' -> image coil dimension
%           'n' -> image contrast dimension
%           'a' -> it means the algorithm will process data along this 
%                  dimension separately (e.g., different MRE phase offsets 
%                  or the slice dimension if you have 2D multi-slice data) 
%   
% (3) combinedMethod: There are 4 different methods.
%     'SOS'  : Sum-of-squares
%     'ACC'  : Adaptive coil combine
%
% (4) Other optional fields:
%     'patchSize'    : The local patch size to estimate signal covariances.
%                      The default is [11,11]. A larger patch will lead to
%                      smoother coil sensitivity maps and longer
%                      computational time.
%     'fov'          : The field of view you want for your output images
%                      and coil sensitivity maps.
%     'slBatch'      : The number of neighboring slices for interference
%                      covariance calculation. The default is 3 for 3D
%                      data, is 1 for 2D data.
%     'accel'        : The acceleration factor for 'ACC', 'BF', 'MVDR'
%                      calculation. The default is 2, which has 2x2
%                      interpolation with 4-fold acceleration.
%     'bgMask'       : You can skip some voxels by adding a binary bgMask. 
%                      1 represents the voxels you want to calculate; 0 
%                      represents the voxels you want to skip calculation.
%                      The size should be (Nx,Ny,Nsl). 
%
% Outputs:
% (1) ccImg: same dimension as img, except that the coil dimension is one.
% (2) coilmaps: same dimension as img, except that the contrast dimension
%     is one.
%
%
% Examples:
%
% This function supports data with different data arrangment. If you have
% data with dimension (Nx,Ny,Nch,Ncon,Nsl) = (256,256,15,6,10), and want to
% do sum-of-squares coil combination. 
%  >> [ccImg] = coilCombine(img, 'xycnz', 'SOS')
% 
% If you want to do adaptive coil combine and has data dimension 
% (Nx,Ny,Nsl,Nch) = (256,256,10,15).
%  >> [ccImg, coilmaps] = coilCombine(img, 'xyzc', 'ACC', 'patchSize', [13,13], 'fov', [128,128])
%
%
% If you have any questions, 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Adjust the "img" dimension to (Nx,Ny,Nz,Ncon,Nch,Na)
order = getCoilCombineDim(dimStr);
[~, order_rev] = sort(order,'ascend');
img = permute(img, order);

[Nx,Ny,Nz,Ncon,Nch,Nrpt] = size(img);

% Initialize opts with default values
opts.patchSize = [11,11];
opts.fov = [Nx,Ny];
opts.noisePatchLoc = [];
if Nz == 1
    opts.slBatch = 1;
else
    opts.slBatch = 3;
end
opts.accel = 2;
opts.bgMask = ones(Nx,Ny,Nz);

% varargin handling (must be option/value pairs)
for k = 1:2:numel(varargin)
    if k==numel(varargin) || ~ischar(varargin{k})
        error('''varargin'' must be option/value pairs.');
    end
    if ~isfield(opts,varargin{k})
        warning('''%s'' is not a valid option.',varargin{k});
    end
    opts.(varargin{k}) = varargin{k+1};
end

% Load to local variables
patchSize = opts.patchSize;
fov = opts.fov;
noisePatchLoc = opts.noisePatchLoc;
slBatch = opts.slBatch;
accel = opts.accel;
bgMask = opts.bgMask;


% Perfrom coil combination
if strcmp(combineMethod, 'SOS')
    
    ccImg = zeros(Nx,Ny,Nz,Ncon,1,Nrpt);
    if Nz ~= 1 % 3D SOS
        for con = 1:Ncon
        for r = 1:Nrpt
            imgTemp = squeeze(img(:,:,:,con,:,r));
            ccImg(:,:,:,con,:,r) = sqrt(sum(abs(imgTemp).^2, 4));
        end
        end
    else % 2D SOS
        for con = 1:Ncon
        for r = 1:Nrpt
            imgTemp = squeeze(img(:,:,1,con,:,r));
            ccImg(:,:,1,con,:,r) = sqrt(sum(abs(imgTemp).^2, 3));
        end
        end
    end

elseif strcmp(combineMethod, 'ACC')

    use_nroi = 0;
    if Nrpt == 1
        [ccImg, coilmaps3D] = cc3D_ACC(img, bgMask, patchSize, fov, noisePatchLoc, slBatch, accel, use_nroi);
        varargout{1} = coilmaps3D;
    else
        ccImg = zeros(fov(1),fov(2),Nz,Ncon,1,Nrpt);
        coilmaps3D = zeros(fov(1),fov(2),Nz,1,Nch,Nrpt);
        for r = 1:Nrpt
            imgTemp = img(:,:,:,:,:,r);
            [ccImg(:,:,:,:,:,r), coilmaps3D(:,:,:,:,:,r)] = cc3D_ACC(imgTemp, bgMask, patchSize, fov, noisePatchLoc, slBatch, accel, use_nroi);
        end
    end

end
    

ccImg = permute(ccImg, order_rev);
if ~strcmp(combineMethod,'SOS')
    coilmaps3D = permute(coilmaps3D, order_rev);
    varargout{1} = coilmaps3D;
end
