function img_dns = RMT_denoising(img,dns_patch_size)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script performs RMT denoising 
%
% Input: (1) ccImg: complex-valued multi-coil multi-echo images. The data
%                   dimension should be [Nx,Ny,Nsl,Nch,Ncon].
%        (2) dnsPatchSize: patch size used for denoising.
%
% Output: (1) img_dns: denoised result.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get data dimensions
[Nx,Ny,Nsl,Nch,Ncon] = size(img);

% Group all local patches
hankel_mat = vol2row(img, dns_patch_size);
hankel_mat = permute(hankel_mat,[3,2,1]);

% The number of local patches
Ninst = size(hankel_mat,3);

% Set up the matrix dimension for one local patch
% beta_ratio is the matrix aspect ratio
m1 = size(hankel_mat,1)*2; % multiply by 2 because it deals with real and imaginary components
n1 = size(hankel_mat,2);
if m1/n1 > 1
    beta_ratio = n1/m1;
else
    beta_ratio = m1/n1;
end

% Perform denoising for all local patches
parfor nn = 1:Ninst

    % Perform SVD
    X = hankel_mat(:,:,nn);
    local_patch = cat(1,real(X),imag(X));
    [u,d,v] = svd(local_patch,"econ");
    
    % Perform singular value shrinkage
    [opt_d,sigma] = optimal_shrinkage(diag(d),beta_ratio,'fro');
    
    % Save denoised matrix
    local_patch_dns = u * diag(opt_d) * v';
    hankel_mat(:,:,nn) = local_patch_dns(1:end/2,:) + 1j*local_patch_dns(end/2+1:end,:);
end

% Reshape back to the original data dimensions
hankel_mat = permute(hankel_mat,[3,2,1]);
img_dns = row2vol(hankel_mat,[Nx,Ny,Nsl],dns_patch_size);
img_dns = reshape(img_dns,[Nx,Ny,Nsl,Nch,Ncon]);
