function ccImgDns = RLLR_denoising(ccImg,dns_patch_size)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script performs RLLR denoising 
%
% Input: (1) ccImg: coil-combined complex-valued multi-echo images. The 
%                   data dimension should be [Nx,Ny,Nsl,Ncon].
%        (2) dnsPatchSize: patch size used for denoising.
%
% Output: (1) ccImgDns: denoised result.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get data dimensions
[Nx,Ny,Nsl,Ncon] = size(ccImg);
px = dns_patch_size(1);
py = dns_patch_size(2);
pz = dns_patch_size(3);

% Group all local patches
hankel_mat = vol2row(ccImg, dns_patch_size);
hankel_mat = permute(hankel_mat,[3,2,1]);

% The number of local patches
Ninst = size(hankel_mat,3);

% Generate Gaussian noise samples
% The median of the last singular values was caluculated
num_simu_instance = 1000;
simu_sval_arr = zeros(1,num_simu_instance);
samp_sigma = 1;
for inst = 1:num_simu_instance
    random_seed_real = 12345 + inst * 7;
    rng(random_seed_real)
    simuReal = samp_sigma * randn(px*py*pz,Ncon);
    random_seed_imag = 23456 + inst * 7;
    rng(random_seed_imag)
    simuImag = samp_sigma * randn(px*py*pz,Ncon);
    simuMat = simuReal + 1j*simuImag;
    [~,s,~] = svd(simuMat);
    simu_sval_arr(inst) = s(Ncon,Ncon);
end
sval_median = median(simu_sval_arr);

% Perform denoising for all local patches
parfor nn = 1:Ninst
    
    % Perform SVD
    local_patch = hankel_mat(:,:,nn);
    [u,d,v] = svd(local_patch,"econ");
    singVal =  diag(d);

    % Estimate noise variance
    estSigma = d(Ncon,Ncon)/sval_median*samp_sigma;
    
    % Perform singular value thresholding
    lam_min = singVal(end)*0.5;
    lam_max = singVal(2)*1.5;
    numRisk = 100;
    riskArr = zeros(1,numRisk);
    searchGrid = linspace(lam_min,lam_max,numRisk);
    count = 1;
    for jj = searchGrid
        % Calculate SURE risk
        riskArr(count) = sure_svt(jj, estSigma, singVal, [px*py*pz,Ncon], 0);
        count = count + 1;
    end
    [minRisk,minRiskIdx] = min(riskArr);
    optimLam = searchGrid(minRiskIdx);
    singValThresh = (singVal - optimLam).*((singVal - optimLam)>0);
    hankel_mat(:,:,nn) = u*diag(singValThresh)*(v');

end

% Reshape data to the original data dimensions
hankel_mat = permute(hankel_mat,[3,2,1]);
ccImgDns = row2vol(hankel_mat,[Nx,Ny,Nsl],dns_patch_size);
