
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script performs RLLR denoising and RMT denoising on multi-echo
% gradient echo datasets acquired at a prototype 0.55T scanner.
%
% We provide two datasets: phantom and in vivo liver. Details of the 
% acquisition parameters can be found in our paper. 
%
% Acknowledgement: We acknowledge the use or adaptation of the following 
% code and toolboxes.
%    a. GRAPPA Reconstruction Tools by Mark Chiew
%       (Link: https://github.com/mchiew/grappa-tools)
%    b. im2row and row2im functions in ESPIRiT toolbox by Michael Lustig
%       (Link: https://people.eecs.berkeley.edu/~mlustig/Software.html)
%    c. sure_svt function [Stein's Unbiased Risk Estimate (SURE) for 
%       Singular Value Thresholding (SVT).] by Emmanuel Candès
%       (Link: https://candes.su.domains/software/SURE/data.html)
%    d. optimal_shrinkage function [Optimal singular value shrinkage] by
%       Matan Gavish and David Donoho
%       (Link: https://purl.stanford.edu/kv623gt2817)
%
%
% (c) by Shu-Fu Shih, August 2024.
% Magnetic Resonance Research Labs - Wu Lab
% University of California, Los Angeles
%
% If you have any questions, please contact:
%    Shu-Fu Shih (sshih@mednet.ucla.edu)
%    Holden Wu (holdenwu@mednet.ucla.edu)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc;

%% Section A: Parameters

% Param1: Dataset to load
dataset_idx = 1; % 1: phantom; 2: liver

% Param 2: Patch size for low-rank matrix construction in x-y-z dimensions
dns_patch_size = [5,5,5]; 

%% Section B: Set up path and parallel imaging toolbox

% Add paths
addpath(pwd)
addpath(genpath("./coil_combination"))
addpath(genpath("./data_loading"))
addpath(genpath("./denoising"))
addpath(genpath("./grappa"))
addpath(genpath("./noise_decorrelation"))
addpath(genpath("./utils"))
addpath(genpath("./example_dataset"))

% Check if parpool can be started
p = gcp('nocreate');
if isempty(p)
    try
        parpool(12);
        disp('Pararrel pool started successfully');
        useParPool = 1;
    catch 
        warning('Failed to start parallel pool');
        useParPool = 0;
    end
end

%% Section C: Experiment

% Load raw data
[kdata,acs,prescan] = data_loading(dataset_idx);
fprintf('Finish loading data.\n');

% Perform noise decorrelation
[kdata,acs] = noise_decor(kdata,acs,prescan);
fprintf('Finish noise decorrelation.\n');

% Reshape the kdata and refscan to Nfe x Npe x Nch x Nsl x Ncon
kdata = permute(kdata, [1,3,2,4,5]);
acs = permute(acs,[1,3,2,4,5]);
[Nkx, Nky, Nch, Nkz, Ncon] = size(kdata);
[~, Nky_acs, ~, ~, ~] = size(acs);

% Perform inverse Fourier Transform along kz
kdataPartRes = iFT_along_specific_dim(kdata,4);
acsPartRes = iFT_along_specific_dim(acs,4);

% Perform GRAPPA reconstruction and calculate g-factor maps
[img, gFactorMap] = perform_GRAPPA(kdataPartRes, acsPartRes(:,:,:,:,1));
img = img./gFactorMap;
img = imrotate(img,90);
img = permute(img,[1,2,4,3,5]);
fprintf('Finish GRAPPA reconstruction.\n');

% Remove oversampling
img = crop_2D_img(img,[Nky,Nky]);

% Method 1: RLLR denoising
% This approach performs coil combination before denoising
fprintf('Start RLLR denoising...(can take several minutes)\n')
tic;
[ccImg, ~] = coil_combine(img, 'xyzcn', 'ACC','patchSize' ,[15,15]); % Coil combination
ccImg = squeeze(ccImg);
ccImg = ccImg * 1e8; % Scaling to avoid numerical errors
ccImgRLLR = RLLR_denoising(ccImg,dns_patch_size); % RLLR denoising
ccImg = ccImg * 1e-8; % Scaling back
ccImgRLLR = ccImgRLLR * 1e-8; % Scaling back
time_1 = toc;
fprintf('Complete RLLR denoising\n')
fprintf('RLLR denoising time: %f sec\n',time_1)

% Method 2: RMT denoising
% This approach performs coil combination after denoising
fprintf('Start RMT denoising...(can take several minutes)\n')
tic;
img_dns = RMT_denoising(img,dns_patch_size); % RMT denoising
[ccImgRMT, ~] = coil_combine(img_dns, 'xyzcn', 'ACC','patchSize',[15,15]); % Coil combination
ccImgRMT = squeeze(ccImgRMT);
time_2 = toc;
fprintf('Complete RMT denoising\n')
fprintf('RMT denoising time: %f sec\n',time_2)

% Display results
sl = 8; % slice to display
ec = 3; % echo to display
scaling = max(abs(ccImg(:)))*0.5;
figure;imshow(abs(cat(2,ccImg(:,:,sl,ec),ccImgRLLR(:,:,sl,ec),ccImgRMT(:,:,sl,ec))),[0,scaling]);
title('Magnitude image: 1) No denoising, 2) RLLR denoising and 3) RMT denoising')
figure;imshow(angle(cat(2,ccImg(:,:,sl,ec),ccImgRLLR(:,:,sl,ec),ccImgRMT(:,:,sl,ec))),[-pi,pi]);
title('Phase image: 1) No denoising, 2) RLLR denoising and 3) RMT denoising')
