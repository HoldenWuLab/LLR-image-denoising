function [kdata_decor,acs_decor] = noise_decor(kdata,acs,prescan)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This NOISE_DECOR function perform noise decorrelation across the
% multi-coil k-space and acs data. This algorithm was implemented based on
% the research work from Pruessmann et al. MRM 2001 (doi: 10.1002/mrm.1241)
%
% Input:(1) kdata:   This is the acquired multi-coil multi-contrast k-space
%                    data with dimensions [Nkx,Nch,Nky,Nkz,Nec].
%       (2) acs:     This is the acquired autocalibration signal data for 
%                    parallel imaging reconstruction. It has the dimensions 
%                    [Nkx,Nch,Nky_acs,Nkz,Nec].
%       (3) prescan: This is the prescan data for noise decorrelation
%                    across coil channels. It has the dimensions 
%                    [Nkx_prescan,Nch,Nky_prescan,Navg].
%
% Output:(1) kdata_decor: Noise-decorrelated k-space data
%        (2) acs_decor: Noise-decorrelated acs data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Nkx,Nch,Nky,Navg] = size(prescan);

prescan = permute(prescan,[2,1,3,4]);
prescan = reshape(prescan,[Nch,Nkx*Nky*Navg]);

% Please see Appendix in Pruessmann et al. MRM 2001
psi = prescan*prescan';
psi = psi./mean(abs(diag(psi)));

% Cholesky decomposition (Equation B4 in Pruessmann et al. MRM 2001)
Lmat = chol(psi,'lower');

% Perform noise decorrelation in k-space data
kdata = permute(kdata, [2,1,3,4,5]);
kdata_decor = pagemldivide(Lmat,kdata);
kdata_decor = permute(kdata_decor, [2,1,3,4,5]);

% Perform noise decorrelation in acs data
acs = permute(acs, [2,1,3,4,5]);
acs_decor = pagemldivide(Lmat,acs);
acs_decor = permute(acs_decor, [2,1,3,4,5]);
