function [result, gfactormap] = perform_GRAPPA(kdata, acs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Inputs:
% (1) kdata: The undersampled kspace data with dimension (Nfe x Npe x Nch x
%            Nsl x Ncon)
% (2) acs: The calibration kspce data with dimension (Nfe x Npe x Nch x
%          Nsl). Only the first contrast was used for ACS.
%
% Outputs:
% (1) result: Reconstructed images (Nfe x Npe x Nch x Nsl x Ncon)
% (2) gfactormap: Coil-by-coil g-factor maps (Nfe x Npe x Nch x Nsl)
%
% Acknowledgements: We used GRAPPA code provided by mchiew on GitHub 
%                  (https://github.com/mchiew/grappa-tools). This function 
%                  provided additional coil-by-coil g-factor maps.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Nfe,Npe,Nch,Nsl,Ncon] = size(kdata);

Raccel = [1,2];
kernel = [5,5];

result = zeros(Nfe,Npe,Nch,Nsl,Ncon);
gfactormap = zeros(Nfe,Npe,Nch,Nsl);

kdata = permute(kdata, [3,1,2,4,5]);
acs = permute(acs, [3,1,2,4]);

idx1 = 1:(Nch+1):(Nch*Nch);
idx2 = 0:1:(Nfe*Npe-1);
idx3 = (idx1 + idx2' * (Nch.^2)).';
idx3 = idx3(:);

parfor sl = 1:Nsl
for con = 1:Ncon
    
    % GRAPPA reconstruction
    [data, ~, w] = grappa_gfactor(kdata(:,:,:,sl,con), acs(:,:,:,sl), eye(Nch), Raccel, kernel);

    % Calculation of coil-by-coil g-factor maps
    q2 = pagemtimes(w,pagectranspose(w));
    q3 = sqrt(q2(idx3));
    q4 = reshape(q3,[Nch,Nfe,Npe]);
    gmap = permute(q4,[2,3,1]);

    % Save GRAPPA reconstrcted images and g-factor maps
    result(:,:,:,sl,con) = permute(data, [2,3,1]);
    gfactormap(:,:,:,sl) = gmap;

end
end
