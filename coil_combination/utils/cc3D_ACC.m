function [cc_img3D, coilmaps3D] = cc3D_ACC(input_img ,bg_mask, sig_patch, fov, nloc3D, sl_batch, accel, use_nroi)

[Nx, Ny, Nsl, Ncon, Nch] = size(input_img);

% Iterate over all pixels
Nx_new = fov(1); 
Ny_new = fov(2);
s_roi_x_len = sig_patch(1);
s_roi_y_len = sig_patch(2);
b = (s_roi_x_len+1)/2;
Qs = s_roi_x_len * s_roi_y_len;

ddx = floor((Nx - Nx_new)/2);
ddy = floor((Ny - Ny_new)/2);

coilmaps3D = zeros(fov(1), fov(2), Nsl, Nch);
cc_img3D = zeros(fov(1), fov(2), Nsl, Ncon);

parfor sl = 1:Nsl

    sl3D_start = sl - (sl_batch-1)/2;
    sl3D_end = sl + (sl_batch-1)/2;

    if (sl3D_start < 1); sl3D_start = 1; end
    if (sl3D_end > Nsl); sl3D_end = Nsl; end

    if use_nroi == 1
        
        roi_num = size(nloc3D,1);
        Ri = zeros(Nch,Nch);

        for ssl = sl3D_start:sl3D_end
            for nn = 1:roi_num
        
                nx1 = nloc3D(nn,1,ssl);
                nx2 = nloc3D(nn,3,ssl);
                ny1 = nloc3D(nn,2,ssl);
                ny2 = nloc3D(nn,4,ssl);
        
                patch_size_x = abs( nx2 - nx1 ) +1;
                patch_size_y = abs( ny2 - ny1 ) +1;
        
                Qi = patch_size_x * patch_size_y;
                for con = 1:Ncon
                    noise_blk = squeeze(input_img(nx1:nx2, ny1:ny2, ssl, con, :));
                    noise_blk = reshape(noise_blk, [Qi, Nch]);
                    Ri = Ri + noise_blk.' * conj(noise_blk);
                end
            end
        end
        Ri = Ri/(roi_num*Ncon*(sl3D_end-sl3D_start+1));
        
        % To improve condition number for robust matrix inversion
        lambda = 0.01;
        [~,S,~] = svd(Ri);
        M = eye(Nch);
        inv_Ri = inv(Ri + lambda*S(1,1)*M);
    else
        inv_Ri = eye(Nch);
    end


    for i = 1:accel:Nx_new
    for j = 1:accel:Ny_new
            
        if bg_mask(i,j,sl) == 1
    
            nx1 = (i+ddx) - (s_roi_x_len-1)/2;
            nx2 = nx1 + s_roi_x_len - 1;
            ny1 = (j+ddy) - (s_roi_y_len-1)/2;
            ny2 = ny1 + s_roi_y_len - 1;
                
            % Extract the patch for signal covariance calculation
            if (nx1>b) && (nx2<=Nx-b) && (ny1>b) && (ny2<=Ny-b)
                    
                Rs = zeros(Nch, Nch);
                for con = 1:Ncon
                    blk_mat = reshape(squeeze(input_img(nx1:nx2, ny1:ny2, sl, con, :)), [Qs, Nch]);
                    Rs = Rs + blk_mat.' * conj(blk_mat);             
                end
                    
            else
                % Solve the boundary issue
                if (nx1 < 1);  nx1 = 1;  end
                if (nx2 > Nx); nx2 = Nx; end
                if (ny1 < 1);  ny1 = 1;  end
                if (ny2 > Ny); ny2 = Ny; end
                    
                Rs = zeros(Nch, Nch);
                for con = 1:Ncon
                    blk_mat = squeeze(input_img(nx1:nx2, ny1:ny2, sl, con, :));
                    blk_mat = reshape(blk_mat, [size(blk_mat,1)*size(blk_mat,2), Nch]);
                    Rs = Rs + blk_mat.' * conj(blk_mat);
                end
            end
                
            [U,~,~] = svd(inv_Ri * Rs);
        
            w = U(:,1);
            coilmaps3D(i,j,sl,:) = w;
            
    
        else
            coilmaps3D(i,j,sl,:) = 0;
            cc_img3D(i,j,sl,:) = 0;
        end
            
    end
    end
end

coilmaps3D = interp2Ddata(coilmaps3D);
coilmaps3D = reshape(coilmaps3D, [fov(1), fov(2), Nsl, 1, Nch]);
%coilmaps3D = coilmaps3D./sum(abs(coilmaps3D),5);

cc_img3D = sum(conj(coilmaps3D).*input_img(ddx+1:ddx+fov(1), ddy+1:ddy+fov(2), :, :, :), 5);

end