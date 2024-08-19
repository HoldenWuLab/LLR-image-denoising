function coilmaps_pad = interp2Ddata(coilmaps)

[Nx, Ny, Nsl, Nch] = size(coilmaps);

flag = 1;
count = 2;
while flag
    if coilmaps(1,count,1,1) == 0
        count = count + 1;
    else
        flag = 0;
    end
end
accel_fac = count - 1;

num = mod(Nx,accel_fac);
if num == 0
    padnum = Nx+1;
elseif num > 1
    padnum = Nx + accel_fac+1-num;
else
    padnum = Nx;
end

coilmaps_pad = zeros(padnum, padnum, Nsl, Nch);
coilmaps_pad(1:Nx,1:Ny,:,:) = coilmaps;

coilmaps_pad(end,:,:,:) = coilmaps_pad(end-accel_fac,:,:,:);
coilmaps_pad(:,end,:,:) = coilmaps_pad(:,end-accel_fac,:,:);

for ii = 2:accel_fac
    tmp1 = coilmaps_pad(1:accel_fac:(end-accel_fac),:,:,:) * (1-1/accel_fac*(ii-1));
    tmp2 = coilmaps_pad((accel_fac+1):accel_fac:end,:,:,:) * (1/accel_fac*(ii-1));
    coilmaps_pad(ii:accel_fac:(end-accel_fac-1+ii),:,:,:) = tmp1 + tmp2;

    tmp1 = coilmaps_pad(:,1:accel_fac:(end-accel_fac),:,:) * (1-1/accel_fac*(ii-1));
    tmp2 = coilmaps_pad(:,(accel_fac+1):accel_fac:end,:,:) * (1/accel_fac*(ii-1));
    coilmaps_pad(:,ii:accel_fac:(end-accel_fac-1+ii),:,:) = tmp1 + tmp2;
end

coilmaps_pad = coilmaps_pad(1:Nx,1:Ny,:,:);

end
