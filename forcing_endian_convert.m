big='b'; % Big endian data
little='l'; % Little endian
accuracy='real*4';

nx=128;
ny=64;
nt=12;
nz=15;

dowrite=0; % Control for writing data or just processing it....

% Bathymetry
fid=fopen('depth_g77.bin','r',big); depth_be=fread(fid,accuracy); fclose(fid);
depth_le=nan(nx,ny);
for i=1:64
    depth_le(:,i)=depth_be(((i-1)*128)+1:i*128);
end
if dowrite==1
fid=fopen('bathy_le.bin','w',little);fwrite(fid,depth_le,accuracy); fclose(fid);
end

% Sea Ice
fid=fopen('fice.bin','r',big); fice_be=fread(fid,accuracy); fclose(fid);
fice_le=nan(nx,ny,nt);
for m=1:12
    for i=1:64
        fice_le(:,i,m)=fice_be((8192*(m-1))+((i-1)*128)+1:(8192*(m-1))+i*128);
    end
end
if dowrite==1
fid=fopen('fice_le.bin','w',little);fwrite(fid,fice_le,accuracy); fclose(fid);
end

% Levitus Salinity Climatology
fid=fopen('lev_clim_salt.bin','r',big); lev_salt_be=fread(fid,accuracy); fclose(fid);
lev_salt_le=nan(nx,ny,nz);
for m=1:15
    for i=1:64
        lev_salt_le(:,i,m)=lev_salt_be((8192*(m-1))+((i-1)*128)+1:(8192*(m-1))+i*128);
    end
end
if dowrite==1
fid=fopen('lev_salt_le.bin','w',little);fwrite(fid,lev_salt_le,accuracy); fclose(fid);
end

% Levitus Potential Temperature Climatology
fid=fopen('lev_clim_temp.bin','r',big); lev_temp_be=fread(fid,accuracy); fclose(fid);
lev_temp_le=nan(nx,ny,nz);
for m=1:15
    for i=1:64
        lev_temp_le(:,i,m)=lev_temp_be((8192*(m-1))+((i-1)*128)+1:(8192*(m-1))+i*128);
    end
end
if dowrite==1
fid=fopen('lev_temp_le.bin','w',little);fwrite(fid,lev_temp_le,accuracy); fclose(fid);
end

% Levitus Monthly SSSalinity Climatology
fid=fopen('lev_monthly_salt.bin','r',big); lev_salt_be=fread(fid,accuracy); fclose(fid);
lev_salt_le=nan(nx,ny,nt);
for m=1:12
    for i=1:64
        lev_salt_le(:,i,m)=lev_salt_be((8192*(m-1))+((i-1)*128)+1:(8192*(m-1))+i*128);
    end
end
if dowrite==1
fid=fopen('lev_monthly_salt_le.bin','w',little);fwrite(fid,lev_salt_le,accuracy); fclose(fid);
end

% Levitus Potential Temperature Climatology
fid=fopen('lev_monthly_temp.bin','r',big); lev_temp_be=fread(fid,accuracy); fclose(fid);
lev_temp_le=nan(nx,ny,nt);
for m=1:12
    for i=1:64
        lev_temp_le(:,i,m)=lev_temp_be((8192*(m-1))+((i-1)*128)+1:(8192*(m-1))+i*128);
    end
end
if dowrite==1
fid=fopen('lev_monthly_temp_le.bin','w',little);fwrite(fid,lev_temp_le,accuracy); fclose(fid);
end

% Net Surface HF Climatology
fid=fopen('shi_qnet.bin','r',big); qnet_be=fread(fid,accuracy); fclose(fid);
qnet_le=nan(nx,ny,nt);
for m=1:12
    for i=1:64
        qnet_le(:,i,m)=qnet_be((8192*(m-1))+((i-1)*128)+1:(8192*(m-1))+i*128);
    end
end
if dowrite==1
fid=fopen('shi_qnet_le.bin','w',little);fwrite(fid,qnet_le,accuracy); fclose(fid);
end

% Evap - Precip - Runoff Climatology
fid=fopen('shi_empmr_year.bin','r',big); emp_be=fread(fid,accuracy); fclose(fid);
emp_le=nan(nx,ny,nt);
for m=1:12
    for i=1:64
        emp_le(:,i,m)=emp_be((8192*(m-1))+((i-1)*128)+1:(8192*(m-1))+i*128);
    end
end
if dowrite==1
fid=fopen('shi_empmr_year_le.bin','w',little);fwrite(fid,emp_le,accuracy); fclose(fid);
end

% Sea surface Silica Climatology
fid=fopen('sillev1.bin','r',big); si_be=fread(fid,accuracy); fclose(fid);
si_le=nan(nx,ny,nt);
for m=1:12
    for i=1:64
        si_le(:,i,m)=si_be((8192*(m-1))+((i-1)*128)+1:(8192*(m-1))+i*128);
    end
end
if dowrite==1
fid=fopen('sillev1_le.bin','w',little);fwrite(fid,si_le,accuracy); fclose(fid);
end

% Trenberth Wiind Speed Climatology
fid=fopen('tren_speed.bin','r',big); sp_be=fread(fid,accuracy); fclose(fid);
sp_le=nan(nx,ny,nt);
for m=1:12
    for i=1:64
        sp_le(:,i,m)=sp_be((8192*(m-1))+((i-1)*128)+1:(8192*(m-1))+i*128);
    end
end
if dowrite==1
fid=fopen('tren_speed_le.bin','w',little);fwrite(fid,sp_le,accuracy); fclose(fid);
end

% Trenberth Zonal Wind Stress Climatology
fid=fopen('tren_taux.bin','r',big); taux_be=fread(fid,accuracy); fclose(fid);
taux_le=nan(nx,ny,nt);
for m=1:12
    for i=1:64
        taux_le(:,i,m)=taux_be((8192*(m-1))+((i-1)*128)+1:(8192*(m-1))+i*128);
    end
end
if dowrite==1
fid=fopen('tren_taux_le.bin','w',little);fwrite(fid,taux_le,accuracy); fclose(fid);
end

% Trenberth Merid. Wind Stress Climatology
fid=fopen('tren_tauy.bin','r',big); tauy_be=fread(fid,accuracy); fclose(fid);
tauy_le=nan(nx,ny,nt);
for m=1:12
    for i=1:64
        tauy_le(:,i,m)=tauy_be((8192*(m-1))+((i-1)*128)+1:(8192*(m-1))+i*128);
    end
end
if dowrite==1
fid=fopen('tren_tauy_le.bin','w',little);fwrite(fid,tauy_le,accuracy); fclose(fid);
end
%clear *le *be m i nx ny nz nt fid big little accuracy