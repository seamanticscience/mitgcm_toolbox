big='b'; % Big endian data
little='l'; % Little endian
accuracy='real*4';

nx=128;
ny=64;
nt=12;
nz=15;

% Dasilva Monthly SST Climatology
fid=fopen('dasilva_sst.bin','r',little); tmp=fread(fid,accuracy); fclose(fid);
da_sst=nan(nx,ny,nt);
for m=1:12
    for i=1:64
        da_sst(:,i,m)=tmp((8192*(m-1))+((i-1)*128)+1:(8192*(m-1))+i*128);
    end
end

% Levitus SST Climatology
fid=fopen('lev_monthly_temp_le.bin','r',little); tmp=fread(fid,accuracy); fclose(fid);
lev_sst=nan(nx,ny,nt);
for m=1:12
    for i=1:64
        lev_sst(:,i,m)=tmp((8192*(m-1))+((i-1)*128)+1:(8192*(m-1))+i*128);
    end
end

% Dasilva Net Surface HF Climatology
fid=fopen('dasilva_qnet.bin','r',little); tmp=fread(fid,accuracy); fclose(fid);
da_qnet=nan(nx,ny,nt);
for m=1:12
    for i=1:64
        da_qnet(:,i,m)=tmp((8192*(m-1))+((i-1)*128)+1:(8192*(m-1))+i*128);
    end
end

% Net Surface HF Climatology
fid=fopen('shi_qnet_le.bin','r',little); tmp=fread(fid,accuracy); fclose(fid);
shi_qnet=nan(nx,ny,nt);
for m=1:12
    for i=1:64
        shi_qnet(:,i,m)=tmp((8192*(m-1))+((i-1)*128)+1:(8192*(m-1))+i*128);
    end
end

% Dasilva Evap - Precip - Runoff Climatology
fid=fopen('dasilva_emp.bin','r',little); tmp=fread(fid,accuracy); fclose(fid);
da_emp=nan(nx,ny,nt);
for m=1:12
    for i=1:64
        da_emp(:,i,m)=tmp((8192*(m-1))+((i-1)*128)+1:(8192*(m-1))+i*128);
    end
end

% Evap - Precip - Runoff Climatology
fid=fopen('shi_empmr_year_le.bin','r',little); tmp=fread(fid,accuracy); fclose(fid);
shi_emp=nan(nx,ny,nt);
for m=1:12
    for i=1:64
        shi_emp(:,i,m)=tmp((8192*(m-1))+((i-1)*128)+1:(8192*(m-1))+i*128);
    end
end

% Dasilva Wiind Speed Climatology
fid=fopen('dasilva_speed.bin','r',little); tmp=fread(fid,accuracy); fclose(fid);
da_sp=nan(nx,ny,nt);
for m=1:12
    for i=1:64
        da_sp(:,i,m)=tmp((8192*(m-1))+((i-1)*128)+1:(8192*(m-1))+i*128);
    end
end

% Trenberth Wiind Speed Climatology
fid=fopen('tren_speed_le.bin','r',little); tmp=fread(fid,accuracy); fclose(fid);
tren_sp=nan(nx,ny,nt);
for m=1:12
    for i=1:64
        tren_sp(:,i,m)=tmp((8192*(m-1))+((i-1)*128)+1:(8192*(m-1))+i*128);
    end
end

% Dasilva Zonal Wind Stress Climatology
fid=fopen('dasilva_taux.bin','r',little); tmp=fread(fid,accuracy); fclose(fid);
da_taux=nan(nx,ny,nt);
for m=1:12
    for i=1:64
        da_taux(:,i,m)=tmp((8192*(m-1))+((i-1)*128)+1:(8192*(m-1))+i*128);
    end
end

% Trenberth Zonal Wind Stress Climatology
fid=fopen('tren_taux_le.bin','r',little); tmp=fread(fid,accuracy); fclose(fid);
tren_taux=nan(nx,ny,nt);
for m=1:12
    for i=1:64
        tren_taux(:,i,m)=tmp((8192*(m-1))+((i-1)*128)+1:(8192*(m-1))+i*128);
    end
end

% Dasilva Merid. Wind Stress Climatology
fid=fopen('dasilva_tauy.bin','r',little); tmp=fread(fid,accuracy); fclose(fid);
da_tauy=nan(nx,ny,nt);
for m=1:12
    for i=1:64
        da_tauy(:,i,m)=tmp((8192*(m-1))+((i-1)*128)+1:(8192*(m-1))+i*128);
    end
end

% Trenberth Merid. Wind Stress Climatology
fid=fopen('tren_tauy_le.bin','r',little); tmp=fread(fid,accuracy); fclose(fid);
tren_tauy=nan(nx,ny,nt);
for m=1:12
    for i=1:64
        tren_tauy(:,i,m)=tmp((8192*(m-1))+((i-1)*128)+1:(8192*(m-1))+i*128);
    end
end

clear tmp m i nx ny nz nt fid little big accuracy
%clear *le *be m i nx ny nz nt fid little little accuracy