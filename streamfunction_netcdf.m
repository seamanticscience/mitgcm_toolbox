function [filename]=streamfunction_netcdf(baro_psi,global_psi,atlantic_psi,pacific_psi,...
                                          baro_eul,global_eul,atlantic_eul,pacific_eul,...
                                          baro_eddys,global_eddys,atlantic_eddys,pacific_eddys,...
                                          X,Y,Z,tavesteps,attributes)

% 
% X=Xp1;
% Y=Yp1;
% Z=gridtmp.zgpsi;
gridtmp.nx=length(X); gridtmp.ny=length(Y); gridtmp.nz=length(Z);

filename='psi_data.nc';
nc_create_empty(filename)
nc_adddim(filename,'Xpsi',gridtmp.nx)
nc_adddim(filename,'Ypsi',gridtmp.ny)
nc_adddim(filename,'Zpsi',gridtmp.nz)
nc_adddim(filename,'T',0)

nc_addvar(filename,struct('Name','Xpsi','Datatype','double','Dimension',{{'Xpsi'}}))
nc_addvar(filename,struct('Name','Ypsi','Datatype','double','Dimension',{{'Ypsi'}}))
nc_addvar(filename,struct('Name','Zpsi','Datatype','double','Dimension',{{'Zpsi'}}))
nc_addvar(filename,struct('Name','T','Datatype','double','Dimension',{{'T'}}))

% Write to NetCDF file
nc_varput(filename,'Xpsi',X);
nc_varput(filename,'Ypsi',Y);
nc_varput(filename,'Zpsi',Z);

nc_addvar(filename,struct('Name','iter','Datatype','int','Dimension',{{'T'}}));
nc_addvar(filename,struct('Name','gbaro_eul','Datatype','float','Dimension',{{'T','Ypsi','Xpsi'}}));
nc_addvar(filename,struct('Name','gbaro_psi','Datatype','float','Dimension',{{'T','Ypsi','Xpsi'}}));
nc_addvar(filename,struct('Name','gbaro_eddy_psi','Datatype','float','Dimension',{{'T','Ypsi','Xpsi'}}));
nc_addvar(filename,struct('Name','global_psi','Datatype','float','Dimension',{{'T','Zpsi','Ypsi'}}))
nc_addvar(filename,struct('Name','global_eul','Datatype','float','Dimension',{{'T','Zpsi','Ypsi'}}))
nc_addvar(filename,struct('Name','atlantic_psi','Datatype','float','Dimension',{{'T','Zpsi','Ypsi'}}));
nc_addvar(filename,struct('Name','pacific_psi','Datatype','float','Dimension',{{'T','Zpsi','Ypsi'}}))
nc_addvar(filename,struct('Name','atlantic_eul','Datatype','float','Dimension',{{'T','Zpsi','Ypsi'}}));
nc_addvar(filename,struct('Name','pacific_eul','Datatype','float','Dimension',{{'T','Zpsi','Ypsi'}}))
nc_addvar(filename,struct('Name','global_eddy_psi','Datatype','float','Dimension',{{'T','Zpsi','Ypsi'}}))
nc_addvar(filename,struct('Name','atlantic_eddy_psi','Datatype','float','Dimension',{{'T','Zpsi','Ypsi'}}));
nc_addvar(filename,struct('Name','pacific_eddy_psi','Datatype','float','Dimension',{{'T','Zpsi','Ypsi'}}))
nc_addvar(filename,struct('Name','gbaro_mask','Datatype','float','Dimension',{{'Ypsi','Xpsi'}}))
nc_addvar(filename,struct('Name','got_mask','Datatype','float','Dimension',{{'Zpsi','Ypsi'}}))
nc_addvar(filename,struct('Name','aot_mask','Datatype','float','Dimension',{{'Zpsi','Ypsi'}}))
nc_addvar(filename,struct('Name','pot_mask','Datatype','float','Dimension',{{'Zpsi','Xpsi'}}))

%nc_padheader(filename,20000);
% nc_attput(filename,nc_global,'title','Streamfunctions');
nc_attput(filename,'gbaro_psi','units','m3 s-1');
nc_attput(filename,'global_psi','units','m3 s-1');
nc_attput(filename,'atlantic_psi','units','m3 s-1');
nc_attput(filename,'pacific_psi','units','m3 s-1');
nc_attput(filename,'gbaro_eul','units','m3 s-1');
nc_attput(filename,'global_eul','units','m3 s-1');
nc_attput(filename,'atlantic_eul','units','m3 s-1');
nc_attput(filename,'pacific_eul','units','m3 s-1');
nc_attput(filename,'gbaro_eddy_psi','units','m3 s-1');
nc_attput(filename,'global_eddy_psi','units','m3 s-1');
nc_attput(filename,'atlantic_eddy_psi','units','m3 s-1');
nc_attput(filename,'pacific_eddy_psi','units','m3 s-1');
% nc_attput(filename,'gm_vbg','units','m s-1');
% nc_attput(filename,'gm_ubg','units','m s-1');
% nc_attput(filename,'gm_wbg','units','m s-1');
nc_attput(filename,'Xpsi','units','degrees');
nc_attput(filename,'Ypsi','units','degrees');
nc_attput(filename,'Xpsi','modulo','');
nc_attput(filename,'T','units',[tavesteps.tuname,'s'])

nc_attput(filename,'gbaro_psi','title','Barotropic Streamfunction (Eularian mean and Eddy Advection)');
nc_attput(filename,'global_psi','title','Global Overturning Streamfunction (Eularian mean and Eddy Advection)');
nc_attput(filename,'atlantic_psi','title','Atlantic Overturning Streamfunction (Eularian mean and Eddy Advection)');
nc_attput(filename,'pacific_psi','title','Pacific Overturning Streamfunction (Eularian mean and Eddy Advection)');
nc_attput(filename,'gbaro_eddy_psi','title','Barotropic Streamfunction (GM Eddy Advection only)');
nc_attput(filename,'global_eddy_psi','title','Global Overturning Streamfunction (GM Eddy Advection only)');
nc_attput(filename,'atlantic_eddy_psi','title','Atlantic Overturning Streamfunction (GM Eddy Advection only)');
nc_attput(filename,'pacific_eddy_psi','title','Pacific Overturning Streamfunction (GM Eddy Advection only)');
nc_attput(filename,'gbaro_eul','title','Barotropic Streamfunction (Eularian mean only)');
nc_attput(filename,'global_eul','title','Global Overturning Streamfunction (Eularian mean only)');
nc_attput(filename,'atlantic_eul','title','Atlantic Overturning Streamfunction (Eularian mean only)');
nc_attput(filename,'pacific_eul','title','Pacific Overturning Streamfunction (Eularian mean only)');
% nc_attput(filename,'gm_vbg','title','Meridional GM Bolus Velocity');
% nc_attput(filename,'gm_ubg','title','Zonal GM Bolus Velocity');
% nc_attput(filename,'gm_wbg','title','Vertical GM Bolus Velocity');

attnames=fieldnames(attributes);
for attnum = 1:length(attnames)
    comm=sprintf('nc_attput(filename,nc_global,''%s'',%s(attributes.(attnames{attnum})));',...
        attnames{attnum},...
        class(attributes.(attnames{attnum}))); % not sure if the class thing works or not...
    eval(comm)
%     disp(['Adding global attributes: ', attnames{attnum}])
end

[~,os]=system('uname');

if ~isempty(strfind(lower(os),'darwin'))
    % only do wait handles on OSX
    wait_handle=waitbar(0,'Writing PSI components to netcdf');
end

for k=1:tavesteps.kmax
    if ~isempty(strfind(lower(os),'darwin'))
        % only do wait handles on OSX
        waitbar(k/tavesteps.kmax,wait_handle)
    end
         
    psi_data=[];
    % Rearrange axes so they are in the right order...
    baro_temp = permute(baro_psi(:,:,k),[3,2,1]);
    global_temp = permute(global_psi(:,:,k),[3,2,1]);
    atlantic_temp = permute(atlantic_psi(:,:,k),[3,2,1]);
    pacific_temp = permute(pacific_psi(:,:,k),[3,2,1]);

    % Masks....
    baro_psi_mask=ones(gridtmp.ny,gridtmp.nx);
    baro_psi_mask(~isnan(baro_temp(1,:,:))) = 0; %=(-1e34);
    nc_varput(filename,'gbaro_mask',baro_psi_mask);

    global_psi_mask=ones(gridtmp.nz,gridtmp.ny);
    global_psi_mask(~isnan(global_temp(1,:,:))) = 0; %(-1e34);
    nc_varput(filename,'got_mask',global_psi_mask);

    atlantic_psi_mask=ones(gridtmp.nz,gridtmp.ny);
    atlantic_psi_mask(~isnan(atlantic_temp(1,:,:))) = 0; %(-1e34);
    nc_varput(filename,'aot_mask',atlantic_psi_mask);

    pacific_psi_mask=ones(gridtmp.nz,gridtmp.ny);
    pacific_psi_mask(~isnan(pacific_temp(1,:,:))) = 0; %(-1e34);
    nc_varput(filename,'pot_mask',pacific_psi_mask);

    baroe_temp = permute(baro_eul(:,:,k),[3,2,1]);
    globale_temp = permute(global_eul(:,:,k),[3,2,1]);
    atlantice_temp = permute(atlantic_eul(:,:,k),[3,2,1]);
    pacifice_temp = permute(pacific_eul(:,:,k),[3,2,1]);
%     vbg_temp = permute(vbg(:,:,:,k),[4,3,2,1]);
%     ubg_temp = permute(ubg(:,:,:,k),[4,3,2,1]);
%     wbg_temp = permute(wbg(:,:,:,k),[4,3,2,1]);

    baroedds_temp = permute(baro_eddys(:,:,k),[3,2,1]);
    globaledds_temp = permute(global_eddys(:,:,k),[3,2,1]);
    atlanticedds_temp = permute(atlantic_eddys(:,:,k),[3,2,1]);
    pacificedds_temp = permute(pacific_eddys(:,:,k),[3,2,1]);
    
    % Data.....-1e34 is ferret's standard "bad" value, it doesnt seem to understand NaN
    baro_temp(isnan(baro_temp))=(-1e34);
    global_temp(isnan(global_temp))=(-1e34);
    atlantic_temp(isnan(atlantic_temp))=(-1e34);
    pacific_temp(isnan(pacific_temp))=(-1e34);

    baroe_temp(isnan(baroe_temp))=(-1e34);
    globale_temp(isnan(globale_temp))=(-1e34);
    atlantice_temp(isnan(atlantice_temp))=(-1e34);
    pacifice_temp(isnan(pacifice_temp))=(-1e34);
%     vbg_temp(isnan(vbg_temp))=(-1e34);
%     ubg_temp(isnan(ubg_temp))=(-1e34);
%     wbg_temp(isnan(wbg_temp))=(-1e34);
    
    baroedds_temp(isnan(baroedds_temp))=(-1e34);
    globaledds_temp(isnan(globaledds_temp))=(-1e34);
    atlanticedds_temp(isnan(atlanticedds_temp))=(-1e34);
    pacificedds_temp(isnan(pacificedds_temp))=(-1e34);

    psi_data.gbaro_psi = squeeze(baro_temp);%clear baro_temp
    psi_data.global_psi = squeeze(global_temp);%clear global_temp
    psi_data.atlantic_psi = squeeze(atlantic_temp);%clear atlantic_temp
    psi_data.pacific_psi = squeeze(pacific_temp);%clear pacific_temp
    psi_data.gbaro_eul = squeeze(baroe_temp);%clear baroe_temp
    psi_data.global_eul = squeeze(globale_temp);%clear globale_temp
    psi_data.atlantic_eul = squeeze(atlantice_temp);%clear atlantice_temp
    psi_data.pacific_eul = squeeze(pacifice_temp);%clear pacifice_temp
    psi_data.gbaro_eddy_psi = squeeze(baroedds_temp);%clear baro_temp
    psi_data.global_eddy_psi = squeeze(globaledds_temp);%clear global_temp
    psi_data.atlantic_eddy_psi = squeeze(atlanticedds_temp);%clear atlantic_temp
    psi_data.pacific_eddy_psi = squeeze(pacificedds_temp);%clear pacific_temp
%     psi_data.gm_vbg = squeeze(vbg_temp);%clear vbg_temp
%     psi_data.gm_ubg = squeeze(ubg_temp);%clear ubg_temp
%     psi_data.gm_wbg = squeeze(wbg_temp);%clear ubg_temp
    psi_data.iter = tavesteps.timesteps(k);
    psi_data.T = tavesteps.tim(k);
    nc_addrecs(filename,psi_data)
end

% Add floating_point missing_value to all variables
eval(['!ncatted -O -a missing_value,,c,f,-1.0e34 ',filename])

if ~isempty(strfind(lower(os),'darwin'))
    % only do wait handles on OSX
    close(wait_handle)
end

% nc=netcdf(filename,'write');
% nc.the_run_name=attributes.global.the_run_name;
% nc.MITgcm_version=attributes.global.MITgcm_version;
% nc.build_user=attributes.global.build_user;
% nc.build_host=attributes.global.build_host;
% nc.build_date=attributes.global.build_date;
% nc.MITgcm_tag_id=attributes.global.MITgcm_tag_id;
% nc.MITgcm_mnc_ver=attributes.global.MITgcm_mnc_ver;
% nc.tile_number='global';
% nc.sNx=attributes.global.sNx;
% nc.sNy=attributes.global.sNy;
% nc.OLx=attributes.global.OLx;
% nc.OLy=attributes.global.OLy;
% nc.nSx=attributes.global.nSx;
% nc.nSy=attributes.global.nSy;
% nc.nPx=attributes.global.nPx;
% nc.nPy=attributes.global.nPy;
% nc.Nx=attributes.global.Nx;
% nc.Ny=attributes.global.Ny;
% nc.Nr=attributes.global.Nr;
% close(nc);

%%clear global_psi atlantic_psi baro_psi pacific_psi ...
%    global_psi_mask atlantic_psi_mask baro_psi_mask pacific_psi_mask...
%    global_temp atlantic_temp baro_temp pacific_temp ...
%    atlosf pacosf globosf glob_baro_psi Mask psi_data
