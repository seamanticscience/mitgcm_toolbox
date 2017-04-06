function [filename]=mit_calc_psi(grd,tavesteps,filename)

if nargin==0 || ~exist('grd','var')
    grd=mit_loadgrid;
end

if nargin==1 || ~exist('tavesteps','var')
    tavesteps=mit_timesteps('tave');
end

if nargin==2 || ~exist('filename','var')
    filename='psi_data.nc';
end

if ~isfield(grd,'pacific_hfacw')
    grd = mit_oceanmasks(grd);
end

nc_create_empty(filename)
nc_adddim(filename,'Xpsi',grd.nx)
nc_adddim(filename,'Ypsi',grd.ny)
nc_adddim(filename,'Zpsi',length(grd.zgpsi))
nc_adddim(filename,'T',0)

nc_addvar(filename,struct('Name','Xpsi','Datatype','double','Dimension',{{'Xpsi'}}))
nc_addvar(filename,struct('Name','Ypsi','Datatype','double','Dimension',{{'Ypsi'}}))
nc_addvar(filename,struct('Name','Zpsi','Datatype','double','Dimension',{{'Zpsi'}}))
nc_addvar(filename,struct('Name','T','Datatype','double','Dimension',{{'T'}}))

% Write to NetCDF file
nc_varput(filename,'Xpsi',grd.long);
nc_varput(filename,'Ypsi',grd.latg);
nc_varput(filename,'Zpsi',grd.zgpsi);

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

% Extract attributes and axis for netcdf creation later....
if ~exist('attributes','var')
    tave = rdmnc(tavesteps.filearr(2:end-1),'Xp1');
    attributes=tave.attributes.global;
end

% Write Global Attributes to file
attnames=fieldnames(attributes);
%attributes.tile_number='global'; % DONT DO THIS!
for attnum = 1:length(attnames)
    comm=sprintf('nc_attput(filename,nc_global,''%s'',%s(attributes.(attnames{attnum})));',...
        attnames{attnum},...
        class(attributes.(attnames{attnum}))); % not sure if the class thing works or not...
    eval(comm)
    %disp(['Adding global attributes: ', comm])
end

[~,os]=system('uname');

if ~isempty(strfind(lower(os),'darwin'))
    % only do wait handles on OSX
    wait_handle=waitbar(0,'Calculating circulation streamfunctions');
end

for k=1:tavesteps.kmax
    if ~isempty(strfind(lower(os),'darwin'))
        % only do wait handles on OSX
        waitbar(k/tavesteps.kmax,wait_handle)
    end
    
    % Load residual velocity diagnostics or calculate from output
    if ~exist(strrep(tavesteps.filearr(2:end-1),'tave','gmDiag'),'file') ...
           && exist(strrep(tavesteps.filearr(2:end-1),'tave','oceDiag'),'file') 
       
       if k==1
           diagsteps=mit_timesteps(mit_getparm('data.diagnostics','filename'));
       end

        if  nc_isvar(strrep(diagsteps.filearr(2:end-1),'surf','gm'),'GM_U_RES') ...
                && nc_isvar(strrep(diagsteps.filearr(2:end-1),'surf','gm'),'GM_V_RES') ...
                && nc_isvar(strrep(diagsteps.filearr(2:end-1),'surf','gm'),'GM_W_RES') ...
                && nc_isvar(strrep(diagsteps.filearr(2:end-1),'surf','oce'),'UVEL') ...
                && nc_isvar(strrep(diagsteps.filearr(2:end-1),'surf','oce'),'VVEL') ...
                && nc_isvar(strrep(diagsteps.filearr(2:end-1),'surf','oce'),'WVEL')
            
            if k==1
                disp(['Loading eddy velocities from ',strrep(diagsteps.filearr(2:end-1),'surf','gm')])
            end
            
            gmdiag=rdmnc(strrep(diagsteps.filearr(2:end-1),'surf','gm'),...
                'GM_V_RES','GM_U_RES','GM_W_RES','GM_U_EDD','GM_V_EDD','GM_W_EDD',...
                diagsteps.timesteps(k));
            
            ocediag=rdmnc(strrep(diagsteps.filearr(2:end-1),'surf','oce'),...
                'UVEL','VVEL','WVEL',tavesteps.timesteps(k));
            
            uvelk=ocediag.UVEL;
            vvelk=ocediag.VVEL;
            wvelk=ocediag.WVEL;
            
            ueddk=gmdiag.GM_U_EDD;
            veddk=gmdiag.GM_V_EDD;
            weddk=gmdiag.GM_W_EDD;
            
            uresk=gmdiag.GM_U_RES;
            vresk=gmdiag.GM_V_RES;
            wresk=gmdiag.GM_W_RES;
        end
    end
    
    % Load the regular data if those arrays have not been filled with the
    % better diagnostics...
    if ~exist('veddk','var') ...
            || ~exist('vresk','var') ...
            || ~exist('ueddk','var') ...
            || ~exist('uresk','var')
               
        tave = rdmnc(tavesteps.filearr(2:end-1),'uVeltave','vVeltave','wVeltave',tavesteps.timesteps(k));
        if ~strcmp(grd.buoyancy,'OCEANIC');
            uvelk = tave.uVeltave(:,:,end:-1:1);
            vvelk = tave.vVeltave(:,:,end:-1:1);
            wvelk = tave.wVeltave(:,:,end:-1:1);
        else
            uvelk = tave.uVeltave;
            vvelk = tave.vVeltave;
            wvelk = tave.wVeltave;
        end
        
        % Calculate GM bolus velocities
        if nc_isvar(strrep(tavesteps.filearr(2:end-1),'tave','gm_tave'),'PsiX')...
                && nc_isvar(strrep(tavesteps.filearr(2:end-1),'tave','gm_tave'),'PsiY')
            if k==1
                disp(['Loading eddy PsiX and PsiY from ',strrep(tavesteps.filearr(2:end-1),'tave','gm_tave')])
            end
            gm = rdmnc(strrep(tavesteps.filearr(2:end-1),'tave','gm_tave'),'PsiX','PsiY',tavesteps.timesteps(k));
            [ueddk,veddk,weddk]=mit_gmvel(gm.PsiX,gm.PsiY,grd.nx,grd.ny,grd.nz,1,grd.dxg,grd.dyg,grd.dz,grd.cmask,grd.umask,grd.vmask,grd.rac);            
        else
            if k==1
                disp(['Loading eddy diffusivities from ',strrep(tavesteps.filearr(2:end-1),'tave','gm_tave')])
            end
            gm = rdmnc(strrep(tavesteps.filearr(2:end-1),'tave','gm_tave'),'Kwx','Kwy','Kwz',tavesteps.timesteps(k));
            [ueddk,veddk,weddk]=mit_gmvel(gm.Kwx,gm.Kwy,grd.nx,grd.ny,grd.nz,1,grd.dxg,grd.dyg,grd.dz,grd.cmask,grd.umask,grd.vmask,grd.rac);
        end
        
        uresk=uvelk+ueddk;
        vresk=vvelk+veddk;
        wresk=wvelk+weddk;
    end
    
    % Integrate the velocities to produce streamfunctions.
    addlayer = 1;
    
    % global overturning
    global_eul = mit_overturning(vvelk,grd.hfacs,grd.dxg,grd.dz,addlayer);
    global_psi = mit_overturning(vresk,grd.hfacs,grd.dxg,grd.dz,addlayer);
    global_eddys = mit_overturning(veddk,grd.hfacs,grd.dxg,grd.dz,addlayer);
    
    % atlantic overturning
    atlantic_eul = mit_overturning(vvelk,grd.atlantic_hfacs,grd.dxg,grd.dz,addlayer);
    atlantic_psi = mit_overturning(vresk,grd.atlantic_hfacs,grd.dxg,grd.dz,addlayer);
    atlantic_eddys = mit_overturning(veddk,grd.atlantic_hfacs,grd.dxg,grd.dz,addlayer);
    
    % pacific overturning
    pacific_eul = mit_overturning(vvelk,grd.pacific_hfacs,grd.dxg,grd.dz,addlayer);
    pacific_psi = mit_overturning(vresk,grd.pacific_hfacs,grd.dxg,grd.dz,addlayer);
    pacific_eddys = mit_overturning(veddk,grd.pacific_hfacs,grd.dxg,grd.dz,addlayer);
    
    % global barotropic stream function
    baro_eul = mit_barostream(uvelk,grd.umask,grd.dyg,grd.dz);
    baro_psi = mit_barostream(uresk,grd.umask,grd.dyg,grd.dz);
    baro_eddys = mit_barostream(ueddk,grd.umask,grd.dyg,grd.dz);
    
    psi_data=[];
    psi_data.gbaro_psi = change(baro_psi','==',NaN,-1e34);
    psi_data.global_psi = change(global_psi','==',NaN,-1e34);
    psi_data.atlantic_psi = change(atlantic_psi','==',NaN,-1e34);
    psi_data.pacific_psi = change(pacific_psi','==',NaN,-1e34);
    psi_data.gbaro_eul = change(baro_eul','==',NaN,-1e34);
    psi_data.global_eul = change(global_eul','==',NaN,-1e34);
    psi_data.atlantic_eul = change(atlantic_eul','==',NaN,-1e34);
    psi_data.pacific_eul = change(pacific_eul','==',NaN,-1e34);
    psi_data.gbaro_eddy_psi = change(baro_eddys','==',NaN,-1e34);
    psi_data.global_eddy_psi = change(global_eddys','==',NaN,-1e34);
    psi_data.atlantic_eddy_psi = change(atlantic_eddys','==',NaN,-1e34);
    psi_data.pacific_eddy_psi = change(pacific_eddys','==',NaN,-1e34);
    psi_data.iter = tavesteps.timesteps(k);
    psi_data.T = tavesteps.tim(k);
    nc_addrecs(filename,psi_data)
    clear uvelk vvelk wvelk ueddk veddk weddk uresk vresk wresk
end

    % Add Masks....
    baro_psi_mask=ones(grd.ny,grd.nx);
    baro_psi_mask(~isnan(baro_psi')) = 0; %=(-1e34);
    nc_varput(filename,'gbaro_mask',baro_psi_mask);
    
    global_psi_mask=ones(length(grd.zgpsi),grd.ny);
    global_psi_mask(~isnan(global_psi')) = 0; %(-1e34);
    nc_varput(filename,'got_mask',global_psi_mask);
    
    atlantic_psi_mask=ones(length(grd.zgpsi),grd.ny);
    atlantic_psi_mask(~isnan(atlantic_psi')) = 0; %(-1e34);
    nc_varput(filename,'aot_mask',atlantic_psi_mask);
    
    pacific_psi_mask=ones(length(grd.zgpsi),grd.ny);
    pacific_psi_mask(~isnan(pacific_psi')) = 0; %(-1e34);
    nc_varput(filename,'pot_mask',pacific_psi_mask);

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
