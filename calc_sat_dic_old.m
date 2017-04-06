function [filename]=calc_sat_dic(grd,tavesteps,filename)
%% Load variables
if nargin==0 || ~exist('grd','var')
    grd=mit_loadgrid;
end

if nargin==1 || ~exist('tavesteps','var')
    tavesteps=mit_timesteps('tave');
end

if nargin==2 || ~exist('filename','var')
    filename='dic_sat.glob.nc';
end
% Change units from the input of mol/m^3 -> mol/kg:
% (1 mol/m^3) x (1 /1024.5 kg m^-3)
permil=1/1024.5 ;
% convert mol/m3 to umol/kg
permumolkg=permil.*1e6;

[atm_co2]=mit_getdicpco2('data',tavesteps);

atm_co2(:,3)=atm_co2(:,3).*1e6; % Convert from atm to uatm

% Silicate is non-prognostic, so load full depth values from hydrography
% and regrid onto model domain
if exist('~/Dropbox_Work/Applications/ferret/woa01_siz.nc','file')
    sifname='~/Dropbox_Work/Applications/ferret/woa01_siz.nc';
elseif exist('/noc/users/jml102/Applications/ferret/woa01_siz.nc','file')
    sifname='/noc/users/jml102/Applications/ferret/woa01_siz.nc';
elseif exist('/home/jml1/Applications/Data/woa01_siz.nc','file')
    sifname='/home/jml1/Applications/Data/woa01_siz.nc';
else
    error('Cannot find Si climatology')
end

sitmp=nc_varget(sifname,'Silicate');
LATITUDE=nc_varget(sifname,'LATITUDE');
LONGITUDE=nc_varget(sifname,'LONGITUDE');
DEPTH=nc_varget(sifname,'DEPTH');

sitmp=permute(sitmp,[3 2 1]);
[lat,long,depth]=meshgrid(LATITUDE,LONGITUDE,DEPTH);
[xc,yc,zc]=meshgrid(grd.lonc,grd.latc,(-grd.zc));
si=interp3(lat,long,depth,inpaint_nans(sitmp,4),yc,xc,zc);
si=permute(si,[2 1 3]); % zc=permute(zc,[2 1 3]);
si=(si./1000).*grd.hfacc; % woa01 si in umol/l convert to mol/m3

% Preallocate arrays
dic_satini=nan(grd.nx,grd.ny,grd.nz);
ph_satini=nan(grd.nx,grd.ny,grd.nz);
revelle_satini=nan(grd.nx,grd.ny,grd.nz);
omegac_satini=nan(grd.nx,grd.ny,grd.nz);
omegaa_satini=nan(grd.nx,grd.ny,grd.nz);
dic_satatm=nan(grd.nx,grd.ny,grd.nz);
ph_satatm=nan(grd.nx,grd.ny,grd.nz);
revelle_satatm=nan(grd.nx,grd.ny,grd.nz);
omegac_satatm=nan(grd.nx,grd.ny,grd.nz);
omegaa_satatm=nan(grd.nx,grd.ny,grd.nz);
dic_sattot=nan(grd.nx,grd.ny,grd.nz);
ph_sattot=nan(grd.nx,grd.ny,grd.nz);
revelle_sattot=nan(grd.nx,grd.ny,grd.nz);
omegac_sattot=nan(grd.nx,grd.ny,grd.nz);
omegaa_sattot=nan(grd.nx,grd.ny,grd.nz);
dic_satalk=nan(grd.nx,grd.ny,grd.nz);
ph_satalk=nan(grd.nx,grd.ny,grd.nz);
revelle_satalk=nan(grd.nx,grd.ny,grd.nz);
omegac_satalk=nan(grd.nx,grd.ny,grd.nz);
omegaa_satalk=nan(grd.nx,grd.ny,grd.nz);
dic_satt=nan(grd.nx,grd.ny,grd.nz);
ph_satt=nan(grd.nx,grd.ny,grd.nz);
revelle_satt=nan(grd.nx,grd.ny,grd.nz);
omegac_satt=nan(grd.nx,grd.ny,grd.nz);
omegaa_satt=nan(grd.nx,grd.ny,grd.nz);
dic_sats=nan(grd.nx,grd.ny,grd.nz);
ph_sats=nan(grd.nx,grd.ny,grd.nz);
revelle_sats=nan(grd.nx,grd.ny,grd.nz);
omegac_sats=nan(grd.nx,grd.ny,grd.nz);
omegaa_sats=nan(grd.nx,grd.ny,grd.nz);

%% Set up Output to NetCDF
% Write out to netcdf file...
%filename='dic_sat.test.glob.nc'; % assign as an input
nc_create_empty(filename,'clobber')
%nc_padheader(filename,20000);
nc_adddim(filename,'X',grd.nx);
nc_adddim(filename,'Y',grd.ny);
nc_adddim(filename,'Z',grd.nz);
nc_adddim(filename,'T',0);

nc_addvar(filename,struct('Name','X','Datatype','double','Dimension',{{'X'}}))
nc_addvar(filename,struct('Name','Y','Datatype','double','Dimension',{{'Y'}}))
nc_addvar(filename,struct('Name','Z','Datatype','double','Dimension',{{'Z'}}))
nc_addvar(filename,struct('Name','T','Datatype','double','Dimension',{{'T'}}))

nc_attput(filename,'X','units','degrees')
nc_attput(filename,'Y','units','degrees')
nc_attput(filename,'Z','units','m')
nc_attput(filename,'T','units',[tavesteps.tuname,'s'])

% Write to NetCDF file
nc_varput(filename,'X',grd.lonc);
nc_varput(filename,'Y',grd.latc);
nc_varput(filename,'Z',-grd.zc);

nc_addvar(filename,struct('Name','iter','Datatype','int','Dimension',{{'T'}}));
nc_addvar(filename,struct('Name','alk_pref','Datatype','float','Dimension',{{'T','Z','Y','X'}}));
nc_attput(filename,'alk_pref','units','mol/m3');
nc_attput(filename,'alk_pref','title','Preformed Alkalinity')
% Not enough precision, should appear as -1.e+34f in ncdumpexi
% nc_attput(filename,'alk_pref','missing_value',-1.0e34)

% wrt control pco2, this is dT(dCsat/dT)+dALK(dCsat/dALK)+dS(dCsat/dS)+dNUTS(dCsat/dNUTS)
nc_addvar(filename,struct('Name','dic_satini','Datatype','float','Dimension',{{'T','Z','Y','X'}}));
nc_addvar(filename,struct('Name','ph_satini','Datatype','float','Dimension',{{'T','Z','Y','X'}}));
nc_addvar(filename,struct('Name','revelle_satini','Datatype','float','Dimension',{{'T','Z','Y','X'}}));
nc_addvar(filename,struct('Name','omegac_satini','Datatype','float','Dimension',{{'T','Z','Y','X'}}));
nc_addvar(filename,struct('Name','omegaa_satini','Datatype','float','Dimension',{{'T','Z','Y','X'}}));

nc_attput(filename,'dic_satini','title','Saturated DIC concentration at control pCO2 (278.05)')
nc_attput(filename,'ph_satini','title','pH at control pCO2 (278.05)')
nc_attput(filename,'revelle_satini','title','Revelle Factor at control pCO2 (278.05)')
nc_attput(filename,'omegac_satini','title','Solubility Ratio for Calcite at control pCO2 (278.05)')
nc_attput(filename,'omegaa_satini','title','Solubility Ratio of Aragonite at control pCO2 (278.05)')

nc_attput(filename,'dic_satini','units','mol/kg')
nc_attput(filename,'ph_satini','units','')
nc_attput(filename,'revelle_satini','units','')
nc_attput(filename,'omegac_satini','units','')
nc_attput(filename,'omegaa_satini','units','')

% nc_attput(filename,'dic_satini','missing_value',-1.0e34)
% nc_attput(filename,'ph_satini','missing_value',-1.0e34)
% nc_attput(filename,'revelle_satini','missing_value',-1.0e34)
% nc_attput(filename,'omegac_satini','missing_value',-1.0e34)
% nc_attput(filename,'omegaa_satini','missing_value',-1.0e34)

% wrt prognostic pco2 only, this is dCO2(dCsat/dCO2)
nc_addvar(filename,struct('Name','dic_satatm','Datatype','float','Dimension',{{'T','Z','Y','X'}}));
nc_addvar(filename,struct('Name','ph_satatm','Datatype','float','Dimension',{{'T','Z','Y','X'}}));
nc_addvar(filename,struct('Name','revelle_satatm','Datatype','float','Dimension',{{'T','Z','Y','X'}}));
nc_addvar(filename,struct('Name','omegac_satatm','Datatype','float','Dimension',{{'T','Z','Y','X'}}));
nc_addvar(filename,struct('Name','omegaa_satatm','Datatype','float','Dimension',{{'T','Z','Y','X'}}));

nc_attput(filename,'dic_satatm','title','Saturated DIC concentration at prognostic atmospheric pCO2')
nc_attput(filename,'ph_satatm','title','pH at prognostic atmospheric pCO2')
nc_attput(filename,'revelle_satatm','title','Revelle Factor at prognostic atmospheric pCO2')
nc_attput(filename,'omegac_satatm','title','Solubility Ratio for Calcite at prognostic atmospheric pCO2')
nc_attput(filename,'omegaa_satatm','title','Solubility Ratio of Aragonite at prognostic atmospheric pCO2')

nc_attput(filename,'dic_satatm','units','mol/kg')
nc_attput(filename,'ph_satatm','units','')
nc_attput(filename,'revelle_satatm','units','')
nc_attput(filename,'omegac_satatm','units','')
nc_attput(filename,'omegaa_satatm','units','')

% nc_attput(filename,'dic_satatm','missing_value',-1.0e34)
% nc_attput(filename,'ph_satatm','missing_value',-1.0e34)
% nc_attput(filename,'revelle_satatm','missing_value',-1.0e34)
% nc_attput(filename,'omegac_satatm','missing_value',-1.0e34)
% nc_attput(filename,'omegaa_satatm','missing_value',-1.0e34)

% wrt prognostic pco2 only, this is dCO2(dCsat/dCO2)+dT(dCsat/dT)+dALK(dCsat/dALK)+dS(dCsat/dS)
nc_addvar(filename,struct('Name','dic_sattot','Datatype','float','Dimension',{{'T','Z','Y','X'}}));
nc_addvar(filename,struct('Name','ph_sattot','Datatype','float','Dimension',{{'T','Z','Y','X'}}));
nc_addvar(filename,struct('Name','revelle_sattot','Datatype','float','Dimension',{{'T','Z','Y','X'}}));
nc_addvar(filename,struct('Name','omegac_sattot','Datatype','float','Dimension',{{'T','Z','Y','X'}}));
nc_addvar(filename,struct('Name','omegaa_sattot','Datatype','float','Dimension',{{'T','Z','Y','X'}}));

nc_attput(filename,'dic_sattot','title','Saturated DIC concentration')
nc_attput(filename,'ph_sattot','title','pH')
nc_attput(filename,'revelle_sattot','title','Revelle Factor')
nc_attput(filename,'omegac_sattot','title','Solubility Ratio for Calcite')
nc_attput(filename,'omegaa_sattot','title','Solubility Ratio of Aragonite')

nc_attput(filename,'dic_sattot','units','mol/kg')
nc_attput(filename,'ph_sattot','units','')
nc_attput(filename,'revelle_sattot','units','')
nc_attput(filename,'omegac_sattot','units','')
nc_attput(filename,'omegaa_sattot','units','')

% nc_attput(filename,'dic_sattot','missing_value',-1.0e34)
% nc_attput(filename,'ph_sattot','missing_value',-1.0e34)
% nc_attput(filename,'revelle_sattot','missing_value',-1.0e34)
% nc_attput(filename,'omegac_sattot','missing_value',-1.0e34)
% nc_attput(filename,'omegaa_sattot','missing_value',-1.0e34)

% wrt alkalinity changes only
nc_addvar(filename,struct('Name','dic_satalk','Datatype','float','Dimension',{{'T','Z','Y','X'}}));
nc_addvar(filename,struct('Name','ph_satalk','Datatype','float','Dimension',{{'T','Z','Y','X'}}));
nc_addvar(filename,struct('Name','revelle_satalk','Datatype','float','Dimension',{{'T','Z','Y','X'}}));
nc_addvar(filename,struct('Name','omegac_satalk','Datatype','float','Dimension',{{'T','Z','Y','X'}}));
nc_addvar(filename,struct('Name','omegaa_satalk','Datatype','float','Dimension',{{'T','Z','Y','X'}}));

nc_attput(filename,'dic_satalk','title','Saturated DIC concentration due to ALK_PREF')
nc_attput(filename,'ph_satalk','title','pH due to ALK_PREF')
nc_attput(filename,'revelle_satalk','title','Revelle Factor due to ALK_PREF')
nc_attput(filename,'omegac_satalk','title','Solubility Ratio for Calcite due to ALK_PREF')
nc_attput(filename,'omegaa_satalk','title','Solubility Ratio of Aragonite due to ALK_PREF')

nc_attput(filename,'dic_satalk','units','mol/kg')
nc_attput(filename,'ph_satalk','units','')
nc_attput(filename,'revelle_satalk','units','')
nc_attput(filename,'omegac_satalk','units','')
nc_attput(filename,'omegaa_satalk','units','')

% nc_attput(filename,'dic_satalk','missing_value',-1.0e34)
% nc_attput(filename,'ph_satalk','missing_value',-1.0e34)
% nc_attput(filename,'revelle_satalk','missing_value',-1.0e34)
% nc_attput(filename,'omegac_satalk','missing_value',-1.0e34)
% nc_attput(filename,'omegaa_satalk','missing_value',-1.0e34)

% wrt T changes only
nc_addvar(filename,struct('Name','dic_satt','Datatype','float','Dimension',{{'T','Z','Y','X'}}));
nc_addvar(filename,struct('Name','ph_satt','Datatype','float','Dimension',{{'T','Z','Y','X'}}));
nc_addvar(filename,struct('Name','revelle_satt','Datatype','float','Dimension',{{'T','Z','Y','X'}}));
nc_addvar(filename,struct('Name','omegac_satt','Datatype','float','Dimension',{{'T','Z','Y','X'}}));
nc_addvar(filename,struct('Name','omegaa_satt','Datatype','float','Dimension',{{'T','Z','Y','X'}}));

nc_attput(filename,'dic_satt','title','Saturated DIC concentration due to T')
nc_attput(filename,'ph_satt','title','pH due to T')
nc_attput(filename,'revelle_satt','title','Revelle Factor due to T')
nc_attput(filename,'omegac_satt','title','Solubility Ratio for Calcite due to T')
nc_attput(filename,'omegaa_satt','title','Solubility Ratio of Aragonite due to T')

nc_attput(filename,'dic_satt','units','mol/kg')
nc_attput(filename,'ph_satt','units','')
nc_attput(filename,'revelle_satt','units','')
nc_attput(filename,'omegac_satt','units','')
nc_attput(filename,'omegaa_satt','units','')

% nc_attput(filename,'dic_satt','missing_value',-1.0e34)
% nc_attput(filename,'ph_satt','missing_value',-1.0e34)
% nc_attput(filename,'revelle_satt','missing_value',-1.0e34)
% nc_attput(filename,'omegac_satt','missing_value',-1.0e34)
% nc_attput(filename,'omegaa_satt','missing_value',-1.0e34)

% wrt S changes only
nc_addvar(filename,struct('Name','dic_sats','Datatype','float','Dimension',{{'T','Z','Y','X'}}));
nc_addvar(filename,struct('Name','ph_sats','Datatype','float','Dimension',{{'T','Z','Y','X'}}));
nc_addvar(filename,struct('Name','revelle_sats','Datatype','float','Dimension',{{'T','Z','Y','X'}}));
nc_addvar(filename,struct('Name','omegac_sats','Datatype','float','Dimension',{{'T','Z','Y','X'}}));
nc_addvar(filename,struct('Name','omegaa_sats','Datatype','float','Dimension',{{'T','Z','Y','X'}}));

nc_attput(filename,'dic_sats','title','Saturated DIC concentration due to S')
nc_attput(filename,'ph_sats','title','pH due to S')
nc_attput(filename,'revelle_sats','title','Revelle Factor due to S')
nc_attput(filename,'omegac_sats','title','Solubility Ratio for Calcite due to S')
nc_attput(filename,'omegaa_sats','title','Solubility Ratio of Aragonite due to S')

nc_attput(filename,'dic_sats','units','mol/kg')
nc_attput(filename,'ph_sats','units','')
nc_attput(filename,'revelle_sats','units','')
nc_attput(filename,'omegac_sats','units','')
nc_attput(filename,'omegaa_sats','units','')

% nc_attput(filename,'dic_sats','missing_value',-1.0e34)
% nc_attput(filename,'ph_sats','missing_value',-1.0e34)
% nc_attput(filename,'revelle_sats','missing_value',-1.0e34)
% nc_attput(filename,'omegac_sats','missing_value',-1.0e34)
% nc_attput(filename,'omegaa_sats','missing_value',-1.0e34)

% Extract attributes and axis for netcdf creation later....
if ~exist('attributes','var')
    tave = rdmnc(tavesteps.filearr(2:end-1),'Xp1');
    attributes=tave.attributes.global;
end

% Write Global Attributes to file
attnames=fieldnames(attributes);
% attributes.tile_number='global'; % DONT DO THIS
for attnum = 1:length(attnames)
    comm=sprintf('nc_attput(filename,nc_global,''%s'',%s(attributes.(attnames{attnum})));',...
        attnames{attnum},...
        class(attributes.(attnames{attnum}))); % not sure if the class thing works or not...
    eval(comm)
    %disp(['Adding global attributes: ', comm])
end

%% Load Control Data
[~,os]=system('uname');
if ~exist('tave.20k.glob.nc','file')
    if ~isempty(strfind(lower(os),'darwin'))
        % MacOSX
        if exist('/Volumes/PhD_Data/','dir')
            filepath='/Volumes/PhD_Data';
        elseif exist('/Volumes/Postdoc_Data/','dir')
            filepath='/Volumes/Postdoc_Data';
        else
            [~,mitpath]=system('echo $mitgcm');
           
            if ~isempty(mitpath)
                filepath=mitpath(1:end-1);
            else
                filepath='./';
            end
            warning('Could not find external drives, so using path: %s',mitpath)
        end
    elseif ~isempty(strfind(lower(os),'linux'))
        % Assume cluster
        [~,mitpath]=system('echo $mitgcm');
        filepath=mitpath(1:end-1);
    else
        error('Could not determine architecture, so cannot continue')
    end
    clear mitpath
% Load respective control runs. Incase
% these are actually the control runs for the individual parameterisations
% then the standard control run is loaded.
    if (~isempty(strfind(lower(pwd),'exf')) || exist('data.exf','file')) && ...
            isempty(strfind(lower(pwd),'cntrl'))
        % external forcing package runs have different starting state
        disp('Loading EXF control run')
        cntrl_tave=rdmnc([filepath,'/controlrun/cntrl_exf_global_1kyr/tave.32k.glob.nc'],'Ttave','Stave');
        cntrl_ptr=rdmnc([filepath,'/controlrun/cntrl_exf_global_1kyr/ptr_tave.32k.glob.nc'],'alk','dic','po4','o2');
    elseif ~isempty(strfind(lower(pwd),'vmhs97')) && ...
            isempty(strfind(lower(pwd),'cntrl'))
        disp('Loading vmhs97 eddy model control run')
        cntrl_tave=rdmnc([filepath,'/perturb_eddies/VMHS97/cntrl_vmhs97/cntrl_vmhs97_global/tave.30k.glob.nc'],'Ttave','Stave',21600000);
        cntrl_ptr=rdmnc([filepath,'/perturb_eddies/VMHS97/cntrl_vmhs97/cntrl_vmhs97_global/ptr_tave.30k.glob.nc'],'alk','dic','po4','o2',21600000);
    elseif ~isempty(strfind(lower(pwd),'hmm2011')) && ...
            isempty(strfind(lower(pwd),'cntrl'))
        disp('Loading hmm2011 eddy model control run')
        cntrl_tave=rdmnc([filepath,'/perturb_eddies/HMM2011/cntrl_hmm2011/cntrl_hmm2011_global/tave.30k.glob.nc'],'Ttave','Stave',21600000);
        cntrl_ptr=rdmnc([filepath,'/perturb_eddies/HMM2011/cntrl_hmm2011/cntrl_hmm2011_global/ptr_tave.30k.glob.nc'],'alk','dic','po4','o2',21600000);
    elseif ~isempty(strfind(lower(pwd),'fmh2005')) && ...
            isempty(strfind(lower(pwd),'cntrl'))
        disp('Loading fmh2005 model control run')
        cntrl_tave=rdmnc([filepath,'/perturb_eddies/FMH2005/cntrl_fmh2005/cntrl_fmh2005_global/tave.30k.glob.nc'],'Ttave','Stave',21600000);
        cntrl_ptr=rdmnc([filepath,'/perturb_eddies/FMH2005/cntrl_fmh2005/cntrl_fmh2005_global/ptr_tave.30k.glob.nc'],'alk','dic','po4','o2',21600000);
    elseif grd.nx==128 && grd.ny==64 && grd.nz==15
        disp('Loading standard model control run')
        cntrl_tave=rdmnc([filepath,'/controlrun/controlrun_anomplotref/tave.20k.glob.nc'],'Ttave','Stave');
        cntrl_ptr=rdmnc([filepath,'/controlrun/controlrun_anomplotref/ptr_tave.20k.glob.nc'],'alk','dic','po4','o2');
    else
        disp('Loading initial binary files')
        cntrl_tave.Ttave=mit_readfield(mit_getparm('data','hydrogThetaFile'),[(grd.nx),(grd.ny),(grd.nz)],'real*4');
        cntrl_tave.Stave=mit_readfield(mit_getparm('data','hydrogSaltFile'),[(grd.nx),(grd.ny),(grd.nz)],'real*4');
        
        %collect ptracer names into an array
        nptr=mit_getparm('data.ptracers','PTRACERS_numInUse');
        ptr_names=cell(1,nptr);
        
        for ip=1:nptr
            ptr_names{ip}=mit_getparm('data.ptracers',['PTRACERS_names(',num2str(ip),')']);
        end
        
        cntrl_ptr.dic=mit_readfield(mit_getparm('data.ptracers',...
            ['PTRACERS_initialFile(',num2str(find(strcmp('dic',ptr_names))),')']),...
            [(grd.nx),(grd.ny),(grd.nz)],'real*4');
        cntrl_ptr.alk=mit_readfield(mit_getparm('data.ptracers',...
            ['PTRACERS_initialFile(',num2str(find(strcmp('alk',ptr_names))),')']),...
            [(grd.nx),(grd.ny),(grd.nz)],'real*4');
        cntrl_ptr.po4=mit_readfield(mit_getparm('data.ptracers',...
            ['PTRACERS_initialFile(',num2str(find(strcmp('po4',ptr_names))),')']),...
            [(grd.nx),(grd.ny),(grd.nz)],'real*4');
        cntrl_ptr.o2 =mit_readfield(mit_getparm('data.ptracers',...
            ['PTRACERS_initialFile(',num2str(find(strcmp('o2',ptr_names))),')']),...
            [(grd.nx),(grd.ny),(grd.nz)],'real*4');
    end
else
    cntrl_tave=rdmnc('tave.20k.glob.nc','Ttave','Stave');
    cntrl_ptr=rdmnc('ptr_tave.20k.glob.nc','alk','dic','po4','o2');    
end

cntrl_s=cntrl_tave.Stave(:,:,:,end).*grd.hfacc;
cntrl_t=cntrl_tave.Ttave(:,:,:,end).*grd.hfacc;

% Use preformed Alkalinity (in umol/kg)
% Define PO=O2-R(O2:P)PO4, where R(O2:P)=-170
cntrl_po=(cntrl_ptr.o2(:,:,:,end)+(170*cntrl_ptr.po4(:,:,:,end))).*grd.hfacc;

% Use the correct controlpreformed alk - if exf control run, use standard
% relationship here, but exf relationship elsewhere.
if (~isempty(strfind(lower(pwd),'exf')) || exist('data.exf','file')) && (isempty(strfind(lower(pwd),'cntrl')))
    % use different relationship for preformed alkalinity using exf
    disp('Using alternate (EXF) relationship for control preformed alkalinity...')
    cntrl_ta=(0.1332+(0.06210.*cntrl_tave.Stave(:,:,:,end))+(0.06737.*cntrl_po(:,:,:,end))).*grd.hfacc.*permumolkg;
else
    disp('use standard (Relaxing) values for control preformed alkalinity...')
    cntrl_ta=(0.1114+(0.062426.*cntrl_tave.Stave(:,:,:,end))+(0.095896.*cntrl_po(:,:,:,end))).*grd.hfacc.*permumolkg;
end

% Use preformed phosphate (in umol/kg)
[cntrl_o2sat]=calc_oxygen_sat(cntrl_t,cntrl_s);
cntrl_aou = cntrl_o2sat-cntrl_ptr.o2(:,:,:,end) ;
cntrl_pbio = cntrl_aou/170 ;
cntrl_pt=(cntrl_ptr.po4(:,:,:,end)-cntrl_pbio).*grd.hfacc.*permumolkg;

% Estimate Preformed Si ROSI=-170/15 (in umol/kg)
% Is constant as not carried by the model - from WOA actually
cntrl_sibio=(cntrl_aou*15)/170;
cntrl_sit=(si-cntrl_sibio).*grd.hfacc.*permumolkg;

% % Do calculation bit in parallel :)
% if license('test','Distrib_Computing_Toolbox') ==1; % if license exists
%     check=license('checkout','Distrib_Computing_Toolbox'); % check to see if there is a spare copy
%     if check==1 && matlabpool('size')==0
%         %checked out uni license for parallel computing toolbox
%         %no existing matlab pool present
%         %therefore start matlabpool
%         matlabpool
%     end
% end

latgrd=repmat(grd.yc,[1,1,grd.nz]).*grd.hfacc;
zgrd=nan(size(grd.hfacc));
for i=1:grd.nz
    zgrd(:,:,i)=grd.hfacc(:,:,i).*grd.zc(i);
end

pdbar=sw_pres(zgrd,latgrd);
    

%% Calculate Carbon Chem Parameters
if ~isempty(strfind(lower(os),'darwin'))
    % only do wait handles on OSX
    wait_handle=waitbar(0,'Calculating DIC components');
end

for tstep=1:tavesteps.kmax;
    %    disp(['Processing timestep ',num2str(timesteps(tstep)),' (',num2str(tstep),')'])
    if ~isempty(strfind(lower(os),'darwin'))
        % only do wait handles on OSX
        waitbar(tstep/tavesteps.kmax,wait_handle)
    end
    
    tave=rdmnc(tavesteps.filearr(2:end-1),'Ttave','Stave',tavesteps.timesteps(tstep));
    ptracer=rdmnc(strrep(tavesteps.filearr(2:end-1),'tave','ptr_tave')...
        ,'alk','dic','po4','o2',tavesteps.timesteps(tstep));
    
    % This command prevents matlab from giving the error
    % "??? Undefined function or method 'ncregister' for input arguments of type 'double'.
    clear mex
    
    % The Follows et al (2006) method used PREFORMED ALKALINITY (A_PRE) - see
    % Williams and Follows (2011) Ch6. Since A_REG=0 at the surface,
    % A_T=A_PRE By definition. But since we are extending this to entire
    % domain we need to account for the carbonate pump. Luckily ALK dominated
    % by FW fluxes so has v.good relationship with Salinity. Use Multiple Linear
    % Regression against Salinity and PO (see Lauderdale thesis, 2010, Ch2)
    % Have checked that PO/Salt (R2=0.97223) gives better results than with
    % Salt alone (R2=0.93822).
 
    s=tave.Stave.*grd.hfacc;
    t=tave.Ttave.*grd.hfacc;
    
    % Use preformed Alkalinity (in umol/kg)
    % Define PO=O2-R(O2:P)PO4, where R(O2:P)=-170
    po=(ptracer.o2+(170*ptracer.po4)).*grd.hfacc;
    
    if ~isempty(strfind(lower(pwd),'exf')) || exist('data.exf','file')
        % use different relationship for preformed alkalinity using exf
        ta=(0.1332+(0.06210.*tave.Stave)+(0.06737.*po)).*grd.hfacc.*permumolkg;
    else
        % use standard values
        ta=(0.1114+(0.062426.*tave.Stave)+(0.095896.*po)).*grd.hfacc.*permumolkg;
    end
    
    % Following Nutrients are used to calculate Carbonate Alkalinity, which
    % is then used to get DIC conc (or Saturated DIC in this
    % case)..probably small effect
    % Use preformed phosphate (in umol/kg)
    [o2sat]=calc_oxygen_sat(t,s);
    aou = o2sat-ptracer.o2 ;
    pbio = aou/170 ;
    pt=(ptracer.po4-pbio).*grd.hfacc.*permumolkg;
    
    % Estimate Preformed Si ROSI=-170/15 (in umol/kg)
    sibio=(aou*15)/170;
    sit=(si-sibio).*grd.hfacc.*permumolkg;
    
    %     % Set up CO2SYS variables to calculate DIC given Preformed ALK and pCO2
    %     par1type =    4; % Data Type is pCO2 (uatm)
    %     par2type =    1; % The first parameter supplied is of type "1", which is "ALK"
    %     pHscale  =    1; % pH scale at which the input pH is reported ("1" means "Total Scale")  - doesn't matter in this example
    %     k1k2c    =    4; % Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("4" means "Mehrbach refit")
    %     kso4c    =    1; % Choice of HSO4- dissociation constants KSO4 ("1" means "Dickson")
    
    
    for j=1:1:grd.ny
        %        parfor k=1:grd.nz
        for k=1:1:grd.nz
            % Calculate Saturated DIC wrt prognostic atm_pco2 at each timestep
            % This is dCO2(dCsat/dCO2)+dT(dCsat/dT)+dALK(dCsat/dALK)+dS(dCsat/dS)
  
            A=CO2SYS(atm_co2(atm_co2(:,1)==tavesteps.tim(tstep),3),ta(:,j,k),4,1,s(:,j,k),t(:,j,k),t(:,j,k),...
                pdbar(:,j,k),pdbar(:,j,k),sit(:,j,k),pt(:,j,k),1,4,1);
            
            dic_sattot(:,j,k)=A(:,2)./permumolkg;
            ph_sattot(:,j,k)=A(:,3);
            revelle_sattot(:,j,k)=A(:,14);
            omegac_sattot(:,j,k)=A(:,15);
            omegaa_sattot(:,j,k)=A(:,16);
            
            % Calculate Saturated DIC wrt atm_pco2 = 278.05 (control)
            % This is dT(dCsat/dT)+dALK(dCsat/dALK)+dS(dCsat/dS)+dNUTS(dCsat/dNUTS)
            A=CO2SYS(atm_co2(1,3),ta(:,j,k),4,1,s(:,j,k),t(:,j,k),t(:,j,k),...
                pdbar(:,j,k),pdbar(:,j,k),sit(:,j,k),pt(:,j,k),1,4,1);
            
            dic_satini(:,j,k)=A(:,2)./permumolkg;
            ph_satini(:,j,k)=A(:,3);
            revelle_satini(:,j,k)=A(:,14);
            omegac_satini(:,j,k)=A(:,15);
            omegaa_satini(:,j,k)=A(:,16);
            
            %             if exist('cntrl_t','var')
            % Only calculate if control values available.
            
            % Calculate Saturated DIC with only pco2 changes.
            % This is dCO2(dCsat/dCO2)
            A=CO2SYS(atm_co2(atm_co2(:,1)==tavesteps.tim(tstep),3),cntrl_ta(:,j,k),4,1,cntrl_s(:,j,k),cntrl_t(:,j,k),cntrl_t(:,j,k),...
                pdbar(:,j,k),pdbar(:,j,k),cntrl_sit(:,j,k),cntrl_pt(:,j,k),1,4,1);
            
            dic_satatm(:,j,k)=A(:,2)./permumolkg;
            ph_satatm(:,j,k)=A(:,3);
            revelle_satatm(:,j,k)=A(:,14);
            omegac_satatm(:,j,k)=A(:,15);
            omegaa_satatm(:,j,k)=A(:,16);
            
            % Calculate Saturated DIC with only ALK.
            % This is dALK(dCsat/dALK)
            A=CO2SYS(atm_co2(1,3),ta(:,j,k),4,1,cntrl_s(:,j,k),cntrl_t(:,j,k),cntrl_t(:,j,k),...
                pdbar(:,j,k),pdbar(:,j,k),cntrl_sit(:,j,k),cntrl_pt(:,j,k),1,4,1);
            
            dic_satalk(:,j,k)=A(:,2)./permumolkg;
            ph_satalk(:,j,k)=A(:,3);
            revelle_satalk(:,j,k)=A(:,14);
            omegac_satalk(:,j,k)=A(:,15);
            omegaa_satalk(:,j,k)=A(:,16);
            
            % Calculate Saturated DIC with only T changes.
            % This is dT(dCsat/dT)
            A=CO2SYS(atm_co2(1,3),cntrl_ta(:,j,k),4,1,cntrl_s(:,j,k),t(:,j,k),t(:,j,k),...
                pdbar(:,j,k),pdbar(:,j,k),cntrl_sit(:,j,k),cntrl_pt(:,j,k),1,4,1);
            
            dic_satt(:,j,k)=A(:,2)./permumolkg;
            ph_satt(:,j,k)=A(:,3);
            revelle_satt(:,j,k)=A(:,14);
            omegac_satt(:,j,k)=A(:,15);
            omegaa_satt(:,j,k)=A(:,16);
            
            % Calculate Saturated DIC with only S changes.
            % This is dS(dCsat/dS)
            A=CO2SYS(atm_co2(1,3),cntrl_ta(:,j,k),4,1,s(:,j,k),cntrl_t(:,j,k),cntrl_t(:,j,k),...
                pdbar(:,j,k),pdbar(:,j,k),cntrl_sit(:,j,k),cntrl_pt(:,j,k),1,4,1);
            
            dic_sats(:,j,k)=A(:,2)./permumolkg;
            ph_sats(:,j,k)=A(:,3);
            revelle_sats(:,j,k)=A(:,14);
            omegac_sats(:,j,k)=A(:,15);
            omegaa_sats(:,j,k)=A(:,16);
        end
    end
    
    ta=ta./permumolkg;
    
    %% Output to NetCDF File
    % Rearrange axes so they are in the right order...
%     ta_temp = permute(ta,[3,2,1]);
%     
%     dic1_temp = permute(dic_satini,[3,2,1]);
%     ph1_temp = permute(ph_satini,[3,2,1]);
%     r1_temp = permute(revelle_satini,[3,2,1]);
%     oc1_temp = permute(omegac_satini,[3,2,1]);
%     oa1_temp = permute(omegaa_satini,[3,2,1]);
%     
%     dic2_temp = permute(dic_satatm,[3,2,1]);
%     ph2_temp = permute(ph_satatm,[3,2,1]);
%     r2_temp = permute(revelle_satatm,[3,2,1]);
%     oc2_temp = permute(omegac_satatm,[3,2,1]);
%     oa2_temp = permute(omegaa_satatm,[3,2,1]);
%     
%     dic3_temp = permute(dic_sattot,[3,2,1]);
%     ph3_temp = permute(ph_sattot,[3,2,1]);
%     r3_temp = permute(revelle_sattot,[3,2,1]);
%     oc3_temp = permute(omegac_sattot,[3,2,1]);
%     oa3_temp = permute(omegaa_sattot,[3,2,1]);
%     
%     dic4_temp = permute(dic_satalk,[3,2,1]);
%     ph4_temp = permute(ph_satalk,[3,2,1]);
%     r4_temp = permute(revelle_satalk,[3,2,1]);
%     oc4_temp = permute(omegac_satalk,[3,2,1]);
%     oa4_temp = permute(omegaa_satalk,[3,2,1]);
%     
%     dic5_temp = permute(dic_satt,[3,2,1]);
%     ph5_temp = permute(ph_satt,[3,2,1]);
%     r5_temp = permute(revelle_satt,[3,2,1]);
%     oc5_temp = permute(omegac_satt,[3,2,1]);
%     oa5_temp = permute(omegaa_satt,[3,2,1]);
%     
%     dic6_temp = permute(dic_sats,[3,2,1]);
%     ph6_temp = permute(ph_sats,[3,2,1]);
%     r6_temp = permute(revelle_sats,[3,2,1]);
%     oc6_temp = permute(omegac_sats,[3,2,1]);
%     oa6_temp = permute(omegaa_sats,[3,2,1]);
    
    % Data.....-1e34 is ferret's standard "bad" value, it doesnt seem to
    % understand NaN
    % ph_temp(isnan(ph_temp))=(-1e34);
%     ta_temp(isnan(ta_temp))=(-1e34);
%     
%     dic1_temp(isnan(dic1_temp))=(-1e34);
%     ph1_temp(isnan(ph1_temp))=(-1e34);
%     r1_temp(isnan(r1_temp))=(-1e34);
%     oc1_temp(isnan(oc1_temp))=(-1e34);
%     oa1_temp(isnan(oa1_temp))=(-1e34);
%     
%     dic2_temp(isnan(dic2_temp))=(-1e34);
%     ph2_temp(isnan(ph2_temp))=(-1e34);
%     r2_temp(isnan(r2_temp))=(-1e34);
%     oc2_temp(isnan(oc2_temp))=(-1e34);
%     oa2_temp(isnan(oa2_temp))=(-1e34);
%     
%     dic3_temp(isnan(dic3_temp))=(-1e34);
%     ph3_temp(isnan(ph3_temp))=(-1e34);
%     r3_temp(isnan(r3_temp))=(-1e34);
%     oc3_temp(isnan(oc3_temp))=(-1e34);
%     oa3_temp(isnan(oa3_temp))=(-1e34);
%     
%     dic4_temp(isnan(dic4_temp))=(-1e34);
%     ph4_temp(isnan(ph4_temp))=(-1e34);
%     r4_temp(isnan(r4_temp))=(-1e34);
%     oc4_temp(isnan(oc4_temp))=(-1e34);
%     oa4_temp(isnan(oa4_temp))=(-1e34);
%     
%     dic5_temp(isnan(dic5_temp))=(-1e34);
%     ph5_temp(isnan(ph5_temp))=(-1e34);
%     r5_temp(isnan(r5_temp))=(-1e34);
%     oc5_temp(isnan(oc5_temp))=(-1e34);
%     oa5_temp(isnan(oa5_temp))=(-1e34);
%     
%     dic6_temp(isnan(dic6_temp))=(-1e34);
%     ph6_temp(isnan(ph6_temp))=(-1e34);
%     r6_temp(isnan(r6_temp))=(-1e34);
%     oc6_temp(isnan(oc6_temp))=(-1e34);
%     oa6_temp(isnan(oa6_temp))=(-1e34);
    
    % write to file
    data_struct.T = tavesteps.tim(tstep);
    data_struct.iter = tavesteps.timesteps(tstep);
    data_struct.alk_pref = change(permute(ta,[3,2,1]),'==',NaN,-1e34);
    
    data_struct.dic_satini = change(permute(dic_satini,[3,2,1]),'==',NaN,-1e34);
    data_struct.ph_satini = change(permute(ph_satini,[3,2,1]),'==',NaN,-1e34);
    data_struct.revelle_satini = change(permute(revelle_satini,[3,2,1]),'==',NaN,-1e34);
    data_struct.omegac_satini = change(permute(omegac_satini,[3,2,1]),'==',NaN,-1e34);
    data_struct.omegaa_satini = change(permute(omegaa_satini,[3,2,1]),'==',NaN,-1e34);
    
    data_struct.dic_satatm = change(permute(dic_satatm,[3,2,1]),'==',NaN,-1e34);
    data_struct.ph_satatm = change(permute(ph_satatm,[3,2,1]),'==',NaN,-1e34);
    data_struct.revelle_satatm = change(permute(revelle_satatm,[3,2,1]),'==',NaN,-1e34);
    data_struct.omegac_satatm = change(permute(omegac_satatm,[3,2,1]),'==',NaN,-1e34);
    data_struct.omegaa_satatm = change(permute(omegaa_satatm,[3,2,1]),'==',NaN,-1e34);
    
    data_struct.dic_sattot = change(permute(dic_sattot,[3,2,1]),'==',NaN,-1e34);
    data_struct.ph_sattot = change(permute(ph_sattot,[3,2,1]),'==',NaN,-1e34);
    data_struct.revelle_sattot = change(permute(revelle_sattot,[3,2,1]),'==',NaN,-1e34);
    data_struct.omegac_sattot = change(permute(omegac_sattot,[3,2,1]),'==',NaN,-1e34);
    data_struct.omegaa_sattot = change(permute(omegaa_sattot,[3,2,1]),'==',NaN,-1e34);
    
    data_struct.dic_satalk = change(permute(dic_satalk,[3,2,1]),'==',NaN,-1e34);
    data_struct.ph_satalk = change(permute(ph_satalk,[3,2,1]),'==',NaN,-1e34);
    data_struct.revelle_satalk = change(permute(ph_satalk,[3,2,1]),'==',NaN,-1e34);
    data_struct.omegac_satalk = change(permute(omegac_satalk,[3,2,1]),'==',NaN,-1e34);
    data_struct.omegaa_satalk = change(permute(omegaa_satalk,[3,2,1]),'==',NaN,-1e34);
    
    data_struct.dic_satt = change(permute(dic_satt,[3,2,1]),'==',NaN,-1e34);
    data_struct.ph_satt = change(permute(ph_satt,[3,2,1]),'==',NaN,-1e34);
    data_struct.revelle_satt = change(permute(revelle_satt,[3,2,1]),'==',NaN,-1e34);
    data_struct.omegac_satt = change(permute(omegac_satt,[3,2,1]),'==',NaN,-1e34);
    data_struct.omegaa_satt = change(permute(omegaa_satt,[3,2,1]),'==',NaN,-1e34);
    
    data_struct.dic_sats = change(permute(dic_sats,[3,2,1]),'==',NaN,-1e34);
    data_struct.ph_sats = change(permute(ph_sats,[3,2,1]),'==',NaN,-1e34);
    data_struct.revelle_sats = change(permute(revelle_sats,[3,2,1]),'==',NaN,-1e34);
    data_struct.omegac_sats = change(permute(omegac_sats,[3,2,1]),'==',NaN,-1e34);
    data_struct.omegaa_sats = change(permute(omegaa_sats,[3,2,1]),'==',NaN,-1e34);
    
    nc_add_recs(filename,data_struct,'T')
end

% Add floating_point missing_value to all variables
eval(['!ncatted -O -a missing_value,,c,f,-1.0e34 ',filename])

if ~isempty(strfind(lower(os),'darwin'))
    % only do wait handles on OSX
    close(wait_handle)
end

% % End parallel matlabs
% if check==1 && matlabpool('size')>0
%     %checked out uni license for parallel computing toolbox
%     %existing matlab pool present
%     %therefore close matlabpool
%     matlabpool close
% end
