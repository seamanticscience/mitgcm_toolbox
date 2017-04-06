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

[atm_co2]=mit_getdicpco2('data');

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

latgrd=repmat(grd.yc,[1,1,grd.nz]).*grd.hfacc;
zgrd=nan(size(grd.hfacc));
for i=1:grd.nz
    zgrd(:,:,i)=grd.hfacc(:,:,i).*grd.zc(i);
end
% Pressure should be in bars by the looks of carbon_chem.F
pdbar=sw_pres(zgrd,latgrd)/10;
pdbar(isnan(pdbar))=0;
grd.cmask(isnan(grd.cmask))=0;

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

nc_attput(filename,'dic_satini','title','Saturated DIC concentration at initial/control pCO2')
nc_attput(filename,'ph_satini','title','pH at initial/control pCO2')

nc_attput(filename,'dic_satini','units','mol/kg')
nc_attput(filename,'ph_satini','units','')

% nc_attput(filename,'dic_satini','missing_value',-1.0e34)
% nc_attput(filename,'ph_satini','missing_value',-1.0e34)

% wrt prognostic pco2 only, this is dCO2(dCsat/dCO2)
nc_addvar(filename,struct('Name','dic_satatm','Datatype','float','Dimension',{{'T','Z','Y','X'}}));
nc_addvar(filename,struct('Name','ph_satatm','Datatype','float','Dimension',{{'T','Z','Y','X'}}));

nc_attput(filename,'dic_satatm','title','Saturated DIC concentration at prognostic atmospheric pCO2')
nc_attput(filename,'ph_satatm','title','pH at prognostic atmospheric pCO2')

nc_attput(filename,'dic_satatm','units','mol/kg')
nc_attput(filename,'ph_satatm','units','')

% nc_attput(filename,'dic_satatm','missing_value',-1.0e34)
% nc_attput(filename,'ph_satatm','missing_value',-1.0e34)

% wrt prognostic pco2 only, this is dCO2(dCsat/dCO2)+dT(dCsat/dT)+dALK(dCsat/dALK)+dS(dCsat/dS)
nc_addvar(filename,struct('Name','dic_sattot','Datatype','float','Dimension',{{'T','Z','Y','X'}}));
nc_addvar(filename,struct('Name','ph_sattot','Datatype','float','Dimension',{{'T','Z','Y','X'}}));

nc_attput(filename,'dic_sattot','title','Saturated DIC concentration')
nc_attput(filename,'ph_sattot','title','pH')

nc_attput(filename,'dic_sattot','units','mol/kg')
nc_attput(filename,'ph_sattot','units','')

% nc_attput(filename,'dic_sattot','missing_value',-1.0e34)
% nc_attput(filename,'ph_sattot','missing_value',-1.0e34)

% wrt alkalinity changes only
nc_addvar(filename,struct('Name','dic_satalk','Datatype','float','Dimension',{{'T','Z','Y','X'}}));
nc_addvar(filename,struct('Name','ph_satalk','Datatype','float','Dimension',{{'T','Z','Y','X'}}));

nc_attput(filename,'dic_satalk','title','Saturated DIC concentration due to ALK_PREF')
nc_attput(filename,'ph_satalk','title','pH due to ALK_PREF')

nc_attput(filename,'dic_satalk','units','mol/kg')
nc_attput(filename,'ph_satalk','units','')

% nc_attput(filename,'dic_satalk','missing_value',-1.0e34)
% nc_attput(filename,'ph_satalk','missing_value',-1.0e34)

% wrt T changes only
nc_addvar(filename,struct('Name','dic_satt','Datatype','float','Dimension',{{'T','Z','Y','X'}}));
nc_addvar(filename,struct('Name','ph_satt','Datatype','float','Dimension',{{'T','Z','Y','X'}}));

nc_attput(filename,'dic_satt','title','Saturated DIC concentration due to T')
nc_attput(filename,'ph_satt','title','pH due to T')

nc_attput(filename,'dic_satt','units','mol/kg')
nc_attput(filename,'ph_satt','units','')

% nc_attput(filename,'dic_satt','missing_value',-1.0e34)
% nc_attput(filename,'ph_satt','missing_value',-1.0e34)

% wrt S changes only
nc_addvar(filename,struct('Name','dic_sats','Datatype','float','Dimension',{{'T','Z','Y','X'}}));
nc_addvar(filename,struct('Name','ph_sats','Datatype','float','Dimension',{{'T','Z','Y','X'}}));

nc_attput(filename,'dic_sats','title','Saturated DIC concentration due to S')
nc_attput(filename,'ph_sats','title','pH due to S')

nc_attput(filename,'dic_sats','units','mol/kg')
nc_attput(filename,'ph_sats','units','')

% nc_attput(filename,'dic_sats','missing_value',-1.0e34)
% nc_attput(filename,'ph_sats','missing_value',-1.0e34)

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
    if exist('data.exf','file') && ...
            strcmpi(mit_getparm('data.pkg','useEXF'),'true') && ...
            isempty(strfind(lower(pwd),'cntrl'))
        % external forcing package runs have different starting state
        disp('Loading EXF control run')
        cntrl_tave=rdmnc([filepath,'/controlrun/cntrl_exf_global_1kyr/tave.32k.glob.nc'],'Ttave','Stave');
        cntrl_ptr=rdmnc([filepath,'/controlrun/cntrl_exf_global_1kyr/ptr_tave.32k.glob.nc'],'alk','dic','po4','o2');
    elseif exist('data.kpp','file') && ...
            strcmpi(mit_getparm('data.pkg','useKPP'),'true') && ...
            mit_getparm('data','nIter0')==18000000;
        disp('Loading KPP control run from 25k')
        cntrl_tave=rdmnc([filepath,'/controlrun/cntrl_kpp_5k/tave.25k.glob.nc'],'Ttave','Stave',18000000);
        cntrl_ptr=rdmnc([filepath,'/controlrun/cntrl_kpp_5k/ptr_tave.25k.glob.nc'],'alk','dic','po4','o2',18000000);
    elseif exist('data.kpp','file') && ...
            strcmpi(mit_getparm('data.pkg','useKPP'),'true') && ...
            mit_getparm('data','nIter0')==21600000;
        disp('Loading KPP control run from 30k')
        cntrl_tave=rdmnc([filepath,'/controlrun/cntrl_kpp_10k/tave.30k.glob.nc'],'Ttave','Stave',21600000);
        cntrl_ptr=rdmnc([filepath,'/controlrun/cntrl_kpp_10k/ptr_tave.30k.glob.nc'],'alk','dic','po4','o2',21600000);
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
    cntrl_ta=(0.1332+(0.06210.*cntrl_tave.Stave(:,:,:,end))+(0.06737.*cntrl_po(:,:,:,end))).*grd.hfacc;
else
    disp('use standard (Relaxing) values for control preformed alkalinity...')
    cntrl_ta=(0.1114+(0.062426.*cntrl_tave.Stave(:,:,:,end))+(0.095896.*cntrl_po(:,:,:,end))).*grd.hfacc;
end

% Use preformed phosphate (in umol/kg)
[cntrl_o2sat]=calc_oxygen_sat(cntrl_t,cntrl_s);
cntrl_aou = cntrl_o2sat-cntrl_ptr.o2(:,:,:,end) ;
cntrl_pbio = cntrl_aou/170 ;
cntrl_pt=(cntrl_ptr.po4(:,:,:,end)-cntrl_pbio).*grd.hfacc;

% Estimate Preformed Si ROSI=-170/15 (in umol/kg)
% Is constant as not carried by the model - from WOA actually
cntrl_sibio=(cntrl_aou*15)/170;
cntrl_sit=(si-cntrl_sibio).*grd.hfacc;

latgrd=repmat(grd.yc,[1,1,grd.nz]).*grd.hfacc;
zgrd=nan(size(grd.hfacc));
for i=1:grd.nz
    zgrd(:,:,i)=grd.hfacc(:,:,i).*grd.zc(i);
end
    
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
        ta=(0.1332+(0.06210.*tave.Stave)+(0.06737.*po)).*grd.hfacc;
    else
        % use standard values
        ta=(0.1114+(0.062426.*tave.Stave)+(0.095896.*po)).*grd.hfacc;
    end
    
    % Following Nutrients are used to calculate Carbonate Alkalinity, which
    % is then used to get DIC conc (or Saturated DIC in this
    % case)..probably small effect
    % Use preformed phosphate (in umol/kg)
    [o2sat]=calc_oxygen_sat(t,s);
    aou = o2sat-ptracer.o2 ;
    pbio = aou/170 ;
    pt=(ptracer.po4-pbio).*grd.hfacc;
    
    % Estimate Preformed Si ROSI=-170/15 (in umol/kg)
    sibio=(aou*15)/170;
    sit=(si-sibio).*grd.hfacc;
    
%% Use MEX version of Mick's efficient fortran CSAT solver

            % Calculate Saturated DIC wrt prognostic atm_pco2 at each timestep
            % This is dCO2(dCsat/dCO2)+dT(dCsat/dT)+dALK(dCsat/dALK)+dS(dCsat/dS)
  			[dic_sattot,~,ph_sattot]=mit_calc_csat(t,s,pdbar,grd.cmask,...
  			            atm_co2(atm_co2(:,1)==tavesteps.tim(tstep),3),ta,pt,sit);
            
            % Calculate Saturated DIC wrt atm_pco2 = 278.05 (control)
            % This is dT(dCsat/dT)+dALK(dCsat/dALK)+dS(dCsat/dS)+dNUTS(dCsat/dNUTS)
  			[dic_satini,~,ph_satini]=mit_calc_csat(t,s,pdbar,grd.cmask,atm_co2(1,3),ta,pt,sit);

            % Calculate Saturated DIC with only pco2 changes.
            % This is dCO2(dCsat/dCO2)
  			[dic_satatm,~,ph_satatm]=mit_calc_csat(cntrl_t,cntrl_s,pdbar,grd.cmask,...
  			            atm_co2(atm_co2(:,1)==tavesteps.tim(tstep),3),cntrl_ta,cntrl_pt,cntrl_sit);
            
            % Calculate Saturated DIC with only ALK.
            % This is dALK(dCsat/dALK)
  			[dic_satalk,~,ph_satalk]=mit_calc_csat(cntrl_t,cntrl_s,pdbar,grd.cmask,...
  			            atm_co2(1,3),ta,cntrl_pt,cntrl_sit);
            
            % Calculate Saturated DIC with only T changes.
            % This is dT(dCsat/dT)
  			[dic_satt,~,ph_satt]=mit_calc_csat(t,cntrl_s,pdbar,grd.cmask,...
  			            atm_co2(1,3),cntrl_ta,cntrl_pt,cntrl_sit);
  			                        
            % Calculate Saturated DIC with only S changes.
            % This is dS(dCsat/dS)
  			[dic_sats,~,ph_sats]=mit_calc_csat(cntrl_t,s,pdbar,grd.cmask,...
  			            atm_co2(1,3),cntrl_ta,cntrl_pt,cntrl_sit);
  			            
%% Output to NetCDF File    
    % write to file
    data_struct.T = tavesteps.tim(tstep);
    data_struct.iter = tavesteps.timesteps(tstep);
    data_struct.alk_pref = change(permute(ta,[3,2,1]),'==',NaN,-1e34);
    
    data_struct.dic_satini = change(permute(dic_satini,[3,2,1]),'==',NaN,-1e34);
    data_struct.ph_satini = change(permute(ph_satini,[3,2,1]),'==',NaN,-1e34);
    
    data_struct.dic_satatm = change(permute(dic_satatm,[3,2,1]),'==',NaN,-1e34);
    data_struct.ph_satatm = change(permute(ph_satatm,[3,2,1]),'==',NaN,-1e34);
    
    data_struct.dic_sattot = change(permute(dic_sattot,[3,2,1]),'==',NaN,-1e34);
    data_struct.ph_sattot = change(permute(ph_sattot,[3,2,1]),'==',NaN,-1e34);
    
    data_struct.dic_satalk = change(permute(dic_satalk,[3,2,1]),'==',NaN,-1e34);
    data_struct.ph_satalk = change(permute(ph_satalk,[3,2,1]),'==',NaN,-1e34);
    
    data_struct.dic_satt = change(permute(dic_satt,[3,2,1]),'==',NaN,-1e34);
    data_struct.ph_satt = change(permute(ph_satt,[3,2,1]),'==',NaN,-1e34);
    
    data_struct.dic_sats = change(permute(dic_sats,[3,2,1]),'==',NaN,-1e34);
    data_struct.ph_sats = change(permute(ph_sats,[3,2,1]),'==',NaN,-1e34);
    
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
