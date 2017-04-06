function [filename]=calc_pco2(grd,tavesteps,filename)
% DIC used in alkalinity estimation therefore creates circular argument
% about use of preformed nutrients, rendering this entirely wrong for
% use to get dic saturation. Use only for in situ values.
%   [filename]=calc_pco2(grd,tavesteps);
%   optional outputs go in this order
%   phnew,pco2,co3,co2s,ff,ffp,k0,k1,k2

%% Load variables
if nargin==0 || ~exist('grd','var')
   grd=mit_loadgrid;
end

if nargin==1 || ~exist('tavesteps','var')
   tavesteps=mit_timesteps('tave');
end

if nargin==2 || ~exist('filename','var')
    filename='ph_data.glob.nc';
end

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
sit=interp3(lat,long,depth,sitmp,yc,xc,zc);
sit=permute(sit,[2 1 3]);  zc=permute(zc,[2 1 3]);
sit=(sit./1000).*grd.hfacc; % woa01 si in umol/l convert to mol/m3
         
%% Set up Output to NetCDF
% Write out to netcdf file...
%filename='ph_data.glob.nc'; % Make an input
nc_create_empty(filename,'clobber')
%nc_padheader(filename,20000);
nc_adddim(filename,'X',grd.nx)
nc_adddim(filename,'Y',grd.ny)
nc_adddim(filename,'Z',grd.nz)
nc_adddim(filename,'T',0)

nc_addvar(filename,struct('Name','X','Datatype','double','Dimension',{{'X'}}))
nc_addvar(filename,struct('Name','Y','Datatype','double','Dimension',{{'Y'}}))
nc_addvar(filename,struct('Name','Z','Datatype','double','Dimension',{{'Z'}}))
nc_addvar(filename,struct('Name','T','Datatype','double','Dimension',{{'T'}}))

nc_attput(filename,'X','units','degrees');
nc_attput(filename,'Y','units','degrees');
nc_attput(filename,'Z','units','m');
nc_attput(filename,'T','units',[tavesteps.tuname,'s'])

% Write to NetCDF file
nc_varput(filename,'X',grd.lonc);
nc_varput(filename,'Y',grd.latc);
nc_varput(filename,'Z',-grd.zc);

nc_addvar(filename,struct('Name','iter','Datatype','int','Dimension',{{'T'}}));
nc_addvar(filename,struct('Name','new_ph_calc','Datatype','float','Dimension',{{'T','Z','Y','X'}}));
nc_addvar(filename,struct('Name','new_pco2_calc','Datatype','float','Dimension',{{'T','Z','Y','X'}}));
nc_addvar(filename,struct('Name','new_co2s_calc','Datatype','float','Dimension',{{'T','Z','Y','X'}}));
nc_addvar(filename,struct('Name','new_co3_calc','Datatype','float','Dimension',{{'T','Z','Y','X'}}));
nc_addvar(filename,struct('Name','new_hco3_calc','Datatype','float','Dimension',{{'T','Z','Y','X'}}));

% nc_attput(filename,nc_global,'title','Streamfunctions');
nc_attput(filename,'new_ph_calc','units','');
nc_attput(filename,'new_pco2_calc','units','atm');
nc_attput(filename,'new_co2s_calc','units','mol/kg');
nc_attput(filename,'new_co3_calc','units','mol/kg');
nc_attput(filename,'new_hco3_calc','units','mol/kg');

nc_attput(filename,'new_ph_calc','title','Full depth pH calculated using MITgcm routine carbon_chem.F');
nc_attput(filename,'new_pco2_calc','title','Full depth pCO2');
nc_attput(filename,'new_co2s_calc','title','Full depth concentration of CO2(aq)+H2CO3(aq)');
nc_attput(filename,'new_co3_calc','title','Full depth concentration of CO3-2 ions');
nc_attput(filename,'new_hco3_calc','title','Full depth concentration of HCO3- ions');

% % Write Global Attributes to file
% Extract attributes and axis for netcdf creation later....
if ~exist('attributes','var')
    tave = rdmnc(tavesteps.filearr(2:end-1),'Xp1');
    attributes=tave.attributes.global;
end

% Write Global Attributes to file
attnames=fieldnames(attributes);
% attributes.tile_number='global'; % DONT DO THIS!
for attnum = 1:length(attnames)
    comm=sprintf('nc_attput(filename,nc_global,''%s'',%s(attributes.(attnames{attnum})));',...
        attnames{attnum},...
        class(attributes.(attnames{attnum}))); % not sure if the class thing works or not...
    eval(comm)
    %disp(['Adding global attributes: ', comm])
end

%% Calculate some constant terms
latgrd=repmat(grd.yc,[1,1,grd.nz]).*grd.hfacc;
zgrd=nan(size(grd.hfacc));
for i=1:grd.nz
    zgrd(:,:,i)=grd.hfacc(:,:,i).*grd.zc(i);
end

pressbarc=sw_pres(zgrd,latgrd)/10;
pressbarc(isnan(pressbarc))=0;
grd.cmask(isnan(grd.cmask))=0;

[~,os]=system('uname');
if ~isempty(strfind(lower(os),'darwin'))
    % only do wait handles on OSX
    wait_handle=waitbar(0,'Calculating global pH and pCO2');
end

%% Calculate Carbon Chem Parameters
for tstep=1:tavesteps.kmax;
    if ~isempty(strfind(lower(os),'darwin'))
        % only do wait handles on OSX
        waitbar(tstep/(tavesteps.kmax),wait_handle)
    end
    
%    disp(['Processing timestep ',num2str(timesteps(l)),' (',num2str(l),')'])
     tave=rdmnc(tavesteps.filearr(2:end-1),'Ttave','Stave',tavesteps.timesteps(tstep));
     ptracer=rdmnc(strrep(tavesteps.filearr(2:end-1),'tave','ptr_tave')...
         ,'alk','dic','po4','o2',tavesteps.timesteps(tstep));
%      dictave=rdmnc(strrep(tavesteps.filearr(2:end-1),'tave','dic_tave')...
%          ,'dic_pH_ave',tavesteps.timesteps(tstep));
    
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

%     s=tave.Stave.*grd.hfacc;
%     t=tave.Ttave.*grd.hfacc;
    
%     % Use preformed phosphate??
%     [o2sat]=calc_oxygen_sat(t,s);
%     aou = o2sat-ptracer.o2 ;
%     pbio = aou/170 ;
%     pt=(ptracer.po4-pbio).*grd.hfacc;
%    pt=ptracer.po4.*grd.hfacc;
    
%     % Estimate Preformed Si ROSI=-170/15
%     sibio=(aou*15)/170;
%     sit=(si-sibio).*grd.hfacc;
    
%     % Use preformed Alkalinity
%     % Define PO=O2-R(O2:P)PO4, where R(O2:P)=-170
%     po=(ptracer.o2+(170*ptracer.po4)).*grd.hfacc;
%     ta=(0.1114+(0.062426.*tave.Stave)+(0.095896.*po)).*grd.hfacc;
%    ta=ptracer.alk.*grd.hfacc;
    
    
    % DIC used in alkalinity estimation therefore creates circular argument
    % about use of preformed nutrients, rendering this entirely wrong for
    % use to get dic saturation. Use for *in situ* only.
%    dic=ptracer.dic.*grd.hfacc;
    
%     phguess=grd.hfacc.*8.179; % pH guess as average preindustrial
%     phguess(:,:,1)=dictave.dic_pH_ave.*grd.hfacc(:,:,1); % update pH guess at the surface
%     
%     [phnew,co2s,pco2,ff,ffp,k0,k1,k2,co3]=carbon_chem(pressbarc,...
%                     tave.Ttave.*grd.hfacc,...
%                     tave.Stave.*grd.hfacc,...
%                     ptracer.dic.*grd.hfacc,...
%                     ptracer.alk.*grd.hfacc,...
%                     ptracer.po4.*grd.hfacc,...
%                     sit,phguess,grd);

      [pco2,co2s,hco3,co3,phnew,~]=mit_calc_pco2(tave.Ttave,tave.Stave,...
          pressbarc,grd.cmask,ptracer.dic,ptracer.alk,ptracer.po4,sit);

%% Output to NetCDF
    % Rearrange axes so they are in the right order...
%     ph_temp = permute(phnew,[3,2,1]);
%     ff_temp = permute(ff,[3,2,1]);
%     ffp_temp = permute(ffp,[3,2,1]);
%     k0_temp = permute(k0,[3,2,1]);
%     k1_temp = permute(k1,[3,2,1]);
%     k2_temp = permute(k2,[3,2,1]);
%     co2s_temp = permute(co2s,[3,2,1]);
%     pco2_temp = permute(pco2,[3,2,1]);
%     co3_temp = permute(co3,[3,2,1]);
% 
%     % Data.....-1e34 is ferret's standard "bad" value, it doesnt seem to
%     % understand NaN
%     ph_temp(isnan(ph_temp))=(-1e34);
%     ff_temp(isnan(ff_temp))=(-1e34);
%     ffp_temp(isnan(ffp_temp))=(-1e34);
%     k0_temp(isnan(k0_temp))=(-1e34);
%     k1_temp(isnan(k1_temp))=(-1e34);
%     k2_temp(isnan(k2_temp))=(-1e34);
%     co2s_temp(isnan(co2s_temp))=(-1e34);
%     pco2_temp(isnan(pco2_temp))=(-1e34);
%     co3_temp(isnan(co3_temp))=(-1e34);

    % write to file
    ph_data.T = tavesteps.tim(tstep);
    ph_data.iter = tavesteps.timesteps(tstep);
    ph_data.new_ph_calc = permute(change(phnew,'==',NaN,-1e34),[3,2,1]);
    ph_data.new_co2s_calc = permute(change(co2s,'==',NaN,-1e34),[3,2,1]);
    ph_data.new_pco2_calc = permute(change(pco2,'==',NaN,-1e34),[3,2,1]);
    ph_data.new_co3_calc = permute(change(co3,'==',NaN,-1e34),[3,2,1]);   
    ph_data.new_hco3_calc = permute(change(hco3,'==',NaN,-1e34),[3,2,1]);

    nc_add_recs(filename,ph_data,'T')
end

% Add floating_point missing_value to all variables
eval(['!ncatted -O -a missing_value,,c,f,-1.0e34 ',filename])

if ~isempty(strfind(lower(os),'darwin'))
    % only do wait handles on OSX
    close(wait_handle)
end

% varargout(1)={phnew};
% varargout(2)={pco2};
% varargout(3)={co3};
% varargout(4)={co2s};
% varargout(5)={ff};
% varargout(6)={ffp};
% varargout(7)={k0};
% varargout(8)={k1};
% varargout(9)={k2};


% %%%%%
% %
% % Phil Goodwin's carbonate calculations
% %
% %%%%%
%
% tas=ta(:,:,1);
% co3s=co3(:,:,1)/permil;
% dics=dic(:,:,1);
% [b,bint,r,rint,stats]=regress(co3s(:),tas(:)-dics(:));
% tas=tas(:);
% dics=dics(:);
% co3s=co3s(:);
% scatter(tas(:)-dics(:),co3s(:),50,'filled','k')
% hold on
% plot([min(tas(:)-dics(:)),max(tas(:)-dics(:))],[min(tas(:)-dics(:)).*b,max(tas(:)-dics(:)).*b],'r--')
% box on
% 
% permil=1/1024.5 ;
% cntrl_po=(cntrl_ptr.o2+(170*cntrl_ptr.po4)).*grd.hfacc;
% [cntrl_o2sat]=calc_oxygen_sat(cntrl_t,cntrl_s);
% cntrl_aou = cntrl_o2sat-cntrl_ptr.o2 ;
% cntrl_pbio = cntrl_aou/170 ;
% 
% cntrl_ta=(0.1114+(0.062426.*cntrl_tave.Stave)+(0.095896.*cntrl_po)).*grd.hfacc;
% cntrl_pt=(cntrl_ptr.po4-cntrl_pbio).*grd.hfacc;
% 
% csoft=(117.*cntrl_aou)/170;
% ccarb=0.5*(ta-cntrl_ta+(16*cntrl_aou)/170); % ta is insitu alkalinity, cntrl_ta is actually preformed alkalinity
% 
% % Have to get saturated C conc
% ncload('~/Dropbox/Applications/ferret/woa01_siz.nc')
% sitmp=permute(change(Silicate,'==',min(Silicate(:)),0),[3 2 1]);
% [lat,long,depth]=meshgrid(LATITUDE,LONGITUDE,DEPTH);
% [xc,yc,zc]=meshgrid(grd.lonc,grd.latc,(-grd.zc));
% si=interp3(lat,long,depth,sitmp,yc,xc,zc);
% si=permute(si,[2 1 3]);  zc=permute(zc,[2 1 3]);
% si=(si./1000).*grd.hfacc; % woa01 si in umol/l convert to mol/m3
% clear DEPTH LATITUDE LONGITUDE Silicate sitmp lat long depth xc yc
% cntrl_sibio=(cntrl_aou*15)/170;
% cntrl_sit=(si-cntrl_sibio).*grd.hfacc.*permumolkg;
% 
% atmos_files=dir('dic_atmos.0014400000'); i=1;
% atm_co2(i,1)=(str2double(atmos_files(i).name(11:end))*deltat)/(60*60*24*360);
% eval(['load ',atmos_files(i).name])
% atm_co2(i,2:length(dic_atmos))=dic_atmos(:,2:end);
% clear dic_atmos atmos_files i
% atm_co2(:,3)=atm_co2(:,3).*1e6; % Convert from atm to uatm
% tstep=1;
% permumolkg=permil.*1e6;
% 
% ncload('dic_sat.glob.nc')
% dic_sattot(find(dic_sattot<-1e33))=NaN;
% dic_satatm(find(dic_satatm<-1e33))=NaN;
% dic_sattot=permute(dic_sattot,[3,2,1]);
% 
% cdis=dic-csoft-ccarb-dic_sattot; % dic_sattot is due to T,S,ALKpre and pCO2
% co3=new_co3_calc/permil; 
%
% % Volume mean reference
% global_mt(k) = nansum(tave.Ttave(:).*grd.volc(:))/nansum(grd.volc(:)); % mean global temperature
% global_co3(k)= nansum(co3(:).*grd.volc(:))/(nansum(grd.volc(:)).*permil); % Global mean PO4 concentration (mol P/m3)
% global_ccarb(k)= nansum(ccarb(:).*grd.volc(:))/(nansum(grd.volc(:))); % Global mean PO4 concentration (mol P/m3)
% global_csoft(k)= nansum(csoft(:).*grd.volc(:))/(nansum(grd.volc(:))); % Global mean PO4 concentration (mol P/m3)
% global_cdis(k)= nansum(cdis(:).*grd.volc(:))/(nansum(grd.volc(:))); % Global mean PO4 concentration (mol P/m3)
% 
% % Reference to 40S/210E/1000m ish
% t0=t(78,17,8);
% csoft0=csoft(78,17,8);
% ccarb0=ccarb(78,17,8);
% cdis0=cdis(78,17,8);
% co30=co3(78,17,8);
%
% co3_calc0=b.*((ccarb-ccarb0)-(csoft-csoft0)-(cdis-cdis0)+0.0083*(t-t0));
% co3_calc_nanmean=b.*((ccarb-nanmean(ccarb(:)))-(csoft-nanmean(csoft(:)))-(cdis-nanmean(cdis(:)))+0.0083*(t-nanmean(t(:))));
% co3_calc_volmean=b.*((ccarb-global_ccarb)-(csoft-global_csoft)-(cdis-global_cdis)+0.0083*(t-global_mt));
% 
% dco30=-co30+(co3(:)/permil);
% dco3_nanmean=(co3-nanmean(co3(:)))/permil;
% dco3_volmean=(co3-global_co3)/permil;
% 
% yfit=polyval([1,0],dco30(:));
% resid_0=co3_calc0(:)-yfit;
% resid_nanmean=co3_calc_nanmean(:)-yfit;
% resid_volmean=co3_calc_volmean(:)-yfit;
% SSresid0=nansum(resid_0.^2);
% SSresid_volmean=nansum(resid_volmean.^2);
% SSresid_nanmean=nansum(resid_nanmean.^2);
% SStotal0=(length(co3_calc0(:))-1)*nanvar(co3_calc0(:));
% SStotal_nanmean=(length(co3_calc_nanmean(:))-1)*nanvar(co3_calc_nanmean(:));
% SStotal_volmean=(length(co3_calc_volmean(:))-1)*nanvar(co3_calc_volmean(:));
% rsq0=1-SSresid0/SStotal0;
% rsq_nanmean=1-SSresid_nanmean/SStotal_nanmean;
% rsq_volmean=1-SSresid_volmean/SStotal_volmean;
% 
% scatter(co3_calc0(:),dco30(:),50,'filled','k')
% hold on
% box on
% plot([-0.15,0.25],[-0.15,0.25],'r--')
