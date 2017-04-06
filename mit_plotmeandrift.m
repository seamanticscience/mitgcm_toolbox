function mit_plotmeandrift(grd,tavesteps,diagsteps)
% mean temperature and salinity drift

% $Header: /u/gcmpack/MITgcm/verification/tutorial_global_oce_latlon/diags_matlab/mit_plotmeandrift.m,v 1.3 2006/08/12 20:25:13 jmc Exp $
% $Name:  $
if nargin==0 || ~exist('grd','var')
    grd=mit_loadgrid;
end

if nargin==1 || ~exist('tavesteps','var')
    tavesteps=mit_timesteps('tave');
end

if nargin==2 || ~exist('diagsteps','var')
    diagsteps=mit_timesteps(mit_getparm('data.diagnostics','filename'));
end

area = grd.rac.*grd.hfacc(:,:,1) ;
r_cp = 117; % Redfield Ratio of Carbon to Phosphate in phytoplankton
% r_fep = ???; % Redfield Ratio of Iron to ??? in Phytoplankton 

% From Tave File
global_mt = zeros(size(1:length(tavesteps.timesteps)));
global_ms = zeros(size(1:length(tavesteps.timesteps)));

for k=1:length(tavesteps.timesteps);
    tave=rdmnc(tavesteps.filearr(2:end-1),'Ttave','Stave',tavesteps.timesteps(k));

    global_mt(k) = nansum(tave.Ttave(:).*grd.volc(:))/nansum(grd.volc(:)); % mean global temperature
    global_ms(k) = nansum(tave.Stave(:).*grd.volc(:))/nansum(grd.volc(:)); % mean global salinity
end

if ~exist('attributes','var')
    attributes=tave.attributes.global;
end

% transport through Drake Passage
tdp = NaN*ones(length(tavesteps.timesteps),1);
tdpz=NaN*ones(grd.nz,length(tavesteps.timesteps));
if grd.nx == 128 && grd.ny == 64
    kx = 105;
    kyg = 5:20;
elseif grd.nx == 360 && grd.ny >= 160 
    kx=294;
    kyg = 5:27;
else
    error('Drake Passage coordinates need to be specified for this configuration')
end
da = grd.dz*grd.dyg(kx,kyg);

for k=1:length(tavesteps.timesteps);
  tave = rdmnc(tavesteps.filearr(2:end-1),'uVeltave','vVeltave',tavesteps.timesteps(k));
  if ~strcmp(grd.buoyancy,'OCEANIC');
    u = tave.uVeltave(:,:,end:-1:1);
    v = tave.vVeltave(:,:,end:-1:1);
  else
    u = tave.uVeltave;
    v = tave.vVeltave;
  end
  
  tdp(k) = sum(nansum(squeeze(u(kx,kyg,:))'.*da))*1e-6; % in Sv
  tdpz(:,k) = nansum(squeeze(u(kx,kyg,:)).*da')'*1e-6; % in Sv
end

clear tave u v

global_alk = zeros(size(1:length(tavesteps.timesteps)));
global_dic = zeros(size(1:length(tavesteps.timesteps)));
global_dop = zeros(size(1:length(tavesteps.timesteps)));
global_po4 = zeros(size(1:length(tavesteps.timesteps)));
total_dic = zeros(size(1:length(tavesteps.timesteps)));
total_dop = zeros(size(1:length(tavesteps.timesteps)));
total_po4 = zeros(size(1:length(tavesteps.timesteps)));
total_flxco2 = zeros(size(1:length(tavesteps.timesteps)));
% on the off chance that Fe isnt included, wont take up too much memory
global_fe=zeros(size(1:length(tavesteps.timesteps)));
total_fe=zeros(size(1:length(tavesteps.timesteps))); 
global_ligand=zeros(size(1:length(tavesteps.timesteps)));

for k=1:length(tavesteps.timesteps);
    if nc_isvar(strrep(tavesteps.filearr(2:end-1),'tave','ptr_tave'),'ligand')
        ptracer=rdmnc(strrep(tavesteps.filearr(2:end-1),'tave','ptr_tave')...
            ,'alk','dic','dop','po4','fe','ligand',tavesteps.timesteps(k));
    
        global_ligand(k)= nansum(ptracer.ligand(:).*grd.volc(:))/nansum(grd.volc(:)); % mean global alkalinity conc (mol eq/m3)
    else
        ptracer=rdmnc(strrep(tavesteps.filearr(2:end-1),'tave','ptr_tave')...
            ,'alk','dic','dop','po4','fe',tavesteps.timesteps(k));
    end
    
    global_alk(k)= nansum(ptracer.alk(:).*grd.volc(:))/nansum(grd.volc(:)); % mean global alkalinity conc (mol eq/m3)
    global_dic(k)= nansum(ptracer.dic(:).*grd.volc(:))/nansum(grd.volc(:)); % mean DIC conc (mol C/m3)
    total_dic(k) = nansum(ptracer.dic(:).*grd.volc(:)); % Total moles of inorganic C in ocean (mol C)
    global_dop(k)= nansum(ptracer.dop(:).*grd.volc(:))/nansum(grd.volc(:)); % mean DOP conc (mol P/m3)
    total_dop(k) = nansum(ptracer.dop(:).*grd.volc(:)); % Total moles of organic P in ocean (mol P)
    global_po4(k)= nansum(ptracer.po4(:).*grd.volc(:))/nansum(grd.volc(:)); % Global mean PO4 concentration (mol P/m3)
    total_po4(k)= nansum(ptracer.po4(:).*grd.volc(:)); % Total moles of inorganic P (mol P)
    if isvar('ptracer.fe');
        global_fe(k)=nansum(ptracer.fe(:).*grd.volc(:))/nansum(grd.volc(:)); % Global Mean Conc of Iron (mol Fe/m3)
        total_fe(k)=nansum(ptracer.fe(:).*grd.volc(:)); % Total Fe in ocean (mol Fe)
    end
end
clear ptracer

global_pco2 = zeros(size(1:length(tavesteps.timesteps)));
global_flxco2 = zeros(size(1:length(tavesteps.timesteps)));
total_bio = zeros(size(1:length(tavesteps.timesteps)));

for k=1:length(tavesteps.timesteps);
        dic=rdmnc(strrep(tavesteps.filearr(2:end-1),'tave','dic_tave')...
            ,'dic_pCO2_ave','dic_BIO_ave','dic_fluxCO2_ave',tavesteps.timesteps(k));

    % sea surface only therefore weighted by area! 
    global_pco2(k)= nansum(dic.dic_pCO2_ave(:).*1e6.*area(:))/nansum(area(:)); % Global mean surface pCO2 (uatm)
    %tbk=(:,:,:,k).*31104000.*r_cp.*12./1e15; % [GtC m-3 yr-1], not [molP m-3 s-1]
    total_bio(k) = nansum(dic.dic_BIO_ave(:).*grd.volc(:).*31104000.*r_cp.*12./1e15); % Total bio production (GtC/yr)
    total_flxco2(k)= nansum(dic.dic_fluxCO2_ave(:).*area(:).*31104000); % Total air-sea flux of CO2 (mol C/yr)
    global_flxco2(k)= nansum(dic.dic_fluxCO2_ave(:).*area(:).*31104000)/nansum(area(:)); % Global mean air-sea CO2 flux
end
clear dic

% From the Surface Diagnostics
% Preallocate:
global_qnet = zeros(length(diagsteps.tim),1);
global_fwflx = zeros(length(diagsteps.tim),1);
global_tflux = zeros(length(diagsteps.tim),1);
global_sflux = zeros(length(diagsteps.tim),1);
% on the off chance that these arent included, wont take up too much memory
global_trelax = zeros(length(diagsteps.tim),1);
global_srelax = zeros(length(diagsteps.tim),1);
global_ocesflux = zeros(length(diagsteps.tim),1);
global_ocefreez = zeros(length(diagsteps.tim),1);
global_surforcs = zeros(length(diagsteps.tim),1);
global_surforct = zeros(length(diagsteps.tim),1);
%global_oceqsw = zeros(length(diagsteps.tim),1);

for k=1:length(diagsteps.tim);
    surfDiag=rdmnc(diagsteps.filearr(2:end-1)...
            ,'oceQnet','oceFWflx','TFLUX','SFLUX','TRELAX','SRELAX'...
            ,'oceSflux','oceFreez','surForcS','surForcT',diagsteps.timesteps(k));
   
   global_qnet(k) = nansum(surfDiag.oceQnet(:).*area(:))/nansum(area(:));
   global_fwflx(k) = nansum(surfDiag.oceFWflx(:).*area(:).*31104000)/nansum(area(:)); % per year, not per s
   global_tflux(k) = nansum(surfDiag.TFLUX(:).*area(:))/nansum(area(:));
   global_sflux(k) = nansum(surfDiag.SFLUX(:).*area(:).*31104000)/nansum(area(:));

   if isvar('surfDiag.TRELAX'); 
        global_trelax(k) = nansum(surfDiag.TRELAX(:).*area(:))/nansum(area(:));
   end
   
   if isvar('surfDiag.oceQsw'); 
        global_oceqsw(k) = nansum(surfDiag.oceQsw(:).*area(:))/nansum(area(:));
   end
   
   if isvar('surfDiag.SRELAX');
       global_srelax(k) = nansum(surfDiag.SRELAX(:).*area(:).*31104000)/nansum(area(:));
   end
   
   if isvar('surfDiag.oceSflux');
       global_ocesflux(k) = nansum(surfDiag.oceSflux(:).*area(:).*31104000)/nansum(area(:));
   end

   if isvar('surfDiag.oceFreez');
       global_ocefreez(k) = nansum(surfDiag.oceFreez(:).*area(:))/nansum(area(:));
   end

   if isvar('surfDiag.surForcS');
       global_surforcs(k) = nansum(surfDiag.surForcS(:).*area(:).*31104000)/nansum(area(:));
   end

   if isvar('surfDiag.surForcT');
       global_surforct(k) = nansum(surfDiag.surForcT(:).*area(:))/nansum(area(:));
   end
end

clear surfDiag atmos_files

%%% IF we use the exf package, load this up
if exist(strrep(diagsteps.filearr(2:end-1),'surf','exf'),'file')
    rhonil=1035;
    global_qnet = zeros(length(diagsteps.tim),1);   % Qnet
    global_oceqsw = zeros(length(diagsteps.tim),1); % Qsw
    global_lwflux = zeros(length(diagsteps.tim),1); % Qlwnet
    global_latent = zeros(length(diagsteps.tim),1); % latent 
    global_sensib = zeros(length(diagsteps.tim),1); % sensible
    
    global_fwflx = zeros(length(diagsteps.tim),1);  % Net fw flux
    global_exfevap = zeros(length(diagsteps.tim),1);% evap
    global_exfprec = zeros(length(diagsteps.tim),1);% precip
    global_exfroff = zeros(length(diagsteps.tim),1);% runoff
    
    for k=1:length(diagsteps.tim);
        exfDiag=rdmnc(strrep(diagsteps.filearr(2:end-1),'surf','exf')...
            ,'EXFempmr','EXFevap','EXFroff','EXFpreci','EXFhl','EXFhs'...
            ,'EXFlwnet','EXFswnet','EXFqnet',diagsteps.timesteps(k));
        
        global_qnet(k) = nansum(-1.*exfDiag.EXFqnet(:).*area(:))/nansum(area(:));
        global_oceqsw(k) = nansum(-1.*exfDiag.EXFswnet(:).*area(:))/nansum(area(:));
        global_lwflux(k) = nansum(-1.*exfDiag.EXFlwnet(:).*area(:))/nansum(area(:));
        global_latent(k) = nansum(exfDiag.EXFhl(:).*area(:))/nansum(area(:));
        global_sensib(k) = nansum(exfDiag.EXFhs(:).*area(:))/nansum(area(:));
        
        global_fwflx(k) = nansum(exfDiag.EXFempmr(:).*area(:).*31104000.*rhonil)/nansum(area(:)); % per year, not per s
        global_exfevap(k) = nansum(exfDiag.EXFevap(:).*area(:).*31104000.*rhonil)/nansum(area(:));
        global_exfprec(k) = nansum(-1.*exfDiag.EXFpreci(:).*area(:).*31104000.*rhonil)/nansum(area(:));
        global_exfroff(k) = nansum(-1.*exfDiag.EXFroff(:).*area(:).*31104000.*rhonil)/nansum(area(:));
    end
end


tpo4=(total_dop+total_po4); % Phosphate conservation (mol P)
tco2=total_dic+(total_dop.*r_cp); % Oceanic Carbon conservation (mol C)

% dicDiags to look for virtual fluxes
global_vflux = zeros(length(diagsteps.tim),1);

for k=1:length(diagsteps.tim);
    if nc_isvar(strrep(diagsteps.filearr(2:end-1),'surf','dic'),'DICVCFLX')
        dicDiag=rdmnc(strrep(diagsteps.filearr(2:end-1),'surf','dic'),'DICVCFLX',diagsteps.timesteps(k));
        global_vflux(k) = nansum((dicDiag.DICVCFLX(:).*area(:).*31104000)/grd.dz(1)); % diagnostic is im mol C/m3/s in surface layer
    elseif nc_isvar(strrep(diagsteps.filearr(2:end-1),'surf','dic_surf'),'DICVCFLX')
        dicDiag=rdmnc(strrep(diagsteps.filearr(2:end-1),'surf','dic_surf'),'DICVCFLX',diagsteps.timesteps(k));
        global_vflux(k) = nansum((dicDiag.DICVCFLX(:).*area(:).*31104000)/grd.dz(1)); % diagnostic is im mol C/m3/s in surface layer
    else
        global_vflux(k) = 0;
    end
end

clear dicDiag

atm_co2=mit_getdicpco2('data');
dicint = mit_getparm('data.dic','DIC_int1');
% % Use DIC_int1 from data.dic to decide which way to handle the atmos_co2
% % files (if any):
% % dic_int1:
% %  0=use default 278.d-6
% %  1=use constant value - dic_pCO2, read in from data.dic
% %  2=read in from file
% %  3=interact with atmospheric box (use dic_pCO2 as initial atmos. value)
% %  4=atmospheric box with specified emissions (bonus code by jml)
% if exist('data.dic','file')
[~,os]=system('uname');
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


if exist([filepath,'/controlrun/control.mat'],'file') ...
        && ( ~strcmp(attributes.the_run_name,'Control Run') ...
        || ~strcmp(attributes.the_run_name,'cntrl') )...
        && min([mit_getparm('data','startTime'),mit_getparm('data','nIter0')])>0;
    if strcmp(tavesteps.timeunit,'yrs')
        timeunit='kyrs';
    end
        
%    if strfind(lower(attributes.the_run_name),'visbeck')
        % Visbeck runs start at 25k and have different starting state
%        disp('Loading Visbeck isopycnal mixing control run')
%        con_file='/Volumes/PhD_Data/controlrun/con_visbeck_variable_k/controlrun_global/control_visk_last1kyr.mat';
%     elseif strfind(attributes.the_run_name,'elax') && strfind(attributes.the_run_name,'ixed')
%         % Fixed Relaxation runs start at 22k and have a slightly different
%         % starting state
%         disp('Loading Fixed Relaxation control run')
%         load /Volumes/PhD_Data/controlrun/altix_fixedrelax/fixedrelaxrun_global/control_fbc_last_1kyr.mat
%else
    if exist('data.exf','file') && ...
            strcmpi(mit_getparm('data.pkg','useEXF'),'true') && ...
            isempty(strfind(lower(pwd),'cntrl'))
        % external forcing package runs have different starting state
        disp('Loading EXF control run')
        con_file=[filepath,'/controlrun/control_exf'];
    elseif exist('data.kpp','file') && ...
            strcmpi(mit_getparm('data.pkg','useKPP'),'true') && ...
            mit_getparm('data','nIter0')==18000000;
        disp('Loading KPP control run from 25k')
        con_file=[filepath,'/controlrun/control_kpp_5k'];
    elseif exist('data.kpp','file') && ...
            strcmpi(mit_getparm('data.pkg','useKPP'),'true') && ...
            mit_getparm('data','nIter0')==21600000;
        disp('Loading KPP control run from 30k')
        con_file=[filepath,'/controlrun/control_kpp_10k'];
    elseif ~isempty(strfind(lower(pwd),'vmhs97')) && ...
        isempty(strfind(lower(pwd),'cntrl'))
        disp('Loading vmhs97 model control run')
        con_file=[filepath,'/perturb_eddies/VMHS97/cntrl_vmhs97/cntrl_vmhs97_global/control_vmhs97'];
    elseif ~isempty(strfind(lower(pwd),'hmm2011')) && ...
        isempty(strfind(lower(pwd),'cntrl'))
        disp('Loading hmm2011 model control run')
        con_file=[filepath,'/perturb_eddies/HMM2011/cntrl_hmm2011/cntrl_hmm2011_global/control_hmm2011'];
    elseif ~isempty(strfind(lower(pwd),'fmh2005')) && ...
        isempty(strfind(lower(pwd),'cntrl'))
        disp('Loading fmh2005 model control run')
        con_file=[filepath,'/perturb_eddies/FMH2005/cntrl_fmh2005/cntrl_fmh2005_global/control_fmh2005'];
    else
        disp('Loading standard model control run')
        con_file=[filepath,'/controlrun/control'];
    end
    
    load(con_file)
    %Back up old variables
    save tmpbackup.mat tavesteps diagsteps tdp global* total*
    if isvar('atm_co2'); save tmpbackup.mat atm_co2 -append; end
    if isvar('tallco2');  save tmpbackup.mat tallco2 -append; end

    if strfind(attributes.the_run_name,'Geotraces')
        % Is this a short geotraces run?
        tave_yrs=[0,tavesteps.tim];
        diag_yrs=[0,diagsteps.tim];
        idx=length(ctim);
        
        if isvar('atm_co2');
            if dicint==4;
                atm_co2=[[catm_co2(idx,:),zeros(length(atm_co2(idx,:)),1)];atm_co2];
            else
                atm_co2=[catm_co2(idx,:);atm_co2];
            end
        end
        atm_co2(1,1)=0; % reset co2 clock
    else
%        size(ctim)
%        size(tavesteps.tim)
        tave_yrs=[ctim(:)'./1000,tavesteps.tim(:)'./1000];
        diag_yrs=[ctim_diag(:)'./1000,diagsteps.tim(:)'./1000];
        idx=':';
        if isvar('atm_co2');
            if dicint==4;
                atm_co2=[[catm_co2(idx,:),zeros(length(atm_co2(idx,:)),1)];atm_co2];
            else
                atm_co2=[catm_co2(idx,:);atm_co2];
            end
        end
    end
    
    global_vflux_cntrl=cglobal_vflux(length(cglobal_vflux));
    
    clear cTDP catm_co2 ctpo4 cdiagsteps ctavesteps ctallco2 cglobal_vflux
    tim=tavesteps.tim;
    tim_diag=diagsteps.tim;
    
    vars=who('c*');
    for ii=1:length(vars)
        if ~strcmp(vars{ii},'con_file')
            eval(['tmp=',vars{ii},'(idx);']);
            eval([vars{ii}(2:end),'=[tmp(:);',vars{ii}(2:end),'(:)];'])
            clear tmp
        end
    end
    
    if isvar('tpo4')
        tmp=(ctotal_dop(idx)+ctotal_po4(idx)); % Phosphate conservation (mol P)
        tpo4=[tmp(:);tpo4(:)];
    end
    
    if isvar('tco2')
        tmp=ctotal_dic+(ctotal_dop.*r_cp);
        tco2=[tmp(:);tco2(:)];
        clear tmp
    end
    
    if isvar('atm_co2')
        % get rid of any repetitive times
        [c,ia]=unique(atm_co2(:,1));
        atm_co2=atm_co2(ia,:);
        clear c ia
    end
    
    if isvar('global_vflux_cntrl');
        diagfreq=mit_getparm('data.diagnostics','frequency(2)');
        % Since this will be compared with DIC+RCP*DOP, which integrates the
        % loss/gain of carbon by the ocean, should integrate the virtual flux.
        global_vflux_anom=cumsum(global_vflux-global_vflux_cntrl).*(diagfreq./31104000); % because of 10 year averages. % control value
        dvflux=global_vflux_anom(end);
    end
    
%     global_ms=[cglobal_ms(idx)',global_ms];
%     global_mt=[cglobal_mt(idx)',global_mt];
%     global_dic=[cglobal_dic(idx)',global_dic];
%     total_dic=[ctotal_dic(idx)',total_dic];
%     global_alk=[cglobal_alk(idx)',global_alk];
%     global_dop=[cglobal_dop(idx)',global_dop];
%     total_dop=[ctotal_dop(idx)',total_dop];
%     global_po4=[cglobal_po4(idx)',global_po4];
%     total_po4=[ctotal_po4(idx)',total_po4];
%     if isvar('global_fe');
%         global_fe=[cglobal_fe(idx)',global_fe];
%         total_fe=[ctotal_fe(idx)',total_fe];
%     end
%     global_pco2=[cglobal_pco2(idx)',global_pco2];
%     total_bio=[ctotal_bio(idx)',total_bio];
%     global_flxco2=[cglobal_flxco2(idx)',global_flxco2];
%     total_flxco2=[ctotal_flxco2(idx)',total_flxco2];
%     global_qnet=[cglobal_qnet(idx);global_qnet];
%     global_fwflx=[cglobal_fwflx(idx);global_fwflx];
%     global_tflux=[cglobal_tflux(idx);global_tflux];
%     global_sflux=[cglobal_sflux(idx);global_sflux];
%     if isvar('global_trelax'); global_trelax=[cglobal_trelax(idx);global_trelax]; end
%     if isvar('global_srelax'); global_srelax=[cglobal_srelax(idx);global_srelax]; end
%     if isvar('global_ocesflux'); global_ocesflux=[cglobal_ocesflux(idx);global_ocesflux]; end
%     if isvar('global_ocefreez'); global_ocefreez=[cglobal_ocefreez(idx);global_ocefreez]; end
%     if isvar('global_surforcs'); global_surforcs=[cglobal_surforcs(idx);global_surforcs]; end
%     if isvar('global_surforct'); global_surforct=[cglobal_surforct(idx);global_surforct]; end
%     if isvar('global_oceqsw');   global_oceqsw=[zeros(size(cglobal_surforct(idx)));global_oceqsw'];end
%     if isvar('tallco2'); tallco2=[ctallco2(idx);tallco2]; end
    clearvars c* vars ii -except con_file

    load(con_file,'ctim','ctim_diag');
    ctim=ctim./1000;
    ctim_diag=ctim_diag./1000;
else
    if strcmp(tavesteps.timeunit,'yrs')
        timeunit='kyrs';
        tave_yrs=tavesteps.tim./1000;
        diag_yrs=diagsteps.tim./1000;
    else
        tave_yrs=tavesteps.tim;
        diag_yrs=diagsteps.tim;
    end
end

if isvar('atm_co2')
    chkptfreq = mit_getparm('data','pChkptFreq');
    tavefreq = mit_getparm('data','taveFreq');
    if (chkptfreq<tavefreq); %more pickups that tave-output
        freqrat=tavefreq/chkptfreq;
        % subsample atm_co2
        temptco2=tco2';
        tempatm_co2=atm_co2(1:freqrat:end,:);
    elseif (chkptfreq>tavefreq); % more tave-output than pickups
        freqrat=chkptfreq/tavefreq;
        % subsample tco2
        temptco2=tco2(:,freqrat:freqrat:end)';
        tempatm_co2=atm_co2;
    else
        % Do nothing, but rename variables
        temptco2=tco2';
        tempatm_co2=atm_co2;
    end
    clear chkptfreq tavefreq
    % Sorted out sampling freq, now sort out length of run. Only provided
    % atm_co2 is equal or less than because it is unlikely that atm_co2
    % will be longer. Similar assumption for any extension to go at the
    % begining of a run. tallco2 is only valid where atmospheric box is ON,
    % therefore extension=NaN and tco2+NaN=NaN.
    if (length(tempatm_co2(:,1)))==length(temptco2) % was -1 in for length(tempatm_co2)
        tallco2=temptco2(:)+tempatm_co2(:,2); % was 2:end for tempatm_co2
    elseif (length(tempatm_co2(:,1)))<length(temptco2)
        ldiff=length(temptco2)-(length(tempatm_co2(:,1)));
        extension=NaN(ldiff,1);
        tallco2=temptco2(:)+[extension;tempatm_co2(:,2)]; % "tco2" is ocean carbon while "extension" and "atm_co2" are the atmospheric C contribution to total system carbon
    elseif (length(tempatm_co2(:,1)))==length(temptco2)+1
        % Here, atm_co2 will sometimes contain two initial atmco2 values,
        % e.g. at 20k and 25k whereas "tim" with start at 20010...hence
        % drop first value
        tallco2=temptco2(:)+tempatm_co2(2:end,2);
    end
end

dTemp = global_mt(end) - global_mt(max(1,end-tavesteps.kmax-1));
dSal = global_ms(end) - global_ms(max(1,end-tavesteps.kmax-1));
dpsi = tdp(end) - tdp(max(1,end-tavesteps.kmax-1));

if exist(strrep(diagsteps.filearr(2:end-1),'surf','exf'),'file') % Use EXF
    if (length(global_exfevap)-diagsteps.kmax-1)>0; then 
        idx=length(global_exfevap)-diagsteps.kmax-1;
    else
        idx=1;
    end
    
    dqnet = global_qnet(end) - global_qnet(idx);
    dqsw = global_oceqsw(end) - global_oceqsw(idx);
    dqlw = global_lwflux(end) - global_lwflux(idx);    
    dlat = global_latent(end) - global_latent(idx);
    dsen = global_sensib(end) - global_sensib(idx);    
    
    dfwf = global_fwflx(end) - global_fwflx(idx);
    devap = global_exfevap(end) - global_exfevap(idx);
    dprec = global_exfprec(end) - global_exfprec(idx);
    droff = global_exfroff(end) - global_exfroff(idx);
else
    dqnet = global_qnet(end) - global_qnet(max(1,end-diagsteps.kmax-1));
    dfwf = global_fwflx(end) - global_fwflx(max(1,end-diagsteps.kmax-1));
    dtf = global_tflux(end) - global_tflux(max(1,end-diagsteps.kmax-1));
    dsf = global_sflux(end) - global_sflux(max(1,end-diagsteps.kmax-1));
    if isvar('global_trelax'); dtr = global_trelax(end) - global_trelax(max(1,end-diagsteps.kmax-1));end
    if isvar('global_srelax'); dsr = global_srelax(end) - global_srelax(max(1,end-diagsteps.kmax-1));end
    if isvar('global_ocesflux'); dosf = global_ocesflux(end) - global_ocesflux(max(1,end-diagsteps.kmax-1));end
    if isvar('global_ocefreez'); dof = global_ocefreez(end) - global_ocefreez(max(1,end-diagsteps.kmax-1));end
    if isvar('global_surforct'); dsft = global_surforct(end) - global_surforct(max(1,end-diagsteps.kmax-1));end
    if isvar('global_surforcs'); dsfs = global_surforcs(end) - global_surforcs(max(1,end-diagsteps.kmax-1));end
    if isvar('global_oceqsw');   dqsw = global_oceqsw(end) - global_oceqsw(max(1,end-diagsteps.kmax-1));end
end
dalk = global_alk(end) - global_alk(max(1,end-tavesteps.kmax-1));
ddic = global_dic(end) - global_dic(max(1,end-tavesteps.kmax-1));

dpo4 = global_po4(end) - global_po4(max(1,end-tavesteps.kmax-1));
dbio = total_bio(end) - total_bio(max(1,end-tavesteps.kmax-1));
if isvar('global_fe');dfe = global_fe(end) - global_fe(max(1,end-tavesteps.kmax-1));end
dpco2 = global_pco2(end) - global_pco2(max(1,end-tavesteps.kmax-1));
dflxco2 = global_flxco2(end) - global_flxco2(max(1,end-tavesteps.kmax-1));
dtflxco2 = total_flxco2(end) - total_flxco2(max(1,end-tavesteps.kmax-1));
if isvar('tallco2'); dtallco2 = tallco2(end) - tallco2(max(1,end-tavesteps.kmax-1)); end
dtpo4 = tpo4(end) - tpo4(max(1,end-tavesteps.kmax-1));
dtco2 = tco2(end) - tco2(max(1,end-tavesteps.kmax-1));

if isvar('atm_co2')
    datmc=atm_co2(end,2)-atm_co2(max(1,end-tavesteps.kmax),2);
    datmp=(atm_co2(end,3)-atm_co2(max(1,end-tavesteps.kmax),3))*1e6;
end

xlim = [(tavesteps.tim(1)/1000)-1,tave_yrs(end)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURE 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
subplot(3,2,1);
[AX,H1,H2] = plotyy(tave_yrs,global_ms,tave_yrs,global_mt);
set(get(AX(1),'Ylabel'),'String','S [psu]')
set(get(AX(2),'Ylabel'),'String','T [degC]')
set(AX,'XLim',xlim,'YTickMode','auto');
set(AX(2),'YColor','r');set(H2,'Color','r');
if isvar('ctim')
   hold on
   ylim=get(gca,'YLim');
   H3=plot([ctim(end),ctim(end)],ylim,'k--');
end
titlevec={['Drift of Mean T (\Delta\theta = ',num2str(dTemp,3),') and Mean S (\DeltaS = ',num2str(dSal,3),')']};  
title(titlevec); % xlabel(['Time [' timeunit ']']);


BX=subplot(3,2,2);
plot(tave_yrs,tdp);set(gca,'XLim',xlim);
titlevec={['Drift of Transport through Drake Passage,(',num2str(tdp(end)),'Sv, \Delta\psi = ',num2str(dpsi,3),')']};
if isvar('ctim')
   hold on
   ylim=get(gca,'YLim');
   plot([ctim(end),ctim(end)],ylim,'k--');
end
title(titlevec); ylabel('Transport [Sv]'); % xlabel(['Time [' timeunit ']']);

suptitle(attributes.the_run_name)
set(AX(:),'Position',[0.13,0.71,0.3346590909090909,0.1610507974588962])
set(BX,'Position',[0.5703409090909091,0.71,0.334659090909091,0.1610507974588962])

if exist(strrep(diagsteps.filearr(2:end-1),'surf','exf'),'file') % Use EXF
    subplot(3,2,[3 5]);
    plot(diag_yrs,global_fwflx,'k',...
        diag_yrs(end-length(global_exfevap)+1:end),global_exfevap,'r',...
        diag_yrs(end-length(global_exfevap)+1:end),global_exfprec,'b',...
        diag_yrs(end-length(global_exfevap)+1:end),global_exfroff,'g');
    if isvar('ctim')
        hold on
        ylim=get(gca,'YLim');
        plot([ctim(end),ctim(end)],ylim,'k--');
    end
    set(gca,'XLim',xlim);xlabel(['Time [' diagsteps.timeunit ']']);ylabel('[kg/m2/yr]');
    title({['Drift of model surface FW forcing: \color{black}\DeltaFWflx = ',num2str(dfwf,3),...
        ', \color{red}\DeltaexfEvap = ',num2str(devap,3)];[
        ', \color{blue}\DeltaexfPrecip = ',num2str(dprec,3),...
        ', \color{green}\DeltaexfRoff = ',num2str(droff,3)]})
    
    subplot(3,2,[4 6]);
    plot(diag_yrs,global_qnet,'k',...
        diag_yrs(end-length(global_exfevap)+1:end),global_oceqsw,'y',...
        diag_yrs(end-length(global_exfevap)+1:end),global_lwflux,'r',...
        diag_yrs(end-length(global_exfevap)+1:end),global_latent,'b',...
        diag_yrs(end-length(global_exfevap)+1:end),global_sensib,'g');
    if isvar('ctim')
        hold on
        ylim=get(gca,'YLim');
        plot([ctim(end),ctim(end)],ylim,'k--');
    end
    set(gca,'XLim',xlim);xlabel(['Time [' diagsteps.timeunit ']']);ylabel('[W m^{-2}]');
    title({['Drift of model surface heat forcing: \color{black}\DeltaQnet = ',num2str(dqnet,3),...
        ', \color{yellow}\DeltaexfSWnet = ',num2str(dqsw,3)];[
        ', \color{red}\DeltaexfLWnet = ',num2str(dqlw,3),...
        ', \color{blue}\Deltaexflatent = ',num2str(dlat,3),...
        ', \color{green}\Deltaexfsensib = ',num2str(dsen,3)]})
else % USE SURFDIAG
    subplot(3,2,[3 5]);
    plot(diag_yrs,global_fwflx,'b',diag_yrs,global_sflux,'g');
    legvec1={'oceFWflux ';'SFLUX'};
    titlevec1=['Drift of model surface FW forcing: \color{blue}\DeltaoceFWflx (kg m^{-2} yr^{-1}) = ',num2str(dfwf,3),','];
    hold on
    titlevec2=['\color{green}\DeltaSFLUX = ',num2str(dsf,3),', '];
    if isvar('global_srelax'); plot(diag_yrs,global_srelax,'c');titlevec2=[titlevec2,'\color{cyan}\DeltaSRELAX = ',num2str(dsr,3)]; end
    if isvar('global_ocesflux'); plot(diag_yrs,global_ocesflux,'m');titlevec2=[titlevec2,', \color{magenta}\DeltaoceSflux = ',num2str(dosf,3)]; end
    if isvar('global_surforcs');
        plot(diag_yrs,global_surforcs,'r','LineWidth',2);titlevec2=[titlevec2,', \color{red}\DeltasurForcS = ',num2str(dsfs,3)];
        plot(diag_yrs,global_fwflx.*35,'b--');titlevec2=[titlevec2,', \color{blue}PmEpR*35'];
    end
    set(gca,'XLim',xlim);xlabel(['Time [' diagsteps.timeunit ']']);ylabel('[g m^{-2} yr^{-1}]')
    if isvar('ctim')
        hold on
        ylim=get(gca,'YLim');
        plot([ctim(end),ctim(end)],ylim,'k--');
    end
    titlevec=vertcat({titlevec1},{titlevec2});title(titlevec);
    
    subplot(3,2,[4 6]);
    plot(diag_yrs,global_qnet,'b',diag_yrs,global_tflux,'g');
    titlevec1=['Drift of model surface heat forcing: \color{blue}\DeltaoceQnet = ',num2str(dqnet,3),', \color{green}\DeltaTFLUX = ',num2str(dtf,3)];
    hold on
    titlevec2=char(0);
    if isvar('global_trelax'); plot(diag_yrs,global_trelax,'c');titlevec2=[titlevec2,'\color{cyan}\DeltaTRELAX = ',num2str(dtr,3)]; end
    if isvar('global_oceqsw'); plot(diag_yrs,global_oceqsw,'r','LineWidth',2);titlevec2=[titlevec2,', \color{red}\DeltaoceQsw = ',num2str(dqsw,3)]; end
    if isvar('global_ocefreez'); plot(diag_yrs,global_ocefreez,'m');titlevec2=[titlevec2,', \color{magenta}\DeltaoceFreez = ',num2str(dof,3)];end
    if isvar('global_surforct'); plot(diag_yrs,global_surforct,'r','LineWidth',2);titlevec2=[titlevec2,', \color{red}\DeltasurForcT = ',num2str(dsft,3)]; end
    set(gca,'XLim',xlim);xlabel(['Time [' diagsteps.timeunit ']']);ylabel('[W m^{-2}]');
    if isvar('ctim')
        hold on
        ylim=get(gca,'YLim');
        plot([ctim(end),ctim(end)],ylim,'k--');
    end
    titlevec=vertcat({titlevec1},{titlevec2(2:end)});title(titlevec);
end

orient landscape
print -dpsc diag_drift1.ps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURE 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
subplot(221);
[AX,H1,H2] = plotyy(tave_yrs,global_alk,tave_yrs,global_dic);
%set(get(AX(1),'Ylabel'),'String','Alk [mol eq m^{-3}]')
%set(get(AX(2),'Ylabel'),'String','DIC [mol C m^{-3}]')

set(AX,'XLim',xlim);
set(AX(1),'YLim',[min(global_alk)-0.0001 max(global_alk)+0.0001],'YTickMode','auto');
set(AX(2),'YLim',[min(global_dic)-0.001 max(global_dic)+0.001],'YColor','r','YTickMode','auto');
set(H2,'Color','r');
if isvar('ctim')
   hold on
   ylim=get(gca,'YLim');
   plot([ctim(end),ctim(end)],ylim,'k--');
end
titlevec={['Drift of \color{blue}Mean Alkalinity [mol eq m^{-3}, \Deltaalk = ',num2str(dalk,3),']'];...
    ['\color{black}and \color{red}Mean DIC [mol C m^{-3}, \DeltaDIC = ',num2str(ddic,3),']']}; 
title(titlevec)
%xlabel(['Time [' timeunit ']']);

subplot(222);
if isvar('global_fe');
    [BX,H1,H2] = plotyy(tave_yrs,global_po4,tave_yrs,global_fe);
    %set(get(AX(1),'Ylabel'),'String','PO_{4}^{-3} [mol P m^{-3}]')
    %set(get(AX(2),'Ylabel'),'String','Fe [mol Fe m^{-3}]')
    
    set(BX,'XLim',xlim,'YTickMode','auto');
    set(BX(2),'YColor','r');
    set(H2,'Color','r');
   
    if nc_isvar(strrep(tavesteps.filearr(2:end-1),'tave','ptr_tave'),'ligand')
        hold on
        H3=plot(tavesteps.tim./1000,global_ligand,'g');
        set(BX(2),'YLimMode','auto','YTickMode','auto');
    end
    
    titlevec={['Drift of \color{blue}Mean Phosphate [mol P m^{-3}, \DeltaPO_{4} = ',num2str(dpo4,3),']'];...
        ['\color{black}and \color{red}Fe [mol Fe m^{-3}, \DeltaFe = ',num2str(dfe,3),']']}; 
else
    plot(tim,global_po4);set(gca,'XLim',xlim);
    ylabel('PO_{4}^{-3} [ mol P m^{-3} yr^{-1}]');
    titlevec=['Drift of \color{blue}Mean Phosphate [mol P m^{-3}, \DeltaPO_{4} = ',num2str(dpo4,3),']'];
end
if isvar('ctim')
   hold on
   ylim=get(gca,'YLim');
   plot([ctim(end),ctim(end)],ylim,'k--');
end
title(titlevec); 
%xlabel(['Time [' timeunit ']']);

subplot(223);
[CX,H1,H2]=plotyy(tave_yrs,tco2,tave_yrs,tpo4);
set(H2,'Color','r')
set(CX(2),'YColor','r','YLim',[(mean(tpo4)-mean(tpo4)./1000),(mean(tpo4)+mean(tpo4)./1000)],'YTickMode','auto','YTickLabelMode','auto')
%set(get(AX(1),'Ylabel'),'String','Carbon Content [mol C]');
%set(get(AX(2),'Ylabel'),'String','DOP+PO_{4} [mol P]','Color','r');
if isvar('tallco2');
   axes(CX(1))
   hold on
   H3=plot(tave_yrs,tallco2,'g');
   set(CX(1),'YLimMode','auto','YTickMode','auto');
   axes(CX(2))
   titlevec={['Drift of \color{blue}Total Oceanic Carbon [mol C, \DeltaC = ',num2str(dtco2,3),']'];...
       ['\color{black}and \color{green}Total Atm+Ocean Carbon [mol C, \DeltaC = ',num2str(dtallco2,3),']'];...
       ['\color{black}and \color{red}Total Phosphorus [mol P, \DeltaP = ',num2str(dpo4,3),']']};
else
    titlevec={['Drift of \color{blue}Total Oceanic Carbon [mol C, \DeltaC = ',num2str(dtco2,3),']'];...
       ['\color{black}and \color{red}Total Phosphorus [mol P, \DeltaP = ',num2str(dpo4,3),']']};
end

if isvar('global_vflux');
	axes(CX(1))
	hold on
	H4=plot(diag_yrs(end-length(global_vflux_anom)+1:end),global_vflux_anom,'b--');
    set(CX(1),'YLimMode','auto','YTickMode','auto');
    axes(CX(2))
    titlevec=vertcat(titlevec,['\color{black}and \color{blue}Virtual Flux of DIC [mol C, \DeltaC = ',num2str(dvflux,3),']']);
end

if isvar('ctim')
   hold on
   ylim=get(gca,'YLim');
   plot([ctim(end),ctim(end)],ylim,'k--');
end
set(CX,'XLim',xlim); 

xlabel(['Time [' tavesteps.timeunit ']']);
title(titlevec);

subplot(224);
plot(tave_yrs,total_bio,'g');set(gca,'XLim',xlim);
if isvar('ctim')
   hold on
   ylim=get(gca,'YLim');
   plot([ctim(end),ctim(end)],ylim,'k--');
   set(gca,'YLim',ylim);
end
title(['Drift of \color{green}Total Primary Production [GtC yr^{-1}, \DeltaPP = ',num2str(dbio,3),']']); 
xlabel(['Time [' tavesteps.timeunit ']']);
%ylabel('PP [GtC yr^{-1}]');

suptitle(attributes.the_run_name)
orient landscape

set(AX(2),'Position',get(AX(1),'Position'));
set(BX(2),'Position',get(BX(1),'Position'));
set(CX(2),'Position',get(CX(1),'Position'));

print -dpsc diag_drift2.ps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURE 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
subplot(221);
[AX,H1,H2]=plotyy(tave_yrs,global_flxco2,tave_yrs,total_flxco2);
set(H2,'Color','r'); set(AX(2),'YColor','r')
%set(get(AX(1),'Ylabel'),'String','CO_{2} [ mol C m^{-2} yr^{-1}]');
%set(get(AX(2),'Ylabel'),'String','CO_{2} [ mol C yr^{-1}]','Color','r');
set(AX,'XLim',xlim);
if isvar('ctim')
   hold on
   ylim=get(AX(2),'YLim');
      plot([ctim(end),ctim(end)],ylim,'k--');
end
titlevec={['Drift of \color{blue}Mean Air-Sea CO_{2} flux [mol C m^{-2} yr^{-1}, \Deltaflx = ',num2str(dflxco2,3),']'];...
    ['\color{black}and \color{red}Molar Equivalent [mol C yr^{-1}, \DeltaCflx = ',num2str(dtflxco2,3),']']}; 
title(titlevec); xlabel(['Time [' tavesteps.timeunit ']']);

subplot(222);
plot(tave_yrs,global_pco2);set(gca,'XLim',xlim)
if isvar('ctim')
   hold on
   ylim=get(gca,'YLim');
   plot([ctim(end),ctim(end)],ylim,'k--');
end
titlevec={['Drift of \color{blue}Mean Oceanic pCO_{2} [\mu atm, \DeltapCO_{2} = ',num2str(dpco2,3),']']}; 
title(titlevec); xlabel(['Time [' tavesteps.timeunit ']']);
%ylabel('pCO_{2} [\mu atm]'); 

suptitle(attributes.the_run_name)

if exist('atm_co2','var')
    subplot(2,2,3:4);

    if exist('ldiff','var')
        extension=[tave_yrs(1),atm_co2(1,2),atm_co2(1,3)];
        atm_co2=[extension;atm_co2];
    end

    if strcmp(timeunit,'kyrs')
        atm_co2(:,1)=atm_co2(:,1)./1000;
    end
        
    titlevec={['Drift of \color{blue}Atmospheric Carbon Concentration [mol, \DeltaCO_2 = ',num2str(datmc,3),']'];...
        ['\color{black}and \color{red}Atmospheric pCO_2 [uatm, \DeltaCO_2 = ',num2str(datmp,3),']']}; 
    
    if dicint==4 && max(atm_co2(:,4))~=0
       [BX,H1,H2,H3]=plotyy(atm_co2(:,1),atm_co2(:,2)./1e16,atm_co2(:,1),atm_co2(:,3)*1e6,atm_co2(:,1),(atm_co2(:,4).*(60*60*24*360))./1e16);
       set(H3,'Color','g'); set(BX(3),'YColor','r');
       set(get(BX(3),'Ylabel'),'String','Atmospheric C emissions (x10^{16} mol/yr)');
       titlevec=vertcat(titlevec,' \color{black}and \color{green}Carbon Emmisions [mol/yr');
    else
       [BX,H1,H2]=plotyy(atm_co2(:,1),atm_co2(:,2)./1e16,atm_co2(:,1),atm_co2(:,3)*1e6);
    end
    
    set(H2,'Color','r'); set(BX(2),'YColor','r');
    set(get(BX(1),'Ylabel'),'String','Atmospheric C conc (x10^{16} moles)');
    set(get(BX(2),'Ylabel'),'String','Atmospheric pCO_2 (uatm)','Color','r');
    set(BX,'XLim',xlim);
    title(titlevec);
    xlabel(['Time [' tavesteps.timeunit ']']);
    if isvar('ctim')
        hold on
        ylim=get(gca,'YLim');
        plot([ctim(end),ctim(end)],ylim,'k--');
    end
end

orient landscape
set(AX(2),'Position',get(AX(1),'Position'));
set(BX(2),'Position',get(BX(1),'Position'));
print -dpsc diag_drift3.ps

% Restore the original variables
if exist([filepath,'/control.mat'],'file') && ( ~strcmp(attributes.the_run_name,'Control Run') || ~strcmp(attributes.the_run_name,'cntrl') );
    if strcmp(tavesteps.timeunit,'kyrs')
        tavesteps.timeunit='yrs';
    end
    load tmpbackup.mat
    delete tmpbackup.mat
end
%close all
save mitloadglobal.mat global_* total_* d* tdp tpo4 *co2 diagsteps tavesteps