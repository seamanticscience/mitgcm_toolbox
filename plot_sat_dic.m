function [ varargout ] = plot_sat_dic(grd)
%plot_sat_dic loads and plots carbon components calculated online with
% MITgcm.

%% Load variables
if nargin==0 || ~exist('grd','var')
    grd=mit_loadgrid;
end

area = grd.rac.*grd.hfacc(:,:,1) ;
grd=mit_oceanmasks(grd);

if ~exist('diagsteps','var') || ~exist('ptrsteps','var')
    diagsteps=mit_timesteps('componentDiag');
    ptrsteps =mit_timesteps('ptr_tave');
end

atm_co2=mit_getdicpco2('data');

%% Do some temporate inventory plots 

kmax=diagsteps.kmax;

global_preg = zeros(kmax,1);
global_po4 = zeros(kmax,1);
global_pstar = zeros(kmax,1);
total_cpre = zeros(kmax,1);
total_csatinit = zeros(kmax,1);
total_cresinit = zeros(kmax,1);
total_csat = zeros(kmax,1);
total_cres = zeros(kmax,1);
total_csatinsitu = nan(kmax,1);
total_cresinsitu = nan(kmax,1);
total_csoft = zeros(kmax,1);
total_ccarb = zeros(kmax,1);
total_cant = zeros(kmax,1);
total_cstar = zeros(kmax,1);
total_cstar_res = zeros(kmax,1);
total_dic = zeros(kmax,1);
total_dop = zeros(kmax,1);

% dicDiags to look for virtual fluxes
global_vflux = zeros(kmax,1);

for k=1:kmax;
    cparts=rdmnc(diagsteps.filearr(2:end-1),'CSAT','CRES','CSATINIT','CRESINIT','CSATSIT','CRESSIT',...
        'CSOFT','CCARB','PREG','AREG',diagsteps.timesteps(k));
    
    ptracers=rdmnc(ptrsteps.filearr(2:end-1),'dic','alk','po4','dop','cpre','apre','ppre',...
        'wm_age','atmpco2',ptrsteps.timesteps(k));
    
    total_dic(k) = nansum(ptracers.dic(:).*grd.volc(:)); % Total moles of inorganic C in ocean (mol C)
    total_dop(k) = nansum(ptracers.dop(:).*grd.volc(:)); % Total moles of organic P in ocean (mol P)
    global_preg(k)= nansum(cparts.PREG(:).*grd.volc(:))/nansum(grd.volc(:));
    global_po4 (k)= nansum(ptracers.po4(:).*grd.volc(:))/nansum(grd.volc(:)); 
    global_pstar(k) = (global_preg(k)./global_po4(k)).*100; % efficiency of soft tissue pump (P*, Ito & Follows, 2005)
    total_csatinit(k) = nansum(cparts.CSATINIT(:).*grd.volc(:)).*12e-15; % Total moles of Csat at initial pCO2 (PgC)
    total_cresinit(k) = nansum(cparts.CRESINIT(:).*grd.volc(:)).*12e-15; % Total moles of Cres at initial pCO2 (PgC)
    total_csat(k)     = nansum(cparts.CSAT(:).*grd.volc(:)).*12e-15; % Total moles of Csat at current pCO2 (PgC)
    total_cres(k)     = nansum(cparts.CRES(:).*grd.volc(:)).*12e-15; % Total moles of Cres at current pCO2 (PgC)
    total_cpre(k)     = nansum(ptracers.cpre(:).*grd.volc(:))*12e-15;
    
    if isvar('cparts.CSATSIT'); 
        total_csatinsitu(k) = nansum(cparts.CSATSIT(:).*grd.volc(:)).*12e-15; % Total moles of Csat at in situ pCO2 (PgC)
        total_cresinsitu(k) = nansum(cparts.CRESSIT(:).*grd.volc(:)).*12e-15; % Total moles of Cres at in situ pCO2 (PgC)
    else
        total_csatinsitu(k) = NaN;
        total_cresinsitu(k) = NaN;
    end
    
    if isvar('ptracers.cpre');
        % Assuming no change in disequilibrium
        cant=ptracers.cpre-cparts.CSATINIT;
        total_cant(k)= nansum(cant(:).*grd.volc(:))*12e-15;
    else
        cant=ptracers.dic.*0;
        total_cant(k)=0;
    end
    
    if isvar('cparts.CSATSIT') && isvar('cparts.CRESSIT');
        cstar     = cparts.CSATSIT-cparts.CSATINIT;
        cstar_res = cparts.CRESINIT+cparts.CRESSIT;
        
        total_cstar(k)     = nansum(cstar(:).*grd.volc(:))*12e-15;
        total_cstar_res(k) = nansum(cstar_res(:).*grd.volc(:))*12e-15;
    else
        cstar=ptracers.dic.*0;
        cstar_res=ptracers.dic.*0;
        total_cstar(k) = 0;
        total_cstar_res(k) = 0;
    end
    
    total_csoft(k) = nansum(cparts.CSOFT(:).*grd.volc(:)).*12e-15; % total moles of Csoft (PgC)
    total_ccarb(k) = nansum(cparts.CCARB(:).*grd.volc(:)).*12e-15; % Total moles of Ccarb (PgC)
    
    if exist(strrep(diagsteps.filearr(2:end-1),'componentDiag','dic_surfDiag'),'file')
       dicDiag=rdmnc(strrep(diagsteps.filearr(2:end-1),'componentDiag','dic_surfDiag'),'DICVCFLX',diagsteps.timesteps(k));
    else
       dicDiag=rdmnc(strrep(diagsteps.filearr(2:end-1),'componentDiag','dicDiag'),'DICVCFLX',diagsteps.timesteps(k));
    end
    
    if isvar('dicDiag.DICVCFLX'); % diagnostic is im mol C/m3/s in surface layer
       global_vflux(k) = nansum((dicDiag.DICVCFLX(:).*area(:).*31104000)/grd.dz(1));
    end    
end

clear dicDiag

diagfreq=mit_getparm('data.diagnostics','frequency(2)');
% % if exist('/Volumes/PhD_Data/controlrun/control_kpp_10k.mat','file') 
%      load /Volumes/PhD_Data/controlrun/control_kpp_10k cglobal_vflux
% %     global_vflux_cntrl=cglobal_vflux(end);
% % elseif exist('/Volumes/Postdoc_Data/controlrun/control_kpp_10k.mat','file')
% %     load /Volumes/Postdoc_Data/controlrun/control_kpp_10k cglobal_vflux
% %     global_vflux_cntrl=cglobal_vflux(end);
% % else
%     global_vflux_cntrl=cglobal_vflux(end);
% %end
% 
% % Since this will be compared with DIC+RCP*DOP, which integrates the
% % loss/gain of carbon by the ocean, should integrate the virtual flux.
 global_vflux_anom=cumsum(global_vflux-global_vflux(1)).*(diagfreq./31104000)*12e-15; % because of 10 year averages. % control value
 dvflux=global_vflux_anom(end);

tco2=(total_dic+(total_dop.*117)).*12e-15; % Oceanic Carbon conservation (PgC)
tallco2=tco2+(atm_co2(2:end,2)*12e-15); % Ocean-atmosphere CO2 content (PgC)

diagyrs=diagsteps.tim;

% This is for the vertical line
ctim=atm_co2(1,1)./1000;

% Prepend control fields if available
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

if exist([filepath,'/controlrun/cntrl_carbon_components.mat'],'file')
    load([filepath,'/controlrun/cntrl_carbon_components.mat'])
    names=fieldnames(cntrl_timeseries);
    for i=1:length(names)-1
        eval([names{i},'=[cntrl_timeseries.',names{i},';',names{i},'];'])
    end
    
    atm_co2=[cntrl_timeseries.atm_co2(1:end-1,:);atm_co2];
end

diagyrs=diagyrs./1000;
xlim = [diagyrs(1),diagyrs(end)];
idx=length(cntrl_timeseries.diagyrs);

datmc=(atm_co2(end,2)-atm_co2(max(1,idx),2))*12e-15;
datmp=(atm_co2(end,3)-atm_co2(max(1,idx),3))*1e6;
dtco2 = tco2(end) - tco2(idx);
dtallco2 = tallco2(end) - tallco2(idx);
dpstar= (global_pstar(end)-global_pstar(idx));
dcsat = total_csat(end) - total_csat(idx);
dcres = total_cres(end) - total_cres(idx);
dcsatinit = total_csatinit(end) - total_csatinit(idx);
dcresinit = total_cresinit(end) - total_cresinit(idx);
dcsatinsitu = total_csatinsitu(end) - total_csatinsitu(idx);
dcresinsitu = total_cresinsitu(end) - total_cresinsitu(idx);
dcant = total_cant(end) - total_cant(idx);
dcstar= total_cstar(end)- total_cstar(idx);
dcstar_res=total_cstar_res(end)- total_cstar_res(idx);
dcsoft = total_csoft(end) - total_csoft(idx);
dccarb = total_ccarb(end) - total_ccarb(idx);

%% Assign output variables

timeseries=struct(...
    'diagyrs',diagyrs,...
    'global_preg', global_preg,...
    'global_po4',global_po4,...
    'global_pstar',global_pstar,...
    'total_cpre',total_cpre,...
    'total_csatinit',total_csatinit,...
    'total_cresinit',total_cresinit,...
    'total_csat',total_csat,...
    'total_cres',total_cres,...
    'total_csatinsitu',total_csatinsitu,...
    'total_cresinsitu',total_cresinsitu,...
    'total_csoft',total_csoft,...
    'total_ccarb',total_ccarb,...
    'total_cant',total_cant,...
    'total_cstar',total_cstar,...
    'total_cstar_res',total_cstar_res,...
    'total_dic',total_dic,...
    'total_dop',total_dop,...
    'global_vflux_anom',global_vflux_anom,...
    'tco2',tco2,...
    'tallco2',tallco2,...
    'atm_co2',atm_co2);

fields=struct(...
    'cant',cant,...
    'cstar',cstar,...
    'cstar_res',cstar_res);

varargout={cparts,ptracers,timeseries,fields};

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURE 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
subplot(221);
[h] = plot(diagyrs,total_csat,'r',diagyrs,total_csatinit,'b');
           
titlevec={['Drift of Saturated C components [PgC]'];...
    [];...
    ['\color{red}\DeltaC_{sat} = ',num2str(round(dcsat),4),'\color{black}, \color{blue}\DeltaC_{sat}^{init} = ',num2str(round(dcsatinit),4)]};

hold on

if dcsatinsitu ~=0;
    plot(diagyrs,total_csatinsitu,'g')
    titlevec(3)={['\color{red}\DeltaC_{sat} = ',num2str(round(dcsat),4),'\color{black}, \color{blue}\DeltaC_{sat}^{init} = ',num2str(round(dcsatinit),4),...
    '\color{black}, and \color{green}\DeltaC_{sat}^{insitu} = ',num2str(round(dcsatinsitu),4)]};
end

% if dcant ~= 0;
%     plot(diagyrs,total_cant,'m')
%     titlevec(3)={['\color{red}\DeltaC_{sat} = ',num2str(round(dcsat),4),'\color{black}, \color{blue}\DeltaC_{sat}^{init} = ',num2str(round(dcsatinit),4),...
%     '\color{black}, \color{green}\DeltaC_{sat}^{insitu} = ',num2str(round(dcsatinsitu),4),'\color{black}, and \color{magenta}\DeltaC_{ant} = ',num2str(round(dcant),4)]};
%     plot(xlim,[0 0],'k--');
% end

if isvar('ctim')
   hold on
   ylim=get(gca,'YLim');
   H3=plot([ctim(end),ctim(end)],ylim,'k--');
end

set(gca,'XLim',xlim);
title(titlevec)
xlabel('Time [kyrs]'); ylabel('Globally integrated C [PgC]') 

subplot(222);
[BX,H1,H2]=plotyy(atm_co2(:,1)./1000,atm_co2(:,2)*12e-15,atm_co2(:,1)./1000,atm_co2(:,3)*1e6);
titlevec={'Drift of total carbon reservoir content [PgC]';...
        ['\color{blue}Atmosphere [mol, \DeltaCO_2 = ',num2str(round(datmc),4),']'];...
        ['\color{black}and \color{red}Atmospheric pCO_2 [uatm, \DeltaCO_2 = ',num2str(round(datmp),4),']']}; 
set(H2,'Color','r'); set(BX(2),'YColor','r');
set(BX(1),'YColor','k');
set(get(BX(1),'Ylabel'),'String','Globally integrated C [PgC]');
set(get(BX(2),'Ylabel'),'String','Atmospheric pCO_2 (uatm)','Color','r');

axes(BX(1))
hold on
H3=plot(diagyrs,tco2,'c');
H4=plot(diagyrs,tallco2,'m');
set(BX,'YLimMode','auto','YTickMode','auto');
axes(BX(2))

titlevec={'Drift of total carbon reservoir content [PgC]';...
        ['\color{red}Atm pCO_2 [\DeltaCO_2 = ',num2str(round(datmp),4),'], \color{blue}Atm [\DeltaCO_2 = ',num2str(round(datmc),4),'], '];...
        ['\color{cyan}Ocean [\DeltaC = ',num2str(round(dtco2),4),'], \color{black}and \color{magenta}Atm+Ocean[\DeltaC = ',num2str(round(dtallco2),4),']']};

if isvar('global_vflux');
	axes(BX(1))
	hold on
	H5=plot(diagyrs,global_vflux_anom,'g--');
    set(BX(1),'YLimMode','auto','YTickMode','auto');
    axes(BX(2))
    titlevec=vertcat(titlevec,['\color{black}and \color{green}Virtual DIC Flux [\DeltaC = ',num2str(round(dvflux),4),']']);
end

set(BX(1),'YLimMode','auto','YTickMode','auto');
set(BX,'XLim',xlim);

if isvar('ctim')
   hold on
   ylim=get(gca,'YLim');
   H6=plot([ctim(end),ctim(end)],ylim,'k--');
end

title(titlevec);
xlabel('Time [kyrs]');

subplot(223);
[h] = plot(diagyrs,total_cres,'r',diagyrs,total_cresinit,'b');
           
titlevec={['\color{red}\DeltaC_{res} = ',num2str(round(dcres),4),'\color{black}, \color{blue}\DeltaC_{res}^{init} = ',num2str(round(dcresinit),4)]};

hold on
plot(xlim,[0 0],'k--');

if dcresinsitu ~=0;
    plot(diagyrs,total_cresinsitu,'g')
    titlevec={['\color{red}\DeltaC_{res} = ',num2str(round(dcres),4),'\color{black}, \color{blue}\DeltaC_{res}^{init} = ',num2str(round(dcresinit),4),...
    '\color{black}, and \color{green}\DeltaC_{res}^{insitu} = ',num2str(round(dcresinsitu),4)]};
end

if isvar('cant') && isvar('cstar');
    plot(diagyrs,total_cant,'c')
    plot(diagyrs,total_cstar,'m')
    plot(diagyrs,total_cstar_res,'y')
        titlevec=vertcat(titlevec,...
        ['\color{cyan}\DeltaC_{ant} = ',num2str(round(dcant),4),'\color{black}, \color{magenta}\DeltaC_{sat}^{*} = ',num2str(round(dcstar),4),...
    '\color{black}, and \color{yellow}\DeltaC_{res}^{*} = ',num2str(round(dcstar_res),4)]);
end

if isvar('ctim')
   hold on
   ylim=get(gca,'YLim');
   H3=plot([ctim(end),ctim(end)],ylim,'k--');
end

set(gca,'XLim',xlim);
title(titlevec)
xlabel('Time [kyrs]'); ylabel('Globally integrated C [PgC]') 

subplot(224);
[CX,H1,H2]=plotyy(diagyrs,total_csoft,diagyrs,global_pstar);
axes(CX(1))
hold on
H3=plot(diagyrs,total_ccarb,'b');
set(CX(1),'YLimMode','auto','YTickMode','auto');
axes(CX(2))
   
titlevec={['Drift of Biogenic C components and Bio efficiency [PgC]'];...
    [];...
    ['\color{red}\DeltaC_{soft} = ',num2str(round(dcsoft),4),'\color{black}, \color{blue}\DeltaC_{carb} = ',num2str(round(dccarb),4),'\color{black}, and \color{green}\DeltaP* = ',num2str(dpstar,2)]}; 
set(H1,'Color','r'); 
set(CX(1),'YColor','k');
set(H2,'Color','g'); set(CX(2),'YColor','k');
set(get(CX(1),'Ylabel'),'String','Globally integrated C [PgC]','Color','k');
set(get(CX(2),'Ylabel'),'String','Biological efficiency [%]','Color','k');
set(CX(2),'YLim',[10 50],'YTickMode','auto')
set(CX,'XLim',xlim);

if isvar('ctim')
   hold on
   ylim=get(gca,'YLim');
   H3=plot([ctim(end),ctim(end)],ylim,'k--');
end

title(titlevec);
xlabel('Time [kyrs]');

suptitle(cparts.attributes.global.the_run_name)

%set(h(2),'Position',get(h(1),'Position'));
%set(BX(2),'Position',get(BX(1),'Position'));
%set(CX(2),'Position',get(CX(1),'Position'));
orient landscape
print -dpsc2 carbon_component_drift.ps
fixpslinestyle('carbon_component_drift.ps')

%% Do some vertical section plots
% Last timestep is already loaded from previous section.

if isvar('cparts.CSATSIT');
    % num of subfigs
    ss=3;
else
    ss=2;
end

%% CSAT at current atm pco2 and initial atm pco2
figure
subplot(ss,2,1)
contourf(grd.latc,-grd.zc,squeeze(nanmean(cparts.CSAT.*grd.atlantic_hfacc(:,:,:),1))'.*1000,[2000:20:2500])
colormap(bluewhitered(length(2000:20:2500)-1))
caxis([2000 2500])
colorbar('YTick',2000:100:2500)
set(gca,'Color','k');
title('Atlantic C_{sat} at current pCO_2 [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')

subplot(ss,2,2)
contourf(grd.latc,-grd.zc,squeeze(nanmean(cparts.CSAT.*grd.pacific_hfacc(:,:,:),1))'.*1000,[2000:20:2500])
colormap(bluewhitered(length(2000:20:2500)-1))
caxis([2000 2500])
colorbar('YTick',2000:100:2500)
set(gca,'Color','k');
title('Pacific C_{sat} at current pCO_2 [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')

subplot(ss,2,3)
contourf(grd.latc,-grd.zc,squeeze(nanmean(cparts.CSATINIT.*grd.atlantic_hfacc(:,:,:),1))'.*1000,[2000:20:2500])
colormap(bluewhitered(length(2000:20:2500)-1))
caxis([2000 2500])
colorbar('YTick',2000:100:2500)
set(gca,'Color','k');
title('Atlantic C_{sat} at initial pCO_2 [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')

subplot(ss,2,4)
contourf(grd.latc,-grd.zc,squeeze(nanmean(cparts.CSATINIT.*grd.pacific_hfacc(:,:,:),1))'.*1000,[2000:20:2500])
colormap(bluewhitered(length(2000:20:2500)-1))
caxis([2000 2500])
colorbar('YTick',2000:100:2500)
set(gca,'Color','k');
title('Pacific C_{sat} at initial pCO_2 [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')

if isvar('cparts.CSATSIT');
    subplot(ss,2,5)
    contourf(grd.latc,-grd.zc,squeeze(nanmean(cparts.CSATSIT.*grd.atlantic_hfacc(:,:,:),1))'.*1000,[2000:20:2500])
    colormap(bluewhitered(length(2000:20:2500)-1))
    caxis([2000 2500])
    colorbar('YTick',2000:100:2500)
    set(gca,'Color','k');
    title('Atlantic C_{sat} at in situ pCO_2 [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')
    
    subplot(ss,2,6)
    contourf(grd.latc,-grd.zc,squeeze(nanmean(cparts.CSATSIT.*grd.pacific_hfacc(:,:,:),1))'.*1000,[2000:20:2500])
    colormap(bluewhitered(length(2000:20:2500)-1))
    caxis([2000 2500])
    colorbar('YTick',2000:100:2500)
    set(gca,'Color','k');
    title('Pacific C_{sat} at in situ pCO_2 [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')
end

orient landscape
print -dpsc2 carbon_component_csat.ps

%% CRES at current atm pco2 and initial atm pco2
figure
subplot(ss,2,1)
contourf(grd.latc,-grd.zc,squeeze(nanmean(cparts.CRES.*grd.atlantic_hfacc(:,:,:),1))'.*1000,[-300:20:300])
colormap(bluewhitered(length(-300:20:300)-1))
caxis([-300 300])
colorbar('YTick',-300:100:300)
set(gca,'Color','k');
title('Atlantic C_{res} at current pCO_2 [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')

subplot(ss,2,2)
contourf(grd.latc,-grd.zc,squeeze(nanmean(cparts.CRES.*grd.pacific_hfacc(:,:,:),1))'.*1000,[-300:20:300])
colormap(bluewhitered(length(-300:20:300)-1))
caxis([-300 300])
colorbar('YTick',-300:100:300)
set(gca,'Color','k');
title('Pacific C_{res} at current pCO_2 [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')

subplot(ss,2,3)
contourf(grd.latc,-grd.zc,squeeze(nanmean(cparts.CRESINIT.*grd.atlantic_hfacc(:,:,:),1))'.*1000,[-300:20:300])
colormap(bluewhitered(length(-300:20:300)-1))
caxis([-300 300])
colorbar('YTick',-300:100:300)
set(gca,'Color','k');
title('Atlantic C_{res} at initial pCO_2 [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')

subplot(ss,2,4)
contourf(grd.latc,-grd.zc,squeeze(nanmean(cparts.CRESINIT.*grd.pacific_hfacc(:,:,:),1))'.*1000,[-300:20:300])
colormap(bluewhitered(length(-300:20:300)-1))
caxis([-300 300])
colorbar('YTick',-300:100:300)
set(gca,'Color','k');
title('Pacific C_{res} at initial pCO_2 [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')

if isvar('cparts.CSATSIT');
    subplot(ss,2,5)
    contourf(grd.latc,-grd.zc,squeeze(nanmean(cparts.CRESSIT.*grd.atlantic_hfacc(:,:,:),1))'.*1000,[-300:20:300])
    colormap(bluewhitered(length(-300:20:300)-1))
    caxis([-300 300])
    colorbar('YTick',-300:100:300)
    set(gca,'Color','k');
    title('Atlantic C_{res} at in situ pCO_2 [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')
    
    subplot(ss,2,6)
    contourf(grd.latc,-grd.zc,squeeze(nanmean(cparts.CRESSIT.*grd.pacific_hfacc(:,:,:),1))'.*1000,[-300:20:300])
    colormap(bluewhitered(length(-300:20:300)-1))
    caxis([-300 300])
    colorbar('YTick',-300:100:300)
    set(gca,'Color','k');
    title('Pacific C_{res} at in situ pCO_2 [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')
end

colormap(bluewhitered(length(-300:20:300)-1))
caxis([-300 300])
colorbar('YTick',-300:100:300)
orient landscape
print -dpsc2 carbon_component_cres.ps

%% CANT methods
figure
subplot(ss,2,1)
contourf(grd.latc,-grd.zc,squeeze(nanmean(cant.*grd.atlantic_hfacc(:,:,:),1))'.*1000,[-300:20:300])
colormap(bluewhitered(length(-300:20:300)-1))
caxis([-300 300])
colorbar('YTick',-300:100:300)
set(gca,'Color','k');
title('Atlantic C_{ant} [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')

subplot(ss,2,2)
contourf(grd.latc,-grd.zc,squeeze(nanmean(cant.*grd.pacific_hfacc(:,:,:),1))'.*1000,[-300:20:300])
colormap(bluewhitered(length(-300:20:300)-1))
caxis([-300 300])
colorbar('YTick',-300:100:300)
set(gca,'Color','k');
title('Pacific C_{ant} [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')

if isvar('cparts.CSATSIT');
    subplot(ss,2,3)
    contourf(grd.latc,-grd.zc,squeeze(nanmean(cstar.*grd.atlantic_hfacc(:,:,:),1))'.*1000,[-300:20:300])
    colormap(bluewhitered(length(-300:20:300)-1))
    caxis([-300 300])
    colorbar('YTick',-300:100:300)
    set(gca,'Color','k');
    title('Atlantic C^{*}_{sat} [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')
    
    subplot(ss,2,4)
    contourf(grd.latc,-grd.zc,squeeze(nanmean(cstar.*grd.pacific_hfacc(:,:,:),1))'.*1000,[-300:20:300])
    colormap(bluewhitered(length(-300:20:300)-1))
    caxis([-300 300])
    colorbar('YTick',-300:100:300)
    set(gca,'Color','k');
    title('Pacific C^{*}_{sat} [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')
    
    subplot(ss,2,5)
    contourf(grd.latc,-grd.zc,squeeze(nanmean(cstar_res.*grd.atlantic_hfacc(:,:,:),1))'.*1000,[-300:20:300])
    colormap(bluewhitered(length(-300:20:300)-1))
    caxis([-300 300])
    colorbar('YTick',-300:100:300)
    set(gca,'Color','k');
    title('Atlantic C^{*}_{res} [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')
    
    subplot(ss,2,6)
    contourf(grd.latc,-grd.zc,squeeze(nanmean(cstar_res.*grd.pacific_hfacc(:,:,:),1))'.*1000,[-300:20:300])
    colormap(bluewhitered(length(-300:20:300)-1))
    caxis([-300 300])
    colorbar('YTick',-300:100:300)
    set(gca,'Color','k');
    title('Pacific C^{*}_{res} [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')
end

colormap(bluewhitered(length(-300:20:300)-1))
caxis([-300 300])
colorbar('YTick',-300:100:300)
orient landscape
print -dpsc2 carbon_component_cant.ps

%% CSOFT and CCARB
figure
subplot(2,2,1)
contourf(grd.latc,-grd.zc,squeeze(nanmean(cparts.CSOFT.*grd.atlantic_hfacc(:,:,:),1))'.*1000,[-300:20:300])
colormap(bluewhitered(length(-300:20:300)-1))
caxis([-300 300])
colorbar('YTick',-300:100:300)
set(gca,'Color','k');
title('Atlantic C_{soft} [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')

subplot(2,2,2)
contourf(grd.latc,-grd.zc,squeeze(nanmean(cparts.CSOFT.*grd.pacific_hfacc(:,:,:),1))'.*1000,[-300:20:300])
colormap(bluewhitered(length(-300:20:300)-1))
caxis([-300 300])
colorbar('YTick',-300:100:300)
set(gca,'Color','k');
title('Pacific C_{soft} [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')

subplot(2,2,3)
contourf(grd.latc,-grd.zc,squeeze(nanmean(cparts.CCARB.*grd.atlantic_hfacc(:,:,:),1))'.*1000,[-300:20:300])
colormap(bluewhitered(length(-300:20:300)-1))
caxis([-300 300])
colorbar('YTick',-300:100:300)
set(gca,'Color','k');
title('Atlantic C_{carb} [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')

subplot(2,2,4)
contourf(grd.latc,-grd.zc,squeeze(nanmean(cparts.CCARB.*grd.pacific_hfacc(:,:,:),1))'.*1000,[-300:20:300])
colormap(bluewhitered(length(-300:20:300)-1))
caxis([-300 300])
colorbar('YTick',-300:100:300)
set(gca,'Color','k');
title('Pacific C_{carb} [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')

colormap(bluewhitered(length(-300:20:300)-1))
caxis([-300 300])
colorbar('YTick',-300:100:300)
orient landscape
print -dpsc2 carbon_component_cbio.ps

%% PPRE and PREG
figure
subplot(2,2,1)
contourf(grd.latc,-grd.zc,squeeze(nanmean(ptracers.ppre.*grd.atlantic_hfacc(:,:,:),1))'.*1000,[0:0.1:2])
colormap(parula(length(0:0.1:2)-1))
caxis([0 2])
colorbar('YTick',0:0.5:2)
set(gca,'Color','k');
title('Atlantic P_{pre} [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')

subplot(2,2,2)
contourf(grd.latc,-grd.zc,squeeze(nanmean(ptracers.ppre.*grd.pacific_hfacc(:,:,:),1))'.*1000,[0:0.1:2])
colormap(parula(length(0:0.1:2)-1))
caxis([0 2])
colorbar('YTick',0:0.5:2)
set(gca,'Color','k');
title('Pacific P_{pre} [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')

subplot(2,2,3)
contourf(grd.latc,-grd.zc,squeeze(nanmean(cparts.PREG.*grd.atlantic_hfacc(:,:,:),1))'.*1000,[0:0.1:2])
colormap(parula(length(0:0.1:2)-1))
caxis([0 2])
colorbar('YTick',0:0.5:2)
set(gca,'Color','k');
title('Atlantic P_{reg} [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')

subplot(2,2,4)
contourf(grd.latc,-grd.zc,squeeze(nanmean(cparts.PREG.*grd.pacific_hfacc(:,:,:),1))'.*1000,[0:0.1:2])
colormap(parula(length(0:0.1:2)-1))
caxis([0 2])
colorbar('YTick',0:0.5:2)
set(gca,'Color','k');
title('Pacific P_{reg} [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')

orient landscape
print -dpsc2 carbon_component_phos.ps

%% CPRE and CREG
figure
subplot(2,2,1)
contourf(grd.latc,-grd.zc,squeeze(nanmean(ptracers.cpre.*grd.atlantic_hfacc(:,:,:),1))'.*1000,[2000:20:2500])
colormap(bluewhitered(length(2000:20:2500)-1))
caxis([2000 2500])
colorbar('YTick',2000:100:2500)
set(gca,'Color','k');
title('Atlantic C_{pre} [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')

subplot(2,2,2)
contourf(grd.latc,-grd.zc,squeeze(nanmean(ptracers.cpre.*grd.pacific_hfacc(:,:,:),1))'.*1000,[2000:20:2500])
colormap(bluewhitered(length(2000:20:2500)-1))
caxis([2000 2500])
colorbar('YTick',2000:100:2500)
set(gca,'Color','k');
title('Pacific C_{pre} [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')

subplot(2,2,3)
contourf(grd.latc,-grd.zc,squeeze(nanmean((ptracers.dic-ptracers.cpre).*grd.atlantic_hfacc(:,:,:),1))'.*1000,[-300:20:300])
colormap(bluewhitered(length(-300:20:300)-1))
caxis([-300 300])
colorbar('YTick',-300:100:300)
set(gca,'Color','k');
title('Atlantic C_{reg} [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')
colormap(bluewhitered(length(-300:20:300)-1))
caxis([-300 300])
colorbar('YTick',-300:100:300)

subplot(2,2,4)
contourf(grd.latc,-grd.zc,squeeze(nanmean((ptracers.dic-ptracers.cpre).*grd.pacific_hfacc(:,:,:),1))'.*1000,[-300:20:300])
colormap(bluewhitered(length(-300:20:300)-1))
caxis([-300 300])
colorbar('YTick',-300:100:300)
set(gca,'Color','k');
title('Pacific C_{reg} [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')
colormap(bluewhitered(length(-300:20:300)-1))
caxis([-300 300])
colorbar('YTick',-300:100:300)

orient landscape
print -dpsc2 carbon_component_cpre.ps

%% APRE and AREG
figure
subplot(2,2,1)
contourf(grd.latc,-grd.zc,squeeze(nanmean(ptracers.apre.*grd.atlantic_hfacc(:,:,:),1))',[2.0:0.01:2.5])
colormap(parula(length(2.0:0.01:2.5)-1))
caxis([2.0 2.4])
colorbar('YTick',2.0:0.1:2.4)
set(gca,'Color','k');
title('Atlantic A_{pre} [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')

subplot(2,2,2)
contourf(grd.latc,-grd.zc,squeeze(nanmean(ptracers.apre.*grd.pacific_hfacc(:,:,:),1))',[2.0:0.01:2.5])
colormap(parula(length(2.0:0.01:2.5)-1))
caxis([2.0 2.4])
colorbar('YTick',2.0:0.1:2.4)
set(gca,'Color','k');
title('Pacific A_{pre} [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')

subplot(2,2,3)
contourf(grd.latc,-grd.zc,squeeze(nanmean(cparts.AREG.*grd.atlantic_hfacc(:,:,:),1))'.*1000,[-10:10:500])
colormap(parula(length(-10:10:500)-1))
caxis([0 500])
colorbar('YTick',0:100:500)
set(gca,'Color','k');
title('Atlantic A_{reg} [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')

subplot(2,2,4)
contourf(grd.latc,-grd.zc,squeeze(nanmean(cparts.AREG.*grd.pacific_hfacc(:,:,:),1))'.*1000,[-10:10:500])
colormap(parula(length(-10:10:500)-1))
caxis([0 500])
colorbar('YTick',0:100:500)
set(gca,'Color','k');
title('Pacific A_{reg} [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')

orient landscape
print -dpsc2 carbon_component_alk.ps

