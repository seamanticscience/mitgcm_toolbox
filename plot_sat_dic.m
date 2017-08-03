function [ varargout ] = plot_sat_dic(grd)
%plot_sat_dic loads and plots carbon components calculated online with
% MITgcm.

%% Load variables
if nargin==0 || ~exist('grd','var')
    grd=mit_loadgrid;
end

surf_area = grd.rac.*grd.hfacc(:,:,1) ;
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
total_trcsat = zeros(kmax,1);
total_trcres = zeros(kmax,1);
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

% variables from the componentDiag file
cvars={};
if nc_isvar(diagsteps.filearr(2:end-1),'CSATINIT'); cvars=[cvars,{'CSATINIT'}]; end
if nc_isvar(diagsteps.filearr(2:end-1),'CRESINIT'); cvars=[cvars,{'CRESINIT'}]; end
if nc_isvar(diagsteps.filearr(2:end-1),'CSATSIT' );  cvars=[cvars,{'CSATSIT' }]; end
if nc_isvar(diagsteps.filearr(2:end-1),'CRESSIT' );  cvars=[cvars,{'CRESSIT' }]; end
if nc_isvar(diagsteps.filearr(2:end-1),'CSAT'    );  cvars=[cvars,{'CSAT'    }]; end
if nc_isvar(diagsteps.filearr(2:end-1),'CRES'    );  cvars=[cvars,{'CRES'    }]; end
if nc_isvar(diagsteps.filearr(2:end-1),'CSOFT'   );  cvars=[cvars,{'CSOFT'   }]; end
if nc_isvar(diagsteps.filearr(2:end-1),'CCARB'   );  cvars=[cvars,{'CCARB'   }]; end
if nc_isvar(diagsteps.filearr(2:end-1),'PREG'    );  cvars=[cvars,{'PREG'    }]; end
if nc_isvar(diagsteps.filearr(2:end-1),'AREG'    );  cvars=[cvars,{'AREG'    }]; end

% variables from the ptr_tave file
pvars={};
if nc_isvar(ptrsteps.filearr(2:end-1),'dic'    ); pvars=[pvars,{'dic'    }]; end
if nc_isvar(ptrsteps.filearr(2:end-1),'alk'    ); pvars=[pvars,{'alk'    }]; end
if nc_isvar(ptrsteps.filearr(2:end-1),'po4'    ); pvars=[pvars,{'po4'    }]; end
if nc_isvar(ptrsteps.filearr(2:end-1),'dop'    ); pvars=[pvars,{'dop'    }]; end
if nc_isvar(ptrsteps.filearr(2:end-1),'cpre'   ); pvars=[pvars,{'cpre'   }]; end
if nc_isvar(ptrsteps.filearr(2:end-1),'apre'   ); pvars=[pvars,{'apre'   }]; end
if nc_isvar(ptrsteps.filearr(2:end-1),'ppre'   ); pvars=[pvars,{'ppre'   }]; end
if nc_isvar(ptrsteps.filearr(2:end-1),'wm_age' ); pvars=[pvars,{'wm_age' }]; end
if nc_isvar(ptrsteps.filearr(2:end-1),'atmpco2'); pvars=[pvars,{'atmpco2'}]; end
if nc_isvar(ptrsteps.filearr(2:end-1),'csat'   ); pvars=[pvars,{'csat'   }]; end
if nc_isvar(ptrsteps.filearr(2:end-1),'cdis'   ); pvars=[pvars,{'cdis'   }]; end

%%    
for k=1:kmax;
    cdiags=rdmnc(diagsteps.filearr(2:end-1),cvars{:},diagsteps.timesteps(k));
    pdiags=rdmnc(ptrsteps.filearr(2:end-1),pvars{:},ptrsteps.timesteps(k));
    
    total_dic(k) = nansum(pdiags.dic(:).*grd.volc(:)); % Total moles of inorganic C in ocean (mol C)
    total_dop(k) = nansum(pdiags.dop(:).*grd.volc(:)); % Total moles of organic P in ocean (mol P)
    global_preg(k)= nansum(cdiags.PREG(:).*grd.volc(:))/nansum(grd.volc(:));
    global_po4 (k)= nansum(pdiags.po4(:).*grd.volc(:))/nansum(grd.volc(:)); 
    global_pstar(k) = (global_preg(k)./global_po4(k)).*100; % efficiency of soft tissue pump (P*, Ito & Follows, 2005)
    
    total_csatinit(k) = nansum(cdiags.CSATINIT(:).*grd.volc(:)).*12e-15; % Total moles of Csat at initial pCO2 (PgC)
    total_cresinit(k) = nansum(cdiags.CRESINIT(:).*grd.volc(:)).*12e-15; % Total moles of Cres at initial pCO2 (PgC)
    
    total_csat(k)     = nansum(cdiags.CSAT(:).*grd.volc(:)).*12e-15; % Total moles of Csat at current pCO2 (PgC)
    total_cres(k)     = nansum(cdiags.CRES(:).*grd.volc(:)).*12e-15; % Total moles of Cres at current pCO2 (PgC)
    
    if isvar('pdiags.csat') && isvar('pdiags.cdis')
        total_trcsat(k)     = nansum(pdiags.csat(:).*grd.volc(:)).*12e-15; % Total moles of the tracer Csat (PgC)
        total_trcres(k)     = nansum(pdiags.cdis(:).*grd.volc(:)).*12e-15; % Total moles of the tracer Cdis (PgC)
    else
        total_trcsat(k) = 0;
        total_trcres(k) = 0;
    end
    
    if isvar('pdiags.cpre');
        total_cpre(k)     = nansum(pdiags.cpre(:).*grd.volc(:))*12e-15;

        % Assuming no change in disequilibrium
        cant=pdiags.cpre-cdiags.CSATINIT;
        total_cant(k)= nansum(cant(:).*grd.volc(:))*12e-15;
        
        % Calculate Cant using actual Csat and Cdis diagnostics
        if isvar('pdiags.csat') && isvar('pdiags.cdis')
           trcant=pdiags.cpre-cdiags.CSATINIT+pdiags.cdis;
           total_trcant(k) = nansum(trcant(:).*grd.volc(:)).*12e-15; % Total moles of the tracer (PgC)
        else
           total_trcant(k) = 0;
           trcant=pdiags.dic.*0;
        end
    else
        total_cpre(k) = 0;
        cant=pdiags.dic.*0;
        total_cant(k)=0;
    end
    
    if isvar('cdiags.CSATSIT') && isvar('cdiags.CRESSIT');
        total_csatinsitu(k) = nansum(cdiags.CSATSIT(:).*grd.volc(:)).*12e-15; % Total moles of Csat at in situ pCO2 (PgC)
        total_cresinsitu(k) = nansum(cdiags.CRESSIT(:).*grd.volc(:)).*12e-15; % Total moles of Cres at in situ pCO2 (PgC)
        
        cstar     = cdiags.CSATSIT-cdiags.CSATINIT;
        cstar_res = cdiags.CRESINIT+cdiags.CRESSIT;
        
        total_cstar(k)     = nansum(cstar(:).*grd.volc(:))*12e-15;
        total_cstar_res(k) = nansum(cstar_res(:).*grd.volc(:))*12e-15;
    else
        total_csatinsitu(k) = 0;
        total_cresinsitu(k) = 0;
        
        cstar=pdiags.dic.*0;
        cstar_res=pdiags.dic.*0;
        total_cstar(k) = 0;
        total_cstar_res(k) = 0;
    end
    
    total_csoft(k) = nansum(cdiags.CSOFT(:).*grd.volc(:)).*12e-15; % total moles of Csoft (PgC)
    total_ccarb(k) = nansum(cdiags.CCARB(:).*grd.volc(:)).*12e-15; % Total moles of Ccarb (PgC)
    
    if exist(strrep(diagsteps.filearr(2:end-1),'componentDiag','dic_surfDiag'),'file')
       dicDiag=rdmnc(strrep(diagsteps.filearr(2:end-1),'componentDiag','dic_surfDiag'),'DICVCFLX',diagsteps.timesteps(k));
    else
       dicDiag=rdmnc(strrep(diagsteps.filearr(2:end-1),'componentDiag','dicDiag'),'DICVCFLX',diagsteps.timesteps(k));
    end
    
    if isvar('dicDiag.DICVCFLX'); % diagnostic is im mol C/m3/s in surface layer
       global_vflux(k) = nansum((dicDiag.DICVCFLX(:).*surf_area(:).*31104000)/grd.dz(1));
    end    
end

%%
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
    elseif exist('/Volumes/My_Passport/','dir')
        filepath='/Volumes/My_Passport';
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
xlim = [floor(diagsteps.tim(1)./1000),ceil(diagsteps.tim(end)/1000)];
idx=length(cntrl_timeseries.diagyrs);

datmc=(atm_co2(end,2)-atm_co2(max(1,idx),2))*12e-15;
datmp=(atm_co2(end,3)-atm_co2(max(1,idx),3))*1e6;
dtco2 = tco2(end) - tco2(idx);
dtallco2 = tallco2(end) - tallco2(idx);
dpstar= (global_pstar(end)-global_pstar(idx));
dcpre=total_cpre(end) - total_cpre(idx);
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

dtrcsat=total_trcsat(end)-total_trcsat(1);
dtrcres=total_trcres(end)-total_trcres(1);
dtrcant=total_trcant(end)-total_trcant(1);

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
    'total_trcsat',total_trcsat,...
    'total_trcres',total_trcres,...
    'total_trcant',total_trcant,...
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
    'cstar_res',cstar_res,...
    'trcant',trcant);

varargout={cdiags,pdiags,timeseries,fields};

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURE 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
subplot(221);
[h] = plot(diagyrs,total_csat,'r',diagyrs,total_csatinit,'b','LineWidth',2);
           
titlevec={['Drift of Saturated C components [PgC]'];...
    [];...
    ['\color{red}\DeltaC_{sat} = ',num2str(round(dcsat),4),'\color{black}, \color{blue}\DeltaC_{sat}^{init} = ',num2str(round(dcsatinit),4)]};

hold on

if dcsatinsitu ~=0;
    plot(diagyrs,total_csatinsitu,'g','LineWidth',2)
    titlevec(3)={['\color{red}\DeltaC_{sat} = ',num2str(round(dcsat),4),'\color{black}, \color{blue}\DeltaC_{sat}^{init} = ',num2str(round(dcsatinit),4),...
    '\color{black}, and \color{green}\DeltaC_{sat}^{insitu} = ',num2str(round(dcsatinsitu),4)]};
end

if dtrcsat ~= 0;
    plot(diagsteps.tim/1000,total_trcsat,'m--','LineWidth',2)
    titlevec(3)={['\color{red}\DeltaC_{sat} = ',num2str(round(dcsat),4),'\color{black}, \color{blue}\DeltaC_{sat}^{init} = ',num2str(round(dcsatinit),4),...
    '\color{black}, \color{green}\DeltaC_{sat}^{insitu} = ',num2str(round(dcsatinsitu),4),'\color{black}, and \color{magenta}\DeltaTrC_{sat} = ',num2str(round(dtrcsat),4)]};
end

% if isvar('ctim')
%    hold on
%    ylim=get(gca,'YLim');
%    H3=plot([ctim(end),ctim(end)],ylim,'k--');
% end

set(gca,'XLim',xlim);
title(titlevec)
xlabel('Time [kyrs]'); ylabel('Globally integrated C [PgC]') 

subplot(222);
[BX,H1,H2]=plotyy(atm_co2(:,1)./1000,atm_co2(:,2)*12e-15,atm_co2(:,1)./1000,atm_co2(:,3)*1e6);
titlevec={'Drift of total carbon reservoir content [PgC]';...
        ['\color{blue}Atmosphere [mol, \DeltaCO_2 = ',num2str(round(datmc),4),']'];...
        ['\color{black}and \color{red}Atmospheric pCO_2 [uatm, \DeltaCO_2 = ',num2str(round(datmp),4),']']}; 
set(H2,'Color','r','LineWidth',2); set(BX(2),'YColor','r');
set(H1,'LineWidth',2);set(BX(1),'YColor','k');
set(get(BX(1),'Ylabel'),'String','Globally integrated C [PgC]');
set(get(BX(2),'Ylabel'),'String','Atmospheric pCO_2 (uatm)','Color','r');

axes(BX(1))
hold on
H3=plot(diagyrs,tco2,'c','LineWidth',2);
H4=plot(diagyrs,tallco2,'m','LineWidth',2);
set(BX,'YLimMode','auto','YTickMode','auto');
axes(BX(2))

titlevec={'Drift of total carbon reservoir content [PgC]';...
        ['\color{red}Atm pCO_2 [\DeltaCO_2 = ',num2str(round(datmp),4),'], \color{blue}Atm [\DeltaCO_2 = ',num2str(round(datmc),4),'], '];...
        ['\color{cyan}Ocean [\DeltaC = ',num2str(round(dtco2),4),'], \color{black}and \color{magenta}Atm+Ocean[\DeltaC = ',num2str(round(dtallco2),4),']']};

% if isvar('global_vflux');
% 	axes(BX(1))
% 	hold on
% 	H5=plot(diagyrs,global_vflux_anom,'g--','LineWidth',2);
%     set(BX(1),'YLimMode','auto','YTickMode','auto');
%     axes(BX(2))
%     titlevec=vertcat(titlevec,['\color{black}and \color{green}Virtual DIC Flux [\DeltaC = ',num2str(round(dvflux),4),']']);
% end

if dcpre ~=0;
    axes(BX(1))
	hold on
    H6=plot(diagyrs,total_cpre,'y','LineWidth',2);
    axes(BX(2))
    titlevec=vertcat(titlevec,['\color{black}and \color{yellow}Preformed Carbon [\DeltaCpre = ',num2str(round(dcpre),4),']']);
end

set(BX(1),'YLimMode','auto','YTickMode','auto');
set(BX,'XLim',xlim);

% if isvar('ctim')
%    hold on
%    ylim=get(gca,'YLim');
%    H6=plot([ctim(end),ctim(end)],ylim,'k--');
% end

title(titlevec);
xlabel('Time [kyrs]');

subplot(223);
[h] = plot(diagyrs,total_cres,'r',diagyrs,total_cresinit,'b','LineWidth',2);
           
titlevec={['\color{red}\DeltaC_{res} = ',num2str(round(dcres),4),'\color{black}, \color{blue}\DeltaC_{res}^{init} = ',num2str(round(dcresinit),4)]};

hold on
plot(xlim,[0 0],'k--');

if dcresinsitu ~=0;
    plot(diagyrs,total_cresinsitu,'g','LineWidth',2)
    titlevec={['\color{red}\DeltaC_{res} = ',num2str(round(dcres),4),'\color{black}, \color{blue}\DeltaC_{res}^{init} = ',num2str(round(dcresinit),4),...
    '\color{black}, and \color{green}\DeltaC_{res}^{insitu} = ',num2str(round(dcresinsitu),4)]};
end

if dtrcres ~=0;
    plot(diagsteps.tim/1000,total_trcres,'k','LineWidth',2)
    titlevec={['\color{red}\DeltaC_{res} = ',num2str(round(dcres),4),'\color{black}, \color{blue}\DeltaC_{res}^{init} = ',num2str(round(dcresinit),4),...
    '\color{black}, and \color{green}\DeltaC_{res}^{insitu} = ',num2str(round(dcresinsitu),4),'\color{black}, and \DeltaTrC_{res} = ',num2str(round(dtrcres),4)]};
end

if isvar('cant') && isvar('cstar');
    plot(diagyrs,total_cant,'c--','LineWidth',1)
    plot(diagsteps.tim/1000,total_trcant,'m--','LineWidth',1)
    plot(diagyrs,total_cstar,'y--','LineWidth',2)
        titlevec=vertcat(titlevec,...
        ['\color{cyan}\DeltaC_{ant} = ',num2str(round(dcant),4),'\color{black}, \color{magenta}\DeltaCantTracer = ',num2str(round(dtrcant),4),...
    '\color{black}, and \color{yellow}\DeltaC^{*} = ',num2str(round(dcstar),4)]);
end

% if isvar('ctim')
%    hold on
%    ylim=get(gca,'YLim');
%    H3=plot([ctim(end),ctim(end)],ylim,'k--');
% end

set(gca,'XLim',xlim);
title(titlevec)
xlabel('Time [kyrs]'); ylabel('Globally integrated C [PgC]') 

subplot(224);
[CX,H1,H2]=plotyy(diagyrs,total_csoft,diagyrs,global_pstar);
axes(CX(1))
hold on
H3=plot(diagyrs,total_ccarb,'b','LineWidth',2);
set(CX(1),'YLimMode','auto','YTickMode','auto');
axes(CX(2))
   
titlevec={['Drift of Biogenic C components and Bio efficiency [PgC]'];...
    [];...
    ['\color{red}\DeltaC_{soft} = ',num2str(round(dcsoft),4),'\color{black}, \color{blue}\DeltaC_{carb} = ',num2str(round(dccarb),4),'\color{black}, and \color{green}\DeltaP* = ',num2str(dpstar,2)]}; 
set(H1,'Color','r','LineWidth',2); 
set(CX(1),'YColor','k');
set(H2,'Color','g','LineWidth',2); set(CX(2),'YColor','k');
set(get(CX(1),'Ylabel'),'String','Globally integrated C [PgC]','Color','k');
set(get(CX(2),'Ylabel'),'String','Biological efficiency [%]','Color','k');
set(CX(2),'YLim',[10 50],'YTickMode','auto')
set(CX,'XLim',xlim);

% if isvar('ctim')
%    hold on
%    ylim=get(gca,'YLim');
%    H3=plot([ctim(end),ctim(end)],ylim,'k--');
% end

title(titlevec);
xlabel('Time [kyrs]');

suptitle(cdiags.attributes.global.the_run_name)

%set(h(2),'Position',get(h(1),'Position'));
%set(BX(2),'Position',get(BX(1),'Position'));
%set(CX(2),'Position',get(CX(1),'Position'));
orient landscape
print -dpsc2 carbon_component_drift.ps
fixpslinestyle('carbon_component_drift.ps')

%% Do some vertical section plots
% Last timestep is already loaded from previous section.

if isvar('cdiags.CSATSIT');
    % num of subfigs
    ss=3;
else
    ss=2;
end

%% CSAT at current atm pco2 and initial atm pco2
figure
subplot(ss,2,1)
contourf(grd.latc,-grd.zc,squeeze(nanmean(cdiags.CSAT.*grd.atlantic_hfacc(:,:,:),1))'.*1000,[2000:20:2500])
colormap(bluewhitered(length(2000:20:2500)-1))
caxis([2000 2500])
colorbar('YTick',2000:100:2500)
set(gca,'Color','k');
title('Atlantic C_{sat} at current pCO_2 [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')

subplot(ss,2,2)
contourf(grd.latc,-grd.zc,squeeze(nanmean(cdiags.CSAT.*grd.pacific_hfacc(:,:,:),1))'.*1000,[2000:20:2500])
colormap(bluewhitered(length(2000:20:2500)-1))
caxis([2000 2500])
colorbar('YTick',2000:100:2500)
set(gca,'Color','k');
title('Pacific C_{sat} at current pCO_2 [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')

subplot(ss,2,3)
contourf(grd.latc,-grd.zc,squeeze(nanmean(cdiags.CSATINIT.*grd.atlantic_hfacc(:,:,:),1))'.*1000,[2000:20:2500])
colormap(bluewhitered(length(2000:20:2500)-1))
caxis([2000 2500])
colorbar('YTick',2000:100:2500)
set(gca,'Color','k');
title('Atlantic C_{sat} at initial pCO_2 [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')

subplot(ss,2,4)
contourf(grd.latc,-grd.zc,squeeze(nanmean(cdiags.CSATINIT.*grd.pacific_hfacc(:,:,:),1))'.*1000,[2000:20:2500])
colormap(bluewhitered(length(2000:20:2500)-1))
caxis([2000 2500])
colorbar('YTick',2000:100:2500)
set(gca,'Color','k');
title('Pacific C_{sat} at initial pCO_2 [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')

if isvar('pdiags.csat')
    subplot(ss,2,5)
    contourf(grd.latc,-grd.zc,squeeze(nanmean(pdiags.csat.*grd.atlantic_hfacc(:,:,:),1))'.*1000,[2000:20:2500])
    colormap(bluewhitered(length(2000:20:2500)-1))
    caxis([2000 2500])
    colorbar('YTick',2000:100:2500)
    set(gca,'Color','k');
    title('Atlantic C_{sat} at in situ pCO_2 (tracer) [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')
    
    subplot(ss,2,6)
    contourf(grd.latc,-grd.zc,squeeze(nanmean(pdiags.csat.*grd.pacific_hfacc(:,:,:),1))'.*1000,[2000:20:2500])
    colormap(bluewhitered(length(2000:20:2500)-1))
    caxis([2000 2500])
    colorbar('YTick',2000:100:2500)
    set(gca,'Color','k');
    title('Pacific C_{sat} at in situ pCO_2 (tracer) [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')   
elseif isvar('cdiags.CSATSIT');
    subplot(ss,2,5)
    contourf(grd.latc,-grd.zc,squeeze(nanmean(cdiags.CSATSIT.*grd.atlantic_hfacc(:,:,:),1))'.*1000,[2000:20:2500])
    colormap(bluewhitered(length(2000:20:2500)-1))
    caxis([2000 2500])
    colorbar('YTick',2000:100:2500)
    set(gca,'Color','k');
    title('Atlantic C_{sat} at in situ pCO_2 (diagnostic) [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')
    
    subplot(ss,2,6)
    contourf(grd.latc,-grd.zc,squeeze(nanmean(cdiags.CSATSIT.*grd.pacific_hfacc(:,:,:),1))'.*1000,[2000:20:2500])
    colormap(bluewhitered(length(2000:20:2500)-1))
    caxis([2000 2500])
    colorbar('YTick',2000:100:2500)
    set(gca,'Color','k');
    title('Pacific C_{sat} at in situ pCO_2 (diagnostic) [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')
end

orient landscape
print -dpsc2 carbon_component_csat.ps

%% CRES at current atm pco2 and initial atm pco2
figure
subplot(ss,2,1)
contourf(grd.latc,-grd.zc,squeeze(nanmean(cdiags.CRES.*grd.atlantic_hfacc(:,:,:),1))'.*1000,[-300:20:300])
colormap(bluewhitered(length(-300:20:300)-1))
caxis([-300 300])
colorbar('YTick',-300:100:300)
set(gca,'Color','k');
title('Atlantic C_{res} at current pCO_2 [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')

subplot(ss,2,2)
contourf(grd.latc,-grd.zc,squeeze(nanmean(cdiags.CRES.*grd.pacific_hfacc(:,:,:),1))'.*1000,[-300:20:300])
colormap(bluewhitered(length(-300:20:300)-1))
caxis([-300 300])
colorbar('YTick',-300:100:300)
set(gca,'Color','k');
title('Pacific C_{res} at current pCO_2 [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')

subplot(ss,2,3)
contourf(grd.latc,-grd.zc,squeeze(nanmean(cdiags.CRESINIT.*grd.atlantic_hfacc(:,:,:),1))'.*1000,[-300:20:300])
colormap(bluewhitered(length(-300:20:300)-1))
caxis([-300 300])
colorbar('YTick',-300:100:300)
set(gca,'Color','k');
title('Atlantic C_{res} at initial pCO_2 [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')

subplot(ss,2,4)
contourf(grd.latc,-grd.zc,squeeze(nanmean(cdiags.CRESINIT.*grd.pacific_hfacc(:,:,:),1))'.*1000,[-300:20:300])
colormap(bluewhitered(length(-300:20:300)-1))
caxis([-300 300])
colorbar('YTick',-300:100:300)
set(gca,'Color','k');
title('Pacific C_{res} at initial pCO_2 [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')

if isvar('pdiags.cdis');
    subplot(ss,2,5)
    contourf(grd.latc,-grd.zc,squeeze(nanmean(pdiags.cdis.*grd.atlantic_hfacc(:,:,:),1))'.*1000,[-300:20:300])
    colormap(bluewhitered(length(-300:20:300)-1))
    caxis([-300 300])
    colorbar('YTick',-300:100:300)
    set(gca,'Color','k');
    title('Atlantic C_{dis} at in situ pCO_2 (tracer) [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')
    
    subplot(ss,2,6)
    contourf(grd.latc,-grd.zc,squeeze(nanmean(pdiags.cdis.*grd.pacific_hfacc(:,:,:),1))'.*1000,[-300:20:300])
    colormap(bluewhitered(length(-300:20:300)-1))
    caxis([-300 300])
    colorbar('YTick',-300:100:300)
    set(gca,'Color','k');
    title('Pacific C_{dis} at in situ pCO_2 (tracer) [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')

elseif isvar('cdiags.CSATSIT');
    subplot(ss,2,5)
    contourf(grd.latc,-grd.zc,squeeze(nanmean(cdiags.CRESSIT.*grd.atlantic_hfacc(:,:,:),1))'.*1000,[-300:20:300])
    colormap(bluewhitered(length(-300:20:300)-1))
    caxis([-300 300])
    colorbar('YTick',-300:100:300)
    set(gca,'Color','k');
    title('Atlantic C_{res} at in situ pCO_2 (diagnostic) [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')
    
    subplot(ss,2,6)
    contourf(grd.latc,-grd.zc,squeeze(nanmean(cdiags.CRESSIT.*grd.pacific_hfacc(:,:,:),1))'.*1000,[-300:20:300])
    colormap(bluewhitered(length(-300:20:300)-1))
    caxis([-300 300])
    colorbar('YTick',-300:100:300)
    set(gca,'Color','k');
    title('Pacific C_{res} at in situ pCO_2 (diagnostic) [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')
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

if isvar('trcant')
    subplot(ss,2,3)
    contourf(grd.latc,-grd.zc,squeeze(nanmean(cstar_res.*grd.atlantic_hfacc(:,:,:),1))'.*1000,[-300:20:300])
    colormap(bluewhitered(length(-300:20:300)-1))
    caxis([-300 300])
    colorbar('YTick',-300:100:300)
    set(gca,'Color','k');
    title('Atlantic C_{ant} (from tracer) [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')
    
    subplot(ss,2,4)
    contourf(grd.latc,-grd.zc,squeeze(nanmean(cstar_res.*grd.pacific_hfacc(:,:,:),1))'.*1000,[-300:20:300])
    colormap(bluewhitered(length(-300:20:300)-1))
    caxis([-300 300])
    colorbar('YTick',-300:100:300)
    set(gca,'Color','k');
    title('Pacific C_{ant} (from tracer) [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')
else
    subplot(ss,2,3)
    contourf(grd.latc,-grd.zc,squeeze(nanmean(cstar_res.*grd.atlantic_hfacc(:,:,:),1))'.*1000,[-300:20:300])
    colormap(bluewhitered(length(-300:20:300)-1))
    caxis([-300 300])
    colorbar('YTick',-300:100:300)
    set(gca,'Color','k');
    title('Atlantic C^{*}_{res} [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')
    
    subplot(ss,2,4)
    contourf(grd.latc,-grd.zc,squeeze(nanmean(cstar_res.*grd.pacific_hfacc(:,:,:),1))'.*1000,[-300:20:300])
    colormap(bluewhitered(length(-300:20:300)-1))
    caxis([-300 300])
    colorbar('YTick',-300:100:300)
    set(gca,'Color','k');
    title('Pacific C^{*}_{res} [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')
end

if isvar('cdiags.CSATSIT');
    subplot(ss,2,5)
    contourf(grd.latc,-grd.zc,squeeze(nanmean(cstar.*grd.atlantic_hfacc(:,:,:),1))'.*1000,[-300:20:300])
    colormap(bluewhitered(length(-300:20:300)-1))
    caxis([-300 300])
    colorbar('YTick',-300:100:300)
    set(gca,'Color','k');
    title('Atlantic C^{*}_{sat} [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')
    
    subplot(ss,2,6)
    contourf(grd.latc,-grd.zc,squeeze(nanmean(cstar.*grd.pacific_hfacc(:,:,:),1))'.*1000,[-300:20:300])
    colormap(bluewhitered(length(-300:20:300)-1))
    caxis([-300 300])
    colorbar('YTick',-300:100:300)
    set(gca,'Color','k');
    title('Pacific C^{*}_{sat} [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')
end

colormap(bluewhitered(length(-300:20:300)-1))
caxis([-300 300])
colorbar('YTick',-300:100:300)
orient landscape
print -dpsc2 carbon_component_cant.ps

%% CSOFT and CCARB
figure
subplot(2,2,1)
contourf(grd.latc,-grd.zc,squeeze(nanmean(cdiags.CSOFT.*grd.atlantic_hfacc(:,:,:),1))'.*1000,[-300:20:300])
colormap(bluewhitered(length(-300:20:300)-1))
caxis([-300 300])
colorbar('YTick',-300:100:300)
set(gca,'Color','k');
title('Atlantic C_{soft} [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')

subplot(2,2,2)
contourf(grd.latc,-grd.zc,squeeze(nanmean(cdiags.CSOFT.*grd.pacific_hfacc(:,:,:),1))'.*1000,[-300:20:300])
colormap(bluewhitered(length(-300:20:300)-1))
caxis([-300 300])
colorbar('YTick',-300:100:300)
set(gca,'Color','k');
title('Pacific C_{soft} [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')

subplot(2,2,3)
contourf(grd.latc,-grd.zc,squeeze(nanmean(cdiags.CCARB.*grd.atlantic_hfacc(:,:,:),1))'.*1000,[-300:20:300])
colormap(bluewhitered(length(-300:20:300)-1))
caxis([-300 300])
colorbar('YTick',-300:100:300)
set(gca,'Color','k');
title('Atlantic C_{carb} [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')

subplot(2,2,4)
contourf(grd.latc,-grd.zc,squeeze(nanmean(cdiags.CCARB.*grd.pacific_hfacc(:,:,:),1))'.*1000,[-300:20:300])
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
contourf(grd.latc,-grd.zc,squeeze(nanmean(pdiags.ppre.*grd.atlantic_hfacc(:,:,:),1))'.*1000,[0:0.1:2])
colormap(parula(length(0:0.1:2)-1))
caxis([0 2])
colorbar('YTick',0:0.5:2)
set(gca,'Color','k');
title('Atlantic P_{pre} [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')

subplot(2,2,2)
contourf(grd.latc,-grd.zc,squeeze(nanmean(pdiags.ppre.*grd.pacific_hfacc(:,:,:),1))'.*1000,[0:0.1:2])
colormap(parula(length(0:0.1:2)-1))
caxis([0 2])
colorbar('YTick',0:0.5:2)
set(gca,'Color','k');
title('Pacific P_{pre} [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')

subplot(2,2,3)
contourf(grd.latc,-grd.zc,squeeze(nanmean(cdiags.PREG.*grd.atlantic_hfacc(:,:,:),1))'.*1000,[0:0.1:2])
colormap(parula(length(0:0.1:2)-1))
caxis([0 2])
colorbar('YTick',0:0.5:2)
set(gca,'Color','k');
title('Atlantic P_{reg} [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')

subplot(2,2,4)
contourf(grd.latc,-grd.zc,squeeze(nanmean(cdiags.PREG.*grd.pacific_hfacc(:,:,:),1))'.*1000,[0:0.1:2])
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
contourf(grd.latc,-grd.zc,squeeze(nanmean(pdiags.cpre.*grd.atlantic_hfacc(:,:,:),1))'.*1000,[2000:20:2500])
colormap(bluewhitered(length(2000:20:2500)-1))
caxis([2000 2500])
colorbar('YTick',2000:100:2500)
set(gca,'Color','k');
title('Atlantic C_{pre} [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')

subplot(2,2,2)
contourf(grd.latc,-grd.zc,squeeze(nanmean(pdiags.cpre.*grd.pacific_hfacc(:,:,:),1))'.*1000,[2000:20:2500])
colormap(bluewhitered(length(2000:20:2500)-1))
caxis([2000 2500])
colorbar('YTick',2000:100:2500)
set(gca,'Color','k');
title('Pacific C_{pre} [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')

subplot(2,2,3)
contourf(grd.latc,-grd.zc,squeeze(nanmean((pdiags.dic-pdiags.cpre).*grd.atlantic_hfacc(:,:,:),1))'.*1000,[-300:20:300])
colormap(bluewhitered(length(-300:20:300)-1))
caxis([-300 300])
colorbar('YTick',-300:100:300)
set(gca,'Color','k');
title('Atlantic C_{reg} [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')
colormap(bluewhitered(length(-300:20:300)-1))
caxis([-300 300])
colorbar('YTick',-300:100:300)

subplot(2,2,4)
contourf(grd.latc,-grd.zc,squeeze(nanmean((pdiags.dic-pdiags.cpre).*grd.pacific_hfacc(:,:,:),1))'.*1000,[-300:20:300])
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
contourf(grd.latc,-grd.zc,squeeze(nanmean(pdiags.apre.*grd.atlantic_hfacc(:,:,:),1))',[2.0:0.01:2.5])
colormap(parula(length(2.0:0.01:2.5)-1))
caxis([2.0 2.4])
colorbar('YTick',2.0:0.1:2.4)
set(gca,'Color','k');
title('Atlantic A_{pre} [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')

subplot(2,2,2)
contourf(grd.latc,-grd.zc,squeeze(nanmean(pdiags.apre.*grd.pacific_hfacc(:,:,:),1))',[2.0:0.01:2.5])
colormap(parula(length(2.0:0.01:2.5)-1))
caxis([2.0 2.4])
colorbar('YTick',2.0:0.1:2.4)
set(gca,'Color','k');
title('Pacific A_{pre} [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')

subplot(2,2,3)
contourf(grd.latc,-grd.zc,squeeze(nanmean(cdiags.AREG.*grd.atlantic_hfacc(:,:,:),1))'.*1000,[-10:10:500])
colormap(parula(length(-10:10:500)-1))
caxis([0 500])
colorbar('YTick',0:100:500)
set(gca,'Color','k');
title('Atlantic A_{reg} [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')

subplot(2,2,4)
contourf(grd.latc,-grd.zc,squeeze(nanmean(cdiags.AREG.*grd.pacific_hfacc(:,:,:),1))'.*1000,[-10:10:500])
colormap(parula(length(-10:10:500)-1))
caxis([0 500])
colorbar('YTick',0:100:500)
set(gca,'Color','k');
title('Pacific A_{reg} [mmol m^{-3}]'); xlabel('Latitude'); ylabel('Depth [m]')

orient landscape
print -dpsc2 carbon_component_alk.ps

