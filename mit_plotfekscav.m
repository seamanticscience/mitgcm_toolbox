% function mit_plotfekscav(grd,tavesteps)
% 
% if nargin==0;
     grd=mit_loadgrid;
     tavesteps=mit_timesteps('tave');
% elseif nargin==1
%     tavesteps=mit_timesteps('tave');
% end
% 
 if ~exist('grd.atlantic_hfacc','var')
     grd=mit_oceanmasks(grd);
 end

if ~exist('dic_tave','var')
    dic_tave=rdmnc(strrep(tavesteps.filearr(2:end-1),'tave','dic_tave'),tavesteps.timesteps(end));
end

% Plot particle concentration = export_flux/sinking_speed
pop_conc=dic_tave.dic_epflux_ave.*((86400*360)/2900); % in mol/m3

figure
subplot(141)
plot(squeeze(nanmean(nanmean(pop_conc.*grd.hfacc*1e6,1),2)),-grd.zc)
xlabel('POP Conc [umol/m3]','FontSize',12)
ylabel('Depth [m]','FontSize',12)
set(gca,'FontSize',12,'YLim',[-5000 0])
subplot(1,4,2:4)
contourf(grd.latc,-grd.zc,squeeze(nanmean(pop_conc.*grd.hfacc*1e6,1))')
colorbar
set(gca,'FontSize',12,'YLim',[-5000 0])
xlabel('Latitude','FontSize',12)
title('Particulate Organic Phosphate (POP) Export Concentration [umol/m3]','FontSize',12)
orient landscape
print -dpsc2 export_flux_pop.ps
fixpslinestyle('export_flux_pop.ps')

% if exist('data.dic','file')
%     scav_surf_min=3.0542695473251e-08 %mit_getparm('data.dic','scav_surf_min')
%     scav_init=mit_getparm('data.dic','scav_init')
%     scav_rat=mit_getparm('data.dic','Kscav')
%     scav_exp=mit_getparm('data.dic','scav_exp')
% else
    % Calculate variable scavenging rate
    scav_surf_min=0.19/(86400*360) % surface scavenging rate (otherwise Kscav goes to zero at surface)
    R_pop2poc = 117*12
    scav_init=0.079./86400 % units of L^(units_scav_exp) mg^-(units_scav_exp) s-1 
    %scav_init=0.125./86400
    scav_ratio=0.2
    
    % scav_rat=0.0035/86400
    % scav_rat = 0.005/d = 1.8/yr max Kscav. Mean Kscav = 0.51
    % scav_rat = 0.0035/d=                   Mean Kscav = 0.323 with mean of 0.19 below 1000m
    % scav_rat = 0.0025/d= 0.9/yr max Kscav. Mean Kscav = 0.26
    % scav_rat = 0.002/d= 0.72/yr max Kscav. Mean Kscav = 0.20
    % scav_rat = 0.001/d= 0.36/yr max Kscav. Mean Kscav = 0.10
    scav_exp=0.58
%end

scav_poc=pop_conc.*R_pop2poc; % converts from the regular model units to weird units

Kscav=(scav_ratio*scav_init*(scav_poc.^scav_exp));   

% Surface scavenging by dust
dust=mit_readfield('mah_somtimind_molfem2s.bin',[128,64,12],'real*4');
wsp_dust=2; % m/s sinking speed
gdustm3=nanmean((dust.*58.845)./(wsp_dust.*0.035),3);
scav_dust=repmat(gdustm3.*((150*1000)./(86400*1000)),[1,1,15]); % m3/g/s

%Kscav(:,:,1)=scav_surf_min;
Kscav(:,:,1)=scav_surf_min;

%Kscav=Kscav+(5e-5./86400)+scav_dust;

figure
subplot(141)
plot(squeeze(nanmean(nanmean(Kscav.*grd.hfacc.*(60*60*24*360),1),2)),-grd.zc)
set(gca,'FontSize',12,'YLim',[-5000 0])
xlabel('Fe Scavenging Rate [yr-1]','FontSize',12)
ylabel('Depth [m]','FontSize',12)
subplot(1,4,2:4)
contourf(grd.latc,-grd.zc,squeeze(nanmean(Kscav.*grd.hfacc*(60*60*24*360),1))',[0:0.05:1])
caxis([0 0.5]);colormap(parula(10));colorbar
set(gca,'FontSize',12,'YLim',[-5000 0])
xlabel('Latitude','FontSize',12)
title('Spatially Variable Fe Scavenging Rate [yr-1]','FontSize',12)
orient landscape
print -dpsc2 variable_kscav_rate.ps
fixpslinestyle('variable_kscav_rate.ps')