function mit_plotfekscav(grd,tavesteps)

if nargin==0;
    grd=mit_loadgrid;
elseif nargin==1
    tavesteps=mit_timesteps('tave');
end

if ~exist('dic_tave','var')
    dic_tave=rdmnc(strrep(tavesteps.filearr(2:end-1),'tave','dic_tave'),tavesteps.timesteps(end));
end

% Plot particle concentration = export_flux/sinking_speed
pop_conc=dic_tave.dic_epflux_ave.*((60*60*24*360)/2900); % in mol/m3

figure(1)
subplot(141)
plot(squeeze(nanmean(nanmean(pop_conc.*grd.atlantic_hfacc*1e6,1),2)),-grd.zc)
xlabel('POP Conc [umol/m3]','FontSize',12)
ylabel('Depth [m]','FontSize',12)
set(gca,'FontSize',12,'YLim',[-5000 0])
subplot(1,4,2:4)
contourf(grd.latc,-grd.zc,squeeze(nanmean(pop_conc.*grd.atlantic_hfacc*1e6,1))')
colorbar
set(gca,'FontSize',12,'YLim',[-5000 0])
xlabel('Latitude','FontSize',12)
title('Atlantic Particulate Organic Phosphate (POP) Export Concentration [umol/m3]','FontSize',12)
orient landscape
print -dpsc2 export_flux_pop.ps
fixpslinestyle('export_flux_pop.ps')

if exist('data.dic','file')
    scav_surf_min=3.0542695473251e-08 %mit_getparm('data.dic','scav_surf_min')
    scav_inter=mit_getparm('data.dic','scav_inter')
    scav_rat=mit_getparm('data.dic','Kscav')
    scav_exp=mit_getparm('data.dic','scav_exp')
else
    % Calculate variable scavenging rate
    scav_surf_min=0.19/(86400*360) % surface scavenging rate (otherwise Kscav goes to zero at surface)
    scav_inter=0.079
    scav_rat=0.0035/86400
    % scav_rat = 0.005/d = 1.8/yr max Kscav. Mean Kscav = 0.51
    % scav_rat = 0.0035/d=                   Mean Kscav = 0.323 with mean of 0.19 below 1000m
    % scav_rat = 0.0025/d= 0.9/yr max Kscav. Mean Kscav = 0.26
    % scav_rat = 0.002/d= 0.72/yr max Kscav. Mean Kscav = 0.20
    % scav_rat = 0.001/d= 0.36/yr max Kscav. Mean Kscav = 0.10
    scav_exp=0.58
end

poc_conc=(pop_conc.*1000)/1.1321e-4; % Steph suggested 1.1321e-4, but Darwin [pop] is mmol m-3
Kscav=scav_inter*scav_rat*(poc_conc.^scav_exp);
Kscav(:,:,1)=scav_surf_min;

figure
subplot(141)
plot(squeeze(nanmean(nanmean(Kscav.*grd.hfacc.*(60*60*24*360),1),2)),-grd.zc)
set(gca,'FontSize',12,'YLim',[-5000 0])
xlabel('Fe Scavenging Rate [yr-1]','FontSize',12)
ylabel('Depth [m]','FontSize',12)
subplot(1,4,2:4)
contourf(grd.latc,-grd.zc,squeeze(nanmean(Kscav.*grd.hfacc*(60*60*24*360),1))')
colorbar
set(gca,'FontSize',12,'YLim',[-5000 0])
xlabel('Latitude','FontSize',12)
title('Spatially Variable Fe Scavenging Rate [yr-1]','FontSize',12)
orient landscape
print -dpsc2 variable_kscav_rate.ps
fixpslinestyle('variable_kscav_rate.ps')