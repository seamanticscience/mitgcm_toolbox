function []=mit_plotsurfacefields(grd,tavesteps,depthlev)
% Plot sea surface figures for t, s, density, fluxes etc

if nargin==0;
    grd=mit_loadgrid;
elseif nargin==1
    tavesteps=mit_timesteps('tave');
end

% Do Atlantic centric rotation
m_proj('mollweide','clong',0);
[nlonc nic] = globe(grd.lonc,0);

% Load data
tave=rdmnc(tavesteps.filearr(2:end-1),'Ttave','Stave','uVeltave','vVeltave',tavesteps.timesteps(end));
surf=rdmnc(strrep(tavesteps.filearr(2:end-1),'tave','surfDiag'),'surForcT','surForcS','TFLUX',tavesteps.timesteps(end));

% Make pressure grid from depths
[x,y]=size(grd.yc);
mit_press=NaN(grd.nx,grd.ny,grd.nz);
for i=1:x
    for j=1:y
        mit_press(i,j,:)=sw_pres(grd.zc,grd.yc(i,j)');
    end
end
%mit_press=permute(mit_press,[2 1 3]);

% Calculate full depth sigma0
mit_sigma0=sw_pden(tave.Stave(nic,:,1),tave.Ttave(nic,:,1),mit_press(nic,:,1),mit_press(nic,:,1).*0).*grd.hfacc(nic,:,1) -1000; 
mit_sigma0=mit_sigma0';

%% Plot comparisons of heat and freshwater fluxes
figure
contourf(nlonc-360,grd.latc,surf.TFLUX(nic,:,1)'.*grd.hfacc(nic,:,1)',-500:20:500);set(gca,'FontSize',14)
canom;cmapa(4);caxis([-150 150]);colorbar('FontSize',14); % format_ticks(gca,'^{\circ}','^{\circ}');
title('Annual mean MITgcm Total Heat Forcing (W/m2, positive for heat gain)','FontSize',14,'FontWeight','bold'); ylabel('Latitude','FontSize',14); xlabel('Longitude','FontSize',14)
orient landscape
print -dpsc surface_heat_flux.ps

scal=(360.*24.*60.*60)/1000;
figure
contourf(nlonc-360,grd.latc,surf.surForcS(nic,:,1)'.*grd.hfacc(nic,:,1)'.*scal,[-150:10:150]);set(gca,'FontSize',14)
canom;cmapa(4);caxis([-80 80]);colorbar('FontSize',14); % format_ticks(gca,'^{\circ}','^{\circ}');
title('Annual mean MITgcm Total Salinity Forcing (Kg/m2/yr, positive increases S)','FontSize',14,'FontWeight','bold'); ylabel('Latitude','FontSize',14); xlabel('Longitude','FontSize',14)
orient landscape
print -dpsc surface_salinity_flux.ps

figure
contourf(nlonc-360,grd.latc,nanmean(tave.Ttave(nic,:,grd.zc<depthlev),3)'.*grd.hfacc(nic,:,1)',-5:2.5:40);set(gca,'FontSize',14)
caxis([-5 35]);colorbar('FontSize',14); % format_ticks(gca,'^{\circ}','^{\circ}');
title('MITgcm SST (degC)','FontSize',14,'FontWeight','bold'); ylabel('Latitude','FontSize',14); xlabel('Longitude','FontSize',14)
orient landscape
print -dpsc surface_temperature.ps

figure
contourf(nlonc-360,grd.latc,nanmean(tave.Stave(nic,:,grd.zc<depthlev),3)'.*grd.hfacc(nic,:,1)',[20:0.25:40]);set(gca,'FontSize',14)
caxis([30 37]);colorbar('FontSize',14); % format_ticks(gca,'^{\circ}','^{\circ}');
title('MITgcm SSS','FontSize',14,'FontWeight','bold'); ylabel('Latitude','FontSize',14); xlabel('Longitude','FontSize',14)
orient landscape
print -dpsc surface_salinity.ps

speed=sqrt(tave.uVeltave.*tave.uVeltave+tave.vVeltave.*tave.vVeltave);
figure
contourf(nlonc-360,grd.latc,nanmean(log10(speed(nic,:,grd.zc<depthlev)),3)'.*grd.hfacc(nic,:,1)',-5:0.25:2)
caxis([-3 0]);colorbar('FontSize',14);
hold on
quiver(nlonc(1:2:end)-360,grd.latc(1:2:end),nanmean(tave.uVeltave(nic(1:2:end),1:2:end,grd.zc<depthlev),3)'.*grd.hfacc(nic(1:2:end),1:2:end,1)',...
     nanmean(tave.vVeltave(nic(1:2:end),1:2:end,grd.zc<depthlev),3)'.*grd.hfacc(nic(1:2:end),1:2:end,1)',10,'k')
title(['MITgcm Velocity in the upper ',num2str(depthlev),'m'],'FontSize',14,'FontWeight','bold'); ylabel('Latitude','FontSize',14); xlabel('Longitude','FontSize',14)
orient landscape
print -dpsc surface_velocity.ps



