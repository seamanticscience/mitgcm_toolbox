function [u]=mit_plotzonalvelocity(grd,tavesteps)
% m-file: mit_plotzonalvelocity.m

% $Header: /u/gcmpack/MITgcm/verification/tutorial_global_oce_latlon/diags_matlab/mit_plotzonalvelocity.m,v 1.3 2006/08/12 20:25:13 jmc Exp $
% $Name:  $
if nargin==0;
    grd=mit_loadgrid;
    grd = mit_oceanmasks(grd);
elseif nargin==1
    tavesteps=mit_timesteps('tave');
end

if ~isfield(grd,'pacific_hfacw')
    grd = mit_oceanmasks(grd);
end

% select timestep
k=tavesteps.kmax;

tave=rdmnc((tavesteps.filearr(2:end-1)),'uVeltave',tavesteps.timesteps(k));

u=tave.uVeltave;
up = tave.uVeltave.*grd.pacific_hfacw;
upzm = mit_zonalmean(up,grd.pacific_hfacw,grd.dxg);
ua = tave.uVeltave.*grd.atlantic_hfacw;
uazm = mit_zonalmean(ua,grd.atlantic_hfacw,grd.dxg);

% caxup = [min(upzm(:)) max(upzm(:))];
 uplev = -1:.01:1;
% if max(abs(caxup)) < .1
%   uplev = .5*uplev;
% end
% if max(abs(caxup)) < .2
%   uplev = .2*uplev;
% end
% 
% caxua = [min(uazm(:)) max(uazm(:))];
 ualev = -1:.025:1;
% if max(abs(caxua)) < .1
%   ualev = .5*ualev;
% end
% if max(abs(caxua)) < .2
%   ualev = .2*ualev;
% end

ixpw = find(grd.long>163,1,'first');
ixpc = find(grd.long>199,1,'first');
ixpe = find(grd.long>255,1,'first');

ixaw = find(grd.long>303,1,'first');
ixac = find(grd.long>331,1,'first');
ixae = find(grd.long>357,1,'first');

zaxis = -grd.zc;
zaxis = -grd.zg;
yaxis = grd.latg;
%zaxis = -[1:grd.nz];
figure('PaperPosition',[0.25 0.621429 8 9.75714])
suptitle(['experiment ' grd.dname ...
	  ', timestep = ' num2str(tavesteps.timesteps(k)) ...
	  ', zonal velocity [m/s]'])
sh(1) = subplot(4,2,1);
contourf(yaxis,zaxis,upzm',uplev);colorbar
title(sh(1),'Pacific Ocean: zonal average') 
sh(3) = subplot(4,2,3);
contourf(yaxis,zaxis,sq(up(ixpw,:,:))',uplev);colorbar
title(sh(3),['section at ',num2str(grd.long(ixpw),3),'^{\circ}E'])
sh(5) = subplot(4,2,5);
contourf(yaxis,zaxis,sq(up(ixpc,:,:))',uplev);colorbar
title(sh(5),['section at ',num2str(grd.long(ixpc),3),'^{\circ}E'])
sh(7) = subplot(4,2,7);
contourf(yaxis,zaxis,sq(up(ixpe,:,:))',uplev);colorbar
title(sh(7),['section at ',num2str(grd.long(ixpe),3),'^{\circ}E'])
sh(2) = subplot(4,2,2);
contourf(yaxis,zaxis,uazm',ualev);colorbar
title(sh(2),'Atlantic Ocean: zonal average') 
sh(4) = subplot(4,2,4);
contourf(yaxis,zaxis,sq(ua(ixaw,:,:))',ualev);colorbar
title(sh(4),['section at ',num2str(grd.long(ixaw),3),'^{\circ}E'])
sh(6) = subplot(4,2,6);
contourf(yaxis,zaxis,sq(ua(ixac,:,:))',ualev);colorbar
title(sh(6),['section at ',num2str(grd.long(ixac),3),'^{\circ}E'])
sh(8) = subplot(4,2,8);
contourf(yaxis,zaxis,sq(ua(ixae,:,:))',ualev);colorbar
title(sh(8),['section at ',num2str(grd.long(ixae),3),'^{\circ}E'])
set(gca,'CLim',[-0.5 0.5])

set(sh,'CLim',[-0.5 0.5])

set(sh,'xlim',[-1 1]*80,'ylim',[-5200 0])
set(sh,'layer','top')
%set(sh(1:2:end),'clim',[uplev(1) uplev(end)])
%set(sh(2:2:end),'clim',[ualev(1) ualev(end)])
%set(gcf,'currentAxes',sh(end-1));colorbar
%set(gcf,'currentAxes',sh(end));colorbar
canom;cmapa(4)
orient portrait
print -dpsc zonal_velocity_sections.ps
