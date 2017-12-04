function [ptr_tave,cntrl_fe]=mit_plotatlanticfe(varargin)

% plots atlantic iron section and anomaly at the last timestep


if nargin~=0
    nyears=varargin{1};
end

ptr_tave=rdmnc('ptr_tave.*.glob.nc','fe');

if ~exist('nyears','var')
    nyears=length(ptr_tave.iters_from_file);
end

grd=mit_loadgrid;
grd=mit_oceanmasks(grd);

atlantic_fe=ptr_tave.fe(:,:,:,nyears).*grd.atlantic_hfacc;

if exist('geotraces1x1_init_fe.bin','file')
    cntrl_fe=mit_readfield('geotraces1x1_init_fe.bin',[360,160,grd.nz],'real*4');
    atlantic_cntrl_fe=cntrl_fe.*grd.atlantic_hfacc;
elseif exist('/Volumes/PhD_Data/controlrun/controlrun_anomplotref/ptr_tave.20k.glob.nc','file')
    cntrl_ptr_tave=rdmnc('/Volumes/PhD_Data/controlrun/controlrun_anomplotref/ptr_tave.20k.glob.nc','fe');
    atlantic_cntrl_fe=cntrl_ptr_tave.fe.*grd.atlantic_hfacc;
else
    atlantic_cntrl_fe=[];
end

figure
contourf(grd.latc,-grd.zc,squeeze(nanmean(atlantic_fe,1))',[0:1e-7:12e-7])
colormap(parula(12));caxis([0 12e-7]); colorbar
set(gca,'FontSize',12)
xlabel('Latitude','Fontsize',12)
ylabel('Depth','Fontsize',12)
title('Zonal average Atlantic section of Fe (mol/m3)','FontSize',14)
orient landscape
print -dpsc2 atlantic_fe_conc.ps
fixpslinestyle('atlantic_fe_conc.ps')

if ~isempty(atlantic_cntrl_fe)
    figure
    contourf(grd.latc,-grd.zc,squeeze(nanmean(atlantic_fe,1)-nanmean(atlantic_cntrl_fe,1))',[-10e-7:1e-7:10e-7])
    caxis([-10e-7 10e-7]); colormap(bluewhitered(20)); colorbar
    set(gca,'FontSize',12)
    xlabel('Latitude','Fontsize',12)
    ylabel('Depth','Fontsize',12)
    title('Zonal average Atlantic section of Fe Anomaly (mol/m3)','FontSize',14)
    orient landscape
    print -dpsc2 anom_atlantic_fe_conc.ps
    fixpslinestyle('anom_atlantic_fe_conc.ps')
end