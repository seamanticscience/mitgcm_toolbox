function [bioac]=mit_plotbio(grd,tavesteps)
% sections of carbon and nutrients

% $Header: /u/gcmpack/MITgcm/verification/tutorial_global_oce_latlon/diags_matlab/mit_plotmeandrift.m,v 1.3 2006/08/12 20:25:13 jmc Exp $
% $Name:  $
if nargin==0 || ~exist('grd','var')
    grd=mit_loadgrid;
    grd=mit_oceanmasks(grd);
end

if nargin==1 || ~exist('tavesteps','var')
    tavesteps=mit_timesteps('tave');
end

dicdiag=rdmnc(strrep(tavesteps.filearr(2:end-1),'tave','dicDiag'),tavesteps.timesteps(tavesteps.kmax));

zgrid=repmat(permute(grd.dz,[3,2,1]),[128,64,1]);
 
dicbioint=117.*nansum(dicdiag.DICBIOA(:,:,:).*zgrid,3).*(60*60*24*360); % mol C/m3/yr

contourf(grd.lonc,grd.latc,dicbioint'.*grd.hfacc(:,:,1)',[0:2:26])
caxis([0 26]);colormap(parula(13));colorbar('YTick',0:2:26);
set(gca,'Xlim',[min(grd.long) max(grd.lonc)],'XTickLabelMode','manual','XTickMode','manual','XTick',[0:120:300,max(grd.long)],'XTickLabel',[0:120:360],'YLim',[min(grd.latc) max(grd.latc)],'YTick',[-80:40:80],'FontSize',16)
title('Biological Community Production [mol C.m-2.yr]','FontSize',16)
orient landscape
print -dpsc gchem_integrated_bioac.ps


bioac=nansum(dicbioint(:).*grd.rac(:)).*12e-15