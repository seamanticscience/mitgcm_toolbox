function []=mit_plotstreamfunctions(grd,tavesteps,filename)

% $Header: /u/gcmpack/MITgcm/verification/tutorial_global_oce_latlon/diags_matlab/mit_plotstreamfunctions.m,v 1.3 2006/08/12 20:25:13 jmc Exp $
% $Name:  $
%Set up matricies
if nargin==0 || ~exist('grd','var')
    grd=mit_loadgrid;
end

if nargin==1 || ~exist('tavesteps','var')
    tavesteps=mit_timesteps('tave');
end

if nargin==2 || ~exist('filename','var')
    filename='psi_data.nc';
end

if ~isfield(grd,'pacific_hfacw')
    grd = mit_oceanmasks(grd);
end

%% Are we looking at a year or longer?
if (mit_getparm('data','nTimeSteps')/(31104000/grd.deltattracer)) == 1
%    disp('Loading one year')
    ktsteps=tavesteps.timesteps;
    ktimsteps=tavesteps.tim;
else
    ktsteps=tavesteps.timesteps(end);
    ktimsteps=tavesteps.tim(end);
%    disp('Loading last timestep')
end

%% Load the last timestep of Circulation files for plotting.
    % Load residual velocity diagnostics or calculate from output
    if exist(strrep(tavesteps.filearr(2:end-1),'tave','gmDiag'),'file') ...
            && exist(strrep(tavesteps.filearr(2:end-1),'tave','oceDiag'),'file')
        
        diagsteps=mit_timesteps(mit_getparm('data.diagnostics','filename'));
        
        if  nc_isvar(strrep(diagsteps.filearr(2:end-1),'surf','gm'),'GM_V_RES')
            
            disp(['Loading eddy velocities from ',strrep(diagsteps.filearr(2:end-1),'surf','gm')])
            
            gmdiag=rdmnc(strrep(diagsteps.filearr(2:end-1),'surf','gm'),'GM_V_RES',ktsteps);          
            vresk=nanmean(gmdiag.GM_V_RES,4);
            
            if nc_isvar(strrep(diagsteps.filearr(2:end-1),'surf','oce'),'THETA') ...
                    && nc_isvar(strrep(diagsteps.filearr(2:end-1),'surf','oce'),'SALT')
                tave=rdmnc(strrep(diagsteps.filearr(2:end-1),'surf','oce'),'THETA','SALT',ktsteps);
                theta=nanmean(tave.THETA,4);
                salt=nanmean(tave.SALT,4);
            else
                tave = rdmnc(tavesteps.filearr(2:end-1),'Ttave','Stave',ktsteps);
                theta=nanmean(tave.Ttave,4);
                salt=nanmean(tave.Stave,4);
            end
        end
    end
    
    % Load the regular data if those arrays have not been filled with the
    % better diagnostics...
    if ~exist('vresk','var')
        disp(['Loading eddy diffusivities from ',strrep(tavesteps.filearr(2:end-1),'tave','gm_tave')])
        
        tave = rdmnc(tavesteps.filearr(2:end-1),'uVeltave','vVeltave','Ttave','Stave',ktsteps);
        if ~strcmp(grd.buoyancy,'OCEANIC');
            vvelk = nanmean(tave.vVeltave(:,:,end:-1:1),4);
        else
            vvelk = nanmean(tave.vVeltave,4);
        end
        
        theta=nanmean(tave.Ttave,4);
        salt=nanmean(tave.Stave,4);
                
        % Calculate GM bolus velocities
        gm = rdmnc(strrep(tavesteps.filearr(2:end-1),'tave','gm_tave'),'Kwx','Kwy','Kwz',ktsteps);
        
        Kwx=nanmean(gm.Kwx,4);
        Kwy=nanmean(gm.Kwy,4);
        
        [~,veddk,~]=mit_gmvel(Kwx,Kwy,grd.nx,grd.ny,grd.nz,1,grd.dxg,grd.dyg,grd.dz,grd.cmask,grd.umask,grd.vmask,grd.rac);
        vresk=vvelk+veddk;
    end
    
    mitpsi=rdmnc(filename,'global_psi','atlantic_psi','pacific_psi','gbaro_psi',ktimsteps);
    mitpsi.global_psi=nanmean(mitpsi.global_psi,3);
    mitpsi.atlantic_psi=nanmean(mitpsi.atlantic_psi,3);
    mitpsi.pacific_psi=nanmean(mitpsi.pacific_psi,3);
    mitpsi.gbaro_psi=nanmean(mitpsi.gbaro_psi,3);
    
    
%% Calculate overturning on density surfaces at last timestep
%mit_press=NaN(grd.nx,grd.ny,grd.nz);
%for i=1:grd.nx
%    for j=1:grd.ny
%        mit_press(i,j,:)=sw_pres(grd.zc,grd.yc(i,j)');
%    end
%end

% Dont forget Ttave is potential temperature already!
P=ini_poly3;
mit_sigma25=dens_poly3(P(find(grd.zc<=2100,1,'last')),theta,salt).*grd.hfacc;
ri=floor(nanmin((mit_sigma25(:)))):0.05:ceil(nanmax((mit_sigma25(:))));

%mit_sigma25=sw_pden(tave.Stave,sw_ptmp(tave.Stave,tave.Ttave,zeros(size(mit_press)),mit_press),mit_press,ones(size(mit_press)).*2500).*grd.hfacc -1000;
%ri=floor(nanmin((mit_sigma25(:)))):0.05:ceil(nanmax((mit_sigma25(:))));

%hfacc_on_density=nan(grd.nx,grd.ny,length(ri));
% for i=1:grd.nx
%     for j=1:grd.ny
%         start=find(ri*100>floor(nanmin(mit_sigma25(i,j,:)*100)),1,'first');
%         stop=find(ri*100>floor(nanmax(mit_sigma25(i,j,:)*100)),1,'first');
%         hfacc_on_density(i,j,start:stop)=1;
%     end
% end
%
dx_on_density=repmat(grd.dxg,[1,1,length(ri)]);
% Pretty sure using griddata is a bad idea, really overestimates low
% latitude isopycnal thickness
% dz_on_density=griddata(xi,yi,change(mit_sigma25,'==',NaN,0),repmat(permute(grd.dz,[3,2,1]),[grd.nx,grd.ny,1]),xo,yo,zo);
dz_on_density=zeros(grd.nx,grd.ny,length(ri));
hfacc_on_density=zeros(grd.nx,grd.ny,length(ri));

for i=1:grd.nx
    for j=1:grd.ny
        if ~isnan(grd.hfacc(i,j,1))
            nzl=sum(isfinite(mit_sigma25(i,j,:)));
            if nzl==1
                dz_on_density(i,j,find(ri>=mit_sigma25(i,j,1),1,'first'))=grd.dz(1);
            elseif issorted(squeeze(mit_sigma25(i,j,1:nzl)))
                zc_on_ri_levels=interp1(squeeze(mit_sigma25(i,j,1:nzl)),grd.zc(1:nzl),ri,'linear','extrap')';
                % Now mask out values that are above the surface (-ve)or below
                % deepest level
                zc_on_ri_levels(zc_on_ri_levels<0)=0;
                zc_on_ri_levels(find(zc_on_ri_levels>=grd.depth(i,j),1,'first'))=grd.depth(i,j);
                zc_on_ri_levels(zc_on_ri_levels>grd.depth(i,j))=0;
                dz_on_density(i,j,:)=[0;diff(zc_on_ri_levels)];
            else
                dz_on_density(i,j,find(ri>=mean(squeeze(mit_sigma25(i,j,1:nzl))),1,'first'))=mean(grd.dz(1:nzl));
            end
        end
    end
end
dz_on_density(dz_on_density<0)=0;
hfacc_on_density=change(+(dz_on_density>0),'==',0,NaN);

%         [xin,yin,zin,tin]=ndgrid(ecco_latg,ecco_long,ecco_zc,1:ntsteps);
%         [xg,yg,zg,tout]=ndgrid(grid.latg,grid.long,grid.zc,1:ntsteps);
%         FS=griddedInterpolant(xin,yin,zin,tin,permute(inpaint_nans(tmpkz,4),[2,1,3,4]),'linear','linear');
%         kppkz=permute(FS(xg,yg,zg,tout),[2,1,3,4]); 
[xi,yi,zi]=meshgrid(grd.lonc,grd.latc,grd.zc);
[xo,yo,zo]=meshgrid(grd.lonc,grd.latc,ri);

vres_on_density=griddata(permute(xi,[2,1,3]),permute(yi,[2,1,3]),change(mit_sigma25,'==',NaN,0),vresk,permute(xo,[2,1,3]),permute(yo,[2,1,3]),permute(zo,[2,1,3]),'linear').*hfacc_on_density;

sigma_psi=mit_sigma_overturning(vres_on_density,hfacc_on_density,dx_on_density,dz_on_density,0);

%% plot stream functions
%
figure 
otlev = [-50:5:50];
global_psi_max=mitpsi.global_psi;
contourf(grd.latg,-grd.zgpsi/1000,global_psi_max'*1e-6,otlev); 
caxis([-1 1].*max(otlev)); colormap(bluewhitered(length(otlev)-1)); colorbar;
hold on
[cs h1] = contour(grd.latg,-grd.zgpsi/1000,global_psi_max'*1e-6,[0 0]); 
%clh1 = clabel(cs);
set(h1,'LineWidth',2,'EdgeColor','k');
set(gca,'XLim',[min(grd.latc) max(grd.latc)],'XTick',[-80:20:80],'FontSize',16)
xlabel('Latitude','FontSize',16);ylabel('Depth [km]','FontSize',16)
hold off
% psimin = min(min(global_psi_max(:,5:end)));
% [iy iz] = find(abs(global_psi_max(:)-psimin)<=1e-4);
% text(grd.latg(iy),-grd.zgpsi(iz), ...
%      ['\leftarrow ' num2str(psimin*1e-6,'%5.1f')], ...
%      'horizontalalignment','left')
title('Global overturning streamfunction [Sv]','FontSize',16)
orient landscape
print -dpsc global_overturning.ps


figure %('PaperOrientation','landscape');
%sh(2) = subplot(2,2,2);
sh(2) = subplot(211);
atlantic_psi_max=mitpsi.atlantic_psi;
pacific_psi_max=mitpsi.pacific_psi;
contourf(grd.latg,-grd.zgpsi/1000,atlantic_psi_max'*1e-6,otlev); 
%caxis([-1 1]*max(abs(atlantic_psi_max(:))*1.e-6)); ch(1) = colorbar;
%max(otlev)=ceil(max(abs([atlantic_psi_max(:);pacific_psi_max(:)])));
caxis([-1 1].*max(otlev)); colormap(bluewhitered(length(otlev)-1)); ch(1) = colorbar('YTick',otlev);
set(ch(1),'Position',[0.917, 0.588, 0.03 ,0.337])
hold on;
[cs h2] = contour(grd.latg,-grd.zgpsi/1000,atlantic_psi_max'*1e-6,[0 0]); 
%clh2 = clabel(cs);
hold off
% psimax = max(atlantic_psi_max(:,5:end));
% %iz = find(abs(atlantic_psi(13,:)-psimax)<=1e-4);
% [iy iz] = find(abs(atlantic_psi_max(:)-psimin)<=1e-4);
% text(grd.latg(iy),-grd.zgpsi(iz), ...
%      [num2str(psimax*1e-6,'%5.1f') ' \rightarrow'], ...
%      'horizontalalignment','right')
% psimin = min(min(atlantic_psi_max(1:35,5:end)));
% [iymin,izmin] = find(abs(atlantic_psi_max(:)-psimin)<=1e-4);
% text(grd.latg(iymin),-grd.zgpsi(izmin), ...
%      [num2str(psimin*1e-6,'%5.1f') ' \rightarrow'], ...
%      'horizontalalignment','right')
title('Atlantic overturning streamfunction [Sv]','FontSize',16)
xlabel('Latitude','FontSize',16);ylabel('Depth [km]','FontSize',16)
set(gca,'XLim',[min(grd.latc) max(grd.latc)],'XTick',[-80:20:80],'FontSize',16)
%sh(3) = subplot(2,2,3);
sh(3) = subplot(212);
pacific_psi_max=mitpsi.pacific_psi;
contourf(grd.latg,-grd.zgpsi/1000,pacific_psi_max'*1e-6,otlev); 
caxis([-1 1].*max(otlev)); colormap(bluewhitered(length(otlev)-1)); ch(2) = colorbar;
%caxis([-1 1]*max(abs(pacific_psi_max(:)))*1.e-6); ch(2) = colorbar;
set(ch(2),'Position',[0.917, 0.114, 0.03 ,0.337])
hold on; 
[cs h3] = contour(grd.latg,-grd.zgpsi/1000,pacific_psi_max'*1e-6,[0 0]); 
set(h1,'LineWidth',2,'LineColor','k')
%clh3 = clabel(cs);
hold off
set(gca,'FontSize',16)
title('Pacific overturning streamfunction [Sv]','FontSize',16)
xlabel('Latitude','FontSize',16);ylabel('Depth [km]','FontSize',16)
set(gca,'XLim',[min(grd.latc) max(grd.latc)],'XTick',[-80:20:80],'FontSize',16)
if ~isempty([h2;h3])
  set([h2;h3],'LineWidth',2,'EdgeColor','k');
end
orient landscape
print -dpsc atlantic_pacific_overturning.ps

% clh = [clh1;clh2;clh3];
% if ~isempty(clh)
%   set(clh(2:2:end),'FontSize',8);
% end

% $$$ [cs h] = contourf(grd.long,grd.latg,baro_psi'*1e-6,20); 
% $$$ if ~isempty(h);
% $$$   set(h,'edgecolor','none'); 
% $$$ end; 
% $$$ axis image; 
% $$$ caxis([-1 1]*max(abs(baro_psi(:)))*1.e-6); colorbar('h');
% $$$ title('global barotropic stream function [Sv]')
bstlev = [-200:20:200];
%sh(4) = subplot(2,2,4);
figure
baro_psi_max=mitpsi.gbaro_psi;
contourf(grd.long,grd.latg,baro_psi_max'*1e-6,bstlev); 
caxis([-1 1].*max(bstlev)); colormap(bluewhitered(length(bstlev)-1)); colorbar('YTick',bstlev);
%caxis([-1 1]*max(abs(mitpsi.gbaro_psi(:))*1.e-6));colorbar
%hold on
%[cs h3] = contour(grd.long,grd.latg,baro_psi_max'*1e-6,[0 0]); 
%[cs h]=contour(grd.long,grd.latg,baro_psi_max'*1e-6,bstlev);
%set(h,'edgecolor','k');
% if ~isempty(h); 
%   clh = clabel(cs,h); 
%   set(clh,'Fontsize',8);
% end
%hold off
%axis xy; % axis image; 
title('Global barotropic stream function [Sv]','FontSize',16)
xlabel('Longitude','FontSize',16);ylabel('Latitude','FontSize',16)
set(gca,'Xlim',[min(grd.long) max(grd.lonc)],'XTickLabelMode','manual','XTickMode','manual','XTick',[0:60:300,max(grd.long)],'XTickLabel',[0:60:360],'YLim',[min(grd.latc) max(grd.latc)],'YTick',[-80:20:80],'FontSize',16);
orient landscape
print -dpsc barotropic_streamfunction.ps
%set(sh,'layer','top')
%suptitle(['timestep = ' num2str(timesteps(k)) ', ' tuname ' = ' num2str(tim(k))])
%%
if size(sigma_psi,2)==size(ri,2)
figure
contourf(grd.latc,ri,sigma_psi'.*1e-6,[-50:5:50])
caxis([-1 1].*max(otlev)); colormap(bluewhitered(length(otlev)-1));colorbar;
set(gca,'ydir','rev','clim',[-50 50])
hold on
plot(grd.latc,(nanmin(mit_sigma25(:,:,1))),'color',[0.5 0.5 0.5],'LineWidth',2)% min surface density
plot(grd.latc,(nanmax(mit_sigma25(:,:,1))),'color',[0.5 0.5 0.5],'LineWidth',2)% max surface density
%contour(grd.latc,ri,sigma_psi'*1e-6,[0 0],'LineWidth',2); 

title('Global overturning streamfunction in density space [Sv]','FontSize',16)
set(gca,'FontSize',16)
orient landscape
print -dpsc global_sigma_overturning.ps
else
    disp('Dimension missmatch in sigma overturning calculation')
end
end
