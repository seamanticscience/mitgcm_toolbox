function   []=mit_plotmeansections(grd,tavesteps)

if nargin==0;
    grd=mit_loadgrid;
    grd = mit_oceanmasks(grd);
elseif nargin==1
    tavesteps=mit_timesteps('tave');
end

if ~isfield(grd,'pacific_hfacw')
    grd = mit_oceanmasks(grd);
end

precision = mit_getparm('data','readBinaryPrec');
if isempty(precision); precision = 32.; end
if precision == 32
    acc = 'real*4';
elseif precision == 64;
    acc = 'real*8';
else
    error('readBinaryPrec contains unknown precision')
end

thetaFile = mit_getparm('data','hydrogThetaFile');
if isempty(thetaFile)
    tdatamean = NaN(size(grd.hfacc));
else
    tdatamean = squeeze(mean(mit_readfield(thetaFile,[grd.nx grd.ny grd.nz],acc),4));
end

saltFile = mit_getparm('data','hydrogSaltFile');
if isempty(saltFile)
    sdatamean = NaN(size(grd.hfacc));
else
    sdatamean = squeeze(mean(mit_readfield(saltFile,[grd.nx grd.ny grd.nz],acc),4));
end

if ~strcmp(grd.buoyancy,'OCEANIC');
    sdatamean = sdatamean(:,:,end:-1:1);
    tdatamean = tdatamean(:,:,end:-1:1);
end

if length(tavesteps.tim) > 1
    rac3d = repmat(grd.rac,[1 1 grd.nz]).*grd.hfacc;
    % horizontal averages of t and s in a hovmoeller-type diagram
    mdt = NaN(grd.nz,length(tavesteps.tim));
    mdtt= mdt;
    mdtstd= mdt;
    mdt_atl = mdt;
    mdtt_atl= mdt;
    mdtstd_atl = mdt;
    mdt_pac = mdt;
    mdtt_pac= mdt;
    mdtstd_pac = mdt;
    mdt_sou = mdt;
    mdtt_sou= mdt;
    mdtstd_sou = mdt;
    mds = mdt;
    mdss= mdt;
    mdsstd = mdt;
    mds_atl = mdt;
    mdss_atl= mdt;
    mdsstd_atl = mdt;
    mds_pac = mdt;
    mdss_pac= mdt;
    mdsstd_pac = mdt;
    mds_sou = mdt;
    mdss_sou= mdt;
    mdsstd_sou = mdt;
    indopac_hfacc = ...
        change(change(grd.pacific_hfacc,'==',NaN,0) ...
        +change(grd.indic_hfacc,'==',NaN,0),'==',0,NaN);
    indopac_hfacc(:,1:12,:) = NaN;
    southern_hfacc = grd.hfacc;
    southern_hfacc(:,13:end,:) = NaN;
    atl_hfacc = grd.atlantic_hfacc;
    atl_hfacc(:,1:12,:) = NaN;
    for k=1:length(tavesteps.tim);
        tave=rdmnc(tavesteps.filearr(2:end-1),'Ttave','Stave',tavesteps.timesteps(k));
        t=tave.Ttave;
        s=tave.Stave;
        dt = t(:,:,:)-tdatamean; % anomalies
        ds = s(:,:,:)-sdatamean; % anomalies
        [mdt(:,k) mdtt(:,k) mdtstd(:,k)]    = mit_globalmean(dt,rac3d);
        [mds(:,k) mdss(:,k) mdsstd(:,k)]    = mit_globalmean(ds,rac3d);
        [mdt_atl(:,k) mdtt_atl(:,k) mdtstd_atl(:,k)]    = mit_globalmean(dt,rac3d.*atl_hfacc);
        [mds_atl(:,k) mdss_atl(:,k) mdsstd_atl(:,k)]    = mit_globalmean(ds,rac3d.*atl_hfacc);
        [mdt_pac(:,k) mdtt_pac(:,k) mdtstd_pac(:,k)]    = mit_globalmean(dt,rac3d.*indopac_hfacc);
        [mds_pac(:,k) mdss_pac(:,k) mdsstd_pac(:,k)]    = mit_globalmean(ds,rac3d.*indopac_hfacc);
        [mdt_sou(:,k) mdtt_sou(:,k) mdtstd_sou(:,k)]    = mit_globalmean(dt,rac3d.*southern_hfacc);
        [mds_sou(:,k) mdss_sou(:,k) mdsstd_sou(:,k)]    = mit_globalmean(ds,rac3d.*southern_hfacc);
    end
    
    if length(tavesteps.tim) > 2
        figure
        clear sh
        sh(1) = subplot(2,2,1);
        contourf(tavesteps.tim,-grd.zc,mdt,20);
        caxis([-1 1]*max(abs(caxis))); colormap(bluewhitered(20));colorbar('v')
        title('T-T_{init} horizontally averaged')
        xlabel(['Time [' tavesteps.timeunit ']']); ylabel('Depth [m]')
        sh(2) = subplot(2,2,2);
        contourf(tavesteps.tim,-grd.zc,mds,20);
        caxis([-1 1]*max(abs(caxis))); colormap(bluewhitered(20));colorbar('v')
        title('S-S_{init} horizontally averaged')
        xlabel(['Time [' tavesteps.timeunit ']']); ylabel('Depth [m]')
        %
        sh(3) = subplot(2,2,3);
        contourf(tavesteps.tim,-grd.zc,mdtt,20);
        caxis([-1 1]*max(abs(caxis))); colormap(bluewhitered(20));colorbar('v')
        title('|T-T_{init}| horizontally averaged');
        xlabel(['Time [' tavesteps.timeunit ']']); ylabel('Depth [m]')
        sh(4) = subplot(2,2,4);
        contourf(tavesteps.tim,-grd.zc,mdss,20);
        caxis([-1 1]*max(abs(caxis))); colormap(bluewhitered(20));colorbar('v')
        title('|S-S_{init}| horizontally averaged')
        xlabel(['Time [' tavesteps.timeunit ']']); ylabel('Depth [m]')
        set(sh,'layer','top');
        % suptitle(['experiment ' dname])
        orient landscape
        print -dpsc ts_hovmoeller_anomalies.ps
    end %if length(tim) > 2
    
    k=length(tavesteps.tim);
    figure
    clear sh
    sh(1) = subplot(2,4,1);
    plot(mdt(:,k),-grd.zc,'b-',mdt(:,k)+mdtstd(:,k),-grd.zc,'r--',mdt(:,k)-mdtstd(:,k),-grd.zc,'r--');
    xlabel('T-T_{lev} [degC]'); ylabel('depth [m]'); title('global')
    legend('horiz. mean','std',3,'Location','best');
    sh(2) = subplot(2,4,2);
    plot(mdt_atl(:,k),-grd.zc,'b-',mdt_atl(:,k)+mdtstd_atl(:,k),-grd.zc,'r--',mdt_atl(:,k)-mdtstd_atl(:,k),-grd.zc,'r--');
    xlabel('T-T_{lev} [degC]'); %ylabel('depth [m]');
    title('atlantic')
    sh(3) = subplot(2,4,3);
    plot(mdt_pac(:,k),-grd.zc,'b-',mdt_pac(:,k)+mdtstd_pac(:,k),-grd.zc,'r--',mdt_pac(:,k)-mdtstd_pac(:,k),-grd.zc,'r--');
    xlabel('T-T_{lev} [degC]'); %ylabel('depth [m]');
    title('indo-pacific')
    sh(4) = subplot(2,4,4);
    plot(mdt_sou(:,k),-grd.zc,'b-',mdt_sou(:,k)+mdtstd_sou(:,k),-grd.zc,'r--',mdt_sou(:,k)-mdtstd_sou(:,k),-grd.zc,'r--');
    xlabel('T-T_{lev} [degC]'); ylabel('depth [m]'); title('southern')
    %
    sh(5) = subplot(2,4,5);
    plot(mds(:,k),-grd.zc,'b-',mds(:,k)+mdsstd(:,k),-grd.zc,'r--',mds(:,k)-mdsstd(:,k),-grd.zc,'r--');
    xlabel('S-S_{lev} [PSU]'); ylabel('depth [m]');
    sh(6) = subplot(2,4,6);
    plot(mds_atl(:,k),-grd.zc,'b-',mds_atl(:,k)+mdsstd_atl(:,k),-grd.zc,'r--',mds_atl(:,k)-mdsstd_atl(:,k),-grd.zc,'r--');
    xlabel('S-S_{lev} [PSU]'); %ylabel('depth [m]');
    sh(7) = subplot(2,4,7);
    plot(mds_pac(:,k),-grd.zc,'b-',mds_pac(:,k)+mdsstd_pac(:,k),-grd.zc,'r--',mds_pac(:,k)-mdsstd_pac(:,k),-grd.zc,'r--');
    xlabel('S-S_{lev} [PSU]'); %ylabel('depth [m]');
    sh(8) = subplot(2,4,8);
    plot(mds_sou(:,k),-grd.zc,'b-',mds_sou(:,k)+mdsstd_sou(:,k),-grd.zc,'r--',mds_sou(:,k)-mdsstd_sou(:,k),-grd.zc,'r--');
    xlabel('S-S_{lev} [PSU]'); ylabel('depth [m]');
    set(sh,'xgrid','on','ygrid','on','XLim',[-5 5],'XTickLabelMode','manual','XTick',[-5:1:5],'XTickLabel',[-5:1:5])
    %   set([sh(4); sh(8)],'YAxisLocation','right')
    %     suptitle(['experiment ' dname ', timestep = ' num2str(tim(k)) ...
    % 	      ', ' tuname ' = ' num2str(tim(k))])
    orient landscape
    print -dpsc ts_anomaly_profiles.ps
end %if length(tim) > 1

% zonal mean of temperature and salinity for different ocean basins
tdzm_glo = mit_zonalmean(tdatamean,grd.hfacc,grd.dxc);
tdzm_atl = mit_zonalmean(tdatamean,grd.atlantic_hfacc,grd.dxc);
tdzm_pac = mit_zonalmean(tdatamean,grd.pacific_hfacc,grd.dxc);
tdzm_ind = mit_zonalmean(tdatamean,grd.indic_hfacc,grd.dxc);
sdzm_glo = mit_zonalmean(sdatamean,grd.hfacc,grd.dxc);
sdzm_atl = mit_zonalmean(sdatamean,grd.atlantic_hfacc,grd.dxc);
sdzm_pac = mit_zonalmean(sdatamean,grd.pacific_hfacc,grd.dxc);
sdzm_ind = mit_zonalmean(sdatamean,grd.indic_hfacc,grd.dxc);

clear sh clh
clh = cell(8,1);
for k = length(tavesteps.tim)
    tave=rdmnc(tavesteps.filearr(2:end-1),'Ttave','Stave',tavesteps.timesteps(k));
    t=tave.Ttave;
    s=tave.Stave;
    tzm_glo = mit_zonalmean(t(:,:,:),grd.hfacc,grd.dxc);
    tzm_atl = mit_zonalmean(t(:,:,:),grd.atlantic_hfacc,grd.dxc);
    tzm_pac = mit_zonalmean(t(:,:,:),grd.pacific_hfacc,grd.dxc);
    tzm_ind = mit_zonalmean(t(:,:,:),grd.indic_hfacc,grd.dxc);
    szm_glo = mit_zonalmean(s(:,:,:),grd.hfacc,grd.dxc);
    szm_atl = mit_zonalmean(s(:,:,:),grd.atlantic_hfacc,grd.dxc);
    szm_pac = mit_zonalmean(s(:,:,:),grd.pacific_hfacc,grd.dxc);
    szm_ind = mit_zonalmean(s(:,:,:),grd.indic_hfacc,grd.dxc);
    caxt = [min(tzm_glo(:)-tdzm_glo(:)) max(tzm_glo(:)-tdzm_glo(:))];
    tlev = -6:.5:6;
    if max(abs(caxt)) < 1;
        tlev = -1:.1:1;
    end
    if max(abs(caxt)) < .1;
        tlev = .1*tlev;
    end
    caxs = [min(szm_glo(:)-sdzm_glo(:)) max(szm_glo(:)-sdzm_glo(:))];
    slev = -5.:.1:5;
    if max(abs(caxs)) < 0.2;
        slev = -.20:.01:.20;
    end
    if max(abs(caxs)) < .02;
        slev = .1*slev;
    end
    
    % Plot global t and s sections
    figure
    sh(1) = subplot(211);
    [cs h] = contourf(grd.latc,-grd.zc,(tzm_glo)',[-5:5:40]);
    colormap(parula(length(-5:5:40)));colorbar;
    title('\theta  [degC]: global ocean')
    set(gca,'clim',[-5 40])
    
    
    sh(2) = subplot(212);
    [cs h] = contourf(grd.latc,-grd.zc,(szm_glo)',[15:0.5:40]);
    colormap(parula(length(15:0.5:40)));colorbar;
    title('S  [PSU]: global ocean')
    set(gca,'clim',[30 40])
    colorbar;
    orient landscape
    print -dpsc ts_zonal_sections.ps
    
    % Plot global deviations from levitus climatology
    figure
    sh(1) = subplot(211);
    [cs h] = contourf(grd.latc,-grd.zc,(tzm_glo-tdzm_glo)',tlev);
    title('\theta-\theta_{lev} [degC]: global ocean')
    set(sh(1),'clim',[tlev(1) tlev(end)]);
    caxis([tlev(1) tlev(end)]);colormap(bluewhitered(length(tlev)));colorbar
    
    sh(2) = subplot(212);
    [cs h] = contourf(grd.latc,-grd.zc,(szm_glo-sdzm_glo)',slev);
    title('S-S_{lev} [PSU]: global ocean')
    set(sh(2),'clim',[slev(1) slev(end)])
    caxis([slev(1) slev(end)]);colormap(bluewhitered(length(slev)));colorbar
    orient landscape
    print -dpsc ts_zonal_anomaly_sections.ps
    
    % Plot basin scale sections of t and s
    figure %('PaperPosition',[0.25 0.368552 8 10.2629])
    sh(3) = subplot(321);
    [cs h] = contourf(grd.latc,-grd.zc,(tzm_atl)',[-5:5:40]);
    colormap(parula(length(-5:5:40)));colorbar;
    title('\theta [degC]: atlantic ocean')
    set(gca,'clim',[-5 40])
    sh(5) = subplot(323);
    [cs h] = contourf(grd.latc,-grd.zc,(tzm_pac)',[-5:5:40]);
    colormap(parula(length(-5:5:40)));colorbar;
    title('\theta [degC]: pacific ocean')
    set(gca,'clim',[-5 40])
    sh(7) = subplot(325);
    [cs h] = contourf(grd.latc,-grd.zc,(tzm_ind)',[-5:5:40]);
    colormap(parula(length(-5:5:40)));colorbar;
    title('\theta [degC]: indian ocean')
    set(gca,'clim',[-5 40])
    
    sh(4) = subplot(322);
    [cs h] = contourf(grd.latc,-grd.zc,(szm_atl)',[15:0.5:40]);
    colormap(parula(length(15:0.5:40)));colorbar;set(gca,'clim',[30 40])
    title('S [PSU]: atlantic ocean')
    sh(6) = subplot(324);
    [cs h] = contourf(grd.latc,-grd.zc,(szm_pac)',[15:0.5:40]);
    colormap(parula(length(15:0.5:40)));colorbar;set(gca,'clim',[30 40])
    title('S [PSU]: pacific ocean')
    sh(8) = subplot(326);
    [cs h] = contourf(grd.latc,-grd.zc,(szm_ind)',[15:0.5:40]);
    title('S [PSU]: indian ocean')
    colormap(parula(length(15:0.5:40)));colorbar;set(gca,'clim',[30 40])
    colorbar;
    orient landscape
    print -dpsc ts_basin_sections.ps
    
    % Plot basin scale deviations from levitus climatology
    figure %('PaperPosition',[0.25 0.368552 8 10.2629])
    sh(3) = subplot(321);
    [cs h] = contourf(grd.latc,-grd.zc,(tzm_atl-tdzm_atl)',tlev);
    caxis([tlev(1) tlev(end)]);colormap(bluewhitered(length(tlev)));colorbar
    title('\theta-\theta_{lev} [degC]: atlantic ocean')
    sh(5) = subplot(323);
    [cs h] = contourf(grd.latc,-grd.zc,(tzm_pac-tdzm_pac)',tlev);
    caxis([tlev(1) tlev(end)]);colormap(bluewhitered(length(tlev)));colorbar
    title('\theta-\theta_{lev} [degC]: pacific ocean')
    sh(7) = subplot(325);
    [cs h] = contourf(grd.latc,-grd.zc,(tzm_ind-tdzm_ind)',tlev);
    caxis([tlev(1) tlev(end)]);colormap(bluewhitered(length(tlev)));colorbar
    title('\theta-\theta_{lev} [degC]: indian ocean')
    set(sh(3:2:end),'clim',[tlev(1) tlev(end)])
    
    sh(4) = subplot(322);
    [cs h] = contourf(grd.latc,-grd.zc,(szm_atl-sdzm_atl)',slev);
    caxis([slev(1) slev(end)]);colormap(bluewhitered(length(slev)));colorbar
    title('S-S_{lev} [PSU]: atlantic ocean')
    sh(6) = subplot(324);
    [cs h] = contourf(grd.latc,-grd.zc,(szm_pac-sdzm_pac)',slev);
    caxis([slev(1) slev(end)]);colormap(bluewhitered(length(slev)));colorbar
    title('S-S_{lev} [PSU]: pacific ocean')
    sh(8) = subplot(326);
    [cs h] = contourf(grd.latc,-grd.zc,(szm_ind-sdzm_ind)',slev);
    caxis([slev(1) slev(end)]);colormap(bluewhitered(length(slev)));colorbar
    title('S-S_{lev} [PSU]: indian ocean')
    set(sh(4:2:end),'clim',[slev(1) slev(end)])
    set(sh,'layer','top')
    orient landscape
    print -dpsc ts_basin_anomaly_sections.ps
    %     if ~isempty(cat(1,clh{:}))
    %       set(cat(1,clh{:}),'fontsize',8);
    %     end
    
    %     suptitle(['experiment ' dname ', timestep = ' num2str(tim(k)) ...
    % 	    ', ' tuname ' = ' num2str(tim(k)) ', zonal averages'])
    drawnow;
    
end
clear clh sh

end
