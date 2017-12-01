function []=mit_plotbiogeochem(grd,tavesteps)
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

nptr=mit_getparm('data.ptracers','PTRACERS_numInUse');
ptr_tave=rdmnc(['ptr_',tavesteps.filearr(2:end-1)],tavesteps.timesteps(tavesteps.kmax));

% dic_bar=1.8:0.02:2.4;
% alk_bar=2.0:0.025:2.5;
po4_bar=-2.5e-4:2.5e-4:3.75e-3;
% dop_bar=-2e-5:2e-5:2.4e-4;
o2_bar =-0.02:0.02:0.42;
fe_bar =-1e-7:1e-7:1.4e-6;
cdis_bar= -32e-3:2e-3:32e-3;
ligand_bar=-10e-7:10e-7:11e-6;
% 
% anom_dic_bar=-30:5:30;
% anom_alk_bar=-30:5:30;
% anom_po4_bar=-0.05:0.005:0.05;
% anom_dop_bar=-2.5e-2:5e-3:2.5e-2;
% anom_o2_bar =-25:5:25;
% anom_fe_bar =-1e-3:1e-4:1e-3;

for i=7 %1:nptr
    ptr_name=mit_getparm('data.ptracers',['PTRACERS_names(',num2str(i),')']);
    ptr=ptr_tave.(ptr_name);
    if exist([ptr_name,'_bar'],'var');
        ptr_bar=eval([ptr_name,'_bar']);
    else
%        ptr_mean    = nansum(ptr(:).*grd.volc(:).*grd.hfacc(:))./nansum(grd.volc(:).*grd.hfacc(:));
        ptr_mean    = nanmean(ptr(:).*grd.hfacc(:));
        ptr_std     = nanstd(ptr(:).*grd.hfacc(:)); % I know this isnt quite right, but best quick and dirty
        ptr_mag     = 10^floor(log10(abs(ptr_mean)));
        ptr_bar_min = floor((ptr_mean-2.*ptr_std)./ptr_mag).*ptr_mag;
        ptr_bar_max = ceil((ptr_mean+2.*ptr_std)./ptr_mag).*ptr_mag;
        ptr_min = floor(nanmin(ptr(:).*grd.hfacc(:))./ptr_mag).*ptr_mag;
        ptr_max = ceil(nanmax(ptr(:).*grd.hfacc(:))./ptr_mag).*ptr_mag;
        
        if (ptr_mag == 1 )
            ptr_mag=ptr_mag/10;
            ptr_bar_min=floor(nanmin(ptr(:).*grd.hfacc(:))./ptr_mag).*ptr_mag;
            ptr_bar_max=ceil(nanmax(ptr(:).*grd.hfacc(:))./ptr_mag).*ptr_mag;
        end
        
        ptr_bar_int=(ptr_bar_max-ptr_bar_min)/((ptr_bar_max-ptr_bar_min)./(ptr_mag));
        ptr_bar=[ptr_min,ptr_bar_min:ptr_bar_int:ptr_bar_max,ptr_max];
        
        if length(ptr_bar)<12
            ptr_bar=[ptr_min,ptr_bar_min:ptr_bar_int/4:ptr_bar_max,ptr_max];
        elseif length(ptr_bar)>50
            ptr_bar_int=ptr_bar_int*10;
%            ptr_min = floor(nanmin(ptr(:).*grd.hfacc(:))./ptr_mag).*ptr_mag;
%            ptr_max = ceil(nanmax(ptr(:).*grd.hfacc(:))./ptr_mag).*ptr_mag;
%            ptr_bar_min=floor(nanmin(ptr(:).*grd.hfacc(:))./ptr_mag).*ptr_mag;
%            ptr_bar_max=ceil(nanmax(ptr(:).*grd.hfacc(:))./ptr_mag).*ptr_mag;
            ptr_bar=[ptr_min,ptr_bar_min:ptr_bar_int:ptr_bar_max,ptr_max];
        end
    end
    
%% surface plot
    figure
    h=axes;
    [~,a]=contourf(grd.lonc,grd.latc,ptr(:,:,1)'.*grd.hfacc(:,:,1)',ptr_bar);
    caxis([ptr_bar(2) ptr_bar(end-1)]);colormap(parula(length(ptr_bar)-3));colorbar('YTick',ptr_bar(2:2:end-1));
    %caxis([min(ptr_bar) max(ptr_bar)]);colorbar;colormap(parula(length(ptr_bar)-1))
    title(['Surface ',ptr_name,' concentration [mol/m3]'],'FontSize',16);xlabel('Longitude','FontSize',16);ylabel('Latitude','FontSize',16)
    set(h,'Xlim',[min(grd.long) max(grd.lonc)],'XTickLabelMode','manual','XTickMode','manual','XTick',[0:120:300,max(grd.long)],'XTickLabel',[0:120:360],...
        'YLim',[min(grd.latc) max(grd.latc)],'YTick',[-80:40:80],'FontSize',16);
    orient landscape
    eval(['print -dpsc gchem_',ptr_name,'_surface.ps'])

%% Atlantic and Pacific Sections
    figure
    h(1)=subplot(211);
    [~,a(1)]=contourf(grd.latc,-grd.zc,squeeze(nanmean(ptr.*grd.atlantic_hfacc,1))',ptr_bar);
    caxis([ptr_bar(2) ptr_bar(end-1)]);colormap(parula(length(ptr_bar)-3));colorbar('YTick',ptr_bar(2:2:end-1));
    %caxis([min(ptr_bar) max(ptr_bar)]);colorbar;colormap(parula(length(ptr_bar)-1))
    title(['Atlantic ',ptr_name,' concentration [mol/m3]'],'FontSize',16);xlabel('Latitude','FontSize',16);ylabel('Depth [km]','FontSize',16)

    h(2)=subplot(212);
    [~,a(2)]=contourf(grd.latc,-grd.zc,squeeze(nanmean(ptr.*grd.pacific_hfacc,1))',ptr_bar);
    caxis([ptr_bar(2) ptr_bar(end-1)]);colormap(parula(length(ptr_bar)-3));colorbar('YTick',ptr_bar(2:2:end-1));
    %caxis([min(ptr_bar) max(ptr_bar)]);colorbar;colormap(parula(length(ptr_bar)-1))
    title(['Pacific ',ptr_name,' concentration [mol/m3]'],'FontSize',16);xlabel('Longitude','FontSize',16);ylabel('Depth [km]','FontSize',16)
    
    set(h,'Xlim',[min(grd.latc) max(grd.latc)],'XTick',[-80:40:80],...
        'YLim',[-5000 min(grd.zc)],'YTickLabelMode','manual','YTickMode','manual',...
        'YTick',[-5000:500:-500,min(grd.zc)],'YTickLabel',[-5:0.5:0],'FontSize',16);
    orient landscape
    eval(['print -dpsc gchem_',ptr_name,'_sections.ps'])


    ptr_init_name=mit_getparm('data.ptracers',['PTRACERS_initialFile(',num2str(i),')']);
    
    if exist(ptr_init_name,'file')
        fid=fopen(ptr_init_name,'r','b');cntrl_ptr=fread(fid,'real*4');fclose(fid);
        cntrl_ptr=reshape(cntrl_ptr,[grd.nx,grd.ny,grd.nz]);
        ptr_anom=(ptr-cntrl_ptr).*1000;
        if exist(['anom_',ptr_name,'_bar'],'var');
            anom_bar=eval(['anom_',ptr_name,'_bar']);
        else
            anom_abs=max([0.5*nanmax(ptr_anom(:).*grd.hfacc(:)),abs(0.5*nanmin(ptr_anom(:).*grd.hfacc(:)))]);
            anom_mag=10^floor(log10(abs(anom_abs)));
            anom_std     = nanstd(ptr_anom(:).*grd.hfacc(:)); % I know this isnt quite right, but best quick and dirty
            anom_bar_max = ceil (( 2.*anom_std)./anom_mag).*anom_mag;
            anom_bar_min = -anom_bar_max;
            anom_max=ceil((max([nanmax(ptr_anom(:).*grd.hfacc(:)),abs(nanmin(ptr_anom(:).*grd.hfacc(:)))]))./anom_mag).*anom_mag;
            anom_min=-anom_max;
            anom_bar_int=(anom_bar_max-anom_bar_min)/((anom_bar_max-anom_bar_min)./anom_mag);
                       
            anom_bar=[anom_min,anom_bar_min:anom_bar_int:anom_bar_max,anom_max];

            if length(anom_bar)<=14
                  anom_bar=[anom_min,anom_bar_min:anom_bar_int/4:anom_bar_max,anom_max];
            elseif length(anom_bar)>50
                  anom_bar=[anom_min,anom_bar_min:anom_bar_int*10:anom_bar_max,anom_max];
            end

%             anom_bar=-anom_mag.*ceil(max([nanmax(ptr_anom(:).*grd.hfacc(:)./anom_mag),abs(nanmin(ptr_anom(:).*grd.hfacc(:)./anom_mag))])):...
%                 (anom_mag/10).*ceil(max([nanmax(ptr_anom(:).*grd.hfacc(:)./anom_mag),abs(nanmin(ptr_anom(:).*grd.hfacc(:)./anom_mag))])):...
%                 anom_mag.*ceil(max([nanmax(ptr_anom(:).*grd.hfacc(:)./anom_mag),abs(nanmin(ptr_anom(:).*grd.hfacc(:)./anom_mag))]));
        end

%% Surface anomaly
        figure
        h=axes;
        [~,a]=contourf(grd.lonc,grd.latc,ptr_anom(:,:,1)'.*grd.hfacc(:,:,1)',anom_bar);
        caxis([anom_bar(2) anom_bar(end-1)]);colormap(bluewhitered(length(anom_bar)-3));colorbar('YTick',anom_bar(2:2:end-1));
        %caxis([min(ptr_bar) max(ptr_bar)]);colorbar;colormap(parula(length(ptr_bar)-1))
        title(['Surface ',ptr_name,' anomaly [mmol/m3]'],'FontSize',16);xlabel('Longitude','FontSize',16);ylabel('Latitude','FontSize',16)
        set(h,'Xlim',[min(grd.long) max(grd.lonc)],'XTickLabelMode','manual','XTickMode','manual','XTick',[0:120:300,max(grd.long)],'XTickLabel',[0:120:360],...
            'YLim',[min(grd.latc) max(grd.latc)],'YTick',[-80:40:80],'FontSize',16);
        orient landscape
        eval(['print -dpsc gchem_',ptr_name,'_surface_anomaly.ps'])

%% Atlantic and Pacific Anomaly Sections
        figure
        h(1)=subplot(211);
        [~,a(1)]=contourf(grd.latc,-grd.zc,squeeze(nanmean(ptr_anom.*grd.atlantic_hfacc,1))',anom_bar);
        caxis([anom_bar(2) anom_bar(end-1)]);colormap(bluewhitered(length(anom_bar)-3));colorbar('YTick',anom_bar(2:2:end-1));
        %caxis([min(ptr_bar) max(ptr_bar)]);colorbar;colormap(parula(length(ptr_bar)-1))
        title(['Atlantic ',ptr_name,' anomaly [mmol/m3]'],'FontSize',16);xlabel('Latitude','FontSize',16);ylabel('Depth [km]','FontSize',16)
        
        h(2)=subplot(212);
        [~,a(2)]=contourf(grd.latc,-grd.zc,squeeze(nanmean(ptr_anom.*grd.pacific_hfacc,1))',anom_bar);
        caxis([anom_bar(2) anom_bar(end-1)]);colormap(bluewhitered(length(anom_bar)-3));colorbar('YTick',anom_bar(2:2:end-1));
        %caxis([min(ptr_bar) max(ptr_bar)]);colorbar;colormap(parula(length(ptr_bar)-1))
        title(['Pacific ',ptr_name,' anomaly [mmol/m3]'],'FontSize',16);xlabel('Longitude','FontSize',16);ylabel('Depth [km]','FontSize',16)
        
        set(h,'Xlim',[min(grd.latc) max(grd.latc)],'XTick',[-80:40:80],...
            'YLim',[-5000 min(grd.zc)],'YTickLabelMode','manual','YTickMode','manual',...
            'YTick',[-5000:500:-500,min(grd.zc)],'YTickLabel',[-5:0.5:0],'FontSize',16);
        orient landscape
        eval(['print -dpsc gchem_',ptr_name,'_section_anomaly.ps'])
    end
end

