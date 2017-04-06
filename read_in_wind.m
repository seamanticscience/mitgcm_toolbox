% Load in taux data file & rearrange
clear all
taux=nan(128,64,12);
fid=fopen('tren_taux.bin','r','b');taux_bin=fread(fid,'real*4');fclose(fid);
for m=1:12
    for i=1:64
        taux(:,i,m)=taux_bin((8192*(m-1))+((i-1)*128)+1:(8192*(m-1))+i*128);
    end
end

% Masking
factor=0.5; % or 0.66
mask=ones(1,64);
mask(1,10:20)=factor; 
smooth_mask=smooth(mask,5,'moving');
mask_grid=repmat(smooth_mask',[128,1]);

figure;plot(mask_grid(1,:)','r','LineWidth',2);hold on
plot([6 6],[0.5 1.75],'k--',[10 10],[0.5 1.75],'k--',[20 20],[0.5 1.75],'k--',[24 24],[0.5 1.75],'k--')
title('Zonally averaged Zonal Wind Stress Perturbation Mask','FontSize',18,'FontWeight','bold')
xlabel('Latitudinal Grid Cell');ylabel('Multiplication Factor')
set(gca,'YLim',[0.5 1.75]); orient landscape
%print -dpsc2 taux_mask.ps

% remove zeros where the land is...
taux(find(taux==0))=NaN;
taux_gmt=[taux(64:128,:,:);taux(1:63,:,:)];

grd=mit_loadgrid;
surf_hfac=repmat(grd.hfacc(:,:,1),[1,1,12]);
taux_zmean=mit_zonalmean(taux,surf_hfac,grd.dxg);
taux_zannmean=nanmean(taux_zmean');
taux_zannmean=taux_zannmean';

figure;plot(grd.yc(1,:),taux_zmean);hold on
plot(grd.yc(1,:),taux_zannmean(:,:),'k','LineWidth',3)
%plot([6 6],[-0.1 0.25],'k--',[10 10],[-0.1 0.25],'k--',[20 20],[-0.1 0.25],'k--',[24 24],[-0.1 0.25],'k--')
plot([-63.281 -63.281],[-0.1 0.25],'k--',[-54.843 -54.843],[-0.1 0.25],'k--') % Drake Passage...i think.
title('Monthly (thin) and Annual (thick) Zonally averaged Zonal Wind Stress','FontSize',18,'FontWeight','bold')
xlabel('Latitudinal Grid Cell');ylabel('Wind Stress (N m^{-2})')
orient landscape
%print -dpsc2 taux_zonal_ave.ps

taux_annmean=mean(taux,3);
taux_meangmt=[taux_annmean(64:128,:);taux_annmean(1:63,:,:)];
figure;contourf(taux_meangmt');colorbar; hold on
plot([0 128],[6 6],'k--',[0 128],[24 24],'k--')
plot([0 128],[10 10],'k--',[0 128],[20 20],'k--','LineWidth',2);
title('Annual Mean Zonal Wind Stress (N m^{-2})','FontSize',18,'FontWeight','bold')
ylabel('Latitudinal Grid Cell');xlabel('Longitudinal Grid Cell')
orient landscape
%print -dpsc2 taux_yr_ave.ps

% TURNS OUT THIS IS NOT HOW I HAVE DONE IT FOR THE PRODUCTION RUNS!!!
go=input('Continue to process and output wind file? (y/n)\n>>','s');
if go=='y';
    % Multiply taux by mask_grid, only inc/dec Westerly wind ~~~~> that
    % way :D
    new_wind=nan(128,64,12);
    for m=1:12
        for y=1:64
            for x=1:128
                if taux(x,y,m)>0;
                    new_wind(x,y,m)=mask_grid(x,y).*taux(x,y,m);
                else
                    new_wind(x,y,m)=taux(x,y,m);
                end
            end
        end
    end
    
    new_annmean=mean(new_wind,3);
    new_meangmt=[new_annmean(64:128,:);new_annmean(1:63,:,:)];
    diff=new_meangmt-taux_meangmt;
    new_wind_gmt=[new_wind(64:128,:,:);new_wind(1:63,:,:)];
    
    figure;contourf(diff');orient landscape;colorbar;hold on
    plot([0 128],[6 6],'k--',[0 128],[24 24],'k--')
    plot([0 128],[10 10],'k--',[0 128],[20 20],'k--','LineWidth',2);
    title('Annual Mean Zonal Wind Stress Perturbation (N m^{-2})','FontSize',18,'FontWeight','bold')
    ylabel('Latitudinal Grid Cell');xlabel('Longitudinal Grid Cell')
 %   print -dpsc2 taux_perturb_ave.ps
   
    diff_monthly=nan(128,64,12);
    figure
    for p=1:12
        diff_monthly(:,:,p)=new_wind_gmt(:,:,p)-taux_gmt(:,:,p);
        eval(['AX',num2str(p),'=subplot(4,3,',num2str(p),');']);
        contourf(diff_monthly(:,:,p)');
        if factor>1;
            set(gca,'CLim',[0 0.2]);
        else
            set(gca,'CLim',[-0.14,0.0]);
        end
        name=['Month ',num2str(p)]; title(name)
    end
    suptitle('Monthly Zonal Wind Stress Field Perturbations (N m^{-2})');
    orient landscape
    colorbar1=colorbar('peer',AX9,[0.915,0.098,0.018,0.733],'Box','on','Location','manual');
  %  print -dpsc2 taux_monthly_perturb.ps
    
    diff_monthly=nan(128,64,12);
    figure
    for p=1:12
        diff_monthly(:,:,p)=new_wind_gmt(:,:,p)./taux_gmt(:,:,p);
        eval(['AX',num2str(p),'=subplot(4,3,',num2str(p),');']);
        
        if factor>1;
            contourf(diff_monthly(:,:,p)',[1:0.5/3:1.5]);
            set(gca,'CLim',[1 1.5]);
            colormap(jet(3))
        else
            contourf(diff_monthly(:,:,p)',[0.5:0.5/3:1]);
            set(gca,'CLim',[0.5 1]);
            colormap(jet(3))
        end
        name=['Month ',num2str(p)]; title(name)
    end
    suptitle('Monthly Zonal Wind Stress Field Perturbation Ratio (N m^{-2})');
    orient landscape
    colorbar1=colorbar('peer',AX9,[0.915,0.098,0.018,0.733],'Box','on','Location','manual');
  %  print -dpsc2 taux_monthly_perturb.ps
  
    new_zmean=mit_zonalmean(new_wind,surf_hfac,grd.dxg);
    new_zannmean=nanmean(new_zmean');
    new_zannmean=new_zannmean';
    
    figure;plot(grd.yc(1,:),new_zmean);hold on
    plot(grd.yc(1,:),new_zannmean(:,:),'k','LineWidth',3)
    if factor>1;
        plot([6 6],[-0.1 0.35],'k--',[10 10],[-0.1 0.35],'k--',[20 20],[-0.1 0.35],'k--',[24 24],[-0.1 0.35],'k--')
    else
        plot([6 6],[-0.1 0.25],'k--',[10 10],[-0.1 0.25],'k--',[20 20],[-0.1 0.25],'k--',[24 24],[-0.1 0.25],'k--')
    end
    title('Monthly (thin) and Annual (thick) Zonally averaged Zonal Wind Stress','FontSize',18,'FontWeight','bold')
    xlabel('Latitudinal Grid Cell');ylabel('Wind Stress (N m^{-2})')
    orient landscape
   % print -dpsc2 taux_perturb_zonal_ave.ps

    
    sure=input('Are you sure you want to save this data? (y/n)\n>>','s');
    if sure=='y'
        fid=fopen('inc_taux_le.bin','w','l');fwrite(fid,new_wind,'real*4'); fclose(fid);
        disp('Save successful');
    end
    
end

