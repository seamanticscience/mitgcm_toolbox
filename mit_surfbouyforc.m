% Read in MITgcm data
% Load grid
grid=mit_loadgrid('/Volumes/Postdoc_Data/controlrun/controlrun_anomplotref/');

file=nix_find('/Users/jml1/MITgcm/MITgcm_63m/perturb_eddies/[a-z]*','-name "tave.2016*"');
%filea=[]; %nix_find('/Volumes/Postdoc_Data/perturb_eddies/[a-z]*','-name "tave.*"');
%file=[{'/Volumes/Postdoc_Data/controlrun/controlrun_anomplotref/tave.20k.glob.nc'};filea;fileb];

clear filea fileb

for ff=1:length(file)
    [fpath,fname,fext]=fileparts(file{ff});
    cd(fpath)
        
    mitgcm_tavefile     = [fname,fext];
    mitgcm_surfdiagfile = [strrep(fname,'tave','surfDiag'),fext];
    
    %% Load data
    tave=rdmnc(mitgcm_tavefile,'Ttave','Stave','uVeltave','vVeltave','T','X','Y','Z');
    surf=rdmnc(mitgcm_surfdiagfile,'TFLUX','oceQnet','TRELAX','SFLUX','oceFWflx','SRELAX','T','X','Y');
    %mitgcm_dicfile=strrep(mitgcm_tavefile,'tave','ptr_tave');
    %ptracers=rdmnc(mitgcm_dicfile,'dic','po4');
    
    if length(tave.T)>1
        if ~exist('tave.20k.glob.nc','file')
            !ln -sf /Volumes/Postdoc_Data/controlrun/controlrun_anomplotref/tave.20k.glob.nc .
            !ln -sf /Volumes/Postdoc_Data/controlrun/controlrun_anomplotref/surfDiag.20k.glob.nc .
        end
        
        cntrl_tave=rdmnc('tave.20k.glob.nc','Ttave','Stave','uVeltave','vVeltave','T','X','Y','Z');
        cntrl_surf=rdmnc('surfDiag.20k.glob.nc','TFLUX','oceQnet','TRELAX','SFLUX','oceFWflx','SRELAX','T','X','Y');
    end
    
    landmask=repmat(grid.hfacc(:,:,1),[1,1,length(tave.T)]);
    
    % transport through Drake Passage
    tdp = NaN*ones(length(tave.T),1);
    tdpz=NaN*ones(grid.nz,length(tave.T));
    kx = 105;
    kyg = 5:20;
    da = grid.dz*grid.dyg(kx,kyg);
    
    for t=1:length(tave.T);
        if ~strcmp(grid.buoyancy,'OCEANIC');
            u = tave.uVeltave(:,:,end:-1:1,t);
            v = tave.vVeltave(:,:,end:-1:1,t);
        else
            u = tave.uVeltave(:,:,:,t);
            v = tave.vVeltave(:,:,:,t);
        end
        
        tdp(t) = sum(nansum(squeeze(u(kx,kyg,:))'.*da))*1e-6; % in Sv
        tdpz(:,t) = nansum(squeeze(u(kx,kyg,:)).*da')'*1e-6; % in Sv
    end
    
    % load barotropic streamfunction data to delimit ACC
    %psi_data=rdmnc('psi_data.nc','gbaro_psi','Xpsi','Ypsi','T');
    gtmp=nc_varget('psi_data.nc','gbaro_psi');
    gT=nc_varget('psi_data.nc','T');
    if length(gT)==1
        gbaro_psi=repmat(permute(gtmp,[3,2,1])./1e6,[1,1,length(tave.T)]);
    else
        gbaro_psi=permute(gtmp,[3,2,1])./1e6;
    end
    clear gtmp
    gbaro_psi=squeeze(gbaro_psi);
    
    psimask=landmask;
    for t=1:length(tave.T)
        psimask(:,25:64,t)=NaN;
        for i=1:grid.nx
            for j=1:grid.ny
                if gbaro_psi(i,j,t)<24; psimask(i,j,t)=NaN; end
                if gbaro_psi(i,j,t)>=floor(tdp(t)); psimask(i,j,t)=NaN; end
            end
        end
    end
    
    % figure
    % contourf(grid.lonc,grid.latc,gbaro_psi(:,:,end)')
    % canom;colormap(bluewhitered(20));colorbar
    % hold on
    % contour(grid.lonc,grid.latc,gbaro_psi(:,:,end)',[floor(tdp(end)) floor(tdp(end))],'g')
    % contour(grid.lonc,grid.latc(1:25),gbaro_psi(:,1:25,end)',[0 0],'g')
    
    % surForcT = oceQnet + TRELAX - oceQsw(=0),
    % surForcS = oceSflux(=0) + SRELAX + EmPmR*35, with (EmPmR = -1*oceFWflx)
    % BUT
    % TFLUX = surForcT + oceQsw(=0) + oceFreez(<<)
    % SFLUX = surForcS
    % so TFLUX and SFLUX collect all buoyancy terms
    
    % nc{'SFLUX'}.description = ncchar(''total salt flux (match salt-content variations), >0 increases salt'');
    % nc{'SFLUX'}.units = ncchar(''g/m^2/s'');
    %
    % SFLUX is virtual salt flux added to the surface layer, but need the quantity of FW instead.
    % EmPmR input forcing is given as m/s and is converted to kg/m2/s by
    % multiplying by rhoConstFresh = 1035 (default) in external_fields_load.F
    % and then converted to salt flux by multiplying by convertFW2Salt = 35
    % (default) in external_forcing_surf.F
    % g/m2/s * m3/kg * kg/g --> m/s or *1/rhoConstFresh*convertFW2Salt
    %
    % nc{'oceFWflx'}.description = ncchar(''net surface Fresh-Water flux into the ocean (+=down), >0 decreases salinity'');
    % nc{'oceFWflx'}.units = ncchar(''kg/m^2/s'');
    % nc{'SRELAX'}.description = ncchar(''surface salinity relaxation, >0 increases salt'');
    % nc{'SRELAX'}.units = ncchar(''g/m^2/s'');
    
    % nc{'TFLUX'}.description = ncchar(''total heat flux (match heat-content variations), >0 increases theta'');
    % nc{'TFLUX'}.units = ncchar(''W/m^2'');
    % nc{'oceQnet'}.description = ncchar(''net surface heat flux into the ocean (+=down), >0 increases theta'');
    % nc{'oceQnet'}.units = ncchar(''W/m^2'');
    % nc{'TRELAX'}.description = ncchar(''surface temperature relaxation, >0 increases theta'');
    % nc{'TRELAX'}.units = ncchar(''W/m^2'');
    
    mit_tsurf=squeeze(tave.Ttave(:,:,1,:)).*landmask;
    mit_ssurf=squeeze(tave.Stave(:,:,1,:)).*landmask;
    % mit_sigma0=sw_pden(mit_s0,mit_t0,zeros(length(tave.Y),length(tave.X),length(tave.T)),zeros(length(tave.Y),length(tave.X),length(tave.T))) - 1000;
    % MITgcm outputs potential temperature as Ttave, so should use sw_dens,
    % where P is zero
    mit_sigmasurf=sw_dens(mit_ssurf,mit_tsurf,ones(length(tave.X),length(tave.Y),length(tave.T)).*grid.zc(1)); % Kg m-3
    mit_ht=squeeze(surf.TFLUX(:,:,1,:)).*landmask; % W m-2
    mit_fw=landmask.*squeeze(surf.SFLUX(:,:,1,:))./(1035*35); %g SALT m-2 s-1 -> m s-1
    
    mit_qn=squeeze(surf.oceQnet(:,:,1,:)).*landmask; % W m-2
    mit_tr=landmask.*squeeze(surf.TRELAX(:,:,1,:)); % W m-2
    mit_nf=landmask.*squeeze(surf.oceFWflx(:,:,1,:))./(-1035); % Kg FW m-2 s-1 -> m s-1
    mit_sr=landmask.*squeeze(surf.SRELAX(:,:,1,:))./(1035*35); % g SALT m-2 s-1 -> m s-1
    
    mit_a = NaN(grid.nx,grid.ny,length(tave.T));
    mit_b = NaN(grid.nx,grid.ny,length(tave.T));
    mit_c = NaN(grid.nx,grid.ny,length(tave.T));
    
    mit_bflux_t = NaN(grid.nx,grid.ny,length(tave.T));
    mit_bflux_s = NaN(grid.nx,grid.ny,length(tave.T));
    mit_bflux   = NaN(grid.nx,grid.ny,length(tave.T));
    
    mit_bflux_tf = NaN(grid.nx,grid.ny,length(tave.T));
    mit_bflux_sf = NaN(grid.nx,grid.ny,length(tave.T));
    mit_bflux_f  = NaN(grid.nx,grid.ny,length(tave.T));
    
    mit_bflux_tr = NaN(grid.nx,grid.ny,length(tave.T));
    mit_bflux_sr = NaN(grid.nx,grid.ny,length(tave.T));
    mit_bflux_r  = NaN(grid.nx,grid.ny,length(tave.T));
    
    for m=1:1:length(tave.T)
        g=9.81; % m s-2
        % You get a slightly different answer with fixed rather than variable
        % alpha, beta and cp's!
        %     mit_a=ones(grid.nx,grid.ny,length(tave.T)).*0.0002;
        %     mit_b=ones(grid.nx,grid.ny,length(tave.T)).*0.00074;
        %     mit_c=ones(grid.nx,grid.ny,length(tave.T)).*3994;
        
        mit_a(:,:,m)=sw_alpha(mit_ssurf(:,:,m),mit_tsurf(:,:,m),grid.zc(1),'ptmp'); % K-1
        mit_b(:,:,m)=sw_beta(mit_ssurf(:,:,m),mit_tsurf(:,:,m),grid.zc(1),'ptmp');  %psu-1
        mit_c(:,:,m)=sw_cp(mit_ssurf(:,:,m),mit_tsurf(:,:,m),grid.zc(1));
        %
        mit_d1= -g.*mit_a(:,:,m)./(mit_c(:,:,m).*mit_sigmasurf(:,:,m));
        mit_d2= g.*mit_b(:,:,m).*mit_ssurf(:,:,m); % already in *s-1 not *yr-1
        %
        mit_bflux_t(:,:,m) =-mit_d1.*mit_ht(:,:,m);
        mit_bflux_s(:,:,m) =-mit_d2.*mit_fw(:,:,m);
        mit_bflux(:,:,m)   = mit_bflux_t(:,:,m) + mit_bflux_s(:,:,m);
        
        mit_bflux_tf(:,:,m) =-mit_d1.*mit_qn(:,:,m);
        mit_bflux_sf(:,:,m) =-mit_d2.*mit_nf(:,:,m);
        mit_bflux_f(:,:,m)  = mit_bflux_tf(:,:,m) + mit_bflux_sf(:,:,m);
        
        mit_bflux_tr(:,:,m) =-mit_d1.*mit_tr(:,:,m);
        mit_bflux_sr(:,:,m) =-mit_d2.*mit_sr(:,:,m);
        mit_bflux_r(:,:,m)  = mit_bflux_tr(:,:,m) + mit_bflux_sr(:,:,m);
    end
    
    if length(tave.T)>1
        cntrl_tsurf=squeeze(cntrl_tave.Ttave(:,:,1)).*grid.hfacc(:,:,1);
        cntrl_ssurf=squeeze(cntrl_tave.Stave(:,:,1)).*grid.hfacc(:,:,1);
        % mit_sigma0=sw_pden(mit_s0,mit_t0,zeros(length(tave.Y),length(tave.X),length(tave.T)),zeros(length(tave.Y),length(tave.X),length(tave.T))) - 1000;
        % MITgcm outputs potential temperature as Ttave, so should use sw_dens,
        % where P is zero
        cntrl_sigmasurf=sw_dens(cntrl_ssurf,cntrl_tsurf,ones(grid.nx,grid.ny).*grid.zc(1)); % Kg m-3
        cntrl_ht=squeeze(cntrl_surf.TFLUX(:,:,1)).*grid.hfacc(:,:,1); % W m-2
        cntrl_fw=grid.hfacc(:,:,1).*squeeze(cntrl_surf.SFLUX(:,:,1))./(1035*35); %g SALT m-2 s-1 -> m s-1
        
        cntrl_qn=squeeze(cntrl_surf.oceQnet(:,:,1)).*grid.hfacc(:,:,1); % W m-2
        cntrl_tr=grid.hfacc(:,:,1).*squeeze(cntrl_surf.TRELAX(:,:,1)); % W m-2
        cntrl_nf=grid.hfacc(:,:,1).*squeeze(cntrl_surf.oceFWflx(:,:,1))./(-1035); % Kg FW m-2 s-1 -> m s-1
        cntrl_sr=grid.hfacc(:,:,1).*squeeze(cntrl_surf.SRELAX(:,:,1))./(1035*35); % g SALT m-2 s-1 -> m s-1
        
        g=9.81; % m s-2
        
        cntrl_a=sw_alpha(cntrl_ssurf,cntrl_tsurf,grid.zc(1),'ptmp'); % K-1
        cntrl_b=sw_beta(cntrl_ssurf,cntrl_tsurf,grid.zc(1),'ptmp');  %psu-1
        cntrl_c=sw_cp(cntrl_ssurf,cntrl_tsurf,grid.zc(1));
        %
        cntrl_d1= -g.*cntrl_a./(cntrl_c.*cntrl_sigmasurf);
        cntrl_d2= g.*cntrl_b.*cntrl_ssurf; % already in *s-1 not *yr-1
        %
        cntrl_bflux_t =-cntrl_d1.*cntrl_ht;
        cntrl_bflux_s =-cntrl_d2.*cntrl_fw;
        cntrl_bflux   = cntrl_bflux_t + cntrl_bflux_s;
        
        cntrl_bflux_tf =-cntrl_d1.*cntrl_qn;
        cntrl_bflux_sf =-cntrl_d2.*cntrl_nf;
        cntrl_bflux_f  = cntrl_bflux_tf + cntrl_bflux_sf;
        
        cntrl_bflux_tr = -cntrl_d1.*cntrl_tr;
        cntrl_bflux_sr = -cntrl_d2.*cntrl_sr;
        cntrl_bflux_r  = cntrl_bflux_tr + cntrl_bflux_sr;
    end
    
    figure
    plot(grid.latc,mit_zonalmean(mit_bflux(:,:,end),grid.hfacc(:,:,1),grid.dxc),'k','Linewidth',2)
    hold on
    plot(grid.latc,mit_zonalmean(mit_bflux_tf(:,:,end),grid.hfacc(:,:,1),grid.dxc),'r','Linewidth',2)
    plot(grid.latc,mit_zonalmean(mit_bflux_sf(:,:,end),grid.hfacc(:,:,1),grid.dxc),'b','Linewidth',2)
    plot(grid.latc,mit_zonalmean(mit_bflux_tr(:,:,end),grid.hfacc(:,:,1),grid.dxc),'m','Linewidth',2)
    plot(grid.latc,mit_zonalmean(mit_bflux_sr(:,:,end),grid.hfacc(:,:,1),grid.dxc),'c','Linewidth',2)
    plot([-80 80],[0 0],'k--')
    if length(tave.T)>1
        plot(grid.latc,mit_zonalmean(cntrl_bflux,grid.hfacc(:,:,1),grid.dxc),'k--','Linewidth',1)
        plot(grid.latc,mit_zonalmean(cntrl_bflux_tf,grid.hfacc(:,:,1),grid.dxc),'r--','Linewidth',1)
        plot(grid.latc,mit_zonalmean(cntrl_bflux_sf,grid.hfacc(:,:,1),grid.dxc),'b--','Linewidth',1)
        plot(grid.latc,mit_zonalmean(cntrl_bflux_tr,grid.hfacc(:,:,1),grid.dxc),'m--','Linewidth',1)
        plot(grid.latc,mit_zonalmean(cntrl_bflux_sr,grid.hfacc(:,:,1),grid.dxc),'c--','Linewidth',1)
    end
    legend('Total Bouyancy Flux','Qnet Buoyancy Flux','FW Buoyancy Flux','TRELAX Buoyancy Flux','SRELAX Buoyancy Flux','Location','SouthWest')
    set(gca,'XLim',[-80 80],'YLim',[-5e-8 5e-8],'FontSize',14)
    title('Surface Buoyancy Fluxes, positive is input of buoyancy to the ocean.','FontSize',16)
    orient landscape
    print -dpsc2 diag_mean_buoyancy_flux.ps
    fixpslinestyle('diag_mean_buoyancy_flux.ps')
    
    acc_bforc=NaN(length(tave.T),1);
    for t=1:length(tave.T)
        bftmp=mit_bflux(:,:,t);
        pmasktmp=psimask(:,:,t);
        acc_bforc(t)=nansum(bftmp(:).*grid.rac(:).*pmasktmp(:))./nansum(grid.rac(:).*pmasktmp(:));
    end
    
    % if length(tave.T)>1
    %    figure
    %    plot(tave.T./(60*60*24*360),acc_bforc)
    %    disp(pwd)
    %    report('Average buoyancy Forcing over ACC circumpolar contours is: %03d m2s-3\n',acc_bforc(end))
    % else
    disp(pwd)
    report('Average buoyancy Forcing over ACC circumpolar contours is: %03d m2s-3\n',acc_bforc(end))
    %end
    clear functions
end