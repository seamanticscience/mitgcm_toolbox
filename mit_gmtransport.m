function [ubg,vbg,wbg]=mit_gmtransport(isAdvForm,gm_psix,gm_psiy,nx,ny,nz,nt,dx,dy,cmask,umask,vmask)
% function [ubg,vbg,zbg]=mit_gmtransport(kwx,kwy,kwz,nx,ny,nz,dz)
%
% Caluculate MITgcm Bolus advective transport for the GM-eddy mixing
% parameterisation. Kwx,Kwy and Kwz are Kgm*Isopycnal_Slope (skew flux) or
% use GM_PSIX and GM_PSIY for the advective form.

ubg=NaN(nx,ny,nz,nt);
vbg=NaN(nx,ny,nz,nt);
wbg=NaN(nx,ny,nz,nt);

%gmadvtest=mit_getparm(data_gmredi,'GM_AdvForm');

%if strcmpi(gmadvtest,'true') % Using advective form
if isAdvForm==1;
% Setting of KGM. Kwx and Kwy are the slope*thickness diffusivity so why are we multiplying by KGM again?
    % This is necessary because the gmredi tensor when the skew flux is used and Kgm=Krho is:
    %    x,    y,   z
    % u: 1,    0,   0
    % v: 0,    1,   0
    % w: 2Sx,  2Sy, |S|^2
    % So to get the slopes to calculate the bolus transport we need to divide by 2.
    % However, the Advective form is used then kwx=Sx and kwy=Sy, hence derived bolus
    % velocities will be half as great as expected.
    kgm=2;
else
    kgm=1;
end  
    for t=1:nt
        for i=1:nx
            for j=1:ny
                for k=1:nz
                    kp1 = min(k+1,nz);
                    
                    if (k<=nz);
                        maskp1 = 1;
                    else
                        maskp1 = 0;
                    end
                    
                    ubg(i,j,k,t) = (kgm/2).*dy(i,j)*(gm_psix(i,j,kp1,t)*maskp1-gm_psix(i,j,k,t))*umask(i,j,k);
                    
                    vbg(i,j,k,t) = (kgm/2).*dx(i,j)*(gm_psiy(i,j,kp1,t)*maskp1-gm_psiy(i,j,k,t))*vmask(i,j,k);
                    
                    if i~=nx
                        ip1=i+1;
                    else
                        ip1=1;
                    end
                    
                    if j~=ny
                        wbg(i,j,k,t) =(...
                             dy(ip1,j)*gm_psix(ip1,j,k,t)...
                            -dy( i ,j)*gm_psix( i ,j,k,t)...
                            +dx(i,j+1)*gm_psiy(i,j+1,k,t)...
                            -dx(i, j )*gm_psiy(i, j ,k,t)...
                            ).*cmask(i,j,k);
                    end
                end
            end
        end
    end
