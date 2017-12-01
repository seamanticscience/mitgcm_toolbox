function [ubg,vbg,wbg]=mit_gmvel(kwx,kwy,nx,ny,nz,nt,dx,dy,dz,cmask,umask,vmask,area,varargin)
% function [ubg,vbg,zbg]=mit_gmvel(kwx,kwy,nx,ny,nz,nt,dx,dy,dz,cmask,umask,vmask,area,varargin)
%
% Caluculate MITgcm Bolus Velocities for the GM-eddy mixing
% parameterisation. Kwx,Kwy and Kwz are Kgm*Isopycnal_Slope
%
% original supplied by stephd (03/10/08) modified by JML for NetCDF (14/10/08)
% Corrected by JML (16/12/13) to actually calculate vertical bolus
%       component rather than just nans
% Modified (2/3/17) to better match MITgcm code by JML (see gmredi_residual_flow.F).

% Original code.
% kgm=1; 
% str2=['GM_Kwy-T']; 
% fid=fopen([run1str,'/',str2,'.',timestr,'.data'],'r',enstr); 
% clear field3, field3=fread(fid,'float32'); 
% field3=reshape(field3,nx,ny,nz); 
% % 
% % calculate bolus velocity at center of grid 
% psiby=field3*kgm/2; 
% for k=1:nz-1, 
%     vbc(:,:,k)=-(psiby(:,:,k)-psiby(:,:,k+1))/dz(k); 
% end % for k 
% 
% vbc(:,:,nz)=(psiby(:,:,nz)-0)/dz(nz); 
% % put back on v location 
% for i=1:nx, 
%    for j=1:ny-1, 
%      vbg(i,j,:)=(vbc(i,j,:)+vbc(i,j+1,:))/2; 
%    end, % for j 
%     vbg(i,ny,:)=(vbc(i,ny,:)+vbc(i,1,:))/2; 
%  end % for i}

% Setting of KGM. Kwx and Kwy are the slope*thickness diffusivity so why are we multiplying by a constant?
% This is necessary because the gmredi tensor when the skew flux is used and Kgm=Krho is:
%    x,    y,   z
% u: 1,    0,   0
% v: 0,    1,   0
% w: 2Sx,  2Sy, |S|^2
% So to get the slopes to calculate the bolus transport we need to divide by 2.
% However, the Advective form is used then kwx=Sx and kwy=Sy, hence derived bolus
% velocities will be half as great as expected.
if ~isempty(varargin)
    if varargin{1} % user indicates GM_ADVFORM is true
        psiby=kwy;
        psibx=kwx;
    else
        psiby=kwy./2;
        psibx=kwx./2;
    end
else
    if strcmpi(mit_getparm('data.gmredi','GM_AdvForm'),'true')
        psiby=kwy;
        psibx=kwx;
    else
        psiby=kwy./2;
        psibx=kwx./2;
    end
end

% calculate (v/y) bolus velocity  
psiby(isnan(psiby))=0;
vbg=zeros(nx,ny,nz,nt);

for t=1:nt
    for k=1:nz-1
        vbg(:,:,k,t)=(psiby(:,:,k+1,t)-psiby(:,:,k,t))/dz(k);
    end
    vbg(:,:,nz,t)=(0-psiby(:,:,nz,t))/dz(nz);
end
vbg=vbg.*repmat(vmask,[1,1,1,nt]);

% % put back on v location 
% vbg=zeros(nx,ny,nz);
% for i=1:nx, 
%    for j=1:ny-1, 
%      vbg(i,j,:)=(vbc(i,j,:)+vbc(i,j+1,:))/2; 
%    end
%    % This is not strictly good, but there is no north/south pole so is zero
%    % anyway
%    vbg(i,ny,:)=(vbc(i,ny,:)+vbc(i,1,:))/2; 
% end
 
% calculate (u/x) bolus velocity 
psibx(isnan(psibx))=0;
ubg=zeros(nx,ny,nz,nt);

for t=1:nt
for k=1:nz-1, 
    ubg(:,:,k,t)=(psibx(:,:,k+1,t)-psibx(:,:,k,t))/dz(k); 
end 
ubg(:,:,nz,t)=(0-psibx(:,:,nz,t))/dz(nz); 
end
ubg=ubg.*repmat(umask,[1,1,1,nt]);

% % put back on u location 
% ubg=zeros(nx,ny,nz);
% for j=1:ny, 
%    for i=1:nx-1, 
%      ubg(i,j,:)=(ubc(i,j,:)+ubc(i+1,j,:))/2; 
%    end, % for i 
%      ubg(nx,j,:)=(ubc(nx,j,:)+ubc(1,j,:))/2; % This is ok because axis is modulo
%  end % for j}

%Calculate d(psiby*dx)
psiby=psiby.*repmat(dx,[1,1,nz,nt]);
vbcy=zeros(nx,ny,nz,nt);
for t=1:nt
for i=1:nx
    for j=1:ny
        for k=1:nz
            if j~=ny
                vbcy(i,j,k,t)=(psiby(i,j+1,k,t)-psiby(i,j,k,t));
            else
                vbcy(i,j,k,t)=(0-psiby(i,j,k,t));
            end
        end
    end
end
end

% Calculate d(psibx*dy)
psibx=psibx.*repmat(dy,[1,1,nz,nt]);
ubcx=zeros(nx,ny,nz,nt);
for t=1:nt
for i=1:nx
    for j=1:ny
        for k=1:nz
            if i~=nx
                ubcx(i,j,k,t)=(psibx(i+1,j,k,t)-psibx(i,j,k,t));
            else
                ubcx(i,j,k,t)=(psibx(1,j,k,t)-psibx(i,j,k,t)); % modulo axis 
            end
        end
    end
end
end

wbg=(vbcy+ubcx)./repmat(area,[1,1,nz,nt]); 
wbg=wbg.*repmat(cmask,[1,1,1,nt]);

%Calculate psiby/dy
% vbcy=zeros(nx,ny,nz);
% for i=1:nx
%     for j=1:ny
%         for k=1:nz
%             if j~=ny
%                 vbcy(i,j,k)=(psiby(i,j+1,k)-psiby(i,j,k))./dy(i,j);
%             else
%                 vbcy(i,j,k)=(0-psiby(i,j,k))./dy(i,j);
%             end
%         end
%     end
% end
% 
% 
% % Calculate psibx/dx
% ubcx=zeros(nx,ny,nz);
% for i=1:nx
%     for j=1:ny
%         for k=1:nz
%             if i~=nx
%                 ubcx(i,j,k)=(psibx(i+1,j,k)-psibx(i,j,k))./dx(i,j);
%             else
%                 ubcx(i,j,k)=(psibx(1,j,k)-psibx(i,j,k))./dx(i,j); % modulo axis 
%             end
%         end
%     end
% end
% 
% wbg=(vbcy+ubcx).*cmask;

% % put back on w location 
% wbg=zeros(nx,ny,nz);
% for j=1:ny, 
%    for k=nz:-1:2, 
%       wbg(:,j,k)=(wbc(:,j,k)+wbc(:,j,k-1))/2; 
%    end, % for k 
%       wbg(:,j,1)=(wbc(:,j,1)+0)/2; % No vertical flow through the surface @z=0,w=0
% end % for j}
