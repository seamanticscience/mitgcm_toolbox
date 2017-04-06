function [psi, psimask] = mit_overturning(v,vmask,dx,dz,addlayer)
%function [psi, psimask] = mit_overturning(v,vmask,dx,dz);
% overturning stream function for mitgcm model, time slabs are handled as 
% cell objects, integration from bottom to top.

% $Header: /u/gcmpack/MITgcm/verification/tutorial_global_oce_latlon/diags_matlab/mit_overturning.m,v 1.3 2006/08/12 20:25:13 jmc Exp $
% $Name:  $

  if nargin < 5
    addlayer = 0;
  elseif (nargin > 5 && nargin < 4)
    error('needs 4 or 5 arguments')
  end
  [nx,ny,nz] = size(vmask);
  for kz=1:nz
    dxdzs(:,:,kz) = (vmask(:,:,kz).*dx)*dz(kz);
  end
  % mask for stream function
  pmask = change(squeeze(nanmean(vmask)),'~=',NaN,1);
  % add another psi-point at the bottom of each column
  for ky=1:ny;
    iz = min(find(isnan(pmask(ky,:))));
    if ~isempty(iz) && iz > 1
      pmask(ky,iz) = 1;
    end
  end

  nt = size(v,4);
  for k=1:nt
    vdxdz = squeeze(v(:,:,:,k)).*dxdzs; %(dxdzs.*vmask);
    zave  = fliplr(squeeze(sum(change(vdxdz,'==',NaN,0),1))); %sum in the x-direction
    psi(:,:,k)   = -fliplr(cumsum(zave,2)).*pmask;
  end
    % add another layer at the bottom; psi is zero in the layer by definition
  if addlayer
    %disp('mit_overturning: adding a layer with psi = 0 at the bottom')
    pmask = [pmask pmask(:,end)];
    psi = cat(2,psi,change(psi(:,end,:),'~=',NaN,0));
  end
  if nargout == 2
    psimask = pmask;
  end
  return




  
