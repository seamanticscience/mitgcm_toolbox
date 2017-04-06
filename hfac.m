function [hfacC,ddz] = hfac(dz,H,hfacmin,dzmin)
%
% Create dz(i,j,k) and hfacC(i,j,k) from dz(k), H(i,j)
%
% e.g.
%   [hfacC,ddz] = hfac(dz,H,hfacmin,dzmin)
% e.g. To create a "p-mask" from above
%   pmask=zeros(size(hfacC)); pmask( find(hfacC>0) )=1;

[nx,ny]=size(H);
nz=prod(size(dz));
N=[nx ny nz];

zf=[0 -cumsum(dz)];

report('Level k = xxx')
for k=1:nz,
 report('\b\b\b%3i',k)
 ddd=zf(k)-H;
 ddd(find(ddd < 0)) = 0;
 ddd(find(ddd > dz(k))) = dz(k);
 ddd(find(ddd < hfacmin*dz(k)/2 & ddd ~= 0)) = 0;
 ddd(find(ddd >= hfacmin*dz(k)/2 & ddd < hfacmin*dz(k))) = hfacmin*dz(k);
 ddd(find(ddd < dzmin/2 & ddd ~= 0)) = 0;
 ddd(find(ddd >= dzmin/2 & ddd < dzmin)) = dzmin;
 ddz(:,:,k)=ddd;
 hfacC(:,:,k)=ddd/dz(k);
end
report('\n')
