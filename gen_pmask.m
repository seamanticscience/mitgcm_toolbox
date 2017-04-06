% This is a script that executes the commands for
% creating a "PMASK" file (notation of MITgcm Classic)
% which is a 3D field containing 1 for ocean and 0 for
% land. This mask is *NOT* necessarily the mask used
% in the model but is representative of the largest
% possible mask that the model might use (cf. hFacMin).
% We need it here in order to guarantee that there
% is valid data in the climatological input data files.
%
% Uses data files: bathymetry.bin HGRID.mat VGRID.mat
% Creates data files: pmask.bin
% Creates arrays: n/a
% Uses m-scripts: hfac
%
% Created  11/11/99 by adcroft@mit.edu
% Modified          by
% Maintained by adcroft@mit.edu, abiastoch@ucsd.edu

clear all
format compact
more off

load FMT.mat
load HGRID.mat nxc nyc
load VGRID.mat dz

hn=rdda('bathy_le.bin',[nxc nyc],1,fmt,Ieee);

% Generate mask and store in pmask.bin
[msk] = hfac(dz,hn,0,0);
msk(find(msk~=0))=1;
wrda('pmask.bin',msk,1,fmt,Ieee);
