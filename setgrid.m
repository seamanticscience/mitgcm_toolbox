% setgrid.m creates the model grid and stores the grid information
% various arrays. Nominally set are:
%    lon,lat spacing (dlon,dlat)
%    level thicknesses (dz)
% from which are generated:
%    tracer  grid coordinates (xc,yc,zc)
%    u-point grid coordinates (xf,yc,zc)
%    v-point grid coordinates (xc,yf,zc)
%    w-point grid coordinates (xc,yc,zf)
%    horizontal grid-lengths along cell boundaries (DXF,DYF) [2d]
%    surface area of grid-cell on sphere (zA) [2d]
%
% It also stores the above information in a MATLAB file GRID.mat
% usage:
%        load GRID.mat                   (loads everything above)
%        load GRID.mat nxc nyc xc yc     (loads just these variables/arrays)
%  
% Uses data files: n/a
% Creates data files: GRID.mat
% Creates arrays: dlon,dlat,dz,xc,yc,zc,xf,yf,zf,zA,DXF,DYF
% Uses m-scripts: grid_sphere
%
% Created  08/15/99 by adcroft@mit.edu
% Modified 11/09/99 by adcroft@mit.edu
% Maintained by adcroft@mit.edu, abiastoch@ucsd.edu

lon_lo=0;
lon_L=128;
dlon=2.8125;
lon_hi=lon_lo+(lon_L*dlon);
lat_lo=-90;
lat_L=64;
dlat=2.8125;
lat_hi=lat_lo+(lat_L*dlat);

nxc=round( (lon_hi-lon_lo)/dlon );
report([sprintf('nxc = %i factors(nxc)=',nxc) sprintf(' %i',factor(nxc)) '\n'])
nyc=round( (lat_hi-lat_lo)/dlat );
report([sprintf('nyc = %i factors(nyc)=',nyc) sprintf(' %i',factor(nyc)) '\n'])

% Convert to 1D arrays
dlon=dlon*ones([nxc 1]);
dlat=dlat*ones([nyc 1]);

[xc,yc,xf,yf,DXF,DYF,zA]=grid_sphere(dlon,dlat,lon_lo,lat_lo);

% Store data in HGRID.mat
save HGRID.mat nxc nyc xc yc xf yf dlon dlat DXF DYF zA ...
               lon_lo lon_hi lat_lo lat_hi

% Vertical grid
dz=[50 70 100 140 190 240 290 340 390 440 490 540 590 640 690];

nzc=prod(size(dz));

zf=-cumsum([0 dz]);
zc=(zf(1:end-1)+zf(2:end))/2;
report([sprintf('nzc = %i dz=',nzc) sprintf(' %4.1f',dz) '\n'])

% Store data in VGRID.mat
save VGRID.mat nzc dz zc zf

% This little ditty copies everything in the above files into GRID.mat
clear
load HGRID.mat
load VGRID.mat
save GRID.mat
