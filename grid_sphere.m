function [xc,yc,xf,yf,DXF,DYF,zA]=grid_sphere(dlon,dlat,varargin)
% Calculate surface area of grid-cell on fixed spherical grid
%
% [xc,yc,xf,yf,DXF,DYF,zA]=grid_sphere(dlon,dlat,lon0,lat0);

% Constants
Aearth=6370e3;

if nargin <=2
 lon0=0;
else
 lon0=varargin{1};
end
if nargin <= 3
 lat0=-90;
else
 lat0=varargin{2};
end
if prod(size(dlon)) == 1
 nx=round( 360/dlon );
 ddlon=dlon*ones([nx 1]);
else
 nx=prod(size(dlon));
 ddlon=reshape( dlon ,[nx 1]);
end
if prod(size(dlat)) == 1
 ny=round( (90-lat0)/dlat );
 ddlat=dlon*ones([1 ny]);
else
 ny=prod(size(dlat));
 ddlat=reshape( dlat ,[1 ny]);
end

% Coordinates
xf=cumsum([lon0 ddlon']');
yf=cumsum([lat0 ddlat]);
xc=(xf(1:end-1)+xf(2:end))/2;
yc=(yf(1:end-1)+yf(2:end))/2;

% Convert to radians
ddlon=ddlon*pi/180;
ddlat=ddlat*pi/180;

% Physical lengths only grid-cell edges
[DXF,DYF]=ndgrid(Aearth*ddlon,Aearth*ddlat);
DXF=Aearth*(ddlon*cos(pi/180*yc));

% Surface area of grid-cell
zA=Aearth^2*( ddlon * ( sin( pi/180*yf(2:end) )-sin( pi/180*yf(1:end-1) ) ) );

xf=xf(1:end-1);
yf=yf(1:end-1);
