function [nc_out] = gluemnc(diags,nIter0)
% gluemnc.m
% written by david wang, ldeo
%
% purpose: use mnc_assembly.m to glue the multi-tile *.*.nc mnc output 
%	   into a single "global" nc file, which is subject to further 
%	   manipulation (e.g., adding/modifying coordiates/attributes)
%	   if necessary. 
%
%	   diags:  diagnostics name
%          nIter0: 10-digit iterate #
% EXAMPLE:
%	   foo = gluemnc('state','0000000000');

% $Header: /u/gcmpack/MITgcm/utils/matlab/gluemnc.m,v 1.2 2007/02/17 23:49:43 jmc Exp $
% $Name:  $

%if nargin ~= 2, error('there have to be two input arguments!'); end

if strcmp(diags,'grid');
    nc_in    = [diags,'.t%03d.nc'];
    nc_inone = [diags,'.t001.nc'];
    nc_out   = [diags,'.glob.nc'];
else
    if isnumeric(nIter0)
        % convert nIter0 to character string of length 10
        itertmp=num2str(nIter0);
        if length(itertmp)<10
           tmp=10-length(itertmp);
           for j=1:tmp
              itertmp=['0',itertmp];
           end
        elseif length(itertmp)>10
           error('nIter0 should have a length of 10')
        end
        iterc=itertmp;
    elseif ~exist('nIter0','var')
        [itern,iterc] = scanforfiles(diags);
        disp('scanning for possible files')
    else
        iterc=nIter0;
    end
% nc_in    = [diags,'.',nIter0,'.t%03d.nc'];
% nc_inone = [diags,'.',nIter0,'.t001.nc'];
% nc_inall = [diags,'.',nIter0,'.t*.nc'];
% nc_out   = [diags,'_',nIter0,'.nc'];
    nc_in    = [diags,'.',iterc,'.t%03d.nc'];
% nc_in    = [diags,'.',iterc,'.t*.nc'];
% diagin   = ls([diags,'.',iterc,'.*']); 
% namelength=length(diags)+19;
% nc_inone = diagin(1,1:namelength);
    nc_inone = [diags,'.',iterc,'.t001.nc'];
    nc_out   = [diags,'.',iterc,'.glob.nc'];
end

% varlist = ncload(nc_inone); % NetCDF toolbox way
tmp=nc_info(nc_inone);        % SNCTOOLS way
varlist={tmp.Dataset.Name};
nvars =   length(varlist);
vars = struct([]);

for i = 1:nvars,
  vars(i).name = char(varlist(i));
end

if strncmp(diags,'pickup',6);
    [nt,nf,exit] = pickup_assembly(nc_in, vars, nc_out);
else
    [nt,nf,exit] = mnc_assembly(nc_in, vars, nc_out);
end

% Attributes handled directly in mnc_assembly
% % reply = input('delete the original tiled files? [y/n] ','s');
% % if isempty(reply), reply = 'y'; end
% % if strcmpi(reply,'y'), delete(nc_inall); end
% if exit == 0; %1=exit with error
%     tmp=rdmnc(nc_inone);
%     attributes=tmp.attributes.global;
%     clear tmp
%     nc=netcdf(nc_out,'write');
%     nc.the_run_name=attributes.the_run_name;
%     nc.MITgcm_version=attributes.MITgcm_version;
%     nc.build_user=attributes.build_user;
%     nc.build_host=attributes.build_host;
%     nc.build_date=attributes.build_date;
%     nc.MITgcm_tag_id=attributes.MITgcm_tag_id;
%     nc.MITgcm_mnc_ver=attributes.MITgcm_mnc_ver;
%     nc.tile_number='global';
%     nc.sNx=attributes.sNx;
%     nc.sNy=attributes.sNy;
%     nc.OLx=attributes.OLx;
%     nc.OLy=attributes.OLy;
%     nc.nSx=attributes.nSx;
%     nc.nSy=attributes.nSy;
%     nc.nPx=attributes.nPx;
%     nc.nPy=attributes.nPy;
%     nc.Nx=attributes.Nx;
%     nc.Ny=attributes.Ny;
%     nc.Nr=attributes.Nr;
%     close(nc);
% end
return


function [itern,iterc] = scanforfiles(diags)

itern=[];
iterc=[];

% for k=1:100 % iteratively find first tile number present
%        exp=sprintf('.*.t%03d.nc',k);
      exp=[diags,'.*.t001.nc'];
%       allfiles=dir(exp) % gets all file details in structural array
      allfiles=ls(exp); % gets all file details in structural array
%       if ~isempty(allfiles); break; end
%end
if isempty(allfiles)
    exp=[diags,'.*.t*.nc'];
    allfiles=ls(exp);
    namelength=length(diags)+19;
    allfiles=allfiles(1,1:namelength);
    if isempty(allfiles)
        error('error finding any tile files')
    end
end
for k=1:size(allfiles,1);
%    hh=allfiles(k).name;
%    itern(k)=str2double( hh(end-17:end-8) );
    itern(k)=str2double(allfiles(1,end-17:end-8));
    itertmp=num2str(itern);
    tmp=10-length(itertmp);
    for j=1:tmp
        itertmp=['0',itertmp];
    end
    iterc=itertmp;
end