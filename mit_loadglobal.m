% m-file: mit_loadglobal.m
% plot various diagnostics of the global 4x4 degree runs

% uses
% mit_loadgrid 
% mit_getparm
% mit_oceanmasks
% mit_timesteps
% mit_loadglobaldata
% mit_globalmovie
% mit_plotstreamfunctions
% mit_plotzonalvelocity
% mit_globalmean
% mit_zonalmean
% mit_meridflux
% mit_plotmeandrift

% $Header: /u/gcmpack/MITgcm/verification/tutorial_global_oce_latlon/diags_matlab/mit_loadglobal.m,v 1.3 2006/08/12 20:25:13 jmc Exp $
% $Name:  $

% read in all grid files, etc. This has to be done at the very beginning!
grd = mit_loadgrid('.');
dname = grd.dname;
% add some masks
grd = mit_oceanmasks(grd);

tavesteps=mit_timesteps('tave');
dlmwrite('times',tavesteps.tim','newline','unix')

% for diagnostics, at the moment they are all the same output frequency, so
% just read the first one in the data.diagnostics file
diagsteps=mit_timesteps(mit_getparm('data.diagnostics','filename'));
dlmwrite('times_diag',diagsteps.tim','newline','unix')

meanfields=1;
disp('reading time-averaged fields')

% Extract attributes and axis for netcdf creation later....
if ~exist('attributes','var')
    tave = rdmnc('tave.*.glob.nc','Xp1','Yp1');
    Xp1=tave.Xp1;
    Yp1=tave.Yp1;
    attributes=tave.attributes.global;
    clear tave state
end
%% Processing
% Calulate transports
if ~exist('psi_data.nc','file')
    disp('Calculating streamfunction data...')
    tic
    [filename]=mit_calc_psi(grd,tavesteps);
    toc
else
    disp('psi_data.nc file exists, skipping')
end

% [global_eul,global_psi,global_eddys,...
%     atlantic_eul,atlantic_psi,atlantic_eddys,...
%     pacific_eul,pacific_psi,pacific_eddys,...
%     baro_eul,baro_eddys,baro_psi]=mit_plotstreamfunctions(grd,tavesteps);

if exist('psi_data.nc','file')
    disp('Plotting MITgcm streamfunction data...')
    tic
    mit_plotstreamfunctions(grd,tavesteps);
    toc
else
    disp('psi_data.nc does not exist, skip plotting')
end

%% Do some biogeochemistry calculations
if ~exist('ph_data.glob.nc','file')
   disp('Calculating global pH...')
   tic
   [filename]=calc_pco2(grd,tavesteps);
   toc
   % optional outputs go in this order
   %phnew,pco2,co3,co2s,ff,ffp,k0,k1,k2
   clear functions % frees up open netcdf file handles
else
    disp('ph_data.glob.nc exists, skipping')
end

if ~exist('dic_sat.glob.nc','file')
   disp('Calculating DIC Components...')
   tic
   [pstar]=calc_sat_dic(grd,tavesteps)
   toc
   clear functions % frees up open netcdf file handles
else
   disp('dic_sat.glob.nc exists, skipping')
end

 if ~isempty(mit_getparm('data.diagnostics','CSAT')) && ...
         exist(strrep(diagsteps.filearr(2:end-1),'surf','component'),'file')
    plot_sat_dic(grd);
 end

%% Plot some other variables
% plot zonal velocity sections
[u]=mit_plotzonalvelocity(grd,tavesteps);

% plot mean t and s sections
mit_plotmeansections(grd,tavesteps);

% plot surface fields
mit_plotsurfacefields(grd,tavesteps,150)

% Plot biogeochemical fields
if exist('data.gchem','file') && strcmpi('true',(mit_getparm('data.gchem','useDIC')))  
   mit_plotbiogeochem(grd,tavesteps);
end

clear mex
%clear some space
clear ubg vbg *_psi *_eul

if length(tavesteps.tim)>1;
    disp('Plotting mean model drift timeseries...')
    tic
    mit_plotmeandrift(grd,tavesteps,diagsteps);
    toc
end

