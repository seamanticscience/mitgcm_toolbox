function [nt,nf,exit] = mnc_assembly(fpat,vars, fout, debugMode)
% Function [nt,nf] = mnc_assembly(fpat,vars, fout)
%
% INPUTS
% fpat   string containing the file pattern
% vars   structure array of variable names
%
% fout   output file pattern (DEF: "all.%05d.nc")
% 'debug'  as an additional text argument will enable extra info to be
% printed to screen
%
% OUTPUTS
% nt     number of usable tiles found
% nf     number of output files written
% exit   success of assembly - 0 is completed, 1 is exited with error
%
% This function "assembles" MNC output.  It finds all the per-tile
% NetCDF files that match the input pattern, does some basic "sanity"
% tests to determine whether the files have compatible sizes, and
% then assembles all of the requested data (all of the variables)
% into one or more "global" NetCDF files.  The global files have
% the following dimension conventions:
%
%   "exch 1": all values are within a global horizontal grid
%             and indicies are (X,Y,Z,T)
%
%   "exch 2": all values are within one of up to six "faces" 
%             of a global cube with indicies (Xf,Yf,F,Z,T)
%
% where "X,Y.Z,T" are global space/time indicies, "Xf,Yf" are local
% per-face spatial indicies, and "F" is a face index.
%
% An example of how to use this script is:
% 
%   vars = struct([]);
%   vars(1).name = 'iter';
%   vars(2).name = 'U';
%   vars(3).name = 'Unk';
%   vars(4).name = 'V';
%   vars(5).name = 'Temp';
%   vars(6).name = 'S';
%   fpat = 'exp0_20041126_0001/state.0000.%06d.nc';
%   [nt,nf] = mnc_assembly(fpat,vars);
%
% and the resutlt is written as "all.00000.nc"

% $Header: /u/gcmpack/MITgcm/utils/matlab/mnc_assembly.m,v 1.6 2007/02/17 23:49:43 jmc Exp $
% $Name:  $
%
% Modified to use SNCTOOLS for output, Dec 2011, Jon Lauderdale + John Evans
%                          for input, Dec 2012, Jon Lauderdale

%=====  Argument checking and defaults  =====

%% Set up vars to run as script for debugging purposes
% iterc='0014400000';
% diags='tave';
% nc_in    = [diags,'.',iterc,'.t%03d.nc'];
% nc_inone = [diags,'.',iterc,'.t001.nc'];
% nc_out   = [diags,'.',iterc,'.glob.nc'];
% tmp=nc_info(nc_inone);        % SNCTOOLS way
% varlist={tmp.Dataset.Name};
% nvars =   length(varlist);
% vars = struct([]);
% 
% for i = 1:nvars,
%   vars(i).name = char(varlist(i));
% end
% 
% fpat=nc_in;
% fout=nc_out;
% 
% clear tmp
%%
% if nargin < 2
%   disp('Error: there must be at least 2 arguments!');
%   return
% end
% 
% if nargin < 3
%   fout = 'all.%05d.nc';
% end
% if nargin < 4
%  fsize = 2.0e+9;
% end


%% =====  Find and open all the matching files  =====


if nargin==2
    debugMode=false;
end

nt      = 0;
nf      = 0;
exit    = 0;
all_ncf = struct([]);

%  Find all of the files
exch2_msg = 0;
tmax = 200;
frdone = 0;
it = 0;
while frdone == 0

  it = it + 1;
  fnm = sprintf(fpat,it);
%  disp(fnm);

  %  Check that the file exists
  if ~exist(fnm,'file') || it >= tmax
%   fid = fopen(fnm, 'r');
%   if fid < 0
%     if it >= tmax
      frdone = 1;
%    end
    continue;
  end
  
% Old NetCDF toolbox way    
%   %  Open the NetCDF file
%   fnc = netcdf(fnm, 'nowrite');
%   if isempty(fnc)
%     continue;
%   end
  fnm_tile_info=nc_info(fnm);
    
%  Check for exch1/exch2 grid
  exch = 1;
%  exch2_myFace = fnc.exch2_myFace(:);
  exch2_test=intersect({fnm_tile_info.Attribute.Name},'exch2_myFace');
  if ~isempty(exch2_test);
%  if ~isempty(exch2_myFace)
    exch = 2;
    if exch2_msg == 0
      exch2_msg = 1;
      disp('  Grid type appears to be: "exch2"');
    end
  end
  
  % use SNC's nc_info function
  n = length(all_ncf) + 1;
  all_ncf(n).name           = fnm;   % No change
  all_ncf(n).nc             = {fnm_tile_info}; %???
  all_ncf(n).exch           = exch;  % No change
  all_ncf(n).tile_number    = fnm_tile_info.Attribute(find(strcmp({fnm_tile_info.Attribute.Name},'tile_number'))).Value;
  all_ncf(n).bi             = fnm_tile_info.Attribute(find(strcmp({fnm_tile_info.Attribute.Name},'bi'))).Value;
  all_ncf(n).bj             = fnm_tile_info.Attribute(find(strcmp({fnm_tile_info.Attribute.Name},'bj'))).Value;
  all_ncf(n).sNx            = fnm_tile_info.Attribute(find(strcmp({fnm_tile_info.Attribute.Name},'sNx'))).Value;
  all_ncf(n).sNy            = fnm_tile_info.Attribute(find(strcmp({fnm_tile_info.Attribute.Name},'sNy'))).Value;
  all_ncf(n).Nx             = fnm_tile_info.Attribute(find(strcmp({fnm_tile_info.Attribute.Name},'Nx'))).Value;
  all_ncf(n).Ny             = fnm_tile_info.Attribute(find(strcmp({fnm_tile_info.Attribute.Name},'Ny'))).Value;
  all_ncf(n).Z              = fnm_tile_info.Attribute(find(strcmp({fnm_tile_info.Attribute.Name},'Nr'))).Value;

  if exch == 2
    all_ncf(n).exch2_myFace = fnm_tile_info.Attribute(find(strcmp({fnm_tile_info.Attribute.Name},'exch2_myFace'))).Value;
    all_ncf(n).exch2_tbasex = fnm_tile_info.Attribute(find(strcmp({fnm_tile_info.Attribute.Name},'exch2_tbasex'))).Value;
    all_ncf(n).exch2_tbasey = fnm_tile_info.Attribute(find(strcmp({fnm_tile_info.Attribute.Name},'exch2_tbasey'))).Value;
  end

  clear fnc
end


%% =====  Do some basic sanity checks =====

%  check for number of files/tiles found
if isempty(all_ncf)
  disp('Error: no tiles found--no need to do any assembly!');
  exit=1;
  return
elseif length(all_ncf) == 1
  disp('Error: one tile found--no need to do any assembly!');
  exit=1;
  return
else
  fprintf('  Found %d files matching the pattern:  "%s"\n', ...
               length(all_ncf), fpat);
end

%  check for consistent "exch" version
if prod(double([all_ncf.exch] == all_ncf(1).exch)) ~= 1
  disp('Error: not all the "exch" types of the files match.');
  exit=1;
  return;
end

%  check for consistent sNx,sNy
if (prod(double([all_ncf.sNx] == all_ncf(1).sNx)) ~= 1) ...
      || (prod(double([all_ncf.sNy] == all_ncf(1).sNy)) ~= 1)
  disp('Error: the "sNx,sNy" values for all the tiles are not');
  disp('   uniform.  Future versions of this function will be');
  disp('   able to handle non-uniform grid sizes but this');
  disp('   feature is not yet implemented.');
  exit=1;
  return;
end

%  check for redundant tiles and "time series" output
if length(all_ncf) ~= length(unique([all_ncf.tile_number]))
  disp('length(all_ncf)')
  length(all_ncf)
  disp('length(unique([all_ncf.tile_number]))')
  length(unique([all_ncf.tile_number]))
  all_ncf(:).tile_number
  
  disp('Error: redundant tiles were found.  Please check that');
  disp('   the file pattern does not specify output spanning');
  disp('   multiple model runs or even multiple time series');
  disp('   within a single model run.  For multi-time-series');
  disp('   data sets, EACH "LEVEL" IN THE OUTPUT SERIES MUST');
  disp('   BE ASSEMBLED SEPARATERLY.');
  exit=1;
  return
end


%% =====  Get the dims/vars associations  =====

mydims = struct('names', {}, 'lens', {});
myvars = struct([]);
clear tncf;
for ivar = 1:length(vars)
  mydim_names = {};
  mydim_sizes = {};
  myatt.names = {};
  myatt.types = {};
  myatt.data  = {};
  myname = vars(ivar).name;
 
  % Look for variable in first tile and then check it is the same in
  % subsequent tiles
  if ~license('test','Distrib_Computing_Toolbox') %|| isempty(gcp('nocreate'));
      disp(['  Looking for variable:   ' myname]);
  end
  
  itile = 1;
  fnm_tile_info = nc_info(all_ncf(itile).name);
  
  % test if the variable exists in first tile, if not, try next name
  % Much more simple using SNCTOOLS and nc_info
  varnamei=find(strcmp({fnm_tile_info.Dataset.Name},myname));
  if isempty(varnamei)
  %if isempty(ncv) 
    warns = ['    Warning: variable "%s" is not defined in "%s"\n' ...
             '      so it will be ignored.\n'];
    fprintf(warns,myname,all_ncf(itile).name);
    continue
  end
  
%   Much more simple using SNCTOOLS and nc_info
   mytype = fnm_tile_info.Dataset(varnamei).Datatype;
   mydim_names=fnm_tile_info.Dataset(varnamei).Dimension;
   mydim_sizes=num2cell(fnm_tile_info.Dataset(varnamei).Size);

%   Much more simple using SNCTOOLS and nc_info
  if ~isempty(fnm_tile_info.Dataset(varnamei).Attribute)
      myatt.names = {fnm_tile_info.Dataset(varnamei).Attribute.Name};
      myatt.types = {fnm_tile_info.Dataset(varnamei).Attribute.Datatype};
      myatt.data = {fnm_tile_info.Dataset(varnamei).Attribute.Value};
  else
      myatt.names={};
      myatt.types={};
      myatt.data={};
  end
  
  %  confirm: vars have same dim names across all files
  ierr = 0;
 
% Using SNCTOOLS and nc_info
  for itile = 2:length(all_ncf)
      fnm_tile_info = nc_info(all_ncf(itile).name);
      varnamei=find(strcmp({fnm_tile_info.Dataset.Name},myname));
      if isempty(varnamei)
          warns = ['    Warning: variable "%s" is not defined in "%s"\n' ...
              '      so it will be ignored.\n'];
          fprintf(warns,myname,all_ncf(itile).name);
          continue
      end
      
      tmpdim_names=fnm_tile_info.Dataset(varnamei).Dimension;
      tmpdim_sizes=num2cell(fnm_tile_info.Dataset(varnamei).Size);
      
      if min(strcmp(mydim_names,tmpdim_names))==0 % Check that the dimensions have same name
          warns = ...
              ['    Warning: variable "%s" is not CONSISTENTLY defined.\n' ...
              '      It has different dimensions in different files so\n' ...
              '      so it will be ignored.\n'];
          fprintf(warns,myname);
          ierr = 1;
          break
      end
      % Not sure why this is necessary as all tiles have same dimensions....?
      %     mydim_sizes{inm} = max([ length(tmpdims{inm}) mydim_sizes{inm} ]);
  end
   
  if ierr == 0 
    %  check: does the variable have a "horizontal" component
    has_horiz   = 0;
    horiz_names = { 'X' 'Y' 'Xp1' 'Yp1' };
    for id = 1:length(mydim_names)
      if ~isempty(intersect(horiz_names,mydim_names{id})) % max(strcmp(horiz_names,mydim_names{id}))
        has_horiz = 1;
      end
    end
    if debugMode
        disp([ ' check for horizontal componant    ' myname ' '...
            sprintf('%d',has_horiz) ]);
    end
    
    imy = length(myvars) + 1;
    myvars(imy).name       = myname;
    myvars(imy).type       = mytype;
    myvars(imy).dim_names  = mydim_names;
    myvars(imy).dim_sizes  = mydim_sizes;
    myvars(imy).atts       = myatt;
    myvars(imy).has_horiz  = has_horiz;
    
    % this is necessary to make it work with Matlab 6.5
    if isempty([mydims.names])
      addl = mydim_names;
    else
      addl = setdiff(mydim_names,[mydims.names]);
    end
    for iaddl = 1:length(addl)
      np1 = length(mydims) + 1;
      mydims(np1).names = addl(iaddl);
      mydims(np1).lens  = mydim_sizes(find(strcmp(addl(iaddl),mydim_names)));
    end
    
    % Can just extract this direct from the file?
    % tmp.names=cellstr({fnm_tile_info.Dimension.Name});
    % tmp.lens={fnm_tile_info.Dimension.Length};    
  end
end

%  For exch == 2, we need to add a "face" dimension
if all_ncf(1).exch == 2
  np1 = length(mydims) + 1;
  mydims(np1).names = { 'iface' };
  mydims(np1).lens  = { length(unique([all_ncf.exch2_myFace])) };
end

% myvars.name
% myvars.dim_names
% myvars.dim_sizes
% myvars(2).dim_names
% myvars(2).dim_names(4)

% mydims
% length(mydims)
% [ mydims.names ]
% [ mydims.lens ]


%% =====  Assemble!  =====


if all_ncf(1).exch == 1

  %  exch "1":
  
% $$$   bi_max = max([all_ncf.bi]);
% $$$   bj_max = max([all_ncf.bj]);
% $$$   Xmax = bi_max * all_ncf(1).sNx;
% $$$   Ymax = bj_max * all_ncf(1).sNy;
  Xmax = all_ncf(1).Nx;
  Ymax = all_ncf(1).Ny;
  % at this point I have to make some assumptions about the domain
  % decomposition 
  bi_max = Xmax/all_ncf(1).sNx;
  bj_max = Ymax/all_ncf(1).sNy;
  tx=Xmax/bi_max; % tx is length of tile x-axis
  ty=Ymax/bj_max; % ty is length of tile y-axis

  itile = 0;
  for bj=1:bj_max
    for bi=1:bi_max
      itile = itile+1;
      all_ncf(itile).bi=bi;
      all_ncf(itile).bj=bj;
    end
  end

  horzdim = struct('names',{},'lens',{});
  horzdim(1).names = { 'X' };  horzdim(1).lens = { Xmax     };
  horzdim(3).names = { 'Y' };  horzdim(3).lens = { Ymax     };
  horzdim(4).names = {'Yp1'};  horzdim(4).lens = { Ymax     }; % Purposely truncate
  horzdim(2).names = {'Xp1'};  horzdim(2).lens = { Xmax     }; % Purposely truncate
%   horzdim(4).names = {'Yp1'};  horzdim(4).lens = { Ymax + 1 }; % Match output from unix routines
%   horzdim(2).names = {'Xp1'};  horzdim(2).lens = { Xmax + 1 }; % Match output from unix routines
  horzdim(5).names = { 'T' };  horzdim(5).lens = { 0        };
  
  iseq = 0;
  foutnm = sprintf(fout, iseq);
  nc_create_empty(foutnm,'clobber');
% add some space to the top of the file to allow attributes to be added
% quickly to gigabyte class files (~20Mb)
%  nc_padheader(foutnm,20000);

% Add Global Attributes to global data file
disp(['Writing global attributes to ',foutnm])

% Add Global Attributes to global data file
% Have to interact with mexnc here, but not sure if this is backwards
% compatible, besides, we already have global attributes in fnm_tile_info
%
% [ncid,status]=mexnc('open',all_ncf(1).name,'nowrite'); if status, error ( mexnc ( 'strerror', status ) ), end
% [ngatts, status] = mexnc('INQ_NATTS',ncid);            if status, error ( mexnc ( 'strerror', status ) ), end
% 
% atts=struct;
% for attnum = 1:ngatts
%     atts(attnum).attname = mexnc('inq_attname',ncid, nc_global, attnum-1);
%     atts(attnum).attval = nc_attget(all_ncf(1).name, nc_global,atts(attnum).attname);
%     if strcmp(atts(attnum).attname,'tile_number'); atts(attnum).attval='global'; end % Edit tilenum to be 'global'
% end
% mexnc('close',ncid);

% Edit tilenum to be 'global' and a character datatype
fnm_tile_info.Attribute(find(strcmp({fnm_tile_info.Attribute.Name},'tile_number'))).Value='global'; 
fnm_tile_info.Attribute(find(strcmp({fnm_tile_info.Attribute.Name},'tile_number'))).Nctype=2;
fnm_tile_info.Attribute(find(strcmp({fnm_tile_info.Attribute.Name},'tile_number'))).Datatype='char';

for attnum = 1:length(fnm_tile_info.Attribute)
    if strcmp(fnm_tile_info.Attribute(attnum).Datatype,'char') % write out character attributes
        comm=sprintf('nc_attput(foutnm,nc_global,''%s'',''%s'');',...
            fnm_tile_info.Attribute(attnum).Name,...
            num2str(fnm_tile_info.Attribute(attnum).Value));
        eval(comm)
    else % write out numerical attributes
        comm=sprintf('nc_attput(foutnm,nc_global,''%s'',%s);',...
            fnm_tile_info.Attribute(attnum).Name,...
            num2str(fnm_tile_info.Attribute(attnum).Value));
        eval(comm)
    end
    % disp(['Adding global attributes: ', comm])
end

%     Define Dimension name/size, e.g.
%     idim-loop comm:  fonc('T') = 0;
%     idim-loop comm:  fonc('X') = 128;
%     idim-loop comm:  fonc('Y') = 64;
%     idim-loop comm:  fonc('Z') = 15;
%     idim-loop comm:  fonc('Xp1') = 128; % Match output from unix routines
%     idim-loop comm:  fonc('Yp1') = 64; % Match output from unix routines
%     idim-loop comm:  fonc('Zl') = 15;
  for idim = 1:length(mydims)
    dname = mydims(idim).names{1};
    ind = find(strcmp(dname,[horzdim.names]));
    if ~isempty(ind)
      dlen = horzdim(ind).lens{1};
    else
      dlen = mydims(idim).lens{1};
    end
    %disp(['dlen = ', dlen])
%    comm = sprintf('fonc(''%s'') = %d;',dname,dlen);
    comm = sprintf('nc_adddim(foutnm,''%s'',%d);',dname,dlen);
    if debugMode
        disp([ 'Define dimensions:  ' comm ]);
    end
    eval(comm);
  end

%     Define Variable dimensions/name/units, e.g.
  for ivar = 1:length(myvars)
    % For each variable, defines it as a type (float/double etc) and
    % appends dimensions
    comm= sprintf('nc_addvar(foutnm,struct(''Name'',''%s'',''Datatype'',''%s'',''Dimension'',{%s}));',...
        myvars(ivar).name,...
        strrep(strrep(myvars(ivar).type,'char','text'),'long','int'),...
        strrep(cell2str(myvars(ivar).dim_names),';',''));
    if debugMode
        disp(['Define variables:  ' comm ]);
    end
    eval(comm);
    for iat = 1:length(myvars(ivar).atts.names)
      % Adds and attributes to the variable definition
      comm=sprintf('nc_attput(foutnm,''%s'',''%s'',''%s'');',...
          myvars(ivar).name,...
          myvars(ivar).atts.names{iat},...
          myvars(ivar).atts.data{iat});
      if debugMode
          disp(['    Add variable attributes:  ' comm ]);
      end
      eval(comm);
    end
  end

%% End Define Mode
  %  Here is where we need to check the output file size and start
  %  another file in the sequence, if necessary.
  
  for ivar = 1:1:length(myvars)
      if ~license('test','Distrib_Computing_Toolbox') %|| isempty(gcp('nocreate'));
          fprintf('  Copying variable:   %s\n',myvars(ivar).name)
      end
      
      ncdata_global=[]; % Temporary array contains accumulated tiles before output
      
      for itile = 1:length(all_ncf)
          
          if (myvars(ivar).has_horiz == 1) || (itile == 1)
              
              clear nct;
% Load variable of particular tile using NetCDF toolbox 
% Made redundant by SNCTOOLS rework
% nct = all_ncf(itile).nc{1};
%              disp(['bi is: ',num2str(all_ncf(itile).bi)])
%              disp(['bj is: ',num2str(all_ncf(itile).bi)])
              xax_offset = (all_ncf(itile).bi - 1)*all_ncf(itile).sNx;
              yax_offset = (all_ncf(itile).bj - 1)*all_ncf(itile).sNy;
              ncglob_pos = ''; % Position of tile in the global array
              nctile_size= ''; % Size of the tile data
              
              for jj = 1:length(myvars(ivar).dim_names)
                  if debugMode
                      disp(['operating on ',myvars(ivar).dim_names{jj}])
                  end
                  doff = 1;
                  if jj > 1
                      ncglob_pos = sprintf('%s,',ncglob_pos);
                      nctile_size=sprintf('%s,',nctile_size);
                  end
                  
                  if strcmp(myvars(ivar).dim_names{jj},'X') == 1
              %        disp(['getting dims for X'])
                      dlen = myvars(ivar).dim_sizes{jj};
                      nctile_size=sprintf('%s%d',nctile_size,dlen);
                      
                      doff = xax_offset + doff;
                      dlen = xax_offset + dlen;
                  elseif strcmp(myvars(ivar).dim_names{jj},'Y') == 1
              %        disp(['getting dims for Y'])
                      dlen = myvars(ivar).dim_sizes{jj};
                      nctile_size=sprintf('%s%d',nctile_size,dlen);
                      
                      doff = yax_offset + doff;
                      dlen = yax_offset + dlen;
                  elseif strcmp(myvars(ivar).dim_names{jj},'Xp1') == 1
               %       disp(['getting dims for Xp1'])
                      dlen = myvars(ivar).dim_sizes{jj};
                      % Test if the "new" dlen will be greater than xmax
                      if (xax_offset + dlen)==Xmax+1;
               %           disp('Adjusting x')
               %           We dont want to read in the whole of this
               %           particular variable from the tile file
                          dlen=Xmax;
                          nctile_size=sprintf('%s%d',nctile_size,tx);
                          doff = xax_offset + doff;
                      else % else carry on as usual
                          nctile_size=sprintf('%s%d',nctile_size,dlen);
                          dlen = xax_offset + dlen;
                          doff = xax_offset + doff;
                      end                   
                  elseif strcmp(myvars(ivar).dim_names{jj},'Yp1') == 1
               %       disp(['getting dims for Yp1'])
                      dlen = myvars(ivar).dim_sizes{jj};
                      % Test if the "new" dlen will be greater than ymax
                      if (yax_offset + dlen)==Ymax+1
               %           disp('Adjusting y') 
               %           We dont want to read in the whole of this
               %           particular variable from the tile file
                          dlen=Ymax;
                          nctile_size=sprintf('%s%d',nctile_size,ty);
                          doff = yax_offset + doff;
                      else % else carry on as usual
                          nctile_size=sprintf('%s%d',nctile_size,dlen);
                          dlen = yax_offset + dlen;
                          doff = yax_offset + doff;
                      end 
                  else
            %          disp(['getting dims for other: ',myvars(ivar).dim_names{jj}])
                      dlen = myvars(ivar).dim_sizes{jj};
                      nctile_size=sprintf('%s%d',nctile_size,dlen);
                  end
% Added section to make matlab and unix gluemncs act the same!
% Xp1 and Yp1 respected in the former but not the latter, both now produce 
% length(X)==length(Xp1) and length(Y)==length(Y) by dropping last value,
% which in Xp1's case was 360, starting from 0 ie duplicate.
% Assume Xmax x Ymax grid & uses bi_max/bj_max
%                   Ymax
%                   doff_test=Ymax-(ty)+1
%                   if doff==doff_test;
%                       dlen
%                       if dlen==Ymax+1 
%                           disp('Adjusting y')
%                           dlen=Ymax;
%                           nctile_size=strrep(nctile_size,num2str(ty+1),num2str(ty))
%                       end
%                   end
%                   
%                   if dlen==Xmax+1;
%                       disp('Adjusting x')
%                       dlen=Xmax;
%                       nctile_size=strrep(nctile_size,num2str(tx+1),num2str(tx))
%                   end
                  
                 ncglob_pos = sprintf('%s%d%s%d',ncglob_pos,doff,':',dlen);
              end
              
              
% This loop produces two critical strings. Examples shown are for 100 timesteps 
%   of a 128x64x15 model with 4 tiles. 
% nctile_size is the size of the tile data loaded using SNCTOOLS nc_varget 
%   and is essentially the "count" input required to load a subset of data.
% nctile_size = 100,15,32,64 
%
% ncglob_pos is the position of the individual tile in the global array
%    used to assemble said global array before dumping to file
% ncglob_pos  = 1:100,1:15,1:32,1:64     % bottom left
% ncglob_pos  = 1:100,1:15,1:32,65:128   % bottom right
% ncglob_pos  = 1:100,1:15,33:64,1:64    % top left
% ncglob_pos  = 1:100,1:15,33:64,65:128  % top right


% Old way, solely based on NetCDF toolbox.           
%            comm = sprintf( ...
%                'fonc{''%s''}(%s) = nct{''%s''}(%s);', ...
%                myvars(ivar).name, ncglob_pos, myvars(ivar).name, diml_in );
            
% Load data into global array using SNCTOOLS.....
            comm = sprintf('ncdata_global(%s) = nc_varget(''%s'',''%s'',[%s],[%s]);', ...
                ncglob_pos, all_ncf(itile).name, myvars(ivar).name,...
                num2str(zeros(1,length(myvars(ivar).dim_names))),... % produce zeros as start point e.g. [0 0 0 0] for 4D var
                nctile_size);
            if debugMode
                disp(['    Create global array:  ' comm ]);
            end
            eval(comm);    
          end
      end
      
      % Produce vector of sizes of data e.g. [nt,nz,ny,nx]
      % Works ok for 2+D data, but for 1D needs special treatment.
      if length(myvars(ivar).dim_names)==1
          vardim_length=length(ncdata_global);
      else
          vardim_length=size(ncdata_global);
      end
      
% .....Then write out using SNCTOOLS.       
      comm=sprintf('nc_varput(foutnm,''%s'',ncdata_global,[%s],[%s]);',...
          myvars(ivar).name,...
          num2str(zeros(1,length(myvars(ivar).dim_names))),... % produce zeros as start point e.g. [0 0 0 0] for 4D var
          num2str(vardim_length)); % vector of global data array e.g. [nt,nz,nx,ny]
      if debugMode
          disp([ 'Write out global data:  ' comm ]);
      end
      eval(comm);
  end

  % Clears/closes any still-open netcdf files
  clear functions

elseif all_ncf(1).exch == 2 % UNTESTED: as I dont have any cube-sphere data available
  warning('The Cube Sphere routine is untested!!')
  
  %  exch "2":
  Xmax = 0;
  Ymax = 0;
  for ii = 1:length(all_ncf)
    Xmax = max(Xmax, (all_ncf(ii).exch2_tbasex + all_ncf(ii).sNx));
    Ymax = max(Ymax, (all_ncf(ii).exch2_tbasey + all_ncf(ii).sNy));
  end
  
  horzdim = struct('names',{},'lens',{});
  horzdim(1).names = { 'X' };  horzdim(1).lens = { Xmax     };
  horzdim(3).names = { 'Y' };  horzdim(3).lens = { Ymax     };
  horzdim(5).names = { 'T' };  horzdim(5).lens = { 0        };
  horzdim(4).names = {'Yp1'};  horzdim(4).lens = { Ymax     };
  horzdim(2).names = {'Xp1'};  horzdim(2).lens = { Xmax     };
%   horzdim(2).names = {'Xp1'};  horzdim(2).lens = { Xmax + 1 }; % Match output from unix routines
%   horzdim(4).names = {'Yp1'};  horzdim(4).lens = { Ymax + 1 }; % Match output from unix routines
  
  iseq = 0;
  foutnm = sprintf(fout, iseq);
%  fonc = netcdf(foutnm,'clobber');  % Should append-or-create!
  nc_create_empty(foutnm,'clobber');

  % Add Global Attributes to global data file
  disp(['Writing global attributes to ',foutnm])
  
  % Add Global Attributes to global data file
  % Have to interact with mexnc here, but not sure if this is backwards
  % compatible, besides, we already have global attributes in fnm_tile_info
  %
  % [ncid,status]=mexnc('open',all_ncf(1).name,'nowrite'); if status, error ( mexnc ( 'strerror', status ) ), end
  % [ngatts, status] = mexnc('INQ_NATTS',ncid);            if status, error ( mexnc ( 'strerror', status ) ), end
  %
  % atts=struct;
  % for attnum = 1:ngatts
  %     atts(attnum).attname = mexnc('inq_attname',ncid, nc_global, attnum-1);
  %     atts(attnum).attval = nc_attget(all_ncf(1).name, nc_global,atts(attnum).attname);
  %     if strcmp(atts(attnum).attname,'tile_number'); atts(attnum).attval='global'; end % Edit tilenum to be 'global'
  % end
  % mexnc('close',ncid);
  
  % Edit tilenum to be 'global' and a character datatype
  fnm_tile_info.Attribute(find(strcmp({fnm_tile_info.Attribute.Name},'tile_number'))).Value='global';
  fnm_tile_info.Attribute(find(strcmp({fnm_tile_info.Attribute.Name},'tile_number'))).Nctype=2;
  fnm_tile_info.Attribute(find(strcmp({fnm_tile_info.Attribute.Name},'tile_number'))).Datatype='char';
  
  for attnum = 1:length(fnm_tile_info.Attribute)
      if strcmp(fnm_tile_info.Attribute(attnum).Datatype,'char') % write out character attributes
          comm=sprintf('nc_attput(foutnm,nc_global,''%s'',''%s'');',...
              fnm_tile_info.Attribute(attnum).Name,...
              num2str(fnm_tile_info.Attribute(attnum).Value));
          eval(comm)
      else % write out numerical attributes
          comm=sprintf('nc_attput(foutnm,nc_global,''%s'',%s);',...
              fnm_tile_info.Attribute(attnum).Name,...
              num2str(fnm_tile_info.Attribute(attnum).Value));
          eval(comm)
      end
      % disp(['Adding global attributes: ', comm])
  end

  for idim = 1:length(mydims)
    dname = mydims(idim).names{1};
    ind = find(strcmp(dname,[horzdim.names]));
    if ~isempty(ind)
      dlen = horzdim(ind).lens{1};
    else
      dlen = mydims(idim).lens{1};
    end
%    comm = sprintf('fonc(''%s'') = %d;',dname,dlen);
    comm = sprintf('nc_adddim(foutnm,''%s'',%d);',dname,dlen);
%    disp([ 'Define dimensions:  ' comm ]);
    eval(comm);
  end

  for ivar = 1:length(myvars)
    % Not sure if this is necessary as I dont have any cube-sphere data available to test
    % and see if the iface dimension also appears in myvars(ivar).dim_names
    % ...presumably not if it is given special treatment
    %comm = sprintf('fonc{''%s''} = nc%s( ',myvars(ivar).name,myvars(ivar).type);
    id = 1;
    %comm = [ comm sprintf('''%s''',myvars(ivar).dim_names{id}) ];
    comm = [ sprintf('''%s''',myvars(ivar).dim_names{id}) ];
    for id = 2:length(myvars(ivar).dim_names)
      dname = myvars(ivar).dim_names{id};
      if (dname(1) == 'Y') && (myvars(ivar).has_horiz == 1)
        comm = [ comm sprintf(',''%s''','iface') ];
      end
      comm = [ comm sprintf(',''%s''',dname) ];
    end
    %comm = [ comm ' );' ];
    exch2_dims=comm; 
    
    comm= sprintf('nc_addvar(foutnm,struct(''Name'',''%s'',''Datatype'',''%s'',''Dimension'',{{%s}});',...
        myvars(ivar).name,...
        strrep(strrep(myvars(ivar).type,'char','text'),'long','int'),...
        exch2_dims);
    eval(comm);

    for iat = 1:length(myvars(ivar).atts.names)
%       comm = sprintf( ...
%           'fonc{''%s''}.%s = nc%s( myvars(ivar).atts.data{iat} );', ...
%           myvars(ivar).name, ...
%           myvars(ivar).atts.names{iat}, ...
%           myvars(ivar).atts.types{iat} );
        comm=sprintf('nc_attput(foutnm,''%s'',''%s'',''%s'');',...
            myvars(ivar).name,...
            myvars(ivar).atts.names{iat},...
            myvars(ivar).atts.data{iat});
%      disp(['    Add variable attributes:  ' comm ]);
        eval(comm);
    end
  end

  %  Here is where we need to check the output file size and start
  %  another file in the sequence, if necessary.
  
  for ivar = 1:length(myvars)
    if ~license('test','Distrib_Computing_Toolbox');
        fprintf('  Copying variable:   %s\n',myvars(ivar).name)
    end
    
    ncdata_global=[]; % Temporary data file contains accumulated tiles before output

    for itile = 1:length(all_ncf)

      if (myvars(ivar).has_horiz == 1) || (itile == 1)
        
        clear nct;
        nct = all_ncf(itile).nc{1};
        xax_offset = all_ncf(itile).exch2_tbasex;
        yax_offset = all_ncf(itile).exch2_tbasey;
        diml_tin = '';
        diml_res = '';
        diml_in  = '';
        ncglob_pos = '';
        if length(myvars(ivar).dim_names) < 2
%           comm = sprintf( ...
%               'fonc{''%s''}(%s%d) = nct{''%s''}(:);', ...
%               myvars(ivar).name, '1:', myvars(ivar).dim_sizes{1}, ...
%               myvars(ivar).name );
           % Accumulate tiles into temporary global file   
           comm = sprintf('ncdata_global(%s) = nct{''%s''}(:);', ...
                ['1:',num2str(myvars(ivar).dim_sizes{1})],myvars(ivar).name);
           % disp([ '    ' comm ]);
           eval(comm);
        else
          for jj = 1:length(myvars(ivar).dim_names)
            doff = 1;
            if jj > 1
              diml_tin = sprintf('%s,',diml_tin);
              diml_res = sprintf('%s,',diml_res);
              diml_in  = sprintf('%s,',diml_in);
              ncglob_pos = sprintf('%s,',ncglob_pos);
            end
            dnam = myvars(ivar).dim_names{jj};
            dlen = myvars(ivar).dim_sizes{jj};
            dlenr = dlen;
            fchar = myvars(ivar).dim_names{jj}(1);
            % disp(['       fchar = ' fchar '  ' myvars(ivar).dim_names{jj}]);
            if strcmp(dnam(1),'X') == 1
              doff = xax_offset + doff;
              dlen = xax_offset + dlen;
            end
            if strcmp(dnam(1),'Y') == 1
              diml_res = sprintf('%s%s',diml_res, '[],');
              diml_in  = sprintf('%s%s',diml_in, ':,');
              ncglob_pos = sprintf('%s%d%s',ncglob_pos,all_ncf(itile).exch2_myFace,',');
              doff = yax_offset + doff;
              dlen = yax_offset + dlen;
            end
            diml_tin = sprintf('%s%s',diml_tin, ':');
            diml_res = sprintf('%s%d',diml_res, dlenr);
            diml_in  = sprintf('%s%s',diml_in, ':');
            ncglob_pos = sprintf('%s%d%s%d',ncglob_pos,doff,':',dlen);
          end
          
          comm = sprintf('ncdata_global(%s) = reshape(nct{''%s''}(%s),%s);',... 
              ncglob_pos,myvars(ivar).name, diml_tin, diml_res);
          eval(comm);
          
%           comm = sprintf('fonc{''%s''}(%s) = tmp(%s);',...
%               myvars(ivar).name,ncglob_pos, diml_in );
%           disp([ '    ' comm ]);
%           eval(comm);
        end
      end
    end
    
    if length(myvars(ivar).dim_names)==1
        vardim_length=length(ncdata_global);
    else
        vardim_length=size(ncdata_global);
    end
    
    comm=sprintf('nc_varput(foutnm,''%s'',ncdata_global,[%s],[%s]);',...
        myvars(ivar).name,...
        num2str(zeros(1,length(myvars(ivar).dim_names))),... % produce zeros as start point e.g. [0 0 0 0] for 4D var
        num2str(vardim_length)); % vector of sizes of data to append e.g. [nt,nz,nx,ny]
    %disp([ 'Write out global data:  ' comm ]);
    eval(comm);
  end
  % end

  % fonc = close(fonc);
end

function string = cell2str(cellstr)
%CELL2STR Convert a cell array of string dimensions to a string.

ncols = size(cellstr,2);
for i=1:ncols-1
    cellstr(:,i) = cellfun(@(x)['''' strrep(x,'''','''''') ''','],...
        cellstr(:,i),'UniformOutput',false);
end
if ncols>0
    cellstr(:,ncols) = cellfun(@(x)['''' strrep(x,'''','''''') ''';'],...
        cellstr(:,ncols),'UniformOutput',false);
end
cellstr = cellstr';
string = ['{' cellstr{:} '}'];
