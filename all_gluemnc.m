function all_gluemnc(issafe)
% all_gluemnc assembles all files from an MITgcm model sequence.
%
% Calls slightly modified gluemnc.m containing function "scanforfiles"
% modified from rdmds.m. mnc_assembly.m remains unchanged.
%
% JML (16/11/2007)

% set default as safe mode
% safe (1) or "clinically insane" (0) timeseries concatenation.
if nargin==0
    issafe=1;
end 

files=dir('*.t001.nc'); % finds different kinds of files
%tmp=rdmnc(files(end).name);
%attributes=tmp.attributes.global;
%clear tmp

% what if tile files have already been cleared?
if isempty(files)
    files=dir('*.glob.nc');
end

for j=1:size(files,1);
    if ~strcmp(files(j).name(1:6),'pickup')
        if strcmp(files(j).name(1:4),'grid')
            varname='grid';
            gname='grid.glob.nc';
        else
            varname=files(j).name(1:end-19);
            gname=[files(j).name(1:end-7),'glob.nc'];
            % Give gluemnc an Iter number to aviod errors
            tmp=files(j).name;
            iter(j,:)=sprintf('%010d',str2double(tmp(end-17:end-8)));
        end
        if ~exist(gname,'file');
            if strcmp(files(j).name(1:4),'grid')
                % dont supply iter because grid files are named "grid.t???.nc"
                [nc_out]=gluemnc(varname);
            else
                [nc_out]=gluemnc(varname,iter(j,:));
            end
            % close any open files
            clear functions
            fclose('all');
        else
            disp(['filename: ',gname,' already exists check manually to confirm it is ok!'])
        end
    end
end

%%   % tidy up leftover links
%!find . -type l | xargs rm
system('find . -type l  -exec unlink {} \;');

%% Now check for timeseries broken across different files
if exist('gluemnc_cat_timeseries.m','file') % prerequisite
    datenow=datestr(now,'yyyymmdd');
    itern=nan(size(iter,1),1);
    for i=1:length(iter);
        itern(i)=str2double(iter(i,:));
    end
    
    itern(isnan(itern))=[]; % remove nans
    if length(unique(itern))~=1
        % create list of file prefixes
        timeser_files=dir(['*.',iter(1,:),'.glob.nc']);
        
        for i=1:length(timeser_files)
            outfnm=[timeser_files(i).name(1:end-19),'.',datenow,'.glob.nc'];
            gluemnc_cat_timeseries(timeser_files(i).name(1:end-19),outfnm,issafe)
            
            % close any open files
            clear functions
            fclose('all');
        end
        
        % tidy up leftover timeseries fragments
        for i=1:length(itern);
            eval(['system(''find . -name "*.',sprintf('%010d',itern(i)),'.glob.nc" -exec rm {} \;'');'])
        end
    end
end

%% check for atmospheric box and whether those files need accumulating
if exist('../input/data','file')
   deltat = mit_getparm('../input/data','deltaTClock'); 
   statesteps=mit_timesteps('state',deltat);
   [fname,atm_co2]=mit_savedicpco2('../input/data',statesteps);
elseif exist('../input/data.orig','file')
   deltat = mit_getparm('../input/data.orig','deltaTClock');
   statesteps=mit_timesteps('state',deltat);
   [fname,atm_co2]=mit_savedicpco2('../input/data.orig',statesteps);
elseif exist('../input/data.001','file')
   deltat = mit_getparm('../input/data.001','deltaTClock');
   statesteps=mit_timesteps('state',deltat);
   [fname,atm_co2]=mit_savedicpco2('../input/data.001',statesteps);
end
% close any open files
clear functions
fclose('all');

% tidy up leftover links again
%!find . -type l | xargs rm
system('find . -type l  -exec unlink {} \;');