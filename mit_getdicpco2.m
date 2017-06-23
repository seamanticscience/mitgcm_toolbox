function [atm_co2]=mit_getdicpco2(data_file,statesteps)
% Use DIC_int1 from data.dic to decide which way to handle the atmos_co2
% files (if any):
% dic_int1:
%  0=use default 278.d-6
%  1=use constant value - dic_pco2, read in from data.dic
%  2=read in from file
%  3=interact with atmospheric box (use dic_pco2 as initial atmos. value)
%  4=atmospheric box with specified emissions (bonus code by jml)
%
% Returns array of values: time, moles of atm carbon and atm pco2.
% JML various dates
if nargin==0
    data_file='data';
end

% get deltat, filearr and list of timesteps
deltat = mit_getparm(data_file,'deltaTClock');
path_string=fileparts(data_file); if ~isempty(path_string); path_string=[path_string,'/']; end
data_dic=[path_string,'data.dic'];
dicint = mit_getparm(data_dic,'DIC_int1');

if nargin<2
    statesteps=mit_timesteps(strrep(data_dic,'data.dic','state'),deltat);
end

if ~exist(data_dic,'file')
    error('Could not find data.dic file')
end

% Check to see if atm_file exists, if so read that in as we have done this
% loop before, else read in pickups or co2atmos.dat or whatever...
atm_file=strrep(strrep(statesteps.filearr,'''state','dic_atmos'),'.glob.nc''','.txt'); % ie dic_atmos.yyyymmdd.txt
if exist(atm_file,'file')
    % dic_atmos files concatonated into one file
    atm_co2 = load(atm_file);
    atm_co2(:,1)=(atm_co2(:,1)*deltat)/(60*60*24*360); % timesteps
    clear dic_atmos atm_file atmos_files
else
    switch dicint
        case 0
            % use default pco2=278.d-6 (atmos mols = 4.9206e+16)
            atm_co2=[statesteps.tim,(ones(length(statesteps.tim),1).*4.9206e+16),...
                (ones(length(statesteps.tim),1).*278.d-6)];
        case 1
            % use constant value read in from data.dic
            [~,str]=system(['grep -i dic_pco2 ',data_dic]);
            eval([lower(str),';']);
            %Does not handle values with exponents in.
            %dic_pco2 = mit_getparm('data.dic','dic_pco2');
            % pco2*total_atmos_moles=mols_co2
            atm_co2=[statesteps.tim,(ones(length(statesteps.tim),1).*(dic_pco2*1.77e20)),...
                (ones(length(statesteps.tim),1).*dic_pco2)];
        case 2
            % read co2 data from file called co2atmos.dat
            % dic_int2=number entries to read
            % dic_int3=start timestep,
            % dic_int4=timestep between file entries
            dicint2 = mit_getparm(data_dic,'DIC_int2');
            dicint3 = mit_getparm(data_dic,'DIC_int3');
            dicint4 = mit_getparm(data_dic,'DIC_int4');
            if exist([fileparts(data_file),'/co2atmos.dat'],'file')
               dic_pco2=load([fileparts(data_file),'/co2atmos.dat']); 
            elseif exist([strrep(fileparts(data_file),'input','build'),'/co2atmos.dat'],'file')
               dic_pco2=load([strrep(fileparts(data_file),'input','build'),'/co2atmos.dat']); 
            else
                error('No co2atmos.dat file could be found')
            end
            atmos_tim=((dicint3):dicint4:(dicint3+dicint4*dicint2))/(dicint4);
            % just get the ones that are within the simulation
            [~,ii]=intersect(atmos_tim',statesteps.tim);
            atm_co2=[atmos_tim(ii)',(dic_pco2(ii)*1.77e20),dic_pco2(ii)]; % pco2*total_atmos_moles=mols_co2
            clear atmos_tim ii
        case {3,4}
            % interact with atmospheric box - restart case only considered
            atm_file=strrep(strrep(statesteps.filearr,'''tave','dic_atmos'),'.glob.nc''','.txt'); % ie dic_atmos.yyyymmdd.txt
            if exist(['dic_atmos.',sprintf('%010d',statesteps.timesteps(1))],'file')
                % dic_atmos files still individual
                % find all atmospheric files
                %atmos_files=dir('dic_atmos.*');
                % Create array of correct size
                eval(['load dic_atmos.',sprintf('%010d',statesteps.timesteps(1))])
                atm_co2=NaN(statesteps.kmax,length(dic_atmos)+1);
                % Load files individually into the array
                for i=1:statesteps.kmax
%                    atm_co2(i,1)=(str2double(atmos_files(i).name(11:end))*deltat)/(60*60*24*360);
                    atm_co2(i,1)=statesteps.tim(i);
                    eval(['load dic_atmos.',sprintf('%010d',statesteps.timesteps(i))])
                    atm_co2(i,2:length(dic_atmos)+1)=dic_atmos; % This is correct!!
                end
                clear dic_atmos atmos_files i
            elseif exist(['pickup_dic_co2atm.',sprintf('%010d',statesteps.timesteps(1)),'.data'],'file')
                % new binary dic_co2atm pickup files post ~63m
                % find all atmospheric files
                %atmos_files=dir('pickup_dic_co2atm.*.data');
                % Create array of correct size
                dic_atmos=mit_readfield(['pickup_dic_co2atm.',sprintf('%010d',statesteps.timesteps(1)),'.data'],[1,1,2],'real*8');
                atm_co2=NaN(statesteps.kmax,length(dic_atmos)+1);
                % Load files individually into the array
                for i=1:statesteps.kmax
%                    atm_co2(i,1)=(str2double(atmos_files(i).name(11:end))*deltat)/(60*60*24*360);
                    atm_co2(i,1)=statesteps.tim(i);
                    dic_atmos=mit_readfield(['pickup_dic_co2atm.',sprintf('%010d',statesteps.timesteps(i)),'.data'],[1,1,2],'real*8');
                    atm_co2(i,2:length(dic_atmos)+1)=squeeze(dic_atmos)'; % This is correct!!
                end
                clear dic_atmos atmos_files i
            elseif exist(atm_file,'file')
                % dic_atmos files concatonated into one file
                atm_co2 = load(atm_file);
                atm_co2(:,1)=(atm_co2(:,1)*deltat)/(60*60*24*360); % timesteps
                clear dic_atmos atm_file atmos_files
            else
                error('Could not detect atmospheric CO2 box output')
            end
        otherwise
            error('Could not detect how you handled atmospheric CO2')
    end
end
