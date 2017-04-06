function [fname,atm_co2]=mit_savedicpco2(data_file,statesteps)

if nargin==0
    data_file='data';
end

if ~exist(data_file,'file')
    error('Could not find "data" file')
end

deltat = mit_getparm(data_file,'deltaTClock');
data_dic=[fileparts(data_file),'/data.dic'];
dicint = mit_getparm(data_dic,'DIC_int1');

if nargin<2
%    statesteps=mit_timesteps('tave',deltat);
    statesteps=mit_timesteps('state',deltat);
end

if dicint >= 3 
    if exist(['dic_atmos.',sprintf('%010d',statesteps.timesteps(2))],'file') ||...
           exist(['pickup_dic_co2atm.',sprintf('%010d',statesteps.timesteps(2)),'.data'],'file') 
    %do nothing
    % else link files into global directory
    elseif exist(['../build/dic_atmos.',sprintf('%010d',statesteps.timesteps(2))],'file')
        !ln -sf ../build/dic_atmos.* .
        % overwrite files actually stored in the input directory
        !ln -sf ../input/dic_atmos.* .
        if exist('../input/oldpickups','dir')
            !ln -sf ../input/oldpickups/dic_atmos.* .
        end
    elseif exist(['../build/pickup_dic_co2atm.',sprintf('%010d',statesteps.timesteps(2)),'.data'],'file')
        !ln -sf ../build/pickup_dic_co2atm.*.data .
        !ln -s ../input/pickup_dic_co2atm.*.data .
        if exist('../input/oldpickups','dir')
            !ln -sf ../input/oldpickups/pickup_dic_co2atm.*.data .
        end
    else
        error('Cannot find atmospheric CO2 files for accumulation')
    end
end
    % Get the atmospheric values from the model output
    [atm_co2]=mit_getdicpco2(data_file,statesteps);
    
    %convert first column of years back to interation number
    atm_co2_i=[(atm_co2(:,1)*(60*60*24*360))/deltat,atm_co2(:,2:end)];
    
    %write to dic_atmos file.
    fname=strrep(strrep(statesteps.filearr,'''state','dic_atmos'),'.glob.nc''','.txt');
    save(fname,'atm_co2_i','-ascii','-double')
    disp([fname,' created successfully'])