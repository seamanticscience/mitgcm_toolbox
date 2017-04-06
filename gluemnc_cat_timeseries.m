function gluemnc_cat_timeseries(var,outfnm,issafe)
% gluemnc_cat_timeseries is called from all_gluemnc (assembles all files
% from an MITgcm model sequence) and glues the timeseries files together.
%
% Result is a file labelled with date rather than iteration number, e.g.
%      tave.20130107.glob.nc
%
%
% JML 20/12/2012

% set default as safe mode
if nargin==0
    error('Requires a var or file prefix to continue')
elseif nargin==1
    issafe=1;
end

if strcmp(var,'pickup') || strcmp(var,'grid')
    disp([var,' cannot be processed with this routine'])
    return
end

% finds different kinds of files based on first timestep file, e.g.
% dic_tave.0000000000.glob.nc
% dic_tave.0000000100.glob.nc
% dic_tave.0000000500.glob.nc
files=dir([var,'.*.glob.nc']);

%outfnm=[var,'.',datestr(now,'yyyymmdd'),'.glob.nc'];
disp(['Concatenating timeseries of ',var,' into ',outfnm]);

%% Accumulate timeseries
if issafe==1
    % Safe way to do it is with mitgcm_nc_cat, which appends the new
    % timesteps from file2 into file 1. Similar to nc_cat but only
    % appends newer timesteps rather than all timesteps, which
    % leads to repetition.
    % Need to cycle through nxfile times appending subsequent files to the first.
%    copyfile(files(1).name,outfnm);
    eval(['!cp ',files(1).name,' ',outfnm]);
    for i=2:length(files) % dont do first file!
        gluemnc_safe_cat(outfnm,files(i).name);
    end
else
    % Alternative way uses John Evan's "clinically insane" nc_cat_a
    % which overwrites the common timesteps in file1 with those from
    % file2. Technically this is more desirable since you will have
    % started file2 from a pickup generated at or near the end of
    % file1, rather than continuing file1 until it runs out and
    % picking up the timeseries n timesteps into file2...
    nc_cat_a({files.name},outfnm,'T')
    
    % Add Global Attributes to data file
    var_info=nc_info(files(1).name);
    for attnum = 1:length(var_info.Attribute)
        if strcmp(var_info.Attribute(attnum).Datatype,'char') % write out character attributes
            comm=sprintf('nc_attput(outfnm,nc_global,''%s'',''%s'');',...
                var_info.Attribute(attnum).Name,...
                num2str(var_info.Attribute(attnum).Value));
            eval(comm)
        else % write out numerical attributes
            comm=sprintf('nc_attput(outfnm,nc_global,''%s'',%s);',...
                var_info.Attribute(attnum).Name,...
                num2str(var_info.Attribute(attnum).Value));
            eval(comm)
        end
        % disp(['Adding global attributes: ', comm])
    end
end

%% Internal Functions
function gluemnc_safe_cat(file1,file2)
    %nc_cat Concatenate netCDF files along unlimited dimension.
    %   nc_cat(file1,file2) concatenates files along the unlimited dimension.
    %   The UNIQUE records of file2 are added to file1.  Variable not defined
    %   along the unlimited dimension are left untouched.
    record_variable = find_record_variable(file1);


    % Verify that the record variable exists in file2.
    if ~nc_isvar(file2,record_variable)
        error('SNCTOOLS:nc_cat:recordVariableNotThere', ...
            'The record variables must have the same name.  Could not find "%s" in %s.', ...
            record_variable, file2);
    end

    % Verify that all the unlimited variables in the 2nd file exist in the first file
    info2 = nc_info(file2);
    for j = 1:numel(info2.Dataset)
        if info2.Dataset(j).Unlimited
            if ~nc_isvar(file1,info2.Dataset(j).Name)
                error('SNCTOOLS:nc_cat:recordVariableMissing', ...
                    'Could not find unlimited variable %s in %s.', ...
                    info2.Dataset(j).Name,file1);
            end
        end
    end
    
    % Load file2 into a structure and append new/unique timesteps to file1
    buff = nc_getbuffer(file2);
    nc_addnewrecs(file1,buff);

    function record_variable = find_record_variable(ncfile)
        %FIND_RECORD_VARIABLE  return name of record variable
        not_found = true;
        info1 = nc_info(ncfile);
        for fj = 1:numel(info1.Dataset)
            if info1.Dataset(fj).Unlimited ...
                    && (numel(info1.Dataset(fj).Dimension) == 1) ...
                    && strcmp(info1.Dataset(fj).Name,info1.Dataset(fj).Dimension{1})
                idx = fj;
                not_found = false;
            end
        end
        if not_found
            error('SNCTOOLS:nc_cat:noRecordVariableFound', ...
                'Could not find a record variable in %s.', ncfile);
        end
        record_variable = info1.Dataset(idx).Name;
    end
end

end

