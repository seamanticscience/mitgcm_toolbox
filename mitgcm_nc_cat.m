function mitgcm_nc_cat(file1,file2)
%nc_cat Concatenate netCDF files along unlimited dimension.
%   nc_cat(file1,file2) concatenates files along the unlimited dimension.
%   The UNIQUE records of file2 are added to file1.  Variable not defined
%   along the unlimited dimension are left untouched.
%
%   Example:
%       nc_create_empty('f1.nc');
%       nc_adddim('f1.nc','time',0);
%       v.Name = 'time';
%       v.Dimension = {'time'};
%       nc_addvar('f1.nc',v);
%       v.Name = 'money';
%       v.Dimension = {'time'};
%       nc_addvar('f1.nc',v);
%       copyfile('f1.nc','f2.nc');
%
%       % Populate the first file.
%       buf = struct('Name','Data');
%       buf(1).Name = 'time';  buf(1).Data = [0 1 2];
%       buf(2).Name = 'money';  buf(2).Data = [0 1000 2000];
%       nc_addnewrecs('f1.nc',buf);
%
%       % Now populate the 2nd file.
%       buf(1).Data = [3 4 5 6];
%       buf(2).Data = [3000 4000 5000 6000];
%       nc_addnewrecs('f2.nc',buf);
%
%       % Now concatenate them.
%       nc_cat('f1.nc','f2.nc');
%       data = nc_varget('f1.nc','money');
%
%
%   See also nc_addnewrecs.

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

% % Load times from each file and decide if the values are unique or not
% % using intersect. Likely to only be 1 but could be more such as
% % diagnostics that have different output frequencies.
% t1=nc_varget(file1,'T');
% t2=nc_varget(file2,'T');
% crossover=intersect(t1,t2); 
% 
% % Only load into memory the unique portions of file2 as is done for jml's
% % unix netcdf file accumulation routine. E.g. for state files this means
% % that the initial state of each startup in file2 is ignored in favour of
% % the final state from file1. "Start" is zero referenced, so for a complete
% % list of unique times length(crossover)=0, and will load from the
% % begining!
% start = length(crossover);
% count = length(t2)-start;
% b = nc_getbuffer(file2,start,count);
% nc_addrecs(file1,b);

% Alternatively could use nc_addnewrecs....
buff = nc_getbuffer(file2);
nc_addnewrecs(file1,buff);

function record_variable = find_record_variable(ncfile)
%FIND_RECORD_VARIABLE  return name of record variable
not_found = true;
info1 = nc_info(ncfile);
for j = 1:numel(info1.Dataset)
    if info1.Dataset(j).Unlimited ...
            && (numel(info1.Dataset(j).Dimension) == 1) ...
            && strcmp(info1.Dataset(j).Name,info1.Dataset(j).Dimension{1})
        idx = j;
        not_found = false;
    end
end
if not_found
    error('SNCTOOLS:nc_cat:noRecordVariableFound', ...
        'Could not find a record variable in %s.', ncfile);
end
record_variable = info1.Dataset(idx).Name;
