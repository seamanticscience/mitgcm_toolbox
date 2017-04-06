
% Bull_gluemnc is a script called during assembly of the parts of a
% parallel model-run sequence. (used to be on the Bull HPC)
%
% Calls slightly modified gluemnc.m containing function "scanforfiles"
% modified from rdmds.m. mnc_assembly.m has been adapted.
% rdmnc.m used to extract attributes
%
% JML (16/11/2007)
% Attempted Parallelisation using PCT - JML (31/70/09)



partdir=dir('*part_*'); % sees how many directories need glueing

if isempty(partdir) %, disp('error finding any part directories');end
   % Not using partdirs so must be a global directory with all files
   % present or linked.
   clear partdir
   partdir.name='.';
end
    
    
for i=1:size(partdir,1);
   cd(partdir(i).name)
   
   files=dir('*.t001.nc'); % finds different kinds of files
%    tmp=rdmnc(files(end).name);
%    attributes=tmp.attributes.global;
%    clear tmp

   % Build cell array of input, output and iteration numbers first
   varname{size(files,1)}=[];
   gname{size(files,1)}=[];
   iter{size(files,1)}=[];
   
   for j=1:size(files,1);
       if ~strcmp(files(j).name(1:6),'pickup') 
            if strcmp(files(j).name(1:4),'grid')
                varname{j}='grid';
                gname{j}='grid.glob.nc';
                iter{j}='0';
            else
                varname{j}=files(j).name(1:end-19);
                gname{j}=[files(j).name(1:end-7),'glob.nc'];
                % Give gluemnc an Iter number (char) to aviod errors
                tmp=files(j).name;
                itertmp=tmp(end-17:end-8);
                if length(itertmp)<10
                    tmp=10-length(itertmp);
                    for z=1:tmp
                        itertmp=['0',iter];
                    end
                end
                iter{j}=itertmp;
            end
       end
   end
     
%    if license('test','Distrib_Computing_Toolbox') ==1; % if license exists
%        check=license('checkout','Distrib_Computing_Toolbox'); % check to see if there is a spare copy
%        if check==1 && matlabpool('size')==0
%            %checked out uni license for parallel computing toolbox
%            %no existing matlab pool present
%            %therefore start matlabpool
%            matlabpool
%        end
%    end
   
   % remove the details of already-glued global files so that they dont
   % get re-glued....have to do this outside the parallel loop.
   for j=1:size(files,1);
       if exist(gname{j},'file');
           disp(['File ',gname{j},' already exists. Removing details so that it is not overwritten.'])
           gname{j}=[];
           varname{j}=[];
           iter{j}=[];
       elseif strcmp(gname{j}(1:7),'monitor');
           disp(['No need to glue monitor files at present and there is only one tile anyway!'])
           gname{j}=[];
           varname{j}=[];
           iter{j}=[];
       end
   end
   
   % Purge gname, varname and iter of empty cells. Not to worry, grid
   % iteration is '0' so that arrays dont get out of sync, but this is
   % never used...
   gname(cellfun('isempty',gname))=[];
   varname(cellfun('isempty',varname))=[];
   iter(cellfun('isempty',iter))=[];
   
   parfor j=1:size(gname,2);
      [nc_out]=gluemnc(varname{j},iter{j});
   end
      
   cd ..
end


% if check==1 && matlabpool('size')>0
%       %checked out uni license for parallel computing toolbox
%       %existing matlab pool present
%       %therefore close matlabpool
%       matlabpool close
% end