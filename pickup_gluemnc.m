
% function result=pickup_gluemnc(new_nTx,new_nTy);
% pickup_gluemnc is a script to glue pickup tiles together - no way to
% split them into smaller parts yet....
%           nx and ny are variables defined in SIZE.h that tell MITgcm
%           how many tiles in x- and y-directions to split the domain.
%
%
%
%
% Calls slightly modified gluemnc.m containing function "scanforfiles"
% modified from rdmds.m. 
% mnc_assembly.m used for standard output, pickup_assembly.m used for pickups 
% rdmnc.m used to extract attributes
%
% JML (16/11/2007)

tiledir=dir('*tile*'); % checks for old tile directories
if ~isempty(tiledir)
    error('There appear to be tile directories present, delete before continuing')
end

old_tiles=dir('pickup.*.nc');
old_no_of_tiles=length(old_tiles);
pickup=rdmnc(old_tiles(1).name);
old_nSx=pickup.attributes.global.nSx;old_nSy=pickup.attributes.global.nSy;
old_nPx=pickup.attributes.global.nPx;old_nPy=pickup.attributes.global.nPy;

old_nTx=old_nSx*old_nPx;
old_nTy=old_nSy*old_nPy;

if old_no_of_tiles~=old_nTx*old_nTy;
    error('Number of starting tiles does not equal nTx*nTy!')
end

new_nTx=2;new_nTy=2; %%%%% EVENTUALLY THIS WILL BE INPUT
new_no_of_tiles=new_nTx*new_nTy;

if new_no_of_tiles > old_no_of_tiles
    error('Cannot increase number of tiles...YET!')
    return
end

% Build tile domain matrix, for example a 16:4 conversion looks like:
%                   S
%
% [     1   |   2   |   3   |   4   ...
%           1               2
%       5   |   6   |   7   |   8   ...
%
%       9   |   10  |   11  |   12  ...
%           3               4
%       13  |   14  |   15  |   16  ]
%
%                   N
%
domain=ones(old_nTy,old_nTx);
for i=1:old_nTy;
    domain(i,:)=((1:old_nTx)+(i-1)*old_nTx);
end

% Number of divisions in x- and y-directions
dx=old_nTx/new_nTx; dy=old_nTy/new_nTy;

count=0;
tiledir=[];

for ytiles=1:dy
    for xtiles=1:dx
        count=count+1;
        uy=ytiles*dy;
        ly=((ytiles-1)*dy)+1; 
        ux=xtiles*dx; 
        lx=((xtiles-1)*dx)+1; 
        
        tiledir(count,:)=sprintf('tile%03d',count); eval(['!mkdir ',tiledir(count,:)]);
        eval(sprintf('tile%03d=domain(ly:uy,lx:ux);',count))
        
        for yfile=ly:uy
            for xfile=lx:ux
                eval(['!mv ',sprintf('pickup*.t%03d.nc',domain(yfile,xfile)),' ',tiledir(count,:),'/.'])
            end
        end
    end
end

% Asssemble tiles
for i=1:new_no_of_tiles
    cd(sprintf('tile%03d',new_no_of_tiles))
        
    for k=1:5000 % Iteratively find first tile number present   
        files=dir(sprintf('*.t%03d.nc',k)); % gets all file details in structural array
        if ~isempty(files); break; end
    end

    tmp=rdmnc(files(1).name);
    attributes=tmp.attributes.global;
    clear tmp
   
    for j=1:size(files,1);
        varname=files(j).name(1:end-19);
        [nc_out]=gluemnc(varname);
    end

    globfiles=dir('*.glob.nc');
   
    if ~isempty(globfiles)
        for l=1:length(globfiles)
            eval(['!mv ',globfiles(l).name,' ../',globfiles(l).name(1:end-7),sprintf('t%03d.nc',i)]);
        end
    else
        error('global files not found')
    end
    cd ..
end