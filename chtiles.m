function chtiles(diags,nIter0,oldtiles,newtiles)
%
% purpose: uses gluemnc to alter the number of tiles which the domain is
%       broken into. 
%
%	   diags:  diagnostics name
%          nIter0: 10-digit iterate #
% EXAMPLE:
%	   foo = chtiles('pickup','0004248000',16,4);
%

for ii=1:newtiles
    newdir=['tile',num2str(newtiles)];
    eval(['!mkdir ',tiledir])
    
    
end
    