function []=check3dmask(datafile)
% check3dmask(datafile)
%
% Load 3D data in "datafile", calculate a missing
% value mask and compare to pmask.bin.
% They should always be the same!!!
%
% Created  11/11/99 by adcroft@mit.edu
% Modified          by
% Maintained by adcroft@mit.edu, abiastoch@ucsd.edu

load HGRID.mat nxc nyc
load VGRID.mat nzc
load FMT.mat

pmaskfile='pmask.bin';

T=rdda(datafile,[nxc nyc nzc],1,fmt,Ieee);

Tmsk=ones([nxc nyc nzc]);
Tmsk( find(T==0) )=0;

Pmsk=rdda(pmaskfile,[nxc nyc nzc],1,fmt,Ieee);

mismatch=find(Pmsk ~= Tmsk);

if isempty(mismatch)
 report('check3dmask("%s") passed  OK\n',datafile)
else
 report('I found %i points where %s mis-matched %s\n',...
   prod(size(mismatch)),datafile,pmaskfile);
 error('check3dmask FAILED!!')
end
