function []=check2dmask(datafile,month)
% check2dmask(datafile,recnum)
%
% Load 2D data from record number recnum 
% in "datafile", calculate a missing
% value mask and compare to pmask.bin.
% They should always be the same!!!
%
% Created  11/11/99 by adcroft@mit.edu
% Modified          by
% Maintained by adcroft@mit.edu, abiastoch@ucsd.edu

load HGRID.mat nxc nyc
load FMT.mat

pmaskfile='pmask.bin';
Pmsk=rdda(pmaskfile,[nxc nyc],1,fmt,Ieee);

T=rdda(datafile,[nxc nyc],month,fmt,Ieee);
Tmsk=ones([nxc nyc]);
Tmsk( find(T==0) )=0;

mismatch=find(Pmsk ~= Tmsk);
%mismatch=find(Pmsk < Tmsk);

if isempty(mismatch)
 report('check2dmask("%s",%i) passed  OK\n',datafile,month)
else
 report('I found %i points where %s mis-matched %s\n',...
   prod(size(mismatch)),datafile,pmaskfile);
 error('check2dmask FAILED!!')
end
