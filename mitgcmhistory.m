function [vals] = mitgcmhistory(file,varargin)
%vals = mitgcmhistory(FILE,EXPR1,...);
%
%Extracts the expressions "expr1","expr2",... from the file "file".
%This assumes output in the standard form defined by the MITgcm
%monitor package and is not a replacement for TEXTREAD.
% 
%e.g.
%&gt;&gt; vals=mitgcmhistory('output.txt','time_secondsf','ke_mean','ke_max');
%&gt;&gt; plot(vals(:,1)/86400,vals(:,2:3));
%
% Written by adcroft@mit.edu, 2001
%$Header:

if nargin lt 2
 error('You must supply a filename and at least one search expression!')
end

tfile=sprintf('/tmp/grepexpr%15.15f',rand);
for k=1:nargin-1;
 try
  eval(['!grep ' varargin{k} ' ' file ' | sed s/.\*=// | sed s/NAN/1.23456789/ &gt; ' tfile])

  % vals(:,k)=textread(tfile,'%f');
  % When output file is from an ongoing integration, one or more of
  % the diagnostics may be missing at the last available time step.
  % The code below accomodates this difference in length.

  if k==1
   vals(:,k)=textread(tfile,'%f');
   lngt=length(vals(:,k));
  else
   tmp=textread(tfile,'%f');
   if abs(length(tmp)-lngt) gt 1
    % try to read one line at a time in order to deal with special case
    % of values like, e.g.: "  -2.9248686233802-321"
    fid=fopen(tfile);
    n=0;
    while(~feof(fid))
     n=n+1;
     tmp2=fgetl(fid);
     val=sscanf(tmp2,'%f');
     if length(val) gt 1
      tmp(n)=0;
     else
      tmp(n)=val;
     end
    end
    fid=fclose(fid);
    tmp=tmp(1:n);
   end
   if abs(length(tmp)-lngt) gt 1
    error(sprintf('An error occured while scanning for: %s',varargin{k}));
   else
    lngt=min(lngt,length(tmp));
    vals(1:lngt,k)=tmp(1:lngt);
   end

  end
  delete(tfile)
 catch
  delete(tfile)
  error(sprintf('An error occured while scanning for: %s',varargin{k}));
 end
end
vals=vals(1:lngt,:);
