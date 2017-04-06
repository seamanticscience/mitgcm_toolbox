function [vals] = grepread(file,varargin)
% VALS = GREPREAD(FILE,EXPR1,...);
%
% Extracts the expressions "expr1","expr2",... from the file "file".
% This assumes output in the standard form defined by the MITgcm
% monitor package and is not a replacement for TEXTREAD.
%
% e.g.
% vals=grepread('output.txt','time_secondsf','ke_mean','ke_max');
% plot(vals(:,1)/86400,vals(:,2:3));

% $Header: /u/gcmpack/MITgcm/utils/matlab/grepread.m,v 1.3 2007/02/17 23:49:43 jmc Exp $
% $Name:  $

if nargin < 2
 error('You must supply a filename and at least one search expression!')
end

tfile=sprintf('grepexpr%15.15f',rand);
for k=1:nargin-1;
    var=strrep(varargin{k},' ','\ ');
%   eval(['!grep ' varargin{k} ' ' file ' | sed s/.\*=// | sed s/NAN/1.23456789/! ' tfile])
    if strcmp(var,'pCo2') || strcmp(var,'atmos C') % Not "current, old diff" style
      disp(['Searching for ' varargin{k}]);
      eval(['!grep ' var ' ' file ' | sed s/.\*pCo2// >' tfile]);
      [tmp1,tmp2]=textread(tfile,'%f%f');
      tmp=[tmp1,tmp2];
    else
      disp(['Searching for ' varargin{k}]);
      eval(['!grep ' var ' ' file ' | sed s/.\*' var '// | sed s/.\*diff// >' tfile]);
      [tmp1,tmp2]=textread(tfile,'%f%f');
      % Calculate diff from current - old
      tmp3=tmp1-tmp2;
      tmp=[tmp1,tmp2,tmp3];
    end
  var=strrep(varargin{k},' ','');
  var=strrep(var,'-','');
  eval(['vals.',var,'=tmp;']);
  clear tmp* var
  delete(tfile)
end
