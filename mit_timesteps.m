function [timestructure]=mit_timesteps(invar,deltat)
% sets up the the variables tname and timesteps for the 4x4 global runs

% $Header: /u/gcmpack/MITgcm/verification/tutorial_global_oce_latlon/diags_matlab/mit_timesteps.m,v 1.3 2006/08/12 20:25:13 jmc Exp $
% $Name:  $


if nargin==1 || isempty(deltat)
  deltat = mit_getparm('data','deltaTClock');
end

clear mex

allfiles=dir([invar,'.*.glob.nc']);
% Find files that have true iterations (eg 20091029), not the control runs (*20k*)
for n=1:length(allfiles)
    itertmp=strfind(allfiles(n).name,'k.glob.nc');
    if isempty(itertmp)
        %                itergood(n)=str2double(allfiles(n).name(end-15:end-8));
        % I know these files are tave.*.glob.nc, so remove the
        % extra bits to leave the file iteration numbers, whether 8
        % digit dates (tave.20112012.glob.nc) or 10 digit model
        % iters (tave.0014400000.glob.nc)
        itergood(n,:)=strrep(strrep(allfiles(n).name,[invar,'.'],''),'.glob.nc','');
    end
end

if strcmp('tave',invar) || ~exist('filearr','var')
    % construct file array
    filearr=[];
    for n=1:size(itergood,1)
        filearr=[filearr,'''',invar,'.',itergood(n,:),'.glob.nc','''',','];
    end
    filearr=filearr(1:end-1); % knock off last comma
end

% Load files, but does not include control run...
fnm_in=rdmnc(filearr(2:end-1),'iter');

timesteps = fnm_in.iter;

nt = length(timesteps);
kt = 1:nt;

% create a time axis
oneday = 3600*24;
onemonth =oneday*30;
oneyr=onemonth*12;
msg_spinup = dir('spinup.tim');

if isempty(msg_spinup)
   tim = deltat*timesteps';
else
  tim = load('spinup.tim');
  itim = find(diff(tim) == 0);
  tim(itim) = [];
end

% create reasonable unit; looks at range rather that magnetude
% if (max(tim(:))-min(tim(:))) == 0 % rare, but does occur assume years
%   tim = tim/oneyr;
%   timeunit = 'yrs';
%   tuname = 'year';
% elseif (max(tim(:))-min(tim(:)))/oneday <= 360
%   ndayiter0=(timesteps(2)-timesteps(1))*deltat/oneday;
%   %nyears=(max(tim)-min(tim)+ndayiter0)/360;
%   tim = (tim/oneday)-((min(tim)/oneday)-ndayiter0);
%   timeunit = 'days';
%   tuname = 'day';
% elseif (max(tim(:))-min(tim(:)))/onemonth <= 120
%   nmonthiter0=(timesteps(2)-timesteps(1))*deltat/onemonth;
%   %nyears=(max(tim)-min(tim)+nmonthiter0)/12;
%   tim = (tim/onemonth)-((min(tim)/onemonth)-nmonthiter0);
%   timeunit = 'months';
%   tuname = 'month';
% else 
  tim = tim/oneyr;
  timeunit = 'yrs';
  tuname = 'year';
%end

kmax=max(kt);

disp(['kmax = ' num2str(kmax) ', diplayed time = ' ...
      num2str(tim(kmax)) ' ' timeunit])
  
timestructure=struct('timesteps',timesteps,...
                        'tim',tim,...
                        'kt',kt,...
                        'kmax',kmax,...
                        'timeunit',timeunit,...
                        'tuname',tuname,...
                        'filearr',filearr);
 
end

