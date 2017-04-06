% m-file: mit_timesteps.m
% sets up the the variables tname and timesteps for the 4x4 global runs

% $Header: /u/gcmpack/MITgcm/verification/tutorial_global_oce_latlon/diags_matlab/mit_timesteps.m,v 1.3 2006/08/12 20:25:13 jmc Exp $
% $Name:  $

% The second of two get timestep retreivals....need to extract timestep
% info from diagnostic files.....

% if strcmp(dname,'baseline.000')
%   timesteps = [0:11]'*36e3; % one hundred years
%   timesteps = [timesteps; [12:30]'*72e3]; % one hundred years (baseline experiment)
% elseif strcmp(dname,'etopo5.000')
%   timesteps = [0:30]'*72e3; % one hundred years (etopo5 experiments)
% else
%   if meanfields
%     [dummy, timesteps] = rdmds('uVeltave',NaN);
%     if isempty(timesteps)
%       meanfields=0;
%       error('There are no averaged fields, use meanfields=0')
%     end
%     timesteps = [0 timesteps];
%   else
%     [dummy, timesteps] = rdmds('U',NaN);
%   end
%   timesteps = timesteps';
%   clear dummy
% end
diags = rdmnc('surf*iag.*.glob.nc','iter');
timesteps_diag = diags.iter;

nt2 = length(timesteps_diag);
kt2 = 1:nt2;
% tname2 = cell(nt2,1);
% for k2 = kt2
%   tname2{k} = sprintf('%010i',timesteps_diag(k2));
% end

%% alternative (unflexible way)
% $$$   timesteps = [0:11]'*36e3; % one hundred years
% $$$   timesteps = [timesteps; [12:30]'*72e3]; % one hundred years (baseline experiment)
% $$$ clear timesteps tname kt nt
% $$$ timesteps = [0:30]'*72e3; % one hundred years (etopo5 experiments)
% $$$ %tincr = 20; timesteps = [0:36]'*20; % one hundred years
% $$$ nt = length(timesteps);
% $$$ kt = 1:nt;
% $$$ tname = cell(nt,1);
% $$$ for k = kt
% $$$   tname{k} = sprintf('%010i',timesteps(k));
% $$$ end

% create a time axis
oneday = 3600*24;
onemonth =oneday*30;
oneyr=onemonth*12;
msg_spinup = dir('spinup.tim');
clear tim_diag
if isempty(msg_spinup)
  tim_diag = deltat*timesteps_diag';
else
  tim_diag = load('spinup.tim');
  itim2 = find(diff(tim_diag) == 0);
  tim_diag(itim2) = [];
end

% create reasonable unit; looks at range rather that magnetude
if (max(tim_diag(:))-min(tim_diag(:))) == 0 % rare, but does occur assume years
  tim_diag = tim_diag/oneyr;
  timeunit_diag = 'yrs';
  tuname_diag = 'year';
elseif (max(tim_diag(:))-min(tim_diag(:)))/oneday <= 360
  ndayiter0_diag=(diags.iter(2)-diags.iter(1))*deltat/oneday;
  %nyears=(max(tim)-min(tim)+ndayiter0)/360;
  tim_diag = (tim_diag/oneday)-((min(tim_diag)/oneday)-ndayiter0_diag);
  timeunit_diag = 'days';
  tuname_diag = 'day';
elseif (max(tim_diag(:))-min(tim_diag(:)))/onemonth <= 120
  nmonthiter0_diag=(diags.iter(2)-diags.iter(1))*deltat/onemonth;
  %nyears=(max(tim)-min(tim)+nmonthiter0)/12;
  tim_diag = (tim_diag/onemonth)-((min(tim_diag)/onemonth)-nmonthiter0_diag);
  timeunit_diag = 'months';
  tuname_diag = 'month';
else 
  tim_diag = tim_diag/oneyr;
  timeunit_diag = 'yrs';
  tuname_diag = 'year';
end


if ~exist('kmax_diag','var')
  kmax_diag=max(kt2);
end
disp(['kmax_diag = ' num2str(kmax_diag) ', diplayed time = ' ...
      num2str(tim_diag(kmax_diag)) ' ' timeunit_diag])

