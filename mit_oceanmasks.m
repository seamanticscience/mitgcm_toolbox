function grid = mit_oceanmasks(grid,plotit)
% function grid = mit_oceanmasks(grid)
% create masks for the individual oceans
% hfacc/s/w have to be available for this
  
% $Header: /u/gcmpack/MITgcm/verification/tutorial_global_oce_latlon/diags_matlab/mit_oceanmasks.m,v 1.3 2006/08/12 20:25:13 jmc Exp $
% $Name:  $

  if grid.nx == 128 && grid.ny == 64 
  % masks for the zonal velocities:
  % hfacw
  
  % Atlantic mask
  grid.atlantic_hfacw = grid.hfacw;
  ixsow = 9:104; % Southern and south Indian
  iysow = 1:36;
  grid.atlantic_hfacw(ixsow,iysow,:) = NaN;
  ixnpw = 1:99; % North Pacific
  iynpw = 36:50;
  grid.atlantic_hfacw(ixnpw,iynpw,:) = NaN;
  ixnnp = 40:85; % The very northern north Pacific
  iynnp = 50:60;
  grid.atlantic_hfacw(ixnnp,iynnp,:) = NaN;
  
  % Pacific mask
  grid.pacific_hfacw = grid.hfacw;
  ixaw = 104:128; % most of the atlantic
  iyaw = 1:64;
  grid.pacific_hfacw(ixaw,iyaw,:) = NaN;
  ixnaw = 99:104; % baffin bay and gulf of Mexico
  iynaw = 36:64;
  grid.pacific_hfacw(ixnaw,iynaw,:) = NaN;
  ixiw = 1:41; % most of the indic
  iyiw = 1:64;
  grid.pacific_hfacw(ixiw,iyiw,:) = NaN;
  ixauw = 41:55;% water south and west of Australia
  iyauw = 1:25;
  grid.pacific_hfacw(ixauw,iyauw,:) = NaN;
  
  % masks for the meridional overturning stream functions:
  % hfacs
  
  % Atlantic mask
  grid.atlantic_hfacs = grid.hfacs;
  ixso = 9:104; % Southern and south Indian
  iyso = 1:36;
  grid.atlantic_hfacs(ixso,iyso,:) = NaN;
  ixnp = 1:97; % North Pacific
  iynp = 36:50;
  grid.atlantic_hfacs(ixnp,iynp,:) = NaN;
  ixnnp = 40:85; % The very northern north Pacific
  iynnp = 50:60;
  grid.atlantic_hfacs(ixnnp,iynnp,:) = NaN;
  ixsoc = 1:128; % complete southern ocean
  iysoc = 1:20; % south of latg(12) = -34 degN;
  grid.atlantic_hfacs(ixsoc,iysoc,:) = NaN;
  
  % Pacific mask
  grid.pacific_hfacs = grid.hfacs;
  ixa = 104:128; % most of the atlantic
  iya = 1:64;
  grid.pacific_hfacs(ixa,iya,:) = NaN;
  ixna = 99:104; % baffin bay and gulf of Mexico
  iyna = 36:64;
  grid.pacific_hfacs(ixna,iyna,:) = NaN;
  ixna2 = 98:100;
  iyna2 = 40:42;
  grid.pacific_hfacs(ixna2,iyna2,:) = NaN;
  ixi = 1:41; % most of the indic
  iyi = 1:64;
  grid.pacific_hfacs(ixi,iyi,:) = NaN;
  ixau = 41:55;% water south and west of Australia
  iyau = 1:25;
  grid.pacific_hfacs(ixau,iyau,:) = NaN;

  % Indian ocean
  grid.indic_hfacs = grid.hfacs;
  ixsoi = 1:128; % southern ocean
  iysoi = 1:20;
  grid.indic_hfacs(ixsoi,iysoi,:) = NaN;
  ixap = [1:10 46:128]; % atlantic and pacific
  iyap = 1:64;  
  grid.indic_hfacs(ixap,iyap,:) = NaN;
  iynps = 55:64; % north polar sea
  grid.indic_hfacs(:,iynps,:) = NaN;
  ixscs = 39:46; % south china sea
  iyscs = 31:45;
  grid.indic_hfacs(ixscs,iyscs,:) = NaN;

  % Southern Ocean
  grid.so_hfacs = grid.hfacs;
  % Mask everything above 36Deg S
  iysom = 20:64;
  grid.so_hfacs(:,iysom,:) = NaN;

  if exist('plotit','var')
    figure;
    spy(isnan(grid.hfacs(:,:,1)'),'kx');
    hold on; 
    spy(~isnan(grid.atlantic_hfacs(:,:,1)')); 
    spy(~isnan(grid.indic_hfacs(:,:,1)'),'g.'); 
    spy(~isnan(grid.so_hfacs(:,:,1)'),'r.'); 
    spy(~isnan(grid.pacific_hfacs(:,:,1)'),'y.');
    title('hfacs')
    axis xy
  end
  
  % masks for zonal averages of C-properties:
  % hfacc
  
  % Atlantic mask
  grid.atlantic_hfacc = grid.hfacc;
  ixso = 10:104; % Southern and south Indian
  iyso = 1:36;
  grid.atlantic_hfacc(ixso,iyso,:) = NaN;
  ixnp = 1:97; % North Pacific
  iynp = 36:50;
  grid.atlantic_hfacc(ixnp,iynp,:) = NaN;
  ixnnp = 40:85; % The very northern north Pacific
  iynnp = 50:60;
  grid.atlantic_hfacc(ixnnp,iynnp,:) = NaN;
%   ixsoc = 1:128; % complete southern ocean
%   iysoc = 1:20; % south of latg(12) = -34 degN;
%   grid.atlantic_hfacc(ixsoc,iysoc,:) = NaN;
  
  % Pacific mask
  grid.pacific_hfacc = grid.hfacc;
  ixa = 105:128; % most of the atlantic
  iya = 1:64;
  grid.pacific_hfacc(ixa,iya,:) = NaN;
  ixna = 98:104; % baffin bay and gulf of Mexico
  iyna = 36:64;
  grid.pacific_hfacc(ixna,iyna,:) = NaN;
  ixna2 = 98:100;
  iyna2 = 40:42;
  grid.pacific_hfacc(ixna2,iyna2,:) = NaN;
  ixi = 1:41; % most of the indic
  iyi = 1:64;
  grid.pacific_hfacc(ixi,iyi,:) = NaN;
%   ixau = 41:55;% water south and west of Australia
%   iyau = 1:25;
%   grid.pacific_hfacc(ixau,iyau,:) = NaN;
  grid.pacific_hfacc(42:45,1:30)=NaN;
  
  % Indian ocean
  grid.indic_hfacc =grid.hfacc;
%   ixsoi = 1:128; % southern ocean
%   iysoi = 1:20;
%   grid.indic_hfacc(ixsoi,iysoi,:) = NaN;
  ixap = [1:9 46:128]; % atlantic and pacific
  iyap = 1:64;  
  grid.indic_hfacc(ixap,iyap,:) = NaN;
  iynps = 55:64; % north polar sea
  grid.indic_hfacc(:,iynps,:) = NaN;
  ixscs = 39:46; % south china sea
  iyscs = 31:45;
  grid.indic_hfacc(ixscs,iyscs,:) = NaN;
  
  if exist('plotit','var')
    figure;
    spy(isnan(grid.hfacc(:,:,1)'),'kx'); axis xy
    hold on; 
    spy(~isnan(grid.atlantic_hfacc(:,:,1)')); 
    spy(~isnan(grid.pacific_hfacc(:,:,1)'),'r.');
    spy(~isnan(grid.indic_hfacc(:,:,1)'),'g.'); 
    axis xy
    title('hfacc');

  end
  elseif grid.nx == 360 && grid.ny >= 160
      % masks for the zonal velocities:
  % hfacw
  
  % Atlantic mask
  grid.atlantic_hfacw = grid.hfacw;
  ixsow = 22:292; % Southern and south Indian
  iysow = 1:89;
  grid.atlantic_hfacw(ixsow,iysow,:) = NaN;
  ixnpw = 1:263; % North Pacific
  iynpw = 89:132;
  grid.atlantic_hfacw(ixnpw,iynpw,:) = NaN;
  ixnnp = 100:250; % The very northern north Pacific
  iynnp = 132:160;
  grid.atlantic_hfacw(ixnnp,iynnp,:) = NaN;
  ixeqp = 263:278; % The Eastern equatorial Pacific
  iyeqp = 89:97;
  grid.atlantic_hfacw(ixeqp,iyeqp,:) = NaN;
  
  % Pacific mask
  grid.pacific_hfacw = grid.hfacw;
  ixaw = 292:360; % most of the atlantic
  iyaw = 1:160;
  grid.pacific_hfacw(ixaw,iyaw,:) = NaN;
  ixnaw = 262:292; % baffin bay and gulf of Mexico
  iynaw = 97:160;
  grid.pacific_hfacw(ixnaw,iynaw,:) = NaN;
  ixsgm = 278:292; % south gulf of Mexico
  iysgm = 90:97;
  grid.pacific_hfacw(ixsgm,iysgm,:) = NaN;
  ixiw = 1:101; % most of the indic
  iyiw = 1:160;
  grid.pacific_hfacw(ixiw,iyiw,:) = NaN;
  ixauw = 101:146;% water south and west of Australia
  iyauw = 1:80;
  grid.pacific_hfacw(ixauw,iyauw,:) = NaN;
  
  
  % masks for the meridional overturning stream functions:
  % hfacs
  
  % Atlantic mask
  grid.atlantic_hfacs = grid.hfacs;
  ixsow = 22:292; % Southern and south Indian
  iysow = 1:89;
  grid.atlantic_hfacs(ixsow,iysow,:) = NaN;
  ixnpw = 1:263; % North Pacific
  iynpw = 89:132;
  grid.atlantic_hfacs(ixnpw,iynpw,:) = NaN;
  ixnnp = 100:250; % The very northern north Pacific
  iynnp = 132:160;
  grid.atlantic_hfacs(ixnnp,iynnp,:) = NaN;
  ixeqp = 263:278; % The Eastern equatorial Pacific
  iyeqp = 89:97;
  grid.atlantic_hfacs(ixeqp,iyeqp,:) = NaN;
  
  % Pacific mask
  grid.pacific_hfacs = grid.hfacs;
  ixaw = 292:360; % most of the atlantic
  iyaw = 1:160;
  grid.pacific_hfacs(ixaw,iyaw,:) = NaN;
  ixnaw = 262:292; % baffin bay and gulf of Mexico
  iynaw = 97:160;
  grid.pacific_hfacs(ixnaw,iynaw,:) = NaN;
  ixsgm = 278:292; % south gulf of Mexico
  iysgm = 90:97;
  grid.pacific_hfacs(ixsgm,iysgm,:) = NaN;
  ixiw = 1:101; % most of the indic
  iyiw = 1:160;
  grid.pacific_hfacs(ixiw,iyiw,:) = NaN;
  ixauw = 101:146;% water south and west of Australia
  iyauw = 1:80;
  grid.pacific_hfacs(ixauw,iyauw,:) = NaN;

  % Indian ocean
  grid.indic_hfacs = grid.hfacs;
  ixsoi = 1:360; % southern ocean
  iysoi = 1:45;
  grid.indic_hfacs(ixsoi,iysoi,:) = NaN;
  ixap = [1:22 146:360]; % atlantic and pacific
  iyap = 45:160;  
  grid.indic_hfacs(ixap,iyap,:) = NaN;
  ixwp = 101:146; % western pacific
  iywp = 80:160;  
  grid.indic_hfacs(ixwp,iywp,:) = NaN;
  ixnps = 22:57; % north polar sea & med
  iynps = 100:160; 
  grid.indic_hfacs(ixnps,iynps,:) = NaN;

  % Southern Ocean
  grid.so_hfacs = grid.hfacs;
  % Mask everything above 36Deg S
  iysom = 45:160;
  grid.so_hfacs(:,iysom,:) = NaN;

  if exist('plotit','var')
    figure;
    spy(isnan(grid.hfacs(:,:,1)'),'kx');
    hold on; 
    spy(~isnan(grid.atlantic_hfacs(:,:,1)')); 
    spy(~isnan(grid.indic_hfacs(:,:,1)'),'g.'); 
    spy(~isnan(grid.so_hfacs(:,:,1)'),'r.'); 
    spy(~isnan(grid.pacific_hfacs(:,:,1)'),'y.');
    title('hfacs')
    axis xy
  end
  
  % masks for zonal averages of C-properties:
  % hfacc
  
  % Atlantic mask
  grid.atlantic_hfacc = grid.hfacc;
  ixsow = 22:292; % Southern and south Indian
  iysow = 1:89;
  grid.atlantic_hfacc(ixsow,iysow,:) = NaN;
  ixnpw = 1:263; % North Pacific
  iynpw = 89:132;
  grid.atlantic_hfacc(ixnpw,iynpw,:) = NaN;
  ixnnp = 100:250; % The very northern north Pacific
  iynnp = 132:160;
  grid.atlantic_hfacc(ixnnp,iynnp,:) = NaN;
  ixeqp = 263:278; % The Eastern equatorial Pacific
  iyeqp = 89:97;
  grid.atlantic_hfacc(ixeqp,iyeqp,:) = NaN;
  
  % Pacific mask
  grid.pacific_hfacc = grid.hfacc;
  ixaw = 292:360; % most of the atlantic
  iyaw = 1:160;
  grid.pacific_hfacc(ixaw,iyaw,:) = NaN;
  ixnaw = 262:292; % baffin bay and gulf of Mexico
  iynaw = 97:160;
  grid.pacific_hfacc(ixnaw,iynaw,:) = NaN;
  ixsgm = 278:292; % south gulf of Mexico
  iysgm = 90:97;
  grid.pacific_hfacc(ixsgm,iysgm,:) = NaN;
  ixiw = 1:101; % most of the indic
  iyiw = 1:160;
  grid.pacific_hfacc(ixiw,iyiw,:) = NaN;
  ixauw = 101:146;% water south and west of Australia
  iyauw = 1:80;
  grid.pacific_hfacc(ixauw,iyauw,:) = NaN;

  % Indian ocean
  grid.indic_hfacc = grid.hfacc;
%   ixsoi = 1:360; % southern ocean
%   iysoi = 1:45;
%   grid.indic_hfacc(ixsoi,iysoi,:) = NaN;
  ixap = [1:22 146:360]; % atlantic and pacific
  iyap = 1:160;  
  grid.indic_hfacc(ixap,iyap,:) = NaN;
  ixwp = 101:146; % western pacific
  iywp = 80:160;  
  grid.indic_hfacc(ixwp,iywp,:) = NaN;
  ixnps = 22:57; % north polar sea & med
  iynps = 100:160; 
  grid.indic_hfacc(ixnps,iynps,:) = NaN;

  if exist('plotit','var')
    figure;
    spy(isnan(grid.hfacc(:,:,1)'),'kx'); axis xy
    hold on; 
    spy(~isnan(grid.atlantic_hfacc(:,:,1)')); 
    spy(~isnan(grid.pacific_hfacc(:,:,1)'),'r.');
    spy(~isnan(grid.indic_hfacc(:,:,1)'),'g.'); 
    axis xy
    title('hfacc');
  end
  else
      error('Basin masks not setup for this configuration')
  end
  return
