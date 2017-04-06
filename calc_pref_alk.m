% Calculate preformed alkalinity by multiple regression to salinity and PO,
% which removes the influence of high-alk waters upwelled in the Southern
% Ocean
grd=mit_loadgrid;
tavesteps=mit_timesteps('tave');

% Load Data
tave=rdmnc('tave.20k.glob.nc','Ttave','Stave');
tave.Stave=change(tave.Stave,'==',0,NaN);
tave.Ttave=change(tave.Ttave,'==',0,NaN);
ptracer=rdmnc('ptr_tave.20k.glob.nc','alk','po4','o2');
ptracer.alk=change(ptracer.alk,'==',0,NaN);
ptracer.po4=change(ptracer.po4,'==',0,NaN);
ptracer.o2=change(ptracer.o2,'==',0,NaN);
hfacc_surf=grd.hfacc(:,:,1);

o2=ptracer.o2(:,:,1);
po4=ptracer.po4(:,:,1);
alk=ptracer.alk(:,:,1);
sal=tave.Stave(:,:,1);
% Define PO
po_surf=o2+(170*po4);
po=ptracer.o2+(170*ptracer.po4);
% Mulitple Regression
[b,bint,r,rint,stats]=regress(alk(:),[hfacc_surf(:),sal(:),po_surf(:)]);

% Define Preformed Alk
alk0_surf=b(1)+(b(2).*sal)+(b(3).*po_surf);
% EDIT: 02/08/2011 ALK at the surface is ALL Preformed by definition.
alk0_surf=alk;

alk0=b(1)+(b(2).*tave.Stave)+(b(3).*po);

% Calculate 1SD uncertainty in Alkalinity
%alkstd=nanstd(alk0(:));

% Calculate 1SD of the difference between surface Alk0 and Alk
dalk_mean=nanmean(alk0_surf(:)-alk(:)); % 3.4205e-16 mol/m3
dalk_std=nanstd(alk0_surf(:)-alk(:)); % 8.4739 umol/kg or 0.0086815 mol/m3

disp(['Coefficients are:'])
b
disp(['With an Rsquared value of:'])
stats(1)

% Usually:
% b =
%        0.1114
%      0.062426
%      0.095896
% 
% With an Rsquared value of:
% 
% ans =
% 
%       0.97223
plot([2.1,2.6],[2.1,2.6],'k-','LineWidth',2)
hold on
scatter(ptracer.alk(:),alk0(:),30,'b','filled')
%hold on
scatter(alk(:),alk0_surf(:),30,'r','filled')
%plot([2.1,2.6],[2.1+dalk_std,2.6+dalk_std],'k--','LineWidth',2)
%plot([2.1,2.6],[2.1-dalk_std,2.6-dalk_std],'k--','LineWidth',2)
set(gca,'XLim',[2.1,2.6],'YLim',[2.1,2.6])
box on
xlabel('Total Alkalinity (mol m^-^3)')
ylabel('Preformed Alkalinity (mol m^-^3)')
