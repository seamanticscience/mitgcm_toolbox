% This is a script that executes a series of commands for
% extracting da Silva surface data, interpolating, and
% adjusting (for conservation).
%
% Uses data files: bathymetry.bin pmask.bin
% Creates data files:  lev_clim_*.bin lev_monthly_*.bin
% Creates arrays: n/a
% Uses m-scripts: setgrid grid_sphere finegrid gen_hxhy xyrecur_*
%                 edtopo hfac
%
% Created  08/15/99 by adcroft@mit.edu
% Modified 11/11/99 by adcroft@mit.edu
% Maintained by adcroft@mit.edu, abiastoch@ucsd.edu

clear all
format compact
more off

load FMT.mat
load HGRID.mat nxc nyc lon_lo lon_hi lat_lo lat_hi xc yc xf yf zA
dd=2;
%fd=0;

msk=rdslice('pmask.bin',[nxc nyc],1,fmt,Ieee);
zmsk=msk;
msk( find(msk==0) )=NaN;
iim1=[nxc 1:nxc-1];
jjm1=[nyc 1:nyc-1];
msku=msk.*msk(iim1,:);
mskv=msk.*msk(:,jjm1);

% JML 08/01/08
% Dasilva data in netcdf format, therefore mods made...no longer uses exract_dasilva_30x30.m
var=ncload('dasilva.cdf');
for ii=1:length(var);
	vartmp=var{ii};
	eval(['test1=size(',vartmp,');']);
	
	if length(test1)==3
		eval([vartmp,'=permute(',vartmp,',[2,3,1]);']); % puts in xyt order
        eval([vartmp,'(find(',vartmp,'==-10e9))=NaN;']) % masks out null value
	end
end

% SST 
for month=1:12,
% [Tx,xh,yh]=extract_dasilva_30x30('taux',month,lon_lo-dd,lon_hi+dd,lat_lo-dd,lat_hi+dd);
 t(:,:,month)=interp2(Y,X,squeeze(sst(:,:,month))',yc,xf);
 t(:,:,month)=xyexpand(t(:,:,month),5).*msk;
end
t( find(t==-1e10) )=0;
wrslice(['dasilva_sst.bin'],t,1,fmt,Ieee);

% WIND SPEED
for month=1:12,
% [Tx,xh,yh]=extract_dasilva_30x30('taux',month,lon_lo-dd,lon_hi+dd,lat_lo-dd,lat_hi+dd);
 w(:,:,month)=interp2(Y,X,squeeze(w3(:,:,month))',yc,xf);
 w(:,:,month)=xyexpand(w(:,:,month),5).*msku;
end
w( find(w==-1e10) )=0;
wrslice(['dasilva_speed.bin'],w,1,fmt,Ieee);

% ZONAL WIND STRESS
for month=1:12,
% [Tx,xh,yh]=extract_dasilva_30x30('taux',month,lon_lo-dd,lon_hi+dd,lat_lo-dd,lat_hi+dd);
 tx(:,:,month)=interp2(Y,X,squeeze(taux(:,:,month))',yc,xf);
 tx(:,:,month)=xyexpand(tx(:,:,month),5).*msku;
end
tx( find(tx==-1e10) )=0;
wrslice(['dasilva_taux.bin'],tx,1,fmt,Ieee);

% MERIIDIONAL WINDSTRESS
for month=1:12,
% [Ty,xh,yh]=extract_dasilva_30x30('tauy',month,lon_lo-dd,lon_hi+dd,lat_lo-dd,lat_hi+dd);
 ty(:,:,month)=interp2(Y,X,squeeze(tauy(:,:,month))',yf,xc);
 ty(:,:,month)=xyexpand(ty(:,:,month),5).*mskv;
end
wrslice(['dasilva_tauy.bin'],ty,1,fmt,Ieee);

% SURFACE NET HEAT FLUX
for month=1:12,
% [Q,xh,yh]=extract_dasilva_30x30('Qnet',month,lon_lo-dd,lon_hi+dd,lat_lo-dd,lat_hi+dd);
 q(:,:,month)=interp2(Y,X,squeeze(netheat(:,:,month))',yc,xc);
 q(:,:,month)=xyexpand(q(:,:,month),5).*msk;
end
q( isnan(q) )=0;
q=-q; % Positive upwards
meanQ=sum(sum( mean(q,3).*zA )) ./ sum(sum( zmsk.*zA ));
q( find(q==0) )=NaN; qc=q-meanQ; qc( isnan(qc) )=0; q(isnan(q))=0;
%report('Adjusted mean Q by %f W/m^2\n',meanQ)
%wrslice(['dasilva_qnet_corrected.bin'],qc,1,fmt,Ieee);
wrslice(['dasilva_qnet.bin'],q,1,fmt,Ieee);

% SURFACE EVAP MINUS PRECIP
for month=1:12,
% [EmP,xh,yh]=extract_dasilva_30x30('EmP',month,lon_lo-dd,lon_hi+dd,lat_lo-dd,lat_hi+dd);
 e(:,:,month)=interp2(Y,X,squeeze(emp(:,:,month))',yc,xc);
 e(:,:,month)=xyexpand(e(:,:,month),5).*msk;
end
e( isnan(e) )=0;
e=e/(1000*3*3600); % Convert from mm/(3 hrs) to m/s
meanEmP=sum(sum( mean(e,3).*zA )) ./ sum(sum( zmsk.*zA ));
e( find(e==0) )=NaN; ec=e-meanEmP; ec( isnan(ec) )=0; e(isnan(e))=0;
%report('Adjusted mean E-P-R by %f cm/yr\n',meanEmP*(100*86400*365))
wrslice(['dasilva_emp.bin'],e,1,fmt,Ieee);
%wrslice(['dasilva_emp_corrected.bin'],ec,1,fmt,Ieee);

% ATMOSPHERIC SPECIFIC HUMIDITY
for month=1:12,
 shair(:,:,month)=interp2(Y,X,squeeze(qair(:,:,month))',yc,xc);
 shair(:,:,month)=xyexpand(shair(:,:,month),5).*msk;
end
shair( isnan(shair) )=0;
meanQair=sum(sum( mean(shair,3).*zA )) ./ sum(sum( zmsk.*zA ));
shair( find(shair==0) )=NaN; shair_corr=shair-meanQair; shair_corr( isnan(shair_corr) )=0;
shair(isnan(shair))=0;
%report('Adjusted mean Specific Humidity (air) by %f cm/yr\n',meanQair)
wrslice(['dasilva_shair.bin'],shair,1,fmt,Ieee);
%wrslice(['dasilva_shair_corrected.bin'],ec,1,fmt,Ieee);

% SEA SURFACE SPECIFIC HUMIDITY
for month=1:12,
 shsea(:,:,month)=interp2(Y,X,squeeze(qsea(:,:,month))',yc,xc);
 shsea(:,:,month)=xyexpand(shsea(:,:,month),5).*msk;
end
shsea( isnan(shsea) )=0;
meanQsea=sum(sum( mean(shsea,3).*zA )) ./ sum(sum( zmsk.*zA ));
shsea( find(shsea==0) )=NaN; shsea_corr=shsea-meanQsea; shsea_corr( isnan(shsea_corr) )=0;
shsea( isnan(shsea) )=0;
%report('Adjusted mean Specific Humidity (sea) by %f cm/yr\n',meanQsea)
wrslice(['dasilva_shsea.bin'],shsea,1,fmt,Ieee);
%wrslice(['dasilva_shsea_corrected.bin'],ec,1,fmt,Ieee);

% SEA SURFACE ATMOSPHERIC TEMPERATURE
for month=1:12,
 s(:,:,month)=interp2(Y,X,squeeze(sat(:,:,month))',yc,xc);
 s(:,:,month)=xyexpand(s(:,:,month),5).*msk;
end
s( isnan(s) )=0;
meansat=sum(sum( mean(s,3).*zA )) ./ sum(sum( zmsk.*zA ));
s( find(s==0) )=NaN; s_corr=s-meansat; s_corr( isnan(s_corr) )=0;s( isnan(s) )=0;
%report('Adjusted mean Surface Air Temp by %f cm/yr\n',meansat)
wrslice(['dasilva_sat.bin'],s,1,fmt,Ieee);
%wrslice(['dasilva_sat_corrected.bin'],ec,1,fmt,Ieee);

