%
% Load grid and control data for carbon, temperature and SSH (for EKE)
%

% Load model domain
grid=mit_loadgrid;
grid.dzc=[grid.zc(1);diff(grid.zc)];

%get ocean basin masks
grid=mit_oceanmasks(grid);

% Geophysical parameters
omega = 2*pi/(3600*23+56*60);
g=9.8;
[X,Y] = meshgrid(grid.long,grid.latg);
f = 2*omega*sind(Y');

% tave=rdmnc([modeldir,'/tave.',num2str(start_iter,'%010d'),'.glob.nc']);
% ptracers=rdmnc([modeldir,'/ptr_tave.',num2str(start_iter,'%010d'),'.glob.nc']);

tave=rdmnc('tave.20150317.glob.nc',21600000);
ptracers=rdmnc('ptr_tave.20150317.glob.nc',21600000);

theta = tave.Ttave;
ssh   = inpaint_nans(tave.ETAtave,4);
dic   = ptracers.dic; 

dsshdx=diff(-ssh,1,1); dsshdx(grid.nx,:)=ssh(grid.nx,:)-ssh(1,:);
dsshdy=diff(-ssh,1,2); dsshdy(:,grid.ny)=NaN(grid.nx,1);

dsshdx=dsshdx./grid.dxg;
dsshdy=dsshdy./grid.dyg;

ug = -dsshdx.*(ones(size(dsshdx)).*g./f);
vg =  dsshdy.*(ones(size(dsshdx)).*g./f);

gv=sqrt(ug.*ug+vg.*vg);
eke=0.5.*(ug.*ug+vg.*vg);