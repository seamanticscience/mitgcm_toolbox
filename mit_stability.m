% Calculate MITgcm stability values to assess new timestep

% deltatmom=mit_getparm('data','deltaTmom');
% deltatclock=mit_getparm('data','deltaTClock');
% deltattracer=mit_getparm('data','deltaTtracer');
% beta=mit_getparm('data','beta');
% f0=mit_getparm('data','f0');
% dphi=mit_getparm('data','delX');
% dphi=dphi(1);
% ah=mit_getparm('data','viscAh');
% az=mit_getparm('data','viscAz');
% kh=mit_getparm('data','diffKhT');
% kz=mit_getparm('data','diffKzT');
% dz=mit_getparm('data','delZ');
% dz=min(dz); % conservative estimates at small dz

roe=6400000; % radius of earth (m)
phi=80; % ref latitude....conservative estimates at high latitude

mlt=(pi*(ah/beta))^(1/3); % Munk Layer thickness (mlt > dxml to resolve properly)
mlt=mlt/1000; % convert to km
dxml=mean([(roe*cosd(45)-roe*cosd(45+dphi)),(roe*cosd(45-dphi)-roe*cosd(45))])/1000; % Resolution at mid-latitude (km)

dx=mean([(roe*cosd(phi)-roe*cosd(phi+dphi)),(roe*cosd(phi-dphi)-roe*cosd(phi))]);
dx2=dx*dx;
dz2=dz*dz;

% Terms for stability of horizontal/vertical Laplacian friction/dissipation
% slnn < 0.3
slah=(4*ah*deltatmom)/dx2;
slaz=(4*az*deltatmom)/dz2;

slkh=(4*kh*deltattracer)/dx2;
slkz=(4*kz*deltattracer)/dz2;

% Term for stability of inertial oscillations
% si < 1 or 0.5
si=(f0^2)*(deltatmom^2);

% Term for stability of Advective CFL for extremem horizontal flow ~2ms-1
% sa < 0.5
uvel=2;
sa=(uvel*deltatmom)/dx;

% Term for stability of internal gravity waves propogating with speed
% ~10ms-1 (sc < 0.5) or ~2ms-1 (sc < 0.25)
cgvel=10;
sc=(cgvel*deltatmom)/dx;

% Produce report
if mlt>dxml
    fprintf(1,'Munk layer width of %0.5f km GREATER THAN resolution at mid-latitude %0.5f km, ensuring frictional boundary layer is well resolved.\n',mlt,dxml)
else 
    fprintf(1,'WARNING: Munk layer width of %0.5f km LESS THAN resolution at mid-latitude %0.5f km, frictional boundary layer IS NOT well resolved.\n',mlt,dxml)
end

if slah<0.3
    fprintf(1,'Horizontal Laplacian Friction stability term (slah) of %0.5f LESS THAN stability criteria of 0.3.\n',slah)
else 
    fprintf(1,'WARNING: Horizontal Laplacian Friction stability term (slah) of %0.5f  GREATER THAN stability criteria of 0.3.\n',slah)
end

if slaz<0.3
    fprintf(1,'Vertical Laplacian Dissipation stability term (slaz) of %0.5f LESS THAN stability criteria of 0.3.\n',slaz)
else 
    fprintf(1,'WARNING: Vertical Laplacian Dissipation stability term (slaz) of %0.5f  GREATER THAN stability criteria of 0.3.\n',slaz)
end

if slkh<0.5
    fprintf(1,'Horizontal Diffusion stability term (slkh) of %0.5f LESS THAN stability criteria of 0.5.\n',slkh)
else 
    fprintf(1,'WARNING: Horizontal Diffusion stability term (slkh) of %0.5f  GREATER THAN stability criteria of 0.5.\n',slkh)
end

if slkz<0.5
    fprintf(1,'Vertical Diffusion stability term (slkz) of %0.5f LESS THAN stability criteria of 0.5.\n',slkz)
else 
    fprintf(1,'WARNING: Vertical Diffusion stability term (slkz) of %0.5f  GREATER THAN stability criteria of 0.5.\n',slkz)
end

if si<1
    fprintf(1,'Term for stability of Inertial Oscillations (si) of %0.5f LESS THAN stability criteria of 1 (or 0.5).\n',si)
else 
    fprintf(1,'WARNING: Term for stability of Inertial Oscillations (si) of %0.5f  GREATER THAN stability criteria of 1 (or 0.5).\n',si)
end

if sa<0.5
    fprintf(1,'Term for stability of Advective CFL with extreme flow (%d ms-1) (sa) of %0.5f LESS THAN stability criteria of 0.5.\n',uvel,sa)
else 
    fprintf(1,'WARNING: Term for stability of Advective CFL with extreme flow (%d ms-1) (sa) of %0.5f GREATER THAN stability criteria of 0.5.\n',uvel,sa)
end

if sc<0.5
    fprintf(1,'Term for stability of Inertial Gravity Waves propergating with extreme flow (%d ms-1) (sc) of %0.5f LESS THAN stability criteria of 0.5 (or 0.25).\n',cgvel,sc)
else 
    fprintf(1,'WARNING: Term for stability of Inertial Gravity Waves propergating with extreme flow (%d ms-1) (sc) of %0.5f GREATER THAN stability criteria of 0.5 (or 0.25).\n',cgvel,sc)
end
