function [phnew,co2s,pco2,ff,ffp,k0,k1,k2,co3]=carbon_chem(pressc,t,s,dic,ta,pt,sit,phguess,grd)
% New efficient pCO2 solver, Mick Follows, Taka Ito, Stephanie Dutkiewicz
% z = depth (m)
% t = temperature (degrees C)
% s = salinity (PSU)
% dic = total inorganic carbon (mol/m^3)
% ta = total alkalinity (eq/m^3)
% pt = inorganic phosphate (mol/^3)
% sit = inorganic silicate (mol/^3)
% phguess = initial guess at pH

% ---------------------------------------------------------------------
% Get CO2 coefficients

% ff=nan(grd.nx,grd.ny,grd.nz);      ffp=nan(grd.nx,grd.ny,grd.nz);
% k0=nan(grd.nx,grd.ny,grd.nz);      k1=nan(grd.nx,grd.ny,grd.nz);
% k2=nan(grd.nx,grd.ny,grd.nz);      hnew=nan(grd.nx,grd.ny,grd.nz);
% phnew=nan(grd.nx,grd.ny,grd.nz);   co2s=nan(grd.nx,grd.ny,grd.nz);
% pco2=nan(grd.nx,grd.ny,grd.nz);    

hguess=nan(grd.nx,grd.ny,grd.nz);

% JML 14/08/11 - Pressure in bars should be positive! Now supplied
% directly from calc_pco2 routine in bars.
% % Convert Depth into pressure in bars
% pressc = ones(grd.nx,grd.ny,grd.nz) - 0.1.*z;

% Change units from the input of mol/m^3 -> mol/kg:
% (1 mol/m^3) x (1 m^3/1024.5 kg)
permil=1/1024.5 ;
% where the ocean's mean surface density is 1024.5 kg/m^3
% Note: mol/kg are actually what the body of this routine uses
% for calculations. Units are reconverted back to mol/m^3 at the
% end of this routine.
% To convert input in mol/m^3 -> mol/kg

[ff,ffp,k0,k1,k2,kb,k1p,k2p,k3p,ksi,kw,ks,kf,bt,st,ft,Ksp_TP_Calc]...
    = carbon_coeffs(t,s,pressc);

for i=1:grd.nx
    for j=1:grd.ny
        for k=1:grd.nz
            if ~isnan(grd.hfacc(i,j,k))
                %mick - new approx method
                %mick - make estimate of htotal (hydrogen ion conc) using
                %mick appromate estimate of CA, carbonate alkalinity
                hguess(i,j,k) = 10.0^(-phguess(i,j,k)); % conc of H+ in mol/kg
                hguess_old=0;
                count=0;
                %for p=1:1 %jon - iterate to get ph in the deep ocean really good
                while (abs(hguess(i,j,k)-hguess_old)>3e-13)
                    %mick - first estimate borate contribution using guess for [H+]
                    bohg = (bt(i,j,k)*kb(i,j,k))/(hguess(i,j,k)+kb(i,j,k));
                    
                    %mick - first estimate of contribution from phosphate
                    %mick based on Dickson and Goyet
                    stuff = (hguess(i,j,k).^3)...
                        + (k1p(i,j,k).*hguess(i,j,k).^2)...
                        + (k1p(i,j,k).*k2p(i,j,k).*hguess(i,j,k))...
                        + (k1p(i,j,k).*k2p(i,j,k).*k3p(i,j,k));
                    h3po4g = (pt(i,j,k).*permil.*(hguess(i,j,k).^3))./stuff;
                    h2po4g = (pt(i,j,k).*permil.*k1p(i,j,k).*(hguess(i,j,k).^2))./stuff;
                    hpo4g = (pt(i,j,k).*permil.*k1p(i,j,k).*k2p(i,j,k).*hguess(i,j,k))./stuff;
                    po4g = (pt(i,j,k).*permil.*k1p(i,j,k).*k2p(i,j,k).*k3p(i,j,k))./stuff;
                    
                    %mick - estimate contribution from silicate
                    %mick based on Dickson and Goyet
                    siooh3g = (sit(i,j,k).*permil.*ksi(i,j,k))./(ksi(i,j,k) + hguess(i,j,k));
                    
                    %mick - now estimate carbonate alkalinity
                    cag = ta(i,j,k).*permil - bohg - (kw(i,j,k)./hguess(i,j,k)) + hguess(i,j,k)...
                        - hpo4g - 2.*po4g + h3po4g - siooh3g - h2po4g;
                    
                    %mick - now evaluate better guess of hydrogen ion conc
                    %mick htotal = [H+], hydrogen ion conc
                    gamm = (dic(i,j,k).*permil)./cag;
                    stuff = ((1.0-gamm).*(1.0-gamm).*k1(i,j,k).*k1(i,j,k))-(4.0.*k1(i,j,k).*k2(i,j,k).*(1.0-2.0.*gamm));
                    hguess_old=hguess(i,j,k);
                    hguess(i,j,k) = 0.5.*(((gamm-1.0).*k1(i,j,k))+sqrt(stuff));
                    count=count+1;
                    if count==100
                        disp('pH finder was terminated at 100 iterations')
                        break
                    end
                end
            end
        end
    end
end

%mick - return update pH to main routine
phnew = -log10(hguess);

%mick - now determine [CO2*]
co2s = dic./(1.0 + (k1./hguess) + (k1.*k2./(hguess.*hguess)));

% NOW EVALUATE CO32-, carbonate ion concentration
% used in determination of calcite compensation depth
% Karsten Friis & Mick - Sep 2004
co3 = (k1.*k2.*dic)./(hguess.*hguess + k1.*hguess + k1.*k2);
% ---------------------------------------------------------------
% surface pCO2 (following Dickson and Goyet, DOE...)
pco2 = co2s./ff;
%pco2(i,j,k) = co2s(i,j,k)/ffp;
% ----------------------------------------------------------------

% Reconvert from mol/kg -> mol/m^3
%	pt=pt./permil;
%	sit=sit./permil;
%	ta=ta./permil;
%	dic=dic./permil;
return