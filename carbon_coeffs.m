function [ff,ffp,ak0,ak1,ak2,akb,ak1p,ak2p,ak3p,aksi,akw,aks,akf,bt,st,ft,Ksparag,Kspcalc] = carbon_coeffs(t,s,pressc)
% C
% C     ./==========================================================\
% C     | SUBROUTINE CARBON_COEFFS                                 |
% C     | determine coefficients for surface carbon chemistry      |
% C     | adapted from OCMIP2:  SUBROUTINE CO2CALC                 |
% C     | mick follows, oct 1999                                   |
% c     | minor changes to tidy, swd aug 2002                      |
% c     | MODIFIED FOR PRESSURE DEPENDENCE                         |
% c     | Karsten Friis and Mick Follows 2004                      |
% C     \==========================================================./
% C INPUT
% C       diclocal = total inorganic carbon (mol./m.^3)
% C             where 1 T = 1 metric ton = 1000 kg
% C       ta  = total alkalinity (eq./m.^3)
% C       pt  = inorganic phosphate (mol./.^3)
% C       sit = inorganic silicate (mol./.^3)
% C       t   = temperature (degrees C)
% C       s   = salinity (PSU)
% C OUTPUT
% C IMPORTANT: Some words about units - (JCO, 4./4./1999)
% c     - Models carry tracers in mol./m.^3 (on a per volume basis)
% c     - Conversely, this routine, which was written by observationalists
% c       (C. Sabine and R. Key), passes input arguments in umol./kg
% c       (i.e., on a per mass basis)
% c     - I have changed things slightly so that input arguments are in mol./m.^3,
% c     - Thus, all input concentrations (diclocal, ta, pt, and st) should be
% c       given in mol./m.^3; output arguments "co2star" and "dco2star"
% c       are likewise be in mol./m.^3.
% c
% c
% c NOW INCLUDES:
% c PRESSURE DEPENDENCE of K1, K2, solubility product of calcite
% c based on Takahashi, GEOSECS Atlantic Report, Vol. 1 (1981)
% c
% C--------------------------------------------------------------------------
% 
% C.....................................................................
% C OCMIP note:
% C Calculate all constants needed to convert between various measured
% C carbon species. References for each equation are noted in the code.
% C Once calculated, the constants are
% C stored and passed in the common block "const". The original version
% C of this code was based on the code by dickson in Version 2 of
% C "Handbook of Methods C for the Analysis of the Various Parameters of
% C the Carbon Dioxide System in Seawater", DOE, 1994 (SOP No. 3, p25-26).
% C....................................................................
% 
% c determine pressure (bar) from depth
% c 1 BAR at z=0m (atmos pressure)
% c use UPPER surface of cell so top layer pressure = 0 bar
% c for surface exchange coeffs
% 
% cmick..............................
% c        write(6,.*)'Klevel ',klevel

%      ff=0.d0;       ffp=0.d0;     ak0 = 0.d0;   ak1 = 0.d0;
%      ak2 = 0.d0;    akb = 0.d0;   ak1p = 0.d0;  ak2p = 0.d0;
%      ak3p = 0.d0;   aksi = 0.d0;  akw = 0.d0;   aks = 0.d0;
%      akf = 0.d0;    bt = 0.d0;    st = 0.d0;    ft = 0.d0;
%      Ksp_TP_Calc = 0.d0;

%       C terms used more than once
        tk = 273.15 + t ;
        tk100 = tk./100.0 ;
        tk1002=tk100.*tk100 ;
        invtk=1.0./tk ;
        dlogtk=log(tk) ;
        is=19.924.*s./(1000.-1.005.*s);
        is2=is.*is ;
        sqrtis=sqrt(is) ;
        s2=s.*s ;
        sqrts=sqrt(s);
        s15=s.^1.5 ;
        scl=s./1.80655 ;

% C------------------------------------------------------------------------
% C f = k0(1-pH2O).*correction term for non-ideality
% C Weiss &amp; Price (1980, Mar. Chem., 8, 347-359; Eq 13 with table 6 values)
           ff = exp(-162.8301 + 218.2968./tk100  + 90.9241.*log(tk100) - ...
               1.47696.*tk1002 + s .* (.025695 - .025225.*tk100 + 0.0049867.*tk1002));
           
% JML From CO2Sys, FugFac contains a pressure term BUT not sure of units
% and no salinity dependence! so not used....
           delta = (57.7 - 0.118.*tk);
           b = -1636.75 + 12.0408.*tk - 0.0327957.*tk.^2 + 3.16528.*0.00001.*tk.^3;
           ffp = exp((b + 2.*delta).*pressc./(83.143.*tk));
% C------------------------------------------------------------------------
% C K0 from Weiss 1974
           ak0 = exp(93.4517./tk100 - 60.2409 + 23.3585 .* log(tk100) + ...
               s .* (0.023517 - 0.023656.*tk100 + 0.0047036.*tk1002));
% C------------------------------------------------------------------------
% C k1 = [H][HCO3]./[H2CO3]
% C k2 = [H][CO3]./[HCO3]
% C Millero p.664 (1995) using Mehrbach et al. data on seawater scale 
           ak1=10.^(-1.*(3670.7.*invtk - 62.008 + 9.7944.*dlogtk - 0.0118 .* s + 0.000116.*s2));
           ak2=10.^(-1.*(1394.7.*invtk + 4.777 - 0.0184.*s + 0.000118.*s2));
% C NOW PRESSURE DEPENDENCE:
% c Following Takahashi (1981) GEOSECS report - quoting Culberson and 
% c Pytkowicz (1968)
% c pressc = pressure in bars
           ak1 = ak1 .* exp( (24.2-0.085 .* t ) .* ((pressc-1.0) ./ (83.143.*tk)) );
% c FIRST GO FOR K2: According to GEOSECS (1982) report
% c          ak2(i,j,bi,bj) = ak2(i,j,bi,bj).*
% c    &amp;             exp( (26.4-0.040.*t).*(pressc-1.0)./(83.143.*tk) )
% c SECOND GO FOR K2: corrected coeff according to CO2sys documentation
% c                   E. Lewis and D. Wallace (1998) ORNL./CDIAC-105
           ak2 = ak2 .* exp( (16.4-0.040.*t) .* ((pressc-1.0) ./ (83.143.*tk)) );
% C------------------------------------------------------------------------
% C kb = [H][BO2]./[HBO2]
% C Millero p.669 (1995) using data from dickson (1990)
           akb=exp((-8966.90 - 2890.53.*sqrts - 77.942.*s + 1.728.*s15 - 0.0996.*s2).*invtk ...
               + (148.0248 + 137.1942.*sqrts + 1.62142.*s) + ...
               (-24.4344 - 25.085.*sqrts - 0.2474.*s) .* dlogtk + 0.053105.*sqrts.*tk);
% C Mick and Karsten - Dec 04
% C ADDING pressure dependence based on Millero (1995), p675
% C with additional info from CO2sys documentation (E. Lewis and 
% C D. Wallace, 1998 - see endnotes for commentary on Millero, 95)
           bigR = 83.145  ;
           dv = -29.48 + 0.1622.*t + 2.608d-3.*t.*t ;
           dk = -2.84d-3 ;
           pfactor = - (dv./(bigR.*tk)).*pressc + (0.5.*dk./(bigR.*tk)).*pressc.*pressc ;
           akb = akb.*exp(pfactor) ;
% C------------------------------------------------------------------------
% C k1p = [H][H2PO4]./[H3PO4]
% C DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)
           ak1p = exp(-4576.752.*invtk + 115.525 - 18.453.*dlogtk +...
               (-106.736.*invtk + 0.69171).*sqrts + (-0.65643.*invtk - 0.01844).*s);
% C------------------------------------------------------------------------
% C k2p = [H][HPO4]./[H2PO4]
% C DOE(1994) eq 7.2.23 with footnote using data from Millero (1974))
           ak2p = exp(-8814.715.*invtk + 172.0883 - 27.927.*dlogtk + (-160.340.*invtk + 1.3566) .* ...
               sqrts + (0.37335.*invtk - 0.05778) .* s);
% C------------------------------------------------------------------------
% C k3p = [H][PO4]./[HPO4]
% C DOE(1994) eq 7.2.26 with footnote using data from Millero (1974)
           ak3p = exp(-3070.75.*invtk - 18.141 + (17.27039.*invtk + 2.81197) .*...
               sqrts + (-44.99486.*invtk - 0.09984) .* s);
% C------------------------------------------------------------------------
% C ksi = [H][SiO(OH)3]./[Si(OH)4]
% C Millero p.671 (1995) using data from Yao and Millero (1995)
           aksi = exp(-8904.2.*invtk + 117.385 - 19.334.*dlogtk + (-458.79.*invtk + 3.5913) .*...
               sqrtis + (188.74.*invtk - 1.5998) .* is + (-12.1652.*invtk + 0.07871) .* is2 +...
               log(1.0-0.001005.*s));
% C------------------------------------------------------------------------
% C kw = [H][OH]
% C Millero p.670 (1995) using composite data
           akw = exp(-13847.26.*invtk + 148.9652 - 23.6521.*dlogtk + ...
               (118.67.*invtk - 5.977 + 1.0495 .* dlogtk) .* sqrts - 0.01615 .* s);
% C------------------------------------------------------------------------
% C ks = [H][SO4]./[HSO4]
% C dickson (1990, J. chem. Thermodynamics 22, 113)
           aks=exp(-4276.1.*invtk + 141.328 - 23.093.*dlogtk + ...
               (-13856.*invtk + 324.57 - 47.986.*dlogtk).*sqrtis + (35474.*invtk -...
               771.54 + 114.723.*dlogtk).*is - 2698.*invtk.*is.^1.5 + 1776.*invtk.*is2 +...
               log(1.0 - 0.001005.*s));
% C------------------------------------------------------------------------
% C kf = [H][F]./[HF]
% C dickson and Riley (1979) -- change pH scale to total
           akf = exp(1590.2.*invtk - 12.641 + 1.525.*sqrtis + log(1.0 - 0.001005.*s) +...
               log(1.0 + (0.1400./96.062).*(scl)./aks));
% C------------------------------------------------------------------------
% C Calculate concentrations for borate, sulfate, and fluoride
% C Uppstrom (1974)
           bt = 0.000232 .* scl./10.811 ;
% C Morris & Riley (1966)
           st = 0.14 .* scl./96.062 ;
% C Riley (1965)
           ft = 0.000067 .* scl./18.9984 ;
% C------------------------------------------------------------------------
% C solubility product for calcite
% C
% c Following Takahashi (1982) GEOSECS handbook
% C NOT SURE THIS IS WORKING???
% C Ingle et al. (1973)
% c          Kspc = ( -34.452 - 39.866.*(s.*.*0.333333)
% c    &amp;                  + 110.21.*log(s) - 7.5752d-6 .* (tk.*.*2.0)
% c    &amp;                  ) .* 1.0d-7
% c with pressure dependence Culberson and Pytkowicz (1968)
% c          xvalue  =  (36-0.20.*t).*(pressc-1.0)./(83.143.*tk)
% c          Ksp_TP_Calc(i,j,bi,bj) = Kspc.*exp(xvalue)
% c
% c
% C Following Mucci (1983) - from Zeebe./Wolf-Gladrow equic.m
         tmpc1 = - 171.9065 - (0.077993.*tk) + (2839.319./tk) + (71.595.*log10(tk)) ;
         tmpc2 = +(-0.77712 + (0.0028426.*tk) + (178.34./tk) ).*sqrts ;
         tmpc3 = -(0.07711.*s) + (0.0041249.*s15) ;
         logKspc = tmpc1 + tmpc2 + tmpc3 ;
         Kspc = 10.0.^logKspc ;
% c        write(6,.*)i,j,k,tmpa1,tmpa2,tmpa3,logkspc,Kspc
% c with pressure dependence Culberson and Pytkowicz (1968)
% c        xvalue  =  (36.0-0.20.*t).*(pressc-1.0)./(83.143.*tk)
% c        Ksp_TP_Calc(i,j,bi,bj) = Kspc.*exp(xvalue)
% 
% c alternative pressure depdendence
% c following Millero (1995) but using info from Appendix A11 of
% c Zeebe and Wolf-Gladrow (2001) book
          dv = -48.76 - 0.5304.*t ;
          dk = -11.76d-3 - 0.3692.*t ;
          xvalue = - (dv./(bigR.*tk)).*pressc + (0.5.*dk./(bigR.*tk)).*pressc.*pressc ;
% c          Ksp_TP_Calc(i,j,bi,bj) = Kspc.*exp(xvalue)

% c alternative pressure dependence from Ingle (1975)
% 
%            zdum   = (pressc.*10.0d0 - 10.0d0)./10.0d0 ;
%            xvalue = ( (48.8d0 - 0.53d0.*t).*zdum + (-0.00588d0 + 0.0001845d0.*t).*zdum.*zdum) ./ ...
%                (188.93d0.*(t + 273.15d0));
           Kspcalc = Kspc.*10.^(xvalue);

% apparent solubility product of aragonite
%
%  Kspa = [Ca2+]T [CO32-]T
%
%  where $[]_T$ refers to the equilibrium total 
% (free + complexed) ion concentration.
%
%  Mucci 1983 mol./kg-soln

          tmpa1 = -171.945-0.077993.*tk+2903.293./tk+71.595.*log10(tk);
          tmpa2 = +(-0.068393+0.0017276.*tk+88.135./tk).*sqrt(s);
          tmpa3 = -0.10018.*s+0.0059415.*s15;
          log10Kspa = tmpa1 + tmpa2 + tmpa3;

          Kspa = 10..^(log10Kspa);

% Pressure dependence from Zeebe and Wolf-Gladrow (2001) book
          dv = -46. - 0.5304.*t ;
          dk = -11.76d-3 - 0.3692.*t ;
          xvalue = - (dv./(bigR.*tk)).*pressc + (0.5.*dk./(bigR.*tk)).*pressc.*pressc ;
          Ksparag = Kspa.*10.^(xvalue);
return