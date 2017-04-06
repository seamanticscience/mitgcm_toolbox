
C !INTERFACE: ==========================================================
       SUBROUTINE CALC_PCO2(iMax,jMax,kMax,mask,pressc,
     I                       theta,salt,dic,po4,si,alk,
     U                       ph,pco2oc,co2aq,hco3,co3,iter)

C !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE CALC_PCO2_APPROX                              |
C     *==========================================================*
C     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     C  New efficient pCO2 solver, Mick Follows         CC
C     C                             Taka Ito             CC
C     C                             Stephanie Dutkiewicz CC
C     C  20 April 2003                                   CC
C     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      Apr 201 fix vapour bug (following Bennington)
C      Oct 201 add CO3 extimation and pass out

C !USES: ===============================================================
      IMPLICIT NONE
C     == Routine arguments ==
C       diclocal = total inorganic carbon (mol/m^3)
C             where 1 T = 1 metric ton = 1000 kg
C       ta  = total alkalinity (eq/m^3)
C       pt  = inorganic phosphate (mol/^3)
C       sit = inorganic silicate (mol/^3)
C       t   = temperature (degrees C)
C       s   = salinity (PSU)
        integer*4 iMax,jMax,kMax,i,j,k
        real*8  t, s, pt, sit, ta, ct
        real*8  ff, bt, st, ft
        real*8  k1, k2
        real*8  k1p, k2p, k3p
        real*8  ks, kb, kw, ksi, kf
        real*8  k0, fugf, Ksp_TP_Calc

C     == Local variables ==
        real*8  phguess
        real*8  cag
        real*8  bohg
        real*8  hguess
        real*8  stuff
        real*8  gamm
        real*8  hnew
        real*8  co2s, hco3loc, co3loc
        real*8  h3po4g, h2po4g, hpo4g, po4g
        real*8  siooh3g
        real*8  fco2

        real*8 mask(iMax,jMax,kMax)
        real*8 pressc(iMax,jMax,kMax)
        real*8 theta(iMax,jMax,kMax)
        real*8 salt(iMax,jMax,kMax)
        real*8 dic (iMax,jMax,kMax)        
        real*8 po4 (iMax,jMax,kMax)
        real*8 si  (iMax,jMax,kMax)
        real*8 alk (iMax,jMax,kMax)
        real*8 ph  (iMax,jMax,kMax)
        real*8 pco2oc(iMax,jMax,kMax)
        real*8 co2aq(iMax,jMax,kMax)        
        real*8 hco3 (iMax,jMax,kMax)        
        real*8 co3(iMax,jMax,kMax)        
        real*8 iter(iMax,jMax,kMax)

        real*8 permil
        parameter(permil = 1.D0/1024.5D0)

        do i=1,iMax
         do j=1,jMax
          do k=1,kMax
           if (mask(i,j,k).gt.0.D0) then
            t  = theta(i,j,k)
            s  = salt(i,j,k) 
c ---------------------------------------------------------------------
C Change units from the input of mol/m^3 -> mol/kg:
c (1 mol/m^3)  x (1 m^3/1024.5 kg)
c where the ocean mean surface density is 1024.5 kg/m^3
c Note: mol/kg are actually what the body of this routine uses
c for calculations.  Units are reconverted back to mol/m^3 at the
c end of this routine.
c To convert input in mol/m^3 -> mol/kg
            pt  = po4(i,j,k)*permil
            sit = si (i,j,k)*permil
            ta  = alk(i,j,k)*permil
            ct  = dic(i,j,k)*permil
            
c Call carbon_coeffs_press_dep to load dissociation constants
            CALL CARBON_COEFFS_PRESSURE_DEP(
     I					 mask(i,j,k),t,s,pressc(i,j,k),
     O					 fugf,ff,k0,k1,k2,kb,k1p,k2p,k3p,
     O                   ksi,kw,ks,kf,bt,st,ft,Ksp_TP_Calc)        

c ---------------------------------------------------------------------
c set first guess and brackets for [H+] solvers
C pH guess as average preindustrial
            ph(i,j,k)=8.179D0
C            hnew = 10.0D0**(-ph(i,j,k))
            hnew   = 0.D0
            phguess= 8.179D0
            hguess = 10.0D0**(-phguess)
            iter(i,j,k)=0.D0
C             WRITE(standardMessageUnit,*) 
C     &      'CALC_CSAT_APPROX: ph guess'

C Iterate from a cold ph start, 3e-11 is the suggestion from Mick's paper
            do while ( abs(hguess - hnew)>1.0D-11)  
             phguess = pH(i,j,k)
cmick - new approx method
cmick - make estimate of htotal (hydrogen ion conc) using
cmick   appromate estimate of CA, carbonate alkalinity
             hguess = 10.0**(-phguess)
cmick - first estimate borate contribution using guess for [H+]
             bohg = bt*kb/(hguess+kb)

cmick - first estimate of contribution from phosphate
cmick based on Dickson and Goyet
             stuff = hguess*hguess*hguess
     &           + (k1p*hguess*hguess)
     &           + (k1p*k2p*hguess)
     &           + (k1p*k2p*k3p)
             h3po4g = (pt*hguess*hguess*hguess) / stuff
             h2po4g = (pt*k1p*hguess*hguess) / stuff
             hpo4g  = (pt*k1p*k2p*hguess) / stuff
             po4g   = (pt*k1p*k2p*k3p) / stuff

cmick - estimate contribution from silicate
cmick based on Dickson and Goyet
             siooh3g = sit*ksi / (ksi + hguess)

cmick - now estimate carbonate alkalinity
             cag = ta - bohg - (kw/hguess) + hguess
     &           - hpo4g - 2.0D0*po4g + h3po4g
     &           - siooh3g

cmick - now evaluate better guess of hydrogen ion conc
cmick   htotal = [H+], hydrogen ion conc
             gamm  = ct/cag
             stuff = (1.0D0-gamm)*(1.0D0-gamm)*k1*k1
     &          - 4.0D0*k1*k2*(1.0D0-2.0D0*gamm)
             hnew  = 0.5D0*( (gamm-1.0D0)*k1 + sqrt(stuff) )

cmick - return update pH to main routine
             ph(i,j,k) = -log10(hnew)

C end the iterating loop
             iter(i,j,k)=iter(i,j,k)+1.D0    
             
             if (iter(i,j,k).eq.50.D0) then
              exit
             endif          
            end do
            
Cmick - determine [CO2*]
            co2s = ct/
     &   (1.0D0 + (k1/hnew) + (k1*k2/(hnew*hnew)))            
C Evaluate HCO3- for c13 air-sea exchange (CO32- also needed)
            hco3loc = k1*ct/
     &         (hnew*hnew + k1*hnew + k1*k2)
c NOW EVALUATE CO32-, carbonate ion concentration
c used in determination of calcite compensation depth
c Karsten Friis & Mick - Sep 2004
            co3loc = k1*k2*ct/
     &         (hnew*hnew + k1*hnew + k1*k2)
c ---------------------------------------------------------------
c surface pCO2 (following Dickson and Goyet, DOE...)
c bug fix by Bennington
C            fco2 = co2s/k0
C            pco2oc(i,j,k) = fco2/fugf

            pco2oc(i,j,k) = co2s/ff
C ----------------------------------------------------------------
c Reconvert from mol/kg -> mol/m^3
            co2aq(i,j,k)=co2s/permil
            hco3 (i,j,k)=hco3loc/permil
            co3  (i,j,k)=co3loc/permil
           else
            pco2oc (i,j,k)=0.D0
            ph   (i,j,k)=0.D0
            co2aq(i,j,k)=0.D0
            hco3 (i,j,k)=0.D0
            co3  (i,j,k)=0.D0
            iter (i,j,k)=0.D0
           endif
          end do
         end do
        end do
         
        RETURN
        END
        
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

CBOP
C !ROUTINE: CALC_CSAT_APPROX

C !INTERFACE: ==========================================================
       SUBROUTINE CALC_CSAT(iMax,jMax,kMax,mask,pressc,
     I                       theta,salt,pco2atm,po4,si,alk,
     U                       ph,csat,co3,iter )

C !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE CALC_CSAT_APPROX                              |
C     *==========================================================*
C      CALCULATE Csat, EQUILIBRIUM DIC FOR GIVEN pCO2, T, Alk etc
C !USES: ===============================================================
      IMPLICIT NONE

C     == Routine arguments ==
C       dic = total inorganic carbon (mol/m^3)
C             where 1 T = 1 metric ton = 1000 kg
C       ta  = total alkalinity (eq/m^3)
C       pt  = inorganic phosphate (mol/^3)
C       sit = inorganic silicate (mol/^3)
C       t   = temperature (degrees C)
C       s   = salinity (PSU)
        integer*4 iMax,jMax,kMax,i,j,k
        real*8 mask(iMax,jMax,kMax)
        real*8 theta(iMax,jMax,kMax)
        real*8 salt(iMax,jMax,kMax)
        real*8 po4 (iMax,jMax,kMax)
        real*8 si  (iMax,jMax,kMax)
        real*8 alk (iMax,jMax,kMax)
        real*8 csat(iMax,jMax,kMax)
        real*8 ph  (iMax,jMax,kMax)
        real*8 co3 (iMax,jMax,kMax)        
        real*8 pressc(iMax,jMax,kMax)
        real*8 iter(iMax,jMax,kMax)
        
        real*8 pco2atm, t, s, pt, ta, sit
        real*8 ff, bt, st, ft
        real*8 k0, k1, k2
        real*8 k1p, k2p, k3p
        real*8 ks, kb, kw, ksi, kf
        real*8 fugf, Ksp_TP_Calc
        real*8 permil
        parameter(permil = 1.D0/1024.5D0)

C     == Local variables ==
        real*8 phguess
        real*8 cag
        real*8 bohg
        real*8 hguess
        real*8 denom
        real*8 stuff, stuff2, stuff3
        real*8 hnew, csatloc, co3loc
        real*8 h3po4g, h2po4g, hpo4g, po4g
        real*8 siooh3g

        do i=1,iMax
         do j=1,jMax
          do k=1,kMax
           if (mask(i,j,k).gt.0.D0) then
            t  = theta(i,j,k)
            s  = salt(i,j,k) 

c ---------------------------------------------------------------------
C Change units from the input of mol/m^3 -> mol/kg:
c (1 mol/m^3)  x (1 m^3/1024.5 kg)
c where the ocean mean surface density is 1024.5 kg/m^3
c Note: mol/kg are actually what the body of this routine uses
c for calculations.  Units are reconverted back to mol/m^3 at the
c end of this routine.
c To convert input in mol/m^3 -> mol/kg

            pt = po4(i,j,k)*permil
            sit= si (i,j,k)*permil
            ta = alk(i,j,k)*permil
                   
c Call carbon_coeffs_press_dep to load dissociation constants
            CALL CARBON_COEFFS_PRESSURE_DEP(
     I					 mask(i,j,k),t,s,pressc(i,j,k),
     O					 fugf,ff,k0,k1,k2,kb,k1p,k2p,k3p,
     O                   ksi,kw,ks,kf,bt,st,ft,Ksp_TP_Calc)        
c ---------------------------------------------------------------------
c set first guess and brackets for [H+] solvers
C pH guess as average preindustrial
            ph(i,j,k)=8.179D0
C            hnew = 10.0D0**(-ph(i,j,k))
            hnew   = 0.D0
            phguess= 8.179D0
            hguess = 10.0D0**(-phguess)
            iter(i,j,k)=0.D0
C             WRITE(standardMessageUnit,*) 
C     &      'CALC_CSAT_APPROX: ph guess'

C Iterate from a cold ph start, 3e-11 is the suggestion from Mick's paper
            do while ( abs(hguess - hnew)>1.0D-11)  
             phguess = ph(i,j,k)
C     - new approx method
C     - make estimate of htotal (hydrogen ion conc) using
C       appromate estimate of CA, carbonate alkalinity
             hguess = 10.0D0**(-phguess)
        
C     - first estimate borate contribution using guess for [H+]
             bohg = (bt*kb)/(hguess+kb)

C     - first estimate of contribution from phosphate
C     based on Dickson and Goyet
             denom  = (hguess*hguess*hguess)
     &           + (k1p*hguess*hguess)
     &           + (k1p*k2p*hguess)
     &           + (k1p*k2p*k3p)
             h3po4g = (pt*hguess*hguess*hguess) / denom
             h2po4g = (pt*k1p*hguess*hguess) / denom
             hpo4g  = (pt*k1p*k2p*hguess) / denom
             po4g   = (pt*k1p*k2p*k3p) / denom

C     - estimate contribution from silicate
C     based on Dickson and Goyet
             siooh3g = (sit*ksi) / (ksi + hguess)

C     - now estimate carbonate alkalinity
             cag = ta - bohg - (kw/hguess) + hguess
     &           - hpo4g - 2.0D0*po4g + h3po4g
     &           - siooh3g

C     - estimate hydrogen ion conc
C        stuff = (k1*fugf*k0*pco2atm)/cag
             stuff = (k1*ff*pco2atm)/cag  
             stuff2 = stuff*(stuff + 8.0D0*k2)  
             stuff3 = sqrt(stuff2)  
             hnew = 0.5D0*(stuff + stuff3)  
C            hguess = hnew  
            
C calc final pH
             ph(i,j,k) = -log10(hnew)
            
C end the iterating loop
             iter(i,j,k)=iter(i,j,k)+1.D0      
                     
             if (iter(i,j,k).eq.50.D0) then
              exit
             endif  
            end do
            
C evaluate csat, equilibrium DIC concentration
            csatloc = (pco2atm*ff)
     &   *(1.0D0 + (k1/hnew) + ((k1*k2)/(hnew*hnew)))

C Output csat in mol/m^3
            csat(i,j,k)=csatloc/permil 
            
c NOW EVALUATE CO32-, carbonate ion concentration
c used in determination of calcite compensation depth
c Karsten Friis & Mick - Sep 2004
            co3loc = k1*k2*csatloc /
     &         (hnew*hnew + k1*hnew + k1*k2)
     
            co3(i,j,k)=co3loc/permil
           else
            csat(i,j,k)=0.D0
            ph  (i,j,k)=0.D0
            co3 (i,j,k)=0.D0
            iter(i,j,k)=0.D0
           endif
          end do
         end do
        end do
       
        RETURN
        END
        
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CCBOP
CC !ROUTINE: CARBON_COEFFS
C
CC !INTERFACE: ==========================================================
C      SUBROUTINE CARBON_COEFFS(
C     I                   ttemp,stemp,
C     I                   bi,bj,iMin,iMax,jMin,jMax,myThid)
C
CC !DESCRIPTION:
CC     *==========================================================*
CC     | SUBROUTINE CARBON_COEFFS                                 |
CC     | determine coefficients for surface carbon chemistry      |
CC     | adapted from OCMIP2:  SUBROUTINE CO2CALC                 |
CC     | mick follows, oct 1999                                   |
CC     | minor changes to tidy, swd aug 2002                      |
CC     *==========================================================*
CC INPUT
CC       diclocal = total inorganic carbon (mol/m^3)
CC             where 1 T = 1 metric ton = 1000 kg
CC       ta  = total alkalinity (eq/m^3)
CC       pt  = inorganic phosphate (mol/^3)
CC       sit = inorganic silicate (mol/^3)
CC       t   = temperature (degrees C)
CC       s   = salinity (PSU)
CC OUTPUT
CC IMPORTANT: Some words about units - (JCO, 4/4/1999)
CC     - Models carry tracers in mol/m^3 (on a per volume basis)
CC     - Conversely, this routine, which was written by observationalists
CC       (C. Sabine and R. Key), passes input arguments in umol/kg
CC       (i.e., on a per mass basis)
CC     - I have changed things slightly so that input arguments are in mol/m^3,
CC     - Thus, all input concentrations (diclocal, ta, pt, and st) should be
CC       given in mol/m^3  output arguments "co2star" and "dco2star"
CC       are likewise be in mol/m^3.
CC
CC Apr 201 fix vapour bug (following Bennington)
CC--------------------------------------------------------------------------
C
CC !USES: ===============================================================
C        IMPLICIT NONE
CC     == GLobal variables ==
C#include "SIZE.h"
C#include "DYNVARS.h"
C#include "EEPARAMS.h"
C#include "PARAMS.h"
C#include "GRID.h"
C#include "FFIELDS.h"
C#include "DIC_VARS.h"
CC     == Routine arguments ==
CC ttemp and stemp are local theta and salt arrays
CC dont really need to pass T and S in, could use theta, salt in
CC common block in DYNVARS.h, but this way keeps subroutine more
CC general
C        real*8  ttemp(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
C        real*8  stemp(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
C        integer*4 bi,bj,iMin,iMax,jMin,jMax
C        integer*4 myThid
CCEOP
C
CC LOCAL VARIABLES
C        real*8  t
C        real*8  s
C        real*8  tk
C        real*8  tk100
C        real*8  tk1002
C        real*8  dlogtk
C        real*8  sqrtis
C        real*8  sqrts
C        real*8  s15
C        real*8  scl
C        real*8  s2
C        real*8  invtk
C        real*8  is
C        real*8  is2
Cc add Bennington
C        real*8  P1atm
C        real*8  Rgas
C        real*8  RT
C        real*8  delta
C        real*8  B1
C        real*8  B
C        integer*4 i
C        integer*4 j
C
CC.....................................................................
CC OCMIP note:
CC Calculate all constants needed to convert between various measured
CC carbon species. References for each equation are noted in the code.
CC Once calculated, the constants are
CC stored and passed in the common block "const". The original version
CC of this code was based on the code by dickson in Version 2 of
CC  Handbook of Methods C for the Analysis of the Various Parameters of
CC the Carbon Dioxide System in Seawater , DOE, 1994 (SOP No. 3, p25-26).
CC....................................................................
C
C        do i=iMin,iMax
C         do j=jMin,jMax
C          if (hFacC(i,j,1,bi,bj).gt.0.D0) then
C           t = ttemp(i,j)
C           s = stemp(i,j)
CC terms used more than once
C           tk = 273.15D0 + t
C           tk100 = tk/100.D0
C           tk1002=tk100*tk100
C           invtk=1.0D0/tk
C           dlogtk=log(tk)
C           is=19.924D0*s/(1000.D0-1.005D0*s)
C           is2=is*is
C           sqrtis=sqrt(is)
C           s2=s*s
C           sqrts=sqrt(s)
C           s15=s**1.5D0
C           scl=s/1.80655D0
CC -----------------------------------------------------------------------
CC added by Val Bennington Nov 2010
CC Fugacity Factor needed for non-ideality in ocean
CC ff used for atmospheric correction for water vapor and pressure
CC Weiss (1974) Marine Chemistry
C           P1atm = 1.01325D0 ! bars
C           Rgas = 83.1451D0 ! bar*cm3/(mol*K)
C           RT = Rgas*tk
C           delta = (57.7D0 - 0.118D0*tk)
C           B1 = -1636.75D0 + 12.0408D0*tk - 0.0327957D0*tk*tk
C           B = B1 + 3.16528D0*tk*tk*tk*(0.00001D0)
C           fugf(i,j,bi,bj) = exp( (B+2.D0*delta) * P1atm / RT)
CC------------------------------------------------------------------------
CC f = k0(1-pH2O)*correction term for non-ideality
CC Weiss & Price (1980, Mar. Chem., 8, 347-359  Eq 13 with table 6 values)
C           ff(i,j,bi,bj) = exp(-162.8301D0 + 218.2968D0/tk100  +
C     &          90.9241D0*log(tk100) - 1.47696D0*tk1002 +
C     &          s * (.025695D0 - .025225D0*tk100 +
C     &          0.0049867D0*tk1002))
CC------------------------------------------------------------------------
CC K0 from Weiss 1974
C           ak0(i,j,bi,bj) = exp(93.4517D0/tk100 - 60.2409D0 +
C     &        23.3585D0 * log(tk100) +
C     &        s * (0.023517D0 - 0.023656D0*tk100 +
C     &        0.0047036D0*tk1002))
CC------------------------------------------------------------------------
CC k1 = [H][HCO3]/[H2CO3]
CC k2 = [H][CO3]/[HCO3]
CC Millero p.664 (1995) using Mehrbach et al. data on seawater scale
C           ak1(i,j,bi,bj)=10.**(-1.D0*(3670.7D0*invtk -
C     &          62.008D0 + 9.7944D0*dlogtk -
C     &          0.0118D0 * s + 0.000116D0*s2))
C           ak2(i,j,bi,bj)=10.**(-1.D0*(1394.7D0*invtk+ 4.777D0-
C     &          0.0184D0*s + 0.000118D0*s2))
CC------------------------------------------------------------------------
CC kb = [H][BO2]/[HBO2]
CC Millero p.669 (1995) using data from dickson (1990)
C           akb(i,j,bi,bj)=exp((-8966.90D0- 2890.53D0*sqrts -
C     &          77.942D0*s + 1.728D0*s15 - 0.0996D0*s2)*invtk +
C     &          (148.0248D0 + 137.1942D0*sqrts + 1.62142D0*s) +
C     &          (-24.4344D0 - 25.085D0*sqrts - 0.2474D0*s) *
C     &          dlogtk + 0.053105D0*sqrts*tk)
CC------------------------------------------------------------------------
CC k1p = [H][H2PO4]/[H3PO4]
CC DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)
C           ak1p(i,j,bi,bj) = exp(-4576.752D0*invtk + 115.525D0 -
C     &          18.453D0*dlogtk +
C     &          (-106.736D0*invtk + 0.69171D0)*sqrts +
C     &          (-0.65643D0*invtk - 0.01844D0)*s)
CC------------------------------------------------------------------------
CC k2p = [H][HPO4]/[H2PO4]
CC DOE(1994) eq 7.2.23 with footnote using data from Millero (1974))
C           ak2p(i,j,bi,bj) = exp(-8814.715D0*invtk + 172.0883D0 -
C     &          27.927D0*dlogtk +
C     &          (-160.340D0*invtk + 1.3566D0) * sqrts +
C     &          (0.37335D0*invtk - 0.05778D0) * s)
CC------------------------------------------------------------------------
CC k3p = [H][PO4]/[HPO4]
CC DOE(1994) eq 7.2.26 with footnote using data from Millero (1974)
C           ak3p(i,j,bi,bj) = exp(-3070.75D0*invtk - 18.141D0 +
C     &          (17.27039D0*invtk + 2.81197D0) *
C     &          sqrts + (-44.99486D0*invtk - 0.09984D0) * s)
CC------------------------------------------------------------------------
CC ksi = [H][SiO(OH)3]/[Si(OH)4]
CC Millero p.671 (1995) using data from Yao and Millero (1995)
C           aksi(i,j,bi,bj) = exp(-8904.2D0*invtk + 117.385D0 -
C     &          19.334D0*dlogtk +
C     &          (-458.79D0*invtk + 3.5913D0) * sqrtis +
C     &          (188.74D0*invtk - 1.5998D0) * is +
C     &          (-12.1652D0*invtk + 0.07871D0) * is2 +
C     &          log(1.0D0-0.001005D0*s))
CC------------------------------------------------------------------------
CC kw = [H][OH]
CC Millero p.670 (1995) using composite data
C           akw(i,j,bi,bj) = exp(-13847.26D0*invtk + 148.9652D0 -
C     &          23.6521D0*dlogtk +
C     &          (118.67D0*invtk - 5.977D0 + 1.0495D0 * dlogtk)
C     &          * sqrts - 0.01615D0 * s)
CC------------------------------------------------------------------------
CC ks = [H][SO4]/[HSO4]
CC dickson (1990, J. chem. Thermodynamics 22, 113)
C           aks(i,j,bi,bj)=exp(-4276.1D0*invtk + 141.328D0 -
C     &          23.093D0*dlogtk +
C     &   (-13856.D0*invtk + 324.57D0 - 47.986D0*dlogtk)*sqrtis+
C     &   (35474.D0*invtk - 771.54D0 + 114.723D0*dlogtk)*is -
C     &          2698.D0*invtk*is**1.5D0 + 1776.D0*invtk*is2 +
C     &          log(1.0D0 - 0.001005D0*s))
CC------------------------------------------------------------------------
CC kf = [H][F]/[HF]
CC dickson and Riley (1979) -- change pH scale to total
C           akf(i,j,bi,bj)=exp(1590.2D0*invtk - 12.641D0 +
C     &   1.525D0*sqrtis + log(1.0D0 - 0.001005D0*s) +
C     &   log(1.0D0 + (0.1400D0/96.062D0)*(scl)/aks(i,j,bi,bj)))
CC------------------------------------------------------------------------
CC Calculate concentrations for borate, sulfate, and fluoride
CC Uppstrom (1974)
C           bt(i,j,bi,bj) = 0.000232D0 * scl/10.811D0
CC Morris & Riley (1966)
C           st(i,j,bi,bj) = 0.14D0 * scl/96.062D0
CC Riley (1965)
C           ft(i,j,bi,bj) = 0.000067D0 * scl/18.9984D0
CC------------------------------------------------------------------------
C         else
Cc add Bennington
C            fugf(i,j,bi,bj)=0.D0
C            ff(i,j,bi,bj)=0.D0
C            ak0(i,j,bi,bj)= 0.D0
C            ak1(i,j,bi,bj)= 0.D0
C            ak2(i,j,bi,bj)= 0.D0
C            akb(i,j,bi,bj)= 0.D0
C            ak1p(i,j,bi,bj) = 0.D0
C            ak2p(i,j,bi,bj) = 0.D0
C            ak3p(i,j,bi,bj) = 0.D0
C            aksi(i,j,bi,bj) = 0.D0
C            akw(i,j,bi,bj) = 0.D0
C            aks(i,j,bi,bj)= 0.D0
C            akf(i,j,bi,bj)= 0.D0
C            bt(i,j,bi,bj) = 0.D0
C            st(i,j,bi,bj) = 0.D0
C            ft(i,j,bi,bj) = 0.D0
C         endif
C         end do
C        end do
C
C        RETURN
C        END
C
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

CBOP
C !ROUTINE: CARBON_COEFFS_PRESSURE_DEP

C !INTERFACE: ==========================================================
      SUBROUTINE CARBON_COEFFS_PRESSURE_DEP(
     I					 mask,t,s,pressc,
     O					 fugf,ff,ak0,ak1,ak2,akb,ak1p,ak2p,ak3p,
     O                   aksi,akw,aks,akf,bt,st,ft,Ksp_TP_Calc)

C !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE CARBON_COEFFS                                 |
C     | determine coefficients for surface carbon chemistry      |
C     | adapted from OCMIP2:  SUBROUTINE CO2CALC                 |
C     | mick follows, oct 1999                                   |
C     | minor changes to tidy, swd aug 2002                      |
C     | MODIFIED FOR PRESSURE DEPENDENCE                         |
C     | Karsten Friis and Mick Follows 2004                      |
C     *==========================================================*
C INPUT
C       diclocal = total inorganic carbon (mol/m^3)
C             where 1 T = 1 metric ton = 1000 kg
C       ta  = total alkalinity (eq/m^3)
C       pt  = inorganic phosphate (mol/^3)
C       sit = inorganic silicate (mol/^3)
C       t   = temperature (degrees C)
C       s   = salinity (PSU)
C OUTPUT
C IMPORTANT: Some words about units - (JCO, 4/4/1999)
C     - Models carry tracers in mol/m^3 (on a per volume basis)
C     - Conversely, this routine, which was written by observationalists
C       (C. Sabine and R. Key), passes input arguments in umol/kg
C       (i.e., on a per mass basis)
C     - I have changed things slightly so that input arguments are in mol/m^3,
C     - Thus, all input concentrations (diclocal, ta, pt, and st) should be
C       given in mol/m^3  output arguments "co2star" and "dco2star"
C       are likewise be in mol/m^3.
C
C Apr 201 fix vapour bug (following Bennington)
C
C NOW INCLUDES:
C PRESSURE DEPENDENCE of K1, K2, solubility product of calcite
C based on Takahashi, GEOSECS Atlantic Report, Vol. 1 (1981)
C
C--------------------------------------------------------------------------

C !USES: ===============================================================
        IMPLICIT NONE
C LOCAL VARIABLES
        real*8  t, s, mask, pressc
        real*8  tk, tk100, tk1002, dlogtk, invtk
        real*8  sqrtis, sqrts, s15, scl, s2, is, is2
        real*8  Ksp_T_Calc, Ksp_TP_Calc
        real*8  xvalue, zdum, tmpa1, tmpa2, tmpa3
        real*8  logKspc
        real*8  dv, dk, pfactor, bigR
c add Bennington
        real*8  P1atm, Rgas, RT, delta, B1, B
        
        real*8 ff, fugf, bt, st, ft
        real*8 ak0, ak1, ak2
        real*8 ak1p, ak2p, ak3p
        real*8 aks, akb, akw, aksi, akf
C.....................................................................
C OCMIP note:
C Calculate all constants needed to convert between various measured
C carbon species. References for each equation are noted in the code.
C Once calculated, the constants are
C stored and passed in the common block "const". The original version
C of this code was based on the code by dickson in Version 2 of
C  Handbook of Methods C for the Analysis of the Various Parameters of
C the Carbon Dioxide System in Seawater , DOE, 1994 (SOP No. 3, p25-26).
C....................................................................

          if (mask.gt.0.D0) then

C terms used more than once
           tk = 273.15D0 + t
           tk100 = tk/100.D0
           tk1002=tk100*tk100
           invtk=1.0D0/tk
           dlogtk=log(tk)
           is=19.924D0*s/(1000.D0-1.005D0*s)
           is2=is*is
           sqrtis=sqrt(is)
           s2=s*s
           sqrts=sqrt(s)
           s15=s**1.5D0
           scl=s/1.80655D0
           
C -----------------------------------------------------------------------
C added by Val Bennington Nov 2010
C Fugacity Factor needed for non-ideality in ocean
C ff used for atmospheric correction for water vapor and pressure
C Weiss (1974) Marine Chemistry
C           P1atm = 1.01325D0 ! bars
           P1atm = pressc
           Rgas = 83.1451D0 ! bar*cm3/(mol*K)
           RT = Rgas*tk
           delta = (57.7D0 - 0.118D0*tk)
           B1 = -1636.75D0 + 12.0408D0*tk - 0.0327957D0*tk*tk
           B = B1 + 3.16528D0*tk*tk*tk*(0.00001D0)
           fugf = exp( (B+2.D0*delta) * P1atm / RT)
C------------------------------------------------------------------------
C f = k0(1-pH2O)*correction term for non-ideality
C Weiss & Price (1980, Mar. Chem., 8, 347-359  Eq 13 with table 6 values)
           ff = exp(-162.8301D0 + 218.2968D0/tk100  +
     &          90.9241D0*log(tk100) - 1.47696D0*tk1002 +
     &          s * (.025695D0 - .025225D0*tk100 +
     &          0.0049867D0*tk1002))
C------------------------------------------------------------------------
C K0 from Weiss 1974
           ak0 = exp(93.4517D0/tk100 - 60.2409D0 +
     &        23.3585D0 * log(tk100) +
     &        s * (0.023517D0 - 0.023656D0*tk100 +
     &        0.0047036D0*tk1002))
C------------------------------------------------------------------------
C k1 = [H][HCO3]/[H2CO3]
C k2 = [H][CO3]/[HCO3]
C Millero p.664 (1995) using Mehrbach et al. data on seawater scale
           ak1=10**(-1D0*(3670.7D0*invtk -
     &          62.008D0 + 9.7944D0*dlogtk -
     &          0.0118D0 * s + 0.000116D0*s2))
           ak2=10**(-1*(1394.7D0*invtk + 4.777D0 -
     &          0.0184D0*s + 0.000118D0*s2))
C NOW PRESSURE DEPENDENCE:
c Following Takahashi (1981) GEOSECS report - quoting Culberson and
c Pytkowicz (1968)
c pressc = pressure in bars
           ak1 = ak1*
     &         exp( (24.2D0-0.085D0*t)*(pressc-1.0D0)/(83.143D0*tk) )
c FIRST GO FOR K2: According to GEOSECS (1982) report
c          ak2 = ak2*
c    &             exp( (26.4-0.040*t)*(pressc-1.0)/(83.143*tk) )
c SECOND GO FOR K2: corrected coeff according to CO2sys documentation
c                   E. Lewis and D. Wallace (1998) ORNL/CDIAC-105
           ak2 = ak2*
     &         exp( (16.4D0-0.040D0*t)*(pressc-1.0D0)/(83.143D0*tk) )
C------------------------------------------------------------------------
C kb = [H][BO2]/[HBO2]
C Millero p.669 (1995) using data from dickson (1990)
           akb=exp((-8966.90D0- 2890.53D0*sqrts -
     &          77.942D0*s + 1.728D0*s15 - 0.0996D0*s2)*invtk +
     &          (148.0248D0 + 137.1942D0*sqrts + 1.62142D0*s) +
     &          (-24.4344D0 - 25.085D0*sqrts - 0.2474D0*s) *
     &          dlogtk + 0.053105D0*sqrts*tk)

C Mick and Karsten - Dec 04
C ADDING pressure dependence based on Millero (1995), p675
C with additional info from CO2sys documentation (E. Lewis and
C D. Wallace, 1998 - see endnotes for commentary on Millero, 95)
           bigR = 83.145D0
           dv = -29.48D0 + 0.1622D0*t + 2.608D-3*t*t
           dk = -2.84D-3
           pfactor = - (dv/(bigR*tk))*pressc
     &               + (0.5D0*dk/(bigR*tk))*pressc*pressc
           akb = akb*exp(pfactor)
C------------------------------------------------------------------------
C k1p = [H][H2PO4]/[H3PO4]
C DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)
           ak1p = exp(-4576.752D0*invtk + 115.525D0 -
     &          18.453D0*dlogtk +
     &          (-106.736D0*invtk + 0.69171D0)*sqrts +
     &          (-0.65643D0*invtk - 0.01844D0)*s)
C------------------------------------------------------------------------
C k2p = [H][HPO4]/[H2PO4]
C DOE(1994) eq 7.2.23 with footnote using data from Millero (1974))
           ak2p = exp(-8814.715D0*invtk + 172.0883D0 -
     &          27.927D0*dlogtk +
     &          (-160.34D00*invtk + 1.3566D0) * sqrts +
     &          (0.37335D0*invtk - 0.05778D0) * s)
C------------------------------------------------------------------------
C k3p = [H][PO4]/[HPO4]
C DOE(1994) eq 7.2.26 with footnote using data from Millero (1974)
           ak3p = exp(-3070.75D0*invtk - 18.141D0 +
     &          (17.27039D0*invtk + 2.81197D0) *
     &          sqrts + (-44.99486D0*invtk - 0.09984D0) * s)
C------------------------------------------------------------------------
C ksi = [H][SiO(OH)3]/[Si(OH)4]
C Millero p.671 (1995) using data from Yao and Millero (1995)
           aksi = exp(-8904.2D0*invtk + 117.385D0 -
     &          19.334D0*dlogtk +
     &          (-458.79D0*invtk + 3.5913D0) * sqrtis +
     &          (188.74D0*invtk - 1.5998D0) * is +
     &          (-12.1652D0*invtk + 0.07871D0) * is2 +
     &          log(1.0D0-0.001005D0*s))
C------------------------------------------------------------------------
C kw = [H][OH]
C Millero p.670 (1995) using composite data
           akw = exp(-13847.26D0*invtk + 148.9652D0 -
     &          23.6521D0*dlogtk +
     &          (118.67D0*invtk - 5.977D0 + 1.0495D0 * dlogtk) *
     &          sqrts - 0.01615D0 * s)
C------------------------------------------------------------------------
C ks = [H][SO4]/[HSO4]
C dickson (1990, J. chem. Thermodynamics 22, 113)
           aks=exp(-4276.1D0*invtk + 141.328D0 -
     &          23.093D0*dlogtk +
     &          (-13856D0*invtk + 324.57D0 - 47.986D0*dlogtk)*sqrtis +
     &          (35474D0*invtk - 771.54D0 + 114.723D0*dlogtk)*is -
     &          2698D0*invtk*is**1.5D0 + 1776D0*invtk*is2 +
     &          log(1.0D0 - 0.001005D0*s))
C------------------------------------------------------------------------
C kf = [H][F]/[HF]
C dickson and Riley (1979) -- change pH scale to total
        akf=exp(1590.2D0*invtk - 12.641D0 + 1.525D0*sqrtis +
     &          log(1.0D0 - 0.001005D0*s) +
     &          log(1.0D0 + (0.1400D0/96.062D0)*(scl)/aks))
C------------------------------------------------------------------------
C Calculate concentrations for borate, sulfate, and fluoride
C Uppstrom (1974)
           bt = 0.000232D0 * scl/10.811D0
C Morris & Riley (1966)
           st = 0.14D0 * scl/96.062D0
C Riley (1965)
           ft = 0.000067D0 * scl/18.9984D0
C------------------------------------------------------------------------
C solubility product for calcite
C
c Following Takahashi (1982) GEOSECS handbook
C NOT SURE THIS IS WORKING???
C Ingle et al. (1973)
c          Ksp_T_Calc = ( -34.452 - 39.866*(s**0.333333)
c    &                  + 110.21*log(s) - 7.5752d-6 * (tk**2.0)
c    &                  ) * 1.0d-7
c with pressure dependence Culberson and Pytkowicz (1968)
c          xvalue  =  (36-0.20*t)*(pressc-1.0)/(83.143*tk)
c          Ksp_TP_Calc = Ksp_T_Calc*exp(xvalue)
c
c
C Following Mucci (1983) - from Zeebe/Wolf-Gladrow equic.m
         tmpa1 = - 171.9065D0 - (0.077993D0*tk) + (2839.319D0/tk)
     &            + (71.595D0*log10(tk))
         tmpa2 = +(-0.77712D0 + (0.0028426D0*tk) + (178.34D0/tk) )*sqrts
         tmpa3 = -(0.07711D0*s) + (0.0041249D0*s15)
         logKspc = tmpa1 + tmpa2 + tmpa3
         Ksp_T_Calc = 10.0**logKspc
c        write(6,*)i,j,k,tmpa1,tmpa2,tmpa3,logkspc,Ksp_T_Calc
c with pressure dependence Culberson and Pytkowicz (1968)
c        xvalue  =  (36.0-0.20*t)*(pressc-1.0)/(83.143*tk)
c        Ksp_TP_Calc = Ksp_T_Calc*exp(xvalue)

c alternative pressure depdendence
c following Millero (1995) but using info from Appendix A11 of
c Zeebe and Wolf-Gladrow (2001) book
c          dv = -48.6 - 0.5304*t
c          dk = -11.76d-3 - 0.3692*t
c          xvalue = - (dv/(bigR*tk))*pressc
c    &               + (0.5*dk/(bigR*tk))*pressc*pressc
c          Ksp_TP_Calc = Ksp_T_Calc*exp(xvalue)

c alternative pressure dependence from Ingle (1975)

           zdum   = (pressc*10.D0 - 10.D0)/10.D0
           xvalue = ( (48.8D0 - 0.53D0*t)*zdum
     &                 + (-0.00588D0 + 0.0001845D0*t)*zdum*zdum)
     &            / (188.93D0*(t + 273.15D0))

           Ksp_TP_Calc = Ksp_T_Calc*10**(xvalue)

C------------------------------------------------------------------------
         else
C add Bennington
            fugf=0.D0
            ff=0.D0
            ak0= 0.D0
            ak1= 0.D0
            ak2= 0.D0
            akb= 0.D0
            ak1p = 0.D0
            ak2p = 0.D0
            ak3p = 0.D0
            aksi = 0.D0
            akw = 0.D0
            aks= 0.D0
            akf= 0.D0
            bt = 0.D0
            st = 0.D0
            ft = 0.D0
            Ksp_TP_Calc = 0.D0
         endif
        return
        end
