function [oxsat]=calc_oxygen_sat(t,s)
% COMPUTE_OXSAT:  compute oxygen saturation from temperature and salinity
%
% Oxygen saturation value is the volume of oxygen gas absorbed from humidity-saturated
% air at a total pressure of one atmosphere, per unit volume of the liquid at the temperature
% of measurement (ml/l)
%
% USAGE: osat = compute_oxsat ( temperature, salinity );
%
% PARAMETERS:
%    Input:
%        temperature:  in degrees celsius
%        salinity: in PSU
%    Output:
%        osat: in mol/m3
%



oA0=  2.00907;
oA1=  3.22014; 
oA2=  4.05010; 
oA3=  4.94457; 
oA4= -2.56847E-1; 
oA5=  3.88767; 
oB0= -6.24523E-3; 
oB1= -7.37614E-3; 
oB2= -1.03410E-2; 
oB3= -8.17083E-3; 
oC0= -4.88682E-7;

aTT = 298.15-t;
aTK = 273.15+t;
aTS = log(aTT./aTK);
aTS2= aTS.*aTS; 
aTS3= aTS2.*aTS;
aTS4= aTS3.*aTS; 
aTS5= aTS4.*aTS;

ocnew= exp(oA0 + oA1.*aTS + oA2.*aTS2 + oA3.*aTS3 + oA4.*aTS4 + oA5.*aTS5...
    + s.*(oB0 + oB1.*aTS + oB2.*aTS2 + oB3.*aTS3) + oC0.*(s.*s));

%Saturation concentration of dissolved O2"/units="mol/m3" 
oxsat = ocnew/22391.6*1000.0;
%Percent saturation of dissolved oxygen"/units="%" oxsat = (new_o2/o2sat)*100
%Apparent Oxygen Utilisation"/units="mol/m3" aou = o2sat-new_o2
return