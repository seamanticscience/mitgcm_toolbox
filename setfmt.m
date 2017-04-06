% setfmt.m specifies the binary precision and IEEE format for
% storing/reading data to files.
%  fmt = 'real*4' or 'real*8'
%  Ieee = 'b', 'l' etc...
%
% Created  11/11/99 by adcroft@mit.edu
% Modified          by
% Maintained by adcroft@mit.edu, abiastoch@ucsd.edu

fmt='real*4';
Ieee='l';

save FMT.mat fmt Ieee
