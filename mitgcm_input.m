% This is a script that executes a series of commands for
% creating climate data input for the MITgcm code.
%
% Created  08/15/99 by adcroft@mit.edu
% Modified 11/11/99 by adcroft@mit.edu
% Maintained by adcroft@mit.edu, abiastoch@ucsd.edu

% Specify binary file formats
setfmt

% Generate grid
setgrid

% Generate/edit topography
%topo

% Generate PMASK file
gen_pmask

% Extract/interpolate/fill Levitus data
%levit

% Extract/interpolate/fill daSilva surface data
dasilva

