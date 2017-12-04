# mitgcm_toolbox
# 
These are all the MATLAB routines that I use to analyse MITgcm output.
Some come packaged with MITgcm in the "utils/matlab" or "verification"
folders and some I have (re)written myself. Probably there are many hard
coded lines that need changing.

To use these routines, you should check for additional dependencies,
such as the Seawater Toolbox (I've not updated to GSW here yet) and the
mexcdf/snctools package for interacting with netcdf output.

A good place to start is "mit_loadglobal" which runs through my output
and produces a bunch of figures, including global residual overturning
circulation (mit_calc_psi, mit_plotstreamfunctions), biogeochemical
fields (mit_plotbiogeochem) and global averages (mit_plotmeandrift).

To make use of the carbon component routines (calc_sat_dic, calc_pco2)
you must use "mex" and a fortran complier to compile the mit_calc_csat,
mit_calc_pco2 and carbon_chem code (ask me how if you need a hand here).

If you find these interesting, please cite:

Lauderdale JM, Naveira Garabato AC, Oliver KIC, Follows MJ & Williams RG
(2013) Wind-driven changes in Southern Ocean residual circulation, ocean
carbon reservoirs and atmospheric CO2. Clim Dyn, 41(7), 2145–2164,
doi:10.1007/s00382-012-1650-3

and

Lauderdale JM, Williams RG, Munday DR & Marshall DP (2016) The impact of
Southern Ocean residual upwelling on atmospheric CO2 on centennial and
millennial timescales. Clim Dyn, 48, 1611-1632,
doi:10.1007/s00382-016-3163-y

And contact me if you find any bugs or need help with this code.
