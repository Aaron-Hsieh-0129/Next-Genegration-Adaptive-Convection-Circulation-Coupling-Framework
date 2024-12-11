DSET ^0.nc
DTYPE netcdf
OPTIONS template
TITLE NetCDF Data for GrADS
UNDEF -9999.0
XDEF 502 LINEAR 0 1   # Update if x-spacing isn’t uniform
YDEF 1 LINEAR 1 1     # Single y-point for 2D model
ZDEF 102 LINEAR 0 200   # z-spacing or levels
TDEF 600 LINEAR 00Z31dec1999 2mn # date

VARS 9
th=>th 102 t,x,z theta 
u=>u 102 t,x,z u 
w=>w 102 t,x,z w 
zeta=>zeta 102 t,x,z zeta 
qv=>qv 102 t,x,z qv 
qc=>qc 102 t,x,z qc
qr=>qr 102 t,x,z qr 
precip=>precip 1 t,x precip
ubarTop=>ubarTop 1 t ubarTop
ENDVARS

