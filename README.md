Oil-Particle-Aggregate module (OPAMOD)
Courtney K Harris ckharris.@vims.edu 
Linlin Cui lcui@vims.edu
Virginia Institute of Marine Sciences
09/12/2019

OPAMOD is based on an existing, population dynamics-based flocculation model (FLOCMOD) in ROMS to account for the formation of OPAs. It was coupled with an oil plume model (https://github.com/DmitryDukh/COAWST-ROMS-OIL) to retrieve oil properties from Eulerian coordinates.

The main part of OPAMOD are two subroutines in ROMS/Nonlinear/Sediment:
-sed_opa.F
-sedopa_mod.h
