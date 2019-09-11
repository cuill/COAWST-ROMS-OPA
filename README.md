Oil-Particle-Aggregate module (OPAMOD)

OPAMOD is based on an existing, population dynamics-based flocculation model (FLOCMOD) in ROMS to account for the formation of OPAs. It was coupled with an oil plume model (https://github.com/DmitryDukh/COAWST-ROMS-OIL) to retrieve oil properties from Eulerian coordinates.

The main part of OPAMOD are two subroutines in ROMS/Nonlinear/Sediment:
-sed_opa.F
-sedopa_mod.h
