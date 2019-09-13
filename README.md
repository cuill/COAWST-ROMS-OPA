Oil-Particle-Aggregate module (OPAMOD)

Courtney K Harris ckharris.@vims.edu 
Linlin Cui lcui@vims.edu
Virginia Institute of Marine Science
09/12/2019

OPAMOD is based on an existing, population dynamics-based flocculation model (FLOCMOD) in ROMS to account for the formation of OPAs. It was coupled with an oil plume model (https://github.com/DmitryDukh/COAWST-ROMS-OIL, July 31, 2019) to retrieve oil properties from Eulerian coordinates.

The main part of OPAMOD are two subroutines in ROMS/Nonlinear/Sediment:
-sed_opa.F
-sedopa_mod.h

At this time, trunk is not updated to the latest version of Dmitry's oil model (as Septempber 4, 2019).

