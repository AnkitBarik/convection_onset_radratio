# Data and scripts for article <br/> "Onset of convection in rotating spherical shells: variations with radius ratio"
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7359204.svg)](https://doi.org/10.5281/zenodo.7359204)
 [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/AnkitBarik/convection_onset_radratio/HEAD?labpath=figures_Bariketal.ipynb)

The data contains five folders for each Ekman number and within each there are 31 folders for each radius ratio
explored. Inside each of these subfolders there are in general four files. The first column of each file is radius.
The radius of the dissipation profile is scaled by the shell thickness L. For the rest of the files, the radius is
scaled by the outer boundary radius:

 - `Ekin_profile.dat` : Profile of kinetic energy with radius, integrated horizontally. Columns are poloidal, toroidal and total kinetic energy.
 - `KE_zavg_s.dat` : Profile of kinetic energy with cylindrical radius, interpolated and averaged in azimuth and along the rotation axis.
 - `dissipation_profile.dat` : Profiles of viscous dissipation with radius, integrated horizontally. Columns are poloidal, toroidal and total dissipation.
 - `temp_rms_profile.dat` : Profile of rms temperature perturbation with radius, integrated horizontally.
 - `modeShape.dat` : Datafile providing the spiral shape of mode in an equatorial slice. Columns are radius, longitude, x, y.

In addition we provide scripts to plot and analyze the data:

 - `dissThick_cmb.py` : Computes boundary layer thicknesses near the outer boundary from dissipation profiles and produces three figures : (1) a schematic of boundary layer thickness estimation, (2) a variation of boundary layer thickness with Ekman number and aspect ratio and (3) fit exponent and coefficient of thickness vs Ekman number at each radius ratio
 - `dissThick_icb.py` : Same as `dissThick_cmb.py` but for boundary layer near the inner boundary
 - `radial_extent.py` : Uses the files `KE_zavg_s.dat` to produce plots of normalized kinetic energy vs scaled cylindrical radius and shows a collapse with theoretical scaling.
 - `plot_diss_cmb.py` : Produces plots of normalized viscous dissipation vs scaled radius near the outer boundary and shows a collapse with the theoretical scaling.
 - `plot_diss_icb.py` : Same as above for the region near the inner boundary. Collapse is shown with respect to the scaling obtained from `dissThick_icb.py`.
 - `modeShape.py` : Plots the quantitative measure of spiralization of modes.

Lastly, in the directory `scripts`, we have provided the scripts `get_KE_profile.py`, `KE_zavg_s.py`, `get_dissipation_profile.py`, `get_temp_profile.py` and `spiral.py` used to obtain the abovementioned radial profiles in the respective order. These require the eigenmode solutions computed using `kore` (not provided with this data due to large size) and use of the kore visualization in the branch `merge_convection`. `KE_zavg_s.py` also uses interpolation code from the MagIC repository : `https://github.com/magic-sph/magic`.
