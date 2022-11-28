# Ocean models
Repository containing model codes for Q-GCM, CABARET and Lagrangian particle model (for CABARET).

## Q-GCM
Original repository is [here](https://github.com/GFDANU/q-gcm). Documentation for the model may be found [here](http://www.q-gcm.org/downloads.html). 

Code is essentially identical to original except very small modificaitons have been made for use on Imperial College HPC server. 

**To compile:**
* module load intel-suite 
* module load netcdf/3.6.3
* make 

**To run:**
* qsub submit.pbs

## CABARET
CABARET 3L QG double-gyre model. See [here](https://www.sciencedirect.com/science/article/abs/pii/S1463500309001267) for description of particular adection scheme used. It is especially good at conserving PV so allows for model to be run at a coarser resolutions.

**To compile:**
* module load intel-suite 
* module load netcdf
* make -f make_cabaret3

Note for CABARET and Lagrangian codes, the more recent versions of netCDF (4.3.3 or later) are used.

**To run:**
* qsub cabaret.pbs

## Lagrangian particles
Advects Lagrangian particles around flow outputted by CABARET. See [here](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/western-boundary-layer-nonlinear-control-of-the-oceanic-gyres/F204931DF6499D3CE94ECF6487FD820F) for example of use.

**To compile:**
* module load intel-suite 
* module load netcdf
* make -f make_advection

**To run**
* qsub advection.pbs

# Post-processing
Some simple example python scripts are included to post-process data. This includes files for:
* EOF decompositon
* (Lagged) SVDs
* Low-pass filtering 
These are relatively simple scripts and may adapted for your requirements. 




