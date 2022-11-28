# Ocean models
Repository containing model codes for Q-GCM, CABARET and Lagrangian particle model (for CABARET).

# Q-GCM
Original repository is [here](https://github.com/GFDANU/q-gcm). Documentation for the model may be found [here](http://www.q-gcm.org/downloads.html). 

Code is essentially identical to original except very small modificaitons have been made for use on Imperial College HPC server. 

**To compile:**
* module load intel-suite 
* module load netcdf/3.6.3
* make 

**To run:**
* qsub submit.pbs

# CABARET
CABARET 3L QG double-gyre model. See [here]{https://doi.org/10.1016/j.ocemod.2009.06.009} for description of particular adection scheme used. It is especially good at conserving PV so allows for model to be run at a coarser resolutions.

**To compile:**
* module load intel-suite 
* module load netcdf
* make -f make_cabaret3

Note for CABARET and Lagrangian codes, the more recent versions of netCDF (4.3.3 or later) are used.

**To run:**
* qsub cabaret.pbs

# Lagrangian particles
Advects Lagrangian particles around flow outputted by CABARET. See [here]{https://doi.org/10.1017/jfm.2021.384} for example of use.

**To compile:**
* module load intel-suite 
* module load netcdf
* make -f make_advection

**To run**
* qsub advection.pbs




