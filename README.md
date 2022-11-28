# Q-GCM, CABARET and advection

Repository containing files for Q-GCM model, original repository is [here](https://github.com/GFDANU/q-gcm). Documentation for the model may be found [here](http://www.q-gcm.org/downloads.html). 

Code is essentially identical to original except very small modificaitons have been made for use on Imperial College HPC server. 

To compile:
* module load intel-suite 
* module load netcdf/3.6.3
* make 

To run:
* qsub submit.pbs
