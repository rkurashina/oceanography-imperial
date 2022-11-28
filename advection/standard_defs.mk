
CFLAGS = -O2 -g -heap-arrays -traceback  -fp-stack-check -mcmodel=large
NETCDF_VERSION=4.3.3
CC = ifort

LIBS = -L${MKL_HOME}/interfaces/fftw3xf -L/apps/netcdf/$(NETCDF_VERSION)/lib -L${MKL_HOME}/interfaces -L${MKL_HOME}/lib/intel64
LINKS = -lfftw3 -lnetcdff -lnetcdf -lmkl_intel_lp64 -lmkl_sequential -lmkl_core

HOME_DIR = /rds/general/user/rk2014/home/WORK/advection
NETCDF_DIR = $(HOME_DIR)/NETCDF
LAGR_DIR = $(HOME_DIR)/LAGR
VAR_DIR = $(HOME_DIR)/VAR
SRC_DIR = $(HOME_DIR)/TRANSPORT
SOLVER_DIR = $(HOME_DIR)/SOLVER
INPUT_DIR = $(HOME_DIR)/INPUT
BUILD_DIR = $(HOME_DIR)/BUILD
MISC_DIR = $(HOME_DIR)/MISC


INCLUDES = -I${MKL_HOME}/include/fftw -I${MKL_HOME}/include/intel64/lp64 -I/apps/netcdf/$(NETCDF_VERSION)/include -I$(BUILD_DIR) -I$(SRC_DIR)
