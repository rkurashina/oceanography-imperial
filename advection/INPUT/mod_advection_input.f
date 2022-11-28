C MODULE THAT CONTAINS INPUT FOR TRANSPORT MODEL

      module mod_advection_input

      implicit none

      integer npoints_sqrt,npoints,isolve,i_pseudo,i_full,i_eddy,regime
     & ,sim,ii,jj,nbins,i_pvbin,jj1
      real*8 H1,H2,Rd,basinscale,visc,U_0,beta
      real*8 release_interval,release_length
      integer release_no
      integer k_save
      real*8 dt

c DYNAMICAL MODEL PARAMETERS

      parameter(ii=513,jj=513,jj1 = 3*jj)
      parameter(basinscale=3840.D5
     & ,H1=1.D5,H2=3.D5
     & ,beta = 2.D-13
     & ,U_0=0.D0
     & ,Rd= 25.D5)            !don't think H1, H2, U_0, Rd matter for uniform binning

C TRANSPORT MODEL PARAMETERS

      parameter(npoints=50) ! NUMBER OF PARTICLES
      parameter(isolve=1 ! if isolve = 0 :bicubic, if isolve = 1 :2Dcubic, 2D CUBIC IS BEST!
     & ,i_eddy=0 ! CALCULATE EDDY ONLY TRAJECTORY
     & ,i_pseudo=0 ! CALCULATE FFE TRAJECTORY
     & ,i_full=1 ! CALCULATE FULL TRAJECTORY
     & ,release_interval = 365 ! RELEASE PARTICLES EVERY release_interval DAYS
     & ,release_length = 365. ! ADVECT PARTICLE FOR release_length DAYS
     & ,release_no = 2! NUMBER OF TIME RELEASES (DO 2 X 2 YEAR RELASES)
     & ,k_save = 1 ! save every k_save days
     & ,dt = 1080 ! TIME STEP
     & ,nbins = 1 ! NUMBER OF BINS
     & ,i_pvbin = 0) ! set = 1 if you wish to bin domain according to PV, otherwise the domain is binned uniformly

      character*(*),parameter :: home_dir =
     & '/rds/general/user/rk2014/home/WORK/advection
     &  /Advection_3_layers/'

       character*(*), parameter :: file_name = 'QG_new.nc' ! NAME OF FILE CONTAINING STREAM FUNCTION DATA
       character*(*), parameter :: ave_name = 'QG_ave.nc' ! NAME OF FILE CONTAINING TIME AVERAGED STREAM FUNCTION DATA

      end module mod_advection_input
