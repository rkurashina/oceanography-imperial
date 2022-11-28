c MODULE THAT CONTAINS INPUT FOR DYNAMICAL MODEL

      module mod_qg_input
      
      implicit none
      
      integer ii,jj,ips,max_spin,max_time,istart
     & ,regime,istart_ave,ii1,jj1,ii2,jj2
      real*8 basinscale,H1,H2,h_max,Rd,visc,visc_bot,U_0
     & ,read_tim,time_out,time_save,cfl 
      
      
      parameter(ips=1,max_spin=0,max_time=20000+max_spin  ! max_time is how long in days you wish to run the code for.
     & ,ii=512,jj=512)                                                   !! Also change ii, jj in solv_ell_mike

      parameter(basinscale=520.D5
     & ,H1=1.D5,H2=3.D5,h_max=1.D2,Rd=25.D5 ! layer depths and rossby deformation radius
     & ,visc=1.D4
     & ,U_0=6.D0 ! velocity in upper layer
     & ,istart=1 ! set to 0 if starting from scratch, otherwise set to 1 if starting from a restart file
     & ,istart_ave = 1 ! set to 0 if you wish to write a new time-averaged file, otherwise set to 1 if starting from a restart file
     & ,regime=2 ! alters the bottom friction to producee different jets - see qg2_dp.f for the bottom friction values
     & ,read_tim=0
     & ,TIME_OUT=10.D0     ! accumulate data every TIME_OUT days
     & ,TIME_SAVE=1.D0	   ! save every time_save days
     & ,CFL=0.4 ) ! courant number, set to be small if spinning up.

      parameter(ii1=ii-1,jj1=jj-1,ii2=ii-2,jj2=jj-2)
     
      character*(*),parameter :: home_dir = 
     & '/work/jp1115/saves/'
     
      character*(*),parameter :: file_name = 'QG_new.nc' ! saves time record of stream function
      character*(*),parameter :: ave_file = 'QG_ave_new.nc' ! saves time averages stream function


      end module mod_qg_input
