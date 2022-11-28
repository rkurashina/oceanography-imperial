c PROGRAM THAT PERFORMS OFFLINE ADVECTION
c RELEASES PARTICLES UNIFORMLY AT DIFFERENT TIMES


c gfortran -o advect.exe mod_1dinterp.f mod_2dcubic.f mod_advection_input.f mod_bicubic.f mod_laplace.f mod_PVbin_advection_constants.f mod_qg2_netcdf.f mod_random.f mod_rk4.f mod_time_interp.f mod_traj_netcdf.f mod_variables.f offline_transport.f offline_transport.f


      program offline_transport

      use mod_advection_input !! IMPORTANT - THIS DEFINES THE INPUT VARIABLES - IN THE INPUT FOLDER
      use mod_PVbin_advection_constants
c      use mod_advection_constants
      use mod_variables
      use mod_qg3_netcdf
      use mod_traj3_netcdf
      use mod_2dcubic
      use mod_bicubic
      use mod_random
      use mod_time_interp3
      use mod_1dinterp
      use mod_laplace
      use mod_rk4

      implicit none

c ---- DEFINE DIMENSIONAL VARIABLES

      uscale=1.
      scale=basinscale/dfloat(ii-1)
      tscale=scale/uscale
      SS=(scale/Rd)**2

      S1=SS/(1.+H1/H2)
      S2=SS-S1
      beta_nondim=beta*scale*scale/uscale

      BETA_NONDIM_U1=BETA_NONDIM+U_0*S1
      BETA_NONDIM_U2=BETA_NONDIM-U_0*S2

        do j = 1,jj1
            y_c(j) = dfloat(j-1)
        enddo

c --- READ TIME AND TIME AVERAGED STREAM FUNCTION ----

      call read_time(file_name,time,t_len)

c ----- CALCULATE THE TRUE TIME ARRAY (WHICH TAKES INTO ACCOUNT THE CHOSEN OFFLINE TIME STEP)

        release_length_secs = release_length*86400 ! time recording in seconds for a release
        max_time_lag = time(t_len)*86400

        t_tot = int(release_length_secs/dt) ! total number of time steps per release
        l_tot_day = int(86400/dt) ! total number of time steps per day

        dt_nondim = dt/tscale

        allocate(time_o(t_tot))
        allocate(time_dim(t_tot))

c ------ CREATE TRAJECTORY DATA FILE ----------------

      if(i_pvbin.eq.1) then

      if (i_eddy .eq. 1) then
        eddy_name = 'eddy_PV_bins_trajectories.nc'
        call create_binned_file(eddy_name,nbins,npoints,release_no)
      endif
      if(i_pseudo .eq. 1) then
         pseudo_name = 'pseudo_PV_bins_trajectories.nc'
         call create_binned_file(pseudo_name,nbins,npoints,release_no)
      endif
      if (i_full .eq. 1) then
          full_name = 'full_PV_bins_trajectories.nc'
          call create_binned_file(full_name,nbins,npoints,release_no)
      endif

      else

      if (i_eddy .eq. 1) then
        eddy_name = 'eddy_uniform_bins_trajectories.nc'
        call create_binned_file(eddy_name,nbins,npoints,release_no)
      endif
      if(i_pseudo .eq. 1) then
         pseudo_name = 'pseudo_uniform_bins_trajectories.nc'
         call create_binned_file(pseudo_name,nbins,npoints,release_no)
      endif
      if (i_full .eq. 1) then
          full_name = 'full_uniform_bins_trajectories.nc'
          call create_binned_file(full_name,nbins,npoints,release_no)            !this is the one we want
      endif

      endif


c ----- Initialise trajectory arrays

        if((i_full.eq.1).or.(i_pseudo.eq.1)) then
        allocate(x1r(nbins,npoints),
     &   x2r(nbins,npoints), x3r(nbins,npoints)
     & ,y1r(nbins,npoints),y2r(nbins,npoints), y3r(nbins,npoints)
     & ,x1r_coord(nbins,npoints)
     & ,x2r_coord(nbins,npoints)
     & ,x3r_coord(nbins,npoints)
     & ,y1r_coord(nbins,npoints)
     & ,y2r_coord(nbins,npoints)
     & ,y3r_coord(nbins,npoints))
        endif

        if (i_eddy.eq.1) then
        allocate(x1r_eddy(nbins,npoints)
     & ,x2r_eddy(nbins,npoints)
     & ,y1r_eddy(nbins,npoints)
     & ,y2r_eddy(nbins,npoints)
     & ,x1r_eddy_coord(nbins,npoints)
     & ,x2r_eddy_coord(nbins,npoints)
     & ,y1r_eddy_coord(nbins,npoints)
     & ,y2r_eddy_coord(nbins,npoints))
        endif

        if (i_pseudo.eq.1) then
        allocate(x1r_pseudo(nbins,npoints)
     & ,x2r_pseudo(nbins,npoints)
     & ,y1r_pseudo(nbins,npoints)
     & ,y2r_pseudo(nbins,npoints)
     & ,x1r_pseudo_coord(nbins,npoints)
     & ,x2r_pseudo_coord(nbins,npoints)
     & ,y1r_pseudo_coord(nbins,npoints)
     & ,y2r_pseudo_coord(nbins,npoints)
     & ,x1r_mean(nbins,npoints)
     & ,y1r_mean(nbins,npoints)
     & ,x2r_mean(nbins,npoints)
     & ,y2r_mean(nbins,npoints))
        endif

        allocate(release_time(release_no) ,nrel(release_no))
        allocate(k_s(release_no))


        do k = 1,release_no
            release_time(k) = (k-1)*release_interval
            nrel(k) = 1 ! counters for writing to file
            k_s(k) = 0
        enddo


c ---- Calculate coefficients for time-averaged stream function
c ---- Needed for spatial interpolation
c ------- ONLY FOR PSEUDO TRAJECTORIES -----------


      if(i_pseudo.eq.1) then


        if (isolve.eq.0) then
            call A_matrix(ii,jj,psi1_av,M1_av)
            call A_matrix(ii,jj,psi2_av,M2_av)


        else if (isolve.eq.1) then

            call cubic_coeff_x(ii,jj,psi1_av
     & ,a1_av,b1_av,c1_av,d1_av)
            call cubic_coeff_x(ii,jj,psi2_av
     & ,a2_av,b2_av,c2_av,d2_av)

        endif

      endif

        iseed = 102

        nrec = 0 ! intiliase counter

        do k = 1,release_no ! loop through different releases

        print*, 'Starting particle advection for release number',k

        do t = 1,t_tot ! define the offline and dimensional time arrays
            time_o(t) = (t-1)*dt_nondim + release_time(k)*86400/tscale
            time_dim(t) = time_o(t)*tscale/86400
        enddo

c ------------ MAIN CYCLE -----------------------------
        do t = 1,t_tot

        if (t.eq.1) then

            if (i_pvbin.eq.1) then

c ------- CONSTRUCT BINS UNIFORM IN PV --------------------
        print*,"this code doesn't do that"
        else

c ---------------------------------------------------------------------
c ---------------------- UNIFORMLY BIN THE DOMAIN AND SEED PARTICLES -----------
c -----------------------------------------------------------------

c bin domain

      write(*,*) 'Binning domain uniformally'

      d_bin = dfloat(jj)/nbins

      do p = 1,nbins
        y_bin(p) = p*d_bin
      enddo

      do p = 1,nbins
      do n = 1,npoints

        x0 = ran1(iseed)*dfloat(ii-1)
        if (p.eq.1) then
        y0 = ran1(iseed)*d_bin
        else
        y0 = (ran1(iseed)*d_bin)+y_bin(p-1)
        endif

        if (i_full.eq.1. .or. i_pseudo.eq.1) then
            x1r(p,n) = x0
            y1r(p,n) = y0
            x2r(p,n) = x0
            y2r(p,n) = y0
            x3r(p,n) = x0
            y3r(p,n) = y0

            x1r_coord(p,n) = 0
            y1r_coord(p,n) = 0
            x2r_coord(p,n) = 0
            y2r_coord(p,n) = 0
            x3r_coord(p,n) = 0
            y3r_coord(p,n) = 0
        endif
        if (i_eddy.eq.1) then
            x1r_eddy(p,n) = x0
            y1r_eddy(p,n) = y0
            x2r_eddy(p,n) = x0
            y2r_eddy(p,n) = y0

            x1r_eddy_coord(p,n) = 0
            y1r_eddy_coord(p,n) = 0
            x2r_eddy_coord(p,n) = 0
            y2r_eddy_coord(p,n) = 0
        endif
        if (i_pseudo .eq. 1) then
            x1r_pseudo(p,n) = x0
            y1r_pseudo(p,n) = y0
            x2r_pseudo(p,n) = x0
            y2r_pseudo(p,n) = y0

            x1r_pseudo_coord(p,n) = 0
            y1r_pseudo_coord(p,n) = 0
            x2r_pseudo_coord(p,n) = 0
            y2r_pseudo_coord(p,n) = 0
        endif

      enddo
      enddo


      do p=1,nbins

        if (i_full.eq.1) then


            x1 = x1r(p,:)
            x2 = x2r(p,:)
            x3 = x3r(p,:)
            y1 = y1r(p,:)
            y2 = y2r(p,:)
            y3 = y3r(p,:)


            x1_coord = x1r_coord(p,:)
            x2_coord = x2r_coord(p,:)
            x3_coord = x3r_coord(p,:)
            y1_coord = y1r_coord(p,:)
            y2_coord = y2r_coord(p,:)
            y3_coord = y3r_coord(p,:)           

            call write_binned_file(full_name,npoints,p,k
     &       ,x1,y1
     &       ,x2,y2
     &       ,x3,y3
     &       ,x1_coord,y1_coord,x2_coord
     &       ,y2_coord,x3_coord,y3_coord
     &       ,release_time(k),nrel(k))
            
            endif

        enddo
        nrel(k) = nrel(k) + 1


        endif

        else ! do if (t .ne. 1), i.e we now advect particles


         time_day = time_o(t)*tscale/86400.
         k_s(k) = k_s(k) + 1

c ----- find entry m in time such that time(m) < time_o(k) <= time(m+1)
        do m = 1,t_len
            if (time(m)-time(1) >= time_dim(t)) then
            k_new = m
                if(k_new .eq. 2)then
                k_new = 3
                else if (k_new.eq.t_len)then
                k_new = k_new-1
                endif
            exit
            endif
        enddo
c ----- do same for k-1
        do m = 1,t_len
            if (time(m)-time(1) >= time_dim(t-1)) then
            k_old = m
                if(k_old .eq. 2)then
                k_old = 3
                elseif(k_old.eq.1)then
                k_old = 3
                else if (k_old.eq.t_len)then
                k_old= k_old-1
                endif
            exit
            endif
        enddo



        time_half = (time_dim(t) + time_dim(t-1))/2

        do m = 1,t_len
            if (time(m)-time(1) >= time_half) then
            k_half = m
                if(k_half .eq. 2)then
                k_half = 3
                else if (k_half.eq.t_len)then
                k_half= k_half-1
                endif
            exit
            endif
        enddo


c ----- cubic interpolate psi in time to find psi
c ----- at time_o(k), time_o(k-1)

c ---- time_o(k)
        call read_psi(file_name,ii,jj,k_new-2,psi1,psi2,psi3)

        time_cubic = time(k_new-2:k_new+1)

        call interp_time(ii,jj,time_cubic,time_dim(t),psi1,psi1_new)
        call interp_time(ii,jj,time_cubic,time_dim(t),psi2,psi2_new)
        call interp_time(ii,jj,time_cubic,time_dim(t),psi3,psi3_new)
        
c ----- time_o(k-1)

        call read_psi(file_name,ii,jj,k_old-2,psi1,psi2,psi3)


        time_cubic = time(k_old-2:k_old+1)

        call interp_time(ii,jj,time_cubic,time_dim(t-1),psi1,psi1_old)
        call interp_time(ii,jj,time_cubic,time_dim(t-1),psi2,psi2_old)
        call interp_time(ii,jj,time_cubic,time_dim(t-1),psi3,psi3_old)

c ----- time_half

        call read_psi(file_name,ii,jj,k_half-2,psi1,psi2,psi3)

        time_cubic = time(k_half-2:k_half+1)

        call interp_time(ii,jj,time_cubic,time_half,psi1,psi1_half)
        call interp_time(ii,jj,time_cubic,time_half,psi2,psi2_half)
        call interp_time(ii,jj,time_cubic,time_half,psi3,psi3_half)


            if (i_eddy.eq.1) then

            ! determine eddy streamfunction

            do i = 1,ii
            do j = 1,jj
                psi1_eddy_old(i,j) = psi1_old(i,j) - psi1_av(i,j)
                psi2_eddy_old(i,j) = psi2_old(i,j) - psi2_av(i,j)


                psi1_eddy_half(i,j) = psi1_half(i,j) - psi1_av(i,j)
                psi2_eddy_half(i,j) = psi2_half(i,j) - psi2_av(i,j)


                psi1_eddy_new(i,j) = psi1_new(i,j) - psi1_av(i,j)
                psi2_eddy_new(i,j) = psi2_new(i,j) - psi2_av(i,j)

            enddo
            enddo

            endif

            !print*,'eddy created'

c --------------- CALCULATE COEFFICIENTS FOR SPATIAL INTERPOLATION -----------

c --------  BICUBIC INTERPOLATION -----------
      if (isolve.eq.0) then

      if(t.eq.2) then

      if ((i_pseudo .eq. 1) .or. (i_full.eq.1)) then

      call A_matrix(ii,jj,psi1_old,M1_old)
      call A_matrix(ii,jj,psi2_old,M2_old)

      endif

      if (i_eddy.eq.1) then

      call A_matrix(ii,jj,psi1_eddy_old,M1_eddy_old)
      call A_matrix(ii,jj,psi2_eddy_old,M2_eddy_old)

      endif

      else

      M1_old = M1_new
      M2_old = M2_new

      if (i_eddy.eq.1) then

      M1_eddy_old = M1_eddy_new
      M2_eddy_old = M2_eddy_new

      endif

      endif

      if ((i_pseudo .eq. 1) .or. (i_full.eq.1)) then


      call A_matrix(ii,jj,psi1_half,M1_half)
      call A_matrix(ii,jj,psi2_half,M2_half)
      call A_matrix(ii,jj,psi1_new,M1_new)
      call A_matrix(ii,jj,psi2_new,M2_new)

      endif


      if (i_eddy.eq.1) then
      call A_matrix(ii,jj,psi1_eddy_half,M1_eddy_half)
      call A_matrix(ii,jj,psi2_eddy_half,M2_eddy_half)
      call A_matrix(ii,jj,psi1_eddy_new,M1_eddy_new)
      call A_matrix(ii,jj,psi2_eddy_new,M2_eddy_new)
      endif

c --------- 2D CUBIC INTERPOLATION


      else if (isolve.eq.1) then

      if (t.eq.2) then

      if ((i_pseudo .eq. 1) .or. (i_full.eq.1)) then

      call cubic_coeff_x(ii,jj,psi1_old
     & ,a1_old,b1_old,c1_old,d1_old)
      call cubic_coeff_x(ii,jj,psi2_old
     & ,a2_old,b2_old,c2_old,d2_old)
      call cubic_coeff_x(ii,jj,psi3_old
     & ,a3_old,b3_old,c3_old,d3_old)

      endif

      if (i_eddy.eq.1) then
      call cubic_coeff_x(ii,jj,psi1_eddy_old
     & ,a1_eddy_old,b1_eddy_old,c1_eddy_old,d1_eddy_old)
      call cubic_coeff_x(ii,jj,psi2_eddy_old
     & ,a2_eddy_old,b2_eddy_old,c2_eddy_old,d2_eddy_old)
      endif

      else

      if ((i_pseudo .eq. 1) .or. (i_full.eq.1)) then

      a1_old = a1_new
      a2_old = a2_new
      a3_old = a3_new
      b1_old = b1_new
      b2_old = b2_new
      b3_old = b3_new
      c1_old = c1_new
      c2_old = c2_new
      c3_old = c3_new
      d1_old = d1_new
      d2_old = d2_new
      d3_old = d3_new

      endif

      if(i_eddy.eq.1) then

      a1_eddy_old = a1_eddy_new
      a2_eddy_old = a2_eddy_new
      b1_eddy_old = b1_eddy_new
      b2_eddy_old = b2_eddy_new
      c1_eddy_old = c1_eddy_new
      c2_eddy_old = c2_eddy_new
      d1_eddy_old = d1_eddy_new
      d2_eddy_old = d2_eddy_new
      endif

      end if

      if ((i_pseudo .eq. 1) .or. (i_full.eq.1)) then


      call cubic_coeff_x(ii,jj,psi1_half
     &,a1_half,b1_half,c1_half,d1_half)
      call cubic_coeff_x(ii,jj,psi2_half
     &,a2_half,b2_half,c2_half,d2_half)
      call cubic_coeff_x(ii,jj,psi3_half
     &,a3_half,b3_half,c3_half,d3_half)
      call cubic_coeff_x(ii,jj,psi1_new
     &,a1_new,b1_new,c1_new,d1_new)
      call cubic_coeff_x(ii,jj,psi2_new
     &,a2_new,b2_new,c2_new,d2_new)
      call cubic_coeff_x(ii,jj,psi3_new
     &,a3_new,b3_new,c3_new,d3_new)

      endif

      if(i_eddy.eq.1) then
      call cubic_coeff_x(ii,jj,psi1_eddy_half
     &,a1_eddy_half,b1_eddy_half,c1_eddy_half,d1_eddy_half)
      call cubic_coeff_x(ii,jj,psi2_eddy_half
     &,a2_eddy_half,b2_eddy_half,c2_eddy_half,d2_eddy_half)
      call cubic_coeff_x(ii,jj,psi1_eddy_new
     &,a1_eddy_new,b1_eddy_new,c1_eddy_new,d1_eddy_new)
      call cubic_coeff_x(ii,jj,psi2_eddy_new
     &,a2_eddy_new,b2_eddy_new,c2_eddy_new,d2_eddy_new)
      endif

      endif



c ------------------ RK4 - FULL ADVECTION ---------------------------


c FULL ADVECTION

            if ((i_full.eq.1).or.(i_pseudo.eq.1)) then
            do p = 1,nbins
            do n = 1,npoints

        if (isolve.eq.0) then

        call rk4_bicubic(ii,jj,x1r(p,n),y1r(p,n)
     &   ,dt_nondim,U_0,M1_old,M1_half,M1_new
     & ,x_diff1,y_diff1)
        call rk4_bicubic(ii,jj,x2r(p,n),y2r(p,n)
     &   ,dt_nondim,dfloat(0),M2_old,M2_half,M2_new
     & ,x_diff2,y_diff2)

        elseif (isolve.eq.1) then


        call rk4_2dcubic(ii,jj,x1r(p,n),y1r(p,n)
     &   ,dt_nondim,U_0,a1_old,b1_old,c1_old
     & ,d1_old,a1_half,b1_half,c1_half,d1_half
     & ,a1_new,b1_new,c1_new,d1_new
     & ,x_diff1,y_diff1)
        call rk4_2dcubic(ii,jj,x2r(p,n),y2r(p,n)
     &   ,dt_nondim,dfloat(0),a2_old,b2_old,c2_old
     & ,d2_old,a2_half,b2_half,c2_half,d2_half
     & ,a2_new,b2_new,c2_new,d2_new
     & ,x_diff2,y_diff2)
        call rk4_2dcubic(ii,jj,x3r(p,n),y3r(p,n)
     &   ,dt_nondim,dfloat(0),a3_old,b3_old,c3_old
     & ,d3_old,a3_half,b3_half,c3_half,d3_half
     & ,a3_new,b3_new,c3_new,d3_new
     & ,x_diff3,y_diff3)

        endif

      if (y_diff1 .gt. 20) then

        print*,'y_diff1 = ',y_diff1
        print*,'psi_old = ',psi1_old
       stop
      endif

c MEAN CONTRIBUTION FOR THE PSEUDO TRAJECTORIES

        if(i_pseudo.eq.1) then

                 if (isolve.eq.0) then

        call rk4_bicubic(ii,jj,x1r(p,n),y1r(p,n),dt_nondim,U_0
     & ,M1_av,M1_av,M1_av
     & ,x_av_diff1,y_av_diff1)
        call rk4_bicubic(ii,jj,x2r(p,n),y2r(p,n),dt_nondim,dfloat(0)
     & ,M2_av,M2_av,M2_av
     & ,x_av_diff2,y_av_diff2)

        elseif (isolve.eq.1) then


        call rk4_2dcubic(ii,jj,x1r(p,n),y1r(p,n),dt_nondim,U_0
     & ,a1_av,b1_av,c1_av
     & ,d1_av,a1_av,b1_av,c1_av,d1_av
     & ,a1_av,b1_av,c1_av,d1_av
     & ,x_av_diff1,y_av_diff1)
        call rk4_2dcubic(ii,jj,x2r(p,n),y2r(p,n),dt_nondim,dfloat(0)
     & ,a2_av,b2_av,c2_av
     & ,d2_av,a2_av,b2_av,c2_av,d2_av
     & ,a2_av,b2_av,c2_av,d2_av
     & ,x_av_diff2,y_av_diff2)

        endif

        ENDIF

        !print*, 'full location'
        !print*, x1r(n,t), y1r(n,t)


        x1r(p,n) = x1r(p,n) + x_diff1
        y1r(p,n) = y1r(p,n) + y_diff1
        x2r(p,n) = x2r(p,n) + x_diff2
        y2r(p,n) = y2r(p,n) + y_diff2
        x3r(p,n) = x3r(p,n) + x_diff3
        y3r(p,n) = y3r(p,n) + y_diff3


        if(x1r(p,n).lt.0)then
            x1r(p,n) = -x1r(p,n)
        end if
        if(x1r(p,n).gt.dfloat(ii-1))then
            x1r(p,n) = 2.D0*dfloat(ii-1) - x1r(p,n)
        endif
        if(y1r(p,n).lt.0.)then
          y1r(p,n)=-y1r(p,n)
        endif
        if(y1r(p,n).gt.dfloat(jj-1))then
          y1r(p,n)=2.D0*dfloat(jj-1) - y1r(p,n)
        endif

        if(x2r(p,n).lt.0)then
            x2r(p,n) = -x2r(p,n)
        end if
        if(x2r(p,n).gt.dfloat(ii-1))then
            x2r(p,n) = 2.D0*dfloat(ii-1) - x2r(p,n)
        endif
        if(y2r(p,n).lt.0.)then
          y2r(p,n)=-y2r(p,n)
        endif
        if(y2r(p,n).gt.dfloat(jj-1))then
          y2r(p,n)=2.D0*dfloat(jj-1) - y2r(p,n)
        endif

        if(x3r(p,n).lt.0)then
            x3r(p,n) = -x3r(p,n)
        end if
        if(x3r(p,n).gt.dfloat(ii-1))then
            x3r(p,n) = 2.D0*dfloat(ii-1) - x3r(p,n)
        endif
        if(y3r(p,n).lt.0.)then
          y3r(p,n)=-y3r(p,n)
        endif
        if(y3r(p,n).gt.dfloat(jj-1))then
          y3r(p,n)=2.D0*dfloat(jj-1) - y3r(p,n)
        endif


            if (i_pseudo.eq.1) then

        x1r_pseudo(p,n) = x1r_pseudo(p,n) + x_diff1 - x_av_diff1
        x2r_pseudo(p,n) = x2r_pseudo(p,n) + x_diff2 - x_av_diff2
        y1r_pseudo(p,n) = y1r_pseudo(p,n) + y_diff1 - y_av_diff1
        y2r_pseudo(p,n) = y2r_pseudo(p,n) + y_diff2 - y_av_diff2

                    if(x1r_pseudo(p,n).lt.0)then
                x1r_pseudo(p,n) = dfloat(ii)+x1r_pseudo(p,n)
                x1r_pseudo_coord(p,n)=x1r_pseudo_coord(p,n)-1
            end if
            if(x1r_pseudo(p,n).gt.dfloat(ii))then
                x1r_pseudo(p,n) = x1r_pseudo(p,n) - dfloat(ii)
                x1r_pseudo_coord(p,n) = x1r_pseudo_coord(p,n)+1
            endif
            if(y1r_pseudo(p,n).lt.0.)then
              y1r_pseudo(p,n)=dfloat(jj)+y1r_pseudo(p,n)
              y1r_pseudo_coord(p,n) = y1r_pseudo_coord(p,n) -1
            endif
            if(y1r_pseudo(p,n).gt.dfloat(jj))then
              y1r_pseudo(p,n)=y1r_pseudo(p,n)-dfloat(jj)
              y1r_pseudo_coord(p,n)= y1r_pseudo_coord(p,n) +1
            endif

            if(x2r_pseudo(p,n).lt.0)then
                x2r_pseudo(p,n) = dfloat(ii)+x2r_pseudo(p,n)
                x2r_pseudo_coord(p,n)=x2r_pseudo_coord(p,n)-1
            end if
            if(x2r_pseudo(p,n).gt.dfloat(ii))then
                x2r_pseudo(p,n) = x2r_pseudo(p,n) - dfloat(ii)
                x2r_pseudo_coord(p,n) = x2r_pseudo_coord(p,n)+1
            endif
            if(y2r_pseudo(p,n).lt.0.)then
              y2r_pseudo(p,n)=dfloat(jj)+y2r_pseudo(p,n)
              y2r_pseudo_coord(p,n) = y2r_pseudo_coord(p,n) -1
            endif
            if(y2r_pseudo(p,n).gt.dfloat(jj))then
              y2r_pseudo(p,n)=y2r_pseudo(p,n)-dfloat(jj)
              y2r_pseudo_coord(p,n)= y2r_pseudo_coord(p,n) +1
            endif

        endif


            enddo
            enddo
            endif


c ----------------- EDDY INDUCED PARTICLE ADVECTION ---------------
        if (i_eddy.eq.1) then

            do p = 1,nbins
            do n = 1,npoints

        if (isolve.eq.0) then

        call rk4_bicubic(ii,jj,x1r_eddy(p,n),y1r_eddy(p,n),dt_nondim
     &   ,dfloat(0)
     & ,M1_eddy_old,M1_eddy_half,M1_eddy_new
     & ,x_diff1,y_diff1)
        call rk4_bicubic(ii,jj,x2r_eddy(p,n),y2r_eddy(p,n),dt_nondim
     &   ,dfloat(0)
     & ,M2_eddy_old,M2_eddy_half,M2_eddy_new
     & ,x_diff2,y_diff2)

        elseif (isolve.eq.1) then


        call rk4_2dcubic(ii,jj,x1r_eddy(p,n),y1r_eddy(p,n)
     &   ,dt_nondim,
     &   dfloat(0)
     & ,a1_eddy_old,b1_eddy_old,c1_eddy_old
     & ,d1_eddy_old,a1_eddy_half,b1_eddy_half,c1_eddy_half,d1_eddy_half
     & ,a1_eddy_new,b1_eddy_new,c1_eddy_new,d1_eddy_new
     & ,x_diff1,y_diff1)

        !print*, 'x1r_eddy,y1r_eddy'
        !print*, x1r_eddy(p,n), y1r_eddy(p,n)
        call rk4_2dcubic(ii,jj,x2r_eddy(p,n),y2r_eddy(p,n),dt_nondim
     &   ,dfloat(0)
     & ,a2_eddy_old,b2_eddy_old,c2_eddy_old
     & ,d2_eddy_old,a2_eddy_half,b2_eddy_half,c2_eddy_half,d2_eddy_half
     & ,a2_eddy_new,b2_eddy_new,c2_eddy_new,d2_eddy_new
     & ,x_diff2,y_diff2)

        endif

        x1r_eddy(p,n) = x1r_eddy(p,n) + x_diff1
        y1r_eddy(p,n) = y1r_eddy(p,n) + y_diff1
        x2r_eddy(p,n) = x2r_eddy(p,n) + x_diff2
        y2r_eddy(p,n) = y2r_eddy(p,n) + y_diff2

            if(x1r_eddy(p,n).lt.0)then
                x1r_eddy(p,n) = dfloat(ii)+x1r_eddy(p,n)
                x1r_eddy_coord(p,n)=x1r_eddy_coord(p,n)-1
            end if
            if(x1r_eddy(p,n).gt.dfloat(ii))then
                x1r_eddy(p,n) = x1r_eddy(p,n) - dfloat(ii)
                x1r_eddy_coord(p,n) = x1r_eddy_coord(p,n)+1
            endif
            if(y1r_eddy(p,n).lt.0.)then
              y1r_eddy(p,n)=dfloat(jj)+y1r_eddy(p,n)
              y1r_eddy_coord(p,n) = y1r_eddy_coord(p,n) -1
            endif
            if(y1r_eddy(p,n).gt.dfloat(jj))then
              y1r_eddy(p,n)=y1r_eddy(p,n)-dfloat(jj)
              y1r_eddy_coord(p,n)= y1r_eddy_coord(p,n) +1
            endif

            if(x2r_eddy(p,n).lt.0)then
                x2r_eddy(p,n) = dfloat(ii)+x2r_eddy(p,n)
                x2r_eddy_coord(p,n)=x2r_eddy_coord(p,n)-1
            end if
            if(x2r_eddy(p,n).gt.dfloat(ii))then
                x2r_eddy(p,n) = x2r_eddy(p,n) - dfloat(ii)
                x2r_eddy_coord(p,n) = x2r_eddy_coord(p,n)+1
            endif
            if(y2r_eddy(p,n).lt.0.)then
              y2r_eddy(p,n)=dfloat(jj)+y2r_eddy(p,n)
              y2r_eddy_coord(p,n) = y2r_eddy_coord(p,n) -1
            endif
            if(y2r_eddy(p,n).gt.dfloat(jj))then
              y2r_eddy(p,n)=y2r_eddy(p,n)-dfloat(jj)
              y2r_eddy_coord(p,n)= y2r_eddy_coord(p,n) +1
            endif


            enddo
            enddo
            endif




            if (k_s(k) .eq. l_tot_day*k_save) then

            k_s(k) = 0

            print *, 'Writing to NETCDF files at time',time_day



            if (i_full.eq.1) then

            do p = 1,nbins

            x1 = x1r(p,:)
            x2 = x2r(p,:)
            x3 = x3r(p,:)
            y1 = y1r(p,:)
            y2 = y2r(p,:)
            y3 = y3r(p,:)

            x1_coord = x1r_coord(p,:)
            x2_coord = x2r_coord(p,:)
            x3_coord = x3r_coord(p,:)
            y1_coord = y1r_coord(p,:)
            y2_coord = y2r_coord(p,:)
            y3_coord = y3r_coord(p,:)


            call write_binned_file(full_name,npoints,p,k
     &       ,x1,y1
     &       ,x2,y2
     &       ,x3,y3
     &       ,x1_coord,y1_coord,x2_coord
     &       ,y2_coord,x3_coord,y3_coord
     &       ,time_day,nrel(k))

            enddo

            endif

            nrel(k) = nrel(k) + 1

            endif


            endif
            enddo


        enddo

      end program offline_transport
