c***********************************************************************
c     Q-GCM Version 1.5.1 : last modified 24/02/2022
c***********************************************************************
c
c     Copyright 2022 Jeff Blundell, Andy Hogg and Bill Dewar.
c     This file is part of Q-GCM.
c
c     Q-GCM is free software: you can redistribute it and/or modify
c     it under the terms of the GNU General Public License as
c     published by the Free Software Foundation, either version 3
c     of the License, or (at your option) any later version.
c
c     Q-GCM is distributed in the hope that it will be useful,
c     but WITHOUT ANY WARRANTY; without even the implied warranty
c     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
c     See the GNU General Public License for more details.
c
c     You should have received a copy of the GNU General Public License
c     along with Q-GCM.  If not, see <http://www.gnu.org/licenses/>.
c
c***********************************************************************
c
      MODULE qocdiag

*     Contains subprograms for computing and outputting (in netCDF format)
*     the time-derivative of ocean vorticity qo, and its component terms.
*     Output interval is noutoc atmospheric timesteps = odiday days
*     (controlled by when qocdiag_out is called from the main program)
*     This version modified 24/02/2022 to additionally output, as qotd2p,
*     qotbet and qotifd, the components of J(qo, po) depending on Del-sqd(po),
*     the planetary vorticity gradient (the "beta" term), and the
*     interface displacements (including bottom topography) respectively.
*     The new terms are evaluated on those boundaries where they are
*     non-zero, albeit using a simpler finite difference formulation that
*     the full Arakawa Jacobian, which is used only at internal p-points.

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: qocdiag_init, qocdiag_out

      integer, PUBLIC, SAVE  :: qocncid
      integer, PRIVATE, SAVE :: dqdt_id, qotjac_id, qotd2p_id, qotbet_id,
     &                          qotifd_id, qt2dif_id, qt4dif_id,
     &                          qotent_id, timop_id
      logical, PRIVATE, SAVE :: qt2dif_out
*     Controls whether the 2nd order diffusive contribution is written
*     out, since this is usually exactly zero by choice of ah2oc

      CONTAINS

c***********************************************************************

      SUBROUTINE qocdiag_init (nsko)

*     Modules
      USE parameters
#ifndef atmos_only
#  ifdef qoc_diag
#    ifdef use_netcdf
      USE occonst
#    endif /* use_netcdf */
#  endif
#endif
      USE timinfo, ONLY : noutstepoc
      USE nc_subs, ONLY : handle_err

      IMPLICIT NONE

#ifdef use_netcdf
      INCLUDE 'netcdf.inc'
#endif /* use_netcdf */
*
*     Subroutine arguments
      integer, INTENT(IN) :: nsko

#ifndef atmos_only
#  ifdef qoc_diag
#    ifdef use_netcdf

*     Local parameters
      character (len=*), parameter :: subnam = 'qocdiag_init'

!!    netCDF variables used locally
      integer ncstat, timopdim, xopdim, yopdim, lodim
      integer axodims, podims(4), xop_id, yop_id, lo_id

!!    Other variables used locally
      double precision xo(nxpo), yo(nypo), zo(nlo)
      integer i, j, k, iwk, mwk

!!    Definition section: define dimensions and variables
!!    ===================================================
!!    Define four dimensions: x, y, z, time

*     Decide if any variables need not be defined
*     We usually run with ah2oc = 0.0 at all levels
*     If the value at any level is nonzero, then
*     output the contribution of this term to dq/dt
      qt2dif_out = .false.
      do k=1,nlo
        if ( ah2oc(k).ne.0.0d0 ) qt2dif_out = .true.
      enddo

!!    Dimension definitions for p-grid output file
!!    --------------------------------------------
      ncstat = nf_def_dim(qocncid, 'time', noutstepoc, timopdim)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!    x and y dimensions for the (subsampled) p-grid
      mwk = mod(nxpo,nsko)
      iwk = min(mwk,1) + (nxpo-mwk)/nsko
      ncstat = nf_def_dim(qocncid, 'xp', iwk, xopdim)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      mwk = mod(nypo,nsko)
      iwk = min(mwk,1) + (nypo-mwk)/nsko
      ncstat = nf_def_dim(qocncid, 'yp', iwk, yopdim)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!    And the z dimension
      ncstat = nf_def_dim(qocncid, 'z', nlo, lodim)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!    Grid variable definitions for p-grid files
!!    ------------------------------------------
!!    Define a one-dimensional variable called
!!    'xp' which stores x gridpoints (p-grid)
      axodims = xopdim
      ncstat = nf_def_var(qocncid, 'xp', NF_FLOAT, 1, axodims, xop_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(qocncid, xop_id, 'units', 2, 'km')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(qocncid, xop_id, 'long_name',
     &            21, 'Ocean X axis (p-grid)')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!    Define a one-dimensional variable called
!!    'yp' which stores y gridpoints (p-grid)
      axodims = yopdim
      ncstat = nf_def_var(qocncid, 'yp', NF_FLOAT, 1, axodims, yop_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(qocncid, yop_id, 'units', 2, 'km')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(qocncid, yop_id, 'long_name',
     &            21, 'Ocean Y axis (p-grid)')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!    Define a one-dimensional variable called
!!    'time' which stores time data in years
      axodims = timopdim
      ncstat=nf_def_var(qocncid, 'time', NF_FLOAT, 1, axodims, timop_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(qocncid, timop_id, 'units', 5, 'years')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat=nf_put_att_text(qocncid,timop_id,'long_name',9,'Time axis')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!    Define a 1-dimensional variable called 'z'
!!    to store mid-layer depths (in km)
      axodims = lodim
      ncstat = nf_def_var(qocncid, 'z', NF_FLOAT, 1, axodims, lo_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(qocncid, lo_id, 'units', 2, 'km')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(qocncid, lo_id, 'long_name',
     &            26, 'Ocean mid-layer depth axis')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!    Define a 4-dimensional variable called 'dqdt'
!!    to store the time derivative of qo
      podims(1) = xopdim
      podims(2) = yopdim
      podims(3) = lodim
      podims(4) = timopdim
      ncstat = nf_def_var(qocncid, 'dqdt', NF_FLOAT, 4, podims, dqdt_id)
      if ( ncstat.ne.NF_NOERR ) then
        print *,' problem defining dqdt in ocean q-diag file'
        call handle_err (ncstat, subnam)
      endif
      ncstat = nf_put_att_text(qocncid, dqdt_id, 'units', 4, 's^-2')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(qocncid, dqdt_id, 'long_name',
     &            23, 'd/dt(oceanic vorticity)')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!    Define a 4-dimensional variable called 'qotjac'
!!    to store advective (Jacobian) component of d(qo)/dt
      podims(1) = xopdim
      podims(2) = yopdim
      podims(3) = lodim
      podims(4) = timopdim
      ncstat = nf_def_var(qocncid, 'qotjac', NF_FLOAT,
     &                    4, podims, qotjac_id)
      if ( ncstat.ne.NF_NOERR ) then
        print *,' problem defining qotjac in ocean q-diag file'
        call handle_err (ncstat, subnam)
      endif
      ncstat = nf_put_att_text(qocncid, qotjac_id, 'units', 4, 's^-2')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(qocncid, qotjac_id, 'long_name',
     &            31, 'advective component of d(qo)/dt')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!    Define a 4-dimensional variable called 'qotd2p' to store the
!!    Del-sqd(po) component of the advective (Jacobian) term in d(qo)/dt
      podims(1) = xopdim
      podims(2) = yopdim
      podims(3) = lodim
      podims(4) = timopdim
      ncstat = nf_def_var(qocncid, 'qotd2p', NF_FLOAT,
     &                    4, podims, qotd2p_id)
      if ( ncstat.ne.NF_NOERR ) then
        print *,' problem defining qotd2p in ocean q-diag file'
        call handle_err (ncstat, subnam)
      endif
      ncstat = nf_put_att_text(qocncid, qotd2p_id, 'units', 4, 's^-2')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(qocncid, qotd2p_id, 'long_name',
     &            33, 'Del-sqd(po) component of Jacobian')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!    Define a 4-dimensional variable called 'qotbet' to store the
!!    beta component of the advective (Jacobian) term in d(qo)/dt
      podims(1) = xopdim
      podims(2) = yopdim
      podims(3) = lodim
      podims(4) = timopdim
      ncstat = nf_def_var(qocncid, 'qotbet', NF_FLOAT,
     &                    4, podims, qotbet_id)
      if ( ncstat.ne.NF_NOERR ) then
        print *,' problem defining qotbet in ocean q-diag file'
        call handle_err (ncstat, subnam)
      endif
      ncstat = nf_put_att_text(qocncid, qotbet_id, 'units', 4, 's^-2')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(qocncid, qotbet_id, 'long_name',
     &            28, '"beta" component of Jacobian')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!    Define a 4-dimensional variable called 'qotifd' to store the interface
!!    displacement component of the advective (Jacobian) term in d(qo)/dt
      podims(1) = xopdim
      podims(2) = yopdim
      podims(3) = lodim
      podims(4) = timopdim
      ncstat = nf_def_var(qocncid, 'qotifd', NF_FLOAT,
     &                    4, podims, qotifd_id)
      if ( ncstat.ne.NF_NOERR ) then
        print *,' problem defining qotifd in ocean q-diag file'
        call handle_err (ncstat, subnam)
      endif
      ncstat = nf_put_att_text(qocncid, qotifd_id, 'units', 4, 's^-2')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(qocncid, qotifd_id, 'long_name',
     &            44, 'interface displacement component of Jacobian')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!    Define a 4-dimensional variable called 'qt2dif'
!!    to store dissipation (2nd order) component of d(qo)/dt
      if ( qt2dif_out ) then
        podims(1) = xopdim
        podims(2) = yopdim
        podims(3) = lodim
        podims(4) = timopdim
        ncstat = nf_def_var(qocncid, 'qt2dif', NF_FLOAT,
     &                      4, podims, qt2dif_id)
        if ( ncstat.ne.NF_NOERR ) then
          print *,' problem defining qt2dif in ocean q-diag file'
          call handle_err (ncstat, subnam)
        endif
        ncstat = nf_put_att_text(qocncid, qt2dif_id, 'units', 4, 's^-2')
        if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
        ncstat = nf_put_att_text(qocncid, qt2dif_id, 'long_name',
     &              45, 'dissipative component of d(qo)/dt (2nd order)')
        if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      endif

!!    Define a 4-dimensional variable called 'qt4dif'
!!    to store dissipation (4th order) component of d(qo)/dt
      podims(1) = xopdim
      podims(2) = yopdim
      podims(3) = lodim
      podims(4) = timopdim
      ncstat = nf_def_var(qocncid, 'qt4dif', NF_FLOAT,
     &                    4, podims, qt4dif_id)
      if ( ncstat.ne.NF_NOERR ) then
        print *,' problem defining qt4dif in ocean q-diag file'
        call handle_err (ncstat, subnam)
      endif
      ncstat = nf_put_att_text(qocncid, qt4dif_id, 'units', 4, 's^-2')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(qocncid, qt4dif_id, 'long_name',
     &            45, 'dissipative component of d(qo)/dt (4th order)')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!    Define a 4-dimensional variable called 'qotent'
!!    to store entrainment (forcing) component of d(qo)/dt
      podims(1) = xopdim
      podims(2) = yopdim
      podims(3) = lodim
      podims(4) = timopdim
      ncstat = nf_def_var(qocncid, 'qotent', NF_FLOAT,
     &                    4, podims, qotent_id)
      if ( ncstat.ne.NF_NOERR ) then
        print *,' problem defining qotent in ocean q-diag file'
        call handle_err (ncstat, subnam)
      endif
      ncstat = nf_put_att_text(qocncid, qotent_id, 'units', 4, 's^-2')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(qocncid, qotent_id, 'long_name',
     &            43, 'entrainment (forcing) component of d(qo)/dt')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!    Leave definition mode: entering data mode
!!    -----------------------------------------
      ncstat = nf_enddef(qocncid)
      if ( ncstat.ne.NF_NOERR ) then
        print *,' qocdiag_init: enddef problem for ocean q-diag file'
        print *,' if problem is "One or more variable sizes'
        print *,' violate format constraints", either coarsen'
        print *,' spatial and/or temporal sampling of the data,'
        print *,' or create the file with "64-bit offset" format'
        call handle_err (ncstat, subnam)
      endif

!!    Calculate x gridpoints and store in 'x' arrays
!!    p-grid points
      mwk = mod(nxpo,nsko)
      iwk = min(mwk,1) + (nxpo-mwk)/nsko
      do i=1,iwk
        xo(i) = 1.0d-3*( xpo(1+(i-1)*nsko) - xpo(1) )
      enddo
      ncstat = nf_put_vara_double(qocncid, xop_id, 1, iwk, xo)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!    Calculate y gridpoints and store in 'y' arrays
!!    p-grid points
      mwk = mod(nypo,nsko)
      iwk = min(mwk,1) + (nypo-mwk)/nsko
      do j=1,iwk
        yo(j) = 1.0d-3*( ypo(1+(j-1)*nsko) - ypo(1) )
      enddo
      ncstat = nf_put_vara_double(qocncid, yop_id, 1, iwk, yo)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!    Convert mid-layer depths into km and store in 'z'
      zo(1) = 0.5d-3*hoc(1)
      do k=2,nlo
        zo(k) = zo(k-1) + 0.5d-3*( hoc(k-1) + hoc(k) )
      enddo
      ncstat = nf_put_var_double(qocncid, lo_id, zo)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
!!!!! put operation seems to lack some arguments!!!
!!!!! This may be because it's a put_var not put_vara!!!

#    endif /* use_netcdf */
#  endif /* qoc_diag */
#endif /* not atmos_only */

      END SUBROUTINE qocdiag_init

c***********************************************************************

      SUBROUTINE qocdiag_out (nsko)

*     Compute oceanic dq/dt and its component terms. The code to do this
*     is borrowed from subroutines qgostep and ocadif in qgosubs.F.
*     An alternative approach would be to be to save the terms in
*     those routines, but this would (a) complicate them, and (b)
*     entail declaring much more storage for transfer. Here the arrays are
*     invoked only when necessary, at the expense of some recomputation.
*     This version also computes some of the components of the Jacobian term.

*     Modules
      USE parameters
#ifndef atmos_only
#  ifdef qoc_diag
#    ifdef use_netcdf
      USE occonst
      USE ocstate
#    endif /* use_netcdf */
#  endif
#endif
      USE timinfo, ONLY : ntdone, noutoc, tyrs
      USE nc_subs, ONLY : handle_err

      IMPLICIT NONE

#ifdef use_netcdf
      INCLUDE 'netcdf.inc'
#endif /* use_netcdf */
*
*     Subroutine arguments
      integer, INTENT(IN) :: nsko

#ifndef atmos_only
#  ifdef qoc_diag
#    ifdef use_netcdf

*     Local parameters
      character (len=*), parameter :: subnam = 'qocdiag_out'

!!    netCDF variables used locally
!!    start4, count4 are for 4-D variables
      integer ncstat,startt,countt,start4(4),count4(4)

!!    Other variables used locally
      integer i,j,k,ipwk,jpwk,mwk
      double precision adfaco,bcfaco,fohfac(nlo),bdrfac,ah2fac,
     &                 ah4fac,d2pm(nxpo,nypo),del4p(nxpo,nypo),
     &                 qt2dif(nxpo,nypo),qt4dif(nxpo,nypo),d6p,
     &                 qotjac(nxpo,nypo),qotd2p(nxpo,nypo),
     &                 qotbet(nxpo,nypo),qotifd(nxpo,nypo),
     &                 qotent(nxpo,nypo),d2pc(nxpo,nypo),
     &                 dqdt(nxpo,nypo),qifd(nxpo,nypo),wrk(nxpo*nypo),
     &                 rdto,pnormf,pnorm, dtangf,pxof0,pyof0,f0Am,f0Ac,f0Ap

!!    Store current time as part of 'time' vector
      startt = ntdone/noutoc + 1
      countt = 1
      ncstat=nf_put_vara_double(qocncid, timop_id, startt, countt, tyrs)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!    Compute subsampling array indices
      mwk = mod(nxpo, nsko)
      ipwk = min(mwk,1) + (nxpo-mwk)/nsko
      mwk = mod(nypo, nsko)
      jpwk = min(mwk,1) + (nypo-mwk)/nsko

*     Compute vorticity time derivative and its components with code from qgosubs.F
*     -----------------------------------------------------------------------------
*     This version also computes the Del-sqd(po), "beta" and
*     interface displacement components of the Jacobian term.
      adfaco = 1.0d0/(12.0d0*dxo*dyo*fnot)
*     Version with nondimensional bccooc
      bcfaco = bccooc*dxom2/( 0.5d0*bccooc + 1.0d0 )
      pnorm = 1.0d0/( dxo*( 0.5d0*bccooc + 1.0d0 ) )
      pnormf = 1.0d0/( fnot*dxo*( 0.5d0*bccooc + 1.0d0 ) )
      dtangf = 0.5d0/( fnot*dxo )
      do k=1,nlo
        fohfac(k) = fnot/hoc(k)
      enddo
      bdrfac = 0.5d0*sign(1.0d0,fnot)*delek/hoc(nlo)
      rdto = 1.0d0/dto

      do k=1,nlo

        ah2fac = ah2oc(k)/fnot
        ah4fac = ah4oc(k)/fnot
*       Set coefficients for interface displacement terms
        if ( k.eq.1 ) then
*         Top layer
          f0Am = 0.0d0
          f0Ac = fnot*amatoc(1,1)
          f0Ap = fnot*amatoc(1,2)
         else if ( k.eq.nlo ) then
*         Bottom layer k = nlo (includes topography)
          f0Am = fnot*amatoc(nlo,nlo-1)
          f0Ac = fnot*amatoc(nlo, nlo )
          f0Ap = 0.0d0
         else
*         Intermediate layers
          f0Am = fnot*amatoc(k,k-1)
          f0Ac = fnot*amatoc(k, k )
          f0Ap = fnot*amatoc(k,k+1)
        endif

*       Set netCDF output counters once only, since they
*       are the same for all fields in any given layer
        start4(1) = 1
        start4(2) = 1
        start4(3) = k
        start4(4) = ntdone/noutoc + 1
        count4(1) = ipwk
        count4(2) = jpwk
        count4(3) = 1
        count4(4) = 1

!$OMP   PARALLEL DEFAULT (NONE)
!$OMP&           PRIVATE (i,j,d6p,pxof0,pyof0)
!$OMP&           SHARED  (k,d2pm,d2pc,bcfaco,pom,dxom2,del4p,qt2dif,
!$OMP&                    qt4dif,qotjac,qotd2p,qotbet,qotifd,qotent,dqdt,
!$OMP&                    po,qo,fohfac,wekpo,yporel,entoc,bdrfac,ah2fac,
!$OMP&                    ah4fac,start4,count4,adfaco,rdto,ipwk,jpwk,
!$OMP&                    wrk,nsko,ncstat,qocncid,dqdt_id,qotjac_id,
!$OMP&                    qotd2p_id,qotbet_id,qotifd_id,qt2dif_id,
!$OMP&                    qt4dif_id,qotent_id,qom,qt2dif_out,pnormf,
!$OMP&                    dtangf,qifd,f0Am,f0Ac,f0Ap,ddynoc,hdxom1)

*       Compute Del-sqd(p) at previous time level (d2pm) for dissipation,
*       and at current time level (d2pc) for Del-sqd(po) contribution to
*       the advective (Jacobian) term in the dq/dt equation.
*       -----------------------------------------------------------------
*       N & S boundaries (including corner points) - mixed BCs
*       No x-derivative term because p is constant along boundary
!$OMP   DO SCHEDULE (STATIC)
        do i=1,nxpo
          d2pm(i,  1 ) = bcfaco*( pom(i,   2  ,k) - pom(i,  1 ,k) )
          d2pm(i,nypo) = bcfaco*( pom(i,nypo-1,k) - pom(i,nypo,k) )
          d2pc(i,  1 ) = bcfaco*( po(i,   2  ,k) - po(i,  1 ,k) )
          d2pc(i,nypo) = bcfaco*( po(i,nypo-1,k) - po(i,nypo,k) )
        enddo
!$OMP   END DO NOWAIT

!$OMP   DO SCHEDULE (STATIC)
        do j=2,nypo-1
*         Western boundary (i=1)
#  ifdef cyclic_ocean
*         Zonally cyclic ocean
          d2pm(1,j) = (  pom(1,j-1,k) + pom(nxpo-1,j,k) + pom(2,j,k)
     &                 + pom(1,j+1,k) - 4.0d0*pom(1,j,k) )*dxom2
          d2pc(1,j) = (  po(1,j-1,k) + po(nxpo-1,j,k) + po(2,j,k)
     &                 + po(1,j+1,k) - 4.0d0*po(1,j,k) )*dxom2
#  else
*         Finite box ocean - mixed BCs
          d2pm( 1 ,j) = bcfaco*( pom(  2  ,j,k) - pom( 1 ,j,k) )
          d2pc( 1 ,j) = bcfaco*( po(  2  ,j,k) - po( 1 ,j,k) )
#  endif /* cyclic_ocean */
*         Inner points (standard case)
          do i=2,nxpo-1
            d2pm(i,j) = (  pom(i,j-1,k) + pom(i-1,j,k) + pom(i+1,j,k)
     &                   + pom(i,j+1,k) - 4.0d0*pom(i,j,k) )*dxom2
            d2pc(i,j) = (  po(i,j-1,k) + po(i-1,j,k) + po(i+1,j,k)
     &                   + po(i,j+1,k) - 4.0d0*po(i,j,k) )*dxom2
          enddo
*         Eastern boundary (i=nxpo)
#  ifdef cyclic_ocean
*         Zonally cyclic ocean
          d2pm(nxpo,j) = d2pm(1,j)
          d2pc(nxpo,j) = d2pc(1,j)
#  else
*         Finite box ocean - mixed BCs
          d2pm(nxpo,j) = bcfaco*( pom(nxpo-1,j,k) - pom(nxpo,j,k) )
          d2pc(nxpo,j) = bcfaco*( po(nxpo-1,j,k) - po(nxpo,j,k) )
#  endif /* cyclic_ocean */
        enddo
!$OMP   END DO

*       Compute Del-4th(p) at previous time level for dissipation
*       ---------------------------------------------------------
*       N & S boundaries (including corners) - mixed BCs
*       No x-derivative term because p is constant along boundary
!$OMP   DO SCHEDULE (STATIC)
        do i=1,nxpo
          del4p(i,  1 ) = bcfaco*( d2pm(i,   2  ) - d2pm(i,  1 ) )
          del4p(i,nypo) = bcfaco*( d2pm(i,nypo-1) - d2pm(i,nypo) )
        enddo
!$OMP   END DO NOWAIT

!$OMP   DO SCHEDULE (STATIC)
        do j=2,nypo-1
*         Western boundary (i=1)
#  ifdef cyclic_ocean
*         Zonally cyclic ocean
          del4p(1,j) = (  d2pm(1,j-1) + d2pm(nxpo-1,j) + d2pm(2,j)
     &                  + d2pm(1,j+1) - 4.0d0*d2pm(1,j) )*dxom2
#  else
*         Finite box ocean - mixed BCs
          del4p( 1 ,j) = bcfaco*( d2pm(  2  ,j) - d2pm( 1 ,j) )
#  endif /* cyclic_ocean */
*         Inner points (standard case)
          do i=2,nxpo-1
            del4p(i,j) = (  d2pm(i,j-1) + d2pm(i-1,j) + d2pm(i+1,j)
     &                    + d2pm(i,j+1) - 4.0d0*d2pm(i,j) )*dxom2
          enddo
*         Eastern boundary (i=nxpo)
#  ifdef cyclic_ocean
*         Zonally cyclic ocean
          del4p(nxpo,j) = del4p(1,j)
#  else
*         Finite box ocean - mixed BCs
          del4p(nxpo,j) = bcfaco*( d2pm(nxpo-1,j) - d2pm(nxpo,j) )
#  endif /* cyclic_ocean */
        enddo
!$OMP   END DO NOWAIT

*       Compute interface displacement component of qo
*       ----------------------------------------------
*       Code is borrowed directly from subroutine qcomp in vorsubs.F
*       This term can be computed at all p-points, including boundaries and corners
*
        if ( k.eq.1 ) then
*         Top layer (k=1)
!$OMP     DO SCHEDULE (STATIC)
          do j=1,nypo
            do i=1,nxpo
              qifd(i,j) = - ( f0Ac*po(i,j,1) + f0Ap*po(i,j,2) )
            enddo
          enddo
!$OMP     END DO NOWAIT
         else if ( k.eq.nlo ) then
*         Bottom layer k = nlo (includes topography)
!$OMP     DO SCHEDULE (STATIC)
          do j=1,nypo
            do i=1,nxpo
              qifd(i,j) = - ( f0Am*po(i,j,nlo-1) + f0Ac*po(i,j,nlo) )
     &                    + ddynoc(i,j)
            enddo
          enddo
!$OMP     END DO NOWAIT
         else
*         Intermediate layers (k=2,nlo-1)
!$OMP     DO SCHEDULE (STATIC)
          do j=1,nypo
            do i=1,nxpo
              qifd(i,j) = - (  f0Am*po(i,j,k-1) + f0Ac*po(i,j,k)
     &                       + f0Ap*po(i,j,k+1) )
            enddo
          enddo
!$OMP     END DO NOWAIT
        endif

*       Zero all the results arrays as a precaution
*       -------------------------------------------
*       This can be removed when the boundary terms are computed
!$OMP   DO SCHEDULE (STATIC)
        do j=1,nypo
          do i=1,nxpo
            qt2dif(i,j) = 0.0d0
            qt4dif(i,j) = 0.0d0
            qotjac(i,j) = 0.0d0
            qotd2p(i,j) = 0.0d0
            qotbet(i,j) = 0.0d0
            qotifd(i,j) = 0.0d0
            qotent(i,j) = 0.0d0
            dqdt(i,j)   = 0.0d0
          enddo
        enddo
!$OMP   END DO

*       Compute advective (Jacobian) and diffusive contributions to dq/dt
*       -----------------------------------------------------------------
*       Compute Del-6th(p) except at N & S boundaries
*       No need to apply N & S boundary conditions
*       Also compute the Del-sqd(po) and "beta" term components of the Jacobian

!$OMP   DO SCHEDULE (STATIC)
        do j=2,nypo-1
*         Western boundary (i=1)
#  ifdef cyclic_ocean
*         Zonally cyclic ocean
          d6p = dxom2*(  del4p(1,j-1) + del4p(nxpo-1,j) + del4p(2,j)
     &                 + del4p(1,j+1) - 4.0d0*del4p(1,j) )
          qt2dif(1,j) = ah2fac*del4p(1,j)
          qt4dif(1,j) = -ah4fac*d6p
          qotjac(1,j) = adfaco*
     & ( (qo( 2 , j ,k)-qo(nxpo-1, j ,k))*(po( 1 ,j+1,k)-po(  1 ,j-1,k))
     &  +(qo( 1 ,j-1,k)-qo(  1 ,j+1,k))*(po( 2 , j ,k)-po(nxpo-1, j ,k))
     &  + qo(   2  , j ,k)*( po(   2  ,j+1,k) - po(   2  ,j-1,k) )
     &  - qo(nxpo-1, j ,k)*( po(nxpo-1,j+1,k) - po(nxpo-1,j-1,k) )
     &  - qo(   1  ,j+1,k)*( po(   2  ,j+1,k) - po(nxpo-1,j+1,k) )
     &  + qo(   1  ,j-1,k)*( po(   2  ,j-1,k) - po(nxpo-1,j-1,k) )
     &  + po(   1  ,j+1,k)*( qo(   2  ,j+1,k) - qo(nxpo-1,j+1,k) )
     &  - po(   1  ,j-1,k)*( qo(   2  ,j-1,k) - qo(nxpo-1,j-1,k) )
     &  - po(   2  , j ,k)*( qo(   2  ,j+1,k) - qo(   2  ,j-1,k) )
     &  + po(nxpo-1, j ,k)*( qo(nxpo-1,j+1,k) - qo(nxpo-1,j-1,k) ) )

*         Compute the (1/f0)*Del-sqd(po) contribution to the Jacobian term
          qotd2p(1,j) = (adfaco/fnot)*
     & ( (d2pc( 2 , j )-d2pc(nxpo-1, j ))*(po( 1 ,j+1,k)-po(  1 ,j-1,k))
     &  +(d2pc( 1 ,j-1)-d2pc(  1 ,j+1))*(po( 2 , j ,k)-po(nxpo-1, j ,k))
     &  + d2pc(   2  , j )*( po(   2  ,j+1,k) - po(   2  ,j-1,k) )
     &  - d2pc(nxpo-1, j )*( po(nxpo-1,j+1,k) - po(nxpo-1,j-1,k) )
     &  - d2pc(   1  ,j+1)*( po(   2  ,j+1,k) - po(nxpo-1,j+1,k) )
     &  + d2pc(   1  ,j-1)*( po(   2  ,j-1,k) - po(nxpo-1,j-1,k) )
     &  + po(   1  ,j+1,k)*( d2pc(   2  ,j+1) - d2pc(nxpo-1,j+1) )
     &  - po(   1  ,j-1,k)*( d2pc(   2  ,j-1) - d2pc(nxpo-1,j-1) )
     &  - po(   2  , j ,k)*( d2pc(   2  ,j+1) - d2pc(   2  ,j-1) )
     &  + po(nxpo-1, j ,k)*( d2pc(nxpo-1,j+1) - d2pc(nxpo-1,j-1) ) )

*         Compute "beta" term contribution to the Jacobian term
          qotbet(1,j) = adfaco*beta*
     & ( (yporel(j-1)-yporel(j+1))*(po( 2 , j ,k)-po(nxpo-1, j ,k))
     &  + yporel( j )*( po(   2  ,j+1,k) - po(   2  ,j-1,k) )
     &  - yporel( j )*( po(nxpo-1,j+1,k) - po(nxpo-1,j-1,k) )
     &  - yporel(j+1)*( po(   2  ,j+1,k) - po(nxpo-1,j+1,k) )
     &  + yporel(j-1)*( po(   2  ,j-1,k) - po(nxpo-1,j-1,k) )
     &  - po(   2  , j ,k)*( yporel(j+1) - yporel(j-1) )
     &  + po(nxpo-1, j ,k)*( yporel(j+1) - yporel(j-1) ) )

*         Compute interface displacement contribution to the Jacobian term
          qotifd(1,j) = adfaco*
     & ( (qifd( 2 , j )-qifd(nxpo-1, j ))*(po( 1 ,j+1,k)-po(  1 ,j-1,k))
     &  +(qifd( 1 ,j-1)-qifd(  1 ,j+1))*(po( 2 , j ,k)-po(nxpo-1, j ,k))
     &  + qifd(   2  , j )*( po(   2  ,j+1,k) - po(   2  ,j-1,k) )
     &  - qifd(nxpo-1, j )*( po(nxpo-1,j+1,k) - po(nxpo-1,j-1,k) )
     &  - qifd(   1  ,j+1)*( po(   2  ,j+1,k) - po(nxpo-1,j+1,k) )
     &  + qifd(   1  ,j-1)*( po(   2  ,j-1,k) - po(nxpo-1,j-1,k) )
     &  + po(   1  ,j+1,k)*( qifd(   2  ,j+1) - qifd(nxpo-1,j+1) )
     &  - po(   1  ,j-1,k)*( qifd(   2  ,j-1) - qifd(nxpo-1,j-1) )
     &  - po(   2  , j ,k)*( qifd(   2  ,j+1) - qifd(   2  ,j-1) )
     &  + po(nxpo-1, j ,k)*( qifd(nxpo-1,j+1) - qifd(nxpo-1,j-1) ) )

*         (layer-dependent) entrainment/forcing term
          if ( k.eq.1 ) then
            qotent(1,j) = fohfac(1)*( wekpo(1,j) - entoc(1,j) )
           else if ( k.eq.2 ) then
            qotent(1,j) = fohfac(2)*entoc(1,j)
           else
            qotent(1,j) = 0.0d0
          endif
*         Add bottom drag contribution
          if ( k.eq.nlo ) then
            qotent(1,j) = qotent(1,j) - bdrfac*d2pm(1,j)
          endif
          dqdt(1,j) =  qotjac(1,j) + qt2dif(1,j)
     &               + qt4dif(1,j) + qotent(1,j)
#  else
*         Finite box ocean
          qt2dif(1,j) = 0.0d0
          qt4dif(1,j) = 0.0d0
          qotjac(1,j) = 0.0d0
          qotent(1,j) = 0.0d0
*         Compute scaled normal derivative (1/f0)*d(po)/dx
*         using the mixed boundary condition (Appendix C).
          pxof0 = pnormf*( po(2,j,k) - po(1,j,k) )
*         On W & E boundaries, Del-sqd term reduces to - (1/f0^2)*d(pxx)/dy*d(po)/dx
          qotd2p(1,j) = 0.5d0*(adfaco/fnot)*
     &   ((d2pc(2,j)-d2pc(1,j))*(po(2,j+1,k)-po(2,j-1,k))
     &   +((d2pc(1,j-1)-d2pc(1,j+1)+d2pc(2,j-1)-
     &      d2pc(2,j+1)))*(po(2,j,k)-po(1,j,k))
     &    + d2pc(2, j)*( po(2,j+1,k)-po(2,j-1,k) )
C     &    - qo(1, j)*( po(1,j+1,k) - po(1,j-1,k) )
     &    - (d2pc(1,j+1)+d2pc(2,j+1))*( po(2,j+1,k) - po(1,j+1,k) )
     &    + (d2pc(1,j-1)+d2pc(2,j-1))*( po(2,j-1,k) - po(1,j-1,k) )
     &    + (po(1,j+1,k)+po(2,j+1,k))*( d2pc(2,j+1) - d2pc(1,j+1) )
     &    - (po(1,j-1,k)+po(2,j-1,k))*( d2pc(2,j-1) - d2pc(1,j-1) )
     &    - po(2, j ,k)*( d2pc(2,j+1) - d2pc(2,j-1) )
     &    + po(1, j ,k)*( d2pc(1,j+1) - d2pc(1,j-1) ) )

          qotifd(1,j) = 0.5d0*adfaco*
     &   ((qifd(2,j)-qifd(1,j))*(po(2,j+1,k)-po(2,j-1,k))
     &   +((qifd(1,j-1)-qifd(1,j+1)+qifd(2,j-1)-
     &      qifd(2,j+1)))*(po(2,j,k)-po(1,j,k))
     &    + qifd(2, j)*( po(2,j+1,k)-po(2,j-1,k) )
C     &    - qo(1, j)*( po(1,j+1,k) - po(1,j-1,k) )
     &    - (qifd(1,j+1)+qifd(2,j+1))*( po(2,j+1,k) - po(1,j+1,k) )
     &    + (qifd(1,j-1)+qifd(2,j-1))*( po(2,j-1,k) - po(1,j-1,k) )
     &    + (po(1,j+1,k)+po(2,j+1,k))*( qifd(2,j+1) - qifd(1,j+1) )
     &    - (po(1,j-1,k)+po(2,j-1,k))*( qifd(2,j-1) - qifd(1,j-1) )
     &    - po(2, j ,k)*( qifd(2,j+1) - qifd(2,j-1) )
     &    + po(1, j ,k)*( qifd(1,j+1) - qifd(1,j-1) ) )

          qotjac(1,j) = 0.5d0*adfaco*
     &   ( (qo(2, j ,k)-qo(1, j ,k))*(po(2,j+1,k)-po(2,j-1,k))
     &    +(qo(1,j-1,k)-qo(1,j+1,k) + qo(2,j-1,k)-qo(2,j+1,k))*
     &     (po(2, j ,k)-po(1, j ,k))
     &    +  qo(2, j ,k)*(  po(2,j+1,k) - po(2,j-1,k))
C     &    - qo(1, j ,k)*(  po(1,j+1,k) - po(1,j-1,k))
     &    - (qo(1,j+1,k)+qo(2,j+1,k))*( po(2,j+1,k) - po(1,j+1,k) )
     &    + (qo(1,j-1,k)+qo(2,j-1,k))*( po(2,j-1,k) - po(1,j-1,k) )
     &    + (po(1,j+1,k)+po(2,j+1,k))*( qo(2,j+1,k) - qo(1,j+1,k) )
     &    - (po(1,j-1,k)+po(2,j-1,k))*( qo(2,j-1,k) - qo(1,j-1,k) )
     &    -  po(2, j ,k)*( qo(2,j+1,k) - qo(2,j-1,k) )
     &    +  po(1, j ,k)*( qo(1,j+1,k) - qo(1,j-1,k) ) )

*         On W & E boundaries, "beta" term = -beta*(px/f0)
          qotbet(1,j) = 0.5d0*adfaco*beta*
c     &   ((qo(2,j,k)-qo(1,j,k))*(po(2,j+1,k)-po(2,j-1,k))
     &   (  ((yporel(j-1)-yporel(j+1)+yporel(j-1)-yporel(j+1)))*
     &       (po(2,j,k)-po(1,j,k))
     &    +   yporel(j)*( po(2,j+1,k) - po(2,j-1,k) )
C     &    -  yporel(j)*( po(1,j+1,k) - po(1,j-1,k) )
     &     - (yporel(j+1)+yporel(j+1))*( po(2,j+1,k) - po(1,j+1,k) )
     &     + (yporel(j-1)+yporel(j-1))*( po(2,j-1,k) - po(1,j-1,k) )
c     &     + (po(1,j+1,k)+po(2,j+1,k))*( yporel(j+1) - yporel(j+1) )
c     &     - (po(1,j-1,k)+po(2,j-1,k))*( yporel(j-1) - yporel(j-1) )
     &     -  po(2, j ,k)*( yporel(j+1) - yporel(j-1) )
     &     +  po(1, j ,k)*( yporel(j+1) - yporel(j-1) ) )

          dqdt(1,j)   = rdto*( qo(1,j,k) - qom(1,j,k) )
#  endif /* cyclic_ocean */
*         Inner points (standard case)
          do i=2,nxpo-1
            d6p = dxom2*(  del4p(i,j-1) + del4p(i-1,j) + del4p(i+1,j)
     &                   + del4p(i,j+1) - 4.0d0*del4p(i,j) )
            qt2dif(i,j) = ah2fac*del4p(i,j)
            qt4dif(i,j) = -ah4fac*d6p
            qotjac(i,j) = adfaco*
     &   ( (qo(i+1, j ,k)-qo(i-1, j ,k))*(po( i ,j+1,k)-po( i ,j-1,k))
     &    +(qo( i ,j-1,k)-qo( i ,j+1,k))*(po(i+1, j ,k)-po(i-1, j ,k))
     &    + qo(i+1, j ,k)*( po(i+1,j+1,k) - po(i+1,j-1,k) )
     &    - qo(i-1, j ,k)*( po(i-1,j+1,k) - po(i-1,j-1,k) )
     &    - qo( i ,j+1,k)*( po(i+1,j+1,k) - po(i-1,j+1,k) )
     &    + qo( i ,j-1,k)*( po(i+1,j-1,k) - po(i-1,j-1,k) )
     &    + po( i ,j+1,k)*( qo(i+1,j+1,k) - qo(i-1,j+1,k) )
     &    - po( i ,j-1,k)*( qo(i+1,j-1,k) - qo(i-1,j-1,k) )
     &    - po(i+1, j ,k)*( qo(i+1,j+1,k) - qo(i+1,j-1,k) )
     &    + po(i-1, j ,k)*( qo(i-1,j+1,k) - qo(i-1,j-1,k) ) )

*           Compute the (1/f0)*Del-sqd(po) contribution to the Jacobian term
            qotd2p(i,j) = (adfaco/fnot)*
     &   ( (d2pc(i+1, j )-d2pc(i-1, j ))*(po( i ,j+1,k)-po( i ,j-1,k))
     &    +(d2pc( i ,j-1)-d2pc( i ,j+1))*(po(i+1, j ,k)-po(i-1, j ,k))
     &    + d2pc(i+1, j )*( po(i+1,j+1,k) - po(i+1,j-1,k) )
     &    - d2pc(i-1, j )*( po(i-1,j+1,k) - po(i-1,j-1,k) )
     &    - d2pc( i ,j+1)*( po(i+1,j+1,k) - po(i-1,j+1,k) )
     &    + d2pc( i ,j-1)*( po(i+1,j-1,k) - po(i-1,j-1,k) )
     &    + po( i ,j+1,k)*( d2pc(i+1,j+1) - d2pc(i-1,j+1) )
     &    - po( i ,j-1,k)*( d2pc(i+1,j-1) - d2pc(i-1,j-1) )
     &    - po(i+1, j ,k)*( d2pc(i+1,j+1) - d2pc(i+1,j-1) )
     &    + po(i-1, j ,k)*( d2pc(i-1,j+1) - d2pc(i-1,j-1) ) )

*           Compute "beta" term contribution to the Jacobian term
            qotbet(i,j) = adfaco*beta*
     &   ( (yporel(j-1)-yporel(j+1))*(po(i+1, j ,k)-po(i-1, j ,k))
     &    + yporel( j )*( po(i+1,j+1,k) - po(i+1,j-1,k) )
     &    - yporel( j )*( po(i-1,j+1,k) - po(i-1,j-1,k) )
     &    - yporel(j+1)*( po(i+1,j+1,k) - po(i-1,j+1,k) )
     &    + yporel(j-1)*( po(i+1,j-1,k) - po(i-1,j-1,k) )
     &    - po(i+1, j ,k)*( yporel(j+1) - yporel(j-1) )
     &    + po(i-1, j ,k)*( yporel(j+1) - yporel(j-1) ) )

*           Compute interface displacement contribution to the Jacobian term
            qotifd(i,j) = adfaco*
     &   ( (qifd(i+1, j )-qifd(i-1, j ))*(po( i ,j+1,k)-po( i ,j-1,k))
     &    +(qifd( i ,j-1)-qifd( i ,j+1))*(po(i+1, j ,k)-po(i-1, j ,k))
     &    + qifd(i+1, j )*( po(i+1,j+1,k) - po(i+1,j-1,k) )
     &    - qifd(i-1, j )*( po(i-1,j+1,k) - po(i-1,j-1,k) )
     &    - qifd( i ,j+1)*( po(i+1,j+1,k) - po(i-1,j+1,k) )
     &    + qifd( i ,j-1)*( po(i+1,j-1,k) - po(i-1,j-1,k) )
     &    + po( i ,j+1,k)*( qifd(i+1,j+1) - qifd(i-1,j+1) )
     &    - po( i ,j-1,k)*( qifd(i+1,j-1) - qifd(i-1,j-1) )
     &    - po(i+1, j ,k)*( qifd(i+1,j+1) - qifd(i+1,j-1) )
     &    + po(i-1, j ,k)*( qifd(i-1,j+1) - qifd(i-1,j-1) ) )

*           (layer-dependent) entrainment/forcing term
            if ( k.eq.1 ) then
              qotent(i,j) = fohfac(1)*( wekpo(i,j) - entoc(i,j) )
             else if ( k.eq.2 ) then
              qotent(i,j) = fohfac(2)*entoc(i,j)
             else
              qotent(i,j) = 0.0d0
            endif
*           Add bottom drag contribution
            if ( k.eq.nlo ) then
              qotent(i,j) = qotent(i,j) - bdrfac*d2pm(i,j)
            endif
            dqdt(i,j) =  qotjac(i,j) + qt2dif(i,j)
     &                 + qt4dif(i,j) + qotent(i,j)
          enddo
*         Eastern boundary (i=nxpo)
#  ifdef cyclic_ocean
*         Zonally cyclic ocean
          qt2dif(nxpo,j) = qt2dif(1,j)
          qt4dif(nxpo,j) = qt4dif(1,j)
          qotjac(nxpo,j) = qotjac(1,j)
          qotd2p(nxpo,j) = qotd2p(1,j)
          qotbet(nxpo,j) = qotbet(1,j)
          qotifd(nxpo,j) = qotifd(1,j)
          qotent(nxpo,j) = qotent(1,j)
          dqdt(nxpo,j)   = dqdt(1,j)
#  else
*         Finite box ocean
          qt2dif(nxpo,j) = 0.0d0
          qt4dif(nxpo,j) = 0.0d0
          qotent(nxpo,j) = 0.0d0
*         Compute scaled normal derivative (1/f0)*d(po)/dx
*         using the mixed boundary condition (Appendix C).
          pxof0 = pnormf*( po(nxpo,j,k) - po(nxpo-1,j,k) )

*         On W & E boundaries, Del-sqd term reduces to - (1/f0^2)*d(pxx)/dy*d(po)/dx
          qotd2p(nxpo,j) = 0.5d0*(adfaco/fnot)*
     &   ((d2pc(nxpo,j)-d2pc(nxpo-1,j))*(po(nxpo-1,j+1,k)-po(nxpo-1,j-1,k))
     &   +(d2pc(nxpo,j-1)  -d2pc(nxpo,j+1)
     &    +d2pc(nxpo-1,j-1)-d2pc(nxpo-1,j+1))
     &    *(po(nxpo, j ,k)-po(nxpo-1, j ,k))
C     &    + d2pc(nxpo, j)*( po(nxpo,j+1,k)-po(nxpo,j-1,k) )
     &    - d2pc(nxpo-1, j)*( po(nxpo-1,j+1,k) - po(nxpo-1,j-1,k) )
     & -(d2pc(nxpo-1,j+1)+d2pc(nxpo,j+1))*(po(nxpo,j+1,k)-po(nxpo-1,j+1,k))
     & +(d2pc(nxpo-1,j-1)+d2pc(nxpo,j-1))*(po(nxpo,j-1,k)-po(nxpo-1,j-1,k))
     & +(po(nxpo-1,j+1,k)+po(nxpo,j+1,k))*(d2pc(nxpo,j+1)-d2pc(nxpo-1,j+1))
     & -(po(nxpo-1,j-1,k)+po(nxpo,j+1,k))*(d2pc(nxpo,j-1)-d2pc(nxpo-1,j-1))
     &    - po(nxpo, j ,k)*(d2pc(nxpo,j+1)-d2pc(nxpo,j-1) )
     &    + po(nxpo-1, j ,k)*(d2pc(nxpo-1,j+1)-d2pc(nxpo-1,j-1) ) )

          qotifd(nxpo,j) = 0.5d0*adfaco*
     &   ((qifd(nxpo,j)-qifd(nxpo-1,j))*(po(nxpo-1,j+1,k)-po(nxpo-1,j-1,k))
     &   +(qifd(nxpo,j-1)-qifd(nxpo,j+1)+qifd(nxpo-1,j-1)-qifd(nxpo-1,j+1))
     &    *(po(nxpo, j ,k)-po(nxpo-1, j ,k))
C     &    + qo(nxpo, j ,k)*( po(nxpo,j+1,k)-po(nxpo,j-1,k) )
     &    - qifd(nxpo-1, j)*( po(nxpo-1,j+1,k) - po(nxpo-1,j-1,k) )
     &    - (qifd(nxpo-1,j+1)+qifd(nxpo,j+1))*(po(nxpo,j+1,k)-po(nxpo-1,j+1,k))
     &    + (qifd(nxpo-1,j-1)+qifd(nxpo,j-1))*(po(nxpo,j-1,k)-po(nxpo-1,j-1,k))
     &    + (po(nxpo-1,j+1,k)+po(nxpo,j+1,k))*(qifd(nxpo,j+1)-qifd(nxpo-1,j+1))
     &    - (po(nxpo-1,j-1,k)+po(nxpo,j+1,k))*(qifd(nxpo,j-1)-qifd(nxpo-1,j-1))
     &    - po(nxpo, j ,k)*( qifd(nxpo,j+1) - qifd(nxpo,j-1) )
     &    + po(nxpo-1, j ,k)*( qifd(nxpo-1,j+1) - qifd(nxpo-1,j-1) ) )
C          qotjac(nxpo,j) = - pxof0*hdxom1*( qo(nxpo,j+1,k) - qo(nxpo,j-1,k) )

          qotjac(nxpo,j) = 0.5d0*adfaco*
     &    ( (qo(nxpo,j,k)-qo(nxpo-1,j,k))*(po(nxpo-1,j+1,k)-po(nxpo-1,j-1,k))
     &     +(qo(nxpo,j-1,k)-qo(nxpo,j+1,k)+qo(nxpo-1,j-1,k)-qo(nxpo-1,j+1,k))
     &    *(po(nxpo, j ,k)-po(nxpo-1, j ,k))
C     &    +qo(nxpo, j ,k)*( po(nxpo,j+1,k)-po(nxpo,j-1,k) )
     &   -  qo(nxpo-1, j ,k)*( po(nxpo-1,j+1,k) - po(nxpo-1,j-1,k) )
     &   - (qo(nxpo-1,j+1,k)+ qo(nxpo,j+1,k))*(po(nxpo,j+1,k)-po(nxpo-1,j+1,k))
     &   + (qo(nxpo-1,j-1,k)+ qo(nxpo,j-1,k))*(po(nxpo,j-1,k)-po(nxpo-1,j-1,k))
     &   + (po(nxpo-1,j+1,k)+ po(nxpo,j+1,k))*(qo(nxpo,j+1,k)-qo(nxpo-1,j+1,k))
     &   - (po(nxpo-1,j-1,k)+ po(nxpo,j+1,k))*(qo(nxpo,j-1,k)-qo(nxpo-1,j-1,k))
     &   -  po( nxpo , j ,k)*(qo( nxpo ,j+1,k)-qo(nxpo,j-1,k) )
     &   +  po(nxpo-1, j ,k)*(qo(nxpo-1,j+1,k)-qo(nxpo-1,j-1,k) ) )

*         On W & E boundaries, "beta" term = -beta*(px/f0)
          qotbet(nxpo,j) = 0.5d0*adfaco*beta*
c     &   ((yporel(j)-yporel(j))*(po(nxpo-1,j+1,k)-po(nxpo-1,j-1,k))
     &   ( (yporel(j-1)-yporel(j+1)+yporel(j-1)-yporel(j+1))*
     &     (po(nxpo, j ,k)- po(nxpo-1, j ,k))
C     &   + yporel( j ) * ( po( nxpo ,j+1,k)-po( nxpo ,j-1,k) )
     &    - yporel( j ) * ( po(nxpo-1,j+1,k)-po(nxpo-1,j-1,k) )
     &    -(yporel(j+1) + yporel(j+1)) * ( po( nxpo ,j+1,k)-po(nxpo-1,j+1,k))
     &    +(yporel(j-1) + yporel(j-1)) * ( po( nxpo ,j-1,k)-po(nxpo-1,j-1,k))
C     &   +(po(nxpo-1,j+1,k)+po(nxpo,j+1,k))*(yporel(j+1)-yporel(j+1))
C     &   -(po(nxpo-1,j-1,k)+po(nxpo,j+1,k))*(yporel(j-1)-yporel(j-1))
     &    - po( nxpo , j ,k)*( yporel(j+1) - yporel(j-1) )
     &    + po(nxpo-1, j ,k)*( yporel(j+1) - yporel(j-1) ) )

          dqdt(nxpo,j)   = rdto*( qo(nxpo,j,k) - qom(nxpo,j,k) )
#  endif /* cyclic_ocean */
        enddo
!$OMP   END DO

*       Add zonal boundary contributions to Del-sqd(po) and interface displacement
*       components of the Jacobian. The "beta" term component vanishes identically.
!$OMP   DO SCHEDULE (STATIC)
        do i=2,nxpo-1
*         On S & N boundaries, Del-sqd term reduces to (1/f0^2)*d(pyy)/dx*d(po)/dy
*         Southern boundary
*         Compute scaled normal derivative (1/f0)*d(po)/dy
*         using the mixed boundary condition (Appendix C).
          pyof0 = pnormf*( po(i,2,k) - po(i,1,k) )
          qotd2p(i,1) = 0.5d0*(adfaco/fnot)*
     &   ( ((d2pc(i+1,1)-d2pc(i-1,1))+d2pc(i+1,2)-d2pc(i-1,2))
     &     *(po(i,2,k)-po(i,1,k))
     &    +(d2pc( i ,1)-d2pc( i ,2))*(po(i+1, 2 ,k)-po(i-1, 2 ,k))
     &    +(d2pc(i+1,1)+d2pc(i+1,2))*( po(i+1,2,k) - po(i+1,1,k) )
     &    - (d2pc(i-1,1)+d2pc(i-1,2))*( po(i-1,2,k) - po(i-1,1,k) )
     &    - d2pc( i ,2)*( po(i+1,2,k) - po(i-1,2,k) )
C     &    + d2pc( i ,1)*( po(i+1,1,k) - po(i-1,1,k) )
     &    + po( i ,2,k)*(d2pc(i+1,2) - d2pc(i-1,2) )
     &    - po( i ,1,k)*(d2pc(i+1,1) - d2pc(i-1,1) )
     &    - (po(i+1,1,k)+po(i+1,2,k))*(d2pc(i+1,2) - d2pc(i+1,1) )
     &    + (po(i-1,1,k)+po(i-1,2,k))*(d2pc(i-1,2) - d2pc(i-1,1) ) )

          qotifd(i,1) = 0.5d0*adfaco*
     &   ( ((qifd(i+1,1)-qifd(i-1,1))+qifd(i+1,2)-qifd(i-1,2))
     &     *(po(i,2,k)-po(i,1,k))
     &    +(qifd( i ,1)-qifd(i,2))*(po(i+1, 2 ,k)-po(i-1, 2 ,k))
     &    + (qifd(i+1,1)+qifd(i+1,2))*( po(i+1,2,k) - po(i+1,1,k) )
     &    - (qifd(i-1,1)+qifd(i-1,2))*( po(i-1,2,k) - po(i-1,1,k) )
     &    - qifd( i ,2)*( po(i+1,2,k) - po(i-1,2,k) )
C     &    + qifd(i,j-1)*( po(i+1,1,k) - po(i-1,1,k) )
     &    + po( i ,2,k)*( qifd(i+1,2) - qifd(i-1,2) )
     &    - po( i ,1,k)*( qifd(i+1,1) - qifd(i-1,1) )
     &    - (po(i+1,1,k)+po(i+1,2,k))*( qifd(i+1,2) - qifd(i+1,1) )
     &    + (po(i-1,1,k)+po(i-1,2,k))*( qifd(i-1,2) - qifd(i-1,1) ) )

          qotjac(i,1) = 0.5d0*adfaco*
     &   ( ((qo(i+1,1,k)- qo(i-1,1,k))+qo(i+1,2,k)-qo(i-1,2,k))
     &    * (po(i,  2,k)- po(i,  1,k))
     &    + (qo(i,  1,k)- qo( i ,2,k))*( po(i+1,2,k) - po(i-1,2,k) )
     &    + (qo(i+1,1,k)+ qo(i+1,2,k))*( po(i+1,2,k) - po(i+1,1,k) )
     &    - (qo(i-1,1,k)+ qo(i-1,2,k))*( po(i-1,2,k) - po(i-1,1,k) )
     &    -  qo( i ,2,k)*(po(i+1,2,k) -  po(i-1,2,k) )
C     &    + qo( i ,1,k)*(po(i+1,1,k) -  po(i-1,1,k) )
     &    +  po( i ,2,k)*(qo(i+1,2,k) -  qo(i-1,2,k) )
     &    -  po( i ,1,k)*(qo(i+1,1,k) -  qo(i-1,1,k) )
     &    - (po(i+1,1,k)+ po(i+1,2,k))*( qo(i+1,2,k) - qo(i+1,1,k) )
     &    + (po(i-1,1,k)+ po(i-1,2,k))*( qo(i-1,2,k) - qo(i-1,1,k) ) )

          qotbet(i,1) = 0.5d0*adfaco*beta*
c     &   ( ((qo(i+1,1,k)- qo(i-1,1,k))+qo(i+1,2,k)-qo(i-1,2,k))
c     &    * (po(i,  2,k)- po(i,  1,k))
     &    ( (yporel(1)- yporel(2))*( po(i+1,2,k) - po(i-1,2,k) )
     &    + (yporel(1)+ yporel(2))*( po(i+1,2,k) - po(i+1,1,k) )
     &    - (yporel(1)+ yporel(2))*( po(i-1,2,k) - po(i-1,1,k) )
     &    -  yporel(2)*(po(i+1,2,k) -  po(i-1,2,k) )
C     &    + yporel(1)*(po(i+1,1,k) -  po(i-1,1,k) )
C     &    +  po( i ,2,k)*(qo(i+1,2,k) -  qo(i-1,2,k) )
C     &    -  po( i ,1,k)*(qo(i+1,1,k) -  qo(i-1,1,k) )
     &    - (po(i+1,1,k)+ po(i+1,2,k))*( yporel(2) - yporel(1) )
     &    + (po(i-1,1,k)+ po(i-1,2,k))*( yporel(2) - yporel(1) ) )
*         Northern boundary
*         Compute scaled normal derivative (1/f0)*d(po)/dy
*         using the mixed boundary condition (Appendix C).
          pyof0 = pnormf*( po(i,nypo,k) - po(i,nypo-1,k) )
          qotd2p(i,nypo) = 0.5d0*(adfaco/fnot)*
     &  ( (d2pc(i+1,nypo)-d2pc(i-1,nypo)+
     &     d2pc(i+1,nypo-1)-d2pc(i-1,nypo-1))
     &    *(po(i,nypo,k)-po(i,nypo-1,k))
     &   +(d2pc(i,nypo-1)-d2pc(i,nypo))*(po(i+1,nypo-1,k)-po(i-1,nypo-1,k))
     &   +(d2pc(i+1,nypo-1)+d2pc(i+1,nypo))*(po(i+1,nypo,k)-po(i+1,nypo-1,k))
     &   -(d2pc(i-1,nypo-1)+d2pc(i-1,nypo))*(po(i-1,nypo,k)-po(i-1,nypo-1,k))
C     &    - qo( i ,nypo,k)*( po(i+1,nypo,k) - po(i-1,nypo,k) )
     &    + d2pc( i ,nypo-1)*( po(i+1,nypo-1,k) - po(i-1,nypo-1,k))
     &    + po( i ,nypo  ,k)*(d2pc(i+1,nypo  ) - d2pc(i-1,nypo  ))
     &    - po( i ,nypo-1,k)*(d2pc(i+1,nypo-1) - d2pc(i-1,nypo-1))
     &  -(po(i+1,nypo-1,k)+po(i+1,nypo,k))*(d2pc(i+1,nypo)-d2pc(i+1,nypo-1))
     &  +(po(i-1,nypo-1,k)+po(i-1,nypo,k))*(d2pc(i-1,nypo)-d2pc(i-1,nypo-1)) )

          qotifd(i,nypo) = 0.5d0*adfaco*
     &  ( (qifd(i+1,nypo)-qifd(i-1,nypo)+qifd(i+1,nypo-1)-qifd(i-1,nypo-1))
     &    *(po(i,nypo,k)-po(i,nypo-1,k))
     &   +(qifd(i,nypo-1)-qifd(i,nypo))*(po(i+1,nypo-1,k)-po(i-1,nypo-1,k))
     &    + (qifd(i+1,nypo-1)+qifd(i+1,nypo))*(po(i+1,nypo,k)-po(i+1,nypo-1,k))
     &    - (qifd(i-1,nypo-1)+qifd(i-1,nypo))*(po(i-1,nypo,k)-po(i-1,nypo-1,k))
C     &    - qifd( i ,nypo)*( po(i+1,nypo,k) - po(i-1,nypo,k) )
     &    + qifd( i ,nypo-1)*( po(i+1,nypo-1,k) - po(i-1,nypo-1,k) )
     &    + po( i ,nypo,k)*( qifd(i+1,nypo) - qifd(i-1,nypo) )
     &    - po( i ,nypo-1,k)*( qifd(i+1,nypo-1) - qifd(i-1,nypo-1) )
     &    - (po(i+1,nypo-1,k)+po(i+1,nypo,k))*(qifd(i+1,nypo)-qifd(i+1,nypo-1))
     &    + (po(i-1,nypo-1,k)+po(i-1,nypo,k))*(qifd(i-1,nypo)-qifd(i-1,nypo-1)))

          qotjac(i,nypo) = 0.5d0*adfaco*
     &  ( (qo(i+1,nypo,k) - qo(i-1, nypo ,k)+qo(i+1,nypo-1,k)-qo(i-1,nypo-1,k))
     &   *(po( i ,nypo,k) - po( i ,nypo-1,k))
     &   +(qo( i ,nypo-1,k)-qo( i ,nypo,k))*(po(i+1,nypo-1,k)-po(i-1,nypo-1,k))
     &   +(qo(i+1,nypo-1,k)+qo(i+1,nypo,k))*(po(i+1, nypo ,k)-po(i+1,nypo-1,k))
     &   -(qo(i-1,nypo-1,k)+qo(i-1,nypo,k))*(po(i-1, nypo ,k)-po(i-1,nypo-1,k))
C     &    - qo( i ,nypo,k)*( po(i+1,nypo,k) - po(i-1,nypo,k) )
     &   + qo( i ,nypo-1,k)*( po(i+1,nypo-1,k) - po(i-1,nypo-1,k) )
     &   + po( i , nypo ,k)*( qo(i+1, nypo ,k) - qo(i-1, nypo ,k) )
     &   - po( i ,nypo-1,k)*( qo(i+1,nypo-1,k) - qo(i-1,nypo-1,k) )
     &   -(po(i+1,nypo-1,k) + po(i+1,nypo,k))*(qo(i+1,nypo,k)-qo(i+1,nypo-1,k))
     &   +(po(i-1,nypo-1,k) + po(i-1,nypo,k))*(qo(i-1,nypo,k)-qo(i-1,nypo-1,k)))

        qotbet(i,nypo) = 0.5d0*adfaco*beta*
C     &  ( (qo(i+1,nypo,k) - qo(i-1, nypo ,k)+qo(i+1,nypo-1,k)-qo(i-1,nypo-1,k))
C     &    (po( i ,nypo,k) - po( i ,nypo-1,k))
     &  ( (yporel(nypo-1)-yporel(nypo))*(po(i+1,nypo-1,k)-po(i-1,nypo-1,k))
     &   +(yporel(nypo-1)+yporel(nypo))*(po(i+1, nypo ,k)-po(i+1,nypo-1,k))
     &   -(yporel(nypo-1)+yporel(nypo))*(po(i-1, nypo ,k)-po(i-1,nypo-1,k))
C     &    - qo( i ,nypo,k)*( po(i+1, nypo ,k) - po(i-1,nypo,k) )
     &   + yporel( nypo-1 )*( po(i+1,nypo-1,k) - po(i-1,nypo-1,k) )
C     &   +po( i , nypo ,k)*( yporel( nypo ) - yporel( nypo ) )
C     &   -po( i ,nypo-1,k)*( yporel(nypo-1) - yporel(nypo-1) )
     &   -(po(i+1,nypo-1,k) + po(i+1,nypo,k))*(yporel(nypo)-yporel(nypo-1))
     &   +(po(i-1,nypo-1,k) + po(i-1,nypo,k))*(yporel(nypo)-yporel(nypo-1)))
*         On S & N boundaries, "beta" term is exactly zero
        enddo
!$OMP   END DO
#  ifdef cyclic_ocean
!$OMP   SINGLE
*       On S & N boundaries, Del-sqd term reduces to (1/f0^2)*d(pyy)/dx*d(po)/dy
*       Special case i = 1
*       Southern boundary
*       Compute scaled normal derivative (1/f0)*d(po)/dy
*       using the mixed boundary condition (Appendix C).
        pyof0 = pnormf*( po(1,2,k) - po(1,1,k) )
        qotd2p(1,1) = pyof0*dtangf*( d2pc(2,1) - d2pc(nxpo-1,1) )
        qotifd(1,1) = pyof0*hdxom1*( qifd(2,1) - qifd(nxpo-1,1) )
        qotjac(1,1) = pyof0*hdxom1*( qo(2,1,k) - qo(nxpo-1,1,k) )
*       Northern boundary
*       Compute scaled normal derivative (1/f0)*d(po)/dy
*       using the mixed boundary condition (Appendix C).
        pyof0 = pnormf*( po(1,nypo,k) - po(1,nypo-1,k) )
        qotd2p(1,nypo) = pyof0*dtangf*( d2pc(2,nypo) - d2pc(nxpo-1,nypo))
        qotifd(1,nypo) = pyof0*hdxom1*( qifd(2,nypo) - qifd(nxpo-1,nypo))
        qotjac(1,nypo) = pyof0*hdxom1*( qo(2,nypo,k) - qo(nxpo-1,nypo,k))
*       Cyclicity in x at i = nxpo
        qotd2p(nxpo,   1) = qotd2p(1,   1)
        qotd2p(nxpo,nypo) = qotd2p(1,nypo)
        qotifd(nxpo,   1) = qotifd(1,   1)
        qotifd(nxpo,nypo) = qotifd(1,nypo)
        qotjac(nxpo,   1) = qotjac(1,   1)
        qotjac(nxpo,nypo) = qotjac(1,nypo)
!$OMP   END SINGLE
#  endif /* cyclic_ocean */

*       Add zonal boundary contributions (to dqdt only), by differencing
*       qo in time. This is slightly displaced in time, but only by dto/2
!$OMP   DO SCHEDULE (STATIC)
        do i=1,nxpo
          dqdt(i,  1 ) = rdto*( qo(i,  1 ,k) - qom(i,  1 ,k) )
          dqdt(i,nypo) = rdto*( qo(i,nypo,k) - qom(i,nypo,k) )
        enddo
!$OMP   END DO

*       Store the (optionally subsampled) results for the current layer
*       ---------------------------------------------------------------
!!      Subsample the dqdt array for the current layer
!!      into contiguous values in the wrk vector
!$OMP   DO SCHEDULE (STATIC)
        do j=1,jpwk
           do i=1,ipwk
              wrk(i+ipwk*(j-1)) = dqdt(1+(i-1)*nsko,1+(j-1)*nsko)
           enddo
        enddo
!$OMP   END DO
!!      Output the subsampled values in wrk to the netCDF file
!$OMP   SINGLE
        ncstat = nf_put_vara_double(qocncid, dqdt_id,
     &                              start4, count4, wrk)
        if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
!$OMP   END SINGLE

!!      Subsample the qotjac array for the current layer
!!      into contiguous values in the wrk vector
!$OMP   DO SCHEDULE (STATIC)
        do j=1,jpwk
           do i=1,ipwk
              wrk(i+ipwk*(j-1)) = qotjac(1+(i-1)*nsko,1+(j-1)*nsko)
           enddo
        enddo
!$OMP   END DO
!!      Output the subsampled values in wrk to the netCDF file
!$OMP   SINGLE
        ncstat = nf_put_vara_double(qocncid, qotjac_id,
     &                              start4, count4, wrk)
        if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
!$OMP   END SINGLE

!!      Subsample the qotd2p array for the current layer
!!      into contiguous values in the wrk vector
!$OMP   DO SCHEDULE (STATIC)
        do j=1,jpwk
           do i=1,ipwk
              wrk(i+ipwk*(j-1)) = qotd2p(1+(i-1)*nsko,1+(j-1)*nsko)
           enddo
        enddo
!$OMP   END DO
!!      Output the subsampled values in wrk to the netCDF file
!$OMP   SINGLE
        ncstat = nf_put_vara_double(qocncid, qotd2p_id,
     &                              start4, count4, wrk)
        if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
!$OMP   END SINGLE

!!      Subsample the qotbet array for the current layer
!!      into contiguous values in the wrk vector
!$OMP   DO SCHEDULE (STATIC)
        do j=1,jpwk
           do i=1,ipwk
              wrk(i+ipwk*(j-1)) = qotbet(1+(i-1)*nsko,1+(j-1)*nsko)
           enddo
        enddo
!$OMP   END DO
!!      Output the subsampled values in wrk to the netCDF file
!$OMP   SINGLE
        ncstat = nf_put_vara_double(qocncid, qotbet_id,
     &                              start4, count4, wrk)
        if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
!$OMP   END SINGLE

!!      Subsample the qotifd array for the current layer
!!      into contiguous values in the wrk vector
!$OMP   DO SCHEDULE (STATIC)
        do j=1,jpwk
           do i=1,ipwk
              wrk(i+ipwk*(j-1)) = qotifd(1+(i-1)*nsko,1+(j-1)*nsko)
           enddo
        enddo
!$OMP   END DO
!!      Output the subsampled values in wrk to the netCDF file
!$OMP   SINGLE
        ncstat = nf_put_vara_double(qocncid, qotifd_id,
     &                              start4, count4, wrk)
        if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
!$OMP   END SINGLE

!!      Subsample the qt2dif array for the current layer
!!      into contiguous values in the wrk vector
        if ( qt2dif_out ) then
!$OMP     DO SCHEDULE (STATIC)
          do j=1,jpwk
             do i=1,ipwk
                wrk(i+ipwk*(j-1)) = qt2dif(1+(i-1)*nsko,1+(j-1)*nsko)
             enddo
          enddo
!$OMP     END DO
!!        Output the subsampled values in wrk to the netCDF file
!$OMP     SINGLE
          ncstat = nf_put_vara_double(qocncid, qt2dif_id,
     &                                start4, count4, wrk)
          if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
!$OMP     END SINGLE
        endif

!!      Subsample the qt4dif array for the current layer
!!      into contiguous values in the wrk vector
!$OMP   DO SCHEDULE (STATIC)
        do j=1,jpwk
           do i=1,ipwk
              wrk(i+ipwk*(j-1)) = qt4dif(1+(i-1)*nsko,1+(j-1)*nsko)
           enddo
        enddo
!$OMP   END DO
!!      Output the subsampled values in wrk to the netCDF file
!$OMP   SINGLE
        ncstat = nf_put_vara_double(qocncid, qt4dif_id,
     &                              start4, count4, wrk)
        if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
!$OMP   END SINGLE

!!      Subsample the qotent array for the current layer
!!      into contiguous values in the wrk vector
!$OMP   DO SCHEDULE (STATIC)
        do j=1,jpwk
           do i=1,ipwk
              wrk(i+ipwk*(j-1)) = qotent(1+(i-1)*nsko,1+(j-1)*nsko)
           enddo
        enddo
!$OMP   END DO
!!      Output the subsampled values in wrk to the netCDF file
!$OMP   SINGLE
        ncstat = nf_put_vara_double(qocncid, qotent_id,
     &                              start4, count4, wrk)
        if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
!$OMP   END SINGLE NOWAIT

!$OMP   END PARALLEL

      enddo

#    endif /* use_netcdf */
#  endif /* qoc_diag */
#endif /* not atmos_only */

*     Detailed notes:
*     ===============
*     dqo/dt = (1/f0)*[ J(qo,po) + ah2oc*Del-4th(pom) + ah4oc*Del-6th(pom) ] + B_arr*ent_vec
*     (7.14), where J(qo,po) = d(qo)/dx*d(po)/dy - d(qo)/dy*d(po)/dx,
*     and qo = (1/f0)*Del-sqd(po) + beta*(y - y0) - f0*A_mat*po_vec + Dtopog  (7.15)
*     where qo and po are the values at the current timestep, and pom is the
*     value at the previous timestep (diffusive terms need to be time-lagged
*     for numerical stability). The vorticity is only time-stepped at internal points.
*     For boundary points, qo is computed using (7.15) once po has been updated by
*     inverting the new qo field and applying suitable constraints. Thus dqo/dt
*     on boundaries has to be inferred by finite-differencing the time levels.

*     Along all boundaries, po is constant, so its derivatives along the boundary
*     all vanish. At the corners of the domain (box_ocean case only) both px and py
*     are zero, so the Jacobian term and its individual components all vanish.
*     S & N boundaries: d(po)/dx = 0, so the Jacobian reduces to d(qo)/dx*d(po)/dy.
*     W & E boundaries: In the cyclic_ocean case there are no physical boundaries,
*     merely a cyclicity condition in x across the ends of the domain. In the
*     box_ocean case, d(po)/dy = 0, so the Jacobian reduces to - d(qo)/dy*d(po)/dx.
*     Note that values on the boundary are evaluated by simple finite differencing,
*     so may not be quite consistent with a reduced form of the Arakawa formulation
*     that should ideally be used. This can be improved if a 6 point formulation
*     of the Arakawa scheme appropriate adjacent to boundaries can be found.
*     All boundary cases involve a factor of (1/f0)*normal derivative of p,
*     so this evaluated, then used for each of the three components below.

*     "Del-sqd(po)" term in Jacobian contribution to dq/dt:
*     -----------------------------------------------------
*     Del-sqd(po) term = (1/f0)*J( (1/f0)*Del-sqd(po), po) = (1/f0^2)*J( Del-sqd(po), po)
*     At internal points, use the usual Arakawa Jacobian formulation, with
*     qo(i,j,k) replaced by d2pc(i,j). Computation loop is within an outer k loop.
*     S & N boundaries: d2(po)/dx2 vanishes identically on these boundaries.
*     So the contribution reduces to (1/f0^2)*d(pyy)/dx*d(po)/dy.
*     W & E boundaries: d2(po)/dy2 vanishes identically on these boundaries.
*     So the contribution reduces to - (1/f0^2)*d(pxx)/dy*d(po)/dx.
*     The normal derivatives can be got from the mixed boundary condition.
*     The third derivatives d(pyy)/dx and d(pxx)/dy can be got by centred differences
*     along the boundaries of the already computed pyy, pxx terms in the d2pc array.

*     "beta" term in Jacobian contribution to dq/dt:
*     ----------------------------------------------
*     beta term = (1/f0)*J( beta*(y - y0), po )
*     The beta term is independent of x, so for the beta term d(qo)/dx is identically
*     zero everywhere. At internal points, evaluate this term using the full Arakawa
*     Jacobian, but simplified by eliminating those terms where qo is differenced in x.
*     S & N boundaries: as previously explained d(qo)/dx is identically zero for the
*     beta term, so this component of the Jacobian vanishes exactly on these boundaries.
*     W & E boundaries: The normal derivative d(po)/dx is given by a mixed boundary
*     condition, implemented as (C.5) and (C.7) (but with x-derivatives), and
*     d(qo)/dy = beta. So the beta term contribution reduces to -(1/f0)*beta*d(po)/dx .

*     "interface displacement" term in Jacobian contribution to dq/dt:
*     ----------------------------------------------------------------
*     interface term = (1/f0)*J( - f0*A_mat*po_vec + Dtopog, po )
*     At internal points, use the usual Arakawa Jacobian formulation, with
*     qo(i,j,k) replaced by qifd(i,j). Computation loop is within an outer k loop.
*     S & N boundaries: the interface displacement component of d(qo)/dx
*     is got by finite differencing qifd along the boundary. Cyclicity in x
*     enables d(qo)/dx to be computed at all points in the cyclic_ocean case.
*     As argued above, the internal points suffice for the box_ocean case,
*     because at the corner points the Jacobian can make no contribution.

      END SUBROUTINE qocdiag_out

c***********************************************************************
*
      END MODULE qocdiag
*
c***********************************************************************
