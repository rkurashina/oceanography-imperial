c
c THIS IS A MULTI-LAYER QG MODEL WITH CABARET
c

C---------------might need to add averaged file------------------

      use MOD_qg3_NETCDF

      implicit none
      integer ips,nn,max_time,ii,jj,ii1,jj1,ii2,jj2,nw,islip,istart
     & ,idata,k_rec,K,MAX_SPIN,KK,K_AV,M,IFMAX,IFMIN,JBL
      real*8 basinscale,D,tau,visc,visc_bot,alpha,weight1,weight2,UMAX
     &,DL1,DL2,TWEIGHT,TWEIGHT_A,TWEIGHT_A1,TWEIGHTU,TWEIGHTU1
     & ,DT_O,DT_XO,DT05,TIME_AV

      parameter(ips=1,nn=3,jj=513,ii=jj,idata=1
     & ,max_time=101,max_spin=0)
      parameter(basinscale=3840.D5,D=4.D5
     & ,tau=.4D0
     & ,visc=20.D4
     & ,visc_bot=4.D-8
     & ,istart=0
     & ,islip=1,alpha=120.D5)

      character*(*),parameter :: file_name = 'QG_new.nc'

      parameter(jj1=jj-1,ii1=ii-1,ii2=ii-2,jj2=jj-2,nw=5*ii2/2+ii
     & ,weight1=0.05D0,weight2=0.90D0)

      real*8 psi(nn,ii,jj),rel(nn,ii,jj),phi(nn,ii,jj)

C CABARET VARIABLES /////////////////////////////
C CELL CENTRES
     & ,ZETA_NEW(NN,II,JJ),ZETA_OLD(NN,II,JJ)
C CELL FACES
     & ,ZETA_FI(NN,II,JJ),ZETA_FJ(NN,II,JJ)
     & ,ZETA_FIN(NN,II,JJ),ZETA_FJN(NN,II,JJ)
     & ,UI(NN,II,JJ),UJ(NN,II,JJ)
     & ,UPI(NN,II,JJ),UPJ(NN,II,JJ)

     &,UIM(NN,II,JJ),UJM(NN,II,JJ),UIMO(NN,II,JJ),UJMO(NN,II,JJ)

     &,FLUX_I(NN,II+1,JJ+1),FLUX_J(NN,II+1,JJ+1)

C CELL EDGES
     &,PSI_EDGE(NN,II,JJ)
C CELL RESIDUALS
     &, RES(NN,II,JJ),RES_VISC(NN,II,JJ),Z_TOT(3)
     &, RES_BETA(NN,II,JJ),RES_BETA_O(NN,II,JJ)

C END OF CABARET VARIABLES ///////////////////////

     & ,eig(nn,ii,jj),wsave(nn,nw),corr(nn,ii,jj)
     & ,eig_aux(ii,jj),wsave_aux(nw),corr_aux(ii,jj)
     & ,eig1(ii,jj),eig2(ii,jj),eig3(ii,jj)
     & ,wsave1(nw),wsave2(nw),wsave3(nw)

     & ,z1_aux(ii,jj),z2_aux(ii,jj),z3_aux(ii,jj)

     & ,h(nn),S(nn,2),SS(nn),theta(nn,nn),omega(nn,nn)
     & ,phi_C(nn)

     & ,force_wind(ii,jj),beta_y(ii),viscosity(ii,jj)

     & ,psi_av(nn,ii,jj)

     & ,phi_aux(nn,ii,jj),psi_aux(nn,ii,jj),z_aux_aux(ii,jj)

      real*8 z0,delta,D_aux
     & ,T,pi,dt,uscale,scale,tscale,beta,cff
     & ,beta_nondim,tau_nondim,visc_nondim,visc_bot_nondim,alpha_nondim
     & ,shift,tilt,asym,A,B,x,y,y_mid
     & ,time_day,time0_day,year
     & ,TIME_OUT,TIME_SAVE,TIME_O,TIME_S
     & ,ekin,epot
     & ,u

      real aux(jj,ii)
     & ,psi1_save(jj,ii),psi2_save(jj,ii),psi3_save(jj,ii)
     & ,unused,FI1,FI2,U1,U2,PSI_N,PSI_O,DTMAX,CFL
     & ,V0,D_F,FMIN,FMAX,DT_DAY, PSI_M,PSI_P,FI_M,FI_P

      integer ngpswk,nn_aux,ii_aux,jj_aux,n,i,j,l
     & ,k_out,k_save,k_o,k_s,k_day,max_time1,max_spin1,l_tot
     & ,i1,j1,n1

      data dt/1080./

      data pi/3.14159265358979323846D0/
      data beta/2.D-13/

      common /ONE/ psi
      common /TWO/ zeta_NEW
      common /THREE/ rel

      DT_O=0.D0
c
c      if(ips.eq.0)then
c        call opngks
c      else
c        call gopks(6,unused)
c        call gopwk(1,2,ngpswk('PS','PORTRAIT','COLOR'))
c        call gacwk(1)
c      endif
c
c--- stratification
c
      open(11,file='stratification.d',form='unformatted')
      read(11) z0,delta,D_aux,nn_aux,jj_aux
      read(11) h,S,SS
      read(11) theta,omega
      close(11)

      if(nn_aux.ne.nn) stop 'Check nn!!!'
      if(jj_aux.ne.jj) stop 'Check jj!!!'
      write(*,*)'z0=',z0,'; delta=',delta,'; D=',D

      do n=1,nn
         write(*,*)n,'; h=',h(n),'; S_n1=',S(n,1),'; S_n2=',S(n,2)
     &              ,'; SS=',SS(n)
      enddo


c ----- read in jet data
c     open(95, file='jet.dat', FORM='unformatted')
c     read(95) y_jet
c close(95)

c      do i = 1, ii
c          write(*, *) y_jet(i)
c      end do

c ----- create new files for storing psi and time averaged psi

      call create_netcdf_file(FILE_NAME,ii,jj,basinscale)
c      call create_ave_file(ave_file,ii,jj,basinscale)

C/// NUMERICAL PARAMETERS FOR OUTPUT
      CFL=0.5

      TIME_DAY=0
      TIME_AV=0
      L=0
      KK=0
      K_AV=0
c
c--- PB counters
c
      TIME_OUT=1.D0     ! accumulate data every TIME_OUT days
      TIME_SAVE=10.D0
C NOTE THAT SINCE THE COUNTERS ARE NORMILISED BY THE ORIGINAL DT,
C WHICH MAY BE VERY DIFFERENT TO THE ACTUAL TIME STEP IN THE CODE,
C THE SAVING INTERVAL IN DAYS MAY NOT CORRESPOND TO K_OUT AND K_SAVE
      TIME_O=0.D0
      TIME_S=0.D0
c
c--- nondimensionalization
c
      uscale=1.
      scale=basinscale/dfloat(ii-1)
      tscale=scale/uscale

      dt=dt/tscale

      beta_nondim=beta*scale*scale/uscale
      tau_nondim=tau*scale/(h(1)*uscale*uscale*dfloat(ii1)) ! DIVIDED BY rho=1
      visc_nondim=visc/(scale*uscale)
      visc_bot_nondim=visc_bot*scale/uscale
      alpha_nondim=scale/alpha
c
c--- counters
c
      l_tot=int(86400./(dt*tscale)+0.001)
      max_time1=max_time*l_tot/2
      max_spin1=max_spin*l_tot/2
c
c--- beta*y
c
      do i=1,ii
         beta_y(i)=beta_nondim*dfloat((ii+1)/2-i)
      enddo
c
c--- windcurl forcing
c
      shift=.0
      tilt =.2
      asym=.9

      A=-1.*asym
      B= 1./asym

      do j=1,jj
         do i=1,ii
            x= dfloat(j-1)/dfloat(jj-1)-.5
            y=-dfloat(i-1)/dfloat(ii-1)+.5

            y_mid=shift+tilt*x

            if(y.lt.y_mid)then
              force_wind(i,j)=A*sin(pi*(.5+y)/(.5+y_mid))
            else
              force_wind(i,j)=B*sin(pi*(y-y_mid)/(.5-y_mid))
            endif
            force_wind(i,j)=tau_nondim*2.*pi*force_wind(i,j)
         enddo
      enddo

      do j=1,jj
         do i=1,ii
            viscosity(i,j)=visc_nondim
         enddo
      enddo
c
c--- initialization for elliptic solvers
c
      call poinit(SS(1),ii2,jj2,ii,eig1,wsave1)
      call poinit(SS(2),ii2,jj2,ii,eig2,wsave2)
      call poinit(SS(3),ii2,jj2,ii,eig3,wsave3)

      do n=1,nn
         call poinit(SS(n),ii2,jj2,ii,eig_aux,wsave_aux)
         do j=1,jj
            do i=1,ii
               eig(n,i,j)=eig_aux(i,j)
            enddo
         enddo
         do n1=1,nw
            wsave(n,n1)=wsave_aux(n1)
         enddo
      enddo
c
c--- correction functions
c
      call correction_functions(nn,ii,jj,nw,eig,wsave,SS
     &                            ,eig_aux,wsave_aux,corr_aux,corr)
c
c--- INITIAL CONDITIONS
c
      do j=1,jj
         do i=1,ii
            do n=1,nn
               psi(n,i,j)=0.
               rel(n,i,j)=0.
            enddo
         enddo
      enddo

      DO J=1,JJ
         DO I=1,II
            DO K=1,NN
               psi_av(k,i,j)=0.

               ZETA_OLD(K,I,J)=0
               ZETA_NEW(K,I,J)=0
               RES_VISC(K,I,J)=0
               RES_BETA(K,I,J)=0
               RES_BETA_O(K,I,J)=0

               UIM(K,I,J)=0
               UJM(K,I,J)=0
               UIMO(K,I,J)=0
               UJMO(K,I,J)=0
               UPI(K,I,J)=0
               UPJ(K,I,J)=0
            END DO
         END DO
      END DO

      DO J=1,JJ
         DO I=2,II
            DO K=1,NN
               ZETA_FI(K,I,J)=0
            END DO
         END DO
      END DO

      DO J=2,JJ
         DO I=1,II
            DO K=1,NN
               ZETA_FJ(K,I,J)=0
            END DO
         END DO
      END DO

      DO J=1,JJ
         DO I=1,II
            DO K=1,NN
               ZETA_FIN(K,I,J)=0
               ZETA_FJN(K,I,J)=0
            END DO
         END DO
      END DO

      DO J=1,JJ+1
         DO I=1,II+1
            DO K=1,NN
               FLUX_I(K,I,J)=0
               FLUX_J(K,I,J)=0
            END DO
         END DO
      END DO

C DEFAULT TIME STEP ////////////

      T=DT
      DL1=DT
      DL2=DT
      TWEIGHT=1.

C ///////// TO START OR TO RESTART IS THE QUESTION

      if(istart.eq.2)then
        write(*,*)'Starting from zeta-restart...'

        if(jj.eq.129) open(11,file='qg_129_start',form='formatted')
        if(jj.eq.257) open(11,file='qg_257_start',form='formatted')
        if(jj.eq.513) open(11,file='qg_513_start',form='formatted')
        if(jj.eq.1025) open(11,file='qg_1025_start',form='formatted')
        if(jj.eq.2049) open(11,file='qg_2049_start',form='formatted')

19    FORMAT(20(3X,F30.15))

        DO J=1,JJ
           DO I=1,II
              DO K=1,NN
                 READ(11,19) ZETA_OLD(K,I,J),RES(K,I,J)
     &                                      ,RES_BETA_O(K,I,J)
              END DO
           END DO
        END DO
        CLOSE(11)

C /// RECONSTRUCT REL AND PSI FROM ZETA /////////////////////////////////////

         do j=2,jj1
            do i=2,ii1
               z1_aux(i,j)= theta(1,1)*ZETA_OLD(1,i,j)
     &                     +theta(1,2)*ZETA_OLD(2,i,j)
     &                     +theta(1,3)*ZETA_OLD(3,i,j)

               z2_aux(i,j)= theta(2,1)*ZETA_OLD(1,i,j)
     &                     +theta(2,2)*ZETA_OLD(2,i,j)
     &                     +theta(2,3)*ZETA_OLD(3,i,j)

               z3_aux(i,j)= theta(3,1)*ZETA_OLD(1,i,j)
     &                     +theta(3,2)*ZETA_OLD(2,i,j)
     &                     +theta(3,3)*ZETA_OLD(3,i,j)
            enddo
         enddo

         call helm(z1_aux(2,2),ii,ii2,jj2,eig1,wsave1)
         call helm(z2_aux(2,2),ii,ii2,jj2,eig2,wsave2)
         call helm(z3_aux(2,2),ii,ii2,jj2,eig3,wsave3)

         do j=1,jj
            do i=1,ii
               phi(1,i,j)=z1_aux(i,j)
               phi(2,i,j)=z2_aux(i,j)
               phi(3,i,j)=z3_aux(i,j)
            enddo
         enddo
c
c--- correction
c
         call mean_correct(nn,ii,jj,phi,phi_C,corr)

         do j=1,jj
            do i=1,ii
               psi(1,i,j)= omega(1,1)*phi(1,i,j)
     &                    +omega(1,2)*phi(2,i,j)
     &                    +omega(1,3)*phi(3,i,j)

               psi(2,i,j)= omega(2,1)*phi(1,i,j)
     &                    +omega(2,2)*phi(2,i,j)
     &                    +omega(2,3)*phi(3,i,j)

               psi(3,i,j)= omega(3,1)*phi(1,i,j)
     &                    +omega(3,2)*phi(2,i,j)
     &                    +omega(3,3)*phi(3,i,j)
            enddo
         enddo
C FIRST DEFINE EDGE VARIABLES FOR PSI

C EDGES
         DO J=2,JJ
            DO I=2,II
               DO K=1,NN
                  PSI_EDGE(K,I,J)=.25D0*(
     &               PSI(K,I,J)+PSI(K,I,J-1)+PSI(K,I-1,J)+PSI(K,I-1,J-1)
     &                                  )
               END DO
            END DO
         END DO

         UMAX=0
         DO J=2,JJ
            DO I=2,II1
               DO K=1,NN
                  UJ(K,I,J)=(PSI_EDGE(K,I+1,J)-PSI_EDGE(K,I,J))
                  IF(ABS(UJ(K,I,J)).GT.UMAX) UMAX=ABS(UJ(K,I,J))
                  UJMO(K,I,J)=UJ(K,I,J)
               END DO
            END DO
         END DO

         DO J=2,JJ1
            DO I=2,II
               DO K=1,NN
                  UI(K,I,J)=-(PSI_EDGE(K,I,J+1)-PSI_EDGE(K,I,J))
                  IF(ABS(UI(K,I,J)).GT.UMAX) UMAX=ABS(UI(K,I,J))
                  UIMO(K,I,J)=UI(K,I,J)
               END DO
            END DO
         END DO
c
c###############################################
c
      elseif(istart.eq.1)then
        write(*,*)'Starting from psi-restart...'

        if(jj.eq.129) open(11,file='qg_129_s',form='unformatted')
        if(jj.eq.257) open(11,file='qg_257_s',form='unformatted')
        if(jj.eq.513) open(11,file='qg_513_s',form='unformatted')
        if(jj.eq.1025) open(11,file='qg_1025_s',form='unformatted')
        if(jj.eq.2049) open(11,file='qg_2049_s',form='unformatted')
        read(11) ii_aux,jj_aux,nn_aux
        read(11) psi
        close(11)

        if(ii_aux.ne.ii) stop 'Check ii!!!'
        if(jj_aux.ne.jj) stop 'Check jj!!!'
        if(nn_aux.ne.nn) stop 'Check nn!!!'

        call rel_from_psi(islip,alpha_nondim,nn,ii,jj,psi,rel)

        do j=1,jj
           do i=1,ii
              do n=1,nn
                 if(n.eq.1)then
                    zeta_OLD(n,i,j)=rel(n,i,j)
     &                              -S(n,2)*(psi(n,i,j)-psi(n+1,i,j))
                 elseif(n.eq.nn)then
                    zeta_OLD(n,i,j)=rel(n,i,j)
     &                              -S(n,1)*(psi(n,i,j)-psi(n-1,i,j))
                 else
                    zeta_OLD(n,i,j)=rel(n,i,j)
     &                              -S(n,1)*(psi(n,i,j)-psi(n-1,i,j))
     &                              -S(n,2)*(psi(n,i,j)-psi(n+1,i,j))
                 endif
                 zeta_NEW(n,i,j)=zeta_OLD(n,i,j)
              enddo
           enddo
        enddo

C FLUX-TYPE VARIABLES

        DO J=1,JJ
           DO I=1,II1
              DO K=1,NN
                 ZETA_FI(K,I+1,J)=
     &                       .5D0*(ZETA_NEW(K,I,J)+ZETA_NEW(K,I+1,J))
              END DO
           END DO
        END DO

        DO J=1,JJ1
           DO I=1,II
              DO K=1,NN
                 ZETA_FJ(K,I,J+1)=
     &                       .5D0*(ZETA_NEW(K,I,J)+ZETA_NEW(K,I,J+1))
              END DO
           END DO
        END DO

C FIRST DEFINE EDGE VARIABLES FOR PSI

C EDGES
        DO J=2,JJ
           DO I=2,II
              DO K=1,NN
                 PSI_EDGE(K,I,J)=.25D0*(
     &              PSI(K,I,J)+PSI(K,I,J-1)+PSI(K,I-1,J)+PSI(K,I-1,J-1)
     &                                 )
              END DO
           END DO
        END DO

        UMAX=0
        DO J=2,JJ
           DO I=2,II1
              DO K=1,NN
                 UJ(K,I,J)= (PSI_EDGE(K,I+1,J)-PSI_EDGE(K,I,J))
                 IF(ABS(UJ(K,I,J)).GT.UMAX) UMAX=ABS(UJ(K,I,J))
                 UJMO(K,I,J)=UJ(K,I,J)
              END DO
           END DO
        END DO

        DO J=2,JJ-1
           DO I=2,II
              DO K=1,NN
                 UI(K,I,J)=-(PSI_EDGE(K,I,J+1)-PSI_EDGE(K,I,J))
                 IF(ABS(UI(K,I,J)).GT.UMAX) UMAX=ABS(UI(K,I,J))
                 UIMO(K,I,J)=UI(K,I,J)
              END DO
           END DO
        END DO

        DO J=2,JJ1
           DO I=2,II1
              DO K=1,NN
                 RES(K,I,J)=-(FLUX_I(K,I+1,J)-FLUX_I(K,I,J)
     &                       +FLUX_J(K,I,J+1)-FLUX_J(K,I,J))
              END DO
           END DO
        END DO

        do j=2,jj1
           do i=2,ii1
              DO K=1,NN
                 RES_BETA_O(K,I,J)=DT*(
     &               .5D0*(UI(K,I+1,J)+UI(K,I,J))*BETA_NONDIM
     &                                )
              enddo
           enddo
        END DO


c
c###############################################
c
      elseif(istart.eq.0)then
        write(*,*)'Starting from the rest...'
      endif
c
c--- initialization of the data files to be recorded
c
      if(idata.eq.1)then
c        open(71,file='psi1.dat'
c     &       ,access='direct',form='unformatted',recl=4*1025*1025)
c     &       ,access='direct',form='unformatted',recl=1025*1025)    ! with ifort/gfort
c        open(72,file='psi2.dat'
c     &       ,access='direct',form='unformatted',recl=4*1025*1025)
c     &       ,access='direct',form='unformatted',recl=1025*1025)    ! with ifort/gfort
c        open(73,file='psi3.dat'
c     &       ,access='direct',form='unformatted',recl=4*1025*1025)
c     &       ,access='direct',form='unformatted',recl=1025*1025)    ! with ifort/gfort
        k_rec=0
      endif

c
c################################# MAIN CYCLE #################################
c
1001   CONTINUE

C DEFINE THE TIME STEP IN DAYS
c      write(*,*)'dt=',dt
      DT_DAY=dt*tscale/86400.
      time_day=time_day+DT_DAY
      k_day=int(time_day+0.001)
      TIME_O=TIME_O+DT_DAY
      TIME_S=TIME_S+DT_DAY
      l=l+1
      KK=KK+1


C COUNTING THE NUMBER OF DAYS
        if(TIME_O.ge.TIME_OUT-.5*DT_DAY)then
           TIME_O=0.
        endif

        if(TIME_S.ge.TIME_SAVE-.5*DT_DAY)then
           TIME_S=0.
        endif

C MAX DT ALLOWABLE
        DTMAX=1.E+6
        DO J=2,JJ1
           DO I=2,II1
              DO K=1,NN
C I
                 cff=CFL/ABS(UI(K,I,J))
                 IF(ABS(UI(K,I,J)).GT.0)THEN
                   IF(DTMAX.GT.cff) DTMAX=cff
                 ELSE
                   DTMAX=DT
                 END IF
C J
                 cff=CFL/ABS(UJ(K,I,J))
                 IF(ABS(UJ(K,I,J)).GT.0)THEN
                   IF(DTMAX.GT.cff) DTMAX=cff
                 ELSE
                   DTMAX=DT
                 END IF
              END DO
           END DO
        END DO

C EXTRAPOLATION FROM THE PREVIOUS TIME LEVELS
        DT_XO=DT_O
        DT_O=DT
        DT=MIN(DTMAX,100*T)
        DL1=DT_XO+DT_O
        DL2=DT_O+DT
        TWEIGHT=DL2/DL1
        TWEIGHTU=DT/DL2
        TWEIGHTU1=1.+TWEIGHTU
        DT05=0.5*DT
        TWEIGHT_A=DT05/DT_O
        TWEIGHT_A1=1.+TWEIGHT_A

C //////CABARET PREDICTOR STEP
        DO J=2,JJ1
           DO I=2,II1
              ZETA_NEW(1,I,J)=ZETA_OLD(1,I,J)+DT05*RES(1,I,J)
              ZETA_NEW(2,I,J)=ZETA_OLD(2,I,J)+DT05*RES(2,I,J)
              ZETA_NEW(3,I,J)=ZETA_OLD(3,I,J)+DT05*RES(3,I,J)
           END DO
        END DO

C/// TIME INEGRATION OF SOURCE TERMS
        cff=DT05*BETA_NONDIM
        do j=2,jj1
           do i=2,ii1
              RES_BETA(1,I,J)=cff*(UI(1,I+1,J)+UI(1,I,J))
              RES_BETA(2,I,J)=cff*(UI(2,I+1,J)+UI(2,I,J))
              RES_BETA(3,I,J)=cff*(UI(3,I+1,J)+UI(3,I,J))
           enddo
        enddo

C SECOND-ORDER IN TIME INTEGRATION, WHICH SUPRESSES SPURIOUS BETA_Y-WAVES;
C NB: INTERSTINGLY, THE SIMPLE OPERATOR ZETA(N+1)-ZETA(N+1/2) = L(ZETA(N+1/2)), AS SUGGESTED BY VM,
C ALSO LEADS TO A STABLE SOLUTION BUT IT IS CONTAMINATED BY SPURIOUS LARGE-SCALE BETA_Y-WAVES

        DO J=2,JJ1
           DO I=2,II1
              ZETA_NEW(1,I,J)=ZETA_NEW(1,I,J)
     &          +TWEIGHT_A1*RES_BETA(1,I,J)-TWEIGHT_A*RES_BETA_O(1,I,J)
              ZETA_NEW(2,I,J)=ZETA_NEW(2,I,J)
     &          +TWEIGHT_A1*RES_BETA(2,I,J)-TWEIGHT_A*RES_BETA_O(2,I,J)
              ZETA_NEW(3,I,J)=ZETA_NEW(3,I,J)
     &          +TWEIGHT_A1*RES_BETA(3,I,J)-TWEIGHT_A*RES_BETA_O(3,I,J)
           END DO
        END DO

c
c--- solve helmholtz problems
c
         do j=2,jj1
            do i=2,ii1
               z1_aux(i,j)= theta(1,1)*ZETA_NEW(1,i,j)
     &                     +theta(1,2)*ZETA_NEW(2,i,j)
     &                     +theta(1,3)*ZETA_NEW(3,i,j)

               z2_aux(i,j)= theta(2,1)*ZETA_NEW(1,i,j)
     &                     +theta(2,2)*ZETA_NEW(2,i,j)
     &                     +theta(2,3)*ZETA_NEW(3,i,j)

               z3_aux(i,j)= theta(3,1)*ZETA_NEW(1,i,j)
     &                     +theta(3,2)*ZETA_NEW(2,i,j)
     &                     +theta(3,3)*ZETA_NEW(3,i,j)
            enddo
         enddo

         call helm(z1_aux(2,2),ii,ii2,jj2,eig1,wsave1)
         call helm(z2_aux(2,2),ii,ii2,jj2,eig2,wsave2)
         call helm(z3_aux(2,2),ii,ii2,jj2,eig3,wsave3)

         do j=1,jj
            do i=1,ii
               phi(1,i,j)=z1_aux(i,j)
               phi(2,i,j)=z2_aux(i,j)
               phi(3,i,j)=z3_aux(i,j)
            enddo
         enddo
c
c--- correction
c
         call mean_correct(nn,ii,jj,phi,phi_C,corr)

         do j=1,jj
            do i=1,ii
               psi(1,i,j)= omega(1,1)*phi(1,i,j)
     &                    +omega(1,2)*phi(2,i,j)
     &                    +omega(1,3)*phi(3,i,j)

               psi(2,i,j)= omega(2,1)*phi(1,i,j)
     &                    +omega(2,2)*phi(2,i,j)
     &                    +omega(2,3)*phi(3,i,j)

               psi(3,i,j)= omega(3,1)*phi(1,i,j)
     &                    +omega(3,2)*phi(2,i,j)
     &                    +omega(3,3)*phi(3,i,j)
            enddo
         enddo
c
c--- relative vorticity
c
         do j=2,jj1
            do i=2,ii1
               rel(1,i,j)=ZETA_NEW(1,i,j)
     &                            +S(1,2)*(
     &                       psi(1,i,j)-psi(2,i,j)
     &                                    )
               rel(nn,i,j)=ZETA_NEW(nn,i,j)
     &                            +S(nn,1)*(
     &                       psi(nn,i,j)-psi(nn-1,i,j)
     &                                    )
               rel(2,i,j)=ZETA_NEW(2,i,j)
     &                            +S(2,1)*(
     &                       psi(2,i,j)-psi(1,i,j)
     &                                    )
     &                            +S(2,2)*(
     &                       psi(2,i,j)-psi(3,i,j)
     &                                    )
            enddo
         enddo
c
c--- boundary condition
c
         call boundary_condition(islip,alpha_nondim,nn,ii,jj,rel,psi)

         DO J=1,JJ,JJ1
            DO I=1,II
               ZETA_NEW(1,i,j)=rel(1,i,j)
     &                            -S(1,2)*(
     &                       psi(1,i,j)-psi(2,i,j)
     &                                    )
               ZETA_NEW(nn,i,j)=rel(nn,i,j)
     &                            -S(nn,1)*(
     &                       psi(nn,i,j)-psi(nn-1,i,j)
     &                                    )
               ZETA_NEW(2,i,j)=rel(2,i,j)
     &                            -S(2,1)*(
     &                       psi(2,i,j)-psi(1,i,j)
     &                                    )
     &                            -S(2,2)*(
     &                       psi(2,i,j)-psi(3,i,j)
     &                                    )
            END DO
         END DO

         DO I=1,II,II1
            DO J=2,JJ1
               ZETA_NEW(1,i,j)=rel(1,i,j)
     &                            -S(1,2)*(
     &                       psi(1,i,j)-psi(2,i,j)
     &                                    )
               ZETA_NEW(nn,i,j)=rel(nn,i,j)
     &                            -S(nn,1)*(
     &                       psi(nn,i,j)-psi(nn-1,i,j)
     &                                    )
               ZETA_NEW(2,i,j)=rel(2,i,j)
     &                            -S(2,1)*(
     &                       psi(2,i,j)-psi(1,i,j)
     &                                    )
     &                            -S(2,2)*(
     &                       psi(2,i,j)-psi(3,i,j)
     &                                    )
            END DO
         END DO

C BULK VISCOSITY

         do j=2,jj1
            do i=2,ii1
               RES_VISC(1,I,J)=DT*(force_wind(i,j)
     & +viscosity(i,j)*(-4.*REL(1,i,j)
     &                  +REL(1,i-1,j)+REL(1,i+1,j)
     &                  +REL(1,i,j-1)+REL(1,i,j+1)
     &                 )
     &                            )

               RES_VISC(NN,I,J)=DT*(
     & +viscosity(i,j)*(-4.*REL(NN,i,j)
     &                  +REL(NN,i-1,j)+REL(NN,i+1,j)
     &                  +REL(NN,i,j-1)+REL(NN,i,j+1)
     &                 )
     &        -visc_bot_nondim*rel(nn,i,j)
     &                             )

               RES_VISC(2,I,J)=DT*(
     & +viscosity(i,j)*(-4.*REL(2,i,j)
     &                  +REL(2,i-1,j)+REL(2,i+1,j)
     &                  +REL(2,i,j-1)+REL(2,i,j+1)
     &                 )
     &                            )
            enddo
         enddo

         DO J=2,JJ1
            DO I=2,II1
               ZETA_NEW(1,I,J)=ZETA_NEW(1,I,J)+RES_VISC(1,I,J)
               ZETA_NEW(2,I,J)=ZETA_NEW(2,I,J)+RES_VISC(2,I,J)
               ZETA_NEW(3,I,J)=ZETA_NEW(3,I,J)+RES_VISC(3,I,J)
            END DO
         END DO

C END OF VISCOSITY


C ///UPDATING CELL-FACE VARIABLES  /////

C FIRST DEFINE EDGE VARIABLES FOR PSI

C EDGES
         DO J=2,JJ
            DO I=2,II1
               UJM(1,I,J)=.25D0*(
     &        PSI(1,I+1,J)+PSI(1,I+1,J-1)-PSI(1,I-1,J)-PSI(1,I-1,J-1)
     &                          )
               UJM(2,I,J)=.25D0*(
     &        PSI(2,I+1,J)+PSI(2,I+1,J-1)-PSI(2,I-1,J)-PSI(2,I-1,J-1)
     &                          )
               UJM(3,I,J)=.25D0*(
     &        PSI(3,I+1,J)+PSI(3,I+1,J-1)-PSI(3,I-1,J)-PSI(3,I-1,J-1)
     &                          )
            END DO
         END DO

         DO J=2,JJ1
            DO I=2,II
               UIM(1,I,J)=-.25D0*(
     &        PSI(1,I,J+1)-PSI(1,I,J-1)+PSI(1,I-1,J+1)-PSI(1,I-1,J-1)
     &                           )
               UIM(2,I,J)=-.25D0*(
     &        PSI(2,I,J+1)-PSI(2,I,J-1)+PSI(2,I-1,J+1)-PSI(2,I-1,J-1)
     &                           )
               UIM(3,I,J)=-.25D0*(
     &        PSI(3,I,J+1)-PSI(3,I,J-1)+PSI(3,I-1,J+1)-PSI(3,I-1,J-1)
     &                           )
            END DO
         END DO

         UMAX=0
         DO J=2,JJ
            DO I=2,II1
               UJ(1,I,J)=TWEIGHTU1*UJM(1,I,J)-TWEIGHTU*UJMO(1,I,J)
               IF(ABS(UJM(1,I,J)).GT.UMAX) UMAX=ABS(UJM(1,I,J))
               UJ(2,I,J)=TWEIGHTU1*UJM(2,I,J)-TWEIGHTU*UJMO(2,I,J)
               IF(ABS(UJM(2,I,J)).GT.UMAX) UMAX=ABS(UJM(2,I,J))
               UJ(3,I,J)=TWEIGHTU1*UJM(3,I,J)-TWEIGHTU*UJMO(3,I,J)
               IF(ABS(UJM(3,I,J)).GT.UMAX) UMAX=ABS(UJM(3,I,J))
            END DO
         END DO

         DO J=2,JJ1
            DO I=2,II
               UI(1,I,J)=TWEIGHTU1*UIM(1,I,J)-TWEIGHTU*UIMO(1,I,J)
               IF(ABS(UIM(1,I,J)).GT.UMAX) UMAX=ABS(UIM(1,I,J))
               UI(2,I,J)=TWEIGHTU1*UIM(2,I,J)-TWEIGHTU*UIMO(2,I,J)
               IF(ABS(UIM(2,I,J)).GT.UMAX) UMAX=ABS(UIM(2,I,J))
               UI(3,I,J)=TWEIGHTU1*UIM(3,I,J)-TWEIGHTU*UIMO(3,I,J)
               IF(ABS(UIM(3,I,J)).GT.UMAX) UMAX=ABS(UIM(3,I,J))
            END DO
         END DO

C CALCULATE FLUX VARIABLES OF THE CABARET SCHEME

C I DIRECTION
      DO J=2,JJ1
         DO I=2,II
            DO K=1,NN
               FI1=ZETA_FI(K,I,J)
               IF(UI(K,I,J).GE.0) THEN
C POSITIVE SPEED
                 IF(I.GT.2) THEN
                   FI2=ZETA_FI(K,I-1,J)
                   V0=UI(K,I,J)+UI(K,I-1,J)
                   PSI_N=ZETA_NEW(K,I-1,J)
                   PSI_O=ZETA_OLD(K,I-1,J)
                   D_F=-(PSI_N-PSI_O)*2.-DT05*V0*(FI1-FI2)
                 ELSE
C BCS
                   FI2=ZETA_NEW(K,I-1,J)
                   PSI_N=ZETA_NEW(K,I-1,J)
                   PSI_O=ZETA_OLD(K,I-1,J)
                   D_F=-(PSI_N-PSI_O)
                 END IF
               ELSE
C NEGATIVE SPEED
                 IF(I.LT.II) THEN
                   FI2=ZETA_FI(K,I+1,J)
                   PSI_N=ZETA_NEW(K,I,J)
                   PSI_O=ZETA_OLD(K,I,J)
                   V0=UI(K,I,J)+UI(K,I+1,J)
                   D_F=-(PSI_N-PSI_O)*2.-DT05*V0*(FI2-FI1)
                 ELSE
C BCS
                   FI2=ZETA_NEW(K,I,J)
                   PSI_N=ZETA_NEW(K,I,J)
                   PSI_O=ZETA_OLD(K,I,J)
                   D_F=-(PSI_N-PSI_O)
                 END IF
               END IF

C NON-LINEAR CORRECTION TO ENFORCE THE MAXIMUM PRINCIPLE
               FMAX=MAX(FI1,FI2,PSI_N)-D_F
               FMIN=MIN(FI1,FI2,PSI_N)-D_F

               ZETA_FIN(K,I,J)=2.*PSI_N-FI2

               IF(ZETA_FIN(K,I,J).GT.FMAX) ZETA_FIN(K,I,J)=FMAX
               IF(ZETA_FIN(K,I,J).LT.FMIN) ZETA_FIN(K,I,J)=FMIN
            END DO
         END DO
      END DO

C J DIRECTION
      DO J=2,JJ
         DO I=2,II1
            DO K=1,NN
               FI1=ZETA_FJ(K,I,J)
               IF(UJ(K,I,J).GE.0) THEN
C POSITIVE SPEED
                 IF(J.GT.2) THEN
                   FI2=ZETA_FJ(K,I,J-1)
                   PSI_N=ZETA_NEW(K,I,J-1)
                   PSI_O=ZETA_OLD(K,I,J-1)
                   V0=UJ(K,I,J)+UJ(K,I,J-1)
                   D_F=-(PSI_N-PSI_O)*2.-DT05*V0*(FI1-FI2)
                 ELSE
C BCS
                   FI2=ZETA_NEW(K,I,J-1)
                   PSI_N=ZETA_NEW(K,I,J-1)
                   PSI_O=ZETA_OLD(K,I,J-1)
                   D_F=-(PSI_N-PSI_O)
                 END IF
               ELSE
C NEGATIVE SPEED
                 IF(J.LT.JJ) THEN
                   FI2=ZETA_FJ(K,I,J+1)
                   PSI_N=ZETA_NEW(K,I,J)
                   PSI_O=ZETA_OLD(K,I,J)
                   V0=UJ(K,I,J)+UJ(K,I,J+1)
                   D_F=-(PSI_N-PSI_O)*2.-DT05*V0*(FI2-FI1)
                 ELSE
C BCS
                   FI2=ZETA_NEW(K,I,J)
                   PSI_N=ZETA_NEW(K,I,J)
                   PSI_O=ZETA_OLD(K,I,J)
                   D_F=-(PSI_N-PSI_O)
                 END IF
               END IF

C NON-LINEAR CORRECTION TO ENFORCE THE MAXIMUM PRINCIPLE
                FMAX=MAX(FI1,FI2,PSI_N)-D_F
                FMIN=MIN(FI1,FI2,PSI_N)-D_F

                ZETA_FJN(K,I,J)=2.*PSI_N-FI2
                IF(ZETA_FJN(K,I,J).GT.FMAX) ZETA_FJN(K,I,J)=FMAX
                IF(ZETA_FJN(K,I,J).LT.FMIN) ZETA_FJN(K,I,J)=FMIN
             END DO
          END DO
      END DO

C UPDATING CELL FACE VARIABLES

      DO J=2,JJ1
         DO I=2,II
            ZETA_FI(1,I,J)=ZETA_FIN(1,I,J)
            ZETA_FI(2,I,J)=ZETA_FIN(2,I,J)
            ZETA_FI(3,I,J)=ZETA_FIN(3,I,J)
         END DO
      END DO

      DO J=2,JJ
         DO I=2,II1
            ZETA_FJ(1,I,J)=ZETA_FJN(1,I,J)
            ZETA_FJ(2,I,J)=ZETA_FJN(2,I,J)
            ZETA_FJ(3,I,J)=ZETA_FJN(3,I,J)
         END DO
      END DO

C COMPUTE FLUXES

      DO J=2,JJ1
         DO I=2,II
            FLUX_I(1,I,J)=UI(1,I,J)*ZETA_FI(1,I,J)
            FLUX_I(2,I,J)=UI(2,I,J)*ZETA_FI(2,I,J)
            FLUX_I(3,I,J)=UI(3,I,J)*ZETA_FI(3,I,J)
C LINEARIZE!
C            FLUX_I(1,I,J)=0
C            FLUX_I(2,I,J)=0
C            FLUX_I(3,I,J)=0
         END DO
      END DO

      DO J=2,JJ
         DO I=2,II1
            FLUX_J(1,I,J)=UJ(1,I,J)*ZETA_FJ(1,I,J)
            FLUX_J(2,I,J)=UJ(2,I,J)*ZETA_FJ(2,I,J)
            FLUX_J(3,I,J)=UJ(3,I,J)*ZETA_FJ(3,I,J)
C            FLUX_J(1,I,J)=0
C            FLUX_J(2,I,J)=0
C            FLUX_J(3,I,J)=0
         END DO
      END DO

C// END OF THE CABARET SOLVER; START POSTPROCESSING /////////////////////

c
c--- time-mean streamfunction
c
        IF(TIME_DAY.GT.MAX_SPIN) THEN
          K_AV=K_AV+1
          TIME_AV=TIME_AV+DT

          do j=1,jj
             do i=1,ii
                psi_av(1,i,j)=psi_av(1,i,j)+dt*psi(1,i,j)
                psi_av(2,i,j)=psi_av(2,i,j)+dt*psi(2,i,j)
                psi_av(3,i,j)=psi_av(3,i,j)+dt*psi(3,i,j)
             enddo
          enddo
        END IF
c
c--- diagnostics
c
        if(TIME_O.eq.0)then
           write(*,*)'Time (days) =',TIME_DAY

           if(idata.eq.1)then
             do j=1,jj
                do i=1,ii
                   psi1_save(j,ii-i+1)=psi(1,i,j)
                   psi2_save(j,ii-i+1)=psi(2,i,j)
                   psi3_save(j,ii-i+1)=psi(3,i,j)
                enddo
             enddo

             k_rec=k_rec+1
c             write(71,rec=k_rec) psi1_save
c             write(72,rec=k_rec) psi2_save
c             write(73,rec=k_rec) psi3_save

           endif

           call energy(nn,ii,jj,D,h,S,psi,ekin,epot)
           write(*,'(A5,F12.6,A7,F12.6)')'epot=',epot,'; ekin=',ekin

           if(idata.eq.1)then
             call write_netcdf(FILE_NAME,psi1_save,psi2_save,psi3_save,
     +     epot,ekin,ekin,ii,jj,TIME_DAY,k_rec)              !probably need to change
           endif

         endif
c
c--- save output
c
        if(TIME_S.eq.0)then
          write(*,*)'Writing to the restart file...'

          if(jj.eq.129) open(11,file='qg_129_f',form='unformatted')
          if(jj.eq.257) open(11,file='qg_257_f',form='unformatted')
          if(jj.eq.513) open(11,file='qg_513_f',form='unformatted')
          if(jj.eq.1025) open(11,file='qg_1025_f',form='unformatted')
          if(jj.eq.2049) open(11,file='qg_2049_f',form='unformatted')
          write(11) ii,jj,nn
          write(11) psi
          close(11)

c          if(jj.eq.129) open(11,file='qg_129_final',form='formatted')
c          if(jj.eq.257) open(11,file='qg_257_final',form='formatted')
c          if(jj.eq.513) open(11,file='qg_513_final',form='formatted')
c          if(jj.eq.1025) open(11,file='qg_1025_final',form='formatted')
c          if(jj.eq.2049) open(11,file='qg_2049_final',form='formatted')
c
c          DO J=1,JJ
c             DO I=1,II
c                DO K=1,NN
c                   WRITE(11,19) ZETA_OLD(K,I,J),RES(K,I,J)
c     &                                         ,RES_BETA_O(K,I,J)
c                END DO
c             END DO
c          END DO
c          CLOSE(11)

C OUTPUTING MEAN STREAMFUNCTIONS

          IF(TIME_DAY.GT.MAX_SPIN) THEN
            if(jj.eq.129) open(11,file='qg_129_av.d',form='unformatted')
            if(jj.eq.257) open(11,file='qg_257_av.d',form='unformatted')
            if(jj.eq.513) open(11,file='qg_513_av.d',form='unformatted')
            if(jj.eq.1025) open(11,file='qg_1025_av.d'
     &                                              ,form='unformatted')
            if(jj.eq.2049) open(11,file='qg_2049_av.d'
     &                                              ,form='unformatted')
            write(11) TIME_AV
            write(11) ii,jj,nn
            write(11) psi_av
            close(11)
          END IF
      endif

C // UPDATING CENTRE-CELL VARIABLES/////

      DO J=2,JJ1
         DO I=2,II1
            RES(1,I,J)=-(FLUX_I(1,I+1,J)-FLUX_I(1,I,J)
     &                        +FLUX_J(1,I,J+1)-FLUX_J(1,I,J))
            RES(2,I,J)=-(FLUX_I(2,I+1,J)-FLUX_I(2,I,J)
     &                        +FLUX_J(2,I,J+1)-FLUX_J(2,I,J))
            RES(3,I,J)=-(FLUX_I(3,I+1,J)-FLUX_I(3,I,J)
     &                        +FLUX_J(3,I,J+1)-FLUX_J(3,I,J))
         END DO
      END DO

      DO J=2,JJ1
         DO I=2,II1
            ZETA_NEW(1,I,J)=ZETA_NEW(1,I,J)+DT05*RES(1,I,J)
            ZETA_NEW(2,I,J)=ZETA_NEW(2,I,J)+DT05*RES(2,I,J)
            ZETA_NEW(3,I,J)=ZETA_NEW(3,I,J)+DT05*RES(3,I,J)
         END DO
      END DO

C UPDATING ZETA_OLD

      DO J=1,JJ
         DO I=1,II
            ZETA_OLD(1,I,J)=ZETA_NEW(1,I,J)
            UIMO(1,I,J)=UIM(1,I,J)
            UJMO(1,I,J)=UJM(1,I,J)
            ZETA_OLD(2,I,J)=ZETA_NEW(2,I,J)
            UIMO(2,I,J)=UIM(2,I,J)
            UJMO(2,I,J)=UJM(2,I,J)
            ZETA_OLD(3,I,J)=ZETA_NEW(3,I,J)
            UIMO(3,I,J)=UIM(3,I,J)
            UJMO(3,I,J)=UJM(3,I,J)
         END DO
      END DO

      DO J=1,JJ
         DO I=1,II
            RES_BETA_O(1,I,J)=RES_BETA(1,I,J)
            RES_BETA_O(2,I,J)=RES_BETA(2,I,J)
            RES_BETA_O(3,I,J)=RES_BETA(3,I,J)
         END DO
      END DO

      IF(TIME_DAY.LT.MAX_TIME+MAX_SPIN) GOTO 1001
c
c#################################################
c
c      if(idata.eq.1)then
c        close(71)
c        close(72)
c        close(73)
c      endif
c
c--- save the results
c
      if(jj.eq.129) open(11,file='qg_129_f',form='unformatted')
      if(jj.eq.257) open(11,file='qg_257_f',form='unformatted')
      if(jj.eq.513) open(11,file='qg_513_f',form='unformatted')
      if(jj.eq.1025) open(11,file='qg_1025_f',form='unformatted')
      if(jj.eq.2049) open(11,file='qg_2049_f',form='unformatted')
      write(11) ii,jj,nn
      write(11) psi
      close(11)

c      if(jj.eq.129) open(11,file='qg_129_final',form='formatted')
c      if(jj.eq.257) open(11,file='qg_257_final',form='formatted')
c      if(jj.eq.513) open(11,file='qg_513_final',form='formatted')
c      if(jj.eq.1025) open(11,file='qg_1025_final',form='formatted')
c      if(jj.eq.2049) open(11,file='qg_2049_final',form='formatted')
c
c      DO J=1,JJ
c         DO I=1,II
c            DO K=1,NN
c               WRITE(11,19) ZETA_OLD(K,I,J),RES(K,I,J)
c     &                                         ,RES_BETA_O(K,I,J)
c            END DO
c         END DO
c      END DO
c      CLOSE(11)
c


c--- time-mean streamfunction
c
      do j=1,jj
         do i=1,ii
            do n=1,nn
               psi_av(n,i,j)=psi_av(n,i,j)/TIME_AV
            enddo
         enddo
      enddo

      if(jj.eq.129) open(11,file='qg_129_av.d',form='unformatted')
      if(jj.eq.257) open(11,file='qg_257_av.d',form='unformatted')
      if(jj.eq.513) open(11,file='qg_513_av.d',form='unformatted')
      if(jj.eq.1025) open(11,file='qg_1025_av.d',form='unformatted')
      if(jj.eq.2049) open(11,file='qg_2049_av.d',form='unformatted')
      write(11) ii,jj,nn
      write(11) psi_av
      close(11)

c      if(ips.eq.0)then
c        call clsgks
c      else
c        call gdawk(1)
c        call gclwk(1)
c        call gclks
c      endif

      stop
      end
C
C============================================================================
C
      subroutine boundary_condition(islip,alpha,nn,ii,jj,rel,psi)
      implicit none
      integer nn,ii,jj,n,i,j,islip
      real*8 rel(nn,ii,jj),psi(nn,ii,jj)
     & ,alpha,cff,alpha_eff,cff_eff

      if(islip.eq.1)then         ! NO-SLIP
        do j=2,jj-1
           rel(1,1 ,j)=-3.5*psi(1,1   ,j)
     &                 +4.0*psi(1,2   ,j)
     &                 -0.5*psi(1,3   ,j)
           rel(1,ii,j)=-3.5*psi(1,ii  ,j)
     &                 +4.0*psi(1,ii-1,j)
     &                 -0.5*psi(1,ii-2,j)

           rel(2,1 ,j)=-3.5*psi(2,1   ,j)
     &                 +4.0*psi(2,2   ,j)
     &                 -0.5*psi(2,3   ,j)
           rel(2,ii,j)=-3.5*psi(2,ii  ,j)
     &                 +4.0*psi(2,ii-1,j)
     &                 -0.5*psi(2,ii-2,j)

           rel(3,1 ,j)=-3.5*psi(3,1   ,j)
     &                 +4.0*psi(3,2   ,j)
     &                 -0.5*psi(3,3   ,j)
           rel(3,ii,j)=-3.5*psi(3,ii  ,j)
     &                 +4.0*psi(3,ii-1,j)
     &                 -0.5*psi(3,ii-2,j)
        enddo
        do i=2,ii-1
           rel(1,i,1 )=-3.5*psi(1,i,1   )
     &                 +4.0*psi(1,i,2   )
     &                 -0.5*psi(1,i,3   )
           rel(1,i,jj)=-3.5*psi(1,i,jj  )
     &                 +4.0*psi(1,i,jj-1)
     &                 -0.5*psi(1,i,jj-2)

           rel(2,i,1 )=-3.5*psi(2,i,1   )
     &                 +4.0*psi(2,i,2   )
     &                 -0.5*psi(2,i,3   )
           rel(2,i,jj)=-3.5*psi(2,i,jj  )
     &                 +4.0*psi(2,i,jj-1)
     &                 -0.5*psi(2,i,jj-2)

           rel(3,i,1 )=-3.5*psi(3,i,1   )
     &                 +4.0*psi(3,i,2   )
     &                 -0.5*psi(3,i,3   )
           rel(3,i,jj)=-3.5*psi(3,i,jj  )
     &                 +4.0*psi(3,i,jj-1)
     &                 -0.5*psi(3,i,jj-2)
        enddo

c-------------------------------- psi_nn-alpha*psi_n=0
        cff=alpha/(alpha+3.D0)
        do j=2,jj-1
           rel(1,1 ,j)=cff*rel(1,1 ,j)
           rel(1,ii,j)=cff*rel(1,ii,j)
           rel(2,1 ,j)=cff*rel(2,1 ,j)
           rel(2,ii,j)=cff*rel(2,ii,j)
           rel(3,1 ,j)=cff*rel(3,1 ,j)
           rel(3,ii,j)=cff*rel(3,ii,j)
        enddo
        do i=2,ii-1
           rel(1,i,1 )=cff*rel(1,i,1 )
           rel(1,i,jj)=cff*rel(1,i,jj)
           rel(2,i,1 )=cff*rel(2,i,1 )
           rel(2,i,jj)=cff*rel(2,i,jj)
           rel(3,i,1 )=cff*rel(3,i,1 )
           rel(3,i,jj)=cff*rel(3,i,jj)
        enddo
      else                       ! FREE-SLIP
        do j=2,jj-1
           rel(1,1,j)=0.
           rel(1,ii,j)=0.
           rel(2,1,j)=0.
           rel(2,ii,j)=0.
           rel(3,1,j)=0.
           rel(3,ii,j)=0.
        enddo
        do i=2,ii-1
           rel(1,i,1)=0.
           rel(1,i,jj)=0.
           rel(2,i,1)=0.
           rel(2,i,jj)=0.
           rel(3,i,1)=0.
           rel(3,i,jj)=0.
        enddo
      endif

      return
      end
c
c-----------------------------------------------------------------------------
c
      subroutine correction_functions(nn,ii,jj,nw,eig,wsave,SS
     &                                ,eig_aux,wsave_aux,corr_aux,corr)
      implicit none
      integer nn,ii,jj,n,i,j,nw,ii2,jj2,n1
      real*8 eig(nn,ii,jj),wsave(nn,nw),SS(nn),corr(nn,ii,jj)
     & ,eig_aux(ii,jj),wsave_aux(nw),corr_aux(ii,jj)
     & ,sum

      ii2=ii-2
      jj2=jj-2

      do 100 n=2,nn
         do j=2,jj-1
            do i=2,ii-1
               corr_aux(i,j)=1.
            enddo
         enddo
         do j=1,jj
            corr_aux(1,j)=0.
            corr_aux(ii,j)=0.
         enddo
         do i=2,ii-1
            corr_aux(i,1)=0.
            corr_aux(i,jj)=0.
         enddo

         do j=1,jj
            do i=1,ii
               eig_aux(i,j)=eig(n,i,j)
            enddo
         enddo
         do n1=1,nw
            wsave_aux(n1)=wsave(n,n1)
         enddo

         call helm(corr_aux(2,2),ii,ii2,jj2,eig_aux,wsave_aux)

         sum=0.
         do j=2,jj-1
            do i=2,ii-1
               sum=sum+corr_aux(i,j)
            enddo
         enddo
         sum=sum/dfloat((ii-1)*(jj-1))

         do j=1,jj
            do i=1,ii
               corr(n,i,j)=(1.D0+SS(n)*corr_aux(i,j))/(1.D0+SS(n)*sum)
            enddo
         enddo
100   continue

      return
      end
c
c-----------------------------------------------------------------------------
c
      subroutine rel_from_psi(islip,alpha,nn,ii,jj,psi,rel)
      implicit none
      integer nn,jj,ii,islip,n,j,i
      real*8 psi(nn,ii,jj),rel(nn,ii,jj)
     & ,alpha

      do j=2,jj-1
         do i=2,ii-1
            do n=1,nn
               rel(n,i,j)=-4.*psi(n,i,j)
     &                     +psi(n,i-1,j)+psi(n,i,j+1)
     &                     +psi(n,i+1,j)+psi(n,i,j-1)
            enddo
         enddo
      enddo

      call boundary_condition(islip,alpha,nn,ii,jj,rel,psi)

      return
      end
c
c----------------------------------------------------------------------
c
      subroutine zeta_from_psi_and_rel(nn,ii,jj,S,psi,rel,zeta)
      implicit none
      integer nn,ii,jj,n,i,j
      real*8 psi(nn,ii,jj),rel(nn,ii,jj),zeta(nn,ii,jj),S(nn,2)

      do j=1,jj
         do i=1,ii
            do n=1,nn
               if(n.eq.1)then
                 zeta(n,i,j)=rel(n,i,j)
     &            -S(n,2)*(psi(n,i,j)-psi(n+1,i,j))
               elseif(n.eq.nn)then
                    zeta(n,i,j)=rel(n,i,j)
     &            -S(n,1)*(psi(n,i,j)-psi(n-1,i,j))
               else
                 zeta(n,i,j)=rel(n,i,j)
     &            -S(n,1)*(psi(n,i,j)-psi(n-1,i,j))
     &            -S(n,2)*(psi(n,i,j)-psi(n+1,i,j))
               endif
            enddo
         enddo
      enddo

      return
      end
c
c----------------------------------------------------------------------
c
      subroutine helm(zeta,ii,m,n,eig,wsave)
C      ==============================================             C
C                                                                 C
C                                                                 C
C       This subroutine uses a Fast Sine Transform algorithm      C
C       to solve  the  HELMHOLTZ  EQUATION:                       C
C                                                                 C
C                 PSI  + PSI   + S*PSI = ZETA                     C
C                    xx     yy                                    C
C                                                                 C
C       using  the Dirichlet's boundary  conditions:              C
C                                                                 C
C         PSI(0,y) = PSI(Lx,y) = PSI(x,0) = PSI(x,Ly) = 0         C
C                                                                 C
C       in a rectangular domain  with dimensions Lx, Ly.          C
C                      ---------------------------                C
C       Note that the the array ZETA inputs the right hand        C
C       side of the equation and  returns the solution, PSI.      C
C                                                                 C
C       WSAVE is an working array; if II2=JJ2 then WSAVE does     C
C       not need to be modified by calling SINTI again.           C
C                    ********************                         C
C       THIS SUBROUTINE CALLS SUBROUTINES SINT AND SINTI,  FROM   C
C       NCAR'S PACKAGE: FFTPACK (LINK WITH LIBRARY MYLIBRY.OLB)   C
c-----------------------------------------------------------------
      implicit none
      integer ii,m,n,n1,m1,k,j
      real*8 zeta(ii,1), eig(ii,1), wsave(1)
     & ,pi,cff

      data pi/3.14159265358979323846D0/

      n1 = n+1
      m1 = m+1
C=================================================================C
C         PART I - COMPUTES TRANSFORM OF ZETA(I,J)                C
C=================================================================C
      cff=0.5D0/float(m1)
      do k=1,n
         do j=1,m
            zeta(j,k)=cff*zeta(j,k)
         enddo
         call sint(m,zeta(1,k),wsave)   ! FORWARD TRANSFORM
      enddo

C--- GAUSSIAN ELIMINATION - RIGHT HAND SIDE
      do k=2,n
         do j=1,m
            zeta(j,k)=zeta(j,k)-zeta(j,k-1)*eig(j,k-1)
         enddo
      enddo

      do j=1,m
         zeta(j,n)=zeta(j,n)*eig(j,n)
      enddo

      do k=n-1,1,-1
         do j=1,m
            zeta(j,k)=(zeta(j,k)-zeta(j,k+1))*eig(j,k)
         enddo
      enddo

C--- BACK TRANSFORM:
      do k=1,n
         call sint(m,zeta(1,k),wsave)
         zeta(m1,k)=0.D0
      enddo

      return
      end
c
c----------------------------------------------------------------------------
c
      subroutine poinit(S,m,n,l,eig,wsave)  ! INITIALIZATION ROUTINE
      implicit none
      integer m,n,l,j,k
      real*8 S,eig(l,n),wsave(*),pi,cff

      data pi/3.14159265358979323846D0/
c
c--- INITIALIZE FFT ROUTINE
c
      call sinti(m,wsave)
c
c--- EIGENVALUES OF TRIDIAGONAL MATRIX
c
      cff=pi/(2.D0*dfloat(m+1))
      do k=1,n
         do j=1,m
            eig(j,k)=-2.D0-S-4.D0*dsin(dfloat(j)*cff)**2
         enddo
      enddo
c
c--- GAUSSIAN ELIMINATION OF TRI-DIAGONAL SYSTEM - LEFT HAND SIDE
c
      do k=2,n
         do j=1,m
            eig(j,k)=eig(j,k)-1.D0/eig(j,k-1)
         enddo
      enddo

      do k=1,n
         do j=1,m
            eig(j,k)=1.D0/eig(j,k)
         enddo
      enddo

      return
      end
c
c----------------------------------------------------------------------------
c
      subroutine sinti(n,wsave)
      implicit none
      integer n,ns2,np1,k
      real*8 wsave(1),pi,dt

      data pi /3.14159265358979D0/
c
      if(n.le.1) return
      ns2=n/2
      np1=n+1
      dt=pi/dfloat(np1)
      do k=1,ns2
         wsave(k)=2.D0*dsin(dfloat(k)*dt)
      enddo
      call rffti(np1,wsave(ns2+1))

      return
      end
c
c----------------------------------------------------------------------------
c
      subroutine sint(n,x,wsave)
      implicit none
      integer n,np1,iw1,iw2,iw3
      real*8 x(1),wsave(1)
c
      np1=n+1
      iw1=n/2+1
      iw2=iw1+np1
      iw3=iw2+np1
      call sint1(n,x,wsave,wsave(iw1),wsave(iw2),wsave(iw3))

      return
      end
c
c----------------------------------------------------------------------------
c
      subroutine sint1(n,war,was,xh,x,ifac)
      implicit none
      integer n,i,np1,ns2,k,kc,modn
      real*8 war(3582),was(3582),x(3582),xh(3582),ifac(3582)
     & ,xhold,sqrt3,t1,t2

      data sqrt3/1.73205080756888D0/
c
      do i=1,n
         xh(i)=war(i)
         war(i)=x(i)
      enddo

      if (n-2) 101,102,103
  101 xh(1)=xh(1)+xh(1)
      go to 106
  102 xhold=sqrt3*(xh(1)+xh(2))
      xh(2)=sqrt3*(xh(1)-xh(2))
      xh(1)=xhold
      go to 106
  103 np1=n+1
      ns2=n/2

      x(1)=0.D0
      do k=1,ns2
         kc=np1-k
         t1=xh(k)-xh(kc)
         t2=was(k)*(xh(k)+xh(kc))
         x(k+1)=t1+t2
         x(kc+1)=t2-t1
      enddo

      modn=mod(n,2)
      if (modn.ne.0) x(ns2+2)=4.D0*xh(ns2+1)
      call rfftf1(np1,x,xh,war,ifac)

      xh(1)=.5D0*x(1)
      do i=3,n,2
         xh(i-1)=-x(i)
         xh(i)=xh(i-2)+x(i-1)
      enddo

      if (modn.ne.0) go to 106
      xh(n)=-x(n+1)
106   continue

      do i=1,n
         x(i)=war(i)
         war(i)=xh(i)
      enddo

      return
      end
c
c----------------------------------------------------------------------------
c
      subroutine rffti(n,wsave)
      implicit none
      integer n
      real*8 wsave(1)
c
      if (n.eq.1) return
      call rffti1(n,wsave(n+1),wsave(2*n+1))

      return
      end
c
c----------------------------------------------------------------------------
c
      subroutine rffti1(n,wa,ifac)
      implicit none
      integer n,nl,nf,j,ntry,nq,nr,i,ib,is,nfm1,l1,k1,ip,ld,l2
     & ,ido,ipm,ii
      real*8 arg,argh,argld,wa(1),ifac(13),ntryh(4)
     & ,tpi,fi

      data ntryh(1),ntryh(2),ntryh(3),ntryh(4)/4.D0,2.D0,3.D0,5.D0/
c
      nl=n
      nf=0
      j=0
  101 j=j+1
      if (j-4) 102,102,103
  102 ntry=ntryh(j)
      go to 104
  103 ntry=ntry+2
  104 nq=nl/ntry
      nr=nl-ntry*nq
      if (nr) 101,105,101
  105 nf=nf+1
      ifac(nf+2)=ntry
      nl=nq
      if (ntry.ne.2) go to 107
      if (nf.eq.1) go to 107

      do i=2,nf
         ib=nf-i+2
         ifac(ib+2)=ifac(ib+1)
      enddo

      ifac(3)=2
  107 if (nl.ne.1) go to 104
      ifac(1)=n
      ifac(2)=nf
      tpi=6.28318530717959D0
      argh=tpi/dfloat(n)
      is=0
      nfm1=nf-1
      l1=1
      if (nfm1.eq.0) return

      do 110 k1=1,nfm1
         ip=ifac(k1+2)
         ld=0
         l2=l1*ip
         ido=n/l2
         ipm=ip-1
         do j=1,ipm
            ld=ld+l1
            i=is
            argld=dfloat(ld)*argh
            fi=0.D0
            do ii=3,ido,2
               i=i+2
               fi=fi+1.D0
               arg=fi*argld
               wa(i-1)=dcos(arg)
               wa(i)=dsin(arg)
            enddo
            is=is+ido
         enddo
         l1=l2
110   continue

      return
      end
c
c----------------------------------------------------------------------------
c
      subroutine rfftf(n,r,wsave)
      implicit none
      integer n
      real*8 r(1),wsave(1)

      if (n.eq.1) return
      call rfftf1(n,r,wsave,wsave(n+1),wsave(2*n+1))

      return
      end
c
c----------------------------------------------------------------------------
c
      subroutine rfftf1(n,c,ch,wa,ifac)
      implicit none
      integer n,nf,na,l2,iw,k1,kh,ip,l1,ido,idl1,ix2,ix3,i
      real*8 ch(1),c(1),wa(1),ifac(13)
c
      nf=ifac(2)
      na=1
      l2=n
      iw=n
      do 111 k1=1,nf
         kh=nf-k1
         ip=ifac(kh+3)
         l1=l2/ip
         ido=n/l2
         idl1=ido*l1
         iw=iw-(ip-1)*ido
         na=1-na
         if (ip.ne.4) go to 102
         ix2=iw+ido
         ix3=ix2+ido
         if (na.ne.0) go to 101
         call radf4(ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
         go to 110
  101    call radf4(ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
         go to 110

  102    if (ip.ne.2) stop 'add "radf..." subroutines!!!'
         IF (na.ne.0) go to 103
         call radf2(ido,l1,c,ch,wa(iw))
         go to 110
  103    call radf2(ido,l1,ch,c,wa(iw))
         go to 110

  110    l2=l1
  111 continue
      IF (na.eq.1) return

      do i=1,n
         c(i)=ch(i)
      enddo

      return
      end
c
c----------------------------------------------------------------------------
c
      subroutine radf2(ido,l1,cc,ch,wa1)
      implicit none
      integer k,l1,ido,idp2,i,ic
      real*8 ch(ido,2,l1),cc(ido,l1,2),wa1(1),tr2,ti2
c
      do k=1,l1
         ch(1,1,k)=cc(1,k,1)+cc(1,k,2)
         ch(ido,2,k)=cc(1,k,1)-cc(1,k,2)
         ch(ido,2,k)=cc(1,k,1)-cc(1,k,2)
      enddo

      if (ido-2) 107,105,102
  102 idp2=ido+2

      do k=1,l1
         do i=3,ido,2
            ic=idp2-i
            tr2=wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            ti2=wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            ch(i,1,k)=cc(i,k,1)+ti2
            ch(ic,2,k)=ti2-cc(i,k,1)
            ch(i-1,1,k)=cc(i-1,k,1)+tr2
            ch(ic-1,2,k)=cc(i-1,k,1)-tr2
         enddo
      enddo
      if (mod(ido,2).eq.1) return
105   continue

      do k=1,l1
         ch(1,2,k)=-cc(ido,k,2)
         ch(ido,1,k)=cc(ido,k,1)
      enddo
107   continue

      return
      end
c
c----------------------------------------------------------------------------
c
      subroutine radf4(ido,l1,cc,ch,wa1,wa2,wa3)
      implicit none
      integer k,l1,ido,idp2,i,ic
      real*8 cc(ido,l1,4),ch(ido,4,l1),wa1(1),wa2(1),wa3(1)
     & ,cr2,ci2,cr3,ci3,cr4,ci4,tr1,tr2,tr3,tr4,ti1,ti2,ti3,ti4
     & ,hsqt2
      data hsqt2/.7071067811865475D0/
c
      do k=1,l1
         tr1=cc(1,k,2)+cc(1,k,4)
         tr2=cc(1,k,1)+cc(1,k,3)
         ch(1,1,k)=tr1+tr2
         ch(ido,4,k)=tr2-tr1
         ch(ido,2,k)=cc(1,k,1)-cc(1,k,3)
         ch(1,3,k)=cc(1,k,4)-cc(1,k,2)
      enddo

      if (ido-2) 107,105,102
  102 idp2=ido+2

      do k=1,l1
         do i=3,ido,2
            ic=idp2-i
            cr2=wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            ci2=wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            cr2=wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            ci2=wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            cr3=wa2(i-2)*cc(i-1,k,3)+wa2(i-1)*cc(i,k,3)
            ci3=wa2(i-2)*cc(i,k,3)-wa2(i-1)*cc(i-1,k,3)
            cr4=wa3(i-2)*cc(i-1,k,4)+wa3(i-1)*cc(i,k,4)
            ci4=wa3(i-2)*cc(i,k,4)-wa3(i-1)*cc(i-1,k,4)
            tr1=cr2+cr4
            tr4=cr4-cr2
            ti1=ci2+ci4
            ti4=ci2-ci4
            ti2=cc(i,k,1)+ci3
            ti3=cc(i,k,1)-ci3
            tr2=cc(i-1,k,1)+cr3
            tr3=cc(i-1,k,1)-cr3
            ch(i-1,1,k)=tr1+tr2
            ch(ic-1,4,k)=tr2-tr1
            ch(i,1,k)=ti1+ti2
            ch(ic,4,k)=ti1-ti2
            ch(i-1,3,k)=ti4+tr3
            ch(ic-1,2,k)=tr3-ti4
            ch(i,3,k)=tr4+ti3
            ch(ic,2,k)=tr4-ti3
         enddo
      enddo
      if (mod(ido,2).eq.1) return
  105 continue

      do k=1,l1
         ti1=-hsqt2*(cc(ido,k,2)+cc(ido,k,4))
         tr1=hsqt2*(cc(ido,k,2)-cc(ido,k,4))
         ch(ido,1,k)=tr1+cc(ido,k,1)
         ch(ido,3,k)=cc(ido,k,1)-tr1
         ch(1,2,k)=ti1-cc(ido,k,3)
         ch(1,4,k)=ti1+cc(ido,k,3)
      enddo
107   continue

      return
      end
c
c--------------------------------------------------------------------------
c
      subroutine mean_correct(nn,ii,jj,phi,phi_C,corr)
      implicit none
      integer nn,ii,jj,n,i,j
      real*8 phi(nn,ii,jj),phi_C(nn),corr(nn,ii,jj)
     & ,sum1,sum2,sum3,cff

      if(nn.ne.3) stop 'Change nn in mean_correct...'
      sum1=0.
      sum2=0.
      sum3=0.
      do j=2,jj-1
         do i=2,ii-1
            sum1=sum1+phi(1,i,j)
            sum2=sum2+phi(2,i,j)
            sum3=sum3+phi(3,i,j)
         enddo
      enddo
      do j=2,jj-1
         sum1=sum1+0.5D0*(phi(1,1,j)+phi(1,ii,j))
         sum2=sum2+0.5D0*(phi(2,1,j)+phi(2,ii,j))
         sum3=sum3+0.5D0*(phi(3,1,j)+phi(3,ii,j))
      enddo
      do i=2,ii-1
         sum1=sum1+0.5D0*(phi(1,i,1)+phi(1,i,jj))
         sum2=sum2+0.5D0*(phi(2,i,1)+phi(2,i,jj))
         sum3=sum3+0.5D0*(phi(3,i,1)+phi(3,i,jj))
      enddo
      sum1=sum1+phi(1,1,1)
      sum2=sum2+phi(2,1,1)  ! 0.25*(* + * + * + *)
      sum3=sum3+phi(3,1,1)
      cff=1./dfloat((ii-1)*(jj-1))
      phi_C(1)=cff*sum1
      phi_C(2)=cff*sum2
      phi_C(3)=cff*sum3

      do j=1,jj
         do i=1,ii
            phi(1,i,j)=phi(1,i,j)-phi_C(1)
            phi(2,i,j)=phi(2,i,j)-phi_C(2)*corr(2,i,j)
            phi(3,i,j)=phi(3,i,j)-phi_C(3)*corr(3,i,j)
         enddo
      enddo

      return
      end
c
c--------------------------------------------------------------------------
c
      subroutine energy(nn,ii,jj,D,h,S,psi,ekin,epot)
      implicit none
      integer nn,ii,jj,n,i,j
      real*8 H(nn),S(nn,2),psi(nn,ii,jj),factor(1000)
     & ,D,ekin,epot,cff

      do n=1,nn-1
         factor(n)=0.25D0*(S(n,2)*h(n)+S(n+1,1)*h(n+1))/D
      enddo

      epot=0.
      do j=2,jj-1
         do i=2,ii-1
            do n=1,nn-1
               epot=epot+factor(n)*(psi(n,i,j)-psi(n+1,i,j))**2
            enddo
         enddo
      enddo

      do j=1,jj,jj-1
         do i=1,ii,ii-1
            do n=1,nn-1
               epot=epot+.5D0*factor(n)*(psi(n,i,j)-psi(n+1,i,j))**2
            enddo
         enddo
      enddo

      ekin=0.
      do j=2,jj
         do i=2,ii
            do n=1,nn
               ekin=ekin+.25*h(n)*(
     &  (psi(n,i-1,j)+psi(n,i,j)-psi(n,i-1,j-1)-psi(n,i,j-1))**2
     & +(psi(n,i,j-1)+psi(n,i,j)-psi(n,i-1,j-1)-psi(n,i-1,j))**2
     &                            )
            enddo
         enddo
      enddo
      ekin=.5D0*ekin/(D*dfloat((ii-1)*(jj-1)))
      epot =epot/dfloat((ii-1)*(jj-1))
c      ekin=1.E-3*ekin
c      epot=1.E-5*epot

      return
      end
