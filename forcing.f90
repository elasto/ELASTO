MODULE forcing

  
  USE precision_tools
  USE wavenumber_tools
  USE datalayout
  USE param
  USE parallel_tools
  USE mpi
  USE communicators_tools
  USE data
  USE transforms_tools
  USE random_tools

  IMPLICIT NONE

  !> Alvelius forcing
  INTEGER,  PRIVATE     :: ialvelius=0
  !> Brandenburg forcing
  INTEGER,  PRIVATE     :: ibrandenburg=0
  !> Kolmogorov forcing
  INTEGER,  PRIVATE     :: ikolmo=0
  !> Linear forcing
  INTEGER,  PRIVATE     :: ilinear=0
  !> Sinus forcing
  INTEGER,  PRIVATE     :: isinus=0
  !> Robert forcing
  INTEGER,  PRIVATE     :: irobert=0
  !> Helicoidal forcing 
  INTEGER,  PRIVATE     :: ihel=0
  !> anisotrop forcing 
  INTEGER,  PRIVATE     :: ianis=0
  !> Helicity but no energy forcing 
  INTEGER,  PRIVATE     :: inoEner=0
  !> For Alvelius forcing
  REAL(WP)              :: A,kf,ka,kb,c
  !> For Alvelius forcing
  REAL(WP)              :: b_anis

  !> For Brandenburg forcing
  REAL(WP)             :: f0, kn,PropForc
  REAL(WP)             :: lxfmin,lxfmax
  REAL(WP)             :: lyfmin,lyfmax
  REAL(WP)             :: lzfmin,lzfmax
  !> For Kolmo forcing
  INTEGER              :: mode_k


CONTAINS

!> Initialise forcing
!! @author Mouloud Kessar
!!    @param[in,out]  iforce    = integer equal to 1 if forcing
!!    @param[in]  AUk           = Storage for the Spectral values of U
!!    @param[in]  AVk           = Storage for the Spectral values of V  
!!    @param[in]  AWk           = Storage for the Spectral values of W  
!!    @param[in]  AVelWN        = Velocity wavenumbers  
!!    @param[in]  mu            = viscosity
!!    @param[in]  spec_rank     = rank inside the spectral communicator
SUBROUTINE forcing_init(iforce, AUk, AVk, AWk, AVelWN, mu, spec_rank)

    implicit none

    INTEGER, INTENT(IN)                     :: spec_rank
    INTEGER, INTENT(INOUT)                  :: iforce
    TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN)   :: AUk, AVk, AWk
    TYPE(WaveNumbers), INTENT(IN)           :: AVelWN
    REAL(WP), INTENT(in)                    :: mu
    REAL(WP)                                :: Pi

    CHARACTER(LEN=32)                       :: forcing
    CHARACTER(LEN=32)                       :: helicity
    CHARACTER(LEN=32)                       :: anisotrop 
    CHARACTER(LEN=32)                       :: noEnergy
 

    CALL parser_read('Forcing Method',forcing)

    Pi = acos(-1.0_WP)
    iforce = 1

    IF (forcing .EQ. 'none') THEN
       iforce = 0
    ELSE IF (forcing .EQ. 'alvelius') THEN
       ialvelius=1
       CALL parser_read('Helicoidal',helicity)
       CALL parser_read('anisotropic',anisotrop)
       IF (anisotrop .EQ. 'yes') THEN
          ianis = 1
          CALL parser_read('anisotropia parameter',b_anis)
       ELSE IF (anisotrop .EQ. 'no') THEN
          ianis = 0
       ELSE
          IF (spec_rank.EQ.0) WRITE(6,'(a)')'[WARNING] In Alvelius forcing: Unknown anisotropic type forcing ['//TRIM(forcing) &
             & //'] : defaulting to no anisotropic.'
          ianis=0
       END IF

       IF (helicity .EQ. 'yes') THEN
          ihel = 1
!           CALL forcing_init_alv(AUk,AVk,AWk,AVelWN,mu,spec_rank)
       ELSE IF (helicity .EQ. 'no') THEN
          ihel = 0
!           CALL forcing_init_alv(AUk,AVk,AWk,AVelWN,mu,spec_rank)
       ELSE
          IF (spec_rank.EQ.0) WRITE(6,'(a)')'[WARNING] In Alvelius forcing: Unknown helicity type forcing ['//TRIM(forcing) &
             & //'] : defaulting to no helicity.'
          ihel=0
       END IF
          CALL forcing_init_alv(AUk,AVk,AWk,AVelWN,mu,spec_rank)
       CALL parser_read('No Energy',noEnergy)

            IF (noEnergy .EQ. 'yes') THEN
               inoEner = 1
            ELSE IF (noEnergy .EQ. 'no') THEN
               inoEner = 0
!        ELSE
!           IF (spec_rank.EQ.0) WRITE(6,'(a)')'[WARNING] In Alvelius forcing: Unknown helicity type forcing ['//TRIM(forcing) &
!              & //'] : defaulting to no helicity.'
!           ihel=0
            END IF

    ELSE IF (forcing .EQ. 'brandenburg') THEN
       CALL parser_read('Helicoidal',helicity)
       ibrandenburg=1
       IF (helicity .EQ. 'yes') THEN
          ihel=1
          CALL parser_read('Peak forcing',kn)
          CALL parser_read('Forcing amplitude',f0)
          CALL parser_read('Forcing area',PropForc)

          lxfmin = Pi*(1.0_WP - PropForc)
          lxfmax = 2*Pi-lxfmin

          lyfmin = lxfmin
          lyfmax = lxfmax
          lzfmin = lxfmin
          lzfmax = lxfmax
          IF (spec_rank.EQ.0) THEN
             WRITE(6,'(a,2(x,f14.6))')'[INFO] internal forcing area(Brand Hel):lxfmin,lxfmax',lxfmin,lxfmax
          END IF
          
       ELSE IF (helicity .EQ. 'no') THEN
          ihel=0
          IF (spec_rank.EQ.0) THEN
             WRITE(6,'(a)')'[ERROR] No helicity brandenburg forcing not yet implemented'
          END IF
          STOP
       ELSE
          IF (spec_rank.EQ.0) WRITE(6,'(a)')'[WARNING] In Brandenburg forcing: Unknown helicity type forcing ['//TRIM(forcing) &
             & //'] : defaulting to helicity.'
          ihel=1
       END IF
    ELSE IF (forcing .EQ.'kolmogorov') THEN
       ikolmo = 1
       CALL parser_read('Additional mode',mode_k)
    ELSE IF (forcing .EQ.'linear') THEN
       CALL parser_read('Forcing amplitude',f0)
       ilinear = 1
    ELSE IF (forcing .EQ.'sinus') THEN
       CALL parser_read('Forcing amplitude',f0)
       isinus = 1
   ELSE IF (forcing .EQ.'robert') THEN
       irobert = 1
     ELSE
       IF(spec_rank.EQ.0) WRITE(6,'(a)')'[WARNING] Unknown forcing method ['//TRIM(forcing) &
          & //'] : defaulting to no forcing.'
       iforce=0
    ENDIF

    !print*,'forcing type :',iforce,ialvelius_hel,ialvelius_nohel,ibrandhel,ibrandenburg
    !print*,'kn,f0',kn,f0

END SUBROUTINE forcing_init


!> Compute forcing
!! @author Mouloud Kessar
!!    @param[in]  AUk           = Storage for the Spectral values of U
!!    @param[in]  AVk           = Storage for the Spectral values of V  
!!    @param[in]  AWk           = Storage for the Spectral values of W  
!!    @param[in]  AVelWN        = Velocity wavenumbers  
!!    @param[in,out]  Anlx      = Nonlinear term in spectral space for U 
!!    @param[in,out]  Anly      = Nonlinear term in spectral space for V 
!!    @param[in,out]  Anlz      = Nonlinear term in spectral space for W 
!!    @param[in]  dt            = time step
!!    @param[in]  spec_rank     = rank inside the spectral communicator
!!    @param[in]  AU, AV, AW    = x, y and z-velocity in physical space, needed for linear forcing
!! @details
!!    Interface for the different forcing method.
SUBROUTINE force_compute(AUk,AVk,AWk,AVelWN,Anlx,Anly,Anlz,dt,spec_rank,AU,AV,AW)

    IMPLICIT NONE 

    TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN)   :: AUk, AVk, AWk
    TYPE(WaveNumbers), INTENT(IN)              :: AVelWN
    REAL(WP), INTENT(IN)                       :: dt
    INTEGER , INTENT(IN)                       :: spec_rank
    TYPE(REAL_DATA_LAYOUT), INTENT(IN)         :: AU, AV, AW
    type(COMPLEX_DATA_LAYOUT), INTENT(inout)   :: Anlx, Anly, Anlz
    REAL(WP),dimension(3)                      :: Na
    REAL(WP),dimension(3)                      :: Nb
    REAL(WP),dimension(3)                      :: kbrand
    REAL(WP)                                   :: phi

    IF (ialvelius.eq.1)  THEN
       if(ianis.EQ.0)then
          if(ihel.EQ.0)then
             CALL force_compute_alv_nohel(AUk,AVk,AWk,AVelWN, &
                                  & Anlx,Anly,Anlz,dt)
             if(inoEner.EQ.1)then
                CALL force_Helicity_no_Energy(AUk,AVk,AWk,AVelWN,Anlx,Anly,Anlz)
             endif
          ELSE
             CALL force_compute_alv_hel(AUk,AVk,AWk,AVelWN, &
                                  & Anlx,Anly,Anlz,dt)
          endif
       else
          if(ihel.EQ.0)then
             CALL force_compute_alv_nohel_anis(AUk,AVk,AWk,AVelWN, &
                                  & Anlx,Anly,Anlz,dt)
          ELSE
             CALL force_compute_alv_hel_anis(AUk,AVk,AWk,AVelWN, &
                                  & Anlx,Anly,Anlz,dt)
          endif
       endif
    ELSE IF (ibrandenburg.eq.1 .and. ihel.eq.1) THEN
       CALL force_brandenburg_init_hel(Na,Nb,kbrand,phi,spec_rank,dt)

       CALL force_brandenburg_hel(Na(1),Nb(1),kbrand,phi,Anlx,AU)
       CALL force_brandenburg_hel(Na(2),Nb(2),kbrand,phi,Anly,AU)
       CALL force_brandenburg_hel(Na(3),Nb(3),kbrand,phi,Anlz,AU)
    ELSE IF (ikolmo.eq.1) THEN
       CALL force_kolmo(Anlx,AU)
    ELSE IF (ilinear.eq.1) THEN
       CALL force_linear(Anlx,Anly,Anlz,AU,AV,AW)
    ELSE IF (isinus.eq.1) THEN
       CALL force_sinus(Anlx,Anly,Anlz,AU)
    ELSE IF (irobert.eq.1) THEN
!        if (spec_rank .EQ.0) print*, 'ça roule Robert!'
       CALL force_Robert_flow(Anlx, Anly, Anlz,AU)
    ELSE

       IF (spec_rank.EQ.0) THEN
             WRITE(6,'(a)')'PROBLEM IN force_compute NO FORCING TO CALL'
             STOP
       END IF
    ENDIF

END SUBROUTINE force_compute


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  Alvelius non-helicoidal forcing          !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> Initialize Alvelius forcing
!! @author Mouloud Kessar
!!    @param[in]  AUk         = Storage for the Spectral values of U
!!    @param[in]  AVk         = Storage for the Spectral values of V  
!!    @param[in]  AWk         = Storage for the Spectral values of W  
!!    @param[in]  AVelWN      = Velocity wavenumbers  
!!    @param[in]  mu          = viscosity
!!    @param[in]  spec_rank   = rank inside the spectral communicator
SUBROUTINE forcing_init_alv(AUk,AVk,AWk,AVelWN,mu,spec_rank)

  IMPLICIT NONE

  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN)   :: AUk, AVk, AWk
  TYPE(WaveNumbers), INTENT(IN)           :: AVelWN
  REAL(WP), INTENT(IN)                    :: mu
  INTEGER, INTENT(IN)                     :: spec_rank
  INTEGER                                 :: i,j,k
  REAL(WP) :: kk
  REAL(WP) :: force_sum
  REAL(WP) :: P
  REAL(WP) :: keta
  REAL(WP) :: eta
  REAL(WP) :: pi
  REAL(WP) :: Rlm
  REAL(WP) :: Pf
  INTEGER :: ik
  pi = acos(-1.0_WP)

  call parser_read('Kmax eta', keta)
  call parser_read('R Lambda',Rlm)
  call parser_read('Power factor',Pf)

  eta = keta*AUk%Lx / (pi*AUk%ny)
  P = mu**3/eta**4
  kf = (2.0_WP*pi)/((3.0_WP/20.0_WP*Rlm*Rlm)**0.75_WP*eta)
  kf = kf/Pf

  if (spec_rank.eq.0) print *,' Forcing wavenumber ', kf

  c = 0.5_WP
  ka = max(1e-10_WP,kf-2.0_WP) !
  kb = kf+1.0_WP

  force_sum = 0.0_WP
  DO k = AUk%zmin,AUk%zmax
     DO j =  AUk%ymin,AUk%ymax
        DO i =  AUk%xmin,AUk%xmax
           kk = sqrt(AVelWN%kx(i)*AVelWN%kx(i) + &
                   & AVelWN%ky(j)*AVelWN%ky(j) + &
                   & AVelWN%kz(k)*AVelWN%kz(k)   )
           ik = 1+idint(kk + 0.5_WP)
          if (kk.ge.ka.and.kk.le.kb) then
              force_sum = force_sum + force_spectrum(kk,kf,c)/(2.0_WP*pi*kk**2)
           ENDIF
        ENDDO
     ENDDO
  ENDDO

  A=DoGlobalSum(spec_rank,force_sum)

  A = P/A
  if (spec_rank.eq.0) print *,' Pre-multiplier for forcing term ', A
  
END SUBROUTINE forcing_init_alv
  

!> Compute Alvelius forcing (helcoidal )
!! @author Mouloud Kessar
!!    @param[in]  AUk           = Storage for the Spectral values of U
!!    @param[in]  AVk           = Storage for the Spectral values of V
!!    @param[in]  AWk           = Storage for the Spectral values of W
!!    @param[in]  AVelWN        = Velocity wavenumbers
!!    @param[in,out]  Anlx      = Nonlinear term in spectral space for U
!!    @param[in,out]  Anly      = Nonlinear term in spectral space for V
!!    @param[in,out]  Anlz      = Nonlinear term in spectral space for W
!!    @param[in]  dt            = time step
SUBROUTINE force_compute_alv_hel(AUk,AVk,AWk,AVelWN,Anlx,Anly,Anlz,dt)

  IMPLICIT NONE

  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN)      :: AUk, AVk, AWk
  TYPE(WaveNumbers), INTENT(IN)              :: AVelWN
  REAL(WP), INTENT(IN)                       :: dt
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT)   :: Anlx, Anly, Anlz
  INTEGER :: i,j,k
  REAL(WP) :: rand,pi,kk
  REAL(WP) :: phi3,psi,theta1,theta2
  REAL(WP) :: e1x,e1y,e1z
  REAL(WP) :: e2x,e2y,e2z
  REAL(WP) :: kk2,ga,gb
  COMPLEX(WP) :: xi1,xi2,af,bf
  REAL(WP) :: rxi1,ixi1,rxi2,ixi2
  REAL(WP) :: num,den,f
  COMPLEX(WP), PARAMETER :: ii=(0.0_WP,1.0_WP)
  pi = acos(-1.0_WP)

  DO k = Anlx%zmin,Anlx%zmax
    DO j =  Anlx%ymin,Anlx%ymax
      DO i =  Anlx%xmin,Anlx%xmax

           kk = sqrt(AVelWN%kx(i)*AVelWN%kx(i) + &
                   & AVelWN%ky(j)*AVelWN%ky(j) + &
                   & AVelWN%kz(k)*AVelWN%kz(k)    )

           kk2 = sqrt( AVelWN%ky(j)*AVelWN%ky(j) + &
                   & AVelWN%kx(i)*AVelWN%kx(i)    )

           if (kk2.gt.1e-10_WP.and.kk.gt.1e-10_WP) then
              if (kk.ge.ka.and.kk.le.kb) then

                 call random_number(rand)
                 phi3 = pi*rand
                 call random_number(rand)
                 psi = 2.0_WP*pi*rand

                 e1x = AVelWN%ky(j)/kk2
                 e1y = -AVelWN%kx(i)/kk2
                 e1z = 0.0_WP

                 e2x = AVelWN%kx(i)*AVelWN%kz(k)/(kk*kk2)
                 e2y = AVelWN%ky(j)*AVelWN%kz(k)/(kk*kk2)
                 e2z = -(kk2)/kk
                 xi1 = AUk%values(i,j,k)*e1x + AVk%values(i,j,k)*e1y + AWk%values(i,j,k)*e1z
                 xi2 = AUk%values(i,j,k)*e2x + AVk%values(i,j,k)*e2y + AWk%values(i,j,k)*e2z
                 rxi1 = real(0.5_WP*(xi1 + conjg(xi1)),WP)
                 rxi2 = real(0.5_WP*(xi2 + conjg(xi2)),WP)
                 ixi1 = real(0.5_WP*(-ii)*(xi1 - conjg(xi1)),WP)
                 ixi2 = real(0.5_WP*(-ii)*(xi2 - conjg(xi2)),WP)
                 ga = sin(2.0_WP*phi3)
                 gb = cos(2.0_WP*phi3)

                    num =  (ga+gb*sin(psi))*rxi1 - gb*cos(psi)*rxi2 + (gb*cos(psi))*ixi1 + (ga+gb*sin(psi))*ixi2
                    den = - gb*cos(psi)*rxi1 + (ga + gb*sin(psi) )*ixi1 - gb*cos(psi)*rxi2 + ( ga + gb*sin(psi))*ixi2

                 IF (num .eq. 0.0_WP) THEN
                    theta1 = 0.0_WP
                 ELSE
                    theta1 = atan(num/den)
                 ENDIF
                 theta2 = psi + theta1

                 f = (A/dt)*force_spectrum(kk,kf,c)/(2.0_WP*pi*kk**2)
                 af = (f**0.5_WP)*exp(ii*theta1)*ga
                 bf = (f**0.5_WP)*exp(ii*theta2)*gb

                    Anlx%values(i,j,k) = Anlx%values(i,j,k) + ( (af+ii*bf)*e1x + ii*(af+ii*bf)*e2x )/(sqrt(2.0))
                    Anly%values(i,j,k) = Anly%values(i,j,k) + ( (af+ii*bf)*e1y + ii*(af+ii*bf)*e2y )/(sqrt(2.0))
                    Anlz%values(i,j,k) = Anlz%values(i,j,k) + ( ii*(af+ii*bf)*e2z )/(sqrt(2.0))

              ENDIF
           ENDIF
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE force_compute_alv_hel

!> Compute Alvelius forcing (no helicoidal)
!! @author Mouloud Kessar
!!    @param[in]  AUk           = Storage for the Spectral values of U
!!    @param[in]  AVk           = Storage for the Spectral values of V
!!    @param[in]  AWk           = Storage for the Spectral values of W
!!    @param[in]  AVelWN        = Velocity wavenumbers
!!    @param[in,out]  Anlx      = Nonlinear term in spectral space for U
!!    @param[in,out]  Anly      = Nonlinear term in spectral space for V
!!    @param[in,out]  Anlz      = Nonlinear term in spectral space for W
!!    @param[in]  dt            = time step
SUBROUTINE force_compute_alv_nohel(AUk,AVk,AWk,AVelWN,Anlx,Anly,Anlz,dt)

  IMPLICIT NONE

  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN)      :: AUk, AVk, AWk
  TYPE(WaveNumbers), INTENT(IN)              :: AVelWN
  REAL(WP), INTENT(IN)                       :: dt
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT)   :: Anlx, Anly, Anlz
  INTEGER :: i,j,k
  REAL(WP) :: rand,pi,kk
  REAL(WP) :: phi3,psi,theta1,theta2
  REAL(WP) :: e1x,e1y,e1z
  REAL(WP) :: e2x,e2y,e2z
  REAL(WP) :: kk2,ga,gb
  COMPLEX(WP) :: xi1,xi2,af,bf
  REAL(WP) :: rxi1,ixi1,rxi2,ixi2
  REAL(WP) :: num,den,f
  COMPLEX(WP), PARAMETER :: ii=(0.0_WP,1.0_WP)
  pi = acos(-1.0_WP)

  DO k = Anlx%zmin,Anlx%zmax
    DO j =  Anlx%ymin,Anlx%ymax
      DO i =  Anlx%xmin,Anlx%xmax

           kk = sqrt(AVelWN%kx(i)*AVelWN%kx(i) + &
                   & AVelWN%ky(j)*AVelWN%ky(j) + &
                   & AVelWN%kz(k)*AVelWN%kz(k)    )

           kk2 = sqrt( AVelWN%ky(j)*AVelWN%ky(j) + &
                   & AVelWN%kx(i)*AVelWN%kx(i)    )

           if (kk2.gt.1e-10_WP.and.kk.gt.1e-10_WP) then
              if (kk.ge.ka.and.kk.le.kb) then

                 call random_number(rand)
                 phi3 = pi*rand
                 call random_number(rand)
                 psi = 2.0_WP*pi*rand

                 e1x = AVelWN%ky(j)/kk2
                 e1y = -AVelWN%kx(i)/kk2
                 e1z = 0.0_WP

                 e2x = AVelWN%kx(i)*AVelWN%kz(k)/(kk*kk2)
                 e2y = AVelWN%ky(j)*AVelWN%kz(k)/(kk*kk2)
                 e2z = -(kk2)/kk
                 xi1 = AUk%values(i,j,k)*e1x + AVk%values(i,j,k)*e1y + AWk%values(i,j,k)*e1z
                 xi2 = AUk%values(i,j,k)*e2x + AVk%values(i,j,k)*e2y + AWk%values(i,j,k)*e2z
                 rxi1 = real(0.5_WP*(xi1 + conjg(xi1)),WP)
                 rxi2 = real(0.5_WP*(xi2 + conjg(xi2)),WP)
                 ixi1 = real(0.5_WP*(-ii)*(xi1 - conjg(xi1)),WP)
                 ixi2 = real(0.5_WP*(-ii)*(xi2 - conjg(xi2)),WP)
                 ga = sin(2.0_WP*phi3)
                 gb = cos(2.0_WP*phi3)

                    num =  ga*rxi1 + gb*(sin(psi)*ixi2 + cos(psi)*rxi2)
                    den = -ga*ixi1 + gb*(sin(psi)*rxi2 - cos(psi)*ixi2)
           
                 IF (num .eq. 0.0_WP) THEN
                    theta1 = 0.0_WP
                 ELSE
                    theta1 = atan(num/den)
                 ENDIF
                 theta2 = psi + theta1

                 f = (A/dt)*force_spectrum(kk,kf,c)/(2.0_WP*pi*kk**2)
                 af = (f**0.5_WP)*exp(ii*theta1)*ga
                 bf = (f**0.5_WP)*exp(ii*theta2)*gb

                    Anlx%values(i,j,k) = Anlx%values(i,j,k) + (af*e1x+bf*e2x)
                    Anly%values(i,j,k) = Anly%values(i,j,k) + (af*e1y+bf*e2y)
                    Anlz%values(i,j,k) = Anlz%values(i,j,k) + (bf*e2z)
 
              ENDIF
           ENDIF
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE force_compute_alv_nohel

!> Compute Alvelius forcing (no helicoidal)
!! @author Mouloud Kessar
!!    @param[in]  AUk           = Storage for the Spectral values of U
!!    @param[in]  AVk           = Storage for the Spectral values of V
!!    @param[in]  AWk           = Storage for the Spectral values of W
!!    @param[in]  AVelWN        = Velocity wavenumbers
!!    @param[in,out]  Anlx      = Nonlinear term in spectral space for U
!!    @param[in,out]  Anly      = Nonlinear term in spectral space for V
!!    @param[in,out]  Anlz      = Nonlinear term in spectral space for W
!!    @param[in]  dt            = time step
SUBROUTINE force_compute_alv_nohel_anis(AUk,AVk,AWk,AVelWN,Anlx,Anly,Anlz,dt)

  IMPLICIT NONE

  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN)      :: AUk, AVk, AWk
  TYPE(WaveNumbers), INTENT(IN)              :: AVelWN
  REAL(WP), INTENT(IN)                       :: dt
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT)   :: Anlx, Anly, Anlz
  INTEGER :: i,j,k
  REAL(WP) :: rand,pi,kk
  REAL(WP) :: phi3,psi,theta1,theta2
  REAL(WP) :: e1x,e1y,e1z
  REAL(WP) :: e2x,e2y,e2z
  REAL(WP) :: kk2,ga,gb
  COMPLEX(WP) :: xi1,xi2,af,bf
  REAL(WP) :: rxi1,ixi1,rxi2,ixi2
  REAL(WP) :: num,den,f
  COMPLEX(WP), PARAMETER :: ii=(0.0_WP,1.0_WP)
  pi = acos(-1.0_WP)

  DO k = Anlx%zmin,Anlx%zmax
    DO j =  Anlx%ymin,Anlx%ymax
      DO i =  Anlx%xmin,Anlx%xmax

           kk = sqrt(AVelWN%kx(i)*AVelWN%kx(i) + &
                   & AVelWN%ky(j)*AVelWN%ky(j) + &
                   & AVelWN%kz(k)*AVelWN%kz(k)    )

           kk2 = sqrt( AVelWN%ky(j)*AVelWN%ky(j) + &
                   & AVelWN%kx(i)*AVelWN%kx(i)    )

           if (kk2.gt.1e-10_WP.and.kk.gt.1e-10_WP) then
              if (kk.ge.ka.and.kk.le.kb) then

                 call random_number(rand)
                 phi3 = pi*rand
                 call random_number(rand)
                 psi = 2.0_WP*pi*rand

                 e1x = AVelWN%ky(j)/kk2
                 e1y = -AVelWN%kx(i)/kk2
                 e1z = 0.0_WP

                 e2x = AVelWN%kx(i)*AVelWN%kz(k)/(kk*kk2)
                 e2y = AVelWN%ky(j)*AVelWN%kz(k)/(kk*kk2)
                 e2z = -(kk2)/kk
                 xi1 = AUk%values(i,j,k)*e1x + AVk%values(i,j,k)*e1y + AWk%values(i,j,k)*e1z
                 xi2 = AUk%values(i,j,k)*e2x + AVk%values(i,j,k)*e2y + AWk%values(i,j,k)*e2z
                 rxi1 = real(0.5_WP*(xi1 + conjg(xi1)),WP)
                 rxi2 = real(0.5_WP*(xi2 + conjg(xi2)),WP)
                 ixi1 = real(0.5_WP*(-ii)*(xi1 - conjg(xi1)),WP)
                 ixi2 = real(0.5_WP*(-ii)*(xi2 - conjg(xi2)),WP)
                 ga = tanh(b_anis*(phi3-pi/2))/tanh(b_anis*(pi/2))
                 gb = sqrt(1-ga*ga)

                    num =  ga*rxi1 + gb*(sin(psi)*ixi2 + cos(psi)*rxi2)
                    den = -ga*ixi1 + gb*(sin(psi)*rxi2 - cos(psi)*ixi2)
           
                 IF (num .eq. 0.0_WP) THEN
                    theta1 = 0.0_WP
                 ELSE
                    theta1 = atan(num/den)
                 ENDIF
                 theta2 = psi + theta1

                 f = (A/dt)*force_spectrum(kk,kf,c)/(2.0_WP*pi*kk**2)
                 af = (f**0.5_WP)*exp(ii*theta1)*ga
                 bf = (f**0.5_WP)*exp(ii*theta2)*gb

                    Anlx%values(i,j,k) = Anlx%values(i,j,k) + (af*e1x+bf*e2x)
                    Anly%values(i,j,k) = Anly%values(i,j,k) + (af*e1y+bf*e2y)
                    Anlz%values(i,j,k) = Anlz%values(i,j,k) + (bf*e2z)
 
              ENDIF
           ENDIF
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE force_compute_alv_nohel_anis

!> Compute Alvelius forcing (helcoidal )
!! @author Mouloud Kessar
!!    @param[in]  AUk           = Storage for the Spectral values of U
!!    @param[in]  AVk           = Storage for the Spectral values of V
!!    @param[in]  AWk           = Storage for the Spectral values of W
!!    @param[in]  AVelWN        = Velocity wavenumbers
!!    @param[in,out]  Anlx      = Nonlinear term in spectral space for U
!!    @param[in,out]  Anly      = Nonlinear term in spectral space for V
!!    @param[in,out]  Anlz      = Nonlinear term in spectral space for W
!!    @param[in]  dt            = time step
SUBROUTINE force_compute_alv_hel_anis(AUk,AVk,AWk,AVelWN,Anlx,Anly,Anlz,dt)

  IMPLICIT NONE

  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN)      :: AUk, AVk, AWk
  TYPE(WaveNumbers), INTENT(IN)              :: AVelWN
  REAL(WP), INTENT(IN)                       :: dt
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT)   :: Anlx, Anly, Anlz
  INTEGER :: i,j,k
  REAL(WP) :: rand,pi,kk
  REAL(WP) :: phi3,psi,theta1,theta2
  REAL(WP) :: e1x,e1y,e1z
  REAL(WP) :: e2x,e2y,e2z
  REAL(WP) :: kk2,ga,gb
  COMPLEX(WP) :: xi1,xi2,af,bf
  REAL(WP) :: rxi1,ixi1,rxi2,ixi2
  REAL(WP) :: num,den,f
  COMPLEX(WP), PARAMETER :: ii=(0.0_WP,1.0_WP)
  pi = acos(-1.0_WP)

  DO k = Anlx%zmin,Anlx%zmax
    DO j =  Anlx%ymin,Anlx%ymax
      DO i =  Anlx%xmin,Anlx%xmax

           kk = sqrt(AVelWN%kx(i)*AVelWN%kx(i) + &
                   & AVelWN%ky(j)*AVelWN%ky(j) + &
                   & AVelWN%kz(k)*AVelWN%kz(k)    )

           kk2 = sqrt( AVelWN%ky(j)*AVelWN%ky(j) + &
                   & AVelWN%kx(i)*AVelWN%kx(i)    )

           if (kk2.gt.1e-10_WP.and.kk.gt.1e-10_WP) then
              if (kk.ge.ka.and.kk.le.kb) then

                 call random_number(rand)
                 phi3 = pi*rand
                 call random_number(rand)
                 psi = 2.0_WP*pi*rand

                 e1x = AVelWN%ky(j)/kk2
                 e1y = -AVelWN%kx(i)/kk2
                 e1z = 0.0_WP

                 e2x = AVelWN%kx(i)*AVelWN%kz(k)/(kk*kk2)
                 e2y = AVelWN%ky(j)*AVelWN%kz(k)/(kk*kk2)
                 e2z = -(kk2)/kk
                 xi1 = AUk%values(i,j,k)*e1x + AVk%values(i,j,k)*e1y + AWk%values(i,j,k)*e1z
                 xi2 = AUk%values(i,j,k)*e2x + AVk%values(i,j,k)*e2y + AWk%values(i,j,k)*e2z
                 rxi1 = real(0.5_WP*(xi1 + conjg(xi1)),WP)
                 rxi2 = real(0.5_WP*(xi2 + conjg(xi2)),WP)
                 ixi1 = real(0.5_WP*(-ii)*(xi1 - conjg(xi1)),WP)
                 ixi2 = real(0.5_WP*(-ii)*(xi2 - conjg(xi2)),WP)
                 ga = tanh(b_anis*(phi3-pi/2))/tanh(b_anis*(pi/2))
                 gb = sqrt(1-ga*ga)

                    num =  (ga+gb*sin(psi))*rxi1 - gb*cos(psi)*rxi2 + (gb*cos(psi))*ixi1 + (ga+gb*sin(psi))*ixi2
                    den = - gb*cos(psi)*rxi1 + (ga + gb*sin(psi) )*ixi1 - gb*cos(psi)*rxi2 + ( ga + gb*sin(psi))*ixi2

                 IF (num .eq. 0.0_WP) THEN
                    theta1 = 0.0_WP
                 ELSE
                    theta1 = atan(num/den)
                 ENDIF
                 theta2 = psi + theta1

                 f = (A/dt)*force_spectrum(kk,kf,c)/(2.0_WP*pi*kk**2)
                 af = (f**0.5_WP)*exp(ii*theta1)*ga
                 bf = (f**0.5_WP)*exp(ii*theta2)*gb

                    Anlx%values(i,j,k) = Anlx%values(i,j,k) + ( (af+ii*bf)*e1x + ii*(af+ii*bf)*e2x )/(sqrt(2.0))
                    Anly%values(i,j,k) = Anly%values(i,j,k) + ( (af+ii*bf)*e1y + ii*(af+ii*bf)*e2y )/(sqrt(2.0))
                    Anlz%values(i,j,k) = Anlz%values(i,j,k) + ( ii*(af+ii*bf)*e2z )/(sqrt(2.0))

              ENDIF
           ENDIF
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE force_compute_alv_hel_anis

FUNCTION force_spectrum(kk,kf,c)

  USE precision_tools
  IMPLICIT NONE
  
  REAL(WP) :: force_spectrum
  REAL(WP) :: kk,kf,c
 
  force_spectrum = exp(-(kk-kf)**2/c)

END FUNCTION force_spectrum


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Brandenburg helicoidal forcing related subroutines  !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> Initialization of Brandenburg forcing at each time step 
!! @author Mouloud Kessar
!!    @param[in,out]  Na            = ??????????????
!!    @param[in,out]  Nb            = ??????????????
!!    @param[in,out]  kbrand        = random real vector 
!!    @param[in,out]  phi           = random phase
!!    @param[in]  spec_rank         = rank inside the spectral communicator
!!    @param[in]  dt                = time step
SUBROUTINE force_brandenburg_init_hel(Na,Nb,kbrand,phi,spec_rank,dt)

  IMPLICIT NONE

  REAL(WP), INTENT(IN)                   :: dt
  INTEGER, INTENT(IN)                    :: spec_rank
  rEAL(WP), INTENT(INOUT)                :: phi
  REAL(WP), DIMENSION(3), INTENT(INOUT)  :: Na, Nb, kbrand

  REAL(WP)                               :: rand,pi,phik,thetak,ex,ey,ez, knn, k2
  REAL(WP)                               :: pvect1,pvect2,pvect3,kpvect1,kpvect2,kpvect3
  REAL(WP)                               :: Na0,Nb0
  INTEGER                                :: ierr

! Pi number
  pi = acos(-1.0_WP)


! arbitrary unit vector
  ex=0.0_WP;
  ey=1.0_WP;
  ez=0.0_WP;

  IF (spec_rank.eq.0) THEN
!    random phase
     CALL random_number(rand)
     phi=pi*(rand-0.5)*2

     CALL random_number(rand)
     phik=pi*(rand-0.5)*2
     CALL random_number(rand)
     thetak=pi*rand

     kbrand(1)=kn*cos(thetak)*sin(phik)
     kbrand(2)=kn*sin(thetak)*sin(phik)
     kbrand(3)=kn*cos(phik)

  END IF


  CALL MPI_BCAST(kbrand,3,MPI_REAL_WP,0,spec_communicator,ierr)
  CALL MPI_BCAST(phi,1,MPI_REAL_WP,0,spec_communicator,ierr)


 ! Compute value constant for the iteration
  knn = sqrt(kbrand(1)**2.0 + kbrand(2)**2.0 + kbrand(3)**2.0)
  k2 = knn*knn

  pvect1=kbrand(2)*ez-kbrand(3)*ey
  pvect2=kbrand(3)*ex-kbrand(1)*ez
  pvect3=kbrand(1)*ey-kbrand(2)*ex

  kpvect1=kbrand(2)*pvect3-kbrand(3)*pvect2
  kpvect2=kbrand(3)*pvect1-kbrand(1)*pvect3
  kpvect3=kbrand(1)*pvect2-kbrand(2)*pvect1

  Na0=   f0*(sqrt(knn/dt))     /   (    &
     &                               2.0*k2*sqrt( &
     &                                            1.0_WP   -    (    &
     &                                                             (   kbrand(1)*ex+kbrand(2)*ey+kbrand(3)*ez    )**2.0   &
     &                                                          )/k2   &
     &                                          ) &
     &                              )
  Nb0=Na0*knn


  Na(1)=Na0*kpvect1
  Na(2)=Na0*kpvect2
  Na(3)=Na0*kpvect3

  Nb(1)=Nb0*pvect1
  Nb(2)=Nb0*pvect2
  Nb(3)=Nb0*pvect3

END SUBROUTINE force_brandenburg_init_hel

!> Computation of forcing term in direction i 
!> Version 2 By P. Begou
!! @author Mouloud Kessar
!!    @param[in]      Na1           = Na component in direction i
!!    @param[in]      Nb1           = Nb component in direction i
!!    @param[in]      kbrand        = random real vector 
!!    @param[in]      phi           = random phase
!!    @param[in,out]  nik           = non linear term in direction i
!!    @param[in]      Ui            = x-velocity in physical space [only used to copy structure, see with Patrick if better is possible]
SUBROUTINE force_brandenburg_hel(Na1,Nb1,kbrand,phi,nik,Ui)

  IMPLICIT NONE

  REAL(WP), INTENT(IN)                       :: Na1, Nb1, phi
  REAL(WP), DIMENSION(3), INTENT(IN)         :: kbrand
  TYPE(REAL_DATA_LAYOUT), INTENT(IN)         :: Ui
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT)   :: nik

  REAL(WP)                                   :: theta, x, y, z
  REAL(WP)                                   :: onx, ony, onz 
  TYPE(REAL_DATA_LAYOUT)                     :: fx
  TYPE(COMPLEX_DATA_LAYOUT)                  :: niktemp
  INTEGER                                    :: i,j,k
  LOGICAL                                    :: res

  IF (.NOT. copyStructOnly(nik,niktemp)) RETURN
  IF (.NOT. copyStructOnly(Ui,fx)) RETURN

! set defaul array value to 0.0
  fx%values= 0.0_WP

! just calculate this out of the loop please...
  onx = fx%lx/fx%nx
  ony = fx%ly/fx%ny
  onz = fx%lz/fx%nz 

  DO k = fx%zmin,fx%zmax
     z = onz*(k-1)
     DO j = fx%ymin,fx%ymax
        y = ony*(j-1)
        DO i = fx%xmin,fx%xmax
           x = onx*(i-1)
           res=((x.GE.lxfmin).AND.(x.LE.lxfmax) .AND. &
              & (y.GE.lyfmin).AND.(y.LE.lyfmax) .AND. &
              & (z.GE.lzfmin).AND.(z.LE.lzfmax))
           IF (res) THEN
              theta = kbrand(1)*x + kbrand(2)*y + kbrand(3)*z + phi
              fx%values(i,j,k) = Na1*cos(theta) + Nb1*sin(theta)
           ENDIF
        END DO
     END DO
  END DO

  CALL ftran(fx,niktemp,res)
  

  nik%values = nik%values+niktemp%values

  CALL deleteDataLayout(fx)
  CALL deleteDataLayout(niktemp)
END SUBROUTINE force_brandenburg_hel
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!   Kolmogorv forcing related subroutines   !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE force_kolmo(nik,Ui)

    IMPLICIT NONE

    TYPE(REAL_DATA_LAYOUT), INTENT(IN)         :: Ui
    TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT)   :: nik

    INTEGER  :: i,j,k
    REAL(WP) :: yy
    TYPE(REAL_DATA_LAYOUT) :: nis
    TYPE(COMPLEX_DATA_LAYOUT) :: niktemp
    LOGICAL     :: res

    IF (.NOT. copyStructOnly(nik,niktemp)) RETURN
    IF (.NOT. copyStructOnly(Ui,nis)) RETURN

    do k = nis%zmin,nis%zmax
       yy = nis%Lz / real(nis%nz)*real(k-1)
       if (mode_k.ne.0) then
          nis%values(:,:,k) = - (sin(yy) + sin(real(mode_k)*yy))  !! Following Rollin et al., JFM, 2011
       else
          nis%values(:,:,k) = - sin(yy)                           !! Following Sarris et al., Phys. Fluids, 2007
       end if
    end do

    CALL ftran(nis,niktemp,res)
    DO k = nik%zmin,nik%zmax
       DO j = nik%ymin,nik%ymax
          DO i = nik%xmin,nik%xmax
             nik%values(i,j,k) = nik%values(i,j,k)+niktemp%values(i,j,k)
          END DO
       END DO
    END DO
    CALL deleteDataLayout(nis)
    CALL deleteDataLayout(niktemp)

END SUBROUTINE force_kolmo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!   Linear forcing (Rosales and Meneveau, PoF, 2005   !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE force_linear(Anlx,Anly,Anlz,AU,AV,AW)

    IMPLICIT NONE

    TYPE(REAL_DATA_LAYOUT), INTENT(IN)         :: AU,AV,AW
    TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT)   :: Anlx, Anly, Anlz

    TYPE(REAL_DATA_LAYOUT) :: compx, compy, compz
    TYPE(COMPLEX_DATA_LAYOUT) :: nkx, nky, nkz
    INTEGER :: i,j,k
    REAL :: kk
    LOGICAL     :: res

    IF (.NOT. copyStructOnly(Anlx,nkx)) RETURN
    IF (.NOT. copyStructOnly(Anlx,nky)) RETURN
    IF (.NOT. copyStructOnly(Anlx,nkz)) RETURN
    IF (.NOT. copyStructOnly(AU,compx)) RETURN
    IF (.NOT. copyStructOnly(AU,compy)) RETURN
    IF (.NOT. copyStructOnly(AU,compz)) RETURN

    compx%values = f0*AU%values
    compy%values = f0*AV%values
    compz%values = f0*AW%values

    CALL ftran(compx,nkx,res)
    CALL ftran(compy,nky,res)
    CALL ftran(compz,nkz,res)

    Anlx%values = Anlx%values + nkx%values
    Anly%values = Anly%values + nky%values
    Anlz%values = Anlz%values + nkz%values

    !!!! To avoid to amplify the error on the first mode
    DO k = Anlx%zmin,Anlx%zmax
        DO j = Anlx%ymin,Anlx%ymax
            DO i = Anlx%xmin,Anlx%xmax
                if (i.eq.1.and.j.eq.1.and.k.eq.1) then
                    Anlx%values(i,j,k) = 0.0_WP
                    Anly%values(i,j,k) = 0.0_WP
                    Anlz%values(i,j,k) = 0.0_WP
                end if
            ENDDO
        ENDDO
    ENDDO

    CALL deleteDataLayout(compx)
    CALL deleteDataLayout(compy)
    CALL deleteDataLayout(compz)
    CALL deleteDataLayout(nkx)
    CALL deleteDataLayout(nky)
    CALL deleteDataLayout(nkz)

END SUBROUTINE force_linear

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!   Sinus forcing (Andel Kareem et al., CMA, 2009)    !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE force_sinus(Anlx,Anly,Anlz,AU)

    IMPLICIT NONE

    TYPE(REAL_DATA_LAYOUT), INTENT(IN)         :: AU
    TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT)   :: Anlx, Anly, Anlz

    TYPE(REAL_DATA_LAYOUT) :: compx, compy, compz
    TYPE(COMPLEX_DATA_LAYOUT) :: nkx, nky, nkz
    INTEGER :: i,j,k
    REAL :: kk
    LOGICAL     :: res

    REAL(WP):: dx, dy, dz
    REAL(WP):: kx, ky, kz, k2
    REAL(WP):: yj,zk,xi
    REAL(WP):: twopi, pi, Gsin

    twopi=2.0_WP*acos(-1.0)
    pi=acos(-1.0)

    kx = twopi/AU%Lx
    ky = twopi/AU%Ly
    kz = twopi/AU%Lz
    k2 = kx**2.0 + ky**2.0 + kz**2.0

    IF (.NOT. copyStructOnly(Anlx,nkx)) RETURN
    IF (.NOT. copyStructOnly(Anlx,nky)) RETURN
    IF (.NOT. copyStructOnly(Anlx,nkz)) RETURN
    IF (.NOT. copyStructOnly(AU,compx)) RETURN
    IF (.NOT. copyStructOnly(AU,compy)) RETURN
    IF (.NOT. copyStructOnly(AU,compz)) RETURN

    dx = AU%Lx/real(AU%nx,WP)
    dy = AU%Ly/real(AU%ny,WP)
    dz = AU%Lz/real(AU%nz,WP)

    DO k = AU%zmin,AU%zmax
       zk = real(k-1,WP)*dz
       DO j = AU%ymin,AU%ymax
          yj = real(j-1,WP)*dy
          DO i = AU%xmin,AU%xmax
             xi = real(i-1,WP)*dx

!             Gsin = sin(xi + yj + zk)
!             compx%values(i,j,k) = 2.0 * f0 * (ky*kz)/k2 * Gsin
!             compy%values(i,j,k) =     - f0 * (kx*kz)/k2 * Gsin 
!             compz%values(i,j,k) =     - f0 * (kx*ky)/k2 * Gsin
!           --> other form because this form does not give the same <u**2.0> for the 3 direction
             compx%values(i,j,k) = f0 * sin(yj-pi)*sin(zk) 
             compy%values(i,j,k) = f0 * sin(xi-pi)*sin(zk)
             compz%values(i,j,k) = f0 * sin(xi-pi)*sin(yj-pi)

            ENDDO
        ENDDO
    ENDDO

    CALL ftran(compx,nkx,res)
    CALL ftran(compy,nky,res)
    CALL ftran(compz,nkz,res)

    Anlx%values = Anlx%values + nkx%values
    Anly%values = Anly%values + nky%values
    Anlz%values = Anlz%values + nkz%values

    !!!! To avoid to amplify the error on the first mode
!    DO k = Anlx%zmin,Anlx%zmax
!        DO j = Anlx%ymin,Anlx%ymax
!            DO i = Anlx%xmin,Anlx%xmax
!                if (i.eq.1.and.j.eq.1.and.k.eq.1) then
!                    Anlx%values(i,j,k) = 0.0_WP
!                    Anly%values(i,j,k) = 0.0_WP
!                    Anlz%values(i,j,k) = 0.0_WP
!                end if
!            ENDDO
!        ENDDO
!    ENDDO

    CALL deleteDataLayout(compx)
    CALL deleteDataLayout(compy)
    CALL deleteDataLayout(compz)
    CALL deleteDataLayout(nkx)
    CALL deleteDataLayout(nky)
    CALL deleteDataLayout(nkz)

END SUBROUTINE force_sinus

!> Add gravity in direction i 
!! @author Guillaume Balarac 
!!    @param[in,out]  nik           = non linear term in direction i in spectral space
!!    @param[in]      grav          = value of the gravity
!!    @param[in]      Ui            = x-velocity in physical space [only used to copy structure, see with Patrick if better is possible]
SUBROUTINE gravity(nik,grav,Ui)

  IMPLICIT NONE

  TYPE(REAL_DATA_LAYOUT), INTENT(IN)         :: Ui
  REAL(WP), INTENT(IN)                       :: grav 
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT)   :: nik

  TYPE(REAL_DATA_LAYOUT)                     :: nis
  TYPE(COMPLEX_DATA_LAYOUT)                  :: niktemp
  INTEGER                                    :: i,j,k
  LOGICAL                                    :: res

  IF (.NOT. copyStructOnly(nik,niktemp)) RETURN
  IF (.NOT. copyStructOnly(Ui,nis)) RETURN
  DO k = nis%zmin,nis%zmax
     DO j = nis%ymin,nis%ymax
        DO i = nis%xmin,nis%xmax
           nis%values(i,j,k)= grav 
        END DO
     END DO
  END DO
  CALL ftran(nis,niktemp,res)
  DO k = nik%zmin,nik%zmax
     DO j = nik%ymin,nik%ymax
        DO i = nik%xmin,nik%xmax
           nik%values(i,j,k) = nik%values(i,j,k)+niktemp%values(i,j,k)
        END DO
     END DO
  END DO
  CALL deleteDataLayout(nis)
  CALL deleteDataLayout(niktemp)

END SUBROUTINE gravity 

SUBROUTINE force_Robert_flow(nikx,niky,nikz,Ui)

    IMPLICIT NONE

    TYPE(REAL_DATA_LAYOUT), INTENT(IN)         :: Ui
    TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT)   :: nikx,niky,nikz

    INTEGER  :: i,j,k
    REAL(WP) :: yy,xx,zz
    TYPE(REAL_DATA_LAYOUT) :: nis
    TYPE(COMPLEX_DATA_LAYOUT) :: niktemp
    LOGICAL     :: res

    IF (.NOT. copyStructOnly(nikx,niktemp)) RETURN
    IF (.NOT. copyStructOnly(Ui,nis)) RETURN
!!!!!Fz
    do i = nis%xmin,nis%xmax
       xx = nis%Lx / real(nis%nx)*real(i-1)
       do j = nis%ymin,nis%ymax
          yy = nis%Ly / real(nis%ny)*real(j-1)
!           if (mode_k.ne.0) then
             nis%values(i,j,:) = sqrt(2.0_WP)*sin(xx)*sin(yy)                           
!           end if
       end do
    end do

    CALL ftran(nis,niktemp,res)
    DO k = nikz%zmin,nikz%zmax
       DO j = nikz%ymin,nikz%ymax
          DO i = nikz%xmin,nikz%xmax
             nikz%values(i,j,k) = nikz%values(i,j,k)+niktemp%values(i,j,k)
          END DO
       END DO
    END DO
!!!!!!Fy
    DO i = nis%xmin,nis%xmax
       xx = nis%Lx / real(nis%nx)*real(i-1)
       DO j = nis%ymin,nis%ymax
          yy = nis%Ly / real(nis%ny)*real(j-1)
!           if (mode_k.ne.0) then
             nis%values(i,j,:) = -cos(xx)*sin(yy)                           
!           end if
       END DO
    END DO

    CALL ftran(nis,niktemp,res)
    DO k = niky%zmin,niky%zmax
       DO j = niky%ymin,niky%ymax
          DO i = niky%xmin,niky%xmax
             niky%values(i,j,k) = niky%values(i,j,k)+niktemp%values(i,j,k)
          END DO
       END DO
    END DO
!!!!!!!Fx
    DO i = nis%xmin,nis%xmax
       xx = nis%Lx / real(nis%nx)*real(i-1)
       DO j = nis%ymin,nis%ymax
          yy = nis%Ly / real(nis%ny)*real(j-1)
!           if (mode_k.ne.0) then
             nis%values(i,j,:) = cos(yy)*sin(xx)                           
!           end if
       END DO
    END DO

    CALL ftran(nis,niktemp,res)
    DO k = nikx%zmin,nikx%zmax
       DO j = nikx%ymin,nikx%ymax
          DO i = nikx%xmin,nikx%xmax
             nikx%values(i,j,k) = nikx%values(i,j,k)+niktemp%values(i,j,k)
          END DO
       END DO
    END DO

    CALL deleteDataLayout(nis)
    CALL deleteDataLayout(niktemp)

END SUBROUTINE force_Robert_flow


SUBROUTINE force_Helicity_no_Energy(AUk,AVk,AWk,AVelWN,Anlx,Anly,Anlz)

  IMPLICIT NONE

  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN)      :: AUk, AVk, AWk
  TYPE(WaveNumbers), INTENT(IN)              :: AVelWN
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT)   :: Anlx, Anly, Anlz
  INTEGER :: i,j,k
  REAL(WP) :: rand,pi,kk
  REAL(WP) :: phi3,psi,theta1,theta2
  REAL(WP) :: eAx,eAy,eAz
  REAL(WP) :: eBx,eBy,eBz
  REAL(WP) :: kk2
  COMPLEX(WP) :: xiA,xiB,af,bf
  REAL(WP) :: A1,A2,UkrB,UkiA !rxiA,ixiA,rxiB,ixiB
  REAL(WP) :: num,den,f
  COMPLEX(WP), PARAMETER :: ii=(0.0_WP,1.0_WP)
  REAL(WP) ::phik,thetak
  REAL(WP) :: SommeEnerFluctuation,SommeSelfCorForcing
  SommeEnerFluctuation =0.0_WP
  SommeSelfCorForcing =0.0_WP

  pi = acos(-1.0_WP)

  DO k = Anlx%zmin,Anlx%zmax
    DO j =  Anlx%ymin,Anlx%ymax
      DO i =  Anlx%xmin,Anlx%xmax

           kk = sqrt(AVelWN%kx(i)*AVelWN%kx(i) + &
                   & AVelWN%ky(j)*AVelWN%ky(j) + &
                   & AVelWN%kz(k)*AVelWN%kz(k)    )

           kk2 = sqrt( AVelWN%ky(j)*AVelWN%ky(j) + &
                   & AVelWN%kx(i)*AVelWN%kx(i)    )

           if (kk2.gt.1e-10_WP.and.kk.gt.1e-10_WP) then
              if (kk.ge.ka.and.kk.le.kb) then


                 CALL random_number(rand)
                 phik=pi*(rand-0.5)*2
                 CALL random_number(rand)
                 thetak=pi*rand
 
                 eAx = cos(thetak)*sin(phik)
                 eAy = sin(thetak)*sin(phik)
                 eAz = cos(phik)


                 CALL random_number(rand)
                 phik=pi*(rand-0.5)*2
                 CALL random_number(rand)
                 thetak=pi*rand

                 eBx = cos(thetak)*sin(phik)
                 eBy = sin(thetak)*sin(phik)
                 eBz = cos(phik)

   !              xiA = AUk%values(i,j,k)*eAx + AVk%values(i,j,k)*eAy + AWk%values(i,j,k)*eAz
   !              xiB = AUk%values(i,j,k)*eBx + AVk%values(i,j,k)*eBy + AWk%values(i,j,k)*eBz


                 !A1=K.(UkR times eA)
                 A1 = AVelWN%kx(i)*((real(0.5_WP*(AVk%values(i,j,k) + conjg(AVk%values(i,j,k))),WP)*eAz) &
                                 & -(real(0.5_WP*(AWk%values(i,j,k) + conjg(AWk%values(i,j,k))),WP)*eAy))&
                  & + AVelWN%ky(j)*((real(0.5_WP*(AWk%values(i,j,k) + conjg(AWk%values(i,j,k))),WP)*eAx) &
                                 & -(real(0.5_WP*(AUk%values(i,j,k) + conjg(AUk%values(i,j,k))),WP)*eAz))&
                  & + AVelWN%kz(k)*((real(0.5_WP*(AUk%values(i,j,k) + conjg(AUk%values(i,j,k))),WP)*eAy) &
                                 & -(real(0.5_WP*(AVk%values(i,j,k) + conjg(AVk%values(i,j,k))),WP)*eAx))
                 !A2=K.(Uki times eb)
                 A2 = AVelWN%kx(i)*((real(0.5_WP*(-ii)*(AVk%values(i,j,k) - conjg(AVk%values(i,j,k))),WP)*eBz) &
                                 & -(real(0.5_WP*(-ii)*(AWk%values(i,j,k) - conjg(AWk%values(i,j,k))),WP)*eBy))&
                  & + AVelWN%ky(j)*((real(0.5_WP*(-ii)*(AWk%values(i,j,k) - conjg(AWk%values(i,j,k))),WP)*eBx) &
                                 & -(real(0.5_WP*(-ii)*(AUk%values(i,j,k) - conjg(AUk%values(i,j,k))),WP)*eBz))&
                  & + AVelWN%kz(k)*((real(0.5_WP*(-ii)*(AUk%values(i,j,k) - conjg(AUk%values(i,j,k))),WP)*eBy) &
                                 & -(real(0.5_WP*(-ii)*(AVk%values(i,j,k) - conjg(AVk%values(i,j,k))),WP)*eBx))


                  
                 

         !        rxiA = eAx*real(0.5_WP*(AUk%values(i,j,k) + conjg(AUk%values(i,j,k))),WP)&
         !         &   + eAy*real(0.5_WP*(AVk%values(i,j,k) + conjg(AVk%values(i,j,k))),WP)&
         !         &   + eAz*real(0.5_WP*(AWk%values(i,j,k) + conjg(AWk%values(i,j,k))),WP)

                 UkrB = eBx*real(0.5_WP*(AUk%values(i,j,k) + conjg(AUk%values(i,j,k))),WP)&
                  &   + eBy*real(0.5_WP*(AVk%values(i,j,k) + conjg(AVk%values(i,j,k))),WP)&
                  &   + eBz*real(0.5_WP*(AWk%values(i,j,k) + conjg(AWk%values(i,j,k))),WP)

                 UkiA = eAx*real(0.5_WP*(-ii)*(AUk%values(i,j,k) - conjg(AUk%values(i,j,k))),WP)&
                  &   + eAy*real(0.5_WP*(-ii)*(AVk%values(i,j,k) - conjg(AVk%values(i,j,k))),WP)&
                  &   + eAz*real(0.5_WP*(-ii)*(AWk%values(i,j,k) - conjg(AWk%values(i,j,k))),WP)

          !       ixiB = eBx*real(0.5_WP*(-ii)*(AUk%values(i,j,k) - conjg(AUk%values(i,j,k))),WP)&
          !        &   + eBy*real(0.5_WP*(-ii)*(AVk%values(i,j,k) - conjg(AVk%values(i,j,k))),WP)&
          !        &   + eBz*real(0.5_WP*(-ii)*(AWk%values(i,j,k) - conjg(AWk%values(i,j,k))),WP)

                    den = UkiA*A2 +UkrB*A1 !-ga*ixi1 + gb*(sin(psi)*rxi2 - cos(psi)*ixi2)
!            
                 IF (den .eq. 0.0_WP) THEN
!                     f = -1.0_WP/(4.0_WP*kk**2)
                    af =0.0_WP 
                    bf =0.0_WP 
                 ELSE
                    f = 0.001_WP/(4.0_WP*kk**2)
                    af = -f*A2/den 
                    bf =  f*A1/den 
                 ENDIF

                    Anlx%values(i,j,k) = Anlx%values(i,j,k) + (AVelWN%ky(j)*(af*eAz+ii*bf*eBz) - AVelWN%kz(k)*(af*eAy+ii*bf*eBy)) ! af*eAx+ii*bf*eBx)
                    Anly%values(i,j,k) = Anly%values(i,j,k) + (AVelWN%kz(k)*(af*eAx+ii*bf*eBx) - AVelWN%kx(i)*(af*eAz+ii*bf*eBz)) !(af*eAy+ii*bf*eBy)
                    Anlz%values(i,j,k) = Anlz%values(i,j,k) + (AVelWN%kx(i)*(af*eAy+ii*bf*eBy) - AVelWN%ky(j)*(af*eAx+ii*bf*eBx)) !(af*eAz+ii*bf*eBz)

    SommeEnerFluctuation = SommeEnerFluctuation +&
                        & real(0.5_WP*( ( AUk%values(i,j,k)*conjg(AVelWN%ky(j)*(af*eAz+ii*bf*eBz) - AVelWN%kz(k)*(af*eAy+ii*bf*eBy))&
                                      & + AVk%values(i,j,k)*conjg(AVelWN%kz(k)*(af*eAx+ii*bf*eBx) - AVelWN%kx(i)*(af*eAz+ii*bf*eBz))&
                                      & + AWk%values(i,j,k)*conjg(AVelWN%kx(i)*(af*eAy+ii*bf*eBy) - AVelWN%ky(j)*(af*eAx+ii*bf*eBx))   )&
                                      & + ( conjg(AUk%values(i,j,k))*(AVelWN%ky(j)*(af*eAz+ii*bf*eBz) - AVelWN%kz(k)*(af*eAy+ii*bf*eBy))&
                                      &   + conjg(AVk%values(i,j,k))*(AVelWN%kz(k)*(af*eAx+ii*bf*eBx) - AVelWN%kx(i)*(af*eAz+ii*bf*eBz))&
                                      &   + conjg(AWk%values(i,j,k))*(AVelWN%kx(i)*(af*eAy+ii*bf*eBy) - AVelWN%ky(j)*(af*eAx+ii*bf*eBx)) ) ),WP) 

 SommeSelfCorForcing = SommeSelfCorForcing +&
                    & real(0.5_WP*( (AVelWN%ky(j)*(af*eAz+ii*bf*eBz) -AVelWN%kz(k)*(af*eAy+ii*bf*eBy))*conjg(AVelWN%ky(j)*(af*eAz+ii*bf*eBz) -AVelWN%kz(k)*(af*eAy+ii*bf*eBy))&
                                 & +(AVelWN%kz(k)*(af*eAx+ii*bf*eBx) -AVelWN%kx(i)*(af*eAz+ii*bf*eBz))*conjg(AVelWN%kz(k)*(af*eAx+ii*bf*eBx) -AVelWN%kx(i)*(af*eAz+ii*bf*eBz))&
                                 & +(AVelWN%kx(i)*(af*eAy+ii*bf*eBy) -AVelWN%ky(j)*(af*eAx+ii*bf*eBx))*conjg(AVelWN%kx(i)*(af*eAy+ii*bf*eBy) -AVelWN%ky(j)*(af*eAx+ii*bf*eBx))   )&    
                       & ,WP)

             ENDIF
           ENDIF
        ENDDO
     ENDDO
  ENDDO
  print*,'integral des fluctuation d énergie,selfcorelation forcing',SommeEnerFluctuation,SommeSelfCorForcing

END SUBROUTINE force_Helicity_no_Energy

END MODULE forcing
