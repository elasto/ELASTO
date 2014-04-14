!> @addtogroup forcing
!! @{
!------------------------------------------------------------------------------
!
! MODULE: forcing_scalar
!
!
! DESCRIPTION: 
!> The module ``forcing_scalar'' gathers all procedure to add forcing term on
!! the scalar field and thus avoid its dissipation and maintaining a statiscally
!! stationary state.
!! @details
!!     This module contain mean gradient forcing (and null) forcing and forcing
!! setup for each scalar field. The implementation allow to avoid a if call.
!! (like number of scalar, presence of magnetic field, ...) but also
!! about I/O and post-process.
!!     More precisly, these values/parameters have (eventually) default value
!! defines in this module and are update by using values defined in the
!! "input" file.
!!     All variables defined in this module are protected : they must be
!! visible by any other modules but must not be changed.
!!
!
!> @author
!! Jean-Baptiste Lagaert, LEGI
!
!------------------------------------------------------------------------------

module forcing_scalar

    use precision_tools
    use wavenumber_tools
    use datalayout
    use random_tools 

    implicit none

    ! ##########################################
    ! ##########      Procedures      ##########
    ! ##########################################
    !
    ! ===== Public procedures =====
    !
    ! ##########################################
    ! ##########      Variables       ##########
    ! ##########################################

    !> Different implementation of the forcing term - mean gradient
    interface forcing_scal_mean
        module procedure forcing_scal_mean_scalar, forcing_scal_mean_source_term
    end interface forcing_scal_mean

    !> Different implementation of the forcing term - forcing in spectral space.
    interface forcing_spec_scal
        module procedure forcing_spec_scal_NL, forcing_spec_scal_no_NL
    end interface forcing_spec_scal
    


contains

!> Add "mean gradient" forcing directly to the scalar
!! @param[in]       spec_rank   = my MPI processe rank in the range [0:nbcpus-1]
!! @param[in,out]   scalar      = scalar field (in real space)
!! @param[in]       Vy          = velocity component along Y (in real space)
!! @param[in]       dt          = time step
!! @author Jean-Baptiste Lagaert
subroutine forcing_scal_mean_scalar(scalar, Vy, dt, spec_rank)

    type(REAL_DATA_LAYOUT), intent(in)      :: Vy
    type(REAL_DATA_LAYOUT), intent(inout)   :: scalar
    real(WP), intent(in)                    :: dt
    integer, intent(in)                     :: spec_rank


    if (.not. (samelayout(scalar,Vy))) then 
        if (spec_rank==0) write(6,'(a)')'[ERROR]Forcing_scal_mean_scalar: internal error, scalar and Vy do not match '
        stop
    end if

    scalar%values = scalar%values + dt*Vy%values

end subroutine forcing_scal_mean_scalar


!> Add "mean gradient" forcing to the source term to apply to the scalar
!! @param[in]       spec_rank   = my MPI processe rank in the range [0:nbcpus-1]
!! @param[in,out]   source_term = source term from the scalar equation (in real space)
!! @param[in]       Vy          = velocity component along Y (in real space)
subroutine forcing_scal_mean_source_term(source_term, Vy, spec_rank)

    type(COMPLEX_DATA_LAYOUT), intent(in)      :: Vy
    type(COMPLEX_DATA_LAYOUT), intent(inout)   :: source_term
    integer, intent(in)                     :: spec_rank


    if (.not. (samelayout(source_term,Vy))) then 
        if (spec_rank==0) write(6,'(a)') "[ERROR] Forcing_scal_mean_source_term: internal error, Vy and source_term do not match"
        stop
    end if

    source_term%values = source_term%values + Vy%values
end subroutine forcing_scal_mean_source_term

!> Compute spectral forcing for scalar 
!! @author Guillaume Balarac
!!    @param[in]  AVelWN        = Scalar wavenumbers  
!!    @param[in,out]  Anlx      = Scalar term
!!    @param[in]  amplitude_sf  = Forcing amplitude
!!    @param[in]  kmin_sf          = lower bound of wavenumber shell of forcing 
!!    @param[in]  kmax_sf          = upper bound of wavenumber shell of forcing 
!!    @param[in]  dt            = time step
SUBROUTINE forcing_spec_scal_NL(AVelWN,Anlx,amplitude_sf,kmin_sf,kmax_sf,dt)

  IMPLICIT NONE
  
  TYPE(WaveNumbers), INTENT(IN)              :: AVelWN
  REAL(WP), INTENT(IN)                       :: dt, kmin_sf, kmax_sf, amplitude_sf
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT)   :: Anlx
  INTEGER :: i,j,k
  REAL(WP) :: rand,pi,kk
  COMPLEX(WP), PARAMETER :: ii=(0.0_WP,1.0_WP)
  pi = acos(-1.0_WP)

  DO k = Anlx%zmin,Anlx%zmax
    DO j =  Anlx%ymin,Anlx%ymax
      DO i =  Anlx%xmin,Anlx%xmax

           kk = sqrt(AVelWN%kx(i)*AVelWN%kx(i) + &
                   & AVelWN%ky(j)*AVelWN%ky(j) + &
                   & AVelWN%kz(k)*AVelWN%kz(k)    )

           IF (kk.ge.kmin_sf.and.kk.le.kmax_sf) then

                 call random_number(rand)
                 !rand = 2.0_WP*(rand - 0.5_WP)
                 
                 Anlx%values(i,j,k) = Anlx%values(i,j,k) + rand*amplitude_sf*dt
                 
           ENDIF
           
        ENDDO
     ENDDO
  ENDDO
  
END SUBROUTINE forcing_spec_scal_NL

!> Compute spectral forcing for scalar 
!! @author Guillaume Balarac
!!    @param[in]  AVelWN        = Scalar wavenumbers  
!!    @param[in,out]  Anlx      = Scalar term
!!    @param[in]  amplitude_sf  = Forcing amplitude
!!    @param[in]  kmin_sf          = lower bound of wavenumber shell of forcing 
!!    @param[in]  kmax_sf          = upper bound of wavenumber shell of forcing 
SUBROUTINE forcing_spec_scal_no_NL(AVelWN,Anlx,amplitude_sf,kmin_sf,kmax_sf)

  TYPE(WaveNumbers), INTENT(IN)              :: AVelWN
  REAL(WP), INTENT(IN)                       :: kmin_sf, kmax_sf, amplitude_sf
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT)   :: Anlx
  INTEGER :: i,j,k
  REAL(WP) :: rand,pi,kk
  COMPLEX(WP), PARAMETER :: ii=(0.0_WP,1.0_WP)
  pi = acos(-1.0_WP)

  DO k = Anlx%zmin,Anlx%zmax
    DO j =  Anlx%ymin,Anlx%ymax
      DO i =  Anlx%xmin,Anlx%xmax

           kk = sqrt(AVelWN%kx(i)*AVelWN%kx(i) + &
                   & AVelWN%ky(j)*AVelWN%ky(j) + &
                   & AVelWN%kz(k)*AVelWN%kz(k)    )

           IF (kk.ge.kmin_sf.and.kk.le.kmax_sf) then

                 call random_number(rand)
                 !rand = 2.0_WP*(rand - 0.5_WP)
                 
                 Anlx%values(i,j,k) = Anlx%values(i,j,k) + rand*amplitude_sf
                 
           ENDIF
           
        ENDDO
     ENDDO
  ENDDO
  
END SUBROUTINE forcing_spec_scal_no_NL

end module forcing_scalar
!! @}
