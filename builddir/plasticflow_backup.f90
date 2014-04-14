!---------------------------------------------------------------------------------
MODULE plasticflow

  USE precision_tools
  USE datalayout
  USE wavenumber_tools

  IMPLICIT NONE

  PUBLIC  :: initplasticflow

CONTAINS

!---------------------------------------------------------------------------------
!> Initialisation with cosinus.
!! @autor = Jean-Baptiste Lagaert, Guillaume Balarac, LEGI
!
!> @param[in,out]   sigma               = longitudinal velocity storage to initialize
!> @param[in]       sigmaWN           = the wavenumbers for U,V and W
!> @param[in]       nbcpus          = the total amount of processes (or cpus)
!> @param[in]       spec_rank       = my MPI processe rank in the range [0:nbcpus-1]
!> @return          ok              = LOGICAL set to .TRUE. if all initialization success.
!---------------------------------------------------------------------------------
FUNCTION initplasticflow(sigma,sigmaWN,nbcpus,spec_rank) result(ok)

    IMPLICIT NONE

    !I/O data
    TYPE(REAL_DATA_LAYOUT), intent(inout)           :: sigma
    TYPE(WaveNumbers), intent(in)                   :: sigmaWN
    INTEGER                                         :: nbcpus, spec_rank
    LOGICAL                                         :: ok

    !Local data
    REAL(WP)        :: zz,  dz                       
    REAL(WP)        :: yy,  dy                       
    REAL(WP)        :: xx,  dx                       
    INTEGER         :: i,j,k
 

    ok = .false.

    ! === Initialization for U
    DO k = sigma%zmin,sigma%zmax
       DO j = sigma%ymin,sigma%ymax
          DO i = sigma%xmin,sigma%xmax
             sigma%values(i,j,k) = real(spec_rank)
          ENDDO
       ENDDO
    ENDDO 
  

    ok = .true.
END FUNCTION initplasticflow

END MODULE plasticflow
