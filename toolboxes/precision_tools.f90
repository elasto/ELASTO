!USEFORTEST toolbox
!USEFORTEST postprocess
!USEFORTEST advec
!USEFORTEST io
!USEFORTEST topo
!USEFORTEST avgcond
!USEFORTEST interpolation
!> @addtogroup toolbox
!! @{
!------------------------------------------------------------------------------
!
! MODULE: precision
!
!> @author
!> Guillaume Balarac, LEGI
!
! DESCRIPTION:
!> The aim of this module is set some parameters to fix the working data
!> representation in the code. It is set to double precision for REAL.
!------------------------------------------------------------------------------

MODULE precision_tools

  IMPLICIT NONE
  INTEGER, PARAMETER  :: SP = kind(1.0)
  INTEGER, PARAMETER  :: DP = kind(1.0d0)
  INTEGER, PARAMETER  :: WP = DP
  REAL(WP), PRIVATE   :: sample_real_at_WP
  REAL(WP), PARAMETER :: MAX_REAL_WP = HUGE(sample_real_at_WP)
  INTEGER, PRIVATE    :: sample_int
  INTEGER, PARAMETER  :: MAX_INTEGER = HUGE(sample_int)
  INTEGER, PARAMETER  :: DI = selected_int_kind(r=12)
  !> the MPI type for REAL exchanges in simple or double precision
  INTEGER, PUBLIC     :: MPI_REAL_WP
  !> the MPI type for COMPLEX exchanges in simple or double precision
  INTEGER, PUBLIC     :: MPI_COMPLEX_WP
  !> the string size 
  INTEGER, PARAMETER  :: str_short  = 8
  INTEGER, PARAMETER  :: str_medium = 64
  INTEGER, PARAMETER  :: str_long   = 4096

END MODULE precision_tools
!> @}
