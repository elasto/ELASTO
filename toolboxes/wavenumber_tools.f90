!USEFORTEST toolbox
!USEFORTEST avgcond 
!USEFORTEST postprocess
!> @addtogroup toolbox 
!! @{
!------------------------------------------------------------------------------
!
! MODULE: toolbox 
!
!> @author
!> Antoine Vollant
!
! DESCRIPTION: 
!> The aim of this module is to provide elementaries piece of code which can be use elsewhere.
!> BEWARE there are no systematic tests in these routines in order to not degrade performances
!> 
!------------------------------------------------------------------------------


module wavenumber_tools

  use precision_tools 

  IMPLICIT NONE

  public 
 
  TYPE wavenumbers
      REAL(WP), DIMENSION(:), POINTER :: kx => NULL()
      REAL(WP), DIMENSION(:), POINTER :: ky => NULL()
      REAL(WP), DIMENSION(:), POINTER :: kz => NULL()
      REAL(WP) :: xkmax,ykmax,zkmax
      REAL(WP) :: xkcut,ykcut,zkcut
  END TYPE wavenumbers

CONTAINS

!------------------------------------------------------------------------------
!> @author 
!> Patrick BEGOU, LEGI
!>
!> @details
!> This function allocate and initialize the wave numbers arrays. It should
!> be called for velocities and for each scalar field.
!> @param [in,out] val the WaveNumber structure to initialize
!> @param [in] spec_rank my MPI processe rank in the range [0:nbcpus-1]
!> @param [in] nx the points number along X axis
!> @param [in] ny the points number along Y axis
!> @param [in] nz the points number along Z axis
!> @return .TRUE. if initialization is successfull. 
!------------------------------------------------------------------------------
FUNCTION initWN(val,spec_rank,nx,ny,nz) result(success)

  !I/O data 
  TYPE(WaveNumbers), INTENT(INOUT)::val
  INTEGER, INTENT(IN) :: spec_rank,nx,ny,nz
  LOGICAL             :: success 
  !Local data 
  INTEGER      :: ires 
  success = .false.
 
!  Allocate memory storage
  IF(ASSOCIATED(val%kx)) DEALLOCATE(val%kx)
  IF(ASSOCIATED(val%ky)) DEALLOCATE(val%ky)
  IF(ASSOCIATED(val%kz)) DEALLOCATE(val%kz)
  ALLOCATE(val%kx(1:(nx/2)+1),val%ky(1:ny), val%kz(1:nz), stat=ires)
  IF (ires .NE.0) THEN
   WRITE(6,'(a,i0,a)')'[ERROR] on process ',spec_rank,&
   & ': not enought memory for waves numbers.'
   RETURN
  ENDIF

  success = .true.

END FUNCTION initWN


!------------------------------------------------------------------------------
!> @author 
!> Patrick BEGOU, LEGI
!
!>
!> @details
!> This function allocate and initialize the wave numbers arrays. It should
!> be called for velocities and for each scalar field.
!> @param [in,out] val the WaveNumber structure to initialize
!> @param [in] spec_rank my MPI processe rank in the range [0:nbcpus-1]
!> @param [in] nx the points number along X axis
!> @param [in] ny the points number along Y axis
!> @param [in] nz the points number along Z axis
!> @param [in] Lx the domain lenght along the X axis
!> @param [in] Ly the domain lenght along the Y axis
!> @param [in] Lz the domain lenght along the Z axis
!> @return .TRUE. if initialization is successfull. 
!------------------------------------------------------------------------------
FUNCTION computeWN(val,spec_rank,Lx,Ly,Lz,nx,ny,nz) result(success)
  
  !I/O data
  TYPE(WaveNumbers), INTENT(INOUT)::val
  INTEGER, INTENT(IN) :: spec_rank,nx,ny,nz
  REAL(WP), INTENT(IN):: Lx,Ly,Lz
  LOGICAL             :: success 
  !Local data 
  INTEGER  :: i
  REAL(WP) :: delt

  success = .FALSE.
  
  IF( .not. ASSOCIATED(val%kx) .and. .not. ASSOCIATED(val%ky) &
    & .and. .not. ASSOCIATED(val%kz) ) THEN
    WRITE(6,'(a,i0,a)')'[ERROR] on process ',spec_rank,&
    & ': waves numbers not initialized.'
    RETURN
  ENDIF
  ! along X for Real <-> complex FFT
  delt=2.0_WP*acos(-1.0_WP)/Lx
  DO i = 1,(nx/2)+1
     val%kx(i)=REAL(i-1,WP)*delt
  ENDDO

  ! along Y for complex <-> complex FFT
  delt=2.0_WP*acos(-1.0_WP)/Ly
  DO i = 1,ny
    IF ( i .gt. ((ny/2)+1) ) THEN
       val%ky(i)=-REAL(ny+1-i,WP)*delt
    ELSE   
       val%ky(i)=REAL(i-1,WP)*delt
    ENDIF
  ENDDO

  ! along Z for complex <-> complex FFT
  delt=2.0_WP*acos(-1.0_WP)/Lz
  DO i = 1,nz
    IF ( i .gt. ((nz/2)+1) ) THEN
       val%kz(i)=-REAL(nz+1-i,WP)*delt
    ELSE   
       val%kz(i)=REAL(i-1,WP)*delt
    ENDIF
  ENDDO

  val%xkmax = real((nx/2),WP)*2.0_WP*acos(-1.0_WP)/Lx ! val%kx(nx/2+1)
  val%ykmax = real((ny/2),WP)*2.0_WP*acos(-1.0_WP)/Ly ! val%ky(ny/2+1)
  val%zkmax = real((nz/2),WP)*2.0_WP*acos(-1.0_WP)/Lz ! val%kz(nz/2+1)
  
  val%xkcut = (2.0_WP/3.0_WP)*val%xkmax
  val%ykcut = (2.0_WP/3.0_WP)*val%ykmax
  val%zkcut = (2.0_WP/3.0_WP)*val%zkmax

  success = .TRUE.

end function computeWN 

!------------------------------------------------------------------------------
!> @author 
!> Patrick BEGOU, LEGI
!>
!> @details
!> This function allocate and initialize the wave numbers arrays. It should
!> be called for velocities and for each scalar field.
!> @param [in,out] val the WaveNumber structure to initialize
!> @param [in] spec_rank my MPI processe rank in the range [0:nbcpus-1]
!> @param [in] nx the points number along X axis
!> @param [in] ny the points number along Y axis
!> @param [in] nz the points number along Z axis
!> @return .TRUE. if initialization is successfull. 
!------------------------------------------------------------------------------
FUNCTION deleteWN(val) result(success)

  IMPLICIT NONE
  
  TYPE(WaveNumbers), INTENT(INOUT)::val
  LOGICAL                         :: success
 
  success = .FALSE.

  DEALLOCATE(val%kx)
  DEALLOCATE(val%ky)
  DEALLOCATE(val%kz)

  success = .TRUE.

END FUNCTION deleteWN


!------------------------------------------------------------------------------
!> @author 
!> Patrick BEGOU, LEGI
!
!>
!> @details
!> This function allocate and initialize the wave numbers arrays. It should
!> be called for velocities and for each scalar field.
!> @param [in,out] val the WaveNumber structure to initialize
!> @param [in] spec_rank my MPI processe rank in the range [0:nbcpus-1]
!> @param [in] nx the points number along X axis
!> @param [in] ny the points number along Y axis
!> @param [in] nz the points number along Z axis
!> @param [in] Lx the domain lenght along the X axis
!> @param [in] Ly the domain lenght along the Y axis
!> @param [in] Lz the domain lenght along the Z axis
!> @return .TRUE. if initialization is successfull. 
!------------------------------------------------------------------------------
FUNCTION wavesNumber_init(val,spec_rank,Lx,Ly,Lz,nx,ny,nz) result(success)

  !I/O data
  TYPE(WaveNumbers), INTENT(INOUT) :: val
  INTEGER, INTENT(IN)              :: spec_rank,nx,ny,nz
  REAL(WP), INTENT(IN)             :: Lx,Ly,Lz
  LOGICAL                          :: success 

  success = .false.
  if (.not. initWN(val,spec_rank,nx,ny,nz)) then 
    write (6,'(a)') '[ERROR] in wavesNumber_init, initWN failed : abort'
    return 
  endif
  if (.not. computeWN(val,spec_rank,Lx,Ly,Lz,nx,ny,nz)) then
    write (6,'(a)') '[ERROR] wavesNumber_init, computeWN failed : abort'
    return
  endif

  success = .true.

END FUNCTION wavesNumber_init



!------------------------------------------------------------------------------
!> @author 
!> Patrick BEGOU, LEGI
!
!>
!> @details
!> This function free the memory used by a variable of type WaveNumbers
!------------------------------------------------------------------------------
SUBROUTINE deleteWavesNumber(val)
  IMPLICIT NONE
  TYPE(WaveNumbers), INTENT(INOUT)::val
  IF (ASSOCIATED(val%kx)) DEALLOCATE(val%kx)
  IF (ASSOCIATED(val%ky)) DEALLOCATE(val%ky)
  IF (ASSOCIATED(val%kz)) DEALLOCATE(val%kz)
  val%kx=>  NULL()
  val%ky=>  NULL()
  val%kz=>  NULL()
 END SUBROUTINE deleteWavesNumber

end module wavenumber_tools 
!> @}
