!USEFORTEST toolbox
!USEFORTEST avgcond 
!USEFORTEST postprocess
!> @addtogroup toolbox
!! @{
!------------------------------------------------------------------------------
!
! MODULE: differential_tools
!
!! DESCRIPTION:
!> This module implements the NS solver for scalars and velocities.
!> @author
!> Guillaume Balarac et Patrick BEGOU, LEGI
!>
!> @details
!! This solver implements 3 choices for dealiasing: none, isotropic and anisotropic.
!! The last one is usefull when the number point in the 3 dimensions are very heterogeneous
!! to avoid loosing to much informations.
!! \image html dealias.png "Influence of dealiasing method on scalar min/max behavior."
!> @author
!> Guillaume Balarac et Patrick BEGOU, LEGI
!------------------------------------------------------------------------------
module differential_tools

    USE precision_tools
    USE datalayout
    USE wavenumber_tools
    USE transforms_tools

    IMPLICIT NONE

   public
   !interface for differential operators in spectral space
   interface computeDivergenceK
     module procedure computeDivergenceKmain
     module procedure computeDivergenceK1
   end interface computeDivergenceK

   interface computeGradientK
     module procedure computeGradientKmain
     module procedure computeGradientK1
   end interface computeGradientK

   interface computeCurlK
     module procedure computeCurlKmain
     module procedure computeCurlK1
     module procedure computeCurlK2
   end interface computeCurlK

   interface computeLaplacianK
     module procedure computeLaplacianScaKmain
     module procedure computeLaplacianVecKmain
     module procedure computeLaplacianVecK1
   end interface computeLaplacianK

   !interface for differential operators in physical space
   interface computeFieldLaplacian
     module procedure computeFieldLaplacianScalar
     module procedure computeFieldLaplacianArray_scaField
     module procedure computeFieldLaplacianArray_vecField
   end interface computeFieldLaplacian

   interface computeFieldCurl
     module procedure computeFieldCurl_scaField
     module procedure computeFieldCurl_vecField
   end interface computeFieldCurl

   interface computeFieldGradient
     module procedure computeFieldGradient_scaField
     module procedure computeFieldGradient_vecField
   end interface computeFieldGradient

   interface computeFieldDivergence
     module procedure computeFieldDivergence_scaField
     module procedure computeFieldDivergence_vecField
   end interface computeFieldDivergence

   interface computeDerivationField
     module procedure computeDerivationField_Phy 
     module procedure computeDerivationField_Spe 
   end interface computeDerivationField

   private foubuff
   private fouBbuff
   private diff2K
   private diffK

   !> work arrays needed for divergence computation of U field
   !TYPE(REAL_DATA_LAYOUT), SAVE     :: phybuff
   TYPE(COMPLEX_DATA_LAYOUT), SAVE:: foubuff
   !> work arrays needed for divergence computation of B field
   !TYPE(REAL_DATA_LAYOUT), SAVE     :: phyBbuff
   TYPE(COMPLEX_DATA_LAYOUT), SAVE:: fouBbuff

contains

!------------------------------------------------------------------------------
!> @author
!> Jean-Baptiste Lagaert
!>
!> @details
!> Allocate buffers if needed.
!> @param[in,out]   Uk      = One of velocity field component in spectral space
!> @param[in,out]   Bk      = One of magnetic field component in spectral space 
!> @return          res     = .TRUE. if all is succcessfull, .FALSE. elsewhere.
!------------------------------------------------------------------------------
function diff_tool_initBuff(Uk, Bk)

   type(COMPLEX_DATA_LAYOUT), intent(inout)            :: Uk
   type(COMPLEX_DATA_LAYOUT), intent(inout), optional  :: Bk
   logical                                     :: diff_tool_initBuff

   diff_tool_initBuff = .false.

   ! IF(spec_rank.EQ.0) WRITE(6,'(a)')'[WARNING] Div(U) not set to 0.0 in init-solver!'
   ! print and set div(U) to 0
   IF (.not. copyStructOnly(Uk,foubuff)) return
   if(present(Bk)) then
       IF (.not. copyStructOnly(Bk,fouBbuff)) return
   end if

   diff_tool_initBuff= .true.

end function diff_tool_initBuff


!------------------------------------------------------------------------------
!> @author
!> Jean-Baptiste Lagaert
!>
!> @details
!> Deallocate buffers if needed.
!> @param[in]   mhd  = logical value to set if mhd is called 
!------------------------------------------------------------------------------
subroutine diff_tool_deleteBuff(mhd)

    logical, intent(in) :: mhd

    call deleteDataLayout(fouBbuff)
    if(mhd) call deleteDataLayout(foubuff)

end subroutine diff_tool_deleteBuff


!------------------------------------------------------------------------------
!> @author
!> Patrick BEGOU, LEGI and Jean-Baptiste Lagaert
!
!>
!> @details
!> Set the div if the input vector to 0.0
!> @param[in]       spec_rank   = my MPI processe rank in the range [0:nbcpus-1]
!> @param[in,out]   vectXk      = vector component along X, in spectral space.
!> @param[in,out]   vectYk      = vector component along Y, in spectral space.
!> @param[in,out]   vectZk      = vector component along Z, in spectral space.
!> @param[in]       WN          = WN numbers associated to the vector
!> @return          res         = .TRUE. if all is succcessfull, .FALSE. elsewhere.
!------------------------------------------------------------------------------
FUNCTION NullifyDiv(vectXk, vectYk, vectZk, WN, spec_rank) result(res)

    !I/O data
    INTEGER,INTENT(IN)                      :: spec_rank
    TYPE(COMPLEX_DATA_LAYOUT),INTENT(INOUT) :: vectXk,vectYk,vectZk
    TYPE(WaveNumbers), INTENT(IN)           :: WN
    LOGICAL                                 :: res

    !Local data
    INTEGER :: i,j,k
    REAL(WP) :: kk
    COMPLEX(WP) :: Pressk
    COMPLEX(WP), PARAMETER :: ii=(0.0_WP,1.0_WP)

    res=.FALSE.
    ! Projection: assert that div(U) == 0
    IF (.NOT. (samelayout(vectXk,vectYk,vectZk))) THEN
        IF (spec_rank.EQ.0) WRITE(6,'(a)')'[ERROR] NullifyDiv: internal error, vectXk, vectYk and vectZk do not match '
        RETURN
    ENDIF
    DO k = vectXk%zmin,vectXk%zmax
        DO j =  vectXk%ymin,vectXk%ymax
            DO i =  vectXk%xmin,vectXk%xmax
            kk =  WN%kx(i)*WN%kx(i) + &
                & WN%ky(j)*WN%ky(j) + &
                & WN%kz(k)*WN%kz(k)
            IF (kk.gt.0.0_WP) THEN
                Pressk =  ( WN%kx(i)*vectXk%values(i,j,k) + &
                    & WN%ky(j)*vectYk%values(i,j,k) + &
                    & WN%kz(k)*vectZk%values(i,j,k) )/kk
                vectXk%values(i,j,k) = vectXk%values(i,j,k) - WN%kx(i)*Pressk
                vectYk%values(i,j,k) = vectYk%values(i,j,k) - WN%ky(j)*Pressk
                vectZk%values(i,j,k) = vectZk%values(i,j,k) - WN%kz(k)*Pressk
            ENDIF
            ENDDO
        ENDDO
    ENDDO
    res=.TRUE.

END FUNCTION NullifyDiv

!------------------------------------------------------------------------------
!> @author
!> Patrick BEGOU, LEGI and Jean-Baptiste Lagaert
!> @details
!> Force the first wave number to zero in spectral space
!> @param[in]       spec_rank   = my MPI processe rank in the range [0:nbcpus-1]
!> @param[in,out]   vectXk      = vector component along X, in spectral space.
!> @param[in,out]   vectYk      = vector component along Y, in spectral space.
!> @param[in,out]   vectZk      = vector component along Z, in spectral space.
!> @param[in]       WN          = WN numbers associated to the vector
!> @return          res         = .TRUE. if all is succcessfull, .FALSE. elsewhere.
!------------------------------------------------------------------------------
FUNCTION NullifyFirstWN(vectXk, vectYk, vectZk, WN, spec_rank) result(res)
    IMPLICIT NONE

    !I/O data
    INTEGER,INTENT(IN)                      :: spec_rank
    TYPE(COMPLEX_DATA_LAYOUT),INTENT(INOUT) :: vectXk,vectYk,vectZk
    TYPE(WaveNumbers), INTENT(IN)           :: WN
    LOGICAL                                 :: res

    !Local data
    INTEGER :: i,j,k
    REAL(WP) :: kk,eps

    res = .false.

    if (vectXk%Lx .ne. vectXk%Ly .or. vectXk%Lx .ne. vectXk%Lz) then
      IF (spec_rank.EQ.0) WRITE(6,'(a)') '[ERROR] Nullify spectral values not implemented  &
                                      & for non cubic box: abort!'
      return
    endif
    eps = acos(-1.0_WP)/vectXk%Lx
    IF (spec_rank.EQ.0) WRITE(6,'(a,1x,f8.5)') '[INFO] Nullify spectral values for wave number &
                                      & under',eps
    do k = vectXk%zmin,vectXk%zmax
        do j =  vectXk%ymin,vectXk%ymax
            do i = vectXk%xmin,vectXk%xmax
                kk = WN%kx(i)*WN%kx(i) + &
                   & WN%ky(j)*WN%ky(j) + &
                   & WN%kz(k)*WN%kz(k)
                kk = sqrt(real(kk,WP))
                if ((kk .lt.eps)) then
                  vectXk%values(i,j,k) = (0.0_WP,0.0_WP)
                  vectYk%values(i,j,k) = (0.0_WP,0.0_WP)
                  vectZk%values(i,j,k) = (0.0_WP,0.0_WP)
                endif
            end do
        end do
    end do

    res = .true.

END FUNCTION NullifyFirstWN



!------------------------------------------------------------------------------
!> @author
!> Patrick BEGOU, LEGI
!>
!> @details
!> Show div(U)
!> @param [in,out] spec_rank my MPI processe rank in the range [0:nbcpus-1]
!> @return .TRUE. if all is succcessfull, .FALSE. elsewhere.
!------------------------------------------------------------------------------
FUNCTION showDivU(Uk,Vk,Wk,VelWN,spec_rank) result(res)


    type(complex_data_layout),intent(in) :: Uk,Vk,Wk 
    type(wavenumbers),intent(in) :: VelWN
    INTEGER,INTENT(IN)  :: spec_rank
    LOGICAL             :: res

    REAL(WP)            :: val

    res=.FALSE.

    ! Projection: assert that div(U) == 0
    IF (.NOT. samelayout(Uk,Vk,Wk,foubuff) ) THEN
         IF (spec_rank.EQ.0) WRITE(6,'(a)')'[ERROR]showDivU: internal error, arrays do not match '
         RETURN
    ENDIF
    val = computeDiv( spec_rank,Uk,Vk,Wk,foubuff,VelWN )
    IF (spec_rank.EQ.0) WRITE(6,'(a,g15.8)')'    [INFO] Divergence max of U field is ',val

    res=.TRUE.

END FUNCTION showDivU


!------------------------------------------------------------------------------
!> @author
!> Patrick BEGOU and Guillaume BALARAC, LEGI
!>
!> @details
!> Show div(B)
!> @param [in,out] spec_rank my MPI processe rank in the range [0:nbcpus-1]
!> @return .TRUE. if all is succcessfull, .FALSE. elsewhere.
!------------------------------------------------------------------------------
FUNCTION showDivB(Bk,BfieldWN,spec_rank) result(res)

    type(complex_data_layout),dimension(:),intent(in) :: Bk
    type(wavenumbers),intent(in) :: BfieldWN 
    INTEGER,INTENT(IN)  :: spec_rank
    LOGICAL             :: res
    REAL(WP)            :: val

    res=.FALSE.

    ! Projection: assert that div(B) == 0
    IF (.NOT. (samelayout(Bk(1),Bk(2),Bk(3),fouBbuff))) THEN
         IF (spec_rank.EQ.0) WRITE(6,'(a)')'[ERROR]showDivB: internal error, arrays do not match '
         RETURN
    ENDIF
    val = computeDiv( spec_rank,Bk(1),Bk(2),Bk(3),fouBbuff,BfieldWN )
    IF (spec_rank.EQ.0) WRITE(6,'(a,g15.8)')'    [INFO] Divergence max of B field is ',val

    res=.TRUE.

END FUNCTION showDivB


!------------------------------------------------------------------------------
!> @author
!> Patrick BEGOU and Guillaume BALARAC, LEGI
!>
!> @details
!> computeDiv
!> @param [in,out] spec_rank my MPI processe rank in the range [0:nbcpus-1]
!> @param [in] in1 first component of array field
!> @param [in] in2 second component of array field
!> @param [in] in3 third component of array field
!> @param [in,out] temp sscalar field which store the divergence in spetral space
!> @param [in] wave wave number associated to the in array field
!> @return .TRUE. if all is succcessfull, .FALSE. elsewhere.
!------------------------------------------------------------------------------
FUNCTION computeDiv(spec_rank,in1,in2,in3,temp,wave) result(res)

    use parallel_tools

  IMPLICIT NONE
  TYPE(COMPLEX_DATA_LAYOUT),INTENT(IN)   :: in1,in2,in3
  TYPE(COMPLEX_DATA_LAYOUT),INTENT(INOUT):: temp
  TYPE(WAVENUMBERS),        INTENT(IN)   :: wave
  INTEGER,                  INTENT(IN)   :: spec_rank
  REAL(WP)                               :: res

  REAL(WP)                               :: sumLoc
  INTEGER                                :: i,j,k

  sumLoc =0.0_WP

  CALL computeDivergenceKmain(wave,in1,in2,in3,temp)
  do k = temp%zmin , temp%zmax
    do j = temp%ymin , temp%ymax
      do i = temp%xmin , temp%xmax
        sumLoc = sumLoc + abs( temp%values(i,j,k) )**2
      enddo
    enddo
  enddo
  res = DoLocalSum(spec_rank, sumLoc)
  res = sqrt( res )

END FUNCTION computeDiv

!---------------------------------------------------------------------------
!> @details 
!> Compute Laplacian of array field in spectral space
!! @author  Antoine Vollant, LEGI
!! @param[in]    wave waves number for field valIn
!! @param[in]    valIn1 first component of scalar field in spectral space
!! @param[in]    valIn2 second component of scalar field in spectral space
!! @param[in]    valIn3 third component of scalar field in spectral space  
!! @param[out]   valOut1 first component of laplacian in spectral space
!! @param[out]   valOut2 second component of laplacian in spectral space
!! @param[out]   valOut3 third component of laplacian in spectral space
!---------------------------------------------------------------------------
 subroutine computeLaplacianVecKmain(wave,valIn1,valIn2,valIn3,valOut1,valOut2,valOut3)


  !IO
  type(complex_data_layout), intent(in)    :: valIn1,valIn2,valIn3
  type(complex_data_layout), intent(inout) :: valOut1,valOut2,valOut3
  type(WaveNumbers), intent(in)        :: wave


#ifdef DEBUGTOOLBOX
    if (.not. sameLayout(valIn1,valOut1,valOut2,valOut3) .or. &
       &.not. sameLayout(valIn1,valIn2,valIn3) ) then
      write(6,'(a)')'[ERROR] in computeCurlKmain data are not the same'
    endif
    if (.not. associated(valIn1%values) .or. &
       &.not. associated(valIn2%values) .or. &
       &.not. associated(valIn3%values) .or. &
       &.not. associated(valOut1%values) .or. &
       &.not. associated(valOut2%values) .or. &
       &.not. associated(valOut3%values) .or. &
       &.not. associated(wave%kx) )  then
      write(6,'(a)')'[ERROR] in computeCurlKmain data are not initialized'
    endif
#endif
  
  call computeLaplacianScaKmain(wave,valIn1,valOut1)
  call computeLaplacianScaKmain(wave,valIn2,valOut2)
  call computeLaplacianScaKmain(wave,valIn3,valOut3)

 end subroutine computeLaplacianVecKmain


!---------------------------------------------------------------------------
!> @details 
!> Compute Laplacian of scalar field in spectral space
!! @author  Antoine Vollant, LEGI
!! @param[in]    wave waves number for field valIn
!! @param[in]    valIn  scalar field in spectral space
!! @param[out]   valOut laplacian in spectral space
!---------------------------------------------------------------------------
 subroutine computeLaplacianScaKmain(wave,valIn,valOut)


  !IO
  type(complex_data_layout), intent(in)    :: valIn
  type(complex_data_layout), intent(inout) :: valOut
  type(waveNumbers), intent(in)        :: wave

  !Internal
  integer                                  :: i,j,k

#ifdef DEBUGTOOLBOX
    if (.not. sameLayout(valIn,valOut) ) then
      write(6,'(a)')'[ERROR] in computeCurlKmain data are not the same'
    endif
    if (.not. associated(valIn%values) .or. &
       &.not. associated(valOut%values) .or. &
       &.not. associated(wave%kx) )  then
      write(6,'(a)')'[ERROR] in computeCurlKmain data are not initialized'
    endif
#endif
  
  do k=valIn%zmin,valIn%zmax 
    do j=valIn%ymin,valIn%ymax 
      do i=valIn%xmin,valIn%xmax 
        valOut%values(i,j,k) = diff2K( wave%kx(i) , valIn%values(i,j,k)) + &
                             & diff2K( wave%ky(j) , valIn%values(i,j,k)) + &
                             & diff2K( wave%kz(k) , valIn%values(i,j,k))
      enddo
    enddo
  enddo

 end subroutine computeLaplacianScaKmain

!---------------------------------------------------------------------------
!> @details 
!> Compute curl in spectral space
!! @author  Antoine Vollant, LEGI
!! @param[in]    wave waves number for field valIn
!! @param[in]    valIn1 first component of scalar field in spectral space
!! @param[in]    valIn2 second component of scalar field in spectral space
!! @param[in]    valIn3 third component of scalar field in spectral space  
!! @param[out]   valOut1 first component of curl in spectral space
!! @param[out]   valOut2 second component of curl in spectral space
!! @param[out]   valOut3 third component of curl in spectral space
!---------------------------------------------------------------------------
 subroutine computeCurlKmain(wave,valIn1,valIn2,valIn3,valOut1,valOut2,valOut3)


  !IO
  type(complex_data_layout), intent(in)    :: valIn1,valIn2,valIn3
  type(complex_data_layout), intent(inout) :: valOut1,valOut2,valOut3
  type(waveNumbers), intent(in)        :: wave

  !Internal
  integer                                  :: i,j,k

#ifdef DEBUGTOOLBOX
    if (.not. sameLayout(valIn1,valOut1,valOut2,valOut3) .or. &
       &.not. sameLayout(valIn1,valIn2,valIn3) ) then
      write(6,'(a)')'[ERROR] in computeCurlKmain data are not the same'
    endif
    if (.not. associated(valIn1%values) .or. &
       &.not. associated(valIn2%values) .or. &
       &.not. associated(valIn3%values) .or. &
       &.not. associated(valOut1%values) .or. &
       &.not. associated(valOut2%values) .or. &
       &.not. associated(valOut3%values) .or. &
       &.not. associated(wave%kx) )  then
      write(6,'(a)')'[ERROR] in computeCurlKmain data are not initialized'
    endif
#endif

  do k=valIn1%zmin,valIn1%zmax 
    do j=valIn1%ymin,valIn1%ymax 
      do i=valIn1%xmin,valIn1%xmax 
        valOut1%values(i,j,k) = diffK( wave%ky(j) , valIn3%values(i,j,k)) - &
                              & diffK( wave%kz(k) , valIn2%values(i,j,k))
        valOut2%values(i,j,k) = diffK( wave%kz(k) , valIn1%values(i,j,k)) - &
                              & diffK( wave%kx(i) , valIn3%values(i,j,k))
        valOut3%values(i,j,k) = diffK( wave%kx(i) , valIn2%values(i,j,k)) - &
                              & diffK( wave%ky(j) , valIn1%values(i,j,k))
      enddo
    enddo
  enddo

 end subroutine computeCurlKmain


!---------------------------------------------------------------------------
!> @details 
!> Compute Gradient in spectral space
!! @author  Antoine Vollant, LEGI
!! @param[in]    wave waves number for field valIn
!! @param[in]    valIn component of scalar field in spectral space
!! @param[out]   valOut1 first component of gradient in spectral space
!! @param[out]   valOut2 second component of gradient in spectral space
!! @param[out]   valOut3 third component of gradient in spectral space
!---------------------------------------------------------------------------
 subroutine computeGradientKmain(wave,valIn,valOut1,valOut2,valOut3)
 

  !IO
  type(complex_data_layout), intent(in)    :: valIn
  type(complex_data_layout), intent(inout) :: valOut1,valOut2,valOut3
  type(waveNumbers), intent(in)        :: wave

  !Internal
  integer                                  :: i,j,k

#ifdef DEBUGTOOLBOX
    if (.not. sameLayout(valIn,valOut1,valOut2,valOut3) ) then
      write(6,'(a)')'[ERROR] in computeGradientKmain data are not the same'
    endif
    if (.not. associated(valIn%values) .or. &
       &.not. associated(valOut1%values) .or. &
       &.not. associated(valOut2%values) .or. &
       &.not. associated(valOut3%values) .or. &
       &.not. associated(wave%kx) )  then
      write(6,'(a)')'[ERROR] in computeGradientKmain data are not initialized'
    endif
#endif

  do k=valIn%zmin,valIn%zmax 
    do j=valIn%ymin,valIn%ymax 
      do i=valIn%xmin,valIn%xmax 
        valOut1%values(i,j,k) = diffK( wave%kx(i) , valIn%values(i,j,k))
        valOut2%values(i,j,k) = diffK( wave%ky(j) , valIn%values(i,j,k))
        valOut3%values(i,j,k) = diffK( wave%kz(k) , valIn%values(i,j,k))
      enddo
    enddo
  enddo
  

 end subroutine computeGradientKmain



!---------------------------------------------------------------------------
!> @details 
!> Compute Divergence in spectral space
!! @author  Antoine Vollant, LEGI
!! @param[in]    wave waves number for field valIn
!! @param[in]    valIn1 first component of scalar field in spectral space
!! @param[in]    valIn2 second component of scalar field in spectral space
!! @param[in]    valIn3 third component of scalar field in spectral space  
!! @param[out]   valOut divergence in spectral space
!---------------------------------------------------------------------------
 subroutine computeDivergenceKmain(wave,valIn1,valIn2,valIn3,valOut)


 !IO
 type(complex_data_layout), intent(in)    :: valIn1,valIn2,valIn3
 type(complex_data_layout), intent(inout) :: valOut
 type(waveNumbers), intent(in)        :: wave

 !Internal
 integer                                  :: i,j,k

#ifdef DEBUGTOOLBOX

    if (.not. sameLayout(valIn1,valIn2,valIn3,valOut) ) then
      write(6,'(a)')'[ERROR] in computeDivergenceKmain data are not the same'
    endif
    if (.not. associated(valIn1%values) .or. &
       &.not. associated(valIn2%values) .or. &
       &.not. associated(valIn3%values) .or. &
       &.not. associated(valOut%values) .or. &
       &.not. associated(wave%kx) )  then
      write(6,'(a)')'[ERROR] in computeDivergenceKmain data are not initialized'
    endif
#endif

  do k=valIn1%zmin,valIn1%zmax 
    do j=valIn1%ymin,valIn1%ymax 
      do i=valIn1%xmin,valIn1%xmax 
        valOut%values(i,j,k) = diffK( wave%kx(i) , valIn1%values(i,j,k)) + &
                             & diffK( wave%ky(j) , valIn2%values(i,j,k)) + &
                             & diffK( wave%kz(k) , valIn3%values(i,j,k))
      enddo
    enddo
  enddo

 end subroutine computeDivergenceKmain


!---------------------------------------------------------------------------
!> @details 
!> Compute differenciation at order 2 in spectral space
!! @author  Antoine Vollant, LEGI
!! @param[in]    wave waves number
!! @param[in]    val is the value to differenciate in spectral space
!! @param[out]   valOut differenciation in spectral space
!---------------------------------------------------------------------------
 function diff2K(wave,val) result(valOut)

  !IO
  complex(WP), intent(in)  :: val
  real(WP), intent(in)     :: wave
  complex(WP)              :: valOut

  valOut = - wave * wave * val

 end function diff2K


!---------------------------------------------------------------------------
!> @details 
!> Compute differenciation at order 1 in spectral space
!! @author  Antoine Vollant, LEGI
!! @param[in]    wave waves number
!! @param[in]    val is the value to differenciate in spectral space
!! @param[out]   valOut differenciation in spectral space
!---------------------------------------------------------------------------
 function diffK(wave,val) result(valOut)
 
  !IO
  complex(WP), intent(in)  :: val
  real(WP), intent(in)     :: wave
  complex(WP)              :: valOut

  !Internal
  complex(WP)              :: ii = (0.0_WP,1.0_WP)

  valOut = ii * wave * val

 end function diffK


!---------------------------------------------------------------------------
!> @details 
!> Compute Divergence in spectral space
!! @author  Antoine Vollant, LEGI
!! @param[in]    wave waves number for field valIn
!! @param[in]    valIn array field in spectral space
!! @param[out]   valOut divergence in spectral space
!---------------------------------------------------------------------------
 subroutine computeDivergenceK1(wave,valIn,valOut)

  !IO
  type(complex_data_layout),dimension(:), intent(in) :: valIn
  type(complex_data_layout), intent(inout)           :: valOut
  type(waveNumbers), intent(in)                      :: wave
 
  call computeDivergenceKmain(wave,valIn(1),valIn(2),valIn(3),valOut)

 end subroutine computeDivergenceK1
 


!---------------------------------------------------------------------------
!> @details 
!> Compute Gradient in spectral space
!! @author  Antoine Vollant, LEGI
!! @param[in]    wave waves number for field valIn
!! @param[in]    valIn component of scalar field in spectral space
!! @param[out]   valOut array field of gradient in spectral space
!---------------------------------------------------------------------------
 subroutine computeGradientK1(wave,valIn,valOut)
 
  !IO
  type(complex_data_layout), intent(in)                :: valIn
  type(complex_data_layout),dimension(:),intent(inout) :: valOut
  type(waveNumbers), intent(in)                        :: wave

  call computeGradientKmain(wave,valIn,valOut(1),valOut(2),valOut(3))

 end subroutine computeGradientK1


!---------------------------------------------------------------------------
!> @details 
!> Compute curl in spectral space
!! @author  Antoine Vollant, LEGI
!! @param[in]    wave waves number for field valIn
!! @param[in]    valIn array field in spectral space
!! @param[out]   valOut array field of curl in spectral space
!---------------------------------------------------------------------------
 subroutine computeCurlK1(wave,valIn,valOut)

  !IO
  type(complex_data_layout),dimension(:),intent(in)    :: valIn
  type(complex_data_layout),dimension(:),intent(inout) :: valOut
  type(waveNumbers), intent(in)        :: wave

#ifdef DEBUGTOOLBOX
    if (.not. sameLayout(valIn(1),valOut(1)) .or. &
       & ( size(valOut,1) .ne. 3) .or. ( size(valIn,1) .ne. 3) ) then
      write(6,'(a)')'[ERROR] in computeCurlK1'
    endif
#endif

  call computeCurlKmain(wave,valIn(1),valIn(2),valIn(3),valOut(1),valOut(2),valOut(3))

 end subroutine computeCurlK1


!---------------------------------------------------------------------------
!> @details 
!> Compute curl in spectral space
!! @author  Antoine Vollant, LEGI
!! @param[in]    wave waves number for field valIn
!! @param[in]    valIn1 first component of scalar field in spectral space
!! @param[in]    valIn2 second component of scalar field in spectral space
!! @param[in]    valIn3 third component of scalar field in spectral space  
!! @param[out]   valOut array field of curl in spectral space
!---------------------------------------------------------------------------
 subroutine computeCurlK2(wave,valIn1,valIn2,valIn3,valOut)

  !IO
  type(complex_data_layout), intent(in)                :: valIn1,valIn2,valIn3
  type(complex_data_layout),dimension(:),intent(inout) :: valOut
  type(waveNumbers), intent(in)                        :: wave

#ifdef DEBUGTOOLBOX
    if (.not. sameLayout(valIn1,valIn2,valIn3,valOut(1)) .or. &
       & ( size(valOut,1) .ne. 3) ) then
      write(6,'(a)')'[ERROR] in computeCurlK2'
    endif
#endif

  call computeCurlKmain(wave,valIn1,valIn2,valIn3,valOut(1),valOut(2),valOut(3))

 end subroutine computeCurlK2


!---------------------------------------------------------------------------
!> @details 
!> Compute Laplacian of array field in spectral space
!! @author  Antoine Vollant, LEGI
!! @param[in]    wave waves number for field valIn
!! @param[in]    valIn array field in spectral space
!! @param[out]   valOut array field of laplacian in spectral space
!---------------------------------------------------------------------------
 subroutine computeLaplacianVecK1(wave,valIn,valOut)

  !IO
  type(complex_data_layout),dimension(:),intent(in)    :: valIn
  type(complex_data_layout),dimension(:),intent(inout) :: valOut
  type(waveNumbers), intent(in)                        :: wave
  
  call computeLaplacianVecKmain(wave,valIn(1),valIn(2),valIn(3),valOut(1),valOut(2),valOut(3))

 end subroutine computeLaplacianVecK1

!-------------------------------------------------------------------------------------------
!> Compute the nth derivative of a field with respect to the ith component
!!
!! @author Antoine Volllant
!!    @param[in]    FieldIn   = field to derivate
!!    @param[inout] FieldOut  = d^nth(FieldIn)/dxi
!!    @param[in]    numComp   = array which contain the components according to derivate
!!                  ex: d^2u/dxdy numComp=(/1,2/)
!!    @param[in]    spec_rank = number of process in cpu pool
!!    @param[in]    WN        = waves number of fieldIn (optinal)
!-------------------------------------------------------------------------------------------
function computeDerivationField_Spe(InK,numComp,FieldOut,spec_rank,WN) result(res)

  use rediscretization_tools

  !I/O data
  type(complex_data_layout),intent(in)  :: InK
  type(real_data_layout),intent(inout)  :: FieldOut
  integer,intent(in),dimension(:)       :: numComp
  type(wavenumbers),intent(in),optional :: WN
  integer,intent(in)                    :: spec_rank

  !Local data
  type(complex_data_layout)             :: tempInK
  integer                               :: i,j,k,l
  integer                               :: degDeriv,nbcpus
  complex(WP)                           :: tempK
  logical                               :: res
  type(wavenumbers)                     :: WNLoc

  res=.false.

  degDeriv=size(numComp,1)
  nbcpus = getnbcpus()
  if (.not. copyStructOnly(inK,tempInk) ) return 
  tempInK%values=Ink%values

  if (  present(WN) ) then
    do l=1,degDeriv
      if ( numComp(l) .eq. 1 ) then
        do k=tempInk%zmin,tempInk%zmax
          do j=tempInk%ymin,tempInk%ymax
            do i=tempInk%xmin,tempInk%xmax
              tempK = diffK( WN%kx(i) , tempInk%values(i,j,k))
              tempInk%values(i,j,k) = tempK
            enddo
          enddo
        enddo
      elseif ( numComp(l) .eq. 2 ) then
        do k=tempInk%zmin,tempInk%zmax
          do j=tempInk%ymin,tempInk%ymax
            do i=tempInk%xmin,tempInk%xmax
              tempK = diffK( WN%ky(j) , tempInk%values(i,j,k))
              tempInk%values(i,j,k) = tempK
            enddo
          enddo
        enddo
      elseif ( numComp(l) .eq. 3 ) then
        do k=tempInk%zmin,tempInk%zmax
          do j=tempInk%ymin,tempInk%ymax
            do i=tempInk%xmin,tempInk%xmax
              tempK = diffK( WN%kz(k) , tempInk%values(i,j,k))
              tempInk%values(i,j,k) = tempK
            enddo
          enddo
        enddo
      else
        write (6,'(a)')'[ERROR] Wrong component !'
      endif
    enddo
  else
    if (.not. initWN(WNLoc,spec_rank,FieldOut%nx,FieldOut%ny,FieldOut%nz)) return
    if (.not. computeWN(WNLoc,spec_rank,FieldOut%Lx,FieldOut%Ly,FieldOut%Lz,FieldOut%nx,&
                    &FieldOut%ny,FieldOut%nz)) return
    do l=1,degDeriv
      if ( numComp(l) .eq. 1 ) then
        do k=tempInk%zmin,tempInk%zmax
          do j=tempInk%ymin,tempInk%ymax
            do i=tempInk%xmin,tempInk%xmax
              tempK = diffK( WNLoc%kx(i) , tempInk%values(i,j,k))
              tempInk%values(i,j,k) = tempK
            enddo
          enddo
        enddo
      elseif ( numComp(l) .eq. 2 ) then
        do k=tempInk%zmin,tempInk%zmax
          do j=tempInk%ymin,tempInk%ymax
            do i=tempInk%xmin,tempInk%xmax
              tempK = diffK( WNLoc%ky(j) , tempInk%values(i,j,k))
              tempInk%values(i,j,k) = tempK
            enddo
          enddo
        enddo
      elseif ( numComp(l) .eq. 3 ) then
        do k=tempInk%zmin,tempInk%zmax
          do j=tempInk%ymin,tempInk%ymax
            do i=tempInk%xmin,tempInk%xmax
              tempK = diffK( WNLoc%kz(k) , tempInk%values(i,j,k))
              tempInk%values(i,j,k) = tempK
            enddo
          enddo
        enddo
      else
        write (6,'(a)')'[ERROR] Wrong component !'
        return
      endif
    enddo
    res = deleteWN(WNLoc)
  endif
  if ( ((FieldOut%nx/2)+1).eq.tempInk%nx .and. FieldOut%ny.eq.tempInk%ny .and.FieldOut%nz.eq.tempInk%nz) then
    call btran(tempInk,FieldOut,res)
    if (.not. res) then
      write(6,'(a)')'[ERROR] btran in computeDerivationField failed !'
      return
    endif
  else
    if (.not. transdiscretization (Ink, FieldOut, FieldOut)) then
      write(6,'(a)')'[ERROR] transdiscretization in computeDerivationField failed !'
      return
    endif
  end if
  call deletedatalayout(tempInK)
  res=.true.

end function computeDerivationField_Spe


!-------------------------------------------------------------------------------------------
!> Compute the nth derivative of a field with respect to the ith component
!!
!! @author Antoine Volllant
!!    @param[in]    FieldIn   = field to derivate
!!    @param[inout] FieldOut  = d^nth(FieldIn)/dxi
!!    @param[in]    numComp   = array which contain the components according to derivate
!!                  ex: d^2u/dxdy numComp=(/1,2/)
!!    @param[in]    spec_rank = number of process in cpu pool
!!    @param[in]    WN        = waves number of fieldIn (optinal)
!-------------------------------------------------------------------------------------------
function computeDerivationField_Phy(FieldIn,numComp,FieldOut,spec_rank,WN) result(res)

  use rediscretization_tools

  !I/O data
  type(real_data_layout),intent(in)     :: FieldIn
  type(real_data_layout),intent(inout)  :: FieldOut
  integer,intent(in),dimension(:)       :: numComp
  type(wavenumbers),intent(in),optional :: WN
  integer,intent(in)                    :: spec_rank

  !Local data
  integer                               :: i,j,k,l
  integer                               :: degDeriv,nbcpus
  type(complex_data_layout)             :: InK
  complex(WP)                           :: tempK
  logical                               :: res
  type(wavenumbers)                     :: WNLoc

  res=.false.

  degDeriv=size(numComp)
  nbcpus = getnbcpus()

  if (.not. initDataLayout("InK",InK,(FieldIn%nx/2)+1,FieldIn%ny,FieldIn%nz, &
     & FieldIn%Lx,FieldIn%Ly,FieldIn%Lz,nbcpus,spec_rank,alongZ)) then
    write(6,'(a,i0,a)')'[ERROR] initDataLayout for TYPE (COMPLEX_DATA_LAYOUT)&
                         & : not enought memory!'
    return
  endif

  call ftran(FieldIn,InK,res)
  if (.not. res) then
    write(6,'(a)')'[ERROR] ftran in computeDerivationField failed !'
    return
  endif
  if (  present(WN) ) then
    do l=1,degDeriv
      if ( numComp(l) .eq. 1 ) then
        do k=InK%zmin,InK%zmax
          do j=InK%ymin,InK%ymax
            do i=InK%xmin,InK%xmax
              tempK = diffK( WN%kx(i) , InK%values(i,j,k))
              InK%values(i,j,k) = tempK
            enddo
          enddo
        enddo
      elseif ( numComp(l) .eq. 2 ) then
        do k=InK%zmin,InK%zmax
          do j=InK%ymin,InK%ymax
            do i=InK%xmin,InK%xmax
              tempK = diffK( WN%ky(j) , InK%values(i,j,k))
              InK%values(i,j,k) = tempK
            enddo
          enddo
        enddo
      elseif ( numComp(l) .eq. 3 ) then
        do k=InK%zmin,InK%zmax
          do j=InK%ymin,InK%ymax
            do i=InK%xmin,InK%xmax
              tempK = diffK( WN%kz(k) , InK%values(i,j,k))
              InK%values(i,j,k) = tempK
            enddo
          enddo
        enddo
      else
        write (6,'(a)')'[ERROR] Wrong component !'
      endif
    enddo
  else
    if (.not. initWN(WNLoc,spec_rank,FieldIn%nx,FieldIn%ny,FieldIn%nz)) return
    if (.not. computeWN(WNLoc,spec_rank,FieldIn%Lx,FieldIn%Ly,FieldIn%Lz,FieldIn%nx,&
                    &FieldIn%ny,FieldIn%nz)) return
    do l=1,degDeriv
      if ( numComp(l) .eq. 1 ) then
        do k=InK%zmin,InK%zmax
          do j=InK%ymin,InK%ymax
            do i=InK%xmin,InK%xmax
              tempK = diffK( WNLoc%kx(i) , InK%values(i,j,k))
              InK%values(i,j,k) = tempK
            enddo
          enddo
        enddo
      elseif ( numComp(l) .eq. 2 ) then
        do k=InK%zmin,InK%zmax
          do j=InK%ymin,InK%ymax
            do i=InK%xmin,InK%xmax
              tempK = diffK( WNLoc%ky(j) , InK%values(i,j,k))
              InK%values(i,j,k) = tempK
            enddo
          enddo
        enddo
      elseif ( numComp(l) .eq. 3 ) then
        do k=InK%zmin,InK%zmax
          do j=InK%ymin,InK%ymax
            do i=InK%xmin,InK%xmax
              tempK = diffK( WNLoc%kz(k) , InK%values(i,j,k))
              InK%values(i,j,k) = tempK
            enddo
          enddo
        enddo
      else
        write (6,'(a)')'[ERROR] Wrong component !'
        return
      endif
    enddo
    res = deleteWN(WNLoc)
  endif
  if (samelayout(FieldIn,FieldOut)) then
    call btran(InK,FieldOut,res)
    if (.not. res) then
      write(6,'(a)')'[ERROR] btran in computeDerivationField failed !'
      return
    endif
  else
    if (.not. transdiscretization (Ink, FieldOut, FieldOut)) then
      write(6,'(a)')'[ERROR] transdiscretization in computeDerivationField failed !'
      return
    endif
  end if

  call deletedatalayout(InK)

  res=.true.

end function computeDerivationField_phy

!------------------------------------------------------------------------------
!> Compute Gradient of a scalar field.
!! @author Antoine Vollant, LEGI
!!    @param[in]    wave            = array wave numbers
!!    @param[in]    inField        = scalar field
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in,out]   outField       = array field of gradient 
!!    @return       success         = logical equals to true is everything is right.
!! @details
!!        This function compute the gradient of a scalar field.
!------------------------------------------------------------------------------
function computeFieldGradient_vecField(wave,inField,outField,nbcpus,spec_rank) result(success)

    ! Input/Output
    integer, intent(in) :: spec_rank,nbcpus
    type(REAL_DATA_LAYOUT), intent(in)  :: inField
    type(WaveNumbers), intent(in) :: wave 
    type(REAL_DATA_LAYOUT),dimension(:),intent(inout) :: outField
    logical :: success

    success = computeFieldGradient_scaField(wave,inField,outField(1),outField(2),outField(3),&
                                           & nbcpus,spec_rank)

end function computeFieldGradient_vecField


!------------------------------------------------------------------------------
!> Compute Gradient of a scalar field.
!! @author Antoine Vollant, LEGI
!!    @param[in]        wave            = array wave numbers
!!    @param[in]        inField         = scalar field 
!!    @param[in]        spec_rank       = mpi rank inside spectral communicator
!!    @param[in,out]    outField1       = first component of gradient 
!!    @param[in,out]    outField2       = second component of gradient
!!    @param[in,out]    outField3       = third component of gradient
!!    @return           success         = logical equals to true is everything is right.
!! @details
!!        This function compute the gradient of a scalar field. 
!------------------------------------------------------------------------------
function computeFieldGradient_scaField(wave,inField,outField1,outField2,outField3,nbcpus,spec_rank) result(success)

    ! Input/Output
    integer, intent(in) :: spec_rank,nbcpus
    type(REAL_DATA_LAYOUT), intent(in)  :: inField
    type(WaveNumbers), intent(in) :: wave 
    type(REAL_DATA_LAYOUT), intent(inout) :: outField1,outField2,outField3
    logical :: success
    ! Other Variables
    type(COMPLEX_DATA_LAYOUT) :: inFieldK,outField1K,outField2K,outField3K
    logical :: res1,res2,res3

    success = .false.

    if ( associated(outField1%values) .and. &
       & associated(outField2%values) .and. &
       & associated(outField3%values) ) then
       if (.not. samelayout(outField1,outField2,outField3,inField) ) then
         write(6,'(a,i0,a)')'[ERROR] Fields in computeFieldGradient are not the same!'
         return
       endif
    endif
    !initialise COMPLEX_DATA_LAYOUT
    if (.not. initDataLayout("inField_Spectral", & 
            & inFieldK,(inField%nx/2)+1,inField%ny,inField%nz, &
            & inField%Lx,inField%Ly,inField%Lz,nbcpus,spec_rank,alongZ)) then
        write(6,'(a,i0,a)')'[ERROR] initDataLayout for TYPE (COMPLEX_DATA_LAYOUT) : not enought memory!'
        return
    endif
    if (.not. copystructonly(inFieldK,outField1K) .or. &
       &.not. copystructonly(inFieldK,outField2K) .or. &
       &.not. copystructonly(inFieldK,outField3K) ) then 
        write(6,'(a,i0,a)')'[ERROR] initDataLayout for TYPE (COMPLEX_DATA_LAYOUT) : not enought memory!'
        return
    endif
    
    !compute FFT
    call ftran(inField,inFieldK,res1)
    if ( .not. res1 ) then
       write(6,'(a,i0,a)')'[ERROR] FFT in computeFieldGradient : failed!'
       return
    endif

    !compute Gradient in spectral space
    call computeGradientK(wave,inFieldK,outField1K,outField2K,outField3K)

    !compute FFT inverse
    call btran(outField1k,outField1,res1)
    call btran(outField2k,outField2,res2)
    call btran(outField3k,outField3,res3)
    if ( .not. res1 .or. .not. res2 .or. .not. res3 ) then
       write(6,'(a,i0,a)')'[ERROR] FFT inverse in computeFieldGradient : failed!'
      return
    endif

    !deallocate COMPLEX_DATA_LAYOUT
    call deleteDataLayout(inFieldK)
    call deleteDataLayout(outField1K)
    call deleteDataLayout(outField2K)
    call deleteDataLayout(outField3K)

    success = .true.

end function computeFieldGradient_scaField


!------------------------------------------------------------------------------
!> Compute divergence of an array field.
!! @author Antoine Vollant, LEGI
!!    @param[in]    wave            = array wave numbers
!!    @param[in]    inField        = array field
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in,out]   outField       = scalar field which contains the divergence
!!    @return       success         = logical equals to true is everything is right.
!! @details
!!        This function compute the divergence of an array field. 
!------------------------------------------------------------------------------
function computeFieldDivergence_vecField(wave,inField,outField,nbcpus,spec_rank) result(success)


    ! Input/Output
    integer, intent(in) :: spec_rank,nbcpus
    type(REAL_DATA_LAYOUT),dimension(:),intent(in)  :: inField
    type(WaveNumbers), intent(in) :: wave 
    type(REAL_DATA_LAYOUT), intent(inout) :: outField
    logical :: success


    success = computeFieldDivergence_scaField(wave,inField(1),inField(2),inField(3),&
                                             &outField,nbcpus,spec_rank)


end function computeFieldDivergence_vecField




!------------------------------------------------------------------------------
!> Compute divergence of an array field.
!! @author Antoine Vollant, LEGI
!!    @param[in]    wave            = array wave numbers
!!    @param[in]    inField1        = scalar field (first component)
!!    @param[in]    inField2        = scalar field (second component)
!!    @param[in]    inField3        = scalar field (third component)
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in,out]   outField       = scalar field which contains the divergence
!!    @return       success         = logical equals to true is everything is right.
!! @details
!!        This function compute the divergence of an array field.
!------------------------------------------------------------------------------
function computeFieldDivergence_scaField(wave,inField1,inField2,inField3,outField,nbcpus,spec_rank) result(success)


    ! Input/Output
    integer, intent(in) :: spec_rank,nbcpus
    type(REAL_DATA_LAYOUT), intent(in)  :: inField1,inField2,inField3
    type(WaveNumbers), intent(in) :: wave
    type(REAL_DATA_LAYOUT), intent(inout) :: outField
    logical :: success
    ! Other Variables
    type(COMPLEX_DATA_LAYOUT) :: inField1K,inField2K,inField3K,outFieldK
    logical :: res1,res2,res3

    success = .false.
    !check if inFields are conforms structure
    if (.not. samelayout(inField1,inField2,inField3) ) then
      write(6,'(a,i0,a)')'[ERROR] Fields in computeFieldDivergence are not the same!'
      return
    endif

    !initialise COMPLEX_DATA_LAYOUT
    if (.not. initDataLayout("inField1_Spectral", & 
            & inField1K,(inField1%nx/2)+1,inField1%ny,inField1%nz, &
            & inField1%Lx,inField1%Ly,inField1%Lz,nbcpus,spec_rank,alongZ)) then
        write(6,'(a,i0,a)')'[ERROR] initDataLayout for TYPE (COMPLEX_DATA_LAYOUT) : not enought memory!'
        return
    endif
    if (.not. copyStructOnly(inField1K,inField2K)) return
    if (.not. copyStructOnly(inField1K,inField3K)) return
    if (.not. copyStructOnly(inField1K,outFieldK)) return
    
    !compute FFT
    call ftran(inField1,inField1K,res1)
    call ftran(inField2,inField2K,res2)
    call ftran(inField3,inField3K,res3)
    if ( .not. res1 .or. .not. res2 .or. .not. res3 ) then
       write(6,'(a,i0,a)')'[ERROR] FFT in computeFieldDivregence : failed!'
       return
    endif

    !compute divergence in spectral space
    call computeDivergenceK(wave,inField1K,inField2K,inField3K,outFieldK)

    !compute FFT inverse
    call btran(outFieldk,outField,res1)
    if ( .not. res1 ) then
       write(6,'(a,i0,a)')'[ERROR] FFT inverse in computeFieldDivergence : failed!'
      return
    endif

    !deallocate COMPLEX_DATA_LAYOUT
    call deleteDataLayout(inField1K)
    call deleteDataLayout(inField2K)
    call deleteDataLayout(inField3K)
    call deleteDataLayout(outFieldK)

    success = .true.

end function computeFieldDivergence_scaField


!-------------------------------------------------------------------------------------------
!> Compute sum[(d(field)/dx_i)**2]
!!
!! @author Jean-Baptiste Lagaert
!!    @param[in]    FieldIn   = field to derivate
!!    @param[inout] FieldOut  = sum(d(FieldIn)/dxi) for i = 1,3
!!    @param[in]    spec_rank = number of process in cpu pool
!!    @param[in]    WN        = waves number of fieldIn (optinal)
!-------------------------------------------------------------------------------------------
function computeSumSquareFirstDerivateField(FieldIn,FieldOut,spec_rank,WN) result(res)

  use rediscretization_tools

  !I/O data
  type(real_data_layout),intent(in)     :: FieldIn
  type(real_data_layout),intent(inout)  :: FieldOut
  type(wavenumbers),intent(in),optional :: WN
  integer,intent(in)                    :: spec_rank

  !Local data
  integer                               :: i,j,k
  integer                               :: nbcpus
  type(complex_data_layout)             :: InKX, InKY, InKZ
  type(real_data_layout)                :: temp
  complex(WP)                           :: tempK
  logical                               :: res
  type(wavenumbers)                     :: WNLoc

  res=.false.

  nbcpus = getnbcpus()

  if (.not. initDataLayout("InKX",InKX,(FieldIn%nx/2)+1,FieldIn%ny,FieldIn%nz, &
     & FieldIn%Lx,FieldIn%Ly,FieldIn%Lz,nbcpus,spec_rank,alongZ)) then
    write(6,'(a,i0,a)')'[ERROR] initDataLayout for TYPE (COMPLEX_DATA_LAYOUT)&
                         & : not enought memory!'
    return
  endif
  if (.not. initDataLayout("InKY",InKY,(FieldIn%nx/2)+1,FieldIn%ny,FieldIn%nz, &
     & FieldIn%Lx,FieldIn%Ly,FieldIn%Lz,nbcpus,spec_rank,alongZ)) then
    write(6,'(a,i0,a)')'[ERROR] initDataLayout for TYPE (COMPLEX_DATA_LAYOUT)&
                         & : not enought memory!'
    return
  endif
  if (.not. initDataLayout("InKZ",InKZ,(FieldIn%nx/2)+1,FieldIn%ny,FieldIn%nz, &
     & FieldIn%Lx,FieldIn%Ly,FieldIn%Lz,nbcpus,spec_rank,alongZ)) then
    write(6,'(a,i0,a)')'[ERROR] initDataLayout for TYPE (COMPLEX_DATA_LAYOUT)&
                         & : not enought memory!'
    return
  endif


  call ftran(FieldIn,InKX,res)
  if (.not. res) then
    write(6,'(a)')'[ERROR] ftran in computeDerivationField failed !'
    return
  endif
  if (  present(WN) ) then
    do k=InKX%zmin,InKX%zmax
      do j=InKX%ymin,InKX%ymax
        do i=InKX%xmin,InKX%xmax
          tempK              = diffK( WN%kx(i) , InKX%values(i,j,k))
          InKY%values(i,j,k) = diffK( WN%ky(j) , InKX%values(i,j,k))
          InKZ%values(i,j,k) = diffK( WN%kz(k) , InKX%values(i,j,k))
          InKX%values(i,j,k) = tempK
        end do
      end do
    end do
  else
    if (.not. initWN(WNLoc,spec_rank,FieldIn%nx,FieldIn%ny,FieldIn%nz)) return
    if (.not. computeWN(WNLoc,spec_rank,FieldIn%Lx,FieldIn%Ly,FieldIn%Lz,FieldIn%nx,&
                    &FieldIn%ny,FieldIn%nz)) return
    do k=InKX%zmin,InKX%zmax
      do j=InKX%ymin,InKX%ymax
        do i=InKX%xmin,InKX%xmax
          tempK              = diffK( WNLoc%kx(i) , InKX%values(i,j,k))
          InKY%values(i,j,k) = diffK( WNLoc%ky(j) , InKX%values(i,j,k))
          InKZ%values(i,j,k) = diffK( WNLoc%kz(k) , InKX%values(i,j,k))
          InKX%values(i,j,k) = tempK
        end do
      end do
    end do
    res = deleteWN(WNLoc)
  endif
  !if (samelayout(FieldIn,FieldOut)) then
  call btran(InKX,FieldOut,res)
  if (.not. res) then
    write(6,'(a)')'[ERROR] btran in computeDerivationField failed !'
    return
  endif
  !else
  !  if (.not. transdiscretization (InK, FieldOut, FieldOut)) then
  !    write(6,'(a)')'[ERROR] transdiscretization in computeDerivationField failed !'
  !    return
  !  endif
  !end if
  call deletedatalayout(InKX)
  if (.not. copyStructOnly(FieldOut,temp)) then
    write(6,'(a,i0,a)')'[ERROR] copyStructOnly : not enought memory!'
    return
  endif
  call btran(InKY,temp,res)
  if (.not. res) then
    write(6,'(a)')'[ERROR] btran in computeSumSqaureFirstDeriv failed !'
    return
  endif
  call deletedatalayout(InKY)
  FieldOut%values = (FieldOut%values**2) + (temp%values**2)
  call btran(InKZ,temp,res)
  if (.not. res) then
    write(6,'(a)')'[ERROR] btran in computeSumSqaureFirstDeriv failed !'
    return
  endif
  call deletedatalayout(InKZ)
  FieldOut%values = FieldOut%values + (temp%values**2)
  call deletedatalayout(temp)

  res=.true.

end function computeSumSquareFirstDerivateField


!------------------------------------------------------------------------------
!> Compute laplacian of an scalar field.
!! @author Antoine Vollant, LEGI
!!    @param[in]    wave            = array wave numbers
!!    @param[in]    inField        = scalar field (first component)
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in,out]   outField       = field out which contains the Laplacian 
!!    @return       success         = logical equals to true is everything is right.
!! @details
!!        This function compute the laplacian of an array field. 
!------------------------------------------------------------------------------
function computeFieldLaplacianScalar(wave,inField,outField,nbcpus,spec_rank) result(success)


    ! Input/Output
    integer, intent(in) :: spec_rank,nbcpus
    type(REAL_DATA_LAYOUT), intent(in)  :: inField
    type(WaveNumbers), intent(in) :: wave 
    type(REAL_DATA_LAYOUT), intent(inout) :: outField
    logical :: success
    ! Other Variables
    type(COMPLEX_DATA_LAYOUT) :: inFieldK,outFieldK
    logical :: res

    success = .false.

    if ( associated(outField%values) ) then
       if (.not. samelayout(outField,inField) ) then
         write(6,'(a,i0,a)')'[ERROR] Fields in computeFieldLaplacian are not the same!'
         return
       endif
    else
      !initialise Real_DATA_LAYOUT as it hasn't been already done !!
      if (.not. associated(outField%values)) then
        if (.not. copyStructOnly(inField,outField,"Laplacian") ) then
          write(6,'(a,i0,a)')'[ERROR] copyStructOnly for TYPE (REAL_DATA_LAYOUT) : not enought memory!'
          return
        endif
      endif
      if (.not. samelayout(outField,inField) ) then
        write(6,'(a,i0,a)')'[ERROR] Fields in computeFieldLaplacian are not the same!'
        return
      endif
    endif
  
    !initialise COMPLEX_DATA_LAYOUT
    if (.not. initDataLayout("inField_Spectral", & 
            & inFieldK,(inField%nx/2)+1,inField%ny,inField%nz, &
            & inField%Lx,inField%Ly,inField%Lz,nbcpus,spec_rank,alongZ)) then
        write(6,'(a,i0,a)')'[ERROR] initDataLayout for TYPE (COMPLEX_DATA_LAYOUT) : not enought memory!'
        return
    endif
    if (.not. copyStructOnly(inFieldK,outFieldK)) then 
        write(6,'(a,i0,a)')'[ERROR] initDataLayout for TYPE (COMPLEX_DATA_LAYOUT) : not enought memory!'
        return
    endif
    
    !compute FFT
    call ftran(inField,inFieldK,res)
    if ( .not. res ) then
       write(6,'(a,i0,a)')'[ERROR] FFT in computeFieldLaplacianArray : failed!'
       return
    endif

    !compute curl in spectral space
    call computeLaplacianK(wave,inFieldK,outFieldK)

    !compute FFT inverse
    call btran(outFieldk,outField,res)
    if ( .not. res ) then
       write(6,'(a,i0,a)')'[ERROR] FFT inverse in computeFieldLaplacianArray : failed!'
      return
    endif

    !deallocate COMPLEX_DATA_LAYOUT
    call deleteDataLayout(inFieldK)
    call deleteDataLayout(outFieldK)

    success = .true.

end function computeFieldLaplacianScalar


!------------------------------------------------------------------------------
!> Compute laplacian of an array field.
!! @author Antoine Vollant, LEGI
!!    @param[in]     wave           = array wave numbers
!!    @param[in]     inField        = array field (first component)
!!    @param[in,out] outField       = array field of laplacian
!!    @param[in]     spec_rank      = mpi rank inside spectral communicator
!!    @return        success        = logical equals to true is everything is right.
!! @details
!!        This function compute the laplacian of an array field. 
!------------------------------------------------------------------------------
function computeFieldLaplacianArray_vecField(wave,inField,outField,nbcpus,&
                                   &spec_rank) result(success)


    integer, intent(in) :: spec_rank,nbcpus
    type(REAL_DATA_LAYOUT),dimension(:),intent(in)  :: inField
    type(WaveNumbers), intent(in) :: wave 
    type(REAL_DATA_LAYOUT),dimension(:),intent(inout) :: outField
    logical :: success


    success = computeFieldLaplacianArray_scaField(wave,inField(1),inField(2),inField(3),&
                                                 &outField(1),outField(2),outField(3),nbcpus,&
                                                 &spec_rank)


end function computeFieldLaplacianArray_vecField



!------------------------------------------------------------------------------
!> Compute laplacian of an array field.
!! @author Antoine Vollant, LEGI
!!    @param[in]    wave            = array wave numbers
!!    @param[in]    inField1        = scalar field (first component)
!!    @param[in]    inField2        = scalar field (second component)
!!    @param[in]    inField3        = scalar field (third component)
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in,out]   outField1       = first component of laplacian
!!    @param[in,out]   outField2       = second component of laplacian
!!    @param[in,out]   outField3       = third component of laplacian
!!    @return       success         = logical equals to true is everything is right.
!! @details
!!        This function compute the laplacian of an array field. 
!------------------------------------------------------------------------------
function computeFieldLaplacianArray_scaField(wave,inField1,inField2,inField3,&
                                   &outField1,outField2,outField3,nbcpus,&
                                   &spec_rank) result(success)


    ! Input/Output
    integer, intent(in) :: spec_rank,nbcpus
    type(REAL_DATA_LAYOUT), intent(in)  :: inField1,inField2,inField3
    type(WaveNumbers), intent(in) :: wave 
    type(REAL_DATA_LAYOUT), intent(inout) :: outField1,outField2,outField3
    logical :: success
    ! Other Variables
    type(COMPLEX_DATA_LAYOUT) :: inField1K,inField2K,inField3K,outField1K,outField2K,outField3K
    logical :: res1,res2,res3

    success = .false.
    !check if inFields are conforms structure
    if (.not. samelayout(inField1,inField2,inField3) ) then
      write(6,'(a,i0,a)')'[ERROR] Fields in computeFieldLaplacian are not the same!'
      return
    endif

    if ( associated(outField1%values) .and. &
       & associated(outField2%values) .and. &
       & associated(outField3%values) ) then
       if (.not. samelayout(outField1,outField2,outField3,inField1) ) then
         write(6,'(a,i0,a)')'[ERROR] Fields in computeFieldLaplacian are not the same!'
         return
       endif
    else
      !initialise Real_DATA_LAYOUT as it hasn't been already done !!
      if (.not. associated(outField1%values)) then
        if (.not. copyStructOnly(inField1,outField1,"rot1") ) then
          write(6,'(a,i0,a)')'[ERROR] copyStructOnly for TYPE (REAL_DATA_LAYOUT) : not enought memory!'
          return
        endif
      endif
      if (.not. associated(outField2%values)) then
        if (.not. copyStructOnly(inField2,outField2,"rot2") ) then 
          write(6,'(a,i0,a)')'[ERROR] copyStructOnly for TYPE (REAL_DATA_LAYOUT) : not enought memory!'
          return
        endif
      endif
      if (.not. associated(outField3%values)) then
        if (.not. copyStructOnly(inField3,outField3,"rot3") ) then 
          write(6,'(a,i0,a)')'[ERROR] copyStructOnly for TYPE (REAL_DATA_LAYOUT) : not enought memory!'
          return
        endif
      endif
      if (.not. samelayout(outField1,outField2,outField3,inField1) ) then
        write(6,'(a,i0,a)')'[ERROR] Fields in computeFieldLaplacian are not the same!'
        return
      endif
    endif
  
    !initialise COMPLEX_DATA_LAYOUT
    if (.not. initDataLayout("inField1_Spectral", & 
            & inField1K,(inField1%nx/2)+1,inField1%ny,inField1%nz, &
            & inField1%Lx,inField1%Ly,inField1%Lz,nbcpus,spec_rank,alongZ)) then
        write(6,'(a,i0,a)')'[ERROR] initDataLayout for TYPE (COMPLEX_DATA_LAYOUT) : not enought memory!'
        return
    endif
    if (.not. initDataLayout("inField2_Spectral", & 
            & inField2K,(inField2%nx)/2+1,inField2%ny,inField2%nz, &
            & inField2%Lx,inField2%Ly,inField2%Lz,nbcpus,spec_rank,alongZ)) then
        write(6,'(a,i0,a)')'[ERROR] initDataLayout for TYPE (COMPLEX_DATA_LAYOUT) : not enought memory!'
       return 
    endif
    if (.not. initDataLayout("inField3_Spectral", & 
            & inField3K,(inField3%nx)/2+1,inField3%ny,inField3%nz, &
            & inField3%Lx,inField3%Ly,inField3%Lz,nbcpus,spec_rank,alongZ)) then
        write(6,'(a,i0,a)')'[ERROR] initDataLayout for TYPE (COMPLEX_DATA_LAYOUT) : not enought memory!'
        return
    endif
    if (.not. initDataLayout("outField1_Spectral", & 
            & outField1K,(inField1%nx)/2+1,inField1%ny,inField1%nz, &
            & inField1%Lx,inField1%Ly,inField1%Lz,nbcpus,spec_rank,alongZ)) then
        write(6,'(a,i0,a)')'[ERROR] initDataLayout for TYPE (COMPLEX_DATA_LAYOUT) : not enought memory!'
        return
    endif
    if (.not. initDataLayout("outField2_Spectral", & 
            & outField2K,(inField2%nx)/2+1,inField2%ny,inField2%nz, &
            & inField2%Lx,inField2%Ly,inField2%Lz,nbcpus,spec_rank,alongZ)) then
        write(6,'(a,i0,a)')'[ERROR] initDataLayout for TYPE (COMPLEX_DATA_LAYOUT) : not enought memory!'
       return 
    endif
    if (.not. initDataLayout("outField3_Spectral", & 
            & outField3K,(inField3%nx)/2+1,inField3%ny,inField3%nz, &
            & inField3%Lx,inField3%Ly,inField3%Lz,nbcpus,spec_rank,alongZ)) then
        write(6,'(a,i0,a)')'[ERROR] initDataLayout for TYPE (COMPLEX_DATA_LAYOUT) : not enought memory!'
        return
    endif
    
    !compute FFT
    call ftran(inField1,inField1K,res1)
    call ftran(inField2,inField2K,res2)
    call ftran(inField3,inField3K,res3)
    if ( .not. res1 .or. .not. res2 .or. .not. res3 ) then
       write(6,'(a,i0,a)')'[ERROR] FFT in computeFieldLaplacianArray : failed!'
       return
    endif

    !compute curl in spectral space
    call computeLaplacianK(wave,inField1K,inField2K,inField3K,outField1K,outField2K,outField3K)

    !compute FFT inverse
    call btran(outField1k,outField1,res1)
    call btran(outField2k,outField2,res2)
    call btran(outField3k,outField3,res3)
    if ( .not. res1 .or. .not. res2 .or. .not. res3 ) then
       write(6,'(a,i0,a)')'[ERROR] FFT inverse in computeFieldLaplacianArray : failed!'
      return
    endif

    !deallocate COMPLEX_DATA_LAYOUT
    call deleteDataLayout(inField1K)
    call deleteDataLayout(inField2K)
    call deleteDataLayout(inField3K)
    call deleteDataLayout(outField1K)
    call deleteDataLayout(outField2K)
    call deleteDataLayout(outField3K)

    success = .true.

end function computeFieldLaplacianArray_scaField



!------------------------------------------------------------------------------
!> Compute Curl of an array field.
!! @author Antoine Vollant, LEGI
!!    @param[in]    wave            = array wave numbers
!!    @param[in]    inField        = scalar field (first component)
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in,out] outField       = first component of curl
!!    @return       success         = logical equals to true is everything is right.
!! @details
!!        This function compute the curl of an array field. 
!------------------------------------------------------------------------------
function computeFieldCurl_vecField(wave,inField,outField,nbcpus,spec_rank) result(success)



    type(REAL_DATA_LAYOUT),dimension(:),intent(in)    :: inField
    type(REAL_DATA_LAYOUT),dimension(:),intent(inout) :: outField
    integer, intent(in) :: spec_rank,nbcpus
    type(WaveNumbers), intent(in) :: wave 
    logical :: success

     success = computeFieldCurl_scaField(wave,inField(1),inField(2),inField(3),&
                                        & outField(1),outField(2),outField(3),nbcpus,&
                                        & spec_rank)

end function computeFieldCurl_vecField




!------------------------------------------------------------------------------
!> Compute Curl of an array field.
!! @author Antoine Vollant, LEGI
!!    @param[in]    wave            = array wave numbers
!!    @param[in]    inField1        = scalar field (first component)
!!    @param[in]    inField2        = scalar field (second component)
!!    @param[in]    inField3        = scalar field (third component)
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in,out]   outField1       = first component of curl
!!    @param[in,out]   outField2       = second component of curl
!!    @param[in,out]   outField3       = third component of curl
!!    @return       success         = logical equals to true is everything is right.
!! @details
!!        This function compute the curl of an array field. 
!------------------------------------------------------------------------------
function computeFieldCurl_scaField(wave,inField1,inField2,inField3,outField1,outField2,outField3,&
                         & nbcpus,spec_rank) result(success)


    ! Input/Output
    integer, intent(in) :: spec_rank,nbcpus
    type(REAL_DATA_LAYOUT), intent(in)  :: inField1,inField2,inField3
    type(WaveNumbers), intent(in) :: wave 
    type(REAL_DATA_LAYOUT), intent(inout) :: outField1,outField2,outField3
    logical :: success
    ! Other Variables
    type(COMPLEX_DATA_LAYOUT) :: inField1K,inField2K,inField3K,outField1K,outField2K,outField3K
    logical :: res1,res2,res3

    success = .false.
    !check if inFields are conforms structure
    if (.not. samelayout(inField1,inField2,inField3) ) then
      write(6,'(a,i0,a)')'[ERROR] Fields in computeFieldCurl are not the same!'
      return
    endif

    if ( associated(outField1%values) .and. &
       & associated(outField2%values) .and. &
       & associated(outField3%values) ) then
       if (.not. samelayout(outField1,outField2,outField3,inField1) ) then
         write(6,'(a,i0,a)')'[ERROR] Fields in computeFieldCurl are not the same!'
         return
       endif
    else
      !initialise Real_DATA_LAYOUT as it hasn't been already done !!
      if (.not. associated(outField1%values)) then
        if (.not. copyStructOnly(inField1,outField1,"rot1") ) then
          write(6,'(a,i0,a)')'[ERROR] copyStructOnly for TYPE (REAL_DATA_LAYOUT) : not enought memory!'
          return
        endif
      endif
      if (.not. associated(outField2%values)) then
        if (.not. copyStructOnly(inField2,outField2,"rot2") ) then 
          write(6,'(a,i0,a)')'[ERROR] copyStructOnly for TYPE (REAL_DATA_LAYOUT) : not enought memory!'
          return
        endif
      endif
      if (.not. associated(outField3%values)) then
        if (.not. copyStructOnly(inField3,outField3,"rot3") ) then 
          write(6,'(a,i0,a)')'[ERROR] copyStructOnly for TYPE (REAL_DATA_LAYOUT) : not enought memory!'
          return
        endif
      endif
    endif
  
    !initialise COMPLEX_DATA_LAYOUT
    if (.not. initDataLayout("inField1_Spectral", & 
            & inField1K,(inField1%nx/2)+1,inField1%ny,inField1%nz, &
            & inField1%Lx,inField1%Ly,inField1%Lz,nbcpus,spec_rank,alongZ)) then
        write(6,'(a,i0,a)')'[ERROR] initDataLayout for TYPE (COMPLEX_DATA_LAYOUT) : not enought memory!'
        return
    endif
    if (.not. initDataLayout("inField2_Spectral", & 
            & inField2K,(inField2%nx)/2+1,inField2%ny,inField2%nz, &
            & inField2%Lx,inField2%Ly,inField2%Lz,nbcpus,spec_rank,alongZ)) then
        write(6,'(a,i0,a)')'[ERROR] initDataLayout for TYPE (COMPLEX_DATA_LAYOUT) : not enought memory!'
       return 
    endif
    if (.not. initDataLayout("inField3_Spectral", & 
            & inField3K,(inField3%nx)/2+1,inField3%ny,inField3%nz, &
            & inField3%Lx,inField3%Ly,inField3%Lz,nbcpus,spec_rank,alongZ)) then
        write(6,'(a,i0,a)')'[ERROR] initDataLayout for TYPE (COMPLEX_DATA_LAYOUT) : not enought memory!'
        return
    endif
    if (.not. initDataLayout("outField1_Spectral", & 
            & outField1K,(inField1%nx)/2+1,inField1%ny,inField1%nz, &
            & inField1%Lx,inField1%Ly,inField1%Lz,nbcpus,spec_rank,alongZ)) then
        write(6,'(a,i0,a)')'[ERROR] initDataLayout for TYPE (COMPLEX_DATA_LAYOUT) : not enought memory!'
        return
    endif
    if (.not. initDataLayout("outField2_Spectral", & 
            & outField2K,(inField2%nx)/2+1,inField2%ny,inField2%nz, &
            & inField2%Lx,inField2%Ly,inField2%Lz,nbcpus,spec_rank,alongZ)) then
        write(6,'(a,i0,a)')'[ERROR] initDataLayout for TYPE (COMPLEX_DATA_LAYOUT) : not enought memory!'
       return 
    endif
    if (.not. initDataLayout("outField3_Spectral", & 
            & outField3K,(inField3%nx)/2+1,inField3%ny,inField3%nz, &
            & inField3%Lx,inField3%Ly,inField3%Lz,nbcpus,spec_rank,alongZ)) then
        write(6,'(a,i0,a)')'[ERROR] initDataLayout for TYPE (COMPLEX_DATA_LAYOUT) : not enought memory!'
        return
    endif
    
    !compute FFT
    call ftran(inField1,inField1K,res1)
    call ftran(inField2,inField2K,res2)
    call ftran(inField3,inField3K,res3)
    if ( .not. res1 .or. .not. res2 .or. .not. res3 ) then
       write(6,'(a,i0,a)')'[ERROR] FFT in computeFieldCurl : failed!'
       return
    endif

    !compute curl in spectral space
    call computeCurlK(wave,inField1K,inField2K,inField3K,outField1K,outField2K,outField3K)

    !compute FFT inverse
    call btran(outField1k,outField1,res1)
    call btran(outField2k,outField2,res2)
    call btran(outField3k,outField3,res3)
    if ( .not. res1 .or. .not. res2 .or. .not. res3 ) then
       write(6,'(a,i0,a)')'[ERROR] FFT inverse in computeFieldCurl : failed!'
      return
    endif

    !deallocate COMPLEX_DATA_LAYOUT
    call deleteDataLayout(inField1K)
    call deleteDataLayout(inField2K)
    call deleteDataLayout(inField3K)
    call deleteDataLayout(outField1K)
    call deleteDataLayout(outField2K)
    call deleteDataLayout(outField3K)

    success = .true.

end function computeFieldCurl_scaField



end module differential_tools
!> @}
