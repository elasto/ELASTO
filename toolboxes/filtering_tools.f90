!> @addtogroup toolbox 
!! @{
!------------------------------------------------------------------------------
!
! MODULE: filtering_tools 
!
!> @author
!> Antoine Vollant
!
! DESCRIPTION: 
!> 
!------------------------------------------------------------------------------


module filtering_tools 

 use precision_tools 
 use wavenumber_tools 
 use transforms_tools 
 use stat_tools
 use parallel_tools
 use datalayout 

 implicit none

 public computeFilter

 interface computeDelta 
   !module procedure computeDelta_anisotrope
   module procedure computeDelta_isotrope
 end interface computeDelta

 interface computeFilter 
   module procedure computeFilter_KtoK
   module procedure computeFilter_KtoP
   module procedure computeFilter_PtoP
   module procedure computeFilter_PtoK
 end interface computeFilter

contains

!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> Compute the filtered quantity of inFieldK
!> @param [in] wave is the waves number of inFieldK 
!> @param [in] filter is the type of filter (see comments next to the case) 
!> @param [in] delta is the cuting lenght 
!> @param [in] inFieldK is the field in Fourier space which is going to be filtered 
!> @param [inout] outFieldK is the filtered field in Fourier space 
!> @return  delta
!------------------------------------------------------------------------------
 subroutine computeFilter_PtoP(wave,delta,inField,outField,filter,Exk)

  !I/O 
  type(real_data_layout),intent(in)             :: inField
  type(real_data_layout),intent(inout)          :: outField
  type(waveNumbers),intent(in)                  :: wave 
  integer,intent(in)                            :: filter  
  real(WP),intent(in)                           :: delta  
  type(complex_data_layout),intent(in),optional :: ExK

  !Local data
  type(complex_data_layout)          :: tempK1
  type(complex_data_layout)          :: tempK2
  logical :: res
  integer :: nbcpus,spec_rank

  if ( present(Exk) ) then
    res = copyStructOnly(Exk,tempK1)  
    res = copyStructOnly(Exk,tempK2)  
  else
    nbcpus    = getnbcpus()
    spec_rank = getnbmycpu()
    res = initDataLayout("tempK1", tempK1,(inField%nx/2)+1,inField%ny,inField%nz, &
                       & inField%Lx,inField%Ly,inField%Lz,nbcpus,spec_rank,alongZ)
    if (.not. res) write (6,'(a)') '[ERROR] initDataLayout in computeFilter_PtoP failed' 
    res = copyStructOnly(tempK1,tempK2)  
    if (.not. res) write (6,'(a)') '[ERROR] copyStructOnly in computeFilter_PtoP failed ' 
  endif
  call ftran(inField,tempK1,res)
  if (.not. res) write (6,'(a)') '[ERROR] ftran in computeFilter_PtoP failed'
  call computeFilter_KtoK(wave,delta,tempK1,tempK2,filter) 
  call btran(tempK2,outField,res)
  if (.not. res) write (6,'(a)') '[ERROR] btran in computeFilter_PtoP failed'
 
  call deletedatalayout(tempK1)
  call deletedatalayout(tempK2)

 end subroutine computeFilter_PtoP

!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> Compute the filtered quantity of inFieldK
!> @param [in] wave is the waves number of inFieldK 
!> @param [in] filter is the type of filter (see comments next to the case) 
!> @param [in] delta is the cuting lenght 
!> @param [in] inFieldK is the field in Fourier space which is going to be filtered 
!> @param [inout] outFieldK is the filtered field in Fourier space 
!> @return  delta
!------------------------------------------------------------------------------
 subroutine computeFilter_KtoK(wave,delta,inFieldK,outFieldK,filter)

  !I/O 
  type(complex_data_layout),intent(in)       :: inFieldK
  type(waveNumbers),intent(in)               :: wave 
  type(complex_data_layout),intent(inout)    :: outFieldK
  integer,intent(in)                         :: filter  
  real(WP),intent(in)                        :: delta  

  !Local data
  integer :: i,j,k

  do k=inFieldK%zmin,inFieldK%zmax
    do j=inFieldK%ymin,inFieldK%ymax
      do i=inFieldK%xmin,inFieldK%xmax
        select case(filter)
          case(1)  !Cut-off      
            outFieldK%values(i,j,k) = inFieldK%values(i,j,k) * &
            & computeCutOffFilterK(wave%kx(i),wave%ky(j),wave%kz(k),delta)
          case(2)  !Box
            outFieldK%values(i,j,k) = inFieldK%values(i,j,k) * &
            & computeBoxFilterK(wave%kx(i),wave%ky(j),wave%kz(k),delta)
          case(3)  !Gaussian
            outFieldK%values(i,j,k) = inFieldK%values(i,j,k) * &
            & computeGaussianFilterK(wave%kx(i),wave%ky(j),wave%kz(k),delta)
          case(4)  !Cut-off + Box
            outFieldK%values(i,j,k) = inFieldK%values(i,j,k) * &
            & computeCutOffFilterK(wave%kx(i),wave%ky(j),wave%kz(k),delta) * &
            & computeBoxFilterK(wave%kx(i),wave%ky(j),wave%kz(k),delta)
          case default
            outFieldK%values(i,j,k) = inFieldK%values(i,j,k) 
            !write(6,'(a)')'[WARNING] computeFilter : Filter unknown !'
         end select
      enddo
    enddo
  enddo

 end subroutine computeFilter_KtoK

!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> Compute the filtered quantity of inFieldK
!> @param [in] wave is the waves number of inFieldK 
!> @param [in] filter is the type of filter (see comments next to the case) 
!> @param [in] delta is the cuting lenght 
!> @param [in] inFieldK is the field in Fourier space which is going to be filtered 
!> @param [inout] outFieldK is the filtered field in Fourier space 
!> @return  delta
!------------------------------------------------------------------------------
 subroutine computeFilter_PtoK(wave,delta,inField,outFieldK,filter,Exk)

  !I/O 
  type(real_data_layout),intent(in)             :: inField
  type(complex_data_layout),intent(inout)          :: outFieldK
  type(waveNumbers),intent(in)                  :: wave 
  integer,intent(in)                            :: filter  
  real(WP),intent(in)                           :: delta  
  type(complex_data_layout),intent(in),optional :: ExK

  !Local data
  type(complex_data_layout)          :: tempK
  logical :: res
  integer :: nbcpus,spec_rank

  if ( present(Exk) ) then
    res = copyStructOnly(Exk,tempK)  
  else
    nbcpus    = getnbcpus()
    spec_rank = getnbmycpu()
    res = initDataLayout("tempK", tempK,(inField%nx/2)+1,inField%ny,inField%nz, &
                       & inField%Lx,inField%Ly,inField%Lz,nbcpus,spec_rank,alongZ)
    if (.not. res) write (6,'(a)') '[ERROR] initDataLayout in computeFilter_PtoP failed' 
  endif
  call ftran(inField,tempK,res)
  if (.not. res) write (6,'(a)') '[ERROR] ftran in computeFilter_PtoP failed'
  call computeFilter_KtoK(wave,delta,tempK,outFieldK,filter) 
 
  call deletedatalayout(tempK)

 end subroutine computeFilter_PtoK

!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> Compute the filtered quantity of inFieldK
!> @param [in] wave is the waves number of inFieldK 
!> @param [in] filter is the type of filter (see comments next to the case) 
!> @param [in] delta is the cuting lenght 
!> @param [in] inFieldK is the field in Fourier space which is going to be filtered 
!> @param [inout] outFieldK is the filtered field in Fourier space 
!> @return  delta
!------------------------------------------------------------------------------
 subroutine computeFilter_KtoP(wave,delta,inFieldK,outField,filter)

  !I/O 
  type(complex_data_layout),intent(in)          :: inFieldK
  type(real_data_layout),intent(inout)          :: outField
  type(waveNumbers),intent(in)                  :: wave 
  integer,intent(in)                            :: filter  
  real(WP),intent(in)                           :: delta  

  !Local data
  type(complex_data_layout)          :: tempK
  logical :: res

  res = copyStructOnly(inFieldK,tempK)  
  if (.not. res) write (6,'(a)') '[ERROR] copyStructOnly in computeFilter_KtoP failed ' 
  call computeFilter_KtoK(wave,delta,inFieldK,tempK,filter) 
  call btran(tempK,outField,res)
  if (.not. res) write (6,'(a)') '[ERROR] btran in computeFilter_KtoP failed'
  call deletedatalayout(tempK)

 end subroutine computeFilter_KtoP

!---------------------------------------------------------------------------
!> @details 
!> Compute cutt-off filter in spectral space
!! @author  Antoine Vollant, LEGI
!! @param[in]    kx wave number
!! @param[in]    ky wave number
!! @param[in]    kz wave number
!! @param[in]    delta is  
!! @param[out]   valOut value of filter 
!---------------------------------------------------------------------------
 function computeCutOffFilterK(kx,ky,kz,delta) result(valOut)

  !IO
  real(WP), intent(in)            :: kx,ky,kz
  real(WP),intent(in)             :: delta
  real(WP)                        :: valOut

  !Local data
  real(WP)                        :: pi
  real(WP)                        :: kc 

  pi = acos(-1.0_WP)
  kc = pi/delta
  if ( sqrt(kx**2.0+ky**2.0+kz**2.0) .le. kc ) then
    valOut = 1.0_WP
  else
    valOut = 0.0_WP
  endif

 end function computeCutOffFilterK


!---------------------------------------------------------------------------
!> @details 
!> Compute box filter in spectral space
!! @author  Antoine Vollant, LEGI
!! @param[in]    kx wave number
!! @param[in]    ky wave number
!! @param[in]    kz wave number
!! @param[in]    delta is  
!! @param[out]   valOut value of filter 
!---------------------------------------------------------------------------
 function computeBoxFilterK(kx,ky,kz,delta) result(valOut)

   
  !IO
  real(WP), intent(in)            :: kx,ky,kz
  real(WP),intent(in)             :: delta
  real(WP)                        :: valOut

!   if(kx.eq.0.and.ky.eq.0.and.kz.eq.0) then 
   if( equalToVal(kx,0.0_WP) .and. &
     & equalToVal(ky,0.0_WP) .and. &
     & equalToVal(kz,0.0_WP) ) then 
     valOut = 1.0_WP
!   elseif( kx.eq.0 .and. ky.eq.0 ) then 
   elseif( equalToVal(kx,0.0_WP) .and. &
         & equalToVal(ky,0.0_WP) ) then 
     valOut = sin( 0.5_WP*delta*kz ) / ( 0.5_WP*delta*kz )
!   elseif( kx.eq.0 .and. kz.eq.0 ) then 
   elseif( equalToVal(kx,0.0_WP) .and. &
         & equalToVal(kz,0.0_WP) ) then 
     valOut = sin( 0.5_WP*delta*ky ) / ( 0.5_WP*delta*ky )
!   elseif( ky.eq.0 .and. kz.eq.0) then 
   elseif( equalToVal(ky,0.0_WP) .and. &
         & equalToVal(kz,0.0_WP) ) then 
     valOut = sin( 0.5_WP*delta*kx ) / ( 0.5_WP*delta*kx )
!   elseif( kx.eq.0 ) then 
   elseif( equalToVal(kx,0.0_WP) ) then 
     valOut = sin( 0.5_WP*delta*ky ) / ( 0.5_WP*delta*ky )* &
            & sin( 0.5_WP*delta*kz ) / ( 0.5_WP*delta*kz )
!   elseif( ky.eq.0 ) then 
   elseif( equalToVal(ky,0.0_WP) ) then 
     valOut = sin( 0.5_WP*delta*kx ) / ( 0.5_WP*delta*kx )* &
            & sin( 0.5_WP*delta*kz ) / ( 0.5_WP*delta*kz )
!   elseif( kz.eq.0) then 
   elseif( equalToVal(kz,0.0_WP) ) then 
     valOut = sin( 0.5_WP*delta*kx ) / ( 0.5_WP*delta*kx )* &
            & sin( 0.5_WP*delta*ky ) / ( 0.5_WP*delta*ky )
   else 
     valOut = sin( 0.5_WP*delta*kx ) / ( 0.5_WP*delta*kx )* &
            & sin( 0.5_WP*delta*ky ) / ( 0.5_WP*delta*ky )* &
            & sin( 0.5_WP*delta*kz ) / ( 0.5_WP*delta*kz )
   endif

 end function computeBoxFilterK


!---------------------------------------------------------------------------
!> @details 
!> Compute Gaussian filter in spectral space 
!! @author  Antoine Vollant, LEGI
!! @param[in]    kx wave number
!! @param[in]    ky wave number
!! @param[in]    kz wave number
!! @param[in]    delta is  
!! @return       valOut = value of filter 
!---------------------------------------------------------------------------
 function computeGaussianFilterK (kx,ky,kz,delta) result(valOut)
   
  !IO
  real(WP),intent(in)             :: kx,ky,kz
  real(WP),intent(in)             :: delta
  real(WP)                        :: valOut

  valOut = exp( - delta**2 * ( kx**2.0 + ky**2.0 + kz**2.0 ) / 24.0_WP )

 end function computeGaussianFilterK

!!------------------------------------------------------------------------------
!!> @author 
!!> Antoine Vollant, LEGI
!!>
!!> @details
!!> Compute the filter size Delta 
!!> @param [in] exField is on of the field which contains informations on the grid 
!!> @return  delta
!!------------------------------------------------------------------------------
! function computeDelta_anisotrope(exField) result (delta)
!
!   !I/O
!   type(real_data_layout),intent(in)   :: exField
!   real(WP),dimension(3)               :: delta  
!   
!   delta(1) = exField%Lx / exField%nx 
!   delta(2) = exField%Ly / exField%ny 
!   delta(3) = exField%Lz / exField%nz 
!
! end function computeDelta_anisotrope
 
!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> Compute the filter size Delta 
!> @param [in] exField is on of the field which contains informations on the grid 
!> @return  delta
!------------------------------------------------------------------------------
 function computeDelta_isotrope(exField) result (delta)

   !I/O
   type(real_data_layout),intent(in)   :: exField
   real(WP)                            :: delta  
   
   delta = ( exField%Lx * exField%Ly * exField%Lz) / &
         & (exField%nx * exField%ny * exField%nz )
   delta = delta**(1.0_WP/3.0_WP)


 end function computeDelta_isotrope

end module filtering_tools
!> @}
