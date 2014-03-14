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


module toolbox

 use precision_tools 
 use parallel_tools
 use datalayout
 use mpilayout_tools
 use wavenumber_tools
 
 implicit none

 public


 interface initWorkArray 
   module procedure initWorkArrayK
   module procedure initWorkArrayP 
 end interface initWorkArray

 interface deleteWorkArray 
   module procedure deleteWorkArrayK
   module procedure deleteWorkArrayP 
 end interface deleteWorkArray

 contains





!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> Init a work array 
!> @param [in] modelK is on of the field which contains informations on the grid 
!> @param [in] nb is the size of the array you wish to initialize 
!> @param [inout] arrayinitK is the array of complex data layout of dim = nb
!------------------------------------------------------------------------------
 function initWorkArrayK(modelK,nb,arrayinitK) result(success)

  !I/O
  type(complex_data_layout), intent(in)                              :: modelK
  integer,intent(in)                                                 :: nb
  type(complex_data_layout), intent(inout),allocatable, dimension(:) :: arrayinitK
  logical ::success 
  
  !Local
  integer :: i,ierr

  success = .false.

  allocate( arrayinitK(nb) ,stat=ierr )
  if (ierr .ne.0) then 
     write(6,'(a)')'[ERROR] initWorkArrayK : not enought memory for allocation'
     return
  endif

  do i=1,nb
    if (.not. copyStructOnly(modelK,arrayinitK(i) )) then
      write (6,'(a)')'[ERROR] initWorkArrayK : not enought memory for allocation'
      return
    endif
  enddo

  success = .true.

 end function initWorkArrayK

!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> Init a work array 
!> @param [in] modelP is on of the field which contains informations on the grid 
!> @param [in] nb is the size of the array you wish to initialize 
!> @param [inout] arrayinitP is the array of real data layout of dim = nb
!------------------------------------------------------------------------------
 function initWorkArrayP(modelP,nb,arrayinitP) result(success)

  !I/O
  type(real_data_layout), intent(in)                              :: modelP
  integer,intent(in)                                              :: nb
  type(real_data_layout), intent(inout),allocatable, dimension(:) :: arrayinitP
  logical ::success 
  
  !Local
  integer :: i,ierr

  success = .false.

  allocate( arrayinitP(nb) ,stat=ierr )
  if (ierr .ne.0) then 
     write(6,'(a)')'[ERROR] initWorkArrayP : not enought memory for allocation'
     return
  endif
  do i=1,nb
    if (.not. copyStructOnly(modelP,arrayinitP(i) )) then
      write (6,'(a)')'[ERROR] initWorkArrayK : not enought memory for allocation'
      return
    endif
  enddo

  success = .true.

 end function initWorkArrayP

!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> delete a work array 
!> @param [inout] arrayinitP is the array of real data layout to delete
!------------------------------------------------------------------------------
 function deleteWorkArrayP(arrayinitP,iMin,iMax) result(success)

  !I/O
  type(real_data_layout), intent(inout), dimension(:) :: arrayinitP
  logical                                             :: success 
  integer,intent(in),optional                         :: iMin,iMax

  !Local
  integer :: i,nb,mini,maxi

  success = .false.

  if (present(iMin) .and. present(iMax) ) then
    mini = iMin 
    maxi = iMax
  else
    nb = size( arrayinitP )
    mini = 1
    maxi = nb
  endif

  do i=mini,maxi
    call deleteDataLayout( arrayinitP(i) )
  enddo

  success = .true.

 end function deleteWorkArrayP

!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> delete a work array 
!> @param [inout] arrayinitK is the array of complex data layout to delete
!------------------------------------------------------------------------------
 function deleteWorkArrayK(arrayinitK,iMin,iMax) result(success)

  !I/O
  type(complex_data_layout), intent(inout), dimension(:) :: arrayinitK
  logical                                                :: success 
  integer,intent(in),optional                            :: iMin,iMax

  !Local
  integer :: i,nb,mini,maxi

  success = .false.

  if (present(iMin) .and. present(iMax) ) then
    mini = iMin 
    maxi = iMax
  else
    nb = size( arrayinitK )
    mini = 1
    maxi = nb
  endif

  do i=mini,maxi
    call deleteDataLayout( arrayinitK(i) )
  enddo

  success = .true.

 end function deleteWorkArrayK


!---------------------------------------------------------------------------
!> @details
!> Compute the number of i^th sample where the targetValue is located in the range
!> maxValue minus minValue uniformly sampled
!! @author  Antoine Vollant, LEGI
!! @param[in]    minValue is the min of range
!! @param[in]    maxValue is the max of range
!! @param[in]    targetValue is the value to locate
!! @param[in]    sampling is the number of sample in total range
!! @return       res is the i^th sample which contains the targetValue [1:sampling]
!---------------------------------------------------------------------------
function computeIndice(minValue,maxValue,targetValue,sampling) result (res)


  !I/O
  real(WP),intent(in)     :: minValue,maxValue,targetValue
  integer ,intent(in)     :: sampling
  integer                 :: res

  !Local data
  real(WP) :: tx,jjfloat,eps
  real(WP) :: newMin,newMax,delta

#ifdef DEBUGTOOLBOX
  if ( targetValue .lt. minValue .or. targetValue .gt. maxValue ) then
    write(6,'(a)')'[ERROR] in computeIndice targetValue is not in range'
  endif
  if (  minValue .gt. maxValue ) then
    write(6,'(a)')'[ERROR] in computeIndice Min and Max are not correct'
  endif
#endif

  !In order to never have targetValue equal to minValue or maxValue
  eps = 1d-5
  delta = maxValue - minValue
  newMin = minValue - eps * delta
  newMax = maxValue + eps * delta
  delta = newMax - newMin

  tx=targetValue - newMin
  tx=tx/delta
  jjfloat=tx*real(sampling,WP) !convert in double precision
  res = ceiling( jjfloat )  !round to the next integer

#ifdef DEBUGTOOLBOX
  if ( res .lt. 0 .or. res .gt. sampling ) then
    write(6,'(a)')'[ERROR] in computeIndice'
  endif
#endif

end function


end module toolbox
!> @}
