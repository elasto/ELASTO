!USEFORTEST avgcond
!> @addtogroup fann 
!! @{
!------------------------------------------------------------------------------
!
! MODULE: fann 
!
!> @author
!> Antoine Vollant
!
! DESCRIPTION: 
!> The aim of this module is to provide fann tools
!> 
!------------------------------------------------------------------------------


module fann


  use datalayout

  implicit none


  contains

!------------------------------------------------------------------------------
!
!
!> @author
!> Antoine Vollant
!> param[in]    fieldIn vector which contains input entries for ann
!> param[in]    nameOfAnnFile is the name of ann file 
!> param[inout] fieldOut vector which contains the out of ann 
!> 
!------------------------------------------------------------------------------
  function computeEOFromANN(fieldIn,fieldOut,nameOfAnnFile) result(success)

    USE, intrinsic::iso_c_binding,only:c_char,c_null_char

    type(real_data_layout),dimension(:),intent(in)    :: fieldIn
    !type(real_data_layout),dimension(:),intent(inout) :: fieldOut
    type(real_data_layout),intent(inout)              :: fieldOut
    character(len=*),intent(in)                       :: nameOfAnnFile
    logical                                           :: success

    real,dimension(:,:),allocatable     :: fieldInBuf
    real,dimension(:,:),allocatable     :: fieldOutBuf
    integer                             :: i,j,k,l,m,ires
    integer                             :: nbInput,numData
    integer                             :: nbOutput
    CHARACTER(kind=c_char,len=128)      :: nomC

    success = .false.
    !if (.not. samelayout(fieldIn(1),fieldOut(1))) return
    if (.not. samelayout(fieldIn(1),fieldOut)) return
    nbInput=size(fieldIn,1)
    !nbOutput=size(fieldOut,1)
    nbOutput=1
    numData=(fieldIn(1)%zmax-fieldIn(1)%zmin+1) &
           &*(fieldIn(1)%ymax-fieldIn(1)%ymin+1) &
           &*(fieldIn(1)%xmax-fieldIn(1)%xmin+1)
    allocate(fieldInBuf(numData,nbInput),stat=ires) 
    if (ires .ne. 0) return
    allocate(fieldOutBuf(numData,nbOutput),stat=ires) 
    if (ires .ne. 0) return
    m=1
    do k =fieldIn(1)%zmin,fieldIn(1)%zmax
      do j =fieldIn(1)%ymin,fieldIn(1)%ymax
        do i =fieldIn(1)%xmin,fieldIn(1)%xmax
          do l=1,nbInput
            fieldInBuf(m,l)=fieldIn(l)%values(i,j,k)
          enddo
          m = m + 1
        enddo
      enddo
    enddo
    !call C routine
    nomC=trim(adjustl(nameOfAnnFile))//c_null_char
#ifdef FANN
    call run_fann_with_scales(fieldInBuf,fieldOutBuf,nomC,numData,nbInput,nbOutput)
#else
    print *,'[WARNING] FANN not compiled !!'
#endif FANN
    m=1
    do k =fieldIn(1)%zmin,fieldIn(1)%zmax
      do j =fieldIn(1)%ymin,fieldIn(1)%ymax
        do i =fieldIn(1)%xmin,fieldIn(1)%xmax
          do l=1,nbOutput
            !fieldOut(l)%values(i,j,k)=fieldOutBuf(m,1)
            fieldOut%values(i,j,k)=fieldOutBuf(m,l)
          enddo
          m = m + 1
        enddo
      enddo
    enddo
    deallocate(fieldInBuf)
    deallocate(fieldOutBuf)
    success = .true.

  end function computeEOFromANN

end module fann
