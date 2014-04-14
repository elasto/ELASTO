!> @addtogroup sub_grid_models
!! @{

module scamodels


 use datalayout
 use toolbox
 use stat_tools 
 use differential_tools 
 use transforms_tools
 use wavenumber_tools 
 use subgrid_tools
 use conditional_mean_tools
 use filtering_tools 
 use physical_values_tools

 implicit none
 private

 public :: DynClarkFabreSca
 public :: DynClarkSca
 public :: DynClarkFabreSca_3d
 public :: DynClarkSca_3d
 public :: DynClarkFabreSca_2d_z
 public :: DynClarkSca_2d_z
 public :: DynWangSca
 public :: gradientSca
 public :: smagorinskySca
 public :: smagorinskySca2
 public :: coefDynSca2
 public :: DynRGMSca
 public :: RGMOSca
 public :: RGMSca_lapack
 public :: RGMScaForTest
 public :: RGMSca_analyt
 public :: DynSmagorinskySca_2d_z
 public :: DynSmagorinskySca_3d
 public :: DynSmagorinskySca
 public :: scaleSimilaritySca
 !public :: sism 

 contains



!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> Sub-grid model based on NDCM one variable
!> @param [in]    U  longitudinal velocity in physical space
!> @param [in]    V  velocity in physical space
!> @param [in]    W  velocity in physical space
!> @param [in]    Uk longitudinal velocity in spectral space
!> @param [in]    Vk velocity in spectral space
!> @param [in]    Wk velocity in spectral space
!> @param [in]    ScalIn scalar in physical space 
!> @param [in]    ScalkIn scalar in spectral space 
!> @param [in]    WaveNum wave numbers for scalar 
!> @param [in]    Filter is the type of test filter 
!> @param [in]    spec_rank is the rank of processus
!> @param [inout] res logical
!> @param [inout] div is the divergence of SGS model
!> @return
!------------------------------------------------------------------------------
 subroutine DynClarkFabreSca(divK,U,V,W,Uk,Vk,Wk,Scal,ScalK,ScalWN,filter,spec_rank,res,type_avg)
  !I/O
  type(complex_data_layout),intent(inout) :: divK
  type(complex_data_layout),intent(in)    :: Uk,Vk,Wk,ScalK
  type(real_data_layout),intent(in)       :: U, V, W, Scal
  type(WaveNumbers), intent(in)           :: ScalWN
  logical, intent(inout)                  :: res
  integer, intent(in)                     :: spec_rank
  integer, intent(in)                     :: filter
  integer, intent(in)                     :: type_avg

  !Local data
  type(real_data_layout),dimension(:),allocatable    :: divSmag
  type(real_data_layout),dimension(:),allocatable    :: divGrad
  type(real_data_layout),dimension(:),allocatable    :: div
  type(real_data_layout),dimension(:),allocatable    :: Ti
  type(complex_data_layout),dimension(:),allocatable :: divGradK
  type(complex_data_layout),dimension(:),allocatable :: divSmagK
  type(complex_data_layout),dimension(:),allocatable :: TiK
  real(WP)                                           :: deltx
  real(WP),dimension(:),allocatable                  :: coefArray
  real(WP)                                           :: coeff
  integer                                            :: ires,k
 
  res = .true.

  !Init des tableaux de travail
  if (.not. initWorkArray(U,1,divSmag) .or. &
     &.not. initWorkArray(U,1,divGrad) .or. &
     &.not. initWorkArray(U,1,div) .or. &
     &.not. initWorkArray(U,3,Ti) .or. &
     &.not. initWorkArray(Uk,1,divSmagK) .or. &
     &.not. initWorkArray(Uk,1,divGradK) .or. &
     &.not. initWorkArray(Uk,3,TiK)) then
    write(6,'(a)')'[ERROR] initWorkArray : not enought memory!'
    return
  endif

  deltx = computeDelta(Scal)
  !Gradient
  call gradientSca(divGradK(1),U,Uk,Vk,Wk,ScalK,Scal,ScalWN,res,1)
  call btran(divGradK(1),divGrad(1),res)
  !Smagorinsky
  call computeSdZ(U,Uk,Vk,Wk,ScalK,ScalWN,Ti(1),Ti(2),Ti(3),res)

  select case (type_avg)
    case(0) !Cas isotrope
      if (.not. coefDynClarkFabreSca(U,V,W,Uk,Vk,Wk,Scal,ScalK,ScalWN,filter,1,spec_rank,&
         & CoeffDynScalar=coeff) ) then
        write (6,'(a)') '[ERROR] in DynClarkFabreSca'
        res = .false.
      endif
      Ti(1)%values = coeff * deltx**2 * Ti(1)%values
      Ti(2)%values = coeff * deltx**2 * Ti(2)%values
      Ti(3)%values = coeff * deltx**2 * Ti(3)%values
    case(1) !Cas anisotrope moyenne dans les plans (X,Y)
      allocate(coefArray(Scal%nz),stat=ires)
      if ( ires .ne. 0 ) then
        write (6,'(a)') '[ERROR] not enought memory to allocate coefArray'
        res = .false.
      endif 
      if (.not. coefDynClarkFabreSca(U,V,W,Uk,Vk,Wk,Scal,ScalK,ScalWN,filter,1,&
          &  spec_rank,CoeffDynArray=coefArray,dir2DAvg=3) ) then
          write (6,'(a)') '[ERROR] Coef dynamic failed !'
        res = .false.
      endif
      do k = U%zmin,U%zmax
         Ti(1)%values(:,:,k) = coefArray(k) * deltx**2 * Ti(1)%values(:,:,k)
         Ti(2)%values(:,:,k) = coefArray(k) * deltx**2 * Ti(2)%values(:,:,k)
         Ti(3)%values(:,:,k) = coefArray(k) * deltx**2 * Ti(3)%values(:,:,k)
      end do
      deallocate(coefArray)
    case(2) !Cas anisotrope moyenne dans les plans (X,Z)
      allocate(coefArray(Scal%ny),stat=ires)
      if ( ires .ne. 0 ) then
        write (6,'(a)') '[ERROR] not enought memory to allocate coeff'
        res = .false.
      endif 
      if (.not. coefDynClarkFabreSca(U,V,W,Uk,Vk,Wk,Scal,ScalK,ScalWN,filter,1,&
          &  spec_rank,CoeffDynArray=coefArray,dir2DAvg=2) ) then
        write (6,'(a)') '[ERROR] in DynSmagorinskySca_3d'
        res = .false.
      endif
      do k = U%ymin,U%ymax
         Ti(1)%values(:,k,:) = coefArray(k) * deltx**2 * Ti(1)%values(:,k,:)
         Ti(2)%values(:,k,:) = coefArray(k) * deltx**2 * Ti(2)%values(:,k,:)
         Ti(3)%values(:,k,:) = coefArray(k) * deltx**2 * Ti(3)%values(:,k,:)
      end do
      deallocate(coefArray)
    case(3) !Cas anisotrope moyenne dans les plans (Y,Z)
      allocate(coefArray(Scal%nx),stat=ires)
      if ( ires .ne. 0 ) then
        write (6,'(a)') '[ERROR] not enought memory to allocate coeff'
        res = .false.
      endif 
      if (.not. coefDynClarkFabreSca(U,V,W,Uk,Vk,Wk,Scal,ScalK,ScalWN,filter,1,&
          &  spec_rank,CoeffDynArray=coefArray,dir2DAvg=1) ) then
        write (6,'(a)') '[ERROR] in DynClarkFabreSca'
        res = .false.
      endif
      do k = U%xmin,U%xmax
         Ti(1)%values(k,:,:) = coefArray(k) * deltx**2 * Ti(1)%values(k,:,:)
         Ti(2)%values(k,:,:) = coefArray(k) * deltx**2 * Ti(2)%values(k,:,:)
         Ti(3)%values(k,:,:) = coefArray(k) * deltx**2 * Ti(3)%values(k,:,:)
      end do
      deallocate(coefArray)
    case default
      write (6,'(a)') '[ERROR] type of averaging unknow in Dynamic Clark Fabre procedure'
      res = .false.
  end select

  call ftran(Ti(1),TiK(1),res)
  if (.not. res ) write(6,'(a)')'[ERROR] in DynClarkFabreSca: compute ftran 1 failed'
  call ftran(Ti(2),TiK(2),res)
  if (.not. res ) write(6,'(a)')'[ERROR] in DynClarkFabreSca: compute ftran 2 failed'
  call ftran(Ti(3),TiK(3),res)
  if (.not. res ) write(6,'(a)')'[ERROR] in DynClarkFabreSca: compute ftran 3 failed'
  call computeDivergenceK(ScalWN,TiK(1),TiK(2),TiK(3),divSmagK(1))
  call btran(divSmagK(1),divSmag(1),res)
  !Mixage
  div(1)%values = divSmag(1)%values + divGrad(1)%values
  call ftran(div(1),divK,res)

  !Désallocation des tableaux de travail
  if ( .not. deleteWorkArray(divGrad) .or. &
  & .not. deleteWorkArray(divGradK) .or. &
  & .not. deleteWorkArray(divSmag) .or. &
  & .not. deleteWorkArray(divSmagK) .or. &
  & .not. deleteWorkArray(div) .or. &
  & .not. deleteWorkArray(Ti) .or. &
  & .not. deleteWorkArray(TiK)) then
    res = .false.
    write(6,'(a)')'[ERROR] in dynClarkFabresca: memory leak'
  endif

 end subroutine DynClarkFabreSca


!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> Sub-grid model based on DCM one variable
!> @param [in]    U  longitudinal velocity in physical space
!> @param [in]    V  velocity in physical space
!> @param [in]    W  velocity in physical space
!> @param [in]    Uk longitudinal velocity in spectral space
!> @param [in]    Vk velocity in spectral space
!> @param [in]    Wk velocity in spectral space
!> @param [in]    tabEOFile contains the EO model
!> @param [in]    Scal scalar in physical space 
!> @param [in]    Scalk scalar in spectral space 
!> @param [in]    ScalWN wave numbers for scalar 
!> @param [in]    spec_rank is the rank of processus
!> @param [inout] res logical
!> @param [inout] div is the divergence of SGS model
!> @return
!------------------------------------------------------------------------------
 subroutine DynClarkSca(divK,U,V,W,Uk,Vk,Wk,Scal,ScalK,ScalWN,filter,spec_rank,res,type_avg)
  !I/O
  type(complex_data_layout),intent(inout) :: divK
  type(complex_data_layout),intent(in)    :: Uk,Vk,Wk,ScalK
  type(real_data_layout),intent(in)       :: U, V, W, Scal
  type(WaveNumbers), intent(in)           :: ScalWN 
  logical, intent(inout)                  :: res
  integer, intent(in)                     :: spec_rank
  integer, intent(in)                     :: type_avg
  integer, intent(in)                     :: filter

  !Local data
  type(real_data_layout),dimension(:),allocatable    :: divSmag
  type(real_data_layout),dimension(:),allocatable    :: divGrad
  type(real_data_layout),dimension(:),allocatable    :: div
  type(real_data_layout),dimension(:),allocatable    :: Ti
  type(complex_data_layout),dimension(:),allocatable :: divGradK
  type(complex_data_layout),dimension(:),allocatable :: divSmagK
  type(complex_data_layout),dimension(:),allocatable :: TiK
  real(WP)                                           :: deltx
  real(WP)                                           :: coeff
  real(WP),dimension(:),allocatable                  :: coefArray
  integer                                            :: ires,k
 
  res = .true.

  !Init des tableaux de travail
  if (.not. initWorkArray(U,1,divSmag) .or. &
     &.not. initWorkArray(U,1,divGrad) .or. &
     &.not. initWorkArray(U,1,div) .or. &
     &.not. initWorkArray(U,3,Ti) .or. &
     &.not. initWorkArray(Uk,1,divSmagK) .or. &
     &.not. initWorkArray(Uk,1,divGradK) .or. &
     &.not. initWorkArray(Uk,3,TiK)) then
    write(6,'(a)')'[ERROR] initWorkArray : not enought memory!'
    return
  endif


  deltx = computeDelta(Scal)
  !Gradient
  call gradientSca(divGradK(1),U,Uk,Vk,Wk,ScalK,Scal,ScalWN,res,1)
  call btran(divGradK(1),divGrad(1),res)
  !Smagorinsky
  call computeSdZ(U,Uk,Vk,Wk,ScalK,ScalWN,Ti(1),Ti(2),Ti(3),res)


  select case (type_avg)
    case(0) !Cas isotrope
      if (.not. coefDynClarkSca(U,V,W,Uk,Vk,Wk,Scal,ScalK,ScalWN,filter,1,spec_rank,&
         & CoeffDynScalar=coeff) ) then
        write (6,'(a)') '[ERROR] in DynClarkFabreSca'
        res = .false.
      endif
      Ti(1)%values = coeff * deltx**2 * Ti(1)%values
      Ti(2)%values = coeff * deltx**2 * Ti(2)%values
      Ti(3)%values = coeff * deltx**2 * Ti(3)%values

    case(1) !Cas anisotrope moyenne dans les plans (X,Y)
      allocate(coefArray(Scal%nz),stat=ires)
      if ( ires .ne. 0 ) then
        write (6,'(a)') '[ERROR] not enought memory to allocate coefArray'
        res = .false.
      endif 
      if (.not. coefDynClarkSca(U,V,W,Uk,Vk,Wk,Scal,ScalK,ScalWN,filter,1,&
          &  spec_rank,CoeffDynArray=coefArray,dir2DAvg=3) ) then
          write (6,'(a)') '[ERROR] in DynClarkSca' 
        res = .false.
      endif
      do k = U%zmin,U%zmax
         Ti(1)%values(:,:,k) = coefArray(k) * deltx**2 * Ti(1)%values(:,:,k)
         Ti(2)%values(:,:,k) = coefArray(k) * deltx**2 * Ti(2)%values(:,:,k)
         Ti(3)%values(:,:,k) = coefArray(k) * deltx**2 * Ti(3)%values(:,:,k)
      end do
      deallocate(coefArray)
    case(2) !Cas anisotrope moyenne dans les plans (X,Z)
      allocate(coefArray(Scal%ny),stat=ires)
      if ( ires .ne. 0 ) then
        write (6,'(a)') '[ERROR] not enought memory to allocate coeff'
        res = .false.
      endif 
      if (.not. coefDynClarkSca(U,V,W,Uk,Vk,Wk,Scal,ScalK,ScalWN,filter,1,&
          &  spec_rank,CoeffDynArray=coefArray,dir2DAvg=2) ) then
        write (6,'(a)') '[ERROR] in DynClarkSca'
        res = .false.
      endif
      do k = U%ymin,U%ymax
         Ti(1)%values(:,k,:) = coefArray(k) * deltx**2 * Ti(1)%values(:,k,:)
         Ti(2)%values(:,k,:) = coefArray(k) * deltx**2 * Ti(2)%values(:,k,:)
         Ti(3)%values(:,k,:) = coefArray(k) * deltx**2 * Ti(3)%values(:,k,:)
      end do
      deallocate(coefArray)
    case(3) !Cas anisotrope moyenne dans les plans (Y,Z)
      allocate(coefArray(Scal%nx),stat=ires)
      if ( ires .ne. 0 ) then
        write (6,'(a)') '[ERROR] not enought memory to allocate coeff'
        res = .false.
      endif 
      if (.not. coefDynClarkSca(U,V,W,Uk,Vk,Wk,Scal,ScalK,ScalWN,filter,1,&
          &  spec_rank,CoeffDynArray=coefArray,dir2DAvg=1) ) then
        write (6,'(a)') '[ERROR] in DynClarkSca'
        res = .false.
      endif
      do k = U%xmin,U%xmax
         Ti(1)%values(k,:,:) = coefArray(k) * deltx**2 * Ti(1)%values(k,:,:)
         Ti(2)%values(k,:,:) = coefArray(k) * deltx**2 * Ti(2)%values(k,:,:)
         Ti(3)%values(k,:,:) = coefArray(k) * deltx**2 * Ti(3)%values(k,:,:)
      end do
      deallocate(coefArray)
    case default
      write (6,'(a)') '[ERROR] type of averaging unknow in Dynamic Clark procedure'
      res = .false.
  end select

  call ftran(Ti(1),TiK(1),res)
  if (.not. res ) write(6,'(a)')'[ERROR] in DynClarkSca : compute ftran 1 failed'
  call ftran(Ti(2),TiK(2),res)
  if (.not. res ) write(6,'(a)')'[ERROR] in DynClarkSca : compute ftran 2 failed'
  call ftran(Ti(3),TiK(3),res)
  if (.not. res ) write(6,'(a)')'[ERROR] in DynClarkSca : compute ftran 3 failed'
  call computeDivergenceK(ScalWN,TiK(1),TiK(2),TiK(3),divSmagK(1))
  call btran(divSmagK(1),divSmag(1),res)
  !Mixage
  div(1)%values = divSmag(1)%values + divGrad(1)%values
  call ftran(div(1),divK,res)
  !Désallocation des tableaux de travail
  if ( .not. deleteWorkArray(divGrad) .or. &
  & .not. deleteWorkArray(divGradK) .or. &
  & .not. deleteWorkArray(divSmag) .or. &
  & .not. deleteWorkArray(divSmagK) .or. &
  & .not. deleteWorkArray(div) .or. &
  & .not. deleteWorkArray(Ti) .or. &
  & .not. deleteWorkArray(TiK)) then
    res = .false.
    write(6,'(a)')'[ERROR] in dynClarksca_3d: memory leak'
  endif

 end subroutine DynClarkSca


!------------------------------------------------------------------------------
!> @author 
!> Guillaume Balarac, LEGI
!>
!> @details
!> Dynamic Wang model
!> @param [in]    U  longitudinal velocity in physical space
!> @param [in]    Uk longitudinal velocity in spectral space
!> @param [in]    Vk velocity in spectral space
!> @param [in]    Wk velocity in spectral space
!> @param [in]    Scalk scalar in spectral space 
!> @param [in]    ScalWN wave numbers for scalar 
!> @param [in]    spec_rank is the rank of processus
!> @param [inout] res logical
!> @param [out]   div is the divergence of SGS model
!------------------------------------------------------------------------------
 subroutine DynWangSca(div,  U,V,W,Uk,VK,Wk,Scal,ScalK,ScalWN,filter,spec_rank,res)
 
   !I/O
  type(complex_data_layout),intent(inout) :: div
  type(complex_data_layout),intent(in)    :: Uk,Vk,Wk,ScalK
  type(real_data_layout),intent(in)       :: U, V, W, Scal
  type(WaveNumbers), intent(in)           :: ScalWN
  logical, intent(inout)                  :: res
  integer, intent(in)                     :: filter, spec_rank

  !Local data
  type(real_data_layout),dimension(:),allocatable    :: temp,temp1, temp2, temp3, temp4, temp5, temp6,Li, Mi, Qi, Pi
  type(real_data_layout),dimension(:),allocatable    :: VelScal_ft, Vel_ft, dZdxi_ft, SmZ, SdZdxi_ft
  type(real_data_layout)                             :: S11, S12, S13, S22, S23, S33
  type(real_data_layout)                             :: S11_ft, S12_ft, S13_ft, S22_ft, S23_ft, S33_ft, Scal_ft
  type(complex_data_layout),dimension(:),allocatable :: tempK, tempKf
  real(WP)                                           :: filtersize, DynCoef, qe_gm, ie_gm
  integer                                            :: i,j,k
 
  res = .true.
  filtersize = computeDelta(Scal)
  
  if ( .not. initWorkArray(U,3,temp) .or. &
  & .not. initWorkArray(U,3,temp1) .or. &
  & .not. initWorkArray(U,3,temp2) .or. &
  & .not. initWorkArray(U,3,temp3) .or. &
  & .not. initWorkArray(U,3,temp4) .or. &
  & .not. initWorkArray(U,3,temp5) .or. &
  & .not. initWorkArray(U,3,temp6) .or. &
  & .not. initWorkArray(U,3,SmZ) .or. &
  & .not. initWorkArray(U,3,Li) .or. &
  & .not. initWorkArray(U,3,Mi) .or. &
  & .not. initWorkArray(U,3,Qi) .or. &
  & .not. initWorkArray(U,3,Pi) .or. &
  & .not. initWorkArray(U,3,VelScal_ft) .or. &
  & .not. initWorkArray(U,3,Vel_ft) .or. &
  & .not. initWorkArray(U,3,dZdxi_ft) .or. &
  & .not. initWorkArray(U,3,SdZdxi_ft) .or. &
  & .not. initWorkArray(Uk,4,tempK) .or. &
  & .not. initWorkArray(Uk,4,tempKf)) then
    write(6,'(a)')'[ERROR] initWorkArray : not enought memory!'
    return
  endif
  
  call computeGradientK(ScalWN,ScalK,tempK(1),tempK(2),tempK(3))
  call btran(tempK(1),temp4(1),res)
  call btran(tempK(2),temp4(2),res)
  call btran(tempK(3),temp4(3),res)
  call computeGradientK(ScalWN,UK,tempK(1),tempK(2),tempK(3))
  call btran(tempK(1),temp1(1),res)
  call btran(tempK(2),temp1(2),res)
  call btran(tempK(3),temp1(3),res)
  call computeGradientK(ScalWN,VK,tempK(1),tempK(2),tempK(3))
  call btran(tempK(1),temp2(1),res)
  call btran(tempK(2),temp2(2),res)
  call btran(tempK(3),temp2(3),res)
  call computeGradientK(ScalWN,WK,tempK(1),tempK(2),tempK(3))
  call btran(tempK(1),temp3(1),res)
  call btran(tempK(2),temp3(2),res)
  call btran(tempK(3),temp3(3),res)
  
  if (.not. (copyStructOnly(U,S11) .and. &
  & copyStructOnly(U,S12) .and. &
  & copyStructOnly(U,S13) .and. &
  & copyStructOnly(U,S22) .and. &
  & copyStructOnly(U,S23) .and. &
  & copyStructOnly(U,S33) .and. &
  & copyStructOnly(U,S11_ft) .and. &
  & copyStructOnly(U,S12_ft) .and. &
  & copyStructOnly(U,S13_ft) .and. &
  & copyStructOnly(U,S22_ft) .and. &
  & copyStructOnly(U,S23_ft) .and. &
  & copyStructOnly(U,S33_ft) .and. &
  & copyStructOnly(U,Scal_ft)) ) then
    write(6,'(a)')'[ERROR] initWorkArray : not enought memory!'
    return
  endif
  
  ! Compute gradients
  call computeGradientK(ScalWN,ScalK,tempK(1),tempK(2),tempK(3))
  call btran(tempK(1),temp4(1),res)
  call btran(tempK(2),temp4(2),res)
  call btran(tempK(3),temp4(3),res)
  call computeGradientK(ScalWN,UK,tempK(1),tempK(2),tempK(3))
  call btran(tempK(1),temp1(1),res)
  call btran(tempK(2),temp1(2),res)
  call btran(tempK(3),temp1(3),res)
  call computeGradientK(ScalWN,VK,tempK(1),tempK(2),tempK(3))
  call btran(tempK(1),temp2(1),res)
  call btran(tempK(2),temp2(2),res)
  call btran(tempK(3),temp2(3),res)
  call computeGradientK(ScalWN,WK,tempK(1),tempK(2),tempK(3))
  call btran(tempK(1),temp3(1),res)
  call btran(tempK(2),temp3(2),res)
  call btran(tempK(3),temp3(3),res)
  
  ! compute Sij = 1/2 * ( duidxj + dujdxi) and Oij = 1/2 * ( duidxj - dujdxi )
  S11%values = temp1(1)%values
  S12%values = 0.5 * ( temp1(2)%values + temp2(1)%values )
  S13%values = 0.5 * ( temp1(3)%values + temp3(1)%values )
  S22%values = temp2(2)%values
  S23%values = 0.5 * ( temp2(3)%values + temp3(2)%values )
  S33%values = temp3(3)%values
  ! compute |Sij| = sqrt(2.0*Sij*Sij)   <-- temp5(1)
  call computeTensorNorme1(S11,S12,S13,S22,S23,S33,temp5(1),res)
  
  ! For dynamic procedure
  ! Compute velocity and scalar test filter
  ! Compute TestFilter(ui) and TestFilter(z)
  call ftran(U,tempK(1),res)
  call ftran(V,tempK(2),res)
  call ftran(W,tempK(3),res)
  call ftran(Scal,tempK(4),res)
  call computeFilter(ScalWN,2.0*filtersize,tempK(1),tempKf(1),filter) !U bar
  call computeFilter(ScalWN,2.0*filtersize,tempK(2),tempKf(2),filter) !V bar
  call computeFilter(ScalWN,2.0*filtersize,tempK(3),tempKf(3),filter) !W bar
  call computeFilter(ScalWN,2.0*filtersize,tempK(4),tempKf(4),filter) !Z bar
  call btran(tempKf(1),Vel_ft(1),res) !U bar
  call btran(tempKf(2),Vel_ft(2),res) !V bar
  call btran(tempKf(3),Vel_ft(3),res) !W bar
  call btran(tempKf(4),Scal_ft,res) !Z bar
  ! Compute Li = TestFilter (uiz) - TestFilter(ui)*TestFilter(z)
  VelScal_ft(1)%values = U%values * Scal%values
  VelScal_ft(2)%values = V%values * Scal%values
  VelScal_ft(3)%values = W%values * Scal%values
  call ftran(VelScal_ft(1),tempK(1),res)
  call ftran(VelScal_ft(2),tempK(2),res)
  call ftran(VelScal_ft(3),tempK(3),res)
  call computeFilter(ScalWN,2.0*filtersize,tempK(1),tempKf(1),filter) !U bar
  call computeFilter(ScalWN,2.0*filtersize,tempK(2),tempKf(2),filter) !V bar
  call computeFilter(ScalWN,2.0*filtersize,tempK(3),tempKf(3),filter) !W bar
  call btran(tempKf(1),VelScal_ft(1),res) !U bar
  call btran(tempKf(2),VelScal_ft(2),res) !V bar
  call btran(tempKf(3),VelScal_ft(3),res) !W bar
  ! Compute Li
  Li(1)%values = VelScal_ft(1)%values - Vel_ft(1)%values * Scal_ft%values
  Li(2)%values = VelScal_ft(2)%values - Vel_ft(2)%values * Scal_ft%values
  Li(3)%values = VelScal_ft(3)%values - Vel_ft(3)%values * Scal_ft%values
  
  ! Smag Part - Mi
  ! Compute Sij_ft by filtering Sij
  call ftran(S11,tempK(1),res)
  call ftran(S12,tempK(2),res)
  call ftran(S13,tempK(3),res)
  call computeFilter(ScalWN,2.0*filtersize,tempK(1),tempKf(1),filter) !U bar
  call computeFilter(ScalWN,2.0*filtersize,tempK(2),tempKf(2),filter) !V bar
  call computeFilter(ScalWN,2.0*filtersize,tempK(3),tempKf(3),filter) !W bar
  call btran(tempKf(1),S11_ft,res) !U bar
  call btran(tempKf(2),S12_ft,res) !V bar
  call btran(tempKf(3),S13_ft,res) !W bar
  call ftran(S22,tempK(1),res)
  call ftran(S23,tempK(2),res)
  call ftran(S33,tempK(3),res)
  call computeFilter(ScalWN,2.0*filtersize,tempK(1),tempKf(1),filter) !U bar
  call computeFilter(ScalWN,2.0*filtersize,tempK(2),tempKf(2),filter) !V bar
  call computeFilter(ScalWN,2.0*filtersize,tempK(3),tempKf(3),filter) !W bar
  call btran(tempKf(1),S22_ft,res) !U bar
  call btran(tempKf(2),S23_ft,res) !V bar
  call btran(tempKf(3),S33_ft,res) !W bar
  ! Compute |Sij_ft|  <-- temp5(2)
  call computeTensorNorme1(S11_ft,S12_ft,S13_ft,S22_ft,S23_ft,S33_ft,temp5(2),res)
  ! Compute dZdxi_ft by filtering of dZdxi (which is in temp4)
  call ftran(temp4(1),tempK(1),res)
  call ftran(temp4(2),tempK(2),res)
  call ftran(temp4(3),tempK(3),res)
  call computeFilter(ScalWN,2.0*filtersize,tempK(1),tempKf(1),filter) !U bar
  call computeFilter(ScalWN,2.0*filtersize,tempK(2),tempKf(2),filter) !V bar
  call computeFilter(ScalWN,2.0*filtersize,tempK(3),tempKf(3),filter) !W bar
  call btran(tempKf(1),dZdxi_ft(1),res) !U bar
  call btran(tempKf(2),dZdxi_ft(2),res) !V bar
  call btran(tempKf(3),dZdxi_ft(3),res) !W bar
  ! Compute SdZdxi_ft
  SdZdxi_ft(1)%values = temp5(1)%values * temp4(1)%values
  SdZdxi_ft(2)%values = temp5(1)%values * temp4(2)%values
  SdZdxi_ft(3)%values = temp5(1)%values * temp4(3)%values
  call ftran(SdZdxi_ft(1),tempK(1),res)
  call ftran(SdZdxi_ft(2),tempK(2),res)
  call ftran(SdZdxi_ft(3),tempK(3),res)
  call computeFilter(ScalWN,2.0*filtersize,tempK(1),tempKf(1),filter) !U bar
  call computeFilter(ScalWN,2.0*filtersize,tempK(2),tempKf(2),filter) !V bar
  call computeFilter(ScalWN,2.0*filtersize,tempK(3),tempKf(3),filter) !W bar
  call btran(tempKf(1),SdZdxi_ft(1),res) !U bar
  call btran(tempKf(2),SdZdxi_ft(2),res) !V bar
  call btran(tempKf(3),SdZdxi_ft(3),res) !W bar
  ! Compute Mi = (2 * filtersize)^2* |Sij_ft|*dZdxi_ft - filtersize^2 * SdZdxi_ft
  Mi(1)%values = (2.0 * filtersize)**2.0 * temp5(2)%values*dZdxi_ft(1)%values - filtersize**2.0 * SdZdxi_ft(1)%values
  Mi(2)%values = (2.0 * filtersize)**2.0 * temp5(2)%values*dZdxi_ft(2)%values - filtersize**2.0 * SdZdxi_ft(2)%values
  Mi(3)%values = (2.0 * filtersize)**2.0 * temp5(2)%values*dZdxi_ft(3)%values - filtersize**2.0 * SdZdxi_ft(3)%values
  
  ! Pang Part
  Qi(1)%values = (2.0*filtersize)**2.0 * ( S11_ft%values*dZdxi_ft(1)%values + S12_ft%values*dZdxi_ft(2)%values + S13_ft%values*dZdxi_ft(3)%values )
  Qi(2)%values = (2.0*filtersize)**2.0 * ( S12_ft%values*dZdxi_ft(1)%values + S22_ft%values*dZdxi_ft(2)%values + S23_ft%values*dZdxi_ft(3)%values )
  Qi(3)%values = (2.0*filtersize)**2.0 * ( S13_ft%values*dZdxi_ft(1)%values + S23_ft%values*dZdxi_ft(2)%values + S33_ft%values*dZdxi_ft(3)%values )
  Pi(1)%values = S11%values*temp4(1)%values + S12%values*temp4(2)%values + S13%values*temp4(3)%values 
  Pi(2)%values = S12%values*temp4(1)%values + S22%values*temp4(2)%values + S23%values*temp4(3)%values 
  Pi(3)%values = S13%values*temp4(1)%values + S23%values*temp4(2)%values + S33%values*temp4(3)%values
  call ftran(Pi(1),tempK(1),res)
  call ftran(Pi(2),tempK(2),res)
  call ftran(Pi(3),tempK(3),res)
  call computeFilter(ScalWN,2.0*filtersize,tempK(1),tempKf(1),filter)
  call computeFilter(ScalWN,2.0*filtersize,tempK(2),tempKf(2),filter)
  call computeFilter(ScalWN,2.0*filtersize,tempK(3),tempKf(3),filter)
  call btran(tempKf(1),Pi(1),res) 
  call btran(tempKf(2),Pi(2),res) 
  call btran(tempKf(3),Pi(3),res)
  Qi(1)%values = Qi(1)%values - (filtersize)**2.0 * Pi(1)%values
  Qi(2)%values = Qi(2)%values - (filtersize)**2.0 * Pi(2)%values
  Qi(3)%values = Qi(3)%values - (filtersize)**2.0 * Pi(3)%values
  
  ! Compute Coef1 (for Smag part)
  temp5(3)%values =  (Li(1)%values*Mi(1)%values + Li(2)%values*Mi(2)%values + Li(3)%values*Mi(3)%values)*(Qi(1)%values*Qi(1)%values + Qi(2)%values*Qi(2)%values + Qi(3)%values*Qi(3)%values) &
               &  -  (Li(1)%values*Qi(1)%values + Li(2)%values*Qi(2)%values + Li(3)%values*Qi(3)%values)*(Mi(1)%values*Qi(1)%values + Mi(2)%values*Qi(2)%values + Mi(3)%values*Qi(3)%values) 
  temp5(2)%values =  (Qi(1)%values*Qi(1)%values + Qi(2)%values*Qi(2)%values + Qi(3)%values*Qi(3)%values)*(Mi(1)%values*Mi(1)%values + Mi(2)%values*Mi(2)%values + Mi(3)%values*Mi(3)%values) &
               &  -  (Mi(1)%values*Qi(1)%values + Mi(2)%values*Qi(2)%values + Mi(3)%values*Qi(3)%values)*(Mi(1)%values*Qi(1)%values + Mi(2)%values*Qi(2)%values + Mi(3)%values*Qi(3)%values)
  !ie_gm = computeFieldAvg( temp5(3) , spec_rank  )
  !qe_gm = computeFieldAvg( temp5(2) , spec_rank  )
  !Coef1 = ie_gm / qe_gm
  !if (spec_rank.eq.0) print*,' Wang model, coef1 =', Coef1
  temp6(1)%values = temp5(3)%values / temp5(2)%values
  ie_gm = computeFieldAvg( temp6(1) , spec_rank  )
  if (spec_rank.eq.0) print*,' Wang model, mean of coef1 - before clipping =', ie_gm
  do k = U%zmin, U%zmax
     do j = U%ymin, U%ymax
        do i = U%xmin, U%xmax
           if (temp6(1)%values(i,j,k).lt.-0.1) temp6(1)%values(i,j,k) = -0.1
           if (temp6(1)%values(i,j,k).gt.0.1) temp6(1)%values(i,j,k) = 0.1
        end do
     end do
  end do
  ie_gm = computeFieldAvg( temp6(1) , spec_rank  )
  if (spec_rank.eq.0) print*,' Wang model, mean of coef1 - after clipping =', ie_gm
  
  ! Compute Coef2 (for other part)
  temp5(3)%values =  (Li(1)%values*Mi(1)%values + Li(2)%values*Mi(2)%values + Li(3)%values*Mi(3)%values)*(Mi(1)%values*Qi(1)%values + Mi(2)%values*Qi(2)%values + Mi(3)%values*Qi(3)%values) &
               &  -  (Li(1)%values*Qi(1)%values + Li(2)%values*Qi(2)%values + Li(3)%values*Qi(3)%values)*(Mi(1)%values*Mi(1)%values + Mi(2)%values*Mi(2)%values + Mi(3)%values*Mi(3)%values) 
  temp5(2)%values =  (Mi(1)%values*Qi(1)%values + Mi(2)%values*Qi(2)%values + Mi(3)%values*Qi(3)%values)*(Qi(1)%values*Mi(1)%values + Qi(2)%values*Mi(2)%values + Qi(3)%values*Mi(3)%values) &
               &  -  (Mi(1)%values*Mi(1)%values + Mi(2)%values*Mi(2)%values + Mi(3)%values*Mi(3)%values)*(Qi(1)%values*Qi(1)%values + Qi(2)%values*Qi(2)%values + Qi(3)%values*Qi(3)%values)
  !ie_gm = computeFieldAvg( temp5(3) , spec_rank  )
  !qe_gm = computeFieldAvg( temp5(2) , spec_rank  )
  !Coef2 = ie_gm / qe_gm
  !if (spec_rank.eq.0) print*,' Wang model, coef2 =', Coef2
  temp6(2)%values = temp5(3)%values / temp5(2)%values
  ie_gm = computeFieldAvg( temp6(2) , spec_rank  )
  if (spec_rank.eq.0) print*,' Wang model, mean of coef2 - before clipping =', ie_gm
  do k = U%zmin, U%zmax
     do j = U%ymin, U%ymax
        do i = U%xmin, U%xmax
           if (temp6(2)%values(i,j,k).lt.-0.1) temp6(2)%values(i,j,k) = -0.1
           if (temp6(2)%values(i,j,k).gt.0.1) temp6(2)%values(i,j,k) = 0.1
        end do
     end do
  end do
  ie_gm = computeFieldAvg( temp6(2) , spec_rank  )
  if (spec_rank.eq.0) print*,' Wang model, mean of coef1 - after clipping =', ie_gm
  
  ! Smag part
  temp1(1)%values = temp6(1)%values *real(filtersize)**2.0 * ( temp5(1)%values * temp4(1)%values )
  temp1(2)%values = temp6(1)%values *real(filtersize)**2.0 * ( temp5(1)%values * temp4(2)%values )
  temp1(3)%values = temp6(1)%values *real(filtersize)**2.0 * ( temp5(1)%values * temp4(3)%values )
  
  ! Pang part
  Pi(1)%values = S11%values*temp4(1)%values + S12%values*temp4(2)%values + S13%values*temp4(3)%values 
  Pi(2)%values = S12%values*temp4(1)%values + S22%values*temp4(2)%values + S23%values*temp4(3)%values 
  Pi(3)%values = S13%values*temp4(1)%values + S23%values*temp4(2)%values + S33%values*temp4(3)%values
  Pi(1)%values = temp6(2)%values * (filtersize)**2.0 * Pi(1)%values
  Pi(2)%values = temp6(2)%values * (filtersize)**2.0 * Pi(2)%values
  Pi(3)%values = temp6(2)%values * (filtersize)**2.0 * Pi(3)%values
  
  ! Model
  temp1(1)%values = temp1(1)%values + Pi(1)%values 
  temp1(2)%values = temp1(2)%values + Pi(2)%values 
  temp1(3)%values = temp1(3)%values + Pi(3)%values
  
  call ftran(temp1(1),tempK(1),res)
  call ftran(temp1(2),tempK(2),res)
  call ftran(temp1(3),tempK(3),res)

  call computeDivergenceK(ScalWN,tempK(1),tempK(2),tempK(3),div)

  call deletedatalayout(S11)
  call deletedatalayout(S12)
  call deletedatalayout(S13)
  call deletedatalayout(S22)
  call deletedatalayout(S23)
  call deletedatalayout(S33)
  
  call deletedatalayout(Scal_ft)
  call deletedatalayout(S11_ft)
  call deletedatalayout(S12_ft)
  call deletedatalayout(S13_ft)
  call deletedatalayout(S22_ft)
  call deletedatalayout(S23_ft)
  call deletedatalayout(S33_ft)
  
  
  res = deleteWorkArray(temp)
  res = deleteWorkArray(temp1)
  res = deleteWorkArray(temp2)
  res = deleteWorkArray(temp3)
  res = deleteWorkArray(temp4)
  res = deleteWorkArray(temp5)
  res = deleteWorkArray(temp6)
  res = deleteWorkArray(tempK)
  res = deleteWorkArray(tempKf)
  res = deleteWorkArray(Li)
  res = deleteWorkArray(Mi)
  res = deleteWorkArray(Qi)
  res = deleteWorkArray(Pi)
  res = deleteWorkArray(dZdxi_ft)
  res = deleteWorkArray(SdZdxi_ft)
  res = deleteWorkArray(VelScal_ft)
  res = deleteWorkArray(Vel_ft)
  res = deleteWorkArray(SmZ)
 
  end subroutine DynWangSca
 
!------------------------------------------------------------------------------
!> @author 
!> Guillaume Balarac, LEGI
!>
!> @details
!>  Dynamic Regularized Gradient model (based on Dynamic Procedure 4 - RGM-VII)
!> @param [in]    U  longitudinal velocity in physical space
!> @param [in]    Uk longitudinal velocity in spectral space
!> @param [in]    Vk velocity in spectral space
!> @param [in]    Wk velocity in spectral space
!> @param [in]    Scalk scalar in spectral space 
!> @param [in]    ScalWN wave numbers for scalar 
!> @param [in]    spec_rank is the rank of processus
!> @param [inout] res logical
!> @param [out]   div is the divergence of SGS model
!------------------------------------------------------------------------------
 subroutine DynRGMSca(div,  U,V,W,Uk,VK,Wk,Scal,ScalK,ScalWN,filter,type_avg,spec_rank,res)




  !I/O
  type(complex_data_layout),intent(inout) :: div
  type(complex_data_layout),intent(in)    :: Uk,Vk,Wk,ScalK
  type(real_data_layout),intent(in)       :: U, V, W, Scal
  type(WaveNumbers), intent(in)           :: ScalWN
  logical, intent(inout)                  :: res
  integer, intent(in)                     :: filter, spec_rank, type_avg

  !Local data
  type(real_data_layout),dimension(:),allocatable    :: temp,temp1, temp2, temp3, temp4, Li, Mi
  type(real_data_layout),dimension(:),allocatable    :: VelScal_ft, Vel_ft, dZdxi_ft, SmZ
  type(real_data_layout)                             :: S11, S12, S13, S22, S23, S33
  type(real_data_layout)                             :: S11_ft, S12_ft, S13_ft, S22_ft, S23_ft, S33_ft, Scal_ft
  type(complex_data_layout),dimension(:),allocatable :: tempK, tempKf
  real(WP)                                           :: filtersize, DynCoef, num, den
  real(WP), dimension(:),allocatable                 :: DynCoef1D, num1D, den1D
  
  !! for eigenvalues computation 
  integer :: a, b, info, i,j, k
  real(WP), dimension(3,3) :: tab,VecP,VL,TvecP
  real(WP), dimension(3,3) :: E1,E2,E3,M1,M2,M3
  real(WP), dimension(3) :: valR,valI
  real(WP) :: tabwork(102)

  res = .true.
  filtersize = computeDelta(U)
  
  res = initWorkArray(U,3,temp)
  res = initWorkArray(U,3,temp1)
  res = initWorkArray(U,3,temp2)
  res = initWorkArray(U,3,temp3)
  res = initWorkArray(U,3,temp4)
  res = initWorkArray(U,3,SmZ)
  res = initWorkArray(U,3,Li)
  res = initWorkArray(U,3,Mi)
  res = initWorkArray(U,3,VelScal_ft)
  res = initWorkArray(U,3,Vel_ft)
  res = initWorkArray(U,3,dZdxi_ft)
  res = initWorkArray(Uk,3,tempK)
  res = initWorkArray(Uk,4,tempKf)
  
  call computeGradientK(ScalWN,ScalK,tempK(1),tempK(2),tempK(3))
  call btran(tempK(1),temp4(1),res)
  call btran(tempK(2),temp4(2),res)
  call btran(tempK(3),temp4(3),res)
  call computeGradientK(ScalWN,UK,tempK(1),tempK(2),tempK(3))
  call btran(tempK(1),temp1(1),res)
  call btran(tempK(2),temp1(2),res)
  call btran(tempK(3),temp1(3),res)
  call computeGradientK(ScalWN,VK,tempK(1),tempK(2),tempK(3))
  call btran(tempK(1),temp2(1),res)
  call btran(tempK(2),temp2(2),res)
  call btran(tempK(3),temp2(3),res)
  call computeGradientK(ScalWN,WK,tempK(1),tempK(2),tempK(3))
  call btran(tempK(1),temp3(1),res)
  call btran(tempK(2),temp3(2),res)
  call btran(tempK(3),temp3(3),res)
  
  res = copyStructOnly(U,S11)
  res = copyStructOnly(U,S12)
  res = copyStructOnly(U,S13)
  res = copyStructOnly(U,S22)
  res = copyStructOnly(U,S23)
  res = copyStructOnly(U,S33)
  
  res = copyStructOnly(U,S11_ft)
  res = copyStructOnly(U,S12_ft)
  res = copyStructOnly(U,S13_ft)
  res = copyStructOnly(U,S22_ft)
  res = copyStructOnly(U,S23_ft)
  res = copyStructOnly(U,S33_ft)
  res = copyStructOnly(U,Scal_ft)
  
  ! compute Sij = 1/2 * ( duidxj + dujdxi)
  S11%values = temp1(1)%values
  S12%values = 0.5 * ( temp1(2)%values + temp2(1)%values )
  S13%values = 0.5 * ( temp1(3)%values + temp3(1)%values )
  S22%values = temp2(2)%values
  S23%values = 0.5 * ( temp2(3)%values + temp3(2)%values )
  S33%values = temp3(3)%values
  
  ! compute Sij-
    E1 = 0
    E2 = 0
    E3 = 0
    E1(1,1) = 1._WP
    E2(2,2) = 1._WP
    E3(3,3) = 1._WP
#ifdef LAPACK
    do k = U%zmin, U%zmax
       do j = U%ymin, U%ymax
          do i = U%xmin, U%xmax
             tab(1,1) = S11%values(i,j,k)
             tab(1,2) = S12%values(i,j,k)
             tab(1,3) = S13%values(i,j,k)
             tab(2,1) = S12%values(i,j,k)
             tab(2,2) = S22%values(i,j,k)
             tab(2,3) = S23%values(i,j,k)
             tab(3,1) = S13%values(i,j,k)
             tab(3,2) = S23%values(i,j,k)
             tab(3,3) = S33%values(i,j,k)
             call DGEEV('N','V',3,tab,3,valR,valI,VL,3,VecP,3,tabwork,102,INFO)
             do a=1,3
               do b=1,3
                 TvecP(a,b)=VecP(b,a)
               enddo
             enddo
             M1 = matmul(E1,TvecP)
             M1 = matmul(vecP,M1)
             M2 = matmul(E2,TvecP)
             M2 = matmul(vecP,M2)
             M3 = matmul(E3,TvecP)
             M3 = matmul(vecP,M3)
             SmZ(1)%values(i,j,k) = ( MIN(valR(1),0.0_WP)* ( M1(1,1)*temp4(1)%values(i,j,k) + M1(1,2)*temp4(2)%values(i,j,k) + M1(1,3)*temp4(3)%values(i,j,k)) &
                                   & + MIN(valR(2),0.0_WP)* ( M2(1,1)*temp4(1)%values(i,j,k) + M2(1,2)*temp4(2)%values(i,j,k) + M2(1,3)*temp4(3)%values(i,j,k)) &
                                   & + MIN(valR(3),0.0_WP)* ( M3(1,1)*temp4(1)%values(i,j,k) + M3(1,2)*temp4(2)%values(i,j,k) + M3(1,3)*temp4(3)%values(i,j,k)))

             SmZ(2)%values(i,j,k) = ( MIN(valR(1),0.0_WP)* ( M1(2,1)*temp4(1)%values(i,j,k) + M1(2,2)*temp4(2)%values(i,j,k) + M1(2,3)*temp4(3)%values(i,j,k)) &
                                   & + MIN(valR(2),0.0_WP)* ( M2(2,1)*temp4(1)%values(i,j,k) + M2(2,2)*temp4(2)%values(i,j,k) + M2(2,3)*temp4(3)%values(i,j,k)) &
                                   & + MIN(valR(3),0.0_WP)* ( M3(2,1)*temp4(1)%values(i,j,k) + M3(2,2)*temp4(2)%values(i,j,k) + M3(2,3)*temp4(3)%values(i,j,k)))

             SmZ(3)%values(i,j,k) = ( MIN(valR(1),0.0_WP)* ( M1(3,1)*temp4(1)%values(i,j,k) + M1(3,2)*temp4(2)%values(i,j,k) + M1(3,3)*temp4(3)%values(i,j,k)) &
                                   & + MIN(valR(2),0.0_WP)* ( M2(3,1)*temp4(1)%values(i,j,k) + M2(3,2)*temp4(2)%values(i,j,k) + M2(3,3)*temp4(3)%values(i,j,k)) &
                                   & + MIN(valR(3),0.0_WP)* ( M3(3,1)*temp4(1)%values(i,j,k) + M3(3,2)*temp4(2)%values(i,j,k) + M3(3,3)*temp4(3)%values(i,j,k)))
          end do
       end do
    end do
#endif

    ! Dynamic procedure 4 (cf. RGM-(VII) ) ::
    ! C * (2.0*filtersize)^2 * Sij^-_ft * dZdxj_ft * dZdxi_ft = Li * dZdxi_ft
    !
    !  C =  <  Li * dZdxi_ft  >  /  <  (2.0*filtersize)^2 * Sij^-_ft * dZdxj_ft * dZdxi_ft  >
    ! 
    ! Li * dZdxi_ft   <-- Li(1)
    ! (2.0*filtersize)^2 * Sij^-_ft * dZdxj_ft * dZdxi_ft    <-- Mi(1)
    !
    !     C = < Li(1) > / < Mi(1) >
    !
    ! Compute Li = TestFilter (uiz) - TestFilter(ui)*TestFilter(z)
    ! Compute TestFilter(ui) and TestFilter(z)
    call computeFilter(ScalWN,2.0*filtersize,Uk,tempKf(1),filter) !U bar
    call computeFilter(ScalWN,2.0*filtersize,VK,tempKf(2),filter) !V bar
    call computeFilter(ScalWN,2.0*filtersize,WK,tempKf(3),filter) !W bar
    call computeFilter(ScalWN,2.0*filtersize,ScalK,tempKf(4),filter) !Z bar
    call btran(tempKf(1),Vel_ft(1),res) !U bar
    call btran(tempKf(2),Vel_ft(2),res) !V bar
    call btran(tempKf(3),Vel_ft(3),res) !W bar
    call btran(tempKf(4),Scal_ft,res) !Z bar
    VelScal_ft(1)%values = U%values * Scal%values
    VelScal_ft(2)%values = V%values * Scal%values
    VelScal_ft(3)%values = W%values * Scal%values
    call ftran(VelScal_ft(1),tempK(1),res)
    call ftran(VelScal_ft(2),tempK(2),res)
    call ftran(VelScal_ft(3),tempK(3),res)
    call computeFilter(ScalWN,2.0*filtersize,tempK(1),tempKf(1),filter) !U bar
    call computeFilter(ScalWN,2.0*filtersize,tempK(2),tempKf(2),filter) !V bar
    call computeFilter(ScalWN,2.0*filtersize,tempK(3),tempKf(3),filter) !W bar
    call btran(tempKf(1),VelScal_ft(1),res) !U bar
    call btran(tempKf(2),VelScal_ft(2),res) !V bar
    call btran(tempKf(3),VelScal_ft(3),res) !W bar
    ! Compute Li
    Li(1)%values = VelScal_ft(1)%values - Vel_ft(1)%values * Scal_ft%values
    Li(2)%values = VelScal_ft(2)%values - Vel_ft(2)%values * Scal_ft%values
    Li(3)%values = VelScal_ft(3)%values - Vel_ft(3)%values * Scal_ft%values
    ! Compuet dZdxi_ft
    call ftran(temp4(1),tempK(1),res)
    call ftran(temp4(2),tempK(2),res)
    call ftran(temp4(3),tempK(3),res)
    call computeFilter(ScalWN,2.0*filtersize,tempK(1),tempKf(1),filter) !U bar
    call computeFilter(ScalWN,2.0*filtersize,tempK(2),tempKf(2),filter) !V bar
    call computeFilter(ScalWN,2.0*filtersize,tempK(3),tempKf(3),filter) !W bar
    call btran(tempKf(1),dZdxi_ft(1),res) !U bar
    call btran(tempKf(2),dZdxi_ft(2),res) !V bar
    call btran(tempKf(3),dZdxi_ft(3),res) !W bar
    ! Li(1) computation  <-- Li * dZdxi_ft 
    Li(1)%values = Li(1)%values*dZdxi_ft(1)%values + Li(2)%values*dZdxi_ft(2)%values + Li(3)%values*dZdxi_ft(3)%values
    
    ! Mi(1) computation
    ! Sij^-_ft * dZdxj_ft <-- temp
    call ftran(S11,tempK(1),res)
    call ftran(S12,tempK(2),res)
    call ftran(S13,tempK(3),res)
    call computeFilter(ScalWN,2.0*filtersize,tempK(1),tempKf(1),filter) !U bar
    call computeFilter(ScalWN,2.0*filtersize,tempK(2),tempKf(2),filter) !V bar
    call computeFilter(ScalWN,2.0*filtersize,tempK(3),tempKf(3),filter) !W bar
    call btran(tempKf(1),S11_ft,res) !U bar
    call btran(tempKf(2),S12_ft,res) !V bar
    call btran(tempKf(3),S13_ft,res) !W bar
    call ftran(S22,tempK(1),res)
    call ftran(S23,tempK(2),res)
    call ftran(S33,tempK(3),res)
    call computeFilter(ScalWN,2.0*filtersize,tempK(1),tempKf(1),filter) !U bar
    call computeFilter(ScalWN,2.0*filtersize,tempK(2),tempKf(2),filter) !V bar
    call computeFilter(ScalWN,2.0*filtersize,tempK(3),tempKf(3),filter) !W bar
    call btran(tempKf(1),S22_ft,res) !U bar
    call btran(tempKf(2),S23_ft,res) !V bar
    call btran(tempKf(3),S33_ft,res) !W bar
    E1 = 0
    E2 = 0
    E3 = 0
    E1(1,1) = 1._WP
    E2(2,2) = 1._WP
    E3(3,3) = 1._WP
#ifdef LAPACK
    do k = U%zmin, U%zmax
       do j = U%ymin, U%ymax
          do i = U%xmin, U%xmax
             tab(1,1) = S11_ft%values(i,j,k)
             tab(1,2) = S12_ft%values(i,j,k)
             tab(1,3) = S13_ft%values(i,j,k)
             tab(2,1) = S12_ft%values(i,j,k)
             tab(2,2) = S22_ft%values(i,j,k)
             tab(2,3) = S23_ft%values(i,j,k)
             tab(3,1) = S13_ft%values(i,j,k)
             tab(3,2) = S23_ft%values(i,j,k)
             tab(3,3) = S33_ft%values(i,j,k)
             call DGEEV('N','V',3,tab,3,valR,valI,VL,3,VecP,3,tabwork,102,INFO)
             do a=1,3
               do b=1,3
                 TvecP(a,b)=VecP(b,a)
               enddo
             enddo
             M1 = matmul(E1,TvecP)
             M1 = matmul(vecP,M1)
             M2 = matmul(E2,TvecP)
             M2 = matmul(vecP,M2)
             M3 = matmul(E3,TvecP)
             M3 = matmul(vecP,M3)
             temp(1)%values(i,j,k) = ( MIN(valR(1),0.0_WP)* ( M1(1,1)*dZdxi_ft(1)%values(i,j,k) + M1(1,2)*dZdxi_ft(2)%values(i,j,k) + M1(1,3)*dZdxi_ft(3)%values(i,j,k)) &
                                   & + MIN(valR(2),0.0_WP)* ( M2(1,1)*dZdxi_ft(1)%values(i,j,k) + M2(1,2)*dZdxi_ft(2)%values(i,j,k) + M2(1,3)*dZdxi_ft(3)%values(i,j,k)) &
                                   & + MIN(valR(3),0.0_WP)* ( M3(1,1)*dZdxi_ft(1)%values(i,j,k) + M3(1,2)*dZdxi_ft(2)%values(i,j,k) + M3(1,3)*dZdxi_ft(3)%values(i,j,k)))

             temp(2)%values(i,j,k) = ( MIN(valR(1),0.0_WP)* ( M1(2,1)*dZdxi_ft(1)%values(i,j,k) + M1(2,2)*dZdxi_ft(2)%values(i,j,k) + M1(2,3)*dZdxi_ft(3)%values(i,j,k)) &
                                   & + MIN(valR(2),0.0_WP)* ( M2(2,1)*dZdxi_ft(1)%values(i,j,k) + M2(2,2)*dZdxi_ft(2)%values(i,j,k) + M2(2,3)*dZdxi_ft(3)%values(i,j,k)) &
                                   & + MIN(valR(3),0.0_WP)* ( M3(2,1)*dZdxi_ft(1)%values(i,j,k) + M3(2,2)*dZdxi_ft(2)%values(i,j,k) + M3(2,3)*dZdxi_ft(3)%values(i,j,k)))

             temp(3)%values(i,j,k) = ( MIN(valR(1),0.0_WP)* ( M1(3,1)*dZdxi_ft(1)%values(i,j,k) + M1(3,2)*dZdxi_ft(2)%values(i,j,k) + M1(3,3)*dZdxi_ft(3)%values(i,j,k)) &
                                   & + MIN(valR(2),0.0_WP)* ( M2(3,1)*dZdxi_ft(1)%values(i,j,k) + M2(3,2)*dZdxi_ft(2)%values(i,j,k) + M2(3,3)*dZdxi_ft(3)%values(i,j,k)) &
                                   & + MIN(valR(3),0.0_WP)* ( M3(3,1)*dZdxi_ft(1)%values(i,j,k) + M3(3,2)*dZdxi_ft(2)%values(i,j,k) + M3(3,3)*dZdxi_ft(3)%values(i,j,k)))

          end do
       end do
    end do
#endif
    
    Mi(1)%values = (2.0*filtersize)**2.0 * ( temp(1)%values * dZdxi_ft(1)%values + &
                                             temp(2)%values * dZdxi_ft(2)%values + &
                                             temp(3)%values * dZdxi_ft(3)%values )
                 
    ! Compute DynCoef = <Li(1)> / <Mi(1)>
    if (type_avg.eq.0) then
    
       num = computeFieldAvg( Li(1) , spec_rank )
       den = computeFieldAvg( Mi(1) , spec_rank )
       DynCoef = num / den
       if (spec_rank.eq.0) print*,'DynCoef Sij^- (procedure 4): ',DynCoef

       temp1(1)%values = DynCoef * filtersize**2.0 * SmZ(1)%values
       temp1(2)%values = DynCoef * filtersize**2.0 * SmZ(2)%values
       temp1(3)%values = DynCoef * filtersize**2.0 * SmZ(3)%values
       
    else if (type_avg.eq.1) then
    
       allocate(num1D(U%Nz))
       allocate(den1D(U%Nz))
       allocate(Dyncoef1D(U%Nz))
    
       if(.not. computeFieldAvg2D(Li(1),num1D,spec_rank)) write(*,'(a)')&
        & '[ERROR] Could not compute AVG for Li'

        if(.not. computeFieldAvg2D(Mi(1),den1D,spec_rank)) write(*,'(a)')&
        & '[ERROR] Could not compute AVG for Mi'

       Dyncoef1D(:) = num1D(:) / den1D(:)
       
       do k = U%zmin, U%zmax
          temp1(1)%values(:,:,k) = Dyncoef1D(k) * filtersize**2.0 * SmZ(1)%values(:,:,k)
          temp1(2)%values(:,:,k) = Dyncoef1D(k) * filtersize**2.0 * SmZ(2)%values(:,:,k)
          temp1(3)%values(:,:,k) = Dyncoef1D(k) * filtersize**2.0 * SmZ(3)%values(:,:,k)
       end do
       
       deallocate(num1D)
       deallocate(den1D)
       deallocate(Dyncoef1D)
       
    else
    
       print*,'[ERROR] AVERGAGING NOT KNOWN FOR DYNAMIC PROCEDURE'
       
    end if
       

  
  temp(1)%values = temp1(1)%values * temp4(1)%values + temp1(2)%values * temp4(2)%values + temp1(3)%values * temp4(3)%values
  !do k = U%zmin, U%zmax
  !   do j = U%ymin, U%ymax
  !      do i = U%xmin, U%xmax
  !         if (temp(1)%values(i,j,k).gt.0.0_WP) then
  !            print*,'DynRGM something wrong', temp(1)%values(i,j,k)
  !         end if
  !      end do
  !   end do
  !end do
  DynCoef = computeFieldAvg( temp(1) , spec_rank )
  if (spec_rank.eq.0) print*,'DynRGM Dissipation =',DynCoef
  
  call ftran(temp1(1),tempK(1),res)
  call ftran(temp1(2),tempK(2),res)
  call ftran(temp1(3),tempK(3),res)

  call computeDivergenceK(ScalWN,tempK(1),tempK(2),tempK(3),div)
  
  call deletedatalayout(S11)
  call deletedatalayout(S12)
  call deletedatalayout(S13)
  call deletedatalayout(S22)
  call deletedatalayout(S23)
  call deletedatalayout(S33)
  
  call deletedatalayout(Scal_ft)
  call deletedatalayout(S11_ft)
  call deletedatalayout(S12_ft)
  call deletedatalayout(S13_ft)
  call deletedatalayout(S22_ft)
  call deletedatalayout(S23_ft)
  call deletedatalayout(S33_ft)
  
  
  res = deleteWorkArray(temp)
  res = deleteWorkArray(temp1)
  res = deleteWorkArray(temp2)
  res = deleteWorkArray(temp3)
  res = deleteWorkArray(temp4)
  res = deleteWorkArray(SmZ)
  res = deleteWorkArray(Li)
  res = deleteWorkArray(Mi)
  res = deleteWorkArray(VelScal_ft)
  res = deleteWorkArray(Vel_ft)
  res = deleteWorkArray(dZdxi_ft)
  res = deleteWorkArray(tempK)
  res = deleteWorkArray(tempKf)
  
  

 end subroutine DynRGMSca
 
 !------------------------------------------------------------------------------
!> @author 
!> Guillaume Balarac, LEGI
!>
!> @details
!>  Regularized Gradient model 
!> @param [in]    U  longitudinal velocity in physical space
!> @param [in]    Uk longitudinal velocity in spectral space
!> @param [in]    Vk velocity in spectral space
!> @param [in]    Wk velocity in spectral space
!> @param [in]    Scalk scalar in spectral space 
!> @param [in]    ScalWN wave numbers for scalar 
!> @param [in]    spec_rank is the rank of processus
!> @param [inout] res logical
!> @param [out]   div is the divergence of SGS model
!------------------------------------------------------------------------------
 subroutine RGMOSca(div,U,Uk,Vk,Wk,ScalK,ScalWN,res,spec_rank,coefRGM)


  !I/O
  type(complex_data_layout),intent(in)    :: Uk,Vk,Wk,ScalK 
  type(complex_data_layout),intent(inout) :: div
  type(real_data_layout),intent(in)       :: U
  type(WaveNumbers), intent(in)           :: ScalWN
  logical, intent(inout)                  :: res
  real(WP), intent(in)                    :: coefRGM
  integer, intent(in)                     :: spec_rank

  !Local data
  type(real_data_layout),dimension(:),allocatable    :: temp,temp1, temp2, temp3, temp4
  type(real_data_layout)                             :: S11, S12, S13, S22, S23, S33, O12, O13, O23
  type(complex_data_layout),dimension(:),allocatable :: tempK
  real(WP)                                           :: Cdeltx2
  
  !! for eigenvalues computation 
  integer :: a, b, info, i,j, k
  real(WP), dimension(3,3) :: tab,VecP,VL,TvecP
  real(WP), dimension(3,3) :: E1,E2,E3,M1,M2,M3
  real(WP), dimension(3) :: valR,valI
  real(WP) :: tabwork(102)

  res = .true.
  Cdeltx2 = coefRGM * computeDelta(U)**2.0
  
  res = initWorkArray(U,3,temp)
  res = initWorkArray(U,3,temp1)
  res = initWorkArray(U,3,temp2)
  res = initWorkArray(U,3,temp3)
  res = initWorkArray(U,3,temp4)
  res = initWorkArray(Uk,3,tempK)
  
  call computeGradientK(ScalWN,ScalK,tempK(1),tempK(2),tempK(3))
  call btran(tempK(1),temp4(1),res)
  call btran(tempK(2),temp4(2),res)
  call btran(tempK(3),temp4(3),res)
  call computeGradientK(ScalWN,UK,tempK(1),tempK(2),tempK(3))
  call btran(tempK(1),temp1(1),res)
  call btran(tempK(2),temp1(2),res)
  call btran(tempK(3),temp1(3),res)
  call computeGradientK(ScalWN,VK,tempK(1),tempK(2),tempK(3))
  call btran(tempK(1),temp2(1),res)
  call btran(tempK(2),temp2(2),res)
  call btran(tempK(3),temp2(3),res)
  call computeGradientK(ScalWN,WK,tempK(1),tempK(2),tempK(3))
  call btran(tempK(1),temp3(1),res)
  call btran(tempK(2),temp3(2),res)
  call btran(tempK(3),temp3(3),res)
  
  res = copyStructOnly(U,S11)
  res = copyStructOnly(U,S12)
  res = copyStructOnly(U,S13)
  res = copyStructOnly(U,S22)
  res = copyStructOnly(U,S23)
  res = copyStructOnly(U,S33)
  res = copyStructOnly(U,O12)
  res = copyStructOnly(U,O13)
  res = copyStructOnly(U,O23)
  
  ! compute Sij = 1/2 * ( duidxj + dujdxi)
  S11%values = temp1(1)%values
  S12%values = 0.5 * ( temp1(2)%values + temp2(1)%values )
  S13%values = 0.5 * ( temp1(3)%values + temp3(1)%values )
  S22%values = temp2(2)%values
  S23%values = 0.5 * ( temp2(3)%values + temp3(2)%values )
  S33%values = temp3(3)%values
  
  O12%values = 0.5 * ( temp1(2)%values - temp2(1)%values )
  O13%values = 0.5 * ( temp1(3)%values - temp3(1)%values )
  O23%values = 0.5 * ( temp2(3)%values - temp3(2)%values )
  
  ! compute Sij-
    E1 = 0
    E2 = 0
    E3 = 0
    E1(1,1) = 1._WP
    E2(2,2) = 1._WP
    E3(3,3) = 1._WP
#ifdef LAPACK
    do k = U%zmin, U%zmax
       do j = U%ymin, U%ymax
          do i = U%xmin, U%xmax
             tab(1,1) = S11%values(i,j,k)
             tab(1,2) = S12%values(i,j,k)
             tab(1,3) = S13%values(i,j,k)
             tab(2,1) = S12%values(i,j,k)
             tab(2,2) = S22%values(i,j,k)
             tab(2,3) = S23%values(i,j,k)
             tab(3,1) = S13%values(i,j,k)
             tab(3,2) = S23%values(i,j,k)
             tab(3,3) = S33%values(i,j,k)
             call DGEEV('N','V',3,tab,3,valR,valI,VL,3,VecP,3,tabwork,102,INFO)
             do a=1,3
               do b=1,3
                 TvecP(a,b)=VecP(b,a)
               enddo
             enddo
             M1 = matmul(E1,TvecP)
             M1 = matmul(vecP,M1)
             M2 = matmul(E2,TvecP)
             M2 = matmul(vecP,M2)
             M3 = matmul(E3,TvecP)
             M3 = matmul(vecP,M3)
             temp(1)%values(i,j,k) = ( MIN(valR(1),0.0_WP)* ( M1(1,1)*temp4(1)%values(i,j,k) + M1(1,2)*temp4(2)%values(i,j,k) + M1(1,3)*temp4(3)%values(i,j,k)) &
                                   & + MIN(valR(2),0.0_WP)* ( M2(1,1)*temp4(1)%values(i,j,k) + M2(1,2)*temp4(2)%values(i,j,k) + M2(1,3)*temp4(3)%values(i,j,k)) &
                                   & + MIN(valR(3),0.0_WP)* ( M3(1,1)*temp4(1)%values(i,j,k) + M3(1,2)*temp4(2)%values(i,j,k) + M3(1,3)*temp4(3)%values(i,j,k)))

             temp(2)%values(i,j,k) = ( MIN(valR(1),0.0_WP)* ( M1(2,1)*temp4(1)%values(i,j,k) + M1(2,2)*temp4(2)%values(i,j,k) + M1(2,3)*temp4(3)%values(i,j,k)) &
                                   & + MIN(valR(2),0.0_WP)* ( M2(2,1)*temp4(1)%values(i,j,k) + M2(2,2)*temp4(2)%values(i,j,k) + M2(2,3)*temp4(3)%values(i,j,k)) &
                                   & + MIN(valR(3),0.0_WP)* ( M3(2,1)*temp4(1)%values(i,j,k) + M3(2,2)*temp4(2)%values(i,j,k) + M3(2,3)*temp4(3)%values(i,j,k)))

             temp(3)%values(i,j,k) = ( MIN(valR(1),0.0_WP)* ( M1(3,1)*temp4(1)%values(i,j,k) + M1(3,2)*temp4(2)%values(i,j,k) + M1(3,3)*temp4(3)%values(i,j,k)) &
                                   & + MIN(valR(2),0.0_WP)* ( M2(3,1)*temp4(1)%values(i,j,k) + M2(3,2)*temp4(2)%values(i,j,k) + M2(3,3)*temp4(3)%values(i,j,k)) &
                                   & + MIN(valR(3),0.0_WP)* ( M3(3,1)*temp4(1)%values(i,j,k) + M3(3,2)*temp4(2)%values(i,j,k) + M3(3,3)*temp4(3)%values(i,j,k)))
          end do
       end do
    end do
#endif

  temp1(1)%values = Cdeltx2 * (temp(1)%values + O12%values*temp4(2)%values + O13%values*temp4(3)%values)
  temp1(2)%values = Cdeltx2 * (temp(2)%values - O12%values*temp4(1)%values + O23%values*temp4(3)%values)
  temp1(3)%values = Cdeltx2 * (temp(3)%values - O13%values*temp4(1)%values - O23%values*temp4(2)%values)
  
  
  temp(1)%values = temp1(1)%values * temp4(1)%values + temp1(2)%values * temp4(2)%values + temp1(3)%values * temp4(3)%values
  do k = U%zmin, U%zmax
     do j = U%ymin, U%ymax
        do i = U%xmin, U%xmax
           if (temp(1)%values(i,j,k).gt.0.0_WP) then
              print*,'RGMO something wrong', temp(1)%values(i,j,k)
           end if
        end do
     end do
  end do
  Cdeltx2 = computeFieldAvg( temp(1) , spec_rank )
  if (spec_rank.eq.0) print*,'RGMO Dissipation =',Cdeltx2
  
  call ftran(temp1(1),tempK(1),res)
  call ftran(temp1(2),tempK(2),res)
  call ftran(temp1(3),tempK(3),res)

  call computeDivergenceK(ScalWN,tempK(1),tempK(2),tempK(3),div)

  call deletedatalayout(S11)
  call deletedatalayout(S12)
  call deletedatalayout(S13)
  call deletedatalayout(S22)
  call deletedatalayout(S23)
  call deletedatalayout(S33)
  call deletedatalayout(O12)
  call deletedatalayout(O13)
  call deletedatalayout(O23)
  
  res = deleteWorkArray(temp)
  res = deleteWorkArray(temp1)
  res = deleteWorkArray(temp2)
  res = deleteWorkArray(temp3)
  res = deleteWorkArray(temp4)
  res = deleteWorkArray(tempK)

 end subroutine RGMOSca


!------------------------------------------------------------------------------
!> @author 
!> Guillaume Balarac, LEGI
!>
!> @details
!>  Regularized Gradient model 
!> @param [in]    U  longitudinal velocity in physical space
!> @param [in]    Uk longitudinal velocity in spectral space
!> @param [in]    Vk velocity in spectral space
!> @param [in]    Wk velocity in spectral space
!> @param [in]    Scalk scalar in spectral space 
!> @param [in]    ScalWN wave numbers for scalar 
!> @param [in]    spec_rank is the rank of processus
!> @param [inout] res logical
!> @param [out]   div is the divergence of SGS model
!------------------------------------------------------------------------------
 subroutine RGMScaForTest(div,U,Uk,Vk,Wk,ScalK,ScalWN,res,spec_rank,coefRGM,Ti)

  use stat_tools , only : equalToVal,computeMSE

  !I/O
  type(complex_data_layout),intent(in)    :: Uk,Vk,Wk,ScalK 
  type(complex_data_layout),intent(inout) :: div
  type(real_data_layout),intent(in)       :: U
  type(real_data_layout),dimension(:),intent(inout),optional :: Ti 
  type(WaveNumbers), intent(in)           :: ScalWN
  logical, intent(inout)                  :: res
  real(WP), intent(in)                 :: coefRGM
  integer, intent (in)                    :: spec_rank

  !Local data
  type(real_data_layout),dimension(:),allocatable    :: temp,temp1, temp2, temp3, temp4
  type(real_data_layout)                             :: temp5,temp6,temp7
  type(real_data_layout)                             :: S11, S12, S13, S22, S23, S33
  type(complex_data_layout),dimension(:),allocatable :: tempK
  real(WP)                                           :: Cdeltx2
  real(WP)                                           :: S11out,S12out,S13out,S22out,S23out,S33out 
  real(WP)                                           :: lam1,lam2,lam3 
  
  !! for eigenvalues computation 
  integer :: a, b, info, i,j, k
  real(WP) :: t1,t2
  real(WP), dimension(3,3) :: tab,VecP,VL,TvecP
  real(WP), dimension(3,3) :: E1,E2,E3,M1,M2,M3
  real(WP), dimension(3) :: valR,valI
  real(WP) :: tabwork(102)

  res = .true.
  Cdeltx2 = coefRGM * computeDelta(U)**2.0
  
  res = initWorkArray(U,3,temp)
  res = initWorkArray(U,3,temp1)
  res = initWorkArray(U,3,temp2)
  res = initWorkArray(U,3,temp3)
  res = initWorkArray(U,3,temp4)
  res = initWorkArray(Uk,3,tempK)
 
  res = copyStructOnly(U,temp5)
  res = copyStructOnly(U,temp6)
  res = copyStructOnly(U,temp7)
 
  call computeGradientK(ScalWN,ScalK,tempK(1),tempK(2),tempK(3))
  call btran(tempK(1),temp4(1),res)
  call btran(tempK(2),temp4(2),res)
  call btran(tempK(3),temp4(3),res)
  call computeGradientK(ScalWN,UK,tempK(1),tempK(2),tempK(3))
  call btran(tempK(1),temp1(1),res)
  call btran(tempK(2),temp1(2),res)
  call btran(tempK(3),temp1(3),res)
  call computeGradientK(ScalWN,VK,tempK(1),tempK(2),tempK(3))
  call btran(tempK(1),temp2(1),res)
  call btran(tempK(2),temp2(2),res)
  call btran(tempK(3),temp2(3),res)
  call computeGradientK(ScalWN,WK,tempK(1),tempK(2),tempK(3))
  call btran(tempK(1),temp3(1),res)
  call btran(tempK(2),temp3(2),res)
  call btran(tempK(3),temp3(3),res)
  
  res = copyStructOnly(U,S11)
  res = copyStructOnly(U,S12)
  res = copyStructOnly(U,S13)
  res = copyStructOnly(U,S22)
  res = copyStructOnly(U,S23)
  res = copyStructOnly(U,S33)
  
  ! compute Sij = 1/2 * ( duidxj + dujdxi)
  S11%values = temp1(1)%values
  S12%values = 0.5 * ( temp1(2)%values + temp2(1)%values )
  S13%values = 0.5 * ( temp1(3)%values + temp3(1)%values )
  S22%values = temp2(2)%values
  S23%values = 0.5 * ( temp2(3)%values + temp3(2)%values )
  S33%values = temp3(3)%values

call cpu_time(t1)  
  ! compute Sij-
    E1 = 0
    E2 = 0
    E3 = 0
    E1(1,1) = 1._WP
    E2(2,2) = 1._WP
    E3(3,3) = 1._WP
#ifdef LAPACK
    do k = U%zmin, U%zmax
       do j = U%ymin, U%ymax
          do i = U%xmin, U%xmax
             tab(1,1) = S11%values(i,j,k)
             tab(1,2) = S12%values(i,j,k)
             tab(1,3) = S13%values(i,j,k)
             tab(2,1) = S12%values(i,j,k)
             tab(2,2) = S22%values(i,j,k)
             tab(2,3) = S23%values(i,j,k)
             tab(3,1) = S13%values(i,j,k)
             tab(3,2) = S23%values(i,j,k)
             tab(3,3) = S33%values(i,j,k)
             call DGEEV('N','V',3,tab,3,valR,valI,VL,3,VecP,3,tabwork,102,INFO)
             do a=1,3
               do b=1,3
                 TvecP(a,b)=VecP(b,a)
               enddo
             enddo
             M1 = matmul(E1,TvecP)
             M1 = matmul(vecP,M1)
             M2 = matmul(E2,TvecP)
             M2 = matmul(vecP,M2)
             M3 = matmul(E3,TvecP)
             M3 = matmul(vecP,M3)
             temp(1)%values(i,j,k) = ( MIN(valR(1),0.0_WP)* ( M1(1,1)*temp4(1)%values(i,j,k) + M1(1,2)*temp4(2)%values(i,j,k) + M1(1,3)*temp4(3)%values(i,j,k)) &
                                   & + MIN(valR(2),0.0_WP)* ( M2(1,1)*temp4(1)%values(i,j,k) + M2(1,2)*temp4(2)%values(i,j,k) + M2(1,3)*temp4(3)%values(i,j,k)) &
                                   & + MIN(valR(3),0.0_WP)* ( M3(1,1)*temp4(1)%values(i,j,k) + M3(1,2)*temp4(2)%values(i,j,k) + M3(1,3)*temp4(3)%values(i,j,k)))

             temp(2)%values(i,j,k) = ( MIN(valR(1),0.0_WP)* ( M1(2,1)*temp4(1)%values(i,j,k) + M1(2,2)*temp4(2)%values(i,j,k) + M1(2,3)*temp4(3)%values(i,j,k)) &
                                   & + MIN(valR(2),0.0_WP)* ( M2(2,1)*temp4(1)%values(i,j,k) + M2(2,2)*temp4(2)%values(i,j,k) + M2(2,3)*temp4(3)%values(i,j,k)) &
                                   & + MIN(valR(3),0.0_WP)* ( M3(2,1)*temp4(1)%values(i,j,k) + M3(2,2)*temp4(2)%values(i,j,k) + M3(2,3)*temp4(3)%values(i,j,k)))

             temp(3)%values(i,j,k) = ( MIN(valR(1),0.0_WP)* ( M1(3,1)*temp4(1)%values(i,j,k) + M1(3,2)*temp4(2)%values(i,j,k) + M1(3,3)*temp4(3)%values(i,j,k)) &
                                   & + MIN(valR(2),0.0_WP)* ( M2(3,1)*temp4(1)%values(i,j,k) + M2(3,2)*temp4(2)%values(i,j,k) + M2(3,3)*temp4(3)%values(i,j,k)) &
                                   & + MIN(valR(3),0.0_WP)* ( M3(3,1)*temp4(1)%values(i,j,k) + M3(3,2)*temp4(2)%values(i,j,k) + M3(3,3)*temp4(3)%values(i,j,k)))

             if (.not. computeSplitingSij(S11%values(i,j,k),S12%values(i,j,k),&
                                         &S13%values(i,j,k),S22%values(i,j,k),&
                                         &S23%values(i,j,k),S33%values(i,j,k),&
                                         &S11out,S12out,S13out,S22out,S23out,S33out,'moins',lam1=lam1,lam2=lam2,lam3=lam3,VPout=TVecP)) then 
               write(6,'(a)') '[ERROR] in RGMSca : new computation of Sij- failed'
               return
             endif
 
             temp5%values(i,j,k) = S11out*temp4(1)%values(i,j,k)+S12out*temp4(2)%values(i,j,k)+S13out*temp4(3)%values(i,j,k)
             temp6%values(i,j,k) = S12out*temp4(1)%values(i,j,k)+S22out*temp4(2)%values(i,j,k)+S23out*temp4(3)%values(i,j,k)
             temp7%values(i,j,k) = S13out*temp4(1)%values(i,j,k)+S23out*temp4(2)%values(i,j,k)+S33out*temp4(3)%values(i,j,k)

             if ( i .eq. (U%xmin+U%xmax)/2 .and. j .eq. (U%ymin+U%ymax)/2 .and. k .eq. (U%zmin+U%zmax)/2 ) then
               write(6,'(a)') ''
               write(6,'(3(a,1x,f10.3,1x))') 'LAPACK lam1=',min(0.0_WP,valR(1)),'lam2=',min(0.0_WP,valR(2)),'lam3=',&
                                             &min(0.0_WP,valR(3))
               write(6,'(3(a,1x,f10.3,1x))') 'ANALYT lam1=',lam1,'lam2=',lam2,'lam3=',lam3
               write(6,'(a)') ''
               write(6,'(a,1x,9(f10.3,1x))') 'LAPACK',VecP
               write(6,'(a,1x,3(f10.3,1x))') 'PLam1 ^ PLam2',VecP(2,1)*VecP(3,2) - VecP(3,1)*VecP(2,2),&
                                                            &VecP(3,1)*VecP(1,2) - VecP(1,1)*VecP(3,2),&
                                                            &VecP(1,1)*VecP(2,2) - VecP(2,1)*VecP(1,2)
               write(6,'(a,1x,9(f10.3,1x))') 'ANALYT',TVecP
               write(6,'(a,1x,3(f10.3,1x))') 'PLam1 ^ PLam2',TVecP(2,1)*TVecP(3,2) - TVecP(3,1)*TVecP(2,2),&
                                                            &TVecP(3,1)*TVecP(1,2) - TVecP(1,1)*TVecP(3,2),&
                                                            &TVecP(1,1)*TVecP(2,2) - TVecP(2,1)*TVecP(1,2)
               write(6,'(a,1x,3(f10.3,1x))') 'LAPACK',temp(1)%values(i,j,k),temp(2)%values(i,j,k),temp(3)%values(i,j,k)
               write(6,'(a,1x,3(f10.3,1x))') 'ANALYT',temp5%values(i,j,k),temp6%values(i,j,k),temp7%values(i,j,k)
             endif

          end do
       end do
    end do
#endif

  if (.not. equalToval (temp(1),temp5) .or. &
     &.not. equalToval (temp(2),temp6) .or. &
     &.not. equalToval (temp(3),temp7) ) then
    write(6,'(a)') '[ERROR] in RGMSca : new computation of Sij- not equal to old one'
    return
  endif

  temp1(1)%values = Cdeltx2 * temp(1)%values
  temp1(2)%values = Cdeltx2 * temp(2)%values
  temp1(3)%values = Cdeltx2 * temp(3)%values

  !Sortie des composantes du vecteur sous maille
  if ( present(Ti)) then
    if ( size(Ti,1) .eq. 3) then
      if ( samelayout(Ti(1),U) ) then
        Ti(1)%values = temp1(1)%values
        Ti(2)%values = temp1(2)%values
        Ti(3)%values = temp1(3)%values
      else
        write(6,'(a,i0,a)')'[ERROR] in computedTidxi : Ti is different from U!'
        return
      endif
    else
      write(6,'(a,i0,a)')'[ERROR] in computedTidxi : Ti is not a vector!'
      return
    endif
  endif
  
  temp(1)%values = temp1(1)%values * temp4(1)%values + temp1(2)%values * temp4(2)%values + temp1(3)%values * temp4(3)%values
  do k = U%zmin, U%zmax
     do j = U%ymin, U%ymax
        do i = U%xmin, U%xmax
           if (temp(1)%values(i,j,k).gt.0.0_WP) then
              print*,'RGM something wrong', temp(1)%values(i,j,k)
           end if
        end do
     end do
  end do
  Cdeltx2 = computeFieldAvg( temp(1) , spec_rank )
!  if (spec_rank.eq.0) print*,'RGM Dissipation =',Cdeltx2
  
  call ftran(temp1(1),tempK(1),res)
  call ftran(temp1(2),tempK(2),res)
  call ftran(temp1(3),tempK(3),res)

  call computeDivergenceK(ScalWN,tempK(1),tempK(2),tempK(3),div)

  call deletedatalayout(S11)
  call deletedatalayout(S12)
  call deletedatalayout(S13)
  call deletedatalayout(S22)
  call deletedatalayout(S23)
  call deletedatalayout(S33)
  
  res = deleteWorkArray(temp)
  res = deleteWorkArray(temp1)
  res = deleteWorkArray(temp2)
  res = deleteWorkArray(temp3)
  res = deleteWorkArray(temp4)
  res = deleteWorkArray(tempK)

  call deletedatalayout(temp5)
  call deletedatalayout(temp6)
  call deletedatalayout(temp7)

 end subroutine RGMScaForTest



!------------------------------------------------------------------------------
!> @author 
!> Guillaume Balarac, LEGI
!>
!> @details
!>  Regularized Gradient model 
!> @param [in]    U  longitudinal velocity in physical space
!> @param [in]    Uk longitudinal velocity in spectral space
!> @param [in]    Vk velocity in spectral space
!> @param [in]    Wk velocity in spectral space
!> @param [in]    Scalk scalar in spectral space 
!> @param [in]    ScalWN wave numbers for scalar 
!> @param [in]    spec_rank is the rank of processus
!> @param [inout] res logical
!> @param [out]   div is the divergence of SGS model
!------------------------------------------------------------------------------
 subroutine RGMSca_lapack(div,U,Uk,Vk,Wk,ScalK,ScalWN,res,spec_rank,coefRGM)

  use stat_tools , only : equalToVal,computeMSE

  !I/O
  type(complex_data_layout),intent(in)    :: Uk,Vk,Wk,ScalK 
  type(complex_data_layout),intent(inout) :: div
  type(real_data_layout),intent(in)       :: U
  type(WaveNumbers), intent(in)           :: ScalWN
  logical, intent(inout)                  :: res
  real(WP), intent(in)                    :: coefRGM
  integer, intent (in)                    :: spec_rank

  !Local data
  type(real_data_layout),dimension(:),allocatable    :: temp,temp1, temp2, temp3, temp4
  type(real_data_layout)                             :: S11, S12, S13, S22, S23, S33
  type(complex_data_layout),dimension(:),allocatable :: tempK
  real(WP)                                           :: Cdeltx2
  
  !! for eigenvalues computation 
  integer :: a, b, info, i,j, k
  real(WP) :: t1,t2
  real(WP), dimension(3,3) :: tab,VecP,VL,TvecP
  real(WP), dimension(3,3) :: E1,E2,E3,M1,M2,M3
  real(WP), dimension(3) :: valR,valI
  real(WP) :: tabwork(102)

  res = .true.
  Cdeltx2 = coefRGM * computeDelta(U)**2.0
  
  res = initWorkArray(U,3,temp)
  res = initWorkArray(U,3,temp1)
  res = initWorkArray(U,3,temp2)
  res = initWorkArray(U,3,temp3)
  res = initWorkArray(U,3,temp4)
  res = initWorkArray(Uk,3,tempK)
 
  call computeGradientK(ScalWN,ScalK,tempK(1),tempK(2),tempK(3))
  call btran(tempK(1),temp4(1),res)
  call btran(tempK(2),temp4(2),res)
  call btran(tempK(3),temp4(3),res)
  call computeGradientK(ScalWN,UK,tempK(1),tempK(2),tempK(3))
  call btran(tempK(1),temp1(1),res)
  call btran(tempK(2),temp1(2),res)
  call btran(tempK(3),temp1(3),res)
  call computeGradientK(ScalWN,VK,tempK(1),tempK(2),tempK(3))
  call btran(tempK(1),temp2(1),res)
  call btran(tempK(2),temp2(2),res)
  call btran(tempK(3),temp2(3),res)
  call computeGradientK(ScalWN,WK,tempK(1),tempK(2),tempK(3))
  call btran(tempK(1),temp3(1),res)
  call btran(tempK(2),temp3(2),res)
  call btran(tempK(3),temp3(3),res)
  
  res = copyStructOnly(U,S11)
  res = copyStructOnly(U,S12)
  res = copyStructOnly(U,S13)
  res = copyStructOnly(U,S22)
  res = copyStructOnly(U,S23)
  res = copyStructOnly(U,S33)
  
  ! compute Sij = 1/2 * ( duidxj + dujdxi)
  S11%values = temp1(1)%values
  S12%values = 0.5 * ( temp1(2)%values + temp2(1)%values )
  S13%values = 0.5 * ( temp1(3)%values + temp3(1)%values )
  S22%values = temp2(2)%values
  S23%values = 0.5 * ( temp2(3)%values + temp3(2)%values )
  S33%values = temp3(3)%values

call cpu_time(t1)  
  ! compute Sij-
    E1 = 0
    E2 = 0
    E3 = 0
    E1(1,1) = 1._WP
    E2(2,2) = 1._WP
    E3(3,3) = 1._WP
#ifdef LAPACK
    do k = U%zmin, U%zmax
       do j = U%ymin, U%ymax
          do i = U%xmin, U%xmax
             tab(1,1) = S11%values(i,j,k)
             tab(1,2) = S12%values(i,j,k)
             tab(1,3) = S13%values(i,j,k)
             tab(2,1) = S12%values(i,j,k)
             tab(2,2) = S22%values(i,j,k)
             tab(2,3) = S23%values(i,j,k)
             tab(3,1) = S13%values(i,j,k)
             tab(3,2) = S23%values(i,j,k)
             tab(3,3) = S33%values(i,j,k)
             call DGEEV('N','V',3,tab,3,valR,valI,VL,3,VecP,3,tabwork,102,INFO)
             do a=1,3
               do b=1,3
                 TvecP(a,b)=VecP(b,a)
               enddo
             enddo
             M1 = matmul(E1,TvecP)
             M1 = matmul(vecP,M1)
             M2 = matmul(E2,TvecP)
             M2 = matmul(vecP,M2)
             M3 = matmul(E3,TvecP)
             M3 = matmul(vecP,M3)
             temp(1)%values(i,j,k) = ( MIN(valR(1),0.0_WP)* ( M1(1,1)*temp4(1)%values(i,j,k) + M1(1,2)*temp4(2)%values(i,j,k) + M1(1,3)*temp4(3)%values(i,j,k)) &
                                   & + MIN(valR(2),0.0_WP)* ( M2(1,1)*temp4(1)%values(i,j,k) + M2(1,2)*temp4(2)%values(i,j,k) + M2(1,3)*temp4(3)%values(i,j,k)) &
                                   & + MIN(valR(3),0.0_WP)* ( M3(1,1)*temp4(1)%values(i,j,k) + M3(1,2)*temp4(2)%values(i,j,k) + M3(1,3)*temp4(3)%values(i,j,k)))

             temp(2)%values(i,j,k) = ( MIN(valR(1),0.0_WP)* ( M1(2,1)*temp4(1)%values(i,j,k) + M1(2,2)*temp4(2)%values(i,j,k) + M1(2,3)*temp4(3)%values(i,j,k)) &
                                   & + MIN(valR(2),0.0_WP)* ( M2(2,1)*temp4(1)%values(i,j,k) + M2(2,2)*temp4(2)%values(i,j,k) + M2(2,3)*temp4(3)%values(i,j,k)) &
                                   & + MIN(valR(3),0.0_WP)* ( M3(2,1)*temp4(1)%values(i,j,k) + M3(2,2)*temp4(2)%values(i,j,k) + M3(2,3)*temp4(3)%values(i,j,k)))

             temp(3)%values(i,j,k) = ( MIN(valR(1),0.0_WP)* ( M1(3,1)*temp4(1)%values(i,j,k) + M1(3,2)*temp4(2)%values(i,j,k) + M1(3,3)*temp4(3)%values(i,j,k)) &
                                   & + MIN(valR(2),0.0_WP)* ( M2(3,1)*temp4(1)%values(i,j,k) + M2(3,2)*temp4(2)%values(i,j,k) + M2(3,3)*temp4(3)%values(i,j,k)) &
                                   & + MIN(valR(3),0.0_WP)* ( M3(3,1)*temp4(1)%values(i,j,k) + M3(3,2)*temp4(2)%values(i,j,k) + M3(3,3)*temp4(3)%values(i,j,k)))
          end do
       end do
    end do
#endif
  call cpu_time(t2)  
  if (spec_rank .eq. 0 ) write (6,'(a,1x,f9.6)') '[INFO LES RGM_lapack] tps de calcul :',t2-t1

  temp1(1)%values = Cdeltx2 * temp(1)%values
  temp1(2)%values = Cdeltx2 * temp(2)%values
  temp1(3)%values = Cdeltx2 * temp(3)%values

  temp(1)%values = temp1(1)%values * temp4(1)%values + temp1(2)%values * temp4(2)%values + temp1(3)%values * temp4(3)%values
  do k = U%zmin, U%zmax
     do j = U%ymin, U%ymax
        do i = U%xmin, U%xmax
           if (temp(1)%values(i,j,k).gt.0.0_WP) then
              print*,'RGM lapack something wrong', temp(1)%values(i,j,k)
           end if
        end do
     end do
  end do
  Cdeltx2 = computeFieldAvg( temp(1) , spec_rank )
  if (spec_rank.eq.0) print*,'RGM lapack Dissipation =',Cdeltx2
  
  call ftran(temp1(1),tempK(1),res)
  call ftran(temp1(2),tempK(2),res)
  call ftran(temp1(3),tempK(3),res)

  call computeDivergenceK(ScalWN,tempK(1),tempK(2),tempK(3),div)

  call deletedatalayout(S11)
  call deletedatalayout(S12)
  call deletedatalayout(S13)
  call deletedatalayout(S22)
  call deletedatalayout(S23)
  call deletedatalayout(S33)
  
  res = deleteWorkArray(temp)
  res = deleteWorkArray(temp1)
  res = deleteWorkArray(temp2)
  res = deleteWorkArray(temp3)
  res = deleteWorkArray(temp4)
  res = deleteWorkArray(tempK)

 end subroutine RGMSca_lapack

!------------------------------------------------------------------------------
!> @author 
!> Guillaume Balarac, LEGI
!>
!> @details
!>  Regularized Gradient model 
!> @param [in]    U  longitudinal velocity in physical space
!> @param [in]    Uk longitudinal velocity in spectral space
!> @param [in]    Vk velocity in spectral space
!> @param [in]    Wk velocity in spectral space
!> @param [in]    Scalk scalar in spectral space 
!> @param [in]    ScalWN wave numbers for scalar 
!> @param [in]    spec_rank is the rank of processus
!> @param [inout] res logical
!> @param [out]   div is the divergence of SGS model
!------------------------------------------------------------------------------
 subroutine RGMSca_analyt(div,U,Uk,Vk,Wk,ScalK,ScalWN,res,spec_rank,coefRGM,Tiout)

  use stat_tools , only : equalToVal,computeMSE

  !I/O
  type(complex_data_layout),intent(in)    :: Uk,Vk,Wk,ScalK 
  type(complex_data_layout),intent(inout) :: div
  type(real_data_layout),dimension(:),intent(inout),optional :: Tiout
  type(real_data_layout),intent(in)       :: U
  type(WaveNumbers), intent(in)           :: ScalWN
  logical, intent(inout)                  :: res
  real(WP), intent(in)                    :: coefRGM
  integer, intent (in)                    :: spec_rank

  !Local data
  type(real_data_layout),dimension(:),allocatable    :: Ti,dZ 
  type(complex_data_layout),dimension(:),allocatable :: tempK
  type(real_data_layout)                             :: S11, S12, S13, S22, S23, S33
  real(WP)                                           :: Cdeltx2,t1,t2
  real(WP)                                           :: S11out,S12out,S13out,S22out,S23out,S33out 
  logical                                            :: success
  integer                                            :: i,j,k
  
  success = .true.
  
  if (.not. (initWorkArray(U,3,Ti) &
     &.and.  initWorkArray(U,3,dZ) &
     &.and.  initWorkArray(Uk,3,tempK) &
     &.and.  copyStructOnly(U,S11) &
     &.and.  copyStructOnly(U,S12) &
     &.and.  copyStructOnly(U,S13) &
     &.and.  copyStructOnly(U,S22) &
     &.and.  copyStructOnly(U,S23) &
     &.and.  copyStructOnly(U,S33))) then
    write(6,'(a)') '[ERROR] in RGM_mod : not enought memory'
    success = .false.
  endif
 
  Cdeltx2 = coefRGM * computeDelta(U)**2
  
  call computeGradientK(ScalWN,ScalK,tempK(1),tempK(2),tempK(3))
  call btran(tempK(1),dZ(1),res)
  if (.not. res ) success = .false. 
  call btran(tempK(2),dZ(2),res)
  if (.not. res ) success = .false. 
  call btran(tempK(3),dZ(3),res)
  if (.not. res ) success = .false. 
  call computeStressTensor(Uk,VK,WK,ScalWN,&
                         &S11,S12,S13,S22,S23,S33,res)
  if (.not. res ) success = .false. 
  call cpu_time(t1)  
  do k = U%zmin, U%zmax
     do j = U%ymin, U%ymax
        do i = U%xmin, U%xmax
           if (.not. computeSplitingSij(S11%values(i,j,k),S12%values(i,j,k),&
                                       &S13%values(i,j,k),S22%values(i,j,k),&
                                       &S23%values(i,j,k),S33%values(i,j,k),&
                                       &S11out,S12out,S13out,S22out,S23out,S33out,'moins')) then 
             write(6,'(a)') '[ERROR] in RGMSca : computeSplitingSij failed'
             success = .false.
           endif
           Ti(1)%values(i,j,k) = S11out*dZ(1)%values(i,j,k)+S12out*dZ(2)%values(i,j,k)+S13out*dZ(3)%values(i,j,k)
           Ti(2)%values(i,j,k) = S12out*dZ(1)%values(i,j,k)+S22out*dZ(2)%values(i,j,k)+S23out*dZ(3)%values(i,j,k)
           Ti(3)%values(i,j,k) = S13out*dZ(1)%values(i,j,k)+S23out*dZ(2)%values(i,j,k)+S33out*dZ(3)%values(i,j,k)
        end do
     end do
  end do
  call cpu_time(t2)  
  if (spec_rank .eq. 0 ) write (6,'(a,1x,f9.6)') '[INFO LES RGM_analyt] tps de calcul :',t2-t1

  Ti(1)%values = Cdeltx2 * Ti(1)%values
  Ti(2)%values = Cdeltx2 * Ti(2)%values
  Ti(3)%values = Cdeltx2 * Ti(3)%values

  do k = U%zmin, U%zmax
     do j = U%ymin, U%ymax
        do i = U%xmin, U%xmax
           S11%values(i,j,k) =  Ti(1)%values(i,j,k) * dZ(1)%values(i,j,k) + &
                              & Ti(2)%values(i,j,k) * dZ(2)%values(i,j,k) + &
                              & Ti(3)%values(i,j,k) * dZ(3)%values(i,j,k)
           if ( S11%values(i,j,k).gt. 0.0_WP ) then
              print*,'RGM analyt something wrong',S11%values(i,j,k)
           end if
        end do
     end do
  end do
  Cdeltx2 = computeFieldAvg( S11 , spec_rank )
  if (spec_rank.eq.0) print*,'RGM analyt Dissipation =',Cdeltx2
  !Sortie des composantes du vecteur sous maille
  if ( present(Tiout)) then
    if ( size(Tiout,1) .eq. 3) then
      if ( samelayout(Tiout(1),U) ) then
        if (spec_rank .eq. 0 ) write (6,'(a)') '[INFO] Tiout...'
        Tiout(1)%values = Ti(1)%values
        Tiout(2)%values = Ti(2)%values
        Tiout(3)%values = Ti(3)%values
      else
        write(6,'(a,i0,a)')'[ERROR] in computedTidxi : Ti is different from U!'
        return
      endif
    else
      write(6,'(a,i0,a)')'[ERROR] in computedTidxi : Ti is not a vector!'
      return
    endif
  endif
  call ftran(Ti(1),tempK(1),res)
  call ftran(Ti(2),tempK(2),res)
  call ftran(Ti(3),tempK(3),res)
  call computeDivergenceK(ScalWN,tempK(1),tempK(2),tempK(3),div)

  call deletedatalayout(S11)
  call deletedatalayout(S12)
  call deletedatalayout(S13)
  call deletedatalayout(S22)
  call deletedatalayout(S23)
  call deletedatalayout(S33)
  res = deleteWorkArray(Ti)
  res = deleteWorkArray(dZ)
  res = deleteWorkArray(tempK)

 end subroutine RGMSca_analyt


!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!>  Gradient model 
!> @param [inout] div is the divergence of SGS model
!> @param [in]    U  longitudinal velocity in physical space (LES quantity)
!> @param [in]    V  vertical velocity in physical space (LES quantity)
!> @param [in]    W  spanwize velocity in physical space (LES quantity)
!> @param [in]    Uk longitudinal velocity in spectral space (LES quantity)
!> @param [in]    Vk velocity in spectral space (LES quantity)
!> @param [in]    Wk velocity in spectral space (LES quantity)
!> @param [in]    Scalk scalar in spectral space (LES quantity)
!> @param [in]    ScalWN wave numbers for scalar (LES quantity)
!> @param [in]    filtersize ( = nd) is the size of filter used in a priori test 
!> @param [in]    newCoef is an optional coef which can be set at the model 
!> @param [inout] res logical
!------------------------------------------------------------------------------
 subroutine ScaleSimilaritySca(div,U,V,W,Scal,Uk,Vk,Wk,ScalK,ScalWN,res,filter,filterSize,newCoef)

  !I/O
  type(complex_data_layout),intent(in)               :: Uk,Vk,Wk,ScalK 
  type(complex_data_layout),intent(inout)            :: div
  type(real_data_layout),intent(in)                  :: U,V,W,Scal
  type(WaveNumbers), intent(in)                      :: ScalWN
  logical, intent(inout)                             :: res
  integer, intent(in)                                :: filterSize
  integer, intent(in)                                :: filter
  real(WP), intent(in),optional                      :: newCoef

  !Local data
  type(real_data_layout),dimension(:),allocatable    :: temp
  type(complex_data_layout),dimension(:),allocatable :: tempK
  real(WP)                                           :: Coef,delta
  real(WP)                                           :: sizeTestFilter

  res = .true.
  delta = computeDelta(Scal)
  sizeTestFilter = 2.0_WP
  if ( present(newCoef) ) then
    Coef = newCoef * (real(filterSize*computeDelta(Scal),WP))**2
  else
    Coef = 1.0_WP/5.0_WP * (real(filterSize*computeDelta(Scal),WP))**2
  endif

 
  if (.not. initWorkArray(U,7,temp)  .or. &
     &.not. initWorkArray(Uk,3,tempK) ) then
    res = .false.
    write(6,'(a)')'[ERROR] in gradientSca : not enought memory'
  endif
 
  !dTidxi = 0.2 * delta**2 * ( (uiBar ZBar)hat - uiBarHat ZbatHat) 
  temp(1)%values = U%values * Scal%values
  temp(2)%values = V%values * Scal%values
  temp(3)%values = W%values * Scal%values
  
  call computeFilter(ScalWN,sizeTestFilter*filtersize*delta,temp(1),temp(5),filter) !UZ bar
  call computeFilter(ScalWN,sizeTestFilter*filtersize*delta,temp(2),temp(6),filter) !UZ bar
  call computeFilter(ScalWN,sizeTestFilter*filtersize*delta,temp(3),temp(7),filter) !UZ bar
  call computeFilter(ScalWN,sizeTestFilter*filtersize*delta,Uk,temp(1),filter) !U bar
  call computeFilter(ScalWN,sizeTestFilter*filtersize*delta,Vk,temp(2),filter) !V bar
  call computeFilter(ScalWN,sizeTestFilter*filtersize*delta,Wk,temp(3),filter) !W bar
  call computeFilter(ScalWN,sizeTestFilter*filtersize*delta,Scalk,temp(4),filter) !Z bar
  
  temp(1)%values = Coef * (temp(5)%values - temp(1)%values * temp(4)%values)
  temp(2)%values = Coef * (temp(6)%values - temp(2)%values * temp(4)%values)
  temp(3)%values = Coef * (temp(7)%values - temp(3)%values * temp(4)%values) 

  call ftran(temp(1),tempK(1),res)
  if (.not. res) write(6,'(a)')'[ERROR] in ScaleSimilaritySca for scalar : compute ftran 1 failed'
  call ftran(temp(2),tempK(2),res)
  if (.not. res) write(6,'(a)')'[ERROR] in ScaleSimilaritySca for scalar : compute ftran 2 failed'
  call ftran(temp(3),tempK(3),res)
  if (.not. res) write(6,'(a)')'[ERROR] in ScaleSimilaritySca for scalar : compute ftran 3 failed'
    
  call computeDivergenceK(ScalWN,tempK(1),tempK(2),tempK(3),div)

  res = deleteWorkArray(temp)
  res = deleteWorkArray(tempK)

 end subroutine ScaleSimilaritySca 



!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!>  Gradient model 
!> @param [inout] div is the divergence of SGS model
!> @param [in]    U  longitudinal velocity in physical space
!> @param [in]    Uk longitudinal velocity in spectral space
!> @param [in]    Vk velocity in spectral space
!> @param [in]    Wk velocity in spectral space
!> @param [in]    Scalk scalar in spectral space 
!> @param [in]    ScalWN wave numbers for scalar 
!> @param [in]    filtersize ( = nd) is the size of filter used in a priori test 
!> @param [in]    newCoef is an optional coef which can be set at the model 
!> @param [inout] res logical
!------------------------------------------------------------------------------
 subroutine GradientSca(div,U,Uk,Vk,Wk,ScalK,Scal,ScalWN,res,filterSize,newCoef)

  !I/O
  type(complex_data_layout),intent(in)               :: Uk,Vk,Wk,ScalK 
  type(complex_data_layout),intent(inout)            :: div
  type(real_data_layout),intent(in)                  :: U,Scal
  type(WaveNumbers), intent(in)                      :: ScalWN
  logical, intent(inout)                             :: res
  integer, intent(in)                                :: filterSize
  real(WP), intent(in),optional                      :: newCoef

  !Local data
  type(real_data_layout),dimension(:),allocatable    :: temp1, temp2, temp3, temp4 
  type(complex_data_layout),dimension(:),allocatable :: tempK
  real(WP)                                           :: Coef
  logical                                            :: res1,res2,res3,res4,res5
  integer   :: spec_rank

  res = .true.
  !spec_rank = getnbmycpu()
  if ( present(newCoef) ) then
    !if (spec_rank .eq. 0) write (6,'(a,1x,f8.5)') '[INFO] Coef for gradient model =',newCoef
    Coef = newCoef * (real(filterSize*computeDelta(Scal),WP))**2.0
  else
    !if (spec_rank .eq. 0) write (6,'(a)') '[INFO] Coef for gradient model = 1/12'
    Coef = 1.0_WP/12.0_WP * (real(filterSize*computeDelta(Scal),WP))**2.0
  endif
 
 
  if (.not. initWorkArray(U,3,temp1)  .or. &
     &.not. initWorkArray(U,3,temp2)  .or. &
     &.not. initWorkArray(U,3,temp3)  .or. &
     &.not. initWorkArray(U,3,temp4)  .or. &
     &.not. initWorkArray(Uk,3,tempK) ) then
    res = .false.
    write(6,'(a)')'[ERROR] in gradientSca : not enought memory'
  endif
  
  call computeGradientK(ScalWN,ScalK,tempK(1),tempK(2),tempK(3))
  call btran(tempK(1),temp4(1),res1)
  call btran(tempK(2),temp4(2),res2)
  call btran(tempK(3),temp4(3),res3)
  if (.not. res1 .or. .not.res2 .or. .not.res3 ) then
    write(6,'(a)')'[ERROR] in gradientSca for scalar : compute btran 1 failed'
    res = .false.
  endif
  call computeGradientK(ScalWN,UK,tempK(1),tempK(2),tempK(3))
  call btran(tempK(1),temp1(1),res1)
  call btran(tempK(2),temp1(2),res2)
  call btran(tempK(3),temp1(3),res3)
  if (.not. res1 .or. .not.res2 .or. .not.res3 ) then
    write(6,'(a)')'[ERROR] in gradientSca for scalar : compute btran 2 failed'
    res = .false.
  endif
  call computeGradientK(ScalWN,VK,tempK(1),tempK(2),tempK(3))
  call btran(tempK(1),temp2(1),res1)
  call btran(tempK(2),temp2(2),res2)
  call btran(tempK(3),temp2(3),res3)
  if (.not. res1 .or. .not.res2 .or. .not.res3 ) then
    write(6,'(a)')'[ERROR] in gradientSca for scalar : compute btran 3 failed'
    res = .false.
  endif
  call computeGradientK(ScalWN,WK,tempK(1),tempK(2),tempK(3))
  call btran(tempK(1),temp3(1),res1)
  call btran(tempK(2),temp3(2),res2)
  call btran(tempK(3),temp3(3),res3)
  if (.not. res1 .or. .not.res2 .or. .not.res3 ) then
    write(6,'(a)')'[ERROR] in gradientSca for scalar : compute btran 4 failed'
    res = .false.
  endif
  temp1(1)%values = Coef * (temp1(1)%values*temp4(1)%values + temp1(2)%values*temp4(2)%values +&
                    & temp1(3)%values*temp4(3)%values)
  temp1(2)%values = Coef * (temp2(1)%values*temp4(1)%values + temp2(2)%values*temp4(2)%values +&
                    & temp2(3)%values*temp4(3)%values)
  temp1(3)%values = Coef * (temp3(1)%values*temp4(1)%values + temp3(2)%values*temp4(2)%values +&
                    & temp3(3)%values*temp4(3)%values)
  call ftran(temp1(1),tempK(1),res1)
  call ftran(temp1(2),tempK(2),res2)
  call ftran(temp1(3),tempK(3),res3)
  if (.not. res1 .or. .not. res2 .or. .not.res3 ) then
    write(6,'(a)')'[ERROR] in gradientSca for scalar : compute ftran 1 failed'
    res = .false.
  endif
  call computeDivergenceK(ScalWN,tempK(1),tempK(2),tempK(3),div)

  res1 = deleteWorkArray(temp1)
  res2 = deleteWorkArray(temp2)
  res3 = deleteWorkArray(temp3)
  res4 = deleteWorkArray(temp4)
  res5 = deleteWorkArray(tempK)
#ifdef DEBUGLES
  if (.not. res1 .or. .not. res2 .or. .not. res3 .or. .not. res4 .or. .not. res5) then
    write(6,'(a)')'[ERROR] in gradientSca for scalar : compute ftran 1 failed'
    res = .false.
  endif
#endif

 end subroutine GradientSca 


!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!>  Smagorinsky model 
!> @param [in]    U  longitudinal velocity in physical space
!> @param [in]    Uk longitudinal velocity in spectral space
!> @param [in]    Vk vertical velocity in spectral space
!> @param [in]    Wk spanwize velocity in spectral space
!> @param [in]    Scalk scalar in spectral space 
!> @param [in]    ScalWN wave numbers for scalar 
!> @param [in]    spec_rank is the rank of processus
!> @param [inout] res logical
!> @param [in]    coeff is the dynamic coef (optinal) 
!> @param [out]   div is the divergence of SGS model
!------------------------------------------------------------------------------
 subroutine SmagorinskySca(div,U,Uk,Vk,Wk,ScalK,ScalWN,res,spec_rank)

  !I/O
  type(complex_data_layout),intent(in)          :: Uk,Vk,Wk,ScalK 
  type(complex_data_layout),intent(inout)       :: div
  type(real_data_layout),intent(in)             :: U
  type(WaveNumbers), intent(in)                 :: ScalWN
  logical, intent(inout)                        :: res
  integer, intent(in)                           :: spec_rank

  !Local data
  type(real_data_layout),dimension(:),allocatable    :: temp 
  type(complex_data_layout),dimension(:),allocatable :: tempK
  real(WP)                                           :: deltx
  real(WP)                                           :: coeff

  res = .true.
  if (.not. initWorkArray(U,3,temp)  .or. &
     &.not. initWorkArray(Uk,3,tempK) ) then
    res = .false.
    write(6,'(a)')'[ERROR] in smagorinskySca : not enought memory'
  endif
  deltx = computeDelta(U)
  coeff = -0.18_WP**2/0.6_WP
  call computeSdZ(U,Uk,Vk,Wk,ScalK,ScalWN,temp(1),temp(2),temp(3),res)
  temp(1)%values = coeff * deltx**2 * temp(1)%values
  temp(2)%values = coeff * deltx**2 * temp(2)%values
  temp(3)%values = coeff * deltx**2 * temp(3)%values

  call ftran(temp(1),tempK(1),res)
  if (.not. res ) write(6,'(a)')'[ERROR] in smagorinsky for scalar : compute ftran 1 failed'
  call ftran(temp(2),tempK(2),res)
  if (.not. res ) write(6,'(a)')'[ERROR] in smagorinsky for scalar : compute ftran 2 failed'
  call ftran(temp(3),tempK(3),res)
  if (.not. res ) write(6,'(a)')'[ERROR] in smagorinsky for scalar : compute ftran 3 failed'
  call computeDivergenceK(ScalWN,tempK(1),tempK(2),tempK(3),div)

  if (.not. deleteWorkArray(temp)  .or. &
     &.not. deleteWorkArray(tempK) ) then
    res = .false.
    write(6,'(a)')'[ERROR] in smagorinskySca : memory leak'
  endif

 end subroutine SmagorinskySca 



!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!>  Smagorinsky model 
!> @param [in]    U  longitudinal velocity in physical space
!> @param [in]    Uk longitudinal velocity in spectral space
!> @param [in]    Vk vertical velocity in spectral space
!> @param [in]    Wk spanwize velocity in spectral space
!> @param [in]    Scalk scalar in spectral space 
!> @param [in]    ScalWN wave numbers for scalar 
!> @param [in]    spec_rank is the rank of processus
!> @param [inout] res logical
!> @param [in]    coeff is the dynamic coef (optinal) 
!> @param [out]   div is the divergence of SGS model
!------------------------------------------------------------------------------
 subroutine DynSmagorinskySca(div,U,V,W,Scal,Uk,Vk,Wk,ScalK,ScalWN,res,filter,spec_rank,type_avg)

  !I/O
  type(complex_data_layout),intent(in)          :: Uk,Vk,Wk,ScalK 
  type(complex_data_layout),intent(inout)       :: div
  type(real_data_layout),intent(in)             :: U,V,W,Scal
  type(WaveNumbers), intent(in)                 :: ScalWN
  logical, intent(inout)                        :: res
  integer, intent(in)                           :: spec_rank
  integer, intent(in)                           :: filter
  integer, intent(in)                           :: type_avg

  !Local data
  type(real_data_layout),dimension(:),allocatable    :: Ti 
  type(complex_data_layout),dimension(:),allocatable :: TiK
  real(WP)                                           :: deltx
  real(WP)                                           :: coeff
  real(WP),dimension(:),allocatable                  :: coefArray
  integer                                            :: ires,k

  res = .true.
  if (.not. initWorkArray(U,3,Ti)  .or. &
     &.not. initWorkArray(Uk,3,TiK) ) then
    res = .false.
    write(6,'(a)')'[ERROR] in smagorinskySca : not enought memory'
  endif
  deltx = computeDelta(U)
  call computeSdZ(U,Uk,Vk,Wk,ScalK,ScalWN,Ti(1),Ti(2),Ti(3),res)

  select case (type_avg)
    case(0) !Cas isotrope
      if (.not. coefDynSmagSca(U,V,W,Uk,Vk,Wk,Scal,ScalK,ScalWN,filter,1,spec_rank,&
         & CoeffDynScalar=coeff) ) then
        write (6,'(a)') '[ERROR] in DynSmagorinskySca_3d'
        res = .false.
      endif
      Ti(1)%values = coeff * deltx**2 * Ti(1)%values
      Ti(2)%values = coeff * deltx**2 * Ti(2)%values
      Ti(3)%values = coeff * deltx**2 * Ti(3)%values
    case(1) !Cas anisotrope moyenne dans les plans (X,Y)
      allocate(coefArray(Scal%nz),stat=ires)
      if ( ires .ne. 0 ) then
        write (6,'(a)') '[ERROR] not enought memory to allocate coeff'
        res = .false.
      endif 
      if (.not. coefDynSmagSca(U,V,W,Uk,Vk,Wk,Scal,ScalK,ScalWN,filter,1,spec_rank,&
         & CoeffDynArray=coefArray,dir2DAvg=3) ) then
        write (6,'(a)') '[ERROR] in DynSmagorinskySca_3d'
        res = .false.
      endif
      do k = U%zmin,U%zmax
         Ti(1)%values(:,:,k) = coefArray(k) * deltx**2.0 * Ti(1)%values(:,:,k)
         Ti(2)%values(:,:,k) = coefArray(k) * deltx**2.0 * Ti(2)%values(:,:,k)
         Ti(3)%values(:,:,k) = coefArray(k) * deltx**2.0 * Ti(3)%values(:,:,k)
      end do
      deallocate(coefArray)
    case(2) !Cas anisotrope moyenne dans les plans (X,Z)
      allocate(coefArray(Scal%ny),stat=ires)
      if ( ires .ne. 0 ) then
        write (6,'(a)') '[ERROR] not enought memory to allocate coeff'
        res = .false.
      endif 
      if (.not. coefDynSmagSca(U,V,W,Uk,Vk,Wk,Scal,ScalK,ScalWN,filter,1,spec_rank,&
         & CoeffDynArray=coefArray,dir2DAvg=2) ) then
        write (6,'(a)') '[ERROR] in DynSmagorinskySca_3d'
        res = .false.
      endif
      do k = U%ymin,U%ymax
         Ti(1)%values(:,k,:) = coefArray(k) * deltx**2.0 * Ti(1)%values(:,k,:)
         Ti(2)%values(:,k,:) = coefArray(k) * deltx**2.0 * Ti(2)%values(:,k,:)
         Ti(3)%values(:,k,:) = coefArray(k) * deltx**2.0 * Ti(3)%values(:,k,:)
      end do
      deallocate(coefArray)
    case(3) !Cas anisotrope moyenne dans les plans (Y,Z)
      allocate(coefArray(Scal%nx),stat=ires)
      if ( ires .ne. 0 ) then
        write (6,'(a)') '[ERROR] not enought memory to allocate coeff'
        res = .false.
      endif 
      if (.not. coefDynSmagSca(U,V,W,Uk,Vk,Wk,Scal,ScalK,ScalWN,filter,1,spec_rank,&
         & CoeffDynArray=coefArray,dir2DAvg=1) ) then
        write (6,'(a)') '[ERROR] in DynSmagorinskySca_3d'
        res = .false.
      endif
      do k = U%xmin,U%xmax
         Ti(1)%values(k,:,:) = coefArray(k) * deltx**2.0 * Ti(1)%values(k,:,:)
         Ti(2)%values(k,:,:) = coefArray(k) * deltx**2.0 * Ti(2)%values(k,:,:)
         Ti(3)%values(k,:,:) = coefArray(k) * deltx**2.0 * Ti(3)%values(k,:,:)
      end do
      deallocate(coefArray)
    case default
      write (6,'(a)') '[ERROR] type of averaging unknow in Dynamic Smagorinsky procedure'
      res = .false.
  end select

  call ftran(Ti(1),TiK(1),res)
  if (.not. res ) write(6,'(a)')'[ERROR] in smagorinsky for scalar : compute ftran 1 failed'
  call ftran(Ti(2),TiK(2),res)
  if (.not. res ) write(6,'(a)')'[ERROR] in smagorinsky for scalar : compute ftran 2 failed'
  call ftran(Ti(3),TiK(3),res)
  if (.not. res ) write(6,'(a)')'[ERROR] in smagorinsky for scalar : compute ftran 3 failed'
  call computeDivergenceK(ScalWN,TiK(1),TiK(2),TiK(3),div)

  if (.not. deleteWorkArray(Ti)  .or. &
     &.not. deleteWorkArray(TiK) ) then
    res = .false.
    write(6,'(a)')'[ERROR] in smagorinskySca : memory leak'
  endif

 end subroutine DynSmagorinskySca



!!!!!!!!!!!!!!!MODELE de TEST!!!!!!!!!!!!!!!!!!!

!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> Sub-grid model based on NDCM one variable
!> @param [in]    U  longitudinal velocity in physical space
!> @param [in]    V  velocity in physical space
!> @param [in]    W  velocity in physical space
!> @param [in]    Uk longitudinal velocity in spectral space
!> @param [in]    Vk velocity in spectral space
!> @param [in]    Wk velocity in spectral space
!> @param [in]    ScalIn scalar in physical space 
!> @param [in]    ScalkIn scalar in spectral space 
!> @param [in]    WaveNum wave numbers for scalar 
!> @param [in]    Filter is the type of test filter 
!> @param [in]    spec_rank is the rank of processus
!> @param [inout] res logical
!> @param [inout] div is the divergence of SGS model
!> @return
!------------------------------------------------------------------------------
 subroutine DynClarkFabreSca_2d_z(divK,U,V,W,Uk,Vk,Wk,ScalIn,ScalKIn,WaveNum,filter,spec_rank,res)
  !I/O
  type(complex_data_layout),intent(inout) :: divK
  type(complex_data_layout),intent(in)    :: Uk,Vk,Wk,ScalKIn
  type(real_data_layout),intent(in)       :: U, V, W, ScalIn
  type(WaveNumbers), intent(in)           :: WaveNum
  logical, intent(inout)                  :: res
  integer, intent(in)                     :: spec_rank
  integer, intent(in)                     :: filter

  !Local data
  type(real_data_layout),dimension(:),allocatable    :: divSmag
  type(real_data_layout),dimension(:),allocatable    :: divGrad
  type(real_data_layout),dimension(:),allocatable    :: div
  type(real_data_layout),dimension(:),allocatable    :: SdZ
  type(complex_data_layout),dimension(:),allocatable :: divGradK
  type(complex_data_layout),dimension(:),allocatable :: divSmagK
  type(complex_data_layout),dimension(:),allocatable :: SdZK
  real(WP)                                           :: delta
  real(WP),dimension(:),allocatable                  :: coefFabre
  integer                                            :: ires,i
 
  res = .true.

  !Init des tableaux de travail
  if (.not. initWorkArray(U,1,divSmag) .or. &
     &.not. initWorkArray(U,1,divGrad) .or. &
     &.not. initWorkArray(U,1,div) .or. &
     &.not. initWorkArray(U,3,SdZ) .or. &
     &.not. initWorkArray(Uk,1,divSmagK) .or. &
     &.not. initWorkArray(Uk,1,divGradK) .or. &
     &.not. initWorkArray(Uk,3,SdZK)) then
    write(6,'(a)')'[ERROR] initWorkArray : not enought memory!'
    return
  endif
  !allocate(coefFabre(ScalIn%nz),stat=ires)
  allocate(coefFabre(ScalIn%nz),stat=ires)
  if ( ires .ne. 0 ) then
    write (6,'(a)') '[ERROR] not enought memory to allocate coefFabre'
    res = .false.
  endif 
  
  delta = computeDelta(ScalIn)

  !CALCUL des TERMES
  !Coefficient dynamique
  if (.not. coefDynClarkFabreSca(U,V,W,Uk,Vk,Wk,ScalIn,ScalKIn,WaveNum,filter,1,&
     &  spec_rank,CoeffDynArray=coefFabre) ) then
     write (6,'(a)') '[ERROR] Coef dynamic failed !'
  endif
  if (.not. coefDynClarkFabreScaTEST(U,V,W,Uk,Vk,Wk,ScalIn,ScalKIn,WaveNum,filter,1,&
     &  spec_rank,CoeffDynArray=coefFabre) ) then
     write (6,'(a)') '[ERROR] Coef dynamic failed !'
  endif
  !Gradient
  call gradientSca(divGradK(1),U,Uk,Vk,Wk,ScalKIn,ScalIn,WaveNum,res,1)
  call btran(divGradK(1),divGrad(1),res)
  !Smagorinsky
  call computeSdZ(U,Uk,Vk,Wk,ScalKIn,WaveNum,SdZ(1),SdZ(2),SdZ(3),res)
  do i=ScalIn%zmin,ScalIn%zmax
    SdZ(1)%values(:,:,i) = coefFabre(i)*delta**2*SdZ(1)%values(:,:,i)
    SdZ(2)%values(:,:,i) = coefFabre(i)*delta**2*SdZ(2)%values(:,:,i)
    SdZ(3)%values(:,:,i) = coefFabre(i)*delta**2*SdZ(3)%values(:,:,i)
  enddo
  call ftran(SdZ(1),SdZK(1),res)
  if (.not. res ) write(6,'(a)')'[ERROR] in smagorinsky for scalar : compute ftran 1 failed'
  call ftran(SdZ(2),SdZK(2),res)
  if (.not. res ) write(6,'(a)')'[ERROR] in smagorinsky for scalar : compute ftran 2 failed'
  call ftran(SdZ(3),SdZK(3),res)
  if (.not. res ) write(6,'(a)')'[ERROR] in smagorinsky for scalar : compute ftran 3 failed'
  call computeDivergenceK(WaveNum,SdZK(1),SdZK(2),SdZK(3),divSmagK(1))
  call btran(divSmagK(1),divSmag(1),res)
  !Mixage
  div(1)%values = divSmag(1)%values + divGrad(1)%values
  call ftran(div(1),divK,res)

  !Désallocation des tableaux de travail
  deallocate(coefFabre)
  if ( .not. deleteWorkArray(divGrad) .or. &
  & .not. deleteWorkArray(divGradK) .or. &
  & .not. deleteWorkArray(divSmag) .or. &
  & .not. deleteWorkArray(divSmagK) .or. &
  & .not. deleteWorkArray(div) .or. &
  & .not. deleteWorkArray(SdZ) .or. &
  & .not. deleteWorkArray(SdZK)) then
    res = .false.
    write(6,'(a)')'[ERROR] in dynClarkFabresca_2d_z: memory leak'
  endif

 end subroutine DynClarkFabreSca_2d_z

!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> Sub-grid model based on NDCM one variable
!> @param [in]    U  longitudinal velocity in physical space
!> @param [in]    V  velocity in physical space
!> @param [in]    W  velocity in physical space
!> @param [in]    Uk longitudinal velocity in spectral space
!> @param [in]    Vk velocity in spectral space
!> @param [in]    Wk velocity in spectral space
!> @param [in]    ScalIn scalar in physical space 
!> @param [in]    ScalkIn scalar in spectral space 
!> @param [in]    WaveNum wave numbers for scalar 
!> @param [in]    Filter is the type of test filter 
!> @param [in]    spec_rank is the rank of processus
!> @param [inout] res logical
!> @param [inout] div is the divergence of SGS model
!> @return
!------------------------------------------------------------------------------
 subroutine DynClarkFabreSca_3d(divK,U,V,W,Uk,Vk,Wk,ScalIn,ScalKIn,WaveNum,filter,spec_rank,res)
  !I/O
  type(complex_data_layout),intent(inout) :: divK
  type(complex_data_layout),intent(in)    :: Uk,Vk,Wk,ScalKIn
  type(real_data_layout),intent(in)       :: U, V, W, ScalIn
  type(WaveNumbers), intent(in)           :: WaveNum
  logical, intent(inout)                  :: res
  integer, intent(in)                     :: spec_rank
  integer, intent(in)                     :: filter

  !Local data
  type(real_data_layout),dimension(:),allocatable    :: divSmag
  type(real_data_layout),dimension(:),allocatable    :: divGrad
  type(real_data_layout),dimension(:),allocatable    :: div
  type(real_data_layout),dimension(:),allocatable    :: SdZ
  type(complex_data_layout),dimension(:),allocatable :: divGradK
  type(complex_data_layout),dimension(:),allocatable :: divSmagK
  type(complex_data_layout),dimension(:),allocatable :: SdZK
  real(WP)                                           :: delta,coefFabre
 
  res = .true.

  !Init des tableaux de travail
  if (.not. initWorkArray(U,1,divSmag) .or. &
     &.not. initWorkArray(U,1,divGrad) .or. &
     &.not. initWorkArray(U,1,div) .or. &
     &.not. initWorkArray(U,3,SdZ) .or. &
     &.not. initWorkArray(Uk,1,divSmagK) .or. &
     &.not. initWorkArray(Uk,1,divGradK) .or. &
     &.not. initWorkArray(Uk,3,SdZK)) then
    write(6,'(a)')'[ERROR] initWorkArray : not enought memory!'
    return
  endif
  delta = computeDelta(ScalIn)

  !CALCUL des TERMES
  !Coefficient dynamique
  if (.not. coefDynClarkFabreSca(U,V,W,Uk,Vk,Wk,ScalIn,ScalKIn,WaveNum,filter,1,&
     &  spec_rank,CoeffDynScalar=coefFabre) ) then
     write (6,'(a)') '[ERROR] Coef dynamic failed !'
  endif
  if (.not. coefDynClarkFabreScaTEST(U,V,W,Uk,Vk,Wk,ScalIn,ScalKIn,WaveNum,filter,1,&
     &  spec_rank,CoeffDynScalar=coefFabre) ) then
     write (6,'(a)') '[ERROR] Coef dynamic failed !'
  endif
  !Gradient
  call gradientSca(divGradK(1),U,Uk,Vk,Wk,ScalKIn,ScalIn,WaveNum,res,1)
  call btran(divGradK(1),divGrad(1),res)
  !Smagorinsky
  call computeSdZ(U,Uk,Vk,Wk,ScalKIn,WaveNum,SdZ(1),SdZ(2),SdZ(3),res)
  SdZ(1)%values = coefFabre*delta**2*SdZ(1)%values
  SdZ(2)%values = coefFabre*delta**2*SdZ(2)%values
  SdZ(3)%values = coefFabre*delta**2*SdZ(3)%values
  call ftran(SdZ(1),SdZK(1),res)
  if (.not. res ) write(6,'(a)')'[ERROR] in smagorinsky for scalar : compute ftran 1 failed'
  call ftran(SdZ(2),SdZK(2),res)
  if (.not. res ) write(6,'(a)')'[ERROR] in smagorinsky for scalar : compute ftran 2 failed'
  call ftran(SdZ(3),SdZK(3),res)
  if (.not. res ) write(6,'(a)')'[ERROR] in smagorinsky for scalar : compute ftran 3 failed'
  call computeDivergenceK(WaveNum,SdZK(1),SdZK(2),SdZK(3),divSmagK(1))
  call btran(divSmagK(1),divSmag(1),res)
  !Mixage
  div(1)%values = divSmag(1)%values + divGrad(1)%values
  call ftran(div(1),divK,res)

  !Désallocation des tableaux de travail
  if ( .not. deleteWorkArray(divGrad) .or. &
  & .not. deleteWorkArray(divGradK) .or. &
  & .not. deleteWorkArray(divSmag) .or. &
  & .not. deleteWorkArray(divSmagK) .or. &
  & .not. deleteWorkArray(div) .or. &
  & .not. deleteWorkArray(SdZ) .or. &
  & .not. deleteWorkArray(SdZK)) then
    res = .false.
    write(6,'(a)')'[ERROR] in dynClarkFabresca_3d: memory leak'
  endif

 end subroutine DynClarkFabreSca_3d




!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> Sub-grid model based on DCM one variable
!> @param [in]    U  longitudinal velocity in physical space
!> @param [in]    V  velocity in physical space
!> @param [in]    W  velocity in physical space
!> @param [in]    Uk longitudinal velocity in spectral space
!> @param [in]    Vk velocity in spectral space
!> @param [in]    Wk velocity in spectral space
!> @param [in]    tabEOFile contains the EO model
!> @param [in]    Scal scalar in physical space 
!> @param [in]    Scalk scalar in spectral space 
!> @param [in]    ScalWN wave numbers for scalar 
!> @param [in]    spec_rank is the rank of processus
!> @param [inout] res logical
!> @param [inout] div is the divergence of SGS model
!> @return
!------------------------------------------------------------------------------
 subroutine DynClarkSca_2d_z(divK,U,V,W,Uk,Vk,Wk,ScalIn,ScalKIn,WaveNum,filter,spec_rank,res)
  !I/O
  type(complex_data_layout),intent(inout) :: divK
  type(complex_data_layout),intent(in)    :: Uk,Vk,Wk,ScalKIn
  type(real_data_layout),intent(in)       :: U, V, W, ScalIn
  type(WaveNumbers), intent(in)           :: WaveNum
  logical, intent(inout)                  :: res
  integer, intent(in)                     :: spec_rank
  integer, intent(in)                     :: filter

  !Local data
  type(real_data_layout),dimension(:),allocatable    :: divSmag
  type(real_data_layout),dimension(:),allocatable    :: divGrad
  type(real_data_layout),dimension(:),allocatable    :: div
  type(real_data_layout),dimension(:),allocatable    :: SdZ
  type(complex_data_layout),dimension(:),allocatable :: divGradK
  type(complex_data_layout),dimension(:),allocatable :: divSmagK
  type(complex_data_layout),dimension(:),allocatable :: SdZK
  real(WP)                                           :: delta
  real(WP),dimension(:),allocatable                  :: coefClark
  integer                                            :: ires,i
 
  res = .true.

  !Init des tableaux de travail
  if (.not. initWorkArray(U,1,divSmag) .or. &
     &.not. initWorkArray(U,1,divGrad) .or. &
     &.not. initWorkArray(U,1,div) .or. &
     &.not. initWorkArray(U,3,SdZ) .or. &
     &.not. initWorkArray(Uk,1,divSmagK) .or. &
     &.not. initWorkArray(Uk,1,divGradK) .or. &
     &.not. initWorkArray(Uk,3,SdZK)) then
    write(6,'(a)')'[ERROR] initWorkArray : not enought memory!'
    return
  endif
  allocate(coefClark(ScalIn%nz),stat=ires)
  if ( ires .ne. 0 ) then
    write (6,'(a)') '[ERROR] not enought memory to allocate coefClark'
    res = .false.
  endif 
  delta = computeDelta(ScalIn)
  !Calcul du coef dynamique de Clark 
  if (.not. coefDynClarkSca(U,V,W,Uk,Vk,Wk,ScalIn,ScalKIn,WaveNum,filter,1,&
     &      spec_rank,CoeffDynArray=coefClark) ) then
    write (6,'(a)') '[INFO] coef Dynamic failed' 
  endif
  if (.not. coefDynClarkScaTEST(U,V,W,Uk,Vk,Wk,ScalIn,ScalKIn,WaveNum,filter,1,&
     &      spec_rank,CoeffDynArray=coefClark) ) then
    write (6,'(a)') '[INFO] coef Dynamic failed' 
  endif
  !CALCUL des TERMES
  !Gradient
  call gradientSca(divGradK(1),U,Uk,Vk,Wk,ScalKIn,ScalIn,WaveNum,res,1)
  call btran(divGradK(1),divGrad(1),res)
  !Smagorinsky
  call computeSdZ(U,Uk,Vk,Wk,ScalKIn,WaveNum,SdZ(1),SdZ(2),SdZ(3),res)
  do i=ScalIn%zmin,ScalIn%zmax
    SdZ(1)%values(:,:,i) = coefClark(i)*delta**2*SdZ(1)%values(:,:,i)
    SdZ(2)%values(:,:,i) = coefClark(i)*delta**2*SdZ(2)%values(:,:,i)
    SdZ(3)%values(:,:,i) = coefClark(i)*delta**2*SdZ(3)%values(:,:,i)
  enddo
  call ftran(SdZ(1),SdZK(1),res)
  if (.not. res ) write(6,'(a)')'[ERROR] in smagorinsky for scalar : compute ftran 1 failed'
  call ftran(SdZ(2),SdZK(2),res)
  if (.not. res ) write(6,'(a)')'[ERROR] in smagorinsky for scalar : compute ftran 2 failed'
  call ftran(SdZ(3),SdZK(3),res)
  if (.not. res ) write(6,'(a)')'[ERROR] in smagorinsky for scalar : compute ftran 3 failed'
  call computeDivergenceK(WaveNum,SdZK(1),SdZK(2),SdZK(3),divSmagK(1))
  call btran(divSmagK(1),divSmag(1),res)
  !Mixage
  div(1)%values = divSmag(1)%values + divGrad(1)%values
  call ftran(div(1),divK,res)
  !Désallocation des tableaux de travail
  deallocate(coefClark)
  if ( .not. deleteWorkArray(divGrad) .or. &
  & .not. deleteWorkArray(divGradK) .or. &
  & .not. deleteWorkArray(divSmag) .or. &
  & .not. deleteWorkArray(divSmagK) .or. &
  & .not. deleteWorkArray(div) .or. &
  & .not. deleteWorkArray(SdZ) .or. &
  & .not. deleteWorkArray(SdZK)) then
    res = .false.
    write(6,'(a)')'[ERROR] in dynClarksca_3d: memory leak'
  endif

 end subroutine DynClarkSca_2d_z
 
!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> Sub-grid model based on DCM one variable
!> @param [in]    U  longitudinal velocity in physical space
!> @param [in]    V  velocity in physical space
!> @param [in]    W  velocity in physical space
!> @param [in]    Uk longitudinal velocity in spectral space
!> @param [in]    Vk velocity in spectral space
!> @param [in]    Wk velocity in spectral space
!> @param [in]    tabEOFile contains the EO model
!> @param [in]    Scal scalar in physical space 
!> @param [in]    Scalk scalar in spectral space 
!> @param [in]    ScalWN wave numbers for scalar 
!> @param [in]    spec_rank is the rank of processus
!> @param [inout] res logical
!> @param [inout] div is the divergence of SGS model
!> @return
!------------------------------------------------------------------------------
 subroutine DynClarkSca_3d(divK,U,V,W,Uk,Vk,Wk,ScalIn,ScalKIn,WaveNum,filter,spec_rank,res)
  !I/O
  type(complex_data_layout),intent(inout) :: divK
  type(complex_data_layout),intent(in)    :: Uk,Vk,Wk,ScalKIn
  type(real_data_layout),intent(in)       :: U, V, W, ScalIn
  type(WaveNumbers), intent(in)           :: WaveNum
  logical, intent(inout)                  :: res
  integer, intent(in)                     :: spec_rank
  integer, intent(in)                     :: filter

  !Local data
  type(real_data_layout),dimension(:),allocatable    :: divSmag
  type(real_data_layout),dimension(:),allocatable    :: divGrad
  type(real_data_layout),dimension(:),allocatable    :: div
  type(real_data_layout),dimension(:),allocatable    :: SdZ
  type(complex_data_layout),dimension(:),allocatable :: divGradK
  type(complex_data_layout),dimension(:),allocatable :: divSmagK
  type(complex_data_layout),dimension(:),allocatable :: SdZK
  real(WP)                                           :: delta,coefClark
 
  res = .true.

  !Init des tableaux de travail
  if (.not. initWorkArray(U,1,divSmag) .or. &
     &.not. initWorkArray(U,1,divGrad) .or. &
     &.not. initWorkArray(U,1,div) .or. &
     &.not. initWorkArray(U,3,SdZ) .or. &
     &.not. initWorkArray(Uk,1,divSmagK) .or. &
     &.not. initWorkArray(Uk,1,divGradK) .or. &
     &.not. initWorkArray(Uk,3,SdZK)) then
    write(6,'(a)')'[ERROR] initWorkArray : not enought memory!'
    return
  endif

  delta = computeDelta(ScalIn)

  !Calcul du coef dynamique de Clark 
  if (.not. coefDynClarkSca(U,V,W,Uk,Vk,Wk,ScalIn,ScalKIn,WaveNum,filter,1,&
     &      spec_rank,CoeffDynScalar=coefClark) ) then
    write (6,'(a)') '[INFO] coef Dynamic failed' 
  endif
  if (.not. coefDynClarkScaTEST(U,V,W,Uk,Vk,Wk,ScalIn,ScalKIn,WaveNum,filter,1,&
     &      spec_rank,CoeffDynScalar=coefClark) ) then
    write (6,'(a)') '[INFO] coef Dynamic test failed' 
  endif
  !CALCUL des TERMES
  !Gradient
  call gradientSca(divGradK(1),U,Uk,Vk,Wk,ScalKIn,ScalIn,WaveNum,res,1)
  if (.not. res ) write(6,'(a)')'[ERROR] in gradient' 
  call btran(divGradK(1),divGrad(1),res)
  if (.not. res ) write(6,'(a)')'[ERROR] in btran' 
  !Smagorinsky
  call computeSdZ(U,Uk,Vk,Wk,ScalKIn,WaveNum,SdZ(1),SdZ(2),SdZ(3),res)
  if (.not. res ) write(6,'(a)')'[ERROR] in computeSdZ' 
  SdZ(1)%values = coefClark*delta**2*SdZ(1)%values
  SdZ(2)%values = coefClark*delta**2*SdZ(2)%values
  SdZ(3)%values = coefClark*delta**2*SdZ(3)%values
  call ftran(SdZ(1),SdZK(1),res)
  if (.not. res ) write(6,'(a)')'[ERROR] in smagorinsky for scalar : compute ftran 1 failed'
  call ftran(SdZ(2),SdZK(2),res)
  if (.not. res ) write(6,'(a)')'[ERROR] in smagorinsky for scalar : compute ftran 2 failed'
  call ftran(SdZ(3),SdZK(3),res)
  if (.not. res ) write(6,'(a)')'[ERROR] in smagorinsky for scalar : compute ftran 3 failed'
  call computeDivergenceK(WaveNum,SdZK(1),SdZK(2),SdZK(3),divSmagK(1))
  call btran(divSmagK(1),divSmag(1),res)
  !Mixage
  div(1)%values = divSmag(1)%values + divGrad(1)%values
  call ftran(div(1),divK,res)
  !Désallocation des tableaux de travail
  if ( .not. deleteWorkArray(divGrad) .or. &
  & .not. deleteWorkArray(divGradK) .or. &
  & .not. deleteWorkArray(divSmag) .or. &
  & .not. deleteWorkArray(divSmagK) .or. &
  & .not. deleteWorkArray(div) .or. &
  & .not. deleteWorkArray(SdZ) .or. &
  & .not. deleteWorkArray(SdZK)) then
    res = .false.
    write(6,'(a)')'[ERROR] in dynClarksca_3d: memory leak'
  endif

 end subroutine DynClarkSca_3d


!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!>  Smagorinsky model 
!> @param [in]    U  longitudinal velocity in physical space
!> @param [in]    Uk longitudinal velocity in spectral space
!> @param [in]    Vk vertical velocity in spectral space
!> @param [in]    Wk spanwize velocity in spectral space
!> @param [in]    Scalk scalar in spectral space 
!> @param [in]    ScalWN wave numbers for scalar 
!> @param [in]    spec_rank is the rank of processus
!> @param [inout] res logical
!> @param [in]    coeff is the dynamic coef (optinal) 
!> @param [out]   div is the divergence of SGS model
!------------------------------------------------------------------------------
 subroutine DynSmagorinskySca_2d_z(div,U,V,W,Scal,Uk,Vk,Wk,ScalK,ScalWN,res,filter,spec_rank)

  !I/O
  type(complex_data_layout),intent(in)          :: Uk,Vk,Wk,ScalK 
  type(complex_data_layout),intent(inout)       :: div
  type(real_data_layout),intent(in)             :: U,V,W,Scal
  type(WaveNumbers), intent(in)                 :: ScalWN
  logical, intent(inout)                        :: res
  integer, intent(in)                           :: spec_rank
  integer, intent(in)                           :: filter

  !Local data
  type(real_data_layout),dimension(:),allocatable    :: temp 
  type(complex_data_layout),dimension(:),allocatable :: tempK
  real(WP)                                           :: deltx
  real(WP),dimension(:),allocatable                  :: coeff
  integer                                            :: ires,k

  res = .true.
  if (.not. initWorkArray(U,3,temp)  .or. &
     &.not. initWorkArray(Uk,3,tempK) ) then
    res = .false.
    write(6,'(a)')'[ERROR] in smagorinskySca : not enought memory'
  endif
  allocate(coeff(Scal%nz),stat=ires)
  if ( ires .ne. 0 ) then
    write (6,'(a)') '[ERROR] not enought memory to allocate coeff'
    res = .false.
  endif 
  deltx = computeDelta(U)
  call computeSdZ(U,Uk,Vk,Wk,ScalK,ScalWN,temp(1),temp(2),temp(3),res)
  if (.not. coefDynSmagSca(U,V,W,Uk,Vk,Wk,Scal,ScalK,ScalWN,filter,1,spec_rank,&
     & CoeffDynArray=coeff) ) then
    write (6,'(a)') '[ERROR] in DynSmagorinskySca_3d'
    res = .false.
  endif
  if (.not. coefDynSmagScaTEST(U,V,W,Uk,Vk,Wk,Scal,ScalK,ScalWN,filter,1,spec_rank,&
     & CoeffDynArray=coeff) ) then
    write (6,'(a)') '[ERROR] in DynSmagorinskySca_3d'
    res = .false.
  endif
  do k = U%zmin,U%zmax
     temp(1)%values(:,:,k) = coeff(k) * deltx**2.0 * temp(1)%values(:,:,k)
     temp(2)%values(:,:,k) = coeff(k) * deltx**2.0 * temp(2)%values(:,:,k)
     temp(3)%values(:,:,k) = coeff(k) * deltx**2.0 * temp(3)%values(:,:,k)
  end do
  call ftran(temp(1),tempK(1),res)
  if (.not. res ) write(6,'(a)')'[ERROR] in smagorinsky for scalar : compute ftran 1 failed'
  call ftran(temp(2),tempK(2),res)
  if (.not. res ) write(6,'(a)')'[ERROR] in smagorinsky for scalar : compute ftran 2 failed'
  call ftran(temp(3),tempK(3),res)
  if (.not. res ) write(6,'(a)')'[ERROR] in smagorinsky for scalar : compute ftran 3 failed'
  call computeDivergenceK(ScalWN,tempK(1),tempK(2),tempK(3),div)

  deallocate(coeff)
  if (.not. deleteWorkArray(temp)  .or. &
     &.not. deleteWorkArray(tempK) ) then
    res = .false.
    write(6,'(a)')'[ERROR] in smagorinskySca : memory leak'
  endif

 end subroutine DynSmagorinskySca_2d_z

!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!>  Smagorinsky model 
!> @param [in]    U  longitudinal velocity in physical space
!> @param [in]    Uk longitudinal velocity in spectral space
!> @param [in]    Vk vertical velocity in spectral space
!> @param [in]    Wk spanwize velocity in spectral space
!> @param [in]    Scalk scalar in spectral space 
!> @param [in]    ScalWN wave numbers for scalar 
!> @param [in]    spec_rank is the rank of processus
!> @param [inout] res logical
!> @param [in]    coeff is the dynamic coef (optinal) 
!> @param [out]   div is the divergence of SGS model
!------------------------------------------------------------------------------
 subroutine DynSmagorinskySca_3d(div,U,V,W,Scal,Uk,Vk,Wk,ScalK,ScalWN,res,filter,spec_rank)

  !I/O
  type(complex_data_layout),intent(in)          :: Uk,Vk,Wk,ScalK 
  type(complex_data_layout),intent(inout)       :: div
  type(real_data_layout),intent(in)             :: U,V,W,Scal
  type(WaveNumbers), intent(in)                 :: ScalWN
  logical, intent(inout)                        :: res
  integer, intent(in)                           :: spec_rank
  integer, intent(in)                           :: filter

  !Local data
  type(real_data_layout),dimension(:),allocatable    :: temp 
  type(complex_data_layout),dimension(:),allocatable :: tempK
  real(WP)                                           :: deltx,coeff
  res = .true.

  if (.not. initWorkArray(U,3,temp)  .or. &
     &.not. initWorkArray(Uk,3,tempK) ) then
    res = .false.
    write(6,'(a)')'[ERROR] in smagorinskySca : not enought memory'
  endif

  deltx = computeDelta(U)
  call computeSdZ(U,Uk,Vk,Wk,ScalK,ScalWN,temp(1),temp(2),temp(3),res)
  if (.not. coefDynSmagSca(U,V,W,Uk,Vk,Wk,Scal,ScalK,ScalWN,filter,1,spec_rank,&
     & CoeffDynScalar=coeff) ) then
    write (6,'(a)') '[ERROR] in DynSmagorinskySca_3d'
    res = .false.
  endif
  if (.not. coefDynSmagScaTEST(U,V,W,Uk,Vk,Wk,Scal,ScalK,ScalWN,filter,1,spec_rank,&
     & CoeffDynScalar=coeff) ) then
    write (6,'(a)') '[ERROR] in DynSmagorinskySca_3d'
    res = .false.
  endif
  temp(1)%values = coeff * deltx**2 * temp(1)%values
  temp(2)%values = coeff * deltx**2 * temp(2)%values
  temp(3)%values = coeff * deltx**2 * temp(3)%values

  call ftran(temp(1),tempK(1),res)
  if (.not. res ) write(6,'(a)')'[ERROR] in smagorinsky for scalar : compute ftran 1 failed'
  call ftran(temp(2),tempK(2),res)
  if (.not. res ) write(6,'(a)')'[ERROR] in smagorinsky for scalar : compute ftran 2 failed'
  call ftran(temp(3),tempK(3),res)
  if (.not. res ) write(6,'(a)')'[ERROR] in smagorinsky for scalar : compute ftran 3 failed'
  call computeDivergenceK(ScalWN,tempK(1),tempK(2),tempK(3),div)

  if (.not. deleteWorkArray(temp)  .or. &
     &.not. deleteWorkArray(tempK) ) then
    res = .false.
    write(6,'(a)')'[ERROR] in DynsmagorinskySca_3d : memory leak'
  endif

 end subroutine DynSmagorinskySca_3d 
 
 
 !!!!!!!!!!!!!!TEST GB 26/10/2012
! ------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!>  Smagorinsky model 
!> @param [in]    U  longitudinal velocity in physical space
!> @param [in]    Uk longitudinal velocity in spectral space
!> @param [in]    Vk vertical velocity in spectral space
!> @param [in]    Wk spanwize velocity in spectral space
!> @param [in]    Scalk scalar in spectral space 
!> @param [in]    ScalWN wave numbers for scalar 
!> @param [in]    spec_rank is the rank of processus
!> @param [inout] res logical
!> @param [in]    coeff is the dynamic coef (optinal) 
!> @param [out]   div is the divergence of SGS model
!------------------------------------------------------------------------------
 subroutine smagorinskySca2(div,U,Uk,Vk,Wk,ScalK,ScalWN,res,spec_rank,dyncoeff,zmin,zmax)

  !I/O
  type(complex_data_layout),intent(in)    :: Uk,Vk,Wk,ScalK 
  type(complex_data_layout),intent(inout) :: div
  type(real_data_layout),intent(in)       :: U
  type(WaveNumbers), intent(in)           :: ScalWN
  logical, intent(inout)                  :: res
  integer, intent(in)                     :: zmin,zmax
  real(WP), intent(in),dimension(zmin:zmax)    :: dyncoeff
  integer, intent(in)                     :: spec_rank

  !Local data
  type(real_data_layout),dimension(:),allocatable    :: temp, temp4 
  type(complex_data_layout),dimension(:),allocatable :: tempK
  real(WP)                                           :: deltx,coeff
  integer  :: i,j,k
  res = .true.

  if (.not. initWorkArray(U,3,temp)  .or. &
     &.not. initWorkArray(Uk,3,tempK) ) then
    res = .false.
    write(6,'(a)')'[ERROR] in smagorinskySca : not enought memory'
  endif

  deltx = computeDelta(U)
  call computeSdZ(U,Uk,Vk,Wk,ScalK,ScalWN,temp(1),temp(2),temp(3),res)

  do k = U%zmin, U%zmax
     temp(1)%values(:,:,k) = dyncoeff(k) * deltx**2.0 * temp(1)%values(:,:,k)
     temp(2)%values(:,:,k) = dyncoeff(k) * deltx**2.0 * temp(2)%values(:,:,k)
     temp(3)%values(:,:,k) = dyncoeff(k) * deltx**2.0 * temp(3)%values(:,:,k)
  end do
  
  res = initWorkArray(U,3,temp4)
  call computeGradientK(ScalWN,ScalK,tempK(1),tempK(2),tempK(3))
  call btran(tempK(1),temp4(1),res)
  call btran(tempK(2),temp4(2),res)
  call btran(tempK(3),temp4(3),res)
  temp4(1)%values = temp(1)%values * temp4(1)%values + temp(2)%values * temp4(2)%values +&
                  & temp(3)%values * temp4(3)%values
  do k = U%zmin, U%zmax
     do j = U%ymin, U%ymax
        do i = U%xmin, U%xmax
           if (temp4(1)%values(i,j,k).gt.0.0_WP) then
              print*,'DynSmag something wrong', temp4(1)%values(i,j,k), dyncoeff(zmin)
           end if
        end do
     end do
  end do
  coeff = computeFieldAvg( temp4(1) , spec_rank )
  if (spec_rank.eq.0) print*,'DynSmag Dissipation =',coeff
  res = deleteWorkArray(temp4)

  call ftran(temp(1),tempK(1),res)
#ifdef DEBUGLES
  if (.not. res ) write(6,'(a)')'[ERROR] in smagorinsky for scalar : compute ftran 1 failed'
#endif
  call ftran(temp(2),tempK(2),res)
#ifdef DEBUGLES
  if (.not. res ) write(6,'(a)')'[ERROR] in smagorinsky for scalar : compute ftran 2 failed'
#endif
  call ftran(temp(3),tempK(3),res)
#ifdef DEBUGLES
  if (.not. res ) write(6,'(a)')'[ERROR] in smagorinsky for scalar : compute ftran 3 failed'
#endif

  call computeDivergenceK(ScalWN,tempK(1),tempK(2),tempK(3),div)

!#ifdef MEMOPTIM
  if (.not. deleteWorkArray(temp)  .or. &
     &.not. deleteWorkArray(tempK) ) then
    res = .false.
    write(6,'(a)')'[ERROR] in smagorinskySca : memory leak'
  endif
!#endif

 end subroutine smagorinskySca2 
 
 
 
 
 
 !------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> Compute the dynamic coef for Smagorinsky dynamic 
!> @param [in] U is the first component of velocity in physical space and filtered
!> @param [in] V is the second component of velocity in physical space and filtered
!> @param [in] W is the third component of velocity in physical space and filtered
!> @param [in] Scal is the scalar in physical space and filtered
!> @param [in] Uk is the first component of velocity in Fourier space and filtered
!> @param [in] Vk is the second component of velocity in Fourier space and filtered
!> @param [in] Wk is the third component of velocity in Fourier space and filtered
!> @param [in] Scalk is the scalar in Fourier space and filtered
!> @param [in] ScalWN is the waves number of Scalar 
!> @param [in] filter is the type of test filter 
!> @param [out] CoeffDyn is dynamic coefficient which is computed by the routine 
!> @param [in] spec_rank is the number of process in the cpu pool
!> @param [out] res is a logical value which is equal to .true. if everything is fine and .false otherwize 
!> @param [out] numer is an optinal field aqual to the numerator of Lilly's term <LiMi> (optional) 
!> @param [out] denom is an optinal field aqual to the denominator of Lilly's term <MiMi> (optional)
!> @param [in] filterSize is the of the first filter in the a priori case (optional) 
!------------------------------------------------------------------------------
 subroutine coefDynSca2(U,V,W,Uk,Vk,Wk,Scal,ScalK,ScalWN,filter,CoeffDyn,zmin,zmax,type_avg,spec_rank,res,numer,denom,filterSize)

  !I/O
  type(complex_data_layout),intent(in)    :: Uk,Vk,Wk,ScalK 
  type(real_data_layout),intent(in)       :: U,V,W,Scal
  type(real_data_layout),intent(inout),optional :: numer,denom
  type(WaveNumbers), intent(in)           :: ScalWN
  integer, intent(in)                     :: spec_rank,filter, type_avg,zmin,zmax
  logical, intent(inout)                  :: res 
  real(WP), intent(inout),dimension(zmin:zmax) :: coeffDyn
  real(WP), intent(in),optional           :: filterSize 

  !Local data
  type(real_data_layout),dimension(:),allocatable     :: temp 
  type(complex_data_layout),dimension(:),allocatable  :: tempK  
  real(WP)                                            :: deltx, coeff
  real(WP)                                            :: sizeTestFilter
  real(WP), dimension(:),allocatable                  :: num1D, den1D
  

  res = .true.

  if (.not. initWorkArray(U,6,temp) .or. &
     &.not. initWorkArray(Uk,6,tempK )) then
    res = .false.
    write(6,'(a)')'[ERROR] in coefDynSca2 for scalar : not enought memory'
  endif

  sizeTestFilter  = 2.0_WP

  !Etape de diffÃ©renciation a priori/ aposteriori
  if ( present(filterSize) ) then
    deltx = filterSize
  else
    deltx = computeDelta(U)
  endif

  !Filtering U
  call computeFilter(ScalWN,sizeTestFilter*deltx,Uk,tempK(1),filter) 
  call computeFilter(ScalWN,sizeTestFilter*deltx,Vk,tempK(2),filter) 
  call computeFilter(ScalWN,sizeTestFilter*deltx,Wk,tempK(3),filter) 
  !Filtering Sca 
  call computeFilter(ScalWN,sizeTestFilter*deltx,ScalK,tempK(4),filter) 

  !Compute dZtilde/dxi
  call computeSdZ(U,tempK(1),tempK(2),tempK(3),tempK(4),ScalWN,temp(1),temp(2),temp(3),res)
  !    computeSdZ(U,UK      ,VK      ,WK      ,ScalK   ,ScalWN,SdZ1   ,SdZ2   ,SdZ3   ,res)

  !Compute dZ/dxi
  call computeSdZ(U,Uk,Vk,Wk,ScalK,ScalWN,temp(4),temp(5),temp(6),res)
  !    computeSdZ(U,UK,VK,WK,ScalK,ScalWN,SdZ1   ,SdZ2   ,SdZ3   ,res)

  call ftran(temp(4),tempK(1),res)
  call ftran(temp(5),tempK(2),res)
  call ftran(temp(6),tempK(3),res)

  call computeFilter(ScalWN,sizeTestFilter*deltx,tempK(1),tempK(4),filter) 
  call computeFilter(ScalWN,sizeTestFilter*deltx,tempK(2),tempK(5),filter) 
  call computeFilter(ScalWN,sizeTestFilter*deltx,tempK(3),tempK(6),filter) 
  
  call btran(tempK(4),temp(4),res)
  call btran(tempK(5),temp(5),res)
  call btran(tempK(6),temp(6),res)

  temp(1)%values = deltx**2.0 *( sizeTestFilter**2.0 * temp(1)%values - temp(4)%values )
  temp(2)%values = deltx**2.0 *( sizeTestFilter**2.0 * temp(2)%values - temp(5)%values )
  temp(3)%values = deltx**2.0 *( sizeTestFilter**2.0 * temp(3)%values - temp(6)%values )
!  temp(1)%values = deltx**2.0 *( sizeTestFilter**2.0 * temp(4)%values - temp(1)%values )
!  temp(2)%values = deltx**2.0 *( sizeTestFilter**2.0 * temp(5)%values - temp(2)%values )
!  temp(3)%values = deltx**2.0 *( sizeTestFilter**2.0 * temp(6)%values - temp(3)%values )

  !compute Germano's term
  call computeGermanoSca(ScalWN,U,V,W,Scal,Uk,Vk,Wk,Scalk,Filter,deltx,temp(4),temp(5),temp(6))
  !    computeGermanoSca(ScalWN,U,V,W,Scal,Uk,Vk,Wk,Scalk,Filter,deltx,S12    ,S13    ,S23)

  !compute Lilly's least squares method
  call computeLilly(temp(1),temp(2),temp(3),temp(4),temp(5),temp(6),Coeff,spec_rank)
  !    computeLilly(S11    ,S22    ,S33    ,S12    ,S13    ,S23    ,CoeffDyn,spec_rank,nbcpus)
  
  if (type_avg.eq.0) then
  
     coeffDyn(:) = coeff
     
  else if (type_avg.eq.1) then
  
     allocate(num1D(zmin:zmax))
     allocate(den1D(zmin:zmax))
     
     temp(4)%values = temp(1)%values*temp(4)%values + temp(2)%values*temp(5)%values + &
                  & temp(3)%values*temp(6)%values
     temp(1)%values = temp(1)%values*temp(1)%values + temp(2)%values*temp(2)%values + &
                  & temp(3)%values*temp(3)%values
     res = computeFieldAvg2D(temp(4),num1D,spec_rank)
     res = computeFieldAvg2D(temp(1),den1D,spec_rank)
     coeffdyn(:) = num1D(:) / den1D(:)
       
     deallocate(num1D)
     deallocate(den1D)
     
  else
    
       print*,'[ERROR] AVERGAGING NOT KNOWN FOR DYNAMIC PROCEDURE'
       
  end if
     
     
     

  if (present(numer) .and. present(denom)) then
     numer%values = temp(1)%values*temp(4)%values + temp(2)%values*temp(5)%values + &
                  & temp(3)%values*temp(6)%values
     denom%values = temp(1)%values*temp(1)%values + temp(2)%values*temp(2)%values + &
                  & temp(3)%values*temp(3)%values
  endif

  if (.not. deleteWorkArray(temp) .or. &
     &.not. deleteWorkArray(tempK ) ) then
    res = .false.
    write(6,'(a)')'[ERROR] in coefDynSca2 for scalar : memory leak is occuring... '
  endif

 end subroutine coefDynSca2
 


end module scamodels
!> @}
