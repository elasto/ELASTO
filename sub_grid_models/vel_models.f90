!> @addtogroup sub_grid_models
!! @{

module velmodels


 use datalayout
 use toolbox
 use stat_tools 
 use differential_tools 
 use transforms_tools
 use wavenumber_tools 
 use subgrid_tools
 use filtering_tools 
 use physical_values_tools

 implicit none

 interface EVSGSMVreman
   module procedure EVSGSMVreman_tij
   module procedure EVSGSMVreman_div
 end interface EVSGSMVreman 

 interface RGMSikmoinsSkj
   module procedure RGMSikmoinsSkj_div
   module procedure RGMSikmoinsSkj_tij
 end interface RGMSikmoinsSkj 

 interface RGMSikOjkOikSjk 
   module procedure RGMSikOjkOikSjk_div
   module procedure RGMSikOjkOikSjk_tij
 end interface RGMSikOjkOikSjk

 interface RGMOikOjk
   module procedure RGMOikOjk_div
   module procedure RGMOikOjk_tij
 end interface RGMOikOjk

 interface RGMSikmoinsSkjSikOjkOikSjk
   module procedure RGMSikmoinsSkjSikOjkOikSjk_div
   module procedure RGMSikmoinsSkjSikOjkOikSjk_tij
 end interface RGMSikmoinsSkjSikOjkOikSjk

 interface RGMSikmoinsSkjOikOjk 
   module procedure RGMSikmoinsSkjOikOjk_div
   module procedure RGMSikmoinsSkjOikOjk_tij
 end interface RGMSikmoinsSkjOikOjk


 contains

!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!>  Smgorinsky model 
!> @param [in]
!> @param [out]
!------------------------------------------------------------------------------
 subroutine smagorinskyVel(div1,div2,div3,U,V,W,Uk,Vk,Wk,VelWN,numVel,spec_rank,res, &
                          &  numVelFilter,type_avg,T11,T12,T13,T22,T23,T33,delta,ite)


  !I/O
  type(complex_data_layout),intent(in)    :: Uk,Vk,Wk 
  type(complex_data_layout),intent(inout) :: div1,div2,div3
  type(real_data_layout),intent(in)       :: U,V,W
  type(WaveNumbers), intent(in)           :: VelWN
  integer, intent(in)                     :: spec_rank, numVel
  logical, intent(inout)                  :: res
  integer, intent(in),optional            :: numVelFilter, type_avg, ite
  type(real_data_layout),intent(inout),optional         :: T11,T12,T13,T22,T23,T33  
  real(WP),intent(in),optional            :: delta
  !Local data
  type(real_data_layout)          :: S11,S12,S13  
  type(real_data_layout)          :: S22,S23,S33  
  type(real_data_layout)          :: S11f,S12f,S13f  
  type(real_data_layout)          :: S22f,S23f,S33f  
  type(real_data_layout)          :: Uf, Vf, Wf
  type(real_data_layout)          :: UUf, UVf, UWf, VVf, VWf, WWf
  type(real_data_layout)          :: normS11f,normS12f,normS13f  
  type(real_data_layout)          :: normS22f,normS23f,normS33f  
  type(real_data_layout)          :: coefff
  type(complex_data_layout)       :: S11K,S12K,S13K  
  type(complex_data_layout)       :: S22K,S23K,S33K
  type(complex_data_layout),dimension(:),allocatable :: tempK, tempKf
  type(real_data_layout)          :: norm, normf
  real(WP), dimension(:),allocatable :: num1D, den1D
  real(WP)                        :: deltx, num, den
  integer                         :: k
  logical                         :: pdffield = .false.

 res = .true.
 if (present(ite)) then
   if ( ite .ne. -1 ) pdffield = .true.
 endif

! call cpu_time(tps1)
  if (.not. copyStructOnly(U,S11) .or. &
     &.not. copyStructOnly(U,S12) .or. &
     &.not. copyStructOnly(U,S13) .or. &
     &.not. copyStructOnly(U,S22) .or. &
     &.not. copyStructOnly(U,S23) .or. &
     &.not. copyStructOnly(U,S33) .or. &
     &.not. copyStructOnly(U,coefff) .or. &
     &.not. copyStructOnly(U,norm) )then
       res = .false.
       write(6,'(a,i0)')'[ERROR] in smagorinsky for vel : not enought memory !',spec_rank
   endif
! call cpu_time(tps2)
!print *,'temps mis pour allocation:',tps2-tps1

  call computeStressTensor(Uk,VK,WK,VelWN,S11,S12,S13,S22,S23,S33,res)
  call computeTensorNorme1(S11,S12,S13,S22,S23,S33,norm,res)
  deltx = computeDelta(U)
  if ((present(T11)).AND.(present(T22)).AND.(present(T33)).AND.&
     &(present(T12)).AND.(present(T13)).AND.(present(T23)) ) then
    if ( .not. present(delta) ) write(6,'(a,i0)')'[ERROR] in smagorinsky for vel, delta is not specified !'
    deltx=delta
  endif
  if (present(delta)) then
    deltx=delta
  endif

  if (numVel.eq.1) then
     ! Static Smagorinsky model
     coefff%values = -2.0_WP * (0.18_WP**2)
  else if (numVel.eq.2) then
     ! Dynamic model --- Formulation
     ! Leonard Terms computation -> (u_i u_j)f - uf_i uf_j
     if ( .not. present(numVelFilter) ) write(6,'(a,i0)')'[ERROR] in dyn smagorinsky for vel, type of filter not specified!'
     if (.not. copyStructOnly(U,Uf) .or. &
        &.not. copyStructOnly(U,Vf) .or. &
        &.not. copyStructOnly(U,Wf) .or. &
        &.not. copyStructOnly(U,UUf) .or. &
        &.not. copyStructOnly(U,UVf) .or. &
        &.not. copyStructOnly(U,UWf) .or. &
        &.not. copyStructOnly(U,VVf) .or. &
        &.not. copyStructOnly(U,VWf) .or. &
        &.not. copyStructOnly(U,WWf) )then
       res = .false.
       write(6,'(a,i0)')'[ERROR] in dyn smagorinsky for vel : not enought memory !',spec_rank
     endif
     res = initWorkArray(Uk,3,tempK)
     res = initWorkArray(Uk,3,tempKf)

     ! uf_i
     call computeFilter(VelWN,2.0*deltx,Uk,tempKf(1),numVelFilter) 
     call computeFilter(VelWN,2.0*deltx,Vk,tempKf(2),numVelFilter) 
     call computeFilter(VelWN,2.0*deltx,Wk,tempKf(3),numVelFilter) 
     call btran(tempKf(1),Uf,res) 
     call btran(tempKf(2),Vf,res) 
     call btran(tempKf(3),Wf,res)

     call computeLijVec(U,V,W,Uf,Vf,Wf&
                     &,Uk,UVf,UWf,VWf,UUf,VVf,WWf,deltx,VelWN,numVelFilter)
 
     ! Compute M_ij
     ! ( norme * Sij )f
      if (.not. copyStructOnly(U,normS11f) .or. &
        &.not. copyStructOnly(U,normS12f) .or. &
        &.not. copyStructOnly(U,normS13f) .or. &
        &.not. copyStructOnly(U,normS22f) .or. &
        &.not. copyStructOnly(U,normS23f) .or. &
        &.not. copyStructOnly(U,normS33f) ) then
       res = .false.
       write(6,'(a,i0)')'[ERROR] in dyn smagorinsky for vel 2 : not enought memory !',spec_rank
     endif
     normS11f%values = norm%values * S11%values
     normS12f%values = norm%values * S12%values
     normS13f%values = norm%values * S13%values
     normS22f%values = norm%values * S22%values
     normS23f%values = norm%values * S23%values
     normS33f%values = norm%values * S33%values
     call ftran(normS11f,tempK(1),res)
     call ftran(normS12f,tempK(2),res)
     call ftran(normS13f,tempK(3),res)
     call computeFilter(VelWN,2.0*deltx,tempK(1),tempKf(1),numVelFilter) 
     call computeFilter(VelWN,2.0*deltx,tempK(2),tempKf(2),numVelFilter) 
     call computeFilter(VelWN,2.0*deltx,tempK(3),tempKf(3),numVelFilter) 
     call btran(tempKf(1),normS11f,res) 
     call btran(tempKf(2),normS12f,res) 
     call btran(tempKf(3),normS13f,res)
     call ftran(normS22f,tempK(1),res)
     call ftran(normS23f,tempK(2),res)
     call ftran(normS33f,tempK(3),res)
     call computeFilter(VelWN,2.0*deltx,tempK(1),tempKf(1),numVelFilter) 
     call computeFilter(VelWN,2.0*deltx,tempK(2),tempKf(2),numVelFilter) 
     call computeFilter(VelWN,2.0*deltx,tempK(3),tempKf(3),numVelFilter) 
     call btran(tempKf(1),normS22f,res) 
     call btran(tempKf(2),normS23f,res) 
     call btran(tempKf(3),normS33f,res)
     ! normef * S11f
     if (.not. copyStructOnly(U,S11f) .or. &
        &.not. copyStructOnly(U,S12f) .or. &
        &.not. copyStructOnly(U,S13f) .or. &
        &.not. copyStructOnly(U,S22f) .or. &
        &.not. copyStructOnly(U,S23f) .or. &
        &.not. copyStructOnly(U,S33f) .or. &
        &.not. copyStructOnly(U,normf) ) then
       res = .false.
       write(6,'(a,i0)')'[ERROR] in dyn smagorinsky for vel 2 : not enought memory !',spec_rank
     endif
     call computeFilter(VelWN,2.0*deltx,Uk,tempKf(1),numVelFilter) 
     call computeFilter(VelWN,2.0*deltx,Vk,tempKf(2),numVelFilter) 
     call computeFilter(VelWN,2.0*deltx,Wk,tempKf(3),numVelFilter) 
     call computeStressTensor(tempKf(1),tempKf(2),tempKf(3),VelWN,S11f,S12f,S13f,S22f,S23f,S33f,res)
     call computeTensorNorme1(S11f,S12f,S13f,S22f,S23f,S33f,normf,res)
     normS11f%values = (2.0 * deltx)**2.0 * normf%values * S11f%values - ( deltx )**2.0 * normS11f%values
     normS12f%values = (2.0 * deltx)**2.0 * normf%values * S12f%values - ( deltx )**2.0 * normS12f%values
     normS13f%values = (2.0 * deltx)**2.0 * normf%values * S13f%values - ( deltx )**2.0 * normS13f%values
     normS22f%values = (2.0 * deltx)**2.0 * normf%values * S22f%values - ( deltx )**2.0 * normS22f%values
     normS23f%values = (2.0 * deltx)**2.0 * normf%values * S23f%values - ( deltx )**2.0 * normS23f%values
     normS33f%values = (2.0 * deltx)**2.0 * normf%values * S33f%values - ( deltx )**2.0 * normS33f%values
     ! L_ij * M_ij
     UUf%values  = UUf%values*normS11f%values + VVf%values*normS22f%values + WWf%values*normS33f%values + &
         &  2.0 * (UVf%values*normS12f%values + UWf%values*normS13f%values + VWf%values*normS23f%values )
     VVf%values  = normS11f%values**2.0 + normS22f%values**2.0 + normS33f%values**2.0 + &
         &  2.0 * (normS12f%values**2.0 + normS13f%values**2.0 + normS23f%values**2.0 )
     
     if (type_avg.eq.0) then
        num = computeFieldAvg( UUf , spec_rank )
        den = computeFieldAvg( VVf , spec_rank )
        coefff%values = num / den
        if (spec_rank.eq.0) print*,'DynCoef Vel : ',num / den
     else if (type_avg.eq.1) then
         allocate(num1D(1:U%nz))
         allocate(den1D(1:U%nz))
         num1D = 0.0
         den1D = 0.0
         if(.not. computeFieldAvg2D(UUf,num1D,spec_rank)) write(*,'(a)')&
                & '[ERROR] Could not compute AVG for UUf'
         if(.not. computeFieldAvg2D(VVf,den1D,spec_rank)) write(*,'(a)')&
                & '[ERROR] Could not compute AVG for VVf'
         do k = U%zmin, U%zmax
            coefff%values(:,:,k) = num1d(k)/den1D(k)
         end do
         deallocate(num1D)
         deallocate(den1D)
     end if
     
     call deleteDataLayout(Uf)
     call deleteDataLayout(Vf)
     call deleteDataLayout(Wf)
     call deleteDataLayout(UUf)
     call deleteDataLayout(UVf)
     call deleteDataLayout(UWf)
     call deleteDataLayout(VVf)
     call deleteDataLayout(VWf)
     call deleteDataLayout(WWf)
     call deleteDataLayout(normS11f)
     call deleteDataLayout(normS12f)
     call deleteDataLayout(normS13f)
     call deleteDataLayout(normS22f)
     call deleteDataLayout(normS23f)
     call deleteDataLayout(normS33f)
     call deleteDataLayout(S11f)
     call deleteDataLayout(S12f)
     call deleteDataLayout(S13f)
     call deleteDataLayout(S22f)
     call deleteDataLayout(S23f)
     call deleteDataLayout(S33f)
     call deleteDataLayout(normf)
     res = deleteWorkArray(tempK)
     res = deleteWorkArray(tempKf)
  end if
     
  
  S11%values = coefff%values * (deltx)**2.0 * norm%values * S11%values
  S12%values = coefff%values * (deltx)**2.0 * norm%values * S12%values
  S13%values = coefff%values * (deltx)**2.0 * norm%values * S13%values
  S22%values = coefff%values * (deltx)**2.0 * norm%values * S22%values
  S23%values = coefff%values * (deltx)**2.0 * norm%values * S23%values
  S33%values = coefff%values * (deltx)**2.0 * norm%values * S33%values

! call cpu_time(tps1)
  if (.not. copyStructOnly(Uk,S11K) .or. &
     &.not. copyStructOnly(Uk,S12K) .or. &
     &.not. copyStructOnly(Uk,S13K) .or. &
     &.not. copyStructOnly(Uk,S22K) .or. &
     &.not. copyStructOnly(Uk,S23K) .or. &
     &.not. copyStructOnly(Uk,S33K) ) then
       res = .false.
       write(6,'(a,i0)')'[ERROR] in smagorinsky for vel : not enought memory !',spec_rank
  endif
! call cpu_time(tps2)
!print *,'temps mis pour allocation:',tps2-tps1

  call ftran(S11,S11K,res)
  call ftran(S12,S12K,res)
  call ftran(S13,S13K,res)
  call ftran(S22,S22K,res)
  call ftran(S23,S23K,res)
  call ftran(S33,S33K,res)
!     num = computeFieldAvg( S11 , spec_rank )
!   if (spec_rank.EQ.0) print*, 'AvgS11:', num
!     num = computeFieldAvg( S12 , spec_rank )
!   if (spec_rank.EQ.0) print*, 'AvgS12:', num
!     num = computeFieldAvg( S13 , spec_rank )
!   if (spec_rank.EQ.0) print*, 'AvgS13:', num
! 
!     num = computeFieldAvg( S22 , spec_rank )
!   if (spec_rank.EQ.0) print*, 'AvgS22:', num
!     num = computeFieldAvg( S23 , spec_rank )
!   if (spec_rank.EQ.0) print*, 'AvgS23:', num
!     num = computeFieldAvg( S33 , spec_rank )
!   if (spec_rank.EQ.0) print*, 'AvgS33:', num
if ((present(T11)).AND.(present(T22)).AND.(present(T33)).AND.&
   &(present(T12)).AND.(present(T13)).AND.(present(T23)) ) then

  T11%values =  S11%values
  T12%values =  S12%values
  T13%values =  S13%values
  T22%values =  S22%values
  T23%values =  S23%values
  T33%values =  S33%values

endif

  call computeDivergenceK(VelWN,S11K,S12K,S13K,div1)
  call computeDivergenceK(VelWN,S12K,S22K,S23K,div2)
  call computeDivergenceK(VelWN,S13K,S23K,S33K,div3)
  if ( pdffield ) then
    call btran(div1,S11,res)
    call btran(div2,S22,res)
    call btran(div3,S33,res)
    S11%name="sgs_Smag_T1"
    S22%name="sgs_Smag_T2"
    S33%name="sgs_Smag_T3"
    if (.not. computeFieldPDF(ite,S11,300,spec_rank,10000)) return
    if (.not. computeFieldPDF(ite,S22,300,spec_rank,10000)) return
    if (.not. computeFieldPDF(ite,S33,300,spec_rank,10000)) return
  endif

  call deleteDataLayout(S11)
  call deleteDataLayout(S12)
  call deleteDataLayout(S13)
  call deleteDataLayout(S22)
  call deleteDataLayout(S23)
  call deleteDataLayout(S33)
  call deleteDataLayout(norm)
  call deleteDataLayout(coefff)
  call deleteDataLayout(S11K)
  call deleteDataLayout(S12K)
  call deleteDataLayout(S13K)
  call deleteDataLayout(S22K)
  call deleteDataLayout(S23K)
  call deleteDataLayout(S33K)

 end subroutine smagorinskyVel

!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!>  Structure function model 
!> @param [in]
!> @param [out]
!------------------------------------------------------------------------------
 subroutine SfVel(div1,div2,div3,U,V,W,Uk,Vk,Wk,VelWN,spec_rank,numVel,res)


  !I/O
  type(complex_data_layout),intent(inout) :: div1,div2,div3
  type(complex_data_layout),intent(in)    :: Uk,Vk,Wk 
  type(real_data_layout),intent(in)       :: U,V,W
  type(WaveNumbers), intent(in)           :: VelWN
  integer, intent(in)                     :: spec_rank, numVel
  logical, intent(inout)                  :: res

  !Local data
  type(real_data_layout)          :: S11,S12,S13  
  type(real_data_layout)          :: S22,S23,S33  
  type(complex_data_layout)       :: S11K,S12K,S13K  
  type(complex_data_layout)       :: S22K,S23K,S33K
  type(real_data_layout)          :: SF
  type(real_data_layout)          :: Ubuf, Vbuf, Wbuf
  type(real_data_layout)          :: Ubuf2, Vbuf2, Wbuf2
  real(WP)                        :: deltx, coeff
  integer                         :: i,j,k,nlap
  integer                         :: ip, im, jm, jp, km, kp

 res = .true.
 if (.not. copyStructOnly(U,Ubuf) .or. &
    &.not. copyStructOnly(U,Vbuf) .or. &
    &.not. copyStructOnly(U,Wbuf) .or. &
    &.not. copyStructOnly(U,SF) )then
      res = .false.
      write(6,'(a,i0)')'[ERROR] in Structure Function model for vel : not enought memory !',spec_rank
 endif
 
  Ubuf = U
  Vbuf = V
  Wbuf = W
  
  if (numVel.eq.4) then
  ! FSF model - Compute Laplacian Filter
     if (.not. copyStructOnly(U,Ubuf2) .or. &
        &.not. copyStructOnly(U,Vbuf2) .or. &
        &.not. copyStructOnly(U,Wbuf2) )then
           res = .false.
           write(6,'(a,i0)')'[ERROR] in Structure Function model for vel : not enought memory !',spec_rank
     endif
     do nlap = 1,3
        ! Along X  --> initial shape
        do i = Ubuf%xmin, Ubuf%xmax
           im = i-1
           ip = i+1
           if (i.eq.1) im = Ubuf%nx
           if (i.eq.Ubuf%nx) ip = 1
           Ubuf2%values (i,:,:) = -6.0_WP*Ubuf%values (i,:,:) + Ubuf%values(ip,:,:) + Ubuf%values(im,:,:)
           Vbuf2%values (i,:,:) = -6.0_WP*Vbuf%values (i,:,:) + Vbuf%values(ip,:,:) + Vbuf%values(im,:,:)
           Wbuf2%values (i,:,:) = -6.0_WP*Wbuf%values (i,:,:) + Wbuf%values(ip,:,:) + Wbuf%values(im,:,:)
        end do    
        ! Along Z 
        IF (.NOT. setAlongZ(spec_rank,Ubuf)) RETURN  
        IF (.NOT. setAlongZ(spec_rank,Vbuf)) RETURN  
        IF (.NOT. setAlongZ(spec_rank,Wbuf)) RETURN  
        IF (.NOT. setAlongZ(spec_rank,Ubuf2)) RETURN  
        IF (.NOT. setAlongZ(spec_rank,Vbuf2)) RETURN  
        IF (.NOT. setAlongZ(spec_rank,Wbuf2)) RETURN  
        do k = Ubuf%zmin, Ubuf%zmax
           km = k-1
           kp = k+1
           if (k.eq.1) km = Ubuf%nz
           if (k.eq.Ubuf%nz) kp = 1
           Ubuf2%values (:,:,k) = Ubuf2%values (:,:,k) + Ubuf%values(:,:,kp) + Ubuf%values(:,:,km)
           Vbuf2%values (:,:,k) = Vbuf2%values (:,:,k) + Vbuf%values(:,:,kp) + Vbuf%values(:,:,km)
           Wbuf2%values (:,:,k) = Wbuf2%values (:,:,k) + Wbuf%values(:,:,kp) + Wbuf%values(:,:,km)
        end do
        ! Along Y
        IF (.NOT. setAlongY(spec_rank,Ubuf)) RETURN  
        IF (.NOT. setAlongY(spec_rank,Vbuf)) RETURN  
        IF (.NOT. setAlongY(spec_rank,Wbuf)) RETURN  
        IF (.NOT. setAlongY(spec_rank,Ubuf2)) RETURN  
        IF (.NOT. setAlongY(spec_rank,Vbuf2)) RETURN  
        IF (.NOT. setAlongY(spec_rank,Wbuf2)) RETURN  
        do j = Ubuf%ymin, Ubuf%ymax
           jm = j-1
           jp = j+1
           if (j.eq.1) jm = Ubuf%ny
           if (j.eq.Ubuf%ny) jp = 1
           Ubuf2%values (:,j,:) = Ubuf2%values (:,j,:) + Ubuf%values(:,jp,:) + Ubuf%values(:,km,:)
           Vbuf2%values (:,j,:) = Vbuf2%values (:,j,:) + Vbuf%values(:,jp,:) + Vbuf%values(:,km,:)
           Wbuf2%values (:,j,:) = Wbuf2%values (:,j,:) + Wbuf%values(:,jp,:) + Wbuf%values(:,km,:)
        end do
        
        Ubuf%values = Ubuf2%values
        Vbuf%values = Vbuf2%values
        Wbuf%values = Wbuf2%values
        ! Along X -> go back in the initial shape
        IF (.NOT. setAlongX(spec_rank,Ubuf)) RETURN  
        IF (.NOT. setAlongX(spec_rank,Vbuf)) RETURN  
        IF (.NOT. setAlongX(spec_rank,Wbuf)) RETURN  
        IF (.NOT. setAlongX(spec_rank,Ubuf2)) RETURN  
        IF (.NOT. setAlongX(spec_rank,Vbuf2)) RETURN  
        IF (.NOT. setAlongX(spec_rank,Wbuf2)) RETURN  
     end do
     call deleteDataLayout(Ubuf2)
     call deleteDataLayout(Vbuf2)
     call deleteDataLayout(Wbuf2)
  end if
  
  ! Structure function computation
  ! Along X  --> initial shape
  do i = Ubuf%xmin, Ubuf%xmax
     im = i-1
     ip = i+1
     if (i.eq.1) im = Ubuf%nx
     if (i.eq.Ubuf%nx) ip = 1
     SF%values(i,:,:) = ( Ubuf%values(i,:,:) - Ubuf%values(im,:,:) )**2.0 + ( Vbuf%values(i,:,:) - Vbuf%values(im,:,:) )**2.0 + ( Wbuf%values(i,:,:) - Wbuf%values(im,:,:) )**2.0 &
                    & + ( Ubuf%values(i,:,:) - Ubuf%values(ip,:,:) )**2.0 + ( Vbuf%values(i,:,:) - Vbuf%values(ip,:,:) )**2.0 + ( Wbuf%values(i,:,:) - Wbuf%values(ip,:,:) )**2.0 
  end do  
  ! Along Z
  IF (.NOT. setAlongZ(spec_rank,SF)) RETURN  
  IF (.NOT. setAlongZ(spec_rank,Ubuf)) RETURN  
  IF (.NOT. setAlongZ(spec_rank,Vbuf)) RETURN  
  IF (.NOT. setAlongZ(spec_rank,Wbuf)) RETURN  
  do k = Ubuf%zmin, Ubuf%zmax
     km = k-1
     kp = k+1
     if (k.eq.1) km = Ubuf%nz
     if (k.eq.Ubuf%nz) kp = 1
     SF%values(:,:,k) = ( Ubuf%values(:,:,k) - Ubuf%values(:,:,km) )**2.0 + ( Vbuf%values(:,:,k) - Vbuf%values(:,:,km) )**2.0 + ( Wbuf%values(:,:,k) - Wbuf%values(:,:,km) )**2.0 &
                    & + ( Ubuf%values(:,:,k) - Ubuf%values(:,:,kp) )**2.0 + ( Vbuf%values(:,:,k) - Vbuf%values(:,:,kp) )**2.0 + ( Wbuf%values(:,:,k) - Wbuf%values(:,:,kp) )**2.0 
  end do  
  ! Along Y
  IF (.NOT. setAlongY(spec_rank,SF)) RETURN  
  IF (.NOT. setAlongY(spec_rank,Ubuf)) RETURN  
  IF (.NOT. setAlongY(spec_rank,Vbuf)) RETURN  
  IF (.NOT. setAlongY(spec_rank,Wbuf)) RETURN  
  do j = Ubuf%ymin, Ubuf%ymax
     jm = j-1
     jp = j+1
     if (j.eq.1) jm = Ubuf%ny
     if (j.eq.Ubuf%ny) jp = 1
     SF%values(:,j,:) = ( Ubuf%values(:,j,:) - Ubuf%values(:,jm,:) )**2.0 + ( Vbuf%values(:,j,:) - Vbuf%values(:,jm,:) )**2.0 + ( Wbuf%values(:,j,:) - Wbuf%values(:,jm,:) )**2.0 &
                    & + ( Ubuf%values(:,j,:) - Ubuf%values(:,jp,:) )**2.0 + ( Vbuf%values(:,j,:) - Vbuf%values(:,jp,:) )**2.0 + ( Wbuf%values(:,j,:) - Wbuf%values(:,jp,:) )**2.0 
  end do  
  ! Along X  --> go back in initial shape and delete buffer
  IF (.NOT. setAlongX(spec_rank,SF)) RETURN 
  call deleteDataLayout(Ubuf)
  call deleteDataLayout(Vbuf)
  call deleteDataLayout(Wbuf)
  SF%values = 1.0_WP / 6.0_WP * SF%values
  
  if (.not. copyStructOnly(U,S11) .or. &
     &.not. copyStructOnly(U,S12) .or. &
     &.not. copyStructOnly(U,S13) .or. &
     &.not. copyStructOnly(U,S22) .or. &
     &.not. copyStructOnly(U,S23) .or. &
     &.not. copyStructOnly(U,S33) ) then
       res = .false.
       write(6,'(a,i0)')'[ERROR] in Structure Function model for vel : not enought memory !',spec_rank
  endif
  call computeStressTensor(Uk,VK,WK,VelWN,S11,S12,S13,S22,S23,S33,res)
  deltx = computeDelta(U)

  if (numVel.eq.4) then
     coeff = -2.0_WP*0.0014_WP*1.4_WP**(-3.0_WP/2.0_WP)
  else
     coeff = -2.0_WP*0.105_WP*1.4_WP**(-3.0_WP/2.0_WP)
  end if
  SF%values = (SF%values)**(0.5_WP)
  S11%values = coeff * SF%values * deltx * S11%values
  S12%values = coeff * SF%values * deltx * S12%values
  S13%values = coeff * SF%values * deltx * S13%values
  S22%values = coeff * SF%values * deltx * S22%values
  S23%values = coeff * SF%values * deltx * S23%values
  S33%values = coeff * SF%values * deltx * S33%values

! call cpu_time(tps1)
  if (.not. copyStructOnly(Uk,S11K) .or. &
     &.not. copyStructOnly(Uk,S12K) .or. &
     &.not. copyStructOnly(Uk,S13K) .or. &
     &.not. copyStructOnly(Uk,S22K) .or. &
     &.not. copyStructOnly(Uk,S23K) .or. &
     &.not. copyStructOnly(Uk,S33K) ) then
       res = .false.
       write(6,'(a,i0)')'[ERROR] in smagorinsky for vel : not enought memory !',spec_rank
  endif
! call cpu_time(tps2)
!print *,'temps mis pour allocation:',tps2-tps1

  call ftran(S11,S11K,res)
  call ftran(S12,S12K,res)
  call ftran(S13,S13K,res)
  call ftran(S22,S22K,res)
  call ftran(S23,S23K,res)
  call ftran(S33,S33K,res)

! call cpu_time(tps1)
  call deleteDataLayout(SF)
  call deleteDataLayout(S11)
  call deleteDataLayout(S12)
  call deleteDataLayout(S13)
  call deleteDataLayout(S22)
  call deleteDataLayout(S23)
  call deleteDataLayout(S33)
! call cpu_time(tps2)
!print *,'temps mis pour deallocation:',tps2-tps1

  call computeDivergenceK(VelWN,S11K,S12K,S13K,div1)
  call computeDivergenceK(VelWN,S12K,S22K,S23K,div2)
  call computeDivergenceK(VelWN,S13K,S23K,S33K,div3)

  call deleteDataLayout(S11K)
  call deleteDataLayout(S12K)
  call deleteDataLayout(S13K)
  call deleteDataLayout(S22K)
  call deleteDataLayout(S23K)
  call deleteDataLayout(S33K)

 end subroutine SfVel 


subroutine computeLijVec(Uin,Vin,Win,Uinf,Vinf,Winf&
                     &,Uk,L12,L13,L23,L11,L22,L33,delta,wave,filter)
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Uin,Vin,Win
!   TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: t12,t13,t23
  
  TYPE(REAL_DATA_LAYOUT), INTENT(INOUT) :: Uinf,Vinf,Winf
  TYPE(REAL_DATA_LAYOUT), INTENT(INOUT) :: L11,L22,L33,L12,L13,L23 
  TYPE(COMPLEX_DATA_LAYOUT) :: Uk

  TYPE(REAL_DATA_LAYOUT) :: tab,tabf
  TYPE(COMPLEX_DATA_LAYOUT) :: tabk,tabkf
  TYPE(WAVENUMBERS), intent(in) :: wave

  REAL(WP), INTENT(IN) :: delta
  logical :: res
  integer :: filter


  IF ( (.NOT. copyStructOnly(Uk,tabk)) .OR. &
      &(.NOT. copyStructOnly(Uk,tabkf)) .OR. &
      &(.NOT. copyStructOnly(Uin,tab)) .OR. &
      &(.NOT. copyStructOnly(Uin,tabf)) ) RETURN 

  tab%values = Uin%values*Uin%values  
     CALL ftran(tab,tabk,res)
     CALL computeFilter(wave,2.0*delta,tabk,tabkf,Filter) 
     CALL btran(tabkf,tabf,res)
  L11%values = tabf%values - (Uinf%values*Uinf%values)

  tab%values = Uin%values*Vin%values  
     CALL ftran(tab,tabk,res)
     CALL computeFilter(wave,2.0*delta,tabk,tabkf,Filter) 
     CALL btran(tabkf,tabf,res)
  L12%values = tabf%values - (Uinf%values*Vinf%values)
  
  tab%values = Uin%values*Win%values  
     CALL ftran(tab,tabk,res)
     CALL computeFilter(wave,2.0*delta,tabk,tabkf,Filter) 
     CALL btran(tabkf,tabf,res)
  L13%values = tabf%values - (Uinf%values*Winf%values)

  tab%values = Vin%values*Vin%values  
     CALL ftran(tab,tabk,res)
     CALL computeFilter(wave,2.0*delta,tabk,tabkf,Filter) 
     CALL btran(tabkf,tabf,res)
  L22%values = tabf%values - (Vinf%values*Vinf%values)

  tab%values = Vin%values*Win%values  
     CALL ftran(tab,tabk,res)
     CALL computeFilter(wave,2.0*delta,tabk,tabkf,Filter) 
     CALL btran(tabkf,tabf,res)
  L23%values = tabf%values - (Vinf%values*Winf%values)

  tab%values = Win%values*Win%values  
     CALL ftran(tab,tabk,res)
     CALL computeFilter(wave,2.0*delta,tabk,tabkf,Filter) 
     CALL btran(tabkf,tabf,res)
  L33%values = tabf%values - (Winf%values*Winf%values)

  call deleteDataLayout(tab)
  call deleteDataLayout(tabf)
  call deleteDataLayout(tabk)
  call deleteDataLayout(tabkf)

end subroutine computeLijVec


!------------------------------------------------------------------------------
!> @author 
!> Mouloud KESSAR, LEGI
!>
!> @details
!>  Computation of modeled Stress tensor for velocity field, with gradient model  
!> this subroutine can be use to compute TijU as well as TijB
!> @param [in/out]  t11,t12,t13,t22,t23,t33
!> @param [in] VecX,VecY,VecZ: componants of the Vector Field
!> @param [in] Wave: wavenumber
!> @param [in] delta: size of filter/ smallest scale resolved in LES
!------------------------------------------------------------------------------
subroutine gradientVelDIV(T1,T2,T3,VecX,VecY,VecZ,Wave,delta,spec_rank,nbcpus)

implicit none

  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: VecX,VecY,VecZ
  INTEGER :: spec_rank,nbcpus
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT) :: T1,T2,T3
  REAL(WP) :: delta
  logical :: res
  TYPE(WaveNumbers), INTENT(IN) :: Wave
  TYPE(REAL_DATA_LAYOUT) :: t11,t12,t13,t22,t23,t33


  IF ((.NOT. copyStructOnly(VecX,t11)) .OR. &
     &(.NOT. copyStructOnly(VecX,t12)) .OR. &
     &(.NOT. copyStructOnly(VecX,t13)) .OR. &
     &(.NOT. copyStructOnly(VecX,t22)) .OR. &
     &(.NOT. copyStructOnly(VecX,t23)) .OR. &
     &(.NOT. copyStructOnly(VecX,t33))    ) RETURN 

  call gradientVel(t11,t12,t13,t22,t23,t33,VecX,VecY,VecZ,Wave,delta,spec_rank,nbcpus)

  call ComputeDivSymTijVec1Vec2(T1,T2,T3,T12,T13,T23,T11,T22,T33,wave,res)

    CALL deleteDataLayout(t11)
    CALL deleteDataLayout(t12)
    CALL deleteDataLayout(t13)
    CALL deleteDataLayout(t22)
    CALL deleteDataLayout(t23)
    CALL deleteDataLayout(t33)

end subroutine gradientVelDIV


!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> Modèle du gradient regularisé 
!> @param [in] VecX velocity component along first direction in spectral space
!> @param [in] VecY velocity component along second direction in spectral space
!> @param [in] VecZ velocity component along third direction in spectral space
!> @param [in] FieldForParam field to get informations about grid (physical space) 
!> @param [in] Wave wave numbers
!> @param [in] spec_rank rank of processor
!> @param [inout] T1 = C delta^2 (Sjk^minus d^2u1/dxjdxk) in spectral space
!> @param [inout] T2 = C delta^2 (Sjk^minus d^2u2/dxjdxk) in spectral space
!> @param [inout] T3 = C delta^2 (Sjk^minus d^2u3/dxjdxk) in spectral space
!------------------------------------------------------------------------------
subroutine RGMDivDynV1(T1,T2,T3,FieldForParam,VecX,VecY,VecZ,Wave,spec_rank,filter,ite)
  
  use stat_tools , only : computeFieldPDF 
  use mpilayout_tools , only : getnbcpus 
 
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: VecX,VecY,VecZ
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT) :: T1,T2,T3
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: FieldForParam 
  TYPE(WaveNumbers), INTENT(IN)      :: Wave
  INTEGER,INTENT(IN)                 :: spec_rank,filter
  INTEGER,INTENT(IN),OPTIONAL        :: ite 

  TYPE(REAL_DATA_LAYOUT) :: M1,M2,M3,L1,L2,L3
  TYPE(REAL_DATA_LAYOUT) :: N1,N2,N3,T1aux,T2aux,T3aux
  REAL(WP) :: delta,Cc,den,num
  logical  :: res,pdffield=.false.

  if (present(ite)) then
    if ( ite .ne. -1 ) pdffield = .true.
  endif

  IF (.NOT. samelayout(VecX,VecY,VecZ)) print *,'[ERROR]'
  IF (.NOT. copyStructOnly(FieldForParam,M1) .OR. &
     &.NOT. copyStructOnly(FieldForParam,M2) .OR. &
     &.NOT. copyStructOnly(FieldForParam,M3) .OR. &
     &.NOT. copyStructOnly(FieldForParam,L1) .OR. &
     &.NOT. copyStructOnly(FieldForParam,L2) .OR. &
     &.NOT. copyStructOnly(FieldForParam,L3) .OR. &
     &.NOT. copyStructOnly(FieldForParam,N1) .OR. &
     &.NOT. copyStructOnly(FieldForParam,N2) .OR. &
     &.NOT. copyStructOnly(FieldForParam,N3) .OR. &
     &.NOT. copyStructOnly(FieldForParam,T1aux) .OR. &
     &.NOT. copyStructOnly(FieldForParam,T2aux) .OR. &
     &.NOT. copyStructOnly(FieldForParam,T3aux)) THEN
    print *,'[ERROR]'
  ENDIF
  
  delta=computeDelta(FieldForParam)

  call computeFilter(Wave,2.0_WP*delta,VecX,N1,filter)
  call computeFilter(Wave,2.0_WP*delta,VecY,N2,filter)
  call computeFilter(Wave,2.0_WP*delta,VecZ,N3,filter)

  !Coef Dyn = <Li Ni>/<Ni Ni>
  ! Ni = delta|^**2 * S|^moins_kj  d^2u|^_i/dxjdxk - delta|**2*(S|moins_kj d^2u|_i/dxjdxk )^
  ! Li = d/dxi(tau_ij)
  call computeT_ijVecA(VecX,VecY,VecZ,N1,N2,N3,Wave,T1aux,T2aux,T3aux,L1,L2,L3,filter,res,2.0_WP*delta)
  call ComputeDivSymTijVec1Vec2(T1,T2,T3,T2aux,T3aux,L2,T1aux,L1,L3,Wave,res)
  call btran(T1,L1,res)
  call btran(T2,L2,res)
  call btran(T3,L3,res)

  call RGMVel(T1,T2,T3,FieldForParam,N1,N2,N3,Wave,spec_rank,'moins',delta=1.0_WP,coef=1.0_WP)
  call btran(T1,N1,res)
  call btran(T2,N2,res)
  call btran(T3,N3,res)
  call RGMVel(T1,T2,T3,FieldForParam,VecX,VecY,VecZ,Wave,spec_rank,'moins',delta=1.0_WP,coef=1.0_WP)
  call btran(T1,M1,res)
  call btran(T2,M2,res)
  call btran(T3,M3,res)
  call computeFilter(Wave,2.0_WP*delta,M1,L1,filter)
  call computeFilter(Wave,2.0_WP*delta,M2,L2,filter)
  call computeFilter(Wave,2.0_WP*delta,M3,L3,filter)
  N1%values = (2.0_WP*delta)**2*N1%values - delta**2* L1%values
  N2%values = (2.0_WP*delta)**2*N2%values - delta**2* L2%values
  N3%values = (2.0_WP*delta)**2*N3%values - delta**2* L3%values

  L1%values = L1%values * N1%values + L2%values * N2%values + L3%values * N3%values
  num = computeFieldAvg(L1,spec_rank)
  L2%values = L1%values**2+ L2%values**2+ L3%values**2
  den = computeFieldAvg(L2,spec_rank)
  Cc = num / den
  if (spec_rank .eq. 0 ) write (6,'(a,1x,g15.8)') '[INFO RGM Vel V2] Coef dyn=',Cc
 
  L1%values = Cc * delta**2* N1%values
  L2%values = Cc * delta**2* N2%values
  L3%values = Cc * delta**2* N3%values

  if ( pdffield ) then
    L1%name="sgs_RGMDynV2_T1"
    L2%name="sgs_RGMDynV2_T2"
    L3%name="sgs_RGMDynV2_T3"
    if (.not. computeFieldPDF(ite,L1,300,spec_rank,10000)) return
    if (.not. computeFieldPDF(ite,L2,300,spec_rank,10000)) return
    if (.not. computeFieldPDF(ite,L3,300,spec_rank,10000)) return
  endif
  call ftran(L1,T1,res)
  call ftran(L2,T2,res)
  call ftran(L3,T3,res)

  CALL deleteDataLayout(L1)
  CALL deleteDataLayout(L2)
  CALL deleteDataLayout(L3)
  CALL deleteDataLayout(M1)
  CALL deleteDataLayout(M2)
  CALL deleteDataLayout(M3)
  CALL deleteDataLayout(N1)
  CALL deleteDataLayout(N2)
  CALL deleteDataLayout(N3)
  CALL deleteDataLayout(T1aux)
  CALL deleteDataLayout(T2aux)
  CALL deleteDataLayout(T3aux)


end subroutine RGMDivDynV1


!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> Modèle du gradient regularisé 
!> @param [in] VecX velocity component along first direction in spectral space
!> @param [in] VecY velocity component along second direction in spectral space
!> @param [in] VecZ velocity component along third direction in spectral space
!> @param [in] FieldForParam field to get informations about grid (physical space) 
!> @param [in] Wave wave numbers
!> @param [in] spec_rank rank of processor
!> @param [inout] T1 = C delta^2 (Sjk^minus d^2u1/dxjdxk) in spectral space
!> @param [inout] T2 = C delta^2 (Sjk^minus d^2u2/dxjdxk) in spectral space
!> @param [inout] T3 = C delta^2 (Sjk^minus d^2u3/dxjdxk) in spectral space
!------------------------------------------------------------------------------
subroutine RGMDivDynV2(T1,T2,T3,FieldForParam,VecX,VecY,VecZ,Wave,spec_rank,filter,ite)
  
  use stat_tools , only : computeFieldPDF 
  use mpilayout_tools , only : getnbcpus 
 
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: VecX,VecY,VecZ
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT) :: T1,T2,T3
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: FieldForParam 
  TYPE(WaveNumbers), INTENT(IN)      :: Wave
  INTEGER,INTENT(IN)                 :: spec_rank,filter
  INTEGER,INTENT(IN),OPTIONAL        :: ite

  TYPE(REAL_DATA_LAYOUT) :: M1,M2,M3,L1,L2,L3,VecXFil,VecYFil,VecZFil
  TYPE(REAL_DATA_LAYOUT) :: N1,N2,N3,T1aux,T2aux,T3aux
  REAL(WP)               :: delta,Cc,den,num
  logical                :: res,pdffield=.false.

  if (present(ite)) then
    if ( ite .ne. -1 ) pdffield = .true.
  endif

  IF (.NOT. samelayout(VecX,VecY,VecZ)) print *,'[ERROR]'
  IF (.NOT. copyStructOnly(FieldForParam,M1) .OR. &
     &.NOT. copyStructOnly(FieldForParam,M2) .OR. &
     &.NOT. copyStructOnly(FieldForParam,M3) .OR. &
     &.NOT. copyStructOnly(FieldForParam,L1) .OR. &
     &.NOT. copyStructOnly(FieldForParam,L2) .OR. &
     &.NOT. copyStructOnly(FieldForParam,L3) .OR. &
     &.NOT. copyStructOnly(FieldForParam,N1) .OR. &
     &.NOT. copyStructOnly(FieldForParam,N2) .OR. &
     &.NOT. copyStructOnly(FieldForParam,N3) .OR. &
     &.NOT. copyStructOnly(FieldForParam,VecXFil) .OR. &
     &.NOT. copyStructOnly(FieldForParam,VecYFil) .OR. &
     &.NOT. copyStructOnly(FieldForParam,VecZFil) .OR. &
     &.NOT. copyStructOnly(FieldForParam,T1aux) .OR. &
     &.NOT. copyStructOnly(FieldForParam,T2aux) .OR. &
     &.NOT. copyStructOnly(FieldForParam,T3aux)) THEN
    print *,'[ERROR]'
  ENDIF
  
  delta=computeDelta(FieldForParam)

  call computeFilter(Wave,2.0_WP*delta,VecX,VecXFil,filter)
  call computeFilter(Wave,2.0_WP*delta,VecY,VecYFil,filter)
  call computeFilter(Wave,2.0_WP*delta,VecZ,VecZFil,filter)

  !Coef Dyn = <(Li-Mi)Ni>/<Ni Ni>
  ! Ni = deltaFt^2 * SjkFt^plus d^2uiFt/dxjdxk  - ( delta^2 Skj^plus d^2ui/dxjdxk )Ft
  ! Mi = deltaFt^2/12*SkjFt^moins d^2uiFt/dxjdxk - ( delta^2/12*Skj^moins d^2ui/dxjdxk )Ft)
  ! Li = d/dxi(tau_ij)

  call RGMVel(T1,T2,T3,FieldForParam,VecXFil,VecYFil,VecZFil,Wave,spec_rank,'plus',delta=2.0_WP*delta,coef=1.0_WP)
  call btran(T1,N1,res)
  call btran(T2,N2,res)
  call btran(T3,N3,res)
  call RGMVel(T1,T2,T3,FieldForParam,VecX,VecY,VecZ,Wave,spec_rank,'plus',delta=1.0_WP,coef=1.0_WP)
  call btran(T1,T1aux,res)
  call btran(T2,T2aux,res)
  call btran(T3,T3aux,res)
  call computeFilter(Wave,2.0_WP*delta,T1aux,M1,filter)
  call computeFilter(Wave,2.0_WP*delta,T2aux,M2,filter)
  call computeFilter(Wave,2.0_WP*delta,T3aux,M3,filter)
  N1%values = N1%values - delta**2* M1%values
  N2%values = N2%values - delta**2* M2%values
  N3%values = N3%values - delta**2* M3%values

  call RGMVel(T1,T2,T3,FieldForParam,VecXFil,VecYFil,VecZFil,Wave,spec_rank,'moins',delta=2.0_WP*delta)
  call btran(T1,M1,res)
  call btran(T2,M2,res)
  call btran(T3,M3,res)
  call RGMVel(T1,T2,T3,FieldForParam,VecX,VecY,VecZ,Wave,spec_rank,'moins',delta=1.0_WP)
  call btran(T1,T1aux,res)
  call btran(T2,T2aux,res)
  call btran(T3,T3aux,res)
  call computeFilter(Wave,2.0_WP*delta,T1aux,L1,filter)
  call computeFilter(Wave,2.0_WP*delta,T2aux,L2,filter)
  call computeFilter(Wave,2.0_WP*delta,T3aux,L3,filter)
  M1%values = M1%values - delta**2* L1%values
  M2%values = M2%values - delta**2* L2%values
  M3%values = M3%values - delta**2* L3%values

  call computeT_ijVecA(VecX,VecY,VecZ,VecXFil,VecYFil,VecZFil,Wave,T1aux,T2aux,T3aux,L1,L2,L3,filter,res,2.0_WP*delta)
  call ComputeDivSymTijVec1Vec2(T1,T2,T3,T2aux,T3aux,L2,T1aux,L1,L3,Wave,res)
  call btran(T1,L1,res)
  call btran(T2,L2,res)
  call btran(T3,L3,res)
  
  L1%values = (L1%values - M1%values) * N1%values + (L2%values - M2%values) * N2%values + (L3%values - M3%values) * N3%values
  num = computeFieldAvg(L1,spec_rank)
  N1%values = N1%values ** 2.0_WP + N2%values ** 2.0_WP + N3%values ** 2.0_WP
  den = computeFieldAvg(N1,spec_rank)
  Cc = num / den
  if (spec_rank .eq. 0 ) write (6,'(a,1x,g15.8)') '[INFO RGM Vel V1] Coef dyn=',Cc
 
  call RGMVel(T1,T2,T3,FieldForParam,VecX,VecY,VecZ,Wave,spec_rank,'plus',delta=delta,coef=1.0_WP)
  call btran(T1,N1,res)
  call btran(T2,N2,res)
  call btran(T3,N3,res)
  call RGMVel(T1,T2,T3,FieldForParam,VecX,VecY,VecZ,Wave,spec_rank,'moins',delta=delta) 
  call btran(T1,M1,res)
  call btran(T2,M2,res)
  call btran(T3,M3,res)

  L1%values = M1%values - Cc * N1%values
  L2%values = M2%values - Cc * N2%values
  L3%values = M3%values - Cc * N3%values

  if ( pdffield ) then
    L1%name="sgs_RGMDynV1_T1"
    L2%name="sgs_RGMDynV1_T2"
    L3%name="sgs_RGMDynV1_T3"
    if (.not. computeFieldPDF(ite,L1,300,spec_rank,10000)) return
    if (.not. computeFieldPDF(ite,L2,300,spec_rank,10000)) return
    if (.not. computeFieldPDF(ite,L3,300,spec_rank,10000)) return
  endif
  call ftran(L1,T1,res)
  call ftran(L2,T2,res)
  call ftran(L3,T3,res)

  CALL deleteDataLayout(L1)
  CALL deleteDataLayout(L2)
  CALL deleteDataLayout(L3)
  CALL deleteDataLayout(M1)
  CALL deleteDataLayout(M2)
  CALL deleteDataLayout(M3)
  CALL deleteDataLayout(N1)
  CALL deleteDataLayout(N2)
  CALL deleteDataLayout(N3)
  CALL deleteDataLayout(T1aux)
  CALL deleteDataLayout(T2aux)
  CALL deleteDataLayout(T3aux)
  CALL deleteDataLayout(VecXFil)
  CALL deleteDataLayout(VecYFil)
  CALL deleteDataLayout(VecZFil)


end subroutine RGMDivDynV2


!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> Modèle du gradient regularisé 
!> @param [in] VecX velocity component along first direction in spectral space
!> @param [in] VecY velocity component along second direction in spectral space
!> @param [in] VecZ velocity component along third direction in spectral space
!> @param [in] FieldForParam field to get informations about grid (physical space) 
!> @param [in] Wave wave numbers
!> @param [in] spec_rank rank of processor
!> @param [inout] T1 = C delta^2 (Sjk^minus d^2u1/dxjdxk) in spectral space
!> @param [inout] T2 = C delta^2 (Sjk^minus d^2u2/dxjdxk) in spectral space
!> @param [inout] T3 = C delta^2 (Sjk^minus d^2u3/dxjdxk) in spectral space
!------------------------------------------------------------------------------
subroutine RGMVel(T1,T2,T3,FieldForParam,VecX,VecY,VecZ,Wave,spec_rank,signe,delta,coef,ite)
  
  use stat_tools , only : computeFieldPDF 
  use mpilayout_tools , only : getnbcpus 
 
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: VecX,VecY,VecZ
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT) :: T1,T2,T3
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: FieldForParam 
  TYPE(WaveNumbers), INTENT(IN)      :: Wave
  REAL(WP),INTENT(IN),OPTIONAL       :: coef,delta
  INTEGER,INTENT(IN)                 :: spec_rank
  CHARACTER(LEN=*),INTENT(IN)        :: signe
  INTEGER,INTENT(IN),OPTIONAL        :: ite

  TYPE(REAL_DATA_LAYOUT) :: S11,S12,S13,S22,S23,S33
  TYPE(REAL_DATA_LAYOUT) :: duidxdx,duidydy,duidzdz,duidxdy,duidxdz,duidydz 
  REAL(WP) :: deltaLoc,coefLoc
  logical :: res,pdffield=.false.

  if (present(ite)) then
    if ( ite .ne. -1 ) pdffield = .true.
  endif

  IF (.NOT. samelayout(VecX,VecY,VecZ)) print *,'[ERROR]'
  IF (.NOT. copyStructOnly(FieldForParam,duidxdx) .OR. &
     &.NOT. copyStructOnly(FieldForParam,duidydy) .OR. &
     &.NOT. copyStructOnly(FieldForParam,duidzdz) .OR. &
     &.NOT. copyStructOnly(FieldForParam,duidxdy) .OR. &
     &.NOT. copyStructOnly(FieldForParam,duidydz) .OR. &
     &.NOT. copyStructOnly(FieldForParam,duidxdz) .OR. &
     &.NOT. copyStructOnly(FieldForParam,S11) .OR. &
     &.NOT. copyStructOnly(FieldForParam,S22) .OR. &
     &.NOT. copyStructOnly(FieldForParam,S33) .OR. &
     &.NOT. copyStructOnly(FieldForParam,S12) .OR. &
     &.NOT. copyStructOnly(FieldForParam,S13) .OR. &
     &.NOT. copyStructOnly(FieldForParam,S23)) THEN
    print *,'[ERROR]'
  ENDIF
  
  if (present(coef)) then
    coefLoc=coef
  else
    coefLoc=1.0_WP/12.0_WP
  endif
  if (present(delta)) then
    deltaLoc=delta
  else
    deltaLoc=computeDelta(FieldForParam)
  endif 
  call ftran(VecX,T1,res)
  call ftran(VecY,T2,res)
  call ftran(VecZ,T3,res)
  !call computeStressTensor(VecX,VecY,VecZ,Wave,S11,S12,S13,S22,S23,S33,res)
  call computeStressTensor(T1,T2,T3,Wave       ,S11,S12,S13,S22,S23,S33,res)
  if (.not. computeSplitingSij(S11,S12,S13,S22,S23,S33,&
                              &duidxdx,duidydy,duidzdz,duidxdy,duidydz,duidxdz,&
                                    &signe)) print *,'[ERROR]' 
  S11%values = duidxdx%values
  S12%values = duidydy%values
  S13%values = duidzdz%values 
  S22%values = duidxdy%values 
  S23%values = duidydz%values 
  S33%values = duidxdz%values 
  !T1
  if (.not. computeDerivationField(VecX,(/1,1/),duidxdx,spec_rank,Wave)) print *,'[ERROR]'
  if (.not. computeDerivationField(VecX,(/2,2/),duidydy,spec_rank,Wave)) print *,'[ERROR]'
  if (.not. computeDerivationField(VecX,(/3,3/),duidzdz,spec_rank,Wave)) print *,'[ERROR]'
  if (.not. computeDerivationField(VecX,(/1,2/),duidxdy,spec_rank,Wave)) print *,'[ERROR]'
  if (.not. computeDerivationField(VecX,(/1,3/),duidxdz,spec_rank,Wave)) print *,'[ERROR]'
  if (.not. computeDerivationField(VecX,(/2,3/),duidydz,spec_rank,Wave)) print *,'[ERROR]'
  duidxdx%values= S11%values*duidxdx%values+ S22%values*duidydy%values+ S33%values*duidzdz%values+&
                 & 2.0_WP*(S12%values*duidxdy%values+ S13%values*duidxdz%values+ S23%values*duidydz%values)
  duidxdx%values= coefLoc*deltaLoc**2*duidxdx%values
  if ( pdffield ) then
    duidxdx%name="sgs_RGMVel_T1"
    S11%name="S11"
    S22%name="S22"
    duidydy%name="duidydy"
    S33%name="S33"
    duidzdz%name="duidzdz"
    S12%name="S12"
    duidxdy%name="duidxdy"
    S13%name="S13"
    duidxdz%name="duidxdz"
    S23%name="S23"
    duidydz%name="duidydz"
    if (.not. computeFieldPDF(ite,duidxdx,300,spec_rank,10000)) return
    if (.not. computeFieldPDF(ite,S11,300,spec_rank,10000)) return
    if (.not. computeFieldPDF(ite,S22,300,spec_rank,10000)) return
    if (.not. computeFieldPDF(ite,S33,300,spec_rank,10000)) return
    if (.not. computeFieldPDF(ite,S12,300,spec_rank,10000)) return
    if (.not. computeFieldPDF(ite,S13,300,spec_rank,10000)) return
    if (.not. computeFieldPDF(ite,S23,300,spec_rank,10000)) return
    if (.not. computeFieldPDF(ite,duidydy,300,spec_rank,10000)) return
    if (.not. computeFieldPDF(ite,duidzdz,300,spec_rank,10000)) return
    if (.not. computeFieldPDF(ite,duidxdy,300,spec_rank,10000)) return
    if (.not. computeFieldPDF(ite,duidxdz,300,spec_rank,10000)) return
    if (.not. computeFieldPDF(ite,duidydz,300,spec_rank,10000)) return
  endif
  CALL ftran(duidxdx,T1,res)
  if (.not. res) print *,'[ERROR]' 
  !T2
  if (.not. computeDerivationField(VecY,(/1,1/),duidxdx,spec_rank,Wave)) print *,'[ERROR]' 
  if (.not. computeDerivationField(VecY,(/2,2/),duidydy,spec_rank,Wave)) print *,'[ERROR]' 
  if (.not. computeDerivationField(VecY,(/3,3/),duidzdz,spec_rank,Wave)) print *,'[ERROR]'
  if (.not. computeDerivationField(VecY,(/1,2/),duidxdy,spec_rank,Wave)) print *,'[ERROR]'
  if (.not. computeDerivationField(VecY,(/1,3/),duidxdz,spec_rank,Wave)) print *,'[ERROR]'
  if (.not. computeDerivationField(VecY,(/2,3/),duidydz,spec_rank,Wave)) print *,'[ERROR]'
  duidxdx%values= S11%values*duidxdx%values+ S22%values*duidydy%values+ S33%values*duidzdz%values+&
                 & 2.0_WP*(S12%values*duidxdy%values+ S13%values*duidxdz%values+ S23%values*duidydz%values)
  duidxdx%values= coefLoc*deltaLoc**2*duidxdx%values
  if ( pdffield ) then
    duidxdx%name="sgs_RGMVel_T2"
    if (.not. computeFieldPDF(ite,duidxdx,300,spec_rank,10000)) return
  endif
  CALL ftran(duidxdx,T2,res)
  if (.not. res) print *,'[ERROR]' 
  !T3
  if (.not. computeDerivationField(VecZ,(/1,1/),duidxdx,spec_rank,Wave)) print *,'[ERROR]' 
  if (.not. computeDerivationField(VecZ,(/2,2/),duidydy,spec_rank,Wave)) print *,'[ERROR]' 
  if (.not. computeDerivationField(VecZ,(/3,3/),duidzdz,spec_rank,Wave)) print *,'[ERROR]'
  if (.not. computeDerivationField(VecZ,(/1,2/),duidxdy,spec_rank,Wave)) print *,'[ERROR]'
  if (.not. computeDerivationField(VecZ,(/1,3/),duidxdz,spec_rank,Wave)) print *,'[ERROR]'
  if (.not. computeDerivationField(VecZ,(/2,3/),duidydz,spec_rank,Wave)) print *,'[ERROR]'
  duidxdx%values= S11%values*duidxdx%values+ S22%values*duidydy%values+ S33%values*duidzdz%values+&
                 & 2.0_WP*(S12%values*duidxdy%values+ S13%values*duidxdz%values+ S23%values*duidydz%values)
  duidxdx%values= coefLoc*deltaLoc**2*duidxdx%values
  if ( pdffield ) then
    duidxdx%name="sgs_RGMVel_T3"
    if (.not. computeFieldPDF(ite,duidxdx,300,spec_rank,10000)) return
  endif
  CALL ftran(duidxdx,T3,res)
  if (.not. res) print *,'[ERROR]' 

  CALL deletedatalayout(duidxdx)
  CALL deletedatalayout(duidydy)
  CALL deletedatalayout(duidzdz)
  CALL deletedatalayout(duidxdy)
  CALL deletedatalayout(duidydz)
  CALL deletedatalayout(duidxdz)
  CALL deletedatalayout(S11)
  CALL deletedatalayout(S22)
  CALL deletedatalayout(S33)
  CALL deletedatalayout(S12)
  CALL deletedatalayout(S13)
  CALL deletedatalayout(S23)

end subroutine RGMVel

!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> @param [in] VecX  velocity component along first direction in spectral space
!> @param [in] VecXk velocity component along first direction in spectral space
!> @param [in] VecYk velocity component along second direction in spectral space
!> @param [in] VecZk velocity component along third direction in spectral space
!> @param [in] 
!> @param [in] Wave wave numbers
!> @param [in] spec_rank rank of processor
!> @param [inout] 
!> @param [inout] 
!> @param [inout] 
!------------------------------------------------------------------------------
subroutine RGMOikOjk_div(T1,T2,T3,VecX,VecXk,VecYk,VecZk,Wave,spec_rank,delta,coef,ite)

    !I/O data
    type(real_data_layout),intent(in)       :: VecX
    type(complex_data_layout),intent(in)    :: VecXk,VecYk,VecZk
    type(complex_data_layout),intent(inout) :: T1,T2,T3
    type(WaveNumbers),intent(in)            :: Wave
    integer,intent(in)                      :: spec_rank
    real(WP),intent(in),optional            :: coef,delta
    integer,intent(in),optional             :: ite

    !Local data
    type(real_data_layout)  :: T11,T12,T13,T22,T23,T33
    type(real_data_layout)  :: F11,F12,F13,F22,F23,F33
    logical                 :: pdffield,res
    real(WP)                :: coefLoc,deltaLoc

    pdffield= .false.
    if (present(ite)) then
      if ( ite .ne. -1 ) pdffield = .true.
    endif
    if (present(coef)) then
      coefLoc = coef
    else
      coefLoc = 1.0_WP/12.0_WP
    endif
    if (present(delta)) then
      deltaLoc = delta
    else
      deltaLoc = computeDelta(VecX)
    endif

    if (.not. copystructonly(VecX,T11) ) return
    if (.not. copystructonly(VecX,T12) ) return
    if (.not. copystructonly(VecX,T13) ) return
    if (.not. copystructonly(VecX,T22) ) return
    if (.not. copystructonly(VecX,T23) ) return
    if (.not. copystructonly(VecX,T33) ) return
    if (spec_rank .eq. 0) write (6,'(a)') '[INFO LES] sgs model for velocity: RGMSikmoinsSkj'

    call RGMOikOjk_tij(T11,T12,T13,T22,T23,T33,VecX,VecXk,VecYk,VecZk,Wave,deltaLoc,coefLoc)
    call ComputeDivSymTijVec1Vec2(T1,T2,T3,T12,T13,T23,T11,T22,T33,Wave,res)
    if (.not. res) write(6,'(a)') '[ERROR] in RGMSikmoinsSkj_div'
    if ( pdffield ) then
      if (.not. copystructonly(VecX,F11) ) return
      if (.not. copystructonly(VecX,F12) ) return
      if (.not. copystructonly(VecX,F13) ) return
      if (.not. copystructonly(VecX,F22) ) return
      if (.not. copystructonly(VecX,F23) ) return
      if (.not. copystructonly(VecX,F33) ) return
      call computeStressTensor(VecXk,VecYk,VecZk,Wave,F11,F12,F13,F22,F23,F33,res)
      T11%values = T11%values*F11%values + T22%values*F22%values + T33%values*F33%values + &
         2.0_WP * (T12%values*F12%values + T13%values*F13%values + T23%values*F23%values )
      T11%name="diss_tijSij"
      coefLoc = computeFieldAvg(T11,spec_rank)
      if (spec_rank .eq. 0) write (6,'(a,1x,g15.8)') '[INFO LES] RGM OikOjk dissipation over the box ',coefLoc
      call btran(T1,T11,res)
      call btran(T2,T22,res)
      call btran(T3,T33,res)
      T11%name="sgs_RGMO1kOjk_T1"
      T22%name="sgs_RGMO1kOjk_T2"
      T33%name="sgs_RGMO1kOjk_T3"
      if (.not. computeFieldPDF(ite,T11,300,spec_rank,10000)) return
      if (.not. computeFieldPDF(ite,T22,300,spec_rank,10000)) return
      if (.not. computeFieldPDF(ite,T33,300,spec_rank,10000)) return
      call deleteDataLayout(F11)
      call deleteDataLayout(F12)
      call deleteDataLayout(F13)
      call deleteDataLayout(F23)
      call deleteDataLayout(F22)
      call deleteDataLayout(F33)
    endif

    call deleteDataLayout(T11)
    call deleteDataLayout(T12)
    call deleteDataLayout(T13)
    call deleteDataLayout(T23)
    call deleteDataLayout(T22)
    call deleteDataLayout(T33)

end subroutine RGMOikOjk_div

!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> Modèle du gradient en prenant en compte l'incompressibilité. Calcul de 3 termes par composantes du
!> tenseur sous maille au lieu de 6. Calcul de d/dxj(tau_ij)=Ti
!> d/dxj(tau_ij) = C delta^2 (duj/dxk d^2ui/dxjdxk) 
!> @param [in] VecX velocity component along first direction in spectral space
!> @param [in] VecY velocity component along second direction in spectral space
!> @param [in] VecZ velocity component along third direction in spectral space
!> @param [in] FieldForParam field to get informations about grid (physical space) 
!> @param [in] Wave wave numbers
!> @param [in] spec_rank rank of processor
!> @param [inout] T1 = eps_1jk eps_jmn C delta^2 duk/dxl d^2un/dxmdxl
!> @param [inout] T2 = eps_2jk eps_jmn C delta^2 duk/dxl d^2un/dxmdxl
!> @param [inout] T3 = eps_3jk eps_jmn C delta^2 duk/dxl d^2un/dxmdxl
!------------------------------------------------------------------------------
subroutine RGMOikOjk_tij(T11,T12,T13,T22,T23,T33,VecX,VecXk,VecYk,VecZk,Wave,delta,coef)
  
    TYPE(REAL_DATA_LAYOUT), INTENT(IN)     :: VecX
    TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN)  :: VecXk,VecYk,VecZk
    TYPE(REAL_DATA_LAYOUT),intent(inout)   :: T11,T12,T13,T22,T23,T33
    TYPE(WaveNumbers), INTENT(IN)          :: Wave
    REAL(WP),INTENT(IN)                    :: delta,coef

    ! Local field
    type(REAL_DATA_LAYOUT) :: O23,O31,O12
    logical :: res

    ! ===== Initialisation =====
    res = .false.
    if (.not. copystructonly(VecX,O23) )return
    if (.not. copystructonly(VecX,O31) ) return
    if (.not. copystructonly(VecX,O12) ) return
    
    call computeVorticity(VecXk,VecYk,VecZk,Wave,O23,O31,O12,res)
    if (.not. res) return
    T11%values =   coef*delta**2*(O12%values**2 + O31%values**2)
    T12%values = - coef*delta**2*O31%values*O23%values
    T13%values = - coef*delta**2*O12%values*O23%values
    T22%values =   coef*delta**2*(O12%values**2 + O23%values**2)
    T23%values = - coef*delta**2*O12%values*O31%values
    T33%values =   coef*delta**2*(O31%values**2 + O23%values**2 )
    ! ===== Free memory and delete temporary fields =====
    call deleteDataLayout(O23)
    call deleteDataLayout(O31)
    call deleteDataLayout(O12)

    res = .true.

end subroutine RGMOikOjk_tij

!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> @param [in] VecX velocity component along first direction in spectral space
!> @param [in] VecY velocity component along second direction in spectral space
!> @param [in] VecZ velocity component along third direction in spectral space
!> @param [in] FieldForParam field to get informations about grid (physical space) 
!> @param [in] Wave wave numbers
!> @param [in] spec_rank rank of processor
!> @param [inout] T1 = Cd delta**2 d/dxj ( S1kmoins Skj)  + delta**2/12 d/dxj ( S1k Ojk + O1k Sjk )
!> @param [inout] T2 = Cd delta**2 d/dxj ( S2kmoins Skj)  + delta**2/12 d/dxj ( S2k Ojk + O2k Sjk )
!> @param [inout] T3 = Cd delta**2 d/dxj ( S3kmoins Skj)  + delta**2/12 d/dxj ( S3k Ojk + O3k Sjk )
!------------------------------------------------------------------------------
subroutine RGMSikmoinsSkjSikOjkOikSjkDyn2Fab(T1,T2,T3,VecX,VecY,VecZ,Wave,spec_rank,filter,delta,ite)

   use stat_tools , only : equalToVal
  
   TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: VecX,VecY,VecZ
   TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT) :: T1,T2,T3
   TYPE(WaveNumbers), INTENT(IN) :: Wave
   INTEGER,INTENT(IN)            :: spec_rank,filter
   REAL(WP),INTENT(IN),OPTIONAL  :: delta
   INTEGER,INTENT(IN),OPTIONAL   :: ite
   LOGICAL                       :: res 

    ! Local field
    type(REAL_DATA_LAYOUT) :: Sm11,Sm12,Sm13,Sm22,Sm23,Sm33
    type(REAL_DATA_LAYOUT) :: S11,S12,S13,S22,S23,S33
    type(REAL_DATA_LAYOUT) :: T11,T12,T13,T22,T23,T33
    type(REAL_DATA_LAYOUT) :: H11,H12,H13,H22,H23,H33
    type(REAL_DATA_LAYOUT) :: O23,O31,O12
    type(COMPLEX_DATA_LAYOUT)::Ufilk,Vfilk,Wfilk 
    integer                  :: nbcpus
!    real(WP) :: Num,Den
    real(WP) :: coefLoc,deltaLoc
    logical :: pdffield

    ! ===== Initialisation =====
    res = .false.
    nbcpus = getnbcpus()
    pdffield=.false.
    if (present(ite)) then
      if ( ite .ne. -1 ) pdffield = .true.
    endif
    if (present(delta)) then
      deltaLoc=delta
    else
      deltaLoc=computeDelta(VecX)
    endif 
    if (.not. initDataLayout("UFilk", UFilk,(VecX%nx/2)+1,VecX%ny,VecX%nz, &
       & VecX%Lx,VecX%Ly,VecX%Lz,nbcpus,spec_rank,alongZ)) then 
      write(6,'(a,i0,a)')'[ERROR] Not able to create Exk'
      return
    endif
    if (.not. copystructonly(VecX,T11) ) return
    if (.not. copystructonly(VecX,T12) ) return
    if (.not. copystructonly(VecX,T13) ) return
    if (.not. copystructonly(VecX,T22) ) return
    if (.not. copystructonly(VecX,T23) ) return
    if (.not. copystructonly(VecX,T33) ) return
    if (.not. copystructonly(VecX,H11) ) return
    if (.not. copystructonly(VecX,H12) ) return
    if (.not. copystructonly(VecX,H13) ) return
    if (.not. copystructonly(VecX,H22) ) return
    if (.not. copystructonly(VecX,H23) ) return
    if (.not. copystructonly(VecX,H33) ) return
    if (.not. copystructonly(VecX,Sm11) ) return
    if (.not. copystructonly(VecX,Sm12) ) return
    if (.not. copystructonly(VecX,Sm13) ) return
    if (.not. copystructonly(VecX,Sm22 ) )return
    if (.not. copystructonly(VecX,Sm23) ) return
    if (.not. copystructonly(VecX,Sm33) ) return
    if (.not. copystructonly(VecX,S11) ) return
    if (.not. copystructonly(VecX,S12) ) return
    if (.not. copystructonly(VecX,S13) ) return
    if (.not. copystructonly(VecX,S22) ) return
    if (.not. copystructonly(VecX,S23) ) return
    if (.not. copystructonly(VecX,S33) ) return
    if (.not. copystructonly(VecX,O23) )return
    if (.not. copystructonly(VecX,O31) ) return
    if (.not. copystructonly(VecX,O12) ) return
    if (.not. copystructonly(Ufilk,Vfilk) ) return
    if (.not. copystructonly(Ufilk,Wfilk) ) return
    

    !Calcul du coef dynamique
    !Cd = < (Sij-Tij)Hij>/<Hij Hij>
    !avec
    ! Sij=(u|_i u|_j)^ - u|^_i u|^_j
    ! Hij = delta|**2/12 S|^moinsik S|^jk
    ! Tij = delta|**2(S|^ik O|^jk + O|^ik S|^jk)
    
    !Calcul de Tij
    call computeFilter(Wave,2.0_WP*deltaLoc,VecX,Ufilk,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,VecY,Vfilk,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,VecZ,Wfilk,filter)
    call computeStressTensor(Ufilk,Vfilk,Wfilk,Wave,S11,S12,S13,S22,S23,S33,res)
    if (.not. computeSplitingSij(S11,S12,S13,S22,S23,S33,Sm11,Sm12,Sm13,Sm22,Sm23,Sm33,'moins')) return
    T11%values = Sm11%values*S11%values+Sm12%values*S12%values+Sm13%values*S13%values
    T12%values = Sm11%values*S12%values+Sm12%values*S22%values+Sm13%values*S23%values
    T13%values = Sm11%values*S13%values+Sm12%values*S23%values+Sm13%values*S33%values
    T22%values = Sm12%values*S12%values+Sm22%values*S22%values+Sm23%values*S23%values
    T23%values = Sm12%values*S13%values+Sm22%values*S23%values+Sm23%values*S33%values
    T33%values = Sm13%values*S13%values+Sm23%values*S23%values+Sm33%values*S33%values
    !Calcul de Hij
    call computeVorticity(Ufilk,Vfilk,Wfilk,Wave,O23,O31,O12,res)
    H11%values = deltaLoc**2 *(2.0_WP * ( S12%values*O12%values - S13%values*O31%values))
    H12%values = deltaLoc**2 *(S13%values*O23%values+S22%values*O12%values-S23%values*O31%values-S11%values*O12%values)
    H13%values = deltaLoc**2 *(S11%values*O31%values-S12%values*O23%values+O12%values*S23%values-O31%values*S33%values)
    H22%values = deltaLoc**2 *(2.0_WP * ( S23%values*O23%values - S12%values*O12%values))
    H23%values = deltaLoc**2 *(S12%values*O31%values-S22%values*O23%values-S13%values*O12%values+S33%values*O23%values)
    H33%values = deltaLoc**2 *(2.0_WP * ( S13%values*O31%values - S23%values*O23%values))
    !Calcul de Sij
    call computeFilter(Wave,2.0_WP*deltaLoc,VecX,Sm11,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,VecY,Sm22,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,VecZ,Sm33,filter)
    call computeT_ijVecA(VecX,VecY,VecZ,Sm11,Sm22,Sm33,Wave,S11,S12,S13,S22,S23,S33,filter,res,2.0_WP*deltaLoc)
    Sm11%values = (S11%values-T11%values)*H11%values + (S22%values-T22%values)*H22%values + (S33%values-T33%values)*H33%values + &
                & 2.0_WP*((S12%values-T12%values)*H12%values + (S13%values-T13%values)*H13%values +(S23%values-T23%values)*H23%values)
    Sm22%values = H11%values**2+H22%values**2+H33%values**2+2.0_WP*(H12%values**2+H13%values**2+H23%values**2) 
    coefLoc = computeFieldAvg(Sm11,spec_rank) / computeFieldAvg(Sm22,spec_rank) 
    if (spec_rank .eq. 0) write (6,'(a,1x,g15.8)') '[INFO LES RGMSikmoinsSkjSikOjkOikSjkDyn2 Fab] Dynanmic coef',coefLoc

    if (.not. computeSplitingSij(S11,S12,S13,S22,S23,S33,Sm11,Sm12,Sm13,Sm22,Sm23,Sm33,'moins')) return
    T11%values = coefLoc*deltaLoc**2*(Sm11%values*S11%values+Sm12%values*S12%values+Sm13%values*S13%values)+&
         deltaLoc**2/12.0_WP*(2.0_WP*(S12%values*O12%values-S13%values*O31%values))
    T12%values = coefLoc*deltaLoc**2*(Sm11%values*S12%values+Sm12%values*S22%values+Sm13%values*S23%values)-&
         deltaLoc**2/12.0_WP*(S11%values*O12%values+S13%values*O23%values+S22%values*O12%values-S23%values*O31%values )
    T13%values = coefLoc*deltaLoc**2*(Sm11%values*S13%values+Sm12%values*S23%values+Sm13%values*S33%values)+&
         deltaLoc**2/12.0_WP*(S11%values*O31%values-S12%values*O23%values+O12%values*S23%values-O31%values*S33%values )
    T22%values = coefLoc*deltaLoc**2*(Sm12%values*S12%values+Sm22%values*S22%values+Sm23%values*S23%values)+&
         deltaLoc**2/12.0_WP*(2.0_WP*(S23%values*O23%values-S12%values*O12%values))
    T23%values = coefLoc*deltaLoc**2*(Sm12%values*S13%values+Sm22%values*S23%values+Sm23%values*S33%values)+&
         deltaLoc**2/12.0_WP*(S12%values*O31%values-S22%values*O23%values-S13%values*O12%values+S33%values*O23%values )
    T33%values = coefLoc*deltaLoc**2*(Sm13%values*S13%values+Sm23%values*S23%values+Sm33%values*S33%values)+&
         deltaLoc**2/12.0_WP*(2.0_WP*(S13%values*O31%values-S23%values*O23%values))
    call ftran(T11,Ufilk,res)
    call ftran(T12,Vfilk,res)
    call ftran(T13,Wfilk,res)
    call computeDivergenceK(Wave,Ufilk,Vfilk,Wfilk,T1)
    call ftran(T22,Ufilk,res)
    call ftran(T23,Wfilk,res)
    call computeDivergenceK(Wave,UfilK,Vfilk,Wfilk,T2)
    call ftran(T13,Vfilk,res)
    call ftran(T33,Ufilk,res)
    call computeDivergenceK(Wave,Vfilk,Ufilk,Wfilk,T3)
    if ( pdffield ) then
      call btran(T1,T11,res)
      call btran(T2,T22,res)
      call btran(T3,T33,res)
      T11%name="sgs_RGMOikOjk_T1"
      T22%name="sgs_RGMOikOjk_T2"
      T33%name="sgs_RGMOikOjk_T3"
      if (.not. computeFieldPDF(ite,T11,300,spec_rank,10000)) return
      if (.not. computeFieldPDF(ite,T22,300,spec_rank,10000)) return
      if (.not. computeFieldPDF(ite,T33,300,spec_rank,10000)) return
    endif

    ! ===== Free memory and delete temporary fields =====
    call deleteDataLayout(H11)
    call deleteDataLayout(H12)
    call deleteDataLayout(H13)
    call deleteDataLayout(H23)
    call deleteDataLayout(H22)
    call deleteDataLayout(H33)
    call deleteDataLayout(T11)
    call deleteDataLayout(T12)
    call deleteDataLayout(T13)
    call deleteDataLayout(T23)
    call deleteDataLayout(T22)
    call deleteDataLayout(T33)
    call deleteDataLayout(Sm11)
    call deleteDataLayout(Sm12)
    call deleteDataLayout(Sm13)
    call deleteDataLayout(Sm23)
    call deleteDataLayout(Sm22)
    call deleteDataLayout(Sm33)
    call deleteDataLayout(S11)
    call deleteDataLayout(S12)
    call deleteDataLayout(S13)
    call deleteDataLayout(S23)
    call deleteDataLayout(S22)
    call deleteDataLayout(S33)
    call deleteDataLayout(O23)
    call deleteDataLayout(O31)
    call deleteDataLayout(O12)
    call deleteDataLayout(Ufilk)
    call deleteDataLayout(Vfilk)
    call deleteDataLayout(Wfilk)

    res = .true.


end subroutine RGMSikmoinsSkjSikOjkOikSjkDyn2Fab

!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> @param [in] VecX velocity component along first direction in spectral space
!> @param [in] VecY velocity component along second direction in spectral space
!> @param [in] VecZ velocity component along third direction in spectral space
!> @param [in] FieldForParam field to get informations about grid (physical space) 
!> @param [in] Wave wave numbers
!> @param [in] spec_rank rank of processor
!> @param [inout] T1 = Cd delta**2 d/dxj ( S1kmoins Skj)  + delta**2/12 d/dxj ( S1k Ojk + O1k Sjk )
!> @param [inout] T2 = Cd delta**2 d/dxj ( S2kmoins Skj)  + delta**2/12 d/dxj ( S2k Ojk + O2k Sjk )
!> @param [inout] T3 = Cd delta**2 d/dxj ( S3kmoins Skj)  + delta**2/12 d/dxj ( S3k Ojk + O3k Sjk )
!------------------------------------------------------------------------------
subroutine RGMSikmoinsSkjSikOjkOikSjkDyn2(T1,T2,T3,VecX,VecY,VecZ,Wave,spec_rank,filter,delta,ite)

   use stat_tools , only : equalToVal
  
   TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: VecX,VecY,VecZ
   TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT) :: T1,T2,T3
   TYPE(WaveNumbers), INTENT(IN) :: Wave
   INTEGER,INTENT(IN)            :: spec_rank,filter
   REAL(WP),INTENT(IN),OPTIONAL  :: delta
   INTEGER,INTENT(IN),OPTIONAL   :: ite
   LOGICAL                       :: res 

    ! Local field
    type(REAL_DATA_LAYOUT) :: Sm11,Sm12,Sm13,Sm22,Sm23,Sm33
    type(REAL_DATA_LAYOUT) :: S11,S12,S13,S22,S23,S33
    type(REAL_DATA_LAYOUT) :: T11,T12,T13,T22,T23,T33
    type(REAL_DATA_LAYOUT) :: L11,L12,L13,L22,L23,L33
    type(REAL_DATA_LAYOUT) :: H11,H12,H13,H22,H23,H33
    type(REAL_DATA_LAYOUT) :: O23,O31,O12
    type(COMPLEX_DATA_LAYOUT)::Ufilk,Vfilk,Wfilk 
    integer                  :: nbcpus
!    real(WP) :: Num,Den
    real(WP) :: coefLoc,deltaLoc
    logical :: pdffield

    ! ===== Initialisation =====
    res = .false.
    nbcpus = getnbcpus()
    pdffield=.false.
    if (present(ite)) then
      if ( ite .ne. -1 ) pdffield = .true.
    endif
    if (present(delta)) then
      deltaLoc=delta
    else
      deltaLoc=computeDelta(VecX)
    endif 
    if (.not. initDataLayout("UFilk", UFilk,(VecX%nx/2)+1,VecX%ny,VecX%nz, &
       & VecX%Lx,VecX%Ly,VecX%Lz,nbcpus,spec_rank,alongZ)) then 
      write(6,'(a,i0,a)')'[ERROR] Not able to create Exk'
      return
    endif
    if (.not. copystructonly(VecX,T11) ) return
    if (.not. copystructonly(VecX,T12) ) return
    if (.not. copystructonly(VecX,T13) ) return
    if (.not. copystructonly(VecX,T22) ) return
    if (.not. copystructonly(VecX,T23) ) return
    if (.not. copystructonly(VecX,T33) ) return
    if (.not. copystructonly(VecX,L11) ) return
    if (.not. copystructonly(VecX,L12) ) return
    if (.not. copystructonly(VecX,L13) ) return
    if (.not. copystructonly(VecX,L22) ) return
    if (.not. copystructonly(VecX,L23) ) return
    if (.not. copystructonly(VecX,L33) ) return
    if (.not. copystructonly(VecX,H11) ) return
    if (.not. copystructonly(VecX,H12) ) return
    if (.not. copystructonly(VecX,H13) ) return
    if (.not. copystructonly(VecX,H22) ) return
    if (.not. copystructonly(VecX,H23) ) return
    if (.not. copystructonly(VecX,H33) ) return
    if (.not. copystructonly(VecX,Sm11) ) return
    if (.not. copystructonly(VecX,Sm12) ) return
    if (.not. copystructonly(VecX,Sm13) ) return
    if (.not. copystructonly(VecX,Sm22 ) )return
    if (.not. copystructonly(VecX,Sm23) ) return
    if (.not. copystructonly(VecX,Sm33) ) return
    if (.not. copystructonly(VecX,S11) ) return
    if (.not. copystructonly(VecX,S12) ) return
    if (.not. copystructonly(VecX,S13) ) return
    if (.not. copystructonly(VecX,S22) ) return
    if (.not. copystructonly(VecX,S23) ) return
    if (.not. copystructonly(VecX,S33) ) return
    if (.not. copystructonly(VecX,O23) )return
    if (.not. copystructonly(VecX,O31) ) return
    if (.not. copystructonly(VecX,O12) ) return
    if (.not. copystructonly(Ufilk,Vfilk) ) return
    if (.not. copystructonly(Ufilk,Wfilk) ) return
    

    !Calcul du coef dynamique
    !Cd = < (Lij-Hij)Tij>/<Tij Tij>
    !avec
    ! Lij=(u|_i u|_j)^ - u|^_i u|^_j
    ! Hij = delta|^**2/12 S|^moinsik S|^jk - delta|**2/12 (S|moinsik S|kj)^
    ! Tij = delta|^**2(S|^ik O|^jk + O|^ik S|^jk) -delta|**2(S|ik O|jk + O|ik S|jk)^
    
    !Calcul de Hij
    call ftran(VecX,Ufilk,res)
    call ftran(VecY,Vfilk,res)
    call ftran(VecZ,Wfilk,res)
    call computeStressTensor(Ufilk,Vfilk,Wfilk,Wave,S11,S12,S13,S22,S23,S33,res)
    if (.not. computeSplitingSij(S11,S12,S13,S22,S23,S33,Sm11,Sm12,Sm13,Sm22,Sm23,Sm33,'moins')) return
    T11%values = Sm11%values*S11%values+Sm12%values*S12%values+Sm13%values*S13%values
    T12%values = Sm11%values*S12%values+Sm12%values*S22%values+Sm13%values*S23%values
    T13%values = Sm11%values*S13%values+Sm12%values*S23%values+Sm13%values*S33%values
    T22%values = Sm12%values*S12%values+Sm22%values*S22%values+Sm23%values*S23%values
    T23%values = Sm12%values*S13%values+Sm22%values*S23%values+Sm23%values*S33%values
    T33%values = Sm13%values*S13%values+Sm23%values*S23%values+Sm33%values*S33%values
    call computeFilter(Wave,2.0_WP*deltaLoc,T11,H11,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,T12,H12,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,T13,H13,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,T22,H22,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,T23,H23,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,T33,H33,filter)

    call computeFilter(Wave,2.0_WP*deltaLoc,VecX,Ufilk,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,VecY,Vfilk,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,VecZ,Wfilk,filter)
    call computeStressTensor(Ufilk,Vfilk,Wfilk,Wave,S11,S12,S13,S22,S23,S33,res)
    if (.not. computeSplitingSij(S11,S12,S13,S22,S23,S33,Sm11,Sm12,Sm13,Sm22,Sm23,Sm33,'moins')) return
    T11%values = Sm11%values*S11%values+Sm12%values*S12%values+Sm13%values*S13%values
    T12%values = Sm11%values*S12%values+Sm12%values*S22%values+Sm13%values*S23%values
    T13%values = Sm11%values*S13%values+Sm12%values*S23%values+Sm13%values*S33%values
    T22%values = Sm12%values*S12%values+Sm22%values*S22%values+Sm23%values*S23%values
    T23%values = Sm12%values*S13%values+Sm22%values*S23%values+Sm23%values*S33%values
    T33%values = Sm13%values*S13%values+Sm23%values*S23%values+Sm33%values*S33%values
    H11%values = (2.0_WP*deltaLoc)**2* T11%values - deltaLoc**2* H11%values
    H12%values = (2.0_WP*deltaLoc)**2* T12%values - deltaLoc**2* H12%values
    H13%values = (2.0_WP*deltaLoc)**2* T13%values - deltaLoc**2* H13%values
    H22%values = (2.0_WP*deltaLoc)**2* T22%values - deltaLoc**2* H22%values
    H23%values = (2.0_WP*deltaLoc)**2* T23%values - deltaLoc**2* H23%values
    H33%values = (2.0_WP*deltaLoc)**2* T33%values - deltaLoc**2* H33%values

    !Calcul de Mij
    ! Tij = delta|^**2(S|^ik O|^jk + O|^ik S|^jk) - delta|**2(S|ik O|jk + O|ik S|jk)^
    call computeVorticity(Ufilk,Vfilk,Wfilk,Wave,O23,O31,O12,res)
    T11%values = 2.0_WP * ( S12%values*O12%values - S13%values*O31%values) 
    T12%values = S13%values*O23%values+S22%values*O12%values-S23%values*O31%values-S11%values*O12%values
    T13%values = S11%values*O31%values-S12%values*O23%values+O12%values*S23%values-O31%values*S33%values
    T22%values = 2.0_WP * ( S23%values*O23%values - S12%values*O12%values)
    T23%values = S12%values*O31%values-S22%values*O23%values-S13%values*O12%values+S33%values*O23%values
    T33%values = 2.0_WP * ( S13%values*O31%values - S23%values*O23%values)
    call ftran(VecX,Ufilk,res)
    call ftran(VecY,Vfilk,res)
    call ftran(VecZ,Wfilk,res)
    call computeStressTensor(Ufilk,Vfilk,Wfilk,Wave,S11,S12,S13,S22,S23,S33,res)
    call computeVorticity(Ufilk,Vfilk,Wfilk,Wave,O23,O31,O12,res)
    Sm11%values = 2.0_WP * ( S12%values*O12%values - S13%values*O31%values) 
    Sm12%values = S13%values*O23%values+S22%values*O12%values-S23%values*O31%values-S11%values*O12%values
    Sm13%values = S11%values*O31%values-S12%values*O23%values+O12%values*S23%values-O31%values*S33%values
    Sm22%values = 2.0_WP * ( S23%values*O23%values - S12%values*O12%values)
    Sm23%values = S12%values*O31%values-S22%values*O23%values-S13%values*O12%values+S33%values*O23%values
    Sm33%values = 2.0_WP * ( S13%values*O31%values - S23%values*O23%values)
    call computeFilter(Wave,2.0_WP*deltaLoc,Sm11,L11,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,Sm12,L12,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,Sm13,L13,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,Sm22,L22,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,Sm23,L23,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,Sm33,L33,filter)
    T11%values = (2.0_WP*deltaLoc)**2 * T11%values - deltaLoc**2 * L11%values
    T12%values = (2.0_WP*deltaLoc)**2 * T12%values - deltaLoc**2 * L12%values
    T13%values = (2.0_WP*deltaLoc)**2 * T13%values - deltaLoc**2 * L13%values
    T22%values = (2.0_WP*deltaLoc)**2 * T22%values - deltaLoc**2 * L22%values
    T23%values = (2.0_WP*deltaLoc)**2 * T23%values - deltaLoc**2 * L23%values
    T33%values = (2.0_WP*deltaLoc)**2 * T33%values - deltaLoc**2 * L33%values

    !Calcul de Lij
    call computeFilter(Wave,2.0_WP*deltaLoc,VecX,Sm11,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,VecY,Sm22,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,VecZ,Sm33,filter)
!        computeT_ijVecA(in1 ,in2 ,in3 ,in1f,in2f,in3f,wave,out11,out12,out13,out22,out23,out33,filter,success,delta)
    call computeT_ijVecA(VecX,VecY,VecZ,Sm11,Sm22,Sm33,Wave,L11  ,L12  ,L13  ,L22  ,L23  ,L33  ,filter,res    ,2.0_WP*deltaLoc)

    Sm11%values = (L11%values-T11%values)*H11%values + (L22%values-T22%values)*H22%values + (L33%values-T33%values)*H33%values + &
                & 2.0_WP*((L12%values-T12%values)*H12%values + (L13%values-T13%values)*H13%values +(L23%values-T23%values)*H23%values)
    Sm22%values = H11%values**2+H22%values**2+H33%values**2+2.0_WP*(H12%values**2+H13%values**2+H23%values**2) 
    coefLoc = computeFieldAvg(Sm11,spec_rank) / computeFieldAvg(Sm22,spec_rank) 
    if (spec_rank .eq. 0) write (6,'(a,1x,g15.8)') '[INFO LES RGMSikmoinsSkjSikOjkOikSjkDyn2] Dynanmic coef',coefLoc

!                                         !out11,out12,out13,out22,out23,out33
!    call computeStressTensor(Ufilk,Vfilk,Wfilk,Wave,S11,S12,S13,S22,S23,S33,res)
    if (.not. computeSplitingSij(S11,S12,S13,S22,S23,S33,Sm11,Sm12,Sm13,Sm22,Sm23,Sm33,'moins')) return
!    call computeVorticity(Ufilk,Vfilk,Wfilk,Wave,O23,O31,O12,res)
    T11%values = coefLoc*deltaLoc**2*(Sm11%values*S11%values+Sm12%values*S12%values+Sm13%values*S13%values)+&
         deltaLoc**2/12.0_WP*(2.0_WP*(S12%values*O12%values-S13%values*O31%values))
    T12%values = coefLoc*deltaLoc**2*(Sm11%values*S12%values+Sm12%values*S22%values+Sm13%values*S23%values)-&
         deltaLoc**2/12.0_WP*(S11%values*O12%values+S13%values*O23%values+S22%values*O12%values-S23%values*O31%values )
    T13%values = coefLoc*deltaLoc**2*(Sm11%values*S13%values+Sm12%values*S23%values+Sm13%values*S33%values)+&
         deltaLoc**2/12.0_WP*(S11%values*O31%values-S12%values*O23%values+O12%values*S23%values-O31%values*S33%values )
    T22%values = coefLoc*deltaLoc**2*(Sm12%values*S12%values+Sm22%values*S22%values+Sm23%values*S23%values)+&
         deltaLoc**2/12.0_WP*(2.0_WP*(S23%values*O23%values-S12%values*O12%values))
    T23%values = coefLoc*deltaLoc**2*(Sm12%values*S13%values+Sm22%values*S23%values+Sm23%values*S33%values)+&
         deltaLoc**2/12.0_WP*(S12%values*O31%values-S22%values*O23%values-S13%values*O12%values+S33%values*O23%values )
    T33%values = coefLoc*deltaLoc**2*(Sm13%values*S13%values+Sm23%values*S23%values+Sm33%values*S33%values)+&
         deltaLoc**2/12.0_WP*(2.0_WP*(S13%values*O31%values-S23%values*O23%values))
    call ftran(T11,Ufilk,res)
    call ftran(T12,Vfilk,res)
    call ftran(T13,Wfilk,res)
    call computeDivergenceK(Wave,Ufilk,Vfilk,Wfilk,T1)
    call ftran(T22,Ufilk,res)
    call ftran(T23,Wfilk,res)
    call computeDivergenceK(Wave,UfilK,Vfilk,Wfilk,T2)
    call ftran(T13,Vfilk,res)
    call ftran(T33,Ufilk,res)
    call computeDivergenceK(Wave,Vfilk,Ufilk,Wfilk,T3)
    if ( pdffield ) then
      call btran(T1,T11,res)
      call btran(T2,T22,res)
      call btran(T3,T33,res)
      T11%name="sgs_RGMOikOjk_T1"
      T22%name="sgs_RGMOikOjk_T2"
      T33%name="sgs_RGMOikOjk_T3"
      if (.not. computeFieldPDF(ite,T11,300,spec_rank,10000)) return
      if (.not. computeFieldPDF(ite,T22,300,spec_rank,10000)) return
      if (.not. computeFieldPDF(ite,T33,300,spec_rank,10000)) return
    endif

    ! ===== Free memory and delete temporary fields =====
    call deleteDataLayout(L11)
    call deleteDataLayout(L12)
    call deleteDataLayout(L13)
    call deleteDataLayout(L23)
    call deleteDataLayout(L22)
    call deleteDataLayout(L33)
    call deleteDataLayout(H11)
    call deleteDataLayout(H12)
    call deleteDataLayout(H13)
    call deleteDataLayout(H23)
    call deleteDataLayout(H22)
    call deleteDataLayout(H33)
    call deleteDataLayout(T11)
    call deleteDataLayout(T12)
    call deleteDataLayout(T13)
    call deleteDataLayout(T23)
    call deleteDataLayout(T22)
    call deleteDataLayout(T33)
    call deleteDataLayout(Sm11)
    call deleteDataLayout(Sm12)
    call deleteDataLayout(Sm13)
    call deleteDataLayout(Sm23)
    call deleteDataLayout(Sm22)
    call deleteDataLayout(Sm33)
    call deleteDataLayout(S11)
    call deleteDataLayout(S12)
    call deleteDataLayout(S13)
    call deleteDataLayout(S23)
    call deleteDataLayout(S22)
    call deleteDataLayout(S33)
    call deleteDataLayout(O23)
    call deleteDataLayout(O31)
    call deleteDataLayout(O12)
    call deleteDataLayout(Ufilk)
    call deleteDataLayout(Vfilk)
    call deleteDataLayout(Wfilk)

    res = .true.


end subroutine RGMSikmoinsSkjSikOjkOikSjkDyn2



!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> @param [in] VecX velocity component along first direction in spectral space
!> @param [in] VecY velocity component along second direction in spectral space
!> @param [in] VecZ velocity component along third direction in spectral space
!> @param [in] FieldForParam field to get informations about grid (physical space) 
!> @param [in] Wave wave numbers
!> @param [in] spec_rank rank of processor
!> @param [inout] T1 =  delta**2/12 d/dxj ( S1kmoins Skj ) + Cd delta**2 d/dxj ( S1k Ojk + O1k Sjk )
!> @param [inout] T2 =  delta**2/12 d/dxj ( S2kmoins Skj ) + Cd delta**2 d/dxj ( S2k Ojk + O2k Sjk )
!> @param [inout] T3 =  delta**2/12 d/dxj ( S3kmoins Skj ) + Cd delta**2 d/dxj ( S3k Ojk + O3k Sjk )
!------------------------------------------------------------------------------
subroutine RGMSikmoinsSkjSikOjkOikSjkDynFab(T1,T2,T3,VecX,VecY,VecZ,Wave,spec_rank,filter,delta,ite)

   use stat_tools , only : equalToVal
  
   TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: VecX,VecY,VecZ
   TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT) :: T1,T2,T3
   TYPE(WaveNumbers), INTENT(IN) :: Wave
   INTEGER,INTENT(IN)            :: spec_rank,filter
   REAL(WP),INTENT(IN),OPTIONAL  :: delta
   INTEGER,INTENT(IN),OPTIONAL   :: ite
   LOGICAL                       :: res 

    ! Local field
    type(REAL_DATA_LAYOUT) :: Sm11,Sm12,Sm13,Sm22,Sm23,Sm33
    type(REAL_DATA_LAYOUT) :: S11,S12,S13,S22,S23,S33
    type(REAL_DATA_LAYOUT) :: T11,T12,T13,T22,T23,T33
    type(REAL_DATA_LAYOUT) :: H11,H12,H13,H22,H23,H33
    type(REAL_DATA_LAYOUT) :: O23,O31,O12
    type(COMPLEX_DATA_LAYOUT)::Ufilk,Vfilk,Wfilk 
    integer                  :: nbcpus
!    real(WP) :: Num,Den
    real(WP) :: coefLoc,deltaLoc
    logical :: pdffield

    ! ===== Initialisation =====
    res = .false.
    nbcpus = getnbcpus()
    pdffield=.false.
    if (present(ite)) then
      if ( ite .ne. -1 ) pdffield = .true.
    endif
    if (present(delta)) then
      deltaLoc=delta
    else
      deltaLoc=computeDelta(VecX)
    endif 
    if (.not. initDataLayout("UFilk", UFilk,(VecX%nx/2)+1,VecX%ny,VecX%nz, &
       & VecX%Lx,VecX%Ly,VecX%Lz,nbcpus,spec_rank,alongZ)) then 
      write(6,'(a,i0,a)')'[ERROR] Not able to create Exk'
      return
    endif
    if (.not. copystructonly(VecX,T11) ) return
    if (.not. copystructonly(VecX,T12) ) return
    if (.not. copystructonly(VecX,T13) ) return
    if (.not. copystructonly(VecX,T22) ) return
    if (.not. copystructonly(VecX,T23) ) return
    if (.not. copystructonly(VecX,T33) ) return
    if (.not. copystructonly(VecX,H11) ) return
    if (.not. copystructonly(VecX,H12) ) return
    if (.not. copystructonly(VecX,H13) ) return
    if (.not. copystructonly(VecX,H22) ) return
    if (.not. copystructonly(VecX,H23) ) return
    if (.not. copystructonly(VecX,H33) ) return
    if (.not. copystructonly(VecX,Sm11) ) return
    if (.not. copystructonly(VecX,Sm12) ) return
    if (.not. copystructonly(VecX,Sm13) ) return
    if (.not. copystructonly(VecX,Sm22 ) )return
    if (.not. copystructonly(VecX,Sm23) ) return
    if (.not. copystructonly(VecX,Sm33) ) return
    if (.not. copystructonly(VecX,S11) ) return
    if (.not. copystructonly(VecX,S12) ) return
    if (.not. copystructonly(VecX,S13) ) return
    if (.not. copystructonly(VecX,S22) ) return
    if (.not. copystructonly(VecX,S23) ) return
    if (.not. copystructonly(VecX,S33) ) return
    if (.not. copystructonly(VecX,O23) )return
    if (.not. copystructonly(VecX,O31) ) return
    if (.not. copystructonly(VecX,O12) ) return
    if (.not. copystructonly(Ufilk,Vfilk) ) return
    if (.not. copystructonly(Ufilk,Wfilk) ) return
    

    !Calcul du coef dynamique
    !Cd = < (Lij-Hij)Tij>/<Tij Tij>
    !avec
    ! Lij=(u|_i u|_j)^ - u|^_i u|^_j
    ! Hij = delta|**2/12 S|^moinsik S|^jk
    ! Tij = delta|**2(S|^ik O|^jk + O|^ik S|^jk)
    
    !Calcul de Hij
    call computeFilter(Wave,2.0_WP*deltaLoc,VecX,Ufilk,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,VecY,Vfilk,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,VecZ,Wfilk,filter)
    call computeStressTensor(Ufilk,Vfilk,Wfilk,Wave,S11,S12,S13,S22,S23,S33,res)
    if (.not. computeSplitingSij(S11,S12,S13,S22,S23,S33,Sm11,Sm12,Sm13,Sm22,Sm23,Sm33,'moins')) return
    H11%values = deltaLoc**2/12.0_WP *(Sm11%values*S11%values+Sm12%values*S12%values+Sm13%values*S13%values)
    H12%values = deltaLoc**2/12.0_WP *(Sm11%values*S12%values+Sm12%values*S22%values+Sm13%values*S23%values)
    H13%values = deltaLoc**2/12.0_WP *(Sm11%values*S13%values+Sm12%values*S23%values+Sm13%values*S33%values)
    H22%values = deltaLoc**2/12.0_WP *(Sm12%values*S12%values+Sm22%values*S22%values+Sm23%values*S23%values)
    H23%values = deltaLoc**2/12.0_WP *(Sm12%values*S13%values+Sm22%values*S23%values+Sm23%values*S33%values)
    H33%values = deltaLoc**2/12.0_WP *(Sm13%values*S13%values+Sm23%values*S23%values+Sm33%values*S33%values)

    !Calcul de Mij
    ! Tij = delta|^**2(S|^ik O|^jk + O|^ik S|^jk) - delta|**2(S|ik O|jk + O|ik S|jk)^
    call computeVorticity(Ufilk,Vfilk,Wfilk,Wave,O23,O31,O12,res)
    T11%values = deltaLoc**2 * (2.0_WP * ( S12%values*O12%values - S13%values*O31%values))
    T12%values = deltaLoc**2 * (S13%values*O23%values+S22%values*O12%values-S23%values*O31%values-S11%values*O12%values)
    T13%values = deltaLoc**2 * (S11%values*O31%values-S12%values*O23%values+O12%values*S23%values-O31%values*S33%values)
    T22%values = deltaLoc**2 * (2.0_WP * ( S23%values*O23%values - S12%values*O12%values))
    T23%values = deltaLoc**2 * (S12%values*O31%values-S22%values*O23%values-S13%values*O12%values+S33%values*O23%values)
    T33%values = deltaLoc**2 * (2.0_WP * ( S13%values*O31%values - S23%values*O23%values))

    !Calcul de Lij
    call computeFilter(Wave,2.0_WP*deltaLoc,VecX,Sm11,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,VecY,Sm22,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,VecZ,Sm33,filter)
    call computeT_ijVecA(VecX,VecY,VecZ,Sm11,Sm22,Sm33,Wave,S11,S12,S13,S22,S23,S33,filter,res,2.0_WP*deltaLoc)
    Sm11%values = (S11%values-H11%values)*T11%values + (S22%values-H22%values)*T22%values + (S33%values-H33%values)*T33%values + &
                & 2.0_WP*((S12%values-H12%values)*T12%values + (S13%values-H13%values)*T13%values +(S23%values-H23%values)*T23%values)
    Sm22%values = T11%values**2+T22%values**2+T33%values**2+2.0_WP*(T12%values**2+T13%values**2+T23%values**2) 
    coefLoc = computeFieldAvg(Sm11,spec_rank) / computeFieldAvg(Sm22,spec_rank) 
    if (spec_rank .eq. 0) write (6,'(a,1x,g15.8)') '[INFO LES RGMSikmoinsSkjSikOjkOikSjkDyn Fab] Dynanmic coef',coefLoc

    if (.not. computeSplitingSij(S11,S12,S13,S22,S23,S33,Sm11,Sm12,Sm13,Sm22,Sm23,Sm33,'moins')) return
    T11%values = deltaLoc**2/12.0_WP*(Sm11%values*S11%values+Sm12%values*S12%values+Sm13%values*S13%values)+&
         coefLoc*deltaLoc**2*(2.0_WP*(S12%values*O12%values-S13%values*O31%values))
    T12%values = deltaLoc**2/12.0_WP*(Sm11%values*S12%values+Sm12%values*S22%values+Sm13%values*S23%values)-&
         coefLoc*deltaLoc**2*(S11%values*O12%values+S13%values*O23%values+S22%values*O12%values-S23%values*O31%values )
    T13%values = deltaLoc**2/12.0_WP*(Sm11%values*S13%values+Sm12%values*S23%values+Sm13%values*S33%values)+&
         coefLoc*deltaLoc**2*(S11%values*O31%values-S12%values*O23%values+O12%values*S23%values-O31%values*S33%values )
    T22%values = deltaLoc**2/12.0_WP*(Sm12%values*S12%values+Sm22%values*S22%values+Sm23%values*S23%values)+&
         coefLoc*deltaLoc**2*(2.0_WP*(S23%values*O23%values-S12%values*O12%values))
    T23%values = deltaLoc**2/12.0_WP*(Sm12%values*S13%values+Sm22%values*S23%values+Sm23%values*S33%values)+&
         coefLoc*deltaLoc**2*(S12%values*O31%values-S22%values*O23%values-S13%values*O12%values+S33%values*O23%values )
    T33%values = deltaLoc**2/12.0_WP*(Sm13%values*S13%values+Sm23%values*S23%values+Sm33%values*S33%values)+&
         coefLoc*deltaLoc**2*(2.0_WP*(S13%values*O31%values-S23%values*O23%values))
    call ftran(T11,Ufilk,res)
    call ftran(T12,Vfilk,res)
    call ftran(T13,Wfilk,res)
    call computeDivergenceK(Wave,Ufilk,Vfilk,Wfilk,T1)
    call ftran(T22,Ufilk,res)
    call ftran(T23,Wfilk,res)
    call computeDivergenceK(Wave,UfilK,Vfilk,Wfilk,T2)
    call ftran(T13,Vfilk,res)
    call ftran(T33,Ufilk,res)
    call computeDivergenceK(Wave,Vfilk,Ufilk,Wfilk,T3)
    if ( pdffield ) then
      call btran(T1,T11,res)
      call btran(T2,T22,res)
      call btran(T3,T33,res)
      T11%name="sgs_RGMOikOjk_T1"
      T22%name="sgs_RGMOikOjk_T2"
      T33%name="sgs_RGMOikOjk_T3"
      if (.not. computeFieldPDF(ite,T11,300,spec_rank,10000)) return
      if (.not. computeFieldPDF(ite,T22,300,spec_rank,10000)) return
      if (.not. computeFieldPDF(ite,T33,300,spec_rank,10000)) return
    endif

    ! ===== Free memory and delete temporary fields =====
    call deleteDataLayout(H11)
    call deleteDataLayout(H12)
    call deleteDataLayout(H13)
    call deleteDataLayout(H23)
    call deleteDataLayout(H22)
    call deleteDataLayout(H33)
    call deleteDataLayout(T11)
    call deleteDataLayout(T12)
    call deleteDataLayout(T13)
    call deleteDataLayout(T23)
    call deleteDataLayout(T22)
    call deleteDataLayout(T33)
    call deleteDataLayout(Sm11)
    call deleteDataLayout(Sm12)
    call deleteDataLayout(Sm13)
    call deleteDataLayout(Sm23)
    call deleteDataLayout(Sm22)
    call deleteDataLayout(Sm33)
    call deleteDataLayout(S11)
    call deleteDataLayout(S12)
    call deleteDataLayout(S13)
    call deleteDataLayout(S23)
    call deleteDataLayout(S22)
    call deleteDataLayout(S33)
    call deleteDataLayout(O23)
    call deleteDataLayout(O31)
    call deleteDataLayout(O12)
    call deleteDataLayout(Ufilk)
    call deleteDataLayout(Vfilk)
    call deleteDataLayout(Wfilk)

    res = .true.

end subroutine RGMSikmoinsSkjSikOjkOikSjkDynFab



!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> @param [in] VecX velocity component along first direction in spectral space
!> @param [in] VecY velocity component along second direction in spectral space
!> @param [in] VecZ velocity component along third direction in spectral space
!> @param [in] FieldForParam field to get informations about grid (physical space) 
!> @param [in] Wave wave numbers
!> @param [in] spec_rank rank of processor
!> @param [inout] T1 =  delta**2/12 d/dxj ( S1kmoins Skj ) + Cd delta**2 d/dxj ( S1k Ojk + O1k Sjk )
!> @param [inout] T2 =  delta**2/12 d/dxj ( S2kmoins Skj ) + Cd delta**2 d/dxj ( S2k Ojk + O2k Sjk )
!> @param [inout] T3 =  delta**2/12 d/dxj ( S3kmoins Skj ) + Cd delta**2 d/dxj ( S3k Ojk + O3k Sjk )
!------------------------------------------------------------------------------
subroutine RGMSikmoinsSkjSikOjkOikSjkDyn(T1,T2,T3,VecX,VecY,VecZ,Wave,spec_rank,filter,delta,ite)

   use stat_tools , only : equalToVal
  
   TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: VecX,VecY,VecZ
   TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT) :: T1,T2,T3
   TYPE(WaveNumbers), INTENT(IN) :: Wave
   INTEGER,INTENT(IN)            :: spec_rank,filter
   REAL(WP),INTENT(IN),OPTIONAL  :: delta
   INTEGER,INTENT(IN),OPTIONAL   :: ite
   LOGICAL                       :: res 

    ! Local field
    type(REAL_DATA_LAYOUT) :: Sm11,Sm12,Sm13,Sm22,Sm23,Sm33
    type(REAL_DATA_LAYOUT) :: S11,S12,S13,S22,S23,S33
    type(REAL_DATA_LAYOUT) :: T11,T12,T13,T22,T23,T33
    type(REAL_DATA_LAYOUT) :: L11,L12,L13,L22,L23,L33
    type(REAL_DATA_LAYOUT) :: H11,H12,H13,H22,H23,H33
    type(REAL_DATA_LAYOUT) :: O23,O31,O12
    type(COMPLEX_DATA_LAYOUT)::Ufilk,Vfilk,Wfilk 
    integer                  :: nbcpus
!    real(WP) :: Num,Den
    real(WP) :: coefLoc,deltaLoc
    logical :: pdffield

    ! ===== Initialisation =====
    res = .false.
    nbcpus = getnbcpus()
    pdffield=.false.
    if (present(ite)) then
      if ( ite .ne. -1 ) pdffield = .true.
    endif
    if (present(delta)) then
      deltaLoc=delta
    else
      deltaLoc=computeDelta(VecX)
    endif 
    if (.not. initDataLayout("UFilk", UFilk,(VecX%nx/2)+1,VecX%ny,VecX%nz, &
       & VecX%Lx,VecX%Ly,VecX%Lz,nbcpus,spec_rank,alongZ)) then 
      write(6,'(a,i0,a)')'[ERROR] Not able to create Exk'
      return
    endif
    if (.not. copystructonly(VecX,T11) ) return
    if (.not. copystructonly(VecX,T12) ) return
    if (.not. copystructonly(VecX,T13) ) return
    if (.not. copystructonly(VecX,T22) ) return
    if (.not. copystructonly(VecX,T23) ) return
    if (.not. copystructonly(VecX,T33) ) return
    if (.not. copystructonly(VecX,L11) ) return
    if (.not. copystructonly(VecX,L12) ) return
    if (.not. copystructonly(VecX,L13) ) return
    if (.not. copystructonly(VecX,L22) ) return
    if (.not. copystructonly(VecX,L23) ) return
    if (.not. copystructonly(VecX,L33) ) return
    if (.not. copystructonly(VecX,H11) ) return
    if (.not. copystructonly(VecX,H12) ) return
    if (.not. copystructonly(VecX,H13) ) return
    if (.not. copystructonly(VecX,H22) ) return
    if (.not. copystructonly(VecX,H23) ) return
    if (.not. copystructonly(VecX,H33) ) return
    if (.not. copystructonly(VecX,Sm11) ) return
    if (.not. copystructonly(VecX,Sm12) ) return
    if (.not. copystructonly(VecX,Sm13) ) return
    if (.not. copystructonly(VecX,Sm22 ) )return
    if (.not. copystructonly(VecX,Sm23) ) return
    if (.not. copystructonly(VecX,Sm33) ) return
    if (.not. copystructonly(VecX,S11) ) return
    if (.not. copystructonly(VecX,S12) ) return
    if (.not. copystructonly(VecX,S13) ) return
    if (.not. copystructonly(VecX,S22) ) return
    if (.not. copystructonly(VecX,S23) ) return
    if (.not. copystructonly(VecX,S33) ) return
    if (.not. copystructonly(VecX,O23) )return
    if (.not. copystructonly(VecX,O31) ) return
    if (.not. copystructonly(VecX,O12) ) return
    if (.not. copystructonly(Ufilk,Vfilk) ) return
    if (.not. copystructonly(Ufilk,Wfilk) ) return
    

    !Calcul du coef dynamique
    !Cd = < (Lij-Hij)Tij>/<Tij Tij>
    !avec
    ! Lij=(u|_i u|_j)^ - u|^_i u|^_j
    ! Hij = delta|^**2/12 S|^moinsik S|^jk - delta|**2/12 (S|moinsik S|kj)^
    ! Tij = delta|^**2(S|^ik O|^jk + O|^ik S|^jk) -delta|**2(S|ik O|jk + O|ik S|jk)^
    
    !Calcul de Hij
    call ftran(VecX,Ufilk,res)
    call ftran(VecY,Vfilk,res)
    call ftran(VecZ,Wfilk,res)
    call computeStressTensor(Ufilk,Vfilk,Wfilk,Wave,S11,S12,S13,S22,S23,S33,res)
    if (.not. computeSplitingSij(S11,S12,S13,S22,S23,S33,Sm11,Sm12,Sm13,Sm22,Sm23,Sm33,'moins')) return
    T11%values = Sm11%values*S11%values+Sm12%values*S12%values+Sm13%values*S13%values
    T12%values = Sm11%values*S12%values+Sm12%values*S22%values+Sm13%values*S23%values
    T13%values = Sm11%values*S13%values+Sm12%values*S23%values+Sm13%values*S33%values
    T22%values = Sm12%values*S12%values+Sm22%values*S22%values+Sm23%values*S23%values
    T23%values = Sm12%values*S13%values+Sm22%values*S23%values+Sm23%values*S33%values
    T33%values = Sm13%values*S13%values+Sm23%values*S23%values+Sm33%values*S33%values
    call computeFilter(Wave,2.0_WP*deltaLoc,T11,H11,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,T12,H12,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,T13,H13,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,T22,H22,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,T23,H23,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,T33,H33,filter)

    call computeFilter(Wave,2.0_WP*deltaLoc,VecX,Ufilk,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,VecY,Vfilk,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,VecZ,Wfilk,filter)
    call computeStressTensor(Ufilk,Vfilk,Wfilk,Wave,S11,S12,S13,S22,S23,S33,res)
    if (.not. computeSplitingSij(S11,S12,S13,S22,S23,S33,Sm11,Sm12,Sm13,Sm22,Sm23,Sm33,'moins')) return
    T11%values = Sm11%values*S11%values+Sm12%values*S12%values+Sm13%values*S13%values
    T12%values = Sm11%values*S12%values+Sm12%values*S22%values+Sm13%values*S23%values
    T13%values = Sm11%values*S13%values+Sm12%values*S23%values+Sm13%values*S33%values
    T22%values = Sm12%values*S12%values+Sm22%values*S22%values+Sm23%values*S23%values
    T23%values = Sm12%values*S13%values+Sm22%values*S23%values+Sm23%values*S33%values
    T33%values = Sm13%values*S13%values+Sm23%values*S23%values+Sm33%values*S33%values
    H11%values = (2.0_WP*deltaLoc)**2/12.0_WP * T11%values - deltaLoc**2/12.0_WP * H11%values
    H12%values = (2.0_WP*deltaLoc)**2/12.0_WP * T12%values - deltaLoc**2/12.0_WP * H12%values
    H13%values = (2.0_WP*deltaLoc)**2/12.0_WP * T13%values - deltaLoc**2/12.0_WP * H13%values
    H22%values = (2.0_WP*deltaLoc)**2/12.0_WP * T22%values - deltaLoc**2/12.0_WP * H22%values
    H23%values = (2.0_WP*deltaLoc)**2/12.0_WP * T23%values - deltaLoc**2/12.0_WP * H23%values
    H33%values = (2.0_WP*deltaLoc)**2/12.0_WP * T33%values - deltaLoc**2/12.0_WP * H33%values

    !Calcul de Mij
    ! Tij = delta|^**2(S|^ik O|^jk + O|^ik S|^jk) - delta|**2(S|ik O|jk + O|ik S|jk)^
    call computeVorticity(Ufilk,Vfilk,Wfilk,Wave,O23,O31,O12,res)
    T11%values = 2.0_WP * ( S12%values*O12%values - S13%values*O31%values) 
    T12%values = S13%values*O23%values+S22%values*O12%values-S23%values*O31%values-S11%values*O12%values
    T13%values = S11%values*O31%values-S12%values*O23%values+O12%values*S23%values-O31%values*S33%values
    T22%values = 2.0_WP * ( S23%values*O23%values - S12%values*O12%values)
    T23%values = S12%values*O31%values-S22%values*O23%values-S13%values*O12%values+S33%values*O23%values
    T33%values = 2.0_WP * ( S13%values*O31%values - S23%values*O23%values)
    call ftran(VecX,Ufilk,res)
    call ftran(VecY,Vfilk,res)
    call ftran(VecZ,Wfilk,res)
    call computeStressTensor(Ufilk,Vfilk,Wfilk,Wave,S11,S12,S13,S22,S23,S33,res)
    call computeVorticity(Ufilk,Vfilk,Wfilk,Wave,O23,O31,O12,res)
    Sm11%values = 2.0_WP * ( S12%values*O12%values - S13%values*O31%values) 
    Sm12%values = S13%values*O23%values+S22%values*O12%values-S23%values*O31%values-S11%values*O12%values
    Sm13%values = S11%values*O31%values-S12%values*O23%values+O12%values*S23%values-O31%values*S33%values
    Sm22%values = 2.0_WP * ( S23%values*O23%values - S12%values*O12%values)
    Sm23%values = S12%values*O31%values-S22%values*O23%values-S13%values*O12%values+S33%values*O23%values
    Sm33%values = 2.0_WP * ( S13%values*O31%values - S23%values*O23%values)
    call computeFilter(Wave,2.0_WP*deltaLoc,Sm11,L11,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,Sm12,L12,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,Sm13,L13,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,Sm22,L22,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,Sm23,L23,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,Sm33,L33,filter)
    T11%values = (2.0_WP*deltaLoc)**2 * T11%values - deltaLoc**2 * L11%values
    T12%values = (2.0_WP*deltaLoc)**2 * T12%values - deltaLoc**2 * L12%values
    T13%values = (2.0_WP*deltaLoc)**2 * T13%values - deltaLoc**2 * L13%values
    T22%values = (2.0_WP*deltaLoc)**2 * T22%values - deltaLoc**2 * L22%values
    T23%values = (2.0_WP*deltaLoc)**2 * T23%values - deltaLoc**2 * L23%values
    T33%values = (2.0_WP*deltaLoc)**2 * T33%values - deltaLoc**2 * L33%values

    !Calcul de Lij
    call computeFilter(Wave,2.0_WP*deltaLoc,VecX,Sm11,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,VecY,Sm22,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,VecZ,Sm33,filter)
!        computeT_ijVecA(in1 ,in2 ,in3 ,in1f,in2f,in3f,wave,out11,out12,out13,out22,out23,out33,filter,success,delta)
    call computeT_ijVecA(VecX,VecY,VecZ,Sm11,Sm22,Sm33,Wave,L11  ,L12  ,L13  ,L22  ,L23  ,L33  ,filter,res    ,2.0_WP*deltaLoc)

    Sm11%values = (L11%values-H11%values)*T11%values + (L22%values-H22%values)*T22%values + (L33%values-H33%values)*T33%values + &
                & 2.0_WP*((L12%values-H12%values)*T12%values + (L13%values-H13%values)*T13%values +(L23%values-H23%values)*T23%values)
    Sm22%values = T11%values**2+T22%values**2+T33%values**2+2.0_WP*(T12%values**2+T13%values**2+T23%values**2) 
    coefLoc = computeFieldAvg(Sm11,spec_rank) / computeFieldAvg(Sm22,spec_rank) 
    if (spec_rank .eq. 0) write (6,'(a,1x,g15.8)') '[INFO LES RGMSikmoinsSkjSikOjkOikSjkDyn] Dynanmic coef',coefLoc

!                                         !out11,out12,out13,out22,out23,out33
!    call computeStressTensor(Ufilk,Vfilk,Wfilk,Wave,S11,S12,S13,S22,S23,S33,res)
    if (.not. computeSplitingSij(S11,S12,S13,S22,S23,S33,Sm11,Sm12,Sm13,Sm22,Sm23,Sm33,'moins')) return
!    call computeVorticity(Ufilk,Vfilk,Wfilk,Wave,O23,O31,O12,res)
    T11%values = deltaLoc**2/12.0_WP*(Sm11%values*S11%values+Sm12%values*S12%values+Sm13%values*S13%values)+&
         coefLoc*deltaLoc**2*(2.0_WP*(S12%values*O12%values-S13%values*O31%values))
    T12%values = deltaLoc**2/12.0_WP*(Sm11%values*S12%values+Sm12%values*S22%values+Sm13%values*S23%values)-&
         coefLoc*deltaLoc**2*(S11%values*O12%values+S13%values*O23%values+S22%values*O12%values-S23%values*O31%values )
    T13%values = deltaLoc**2/12.0_WP*(Sm11%values*S13%values+Sm12%values*S23%values+Sm13%values*S33%values)+&
         coefLoc*deltaLoc**2*(S11%values*O31%values-S12%values*O23%values+O12%values*S23%values-O31%values*S33%values )
    T22%values = deltaLoc**2/12.0_WP*(Sm12%values*S12%values+Sm22%values*S22%values+Sm23%values*S23%values)+&
         coefLoc*deltaLoc**2*(2.0_WP*(S23%values*O23%values-S12%values*O12%values))
    T23%values = deltaLoc**2/12.0_WP*(Sm12%values*S13%values+Sm22%values*S23%values+Sm23%values*S33%values)+&
         coefLoc*deltaLoc**2*(S12%values*O31%values-S22%values*O23%values-S13%values*O12%values+S33%values*O23%values )
    T33%values = deltaLoc**2/12.0_WP*(Sm13%values*S13%values+Sm23%values*S23%values+Sm33%values*S33%values)+&
         coefLoc*deltaLoc**2*(2.0_WP*(S13%values*O31%values-S23%values*O23%values))
    call ftran(T11,Ufilk,res)
    call ftran(T12,Vfilk,res)
    call ftran(T13,Wfilk,res)
    call computeDivergenceK(Wave,Ufilk,Vfilk,Wfilk,T1)
    call ftran(T22,Ufilk,res)
    call ftran(T23,Wfilk,res)
    call computeDivergenceK(Wave,UfilK,Vfilk,Wfilk,T2)
    call ftran(T13,Vfilk,res)
    call ftran(T33,Ufilk,res)
    call computeDivergenceK(Wave,Vfilk,Ufilk,Wfilk,T3)
    if ( pdffield ) then
      call btran(T1,T11,res)
      call btran(T2,T22,res)
      call btran(T3,T33,res)
      T11%name="sgs_RGMOikOjk_T1"
      T22%name="sgs_RGMOikOjk_T2"
      T33%name="sgs_RGMOikOjk_T3"
      if (.not. computeFieldPDF(ite,T11,300,spec_rank,10000)) return
      if (.not. computeFieldPDF(ite,T22,300,spec_rank,10000)) return
      if (.not. computeFieldPDF(ite,T33,300,spec_rank,10000)) return
    endif

    ! ===== Free memory and delete temporary fields =====
    call deleteDataLayout(L11)
    call deleteDataLayout(L12)
    call deleteDataLayout(L13)
    call deleteDataLayout(L23)
    call deleteDataLayout(L22)
    call deleteDataLayout(L33)
    call deleteDataLayout(H11)
    call deleteDataLayout(H12)
    call deleteDataLayout(H13)
    call deleteDataLayout(H23)
    call deleteDataLayout(H22)
    call deleteDataLayout(H33)
    call deleteDataLayout(T11)
    call deleteDataLayout(T12)
    call deleteDataLayout(T13)
    call deleteDataLayout(T23)
    call deleteDataLayout(T22)
    call deleteDataLayout(T33)
    call deleteDataLayout(Sm11)
    call deleteDataLayout(Sm12)
    call deleteDataLayout(Sm13)
    call deleteDataLayout(Sm23)
    call deleteDataLayout(Sm22)
    call deleteDataLayout(Sm33)
    call deleteDataLayout(S11)
    call deleteDataLayout(S12)
    call deleteDataLayout(S13)
    call deleteDataLayout(S23)
    call deleteDataLayout(S22)
    call deleteDataLayout(S33)
    call deleteDataLayout(O23)
    call deleteDataLayout(O31)
    call deleteDataLayout(O12)
    call deleteDataLayout(Ufilk)
    call deleteDataLayout(Vfilk)
    call deleteDataLayout(Wfilk)

    res = .true.

end subroutine RGMSikmoinsSkjSikOjkOikSjkDyn

!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> @param [in] VecX  velocity component along first direction in spectral space
!> @param [in] VecXk velocity component along first direction in spectral space
!> @param [in] VecYk velocity component along second direction in spectral space
!> @param [in] VecZk velocity component along third direction in spectral space
!> @param [in] 
!> @param [in] Wave wave numbers
!> @param [in] spec_rank rank of processor
!> @param [inout] 
!> @param [inout] 
!> @param [inout] 
!------------------------------------------------------------------------------
subroutine RGMSikmoinsSkjSikOjkOikSjk_div(T1,T2,T3,VecX,VecXk,VecYk,VecZk,Wave,spec_rank,delta,coef1,coef2,ite)

    !I/O data
    type(real_data_layout),intent(in)       :: VecX
    type(complex_data_layout),intent(in)    :: VecXk,VecYk,VecZk
    type(complex_data_layout),intent(inout) :: T1,T2,T3
    type(WaveNumbers),intent(in)            :: Wave
    integer,intent(in)                      :: spec_rank
    real(WP),intent(in),optional            :: coef1,coef2,delta
    integer,intent(in),optional             :: ite

    !Local data
    type(real_data_layout)  :: T11,T12,T13,T22,T23,T33
    type(real_data_layout)  :: F11,F12,F13,F22,F23,F33
    logical                 :: pdffield,res
    real(WP)                :: coefLoc1,coefLoc2,deltaLoc

    pdffield= .false.
    if (present(ite)) then
      if ( ite .ne. -1 ) pdffield = .true.
    endif
    if (present(coef1)) then
      coefLoc1 = coef1
    else
      coefLoc1 = 1.0_WP/12.0_WP
    endif
    if (present(coef2)) then
      coefLoc2 = coef2
    else
      coefLoc2 = 1.0_WP/12.0_WP
    endif
    if (present(delta)) then
      deltaLoc = delta
    else
      deltaLoc = computeDelta(VecX)
    endif

    if (.not. copystructonly(VecX,T11) ) return
    if (.not. copystructonly(VecX,T12) ) return
    if (.not. copystructonly(VecX,T13) ) return
    if (.not. copystructonly(VecX,T22) ) return
    if (.not. copystructonly(VecX,T23) ) return
    if (.not. copystructonly(VecX,T33) ) return
    if (spec_rank .eq. 0) write (6,'(a)') '[INFO LES] sgs model for velocity: RGMSikmoinsSkj'

    call RGMSikmoinsSkjSikOjkOikSjk_tij(T11,T12,T13,T22,T23,T33,VecX,VecXk,VecYk,VecZk,Wave,spec_rank,deltaLoc,coefLoc1,coefLoc2)
    call ComputeDivSymTijVec1Vec2(T1,T2,T3,T12,T13,T23,T11,T22,T33,Wave,res)
    if (.not. res) write(6,'(a)') '[ERROR] in RGMSikmoinsSkj_div'
    if ( pdffield ) then
      if (.not. copystructonly(VecX,F11) ) return
      if (.not. copystructonly(VecX,F12) ) return
      if (.not. copystructonly(VecX,F13) ) return
      if (.not. copystructonly(VecX,F22) ) return
      if (.not. copystructonly(VecX,F23) ) return
      if (.not. copystructonly(VecX,F33) ) return
      call computeStressTensor(VecXk,VecYk,VecZk,Wave,F11,F12,F13,F22,F23,F33,res)
      T11%values = T11%values*F11%values + T22%values*F22%values + T33%values*F33%values + &
         2.0_WP * (T12%values*F12%values + T13%values*F13%values + T23%values*F23%values )
      T11%name="diss_tijSij"
      coefLoc1 = computeFieldAvg(T11,spec_rank)
      if (spec_rank .eq. 0) write (6,'(a,1x,g15.8)') '[INFO LES] RGM SikOjk+OikSjk dissipation over the box ',coefLoc1
      call btran(T1,T11,res)
      call btran(T2,T22,res)
      call btran(T3,T33,res)
      T11%name="sgs_RGMS1kOjk+O1kSjk_T1"
      T22%name="sgs_RGMS1kOjk+O2kSjk_T2"
      T33%name="sgs_RGMS1kOjk+O3kSjk_T3"
      if (.not. computeFieldPDF(ite,T11,300,spec_rank,10000)) return
      if (.not. computeFieldPDF(ite,T22,300,spec_rank,10000)) return
      if (.not. computeFieldPDF(ite,T33,300,spec_rank,10000)) return
      call deleteDataLayout(F11)
      call deleteDataLayout(F12)
      call deleteDataLayout(F13)
      call deleteDataLayout(F23)
      call deleteDataLayout(F22)
      call deleteDataLayout(F33)
    endif

    call deleteDataLayout(T11)
    call deleteDataLayout(T12)
    call deleteDataLayout(T13)
    call deleteDataLayout(T23)
    call deleteDataLayout(T22)
    call deleteDataLayout(T33)

end subroutine RGMSikmoinsSkjSikOjkOikSjk_div 


!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> @param [in] VecX velocity component along first direction in spectral space
!> @param [in] VecXk velocity component along first direction in spectral space
!> @param [in] VecYk velocity component along second direction in spectral space
!> @param [in] VecZk velocity component along third direction in spectral space
!> @param [in] Wave wave numbers
!> @param [in] spec_rank rank of processor
!> @param [in] coef1 
!> @param [in] coef2
!> @param [inout] T11 = delta**2 ( coef1*S1k-S1k + coef2*(S1kO1k+O1kS1k) )
!> @param [inout] T22 = delta**2 ( coef1*S2k-S2k + coef2*(S2kO2k+O2kS2k) )
!> @param [inout] T33 = delta**2 ( coef1*S3k-S3k + coef2*(S3kO3k+O3kS3k) )
!> @param [inout] T12 = delta**2 ( coef1*S1k-S2k + coef2*(S1kO2k+O1kS2k) )
!> @param [inout] T13 = delta**2 ( coef1*S1k-S3k + coef2*(S1kO3k+O1kS3k) )
!> @param [inout] T23 = delta**2 ( coef1*S2k-S3k + coef2*(S2kO3k+O2kS3k) )
!------------------------------------------------------------------------------
subroutine RGMSikmoinsSkjSikOjkOikSjk_tij(T11,T12,T13,T22,T23,T33,VecX,VecXk,VecYk,VecZk,Wave,spec_rank,delta,coef1,coef2)

    TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: VecX
    TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: VecXk,VecYk,VecZk
    TYPE(REAL_DATA_LAYOUT), INTENT(INOUT) :: T11,T12,T13,T22,T23,T33
    TYPE(WaveNumbers), INTENT(IN) :: Wave
    INTEGER,INTENT(IN)            :: spec_rank
    REAL(WP),INTENT(IN)           :: delta,coef1,coef2

    ! Local field
    TYPE(REAL_DATA_LAYOUT)   :: S11,S12,S13,S22,S23,S33
    LOGICAL                  :: res

    ! ===== Initialisation =====
    res = .false.
    if (.not. copystructonly(VecX,S11) ) return
    if (.not. copystructonly(VecX,S12) ) return
    if (.not. copystructonly(VecX,S13) ) return
    if (.not. copystructonly(VecX,S22) ) return
    if (.not. copystructonly(VecX,S23) ) return
    if (.not. copystructonly(VecX,S33) ) return
    
    call RGMSikOjkOikSjk_tij(S11,S12,S13,S22,S23,S33,VecX,VecXk,VecYk,VecZk,Wave,spec_rank,1.0_WP,1.0_WP)
    call RGMSikmoinsSkj_tij (T11,T12,T13,T22,T23,T33,VecX,VecXk,VecYk,VecZk,Wave,spec_rank,1.0_WP,1.0_WP)
    T11%values = delta**2*(coef1*T11%values + coef2*S11%values )
    T12%values = delta**2*(coef1*T12%values + coef2*S12%values )
    T13%values = delta**2*(coef1*T13%values + coef2*S23%values )
    T22%values = delta**2*(coef1*T22%values + coef2*S22%values )
    T23%values = delta**2*(coef1*T23%values + coef2*S23%values )
    T33%values = delta**2*(coef1*T33%values + coef2*S33%values )

    ! ===== Free memory and delete temporary fields =====
    call deleteDataLayout(S11)
    call deleteDataLayout(S12)
    call deleteDataLayout(S13)
    call deleteDataLayout(S23)
    call deleteDataLayout(S22)
    call deleteDataLayout(S33)

    res = .true.

end subroutine RGMSikmoinsSkjSikOjkOikSjk_tij

!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> @param [in] VecX  velocity component along first direction in spectral space
!> @param [in] VecXk velocity component along first direction in spectral space
!> @param [in] VecYk velocity component along second direction in spectral space
!> @param [in] VecZk velocity component along third direction in spectral space
!> @param [in] 
!> @param [in] Wave wave numbers
!> @param [in] spec_rank rank of processor
!> @param [inout] 
!> @param [inout] 
!> @param [inout] 
!------------------------------------------------------------------------------
subroutine RGMSikOjkOikSjk_div(T1,T2,T3,VecX,VecXk,VecYk,VecZk,Wave,spec_rank,delta,coef,ite)

    !I/O data
    type(real_data_layout),intent(in)       :: VecX
    type(complex_data_layout),intent(in)    :: VecXk,VecYk,VecZk
    type(complex_data_layout),intent(inout) :: T1,T2,T3
    type(WaveNumbers),intent(in)            :: Wave
    integer,intent(in)                      :: spec_rank
    real(WP),intent(in),optional            :: coef,delta
    integer,intent(in),optional             :: ite

    !Local data
    type(real_data_layout)  :: T11,T12,T13,T22,T23,T33
    type(real_data_layout)  :: F11,F12,F13,F22,F23,F33
    logical                 :: pdffield,res
    real(WP)                :: coefLoc,deltaLoc

    pdffield= .false.
    if (present(ite)) then
      if ( ite .ne. -1 ) pdffield = .true.
    endif
    if (present(coef)) then
      coefLoc = coef
    else
      coefLoc = 1.0_WP/12.0_WP
    endif
    if (present(delta)) then
      deltaLoc = delta
    else
      deltaLoc = computeDelta(VecX)
    endif

    if (.not. copystructonly(VecX,T11) ) return
    if (.not. copystructonly(VecX,T12) ) return
    if (.not. copystructonly(VecX,T13) ) return
    if (.not. copystructonly(VecX,T22) ) return
    if (.not. copystructonly(VecX,T23) ) return
    if (.not. copystructonly(VecX,T33) ) return
    if (spec_rank .eq. 0) write (6,'(a)') '[INFO LES] sgs model for velocity: RGMSikmoinsSkj'

    call RGMSikOjkOikSjk_tij(T11,T12,T13,T22,T23,T33,VecX,VecXk,VecYk,VecZk,Wave,spec_rank,deltaLoc,coefLoc)
    call ComputeDivSymTijVec1Vec2(T1,T2,T3,T12,T13,T23,T11,T22,T33,Wave,res)
    if (.not. res) write(6,'(a)') '[ERROR] in RGMSikmoinsSkj_div'
    if ( pdffield ) then
      if (.not. copystructonly(VecX,F11) ) return
      if (.not. copystructonly(VecX,F12) ) return
      if (.not. copystructonly(VecX,F13) ) return
      if (.not. copystructonly(VecX,F22) ) return
      if (.not. copystructonly(VecX,F23) ) return
      if (.not. copystructonly(VecX,F33) ) return
      call computeStressTensor(VecXk,VecYk,VecZk,Wave,F11,F12,F13,F22,F23,F33,res)
      T11%values = T11%values*F11%values + T22%values*F22%values + T33%values*F33%values + &
         2.0_WP * (T12%values*F12%values + T13%values*F13%values + T23%values*F23%values )
      T11%name="diss_tijSij"
      coefLoc = computeFieldAvg(T11,spec_rank)
      if (spec_rank .eq. 0) write (6,'(a,1x,g15.8)') '[INFO LES] RGM SikOjk+OikSjk dissipation over the box ',coefLoc
      call btran(T1,T11,res)
      call btran(T2,T22,res)
      call btran(T3,T33,res)
      T11%name="sgs_RGMS1kOjk+O1kSjk_T1"
      T22%name="sgs_RGMS1kOjk+O2kSjk_T2"
      T33%name="sgs_RGMS1kOjk+O3kSjk_T3"
      if (.not. computeFieldPDF(ite,T11,300,spec_rank,10000)) return
      if (.not. computeFieldPDF(ite,T22,300,spec_rank,10000)) return
      if (.not. computeFieldPDF(ite,T33,300,spec_rank,10000)) return
      call deleteDataLayout(F11)
      call deleteDataLayout(F12)
      call deleteDataLayout(F13)
      call deleteDataLayout(F23)
      call deleteDataLayout(F22)
      call deleteDataLayout(F33)
    endif

    call deleteDataLayout(T11)
    call deleteDataLayout(T12)
    call deleteDataLayout(T13)
    call deleteDataLayout(T23)
    call deleteDataLayout(T22)
    call deleteDataLayout(T33)

end subroutine RGMSikOjkOikSjk_div
!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> @param [in] VecX velocity component along first direction in spectral space
!> @param [in] VecXk velocity component along first direction in spectral space
!> @param [in] VecYk velocity component along second direction in spectral space
!> @param [in] VecZk velocity component along third direction in spectral space
!> @param [in] FieldForParam field to get informations about grid (physical space) 
!> @param [in] Wave wave numbers
!> @param [in] spec_rank rank of processor
!> @param [in] delta is the lenght of filter
!> @param [in] coef is the coefficient in front of the model (optional) 
!> @param [inout] T11 = coef * delta**2 ( S1kS1k + O1kS1k )
!> @param [inout] T22 = coef * delta**2 ( S2kS2k + O2kS2k )
!> @param [inout] T33 = coef * delta**2 ( S3kS3k + O3kS3k )
!> @param [inout] T12 = coef * delta**2 ( S1kS2k + O1kS2k )
!> @param [inout] T13 = coef * delta**2 ( S2kS3k + O2kS3k )
!> @param [inout] T23 = coef * delta**2 ( S3kS3k + O3kS3k )
!------------------------------------------------------------------------------
subroutine RGMSikOjkOikSjk_tij(T11,T12,T13,T22,T23,T33,VecX,VecXk,VecYk,VecZk,Wave,spec_rank,delta,coef)

    TYPE(REAL_DATA_LAYOUT), INTENT(IN)     :: VecX
    TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN)  :: VecXk,VecYk,VecZk
    TYPE(REAL_DATA_LAYOUT),intent(inout)   :: T11,T12,T13,T22,T23,T33
    TYPE(WaveNumbers), INTENT(IN)          :: Wave
    INTEGER,INTENT(IN)                     :: spec_rank
    REAL(WP),INTENT(IN)                    :: delta,coef

    ! Local field
    type(REAL_DATA_LAYOUT) :: S11,S12,S13,S22,S23,S33
    type(REAL_DATA_LAYOUT) :: O23,O31,O12
    logical                :: res
  
    ! ===== Initialisation =====
    res = .false.
    if (.not. copystructonly(VecX,S11) ) return
    if (.not. copystructonly(VecX,S12) ) return
    if (.not. copystructonly(VecX,S13) ) return
    if (.not. copystructonly(VecX,S22) ) return
    if (.not. copystructonly(VecX,S23) ) return
    if (.not. copystructonly(VecX,S33) ) return
    if (.not. copystructonly(VecX,O23) )return
    if (.not. copystructonly(VecX,O31) ) return
    if (.not. copystructonly(VecX,O12) ) return
    
    call computeStressTensor(VecXk,VecYk,VecZk,Wave,S11,S12,S13,S22,S23,S33,res)
    call computeVorticity(VecXk,VecYk,VecZk,Wave,O23,O31,O12,res)
    T11%values = 2.0_WP*(S12%values*O12%values - S13%values*O31%values)
    T11%values = coef*delta**2 * T11%values 
    T22%values = 2.0_WP*(S23%values*O23%values - S12%values*O12%values)
    T22%values = coef*delta**2 * T22%values 
    T33%values = 2.0_WP*(S13%values*O31%values - S23%values*O23%values)
    T33%values = coef*delta**2 * T33%values
    T12%values = S13%values*O23%values-S11%values*O12%values+S22%values*O12%values-S23%values*O31%values 
    T12%values = coef*delta**2 * T12%values
    T13%values = S11%values*O31%values-S12%values*O23%values+O12%values*S23%values-O31%values*S33%values
    T13%values = coef*delta**2 * T13%values
    T23%values = S12%values*O31%values-S22%values*O23%values-S13%values*O12%values+S33%values*O23%values 
    T23%values = coef*delta**2 * T23%values

    ! ===== Free memory and delete temporary fields =====
    call deleteDataLayout(S11)
    call deleteDataLayout(S12)
    call deleteDataLayout(S13)
    call deleteDataLayout(S23)
    call deleteDataLayout(S22)
    call deleteDataLayout(S33)
    call deleteDataLayout(O23)
    call deleteDataLayout(O31)
    call deleteDataLayout(O12)

    res = .true.

end subroutine RGMSikOjkOikSjk_tij

!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> Modèle du gradient en prenant en compte l'incompressibilité. Calcul de 3 termes par composantes du
!> tenseur sous maille au lieu de 6. Calcul de d/dxj(tau_ij)=Ti
!> d/dxj(tau_ij) = C delta^2 (duj/dxk d^2ui/dxjdxk) 
!> @param [in] VecX velocity component along first direction in spectral space
!> @param [in] VecY velocity component along second direction in spectral space
!> @param [in] VecZ velocity component along third direction in spectral space
!> @param [in] FieldForParam field to get informations about grid (physical space) 
!> @param [in] Wave wave numbers
!> @param [in] spec_rank rank of processor
!> @param [inout] T1 = eps_1jk eps_jmn C delta^2 duk/dxl d^2un/dxmdxl
!> @param [inout] T2 = eps_2jk eps_jmn C delta^2 duk/dxl d^2un/dxmdxl
!> @param [inout] T3 = eps_3jk eps_jmn C delta^2 duk/dxl d^2un/dxmdxl
!------------------------------------------------------------------------------
!subroutine RGMSikmoinsSkjDynFab(T1,T2,T3,VecX,VecY,VecZ,Wave,spec_rank,filter,delta,ite)
subroutine RGMSikmoinsSkjDynFab(T1,T2,T3,VecX,VecXk,VecYk,VecZk,Wave,spec_rank,filter,delta,ite)

  
   TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: VecX
   TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: VecXk,VecYk,VecZk
   TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT) :: T1,T2,T3
   TYPE(WaveNumbers), INTENT(IN) :: Wave
   INTEGER,INTENT(IN)            :: spec_rank
   INTEGER,INTENT(IN)            :: filter 
   REAL(WP),INTENT(IN),OPTIONAL  :: delta
   INTEGER,INTENT(IN),OPTIONAL   :: ite

    ! Local field
    type(REAL_DATA_LAYOUT) :: Sm11,Sm12,Sm13,Sm22,Sm23,Sm33
    type(REAL_DATA_LAYOUT) :: S11,S12,S13,S22,S23,S33
    type(REAL_DATA_LAYOUT) :: T11,T12,T13,T22,T23,T33
    type(COMPLEX_DATA_LAYOUT)::Ufilk,Vfilk,Wfilk 
    integer                  :: nbcpus
    real(WP) :: coefLoc,deltaLoc
    logical :: res,pdffield

    ! ===== Initialisation =====
    res = .false.
    pdffield=.false.
    nbcpus = getnbcpus()
    if (present(ite)) then
      if ( ite .ne. -1 ) pdffield = .true.
    endif
    if (present(delta)) then
      deltaLoc=delta
    else
      deltaLoc=computeDelta(VecX)
    endif 
    if (.not. initDataLayout("UFilk", UFilk,(VecX%nx/2)+1,VecX%ny,VecX%nz, &
       & VecX%Lx,VecX%Ly,VecX%Lz,nbcpus,spec_rank,alongZ)) then 
      write(6,'(a,i0,a)')'[ERROR] Not able to create Exk'
      return
    endif
    if (.not. copystructonly(VecX,T11) ) return
    if (.not. copystructonly(VecX,T12) ) return
    if (.not. copystructonly(VecX,T13) ) return
    if (.not. copystructonly(VecX,T22) ) return
    if (.not. copystructonly(VecX,T23) ) return
    if (.not. copystructonly(VecX,T33) ) return
    if (.not. copystructonly(VecX,S11) ) return
    if (.not. copystructonly(VecX,S12) ) return
    if (.not. copystructonly(VecX,S13) ) return
    if (.not. copystructonly(VecX,S22) ) return
    if (.not. copystructonly(VecX,S23) ) return
    if (.not. copystructonly(VecX,S33) ) return
    if (.not. copystructonly(VecX,Sm11) ) return
    if (.not. copystructonly(VecX,Sm12) ) return
    if (.not. copystructonly(VecX,Sm13) ) return
    if (.not. copystructonly(VecX,Sm22 ) )return
    if (.not. copystructonly(VecX,Sm23) ) return
    if (.not. copystructonly(VecX,Sm33) ) return
    if (.not. copystructonly(Ufilk,Vfilk) ) return
    if (.not. copystructonly(Ufilk,Wfilk) ) return
   

    !Calcul du coef dynamique
    !Cd = < (Lij Hij>/<Hij Hij>
    !avec
    ! Lij=(u|_i u|_j)^ - u|^_i u|^_j
    ! Hij = delta|**2/12 S|^moinsik S|^jk 

    !Calcul de Hij
    call computeFilter(Wave,2.0_WP*deltaLoc,VecXk,Ufilk,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,VecYk,Vfilk,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,VecZk,Wfilk,filter)
    call computeStressTensor(Ufilk,Vfilk,Wfilk,Wave,S11,S12,S13,S22,S23,S33,res)
    if (.not. computeSplitingSij(S11,S12,S13,S22,S23,S33,Sm11,Sm12,Sm13,Sm22,Sm23,Sm33,'moins')) return
    T11%values = (2.0_WP*deltaLoc)**2 *(Sm11%values*S11%values+Sm12%values*S12%values+Sm13%values*S13%values)
    T12%values = (2.0_WP*deltaLoc)**2 *(Sm11%values*S12%values+Sm12%values*S22%values+Sm13%values*S23%values)
    T13%values = (2.0_WP*deltaLoc)**2 *(Sm11%values*S13%values+Sm12%values*S23%values+Sm13%values*S33%values)
    T22%values = (2.0_WP*deltaLoc)**2 *(Sm12%values*S12%values+Sm22%values*S22%values+Sm23%values*S23%values)
    T23%values = (2.0_WP*deltaLoc)**2 *(Sm12%values*S13%values+Sm22%values*S23%values+Sm23%values*S33%values)
    T33%values = (2.0_WP*deltaLoc)**2 *(Sm13%values*S13%values+Sm23%values*S23%values+Sm33%values*S33%values)
    !Calcul de Lij
    call computeFilter(Wave,2.0_WP*deltaLoc,VecXk,Sm11,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,VecYk,Sm22,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,VecZk,Sm33,filter)
    !call computeT_ijVecA(VecX,VecY,VecZ,Sm11,Sm22,Sm33,Wave,S11,S12,S13,S22,S23,S33,filter,res,2.0_WP*deltaLoc)
    call btran(VecYk,Sm13,res)
    call btran(VecZk,Sm23,res)
    call computeT_ijVecA(VecX,Sm13,Sm23,Sm11,Sm22,Sm33,Wave,S11,S12,S13,S22,S23,S33,filter,res,2.0_WP*deltaLoc)
    Sm11%values = T11%values*S11%values + T22%values*S22%values + T33%values*S33%values + &
                & 2.0_WP*(T12%values*S12%values + T13%values*S13%values+T23%values*S23%values)
    Sm22%values = T11%values**2+T22%values**2+T33%values**2+2.0_WP*(T12%values**2+T13%values**2+T23%values**2) 
    coefLoc = computeFieldAvg(Sm11,spec_rank) / computeFieldAvg(Sm22,spec_rank) 
    if (spec_rank .eq. 0) write (6,'(a,1x,g15.8)') '[INFO LES RGM SikmoinsSkjDyn Fab] Dynanmic coef',coefLoc

    call computeStressTensor(VecXk,VecYk,VecZk,Wave,S11,S12,S13,S22,S23,S33,res)
    if (.not. computeSplitingSij(S11,S12,S13,S22,S23,S33,Sm11,Sm12,Sm13,Sm22,Sm23,Sm33,'moins')) return
    T11%values = coefLoc*deltaLoc**2*(Sm11%values*S11%values+Sm12%values*S12%values+Sm13%values*S13%values)
    T12%values = coefLoc*deltaLoc**2*(Sm11%values*S12%values+Sm12%values*S22%values+Sm13%values*S23%values)
    T13%values = coefLoc*deltaLoc**2*(Sm11%values*S13%values+Sm12%values*S23%values+Sm13%values*S33%values)
    T22%values = coefLoc*deltaLoc**2*(Sm12%values*S12%values+Sm22%values*S22%values+Sm23%values*S23%values)
    T23%values = coefLoc*deltaLoc**2*(Sm12%values*S13%values+Sm22%values*S23%values+Sm23%values*S33%values)
    T33%values = coefLoc*deltaLoc**2*(Sm13%values*S13%values+Sm23%values*S23%values+Sm33%values*S33%values)

    call ftran(T11,Ufilk,res)
    call ftran(T12,Vfilk,res)
    call ftran(T13,Wfilk,res)
    call computeDivergenceK(Wave,Ufilk,Vfilk,Wfilk,T1)
    call ftran(T23,Ufilk,res)
    call ftran(T22,Wfilk,res)
    call computeDivergenceK(Wave,VfilK,Wfilk,Ufilk,T2)
    call ftran(T13,Vfilk,res)
    call ftran(T33,Wfilk,res)
    call computeDivergenceK(Wave,Vfilk,Ufilk,Wfilk,T3)

    if ( pdffield ) then
      call btran(T1,T11,res)
      call btran(T2,T22,res)
      call btran(T3,T33,res)
      T11%name="sgs_RGMSikOjkOikSjk_T1"
      T22%name="sgs_RGMSikOjkOikSjk_T2"
      T33%name="sgs_RGMSikOjkOikSjk_T3"
      if (.not. computeFieldPDF(ite,T11,300,spec_rank,10000)) return
      if (.not. computeFieldPDF(ite,T22,300,spec_rank,10000)) return
      if (.not. computeFieldPDF(ite,T33,300,spec_rank,10000)) return
      call btran(VecXk,Sm12,res)
      call btran(VecYk,Sm13,res)
      call btran(VecZk,Sm23,res)
      T11%values = Sm12%values*T11%values + Sm13%values*T22%values + Sm23%values*T33%values
      coefLoc = - computeFieldAvg(T11,spec_rank)
      if (spec_rank .eq. 0) write (6,'(a,1x,g15.8)') '[INFO LES] RGM SikmoinsSkj dissipation over the box ',coefLoc
    endif

    ! ===== Free memory and delete temporary fields =====
    call deleteDataLayout(T11)
    call deleteDataLayout(T12)
    call deleteDataLayout(T13)
    call deleteDataLayout(T23)
    call deleteDataLayout(T22)
    call deleteDataLayout(T33)
    call deleteDataLayout(S11)
    call deleteDataLayout(S12)
    call deleteDataLayout(S13)
    call deleteDataLayout(S23)
    call deleteDataLayout(S22)
    call deleteDataLayout(S33)
    call deleteDataLayout(Sm11)
    call deleteDataLayout(Sm12)
    call deleteDataLayout(Sm13)
    call deleteDataLayout(Sm23)
    call deleteDataLayout(Sm22)
    call deleteDataLayout(Sm33)
    call deleteDataLayout(Ufilk)
    call deleteDataLayout(Vfilk)
    call deleteDataLayout(Wfilk)

    res = .true.

end subroutine RGMSikmoinsSkjDynFab


!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> Modèle du gradient en prenant en compte l'incompressibilité. Calcul de 3 termes par composantes du
!> tenseur sous maille au lieu de 6. Calcul de d/dxj(tau_ij)=Ti
!> d/dxj(tau_ij) = C delta^2 (duj/dxk d^2ui/dxjdxk) 
!> @param [in] VecX velocity component along first direction in spectral space
!> @param [in] VecY velocity component along second direction in spectral space
!> @param [in] VecZ velocity component along third direction in spectral space
!> @param [in] FieldForParam field to get informations about grid (physical space) 
!> @param [in] Wave wave numbers
!> @param [in] spec_rank rank of processor
!> @param [inout] T1 = d/dxj (Cd delta|**2 S|1k- Sjk) Cd = < (eps_nu + TijSij|^>/<delta**2 ( (Sik-| Sjk| Sij|)^ - (Sik-| Sjk|)^ Sij|^)> 
!> @param [inout] T2 = d/dxj (Cd delta|**2 S|2k- Sjk)
!> @param [inout] T3 = d/dxj (Cd delta|**2 S|3k- Sjk)
!------------------------------------------------------------------------------
subroutine RGMSikmoinsSkjDynYouMoin(T1,T2,T3,VecX,VecY,VecZ,VecXk,VecYk,VecZk,Wave,spec_rank,filter,delta,ite)

  
   TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: VecX,VecY,VecZ
   TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: VecXk,VecYk,VecZk
   TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT) :: T1,T2,T3
   TYPE(WaveNumbers), INTENT(IN) :: Wave
   INTEGER,INTENT(IN)            :: spec_rank
   INTEGER,INTENT(IN)            :: filter 
   REAL(WP),INTENT(IN),OPTIONAL  :: delta
   INTEGER,INTENT(IN),OPTIONAL   :: ite

    ! Local field
    type(REAL_DATA_LAYOUT) :: Sm11,Sm12,Sm13,Sm22,Sm23,Sm33
    type(REAL_DATA_LAYOUT) :: S11,S12,S13,S22,S23,S33
    type(REAL_DATA_LAYOUT) :: H11,H12,H13,H22,H23,H33
    type(REAL_DATA_LAYOUT) :: T11,T12,T13,T22,T23,T33,temp
    type(COMPLEX_DATA_LAYOUT)::Ufilk,Vfilk,Wfilk 
    real(WP) :: coefLoc,deltaLoc,avgtijsijbar,avg_epsNu,nu
    logical :: res,pdffield
    real(wp) :: SikmoinsSjkSijSijtoutchap,SikmoinsSjkSijchapSijchap

    ! ===== Initialisation =====
    res = .false.
    pdffield=.false.
    if (present(ite)) then
      if ( ite .ne. -1 ) pdffield = .true.
    endif
    if (present(delta)) then
      deltaLoc=delta
    else
      deltaLoc=computeDelta(VecX)
    endif 
    CALL parser_read('Viscosity',nu)
    if (.not. copystructonly(VecXk,Ufilk)) return
    if (.not. copystructonly(VecX,T11) ) return
    if (.not. copystructonly(VecX,T12) ) return
    if (.not. copystructonly(VecX,T13) ) return
    if (.not. copystructonly(VecX,T22) ) return
    if (.not. copystructonly(VecX,T23) ) return
    if (.not. copystructonly(VecX,T33) ) return
    if (.not. copystructonly(VecX,S11) ) return
    if (.not. copystructonly(VecX,S12) ) return
    if (.not. copystructonly(VecX,S13) ) return
    if (.not. copystructonly(VecX,S22) ) return
    if (.not. copystructonly(VecX,S23) ) return
    if (.not. copystructonly(VecX,S33) ) return
    if (.not. copystructonly(VecX,temp) ) return
    if (.not. copystructonly(VecX,H11) ) return
    if (.not. copystructonly(VecX,H12) ) return
    if (.not. copystructonly(VecX,H13) ) return
    if (.not. copystructonly(VecX,H22) ) return
    if (.not. copystructonly(VecX,H23) ) return
    if (.not. copystructonly(VecX,H33) ) return
    if (.not. copystructonly(VecX,Sm11) ) return
    if (.not. copystructonly(VecX,Sm12) ) return
    if (.not. copystructonly(VecX,Sm13) ) return
    if (.not. copystructonly(VecX,Sm22 ) )return
    if (.not. copystructonly(VecX,Sm23) ) return
    if (.not. copystructonly(VecX,Sm33) ) return
    if (.not. copystructonly(Ufilk,Vfilk) ) return
    if (.not. copystructonly(Ufilk,Wfilk) ) return
   
    !Calcul du coef dynamique
    !Cd = < (eps_nu + TijSij|^>/<delta**2 ( (Sik-| Sjk| Sij|)^ - (Sik-| Sjk|)^ Sij|^)>
    !avec
    ! Tij=(u|_i u|_j)^ - u|^_i u|^_j
    
    if (.not. computeFieldPDF(1,VecX,300,spec_rank,10000)) return
    if (.not. computeFieldPDF(1,VecY,300,spec_rank,10000)) return
    if (.not. computeFieldPDF(1,VecZ,300,spec_rank,10000)) return

    !Calcul de Hij
    call computeStressTensor(VecXk,VecYk,VecZk,Wave,S11,S12,S13,S22,S23,S33,res)
    call RGMSikmoinsSkj_tij(T11,T12,T13,T22,T23,T33,VecX,VecXk,VecYk,VecZk,Wave,spec_rank,1.0_WP,1.0_WP)
    temp%values = T11%values*S11%values + T22%values*S22%values + T33%values*S33%values + &
        2.0_WP *(T12%values*S12%values + T13%values*S13%values + T23%values*S23%values) !  ->  Sik-| Sjk| Sij|
    call computeFilter(Wave,2.0_WP*deltaLoc,T11,H11,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,T12,H12,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,T13,H13,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,T22,H22,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,T23,H23,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,T33,H33,filter)  ! -> Hij = (S|ik- S|jk)^
    call computeFilter(Wave,2.0_WP*deltaLoc,VecXk,Ufilk,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,VecYk,Vfilk,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,VecZk,Wfilk,filter)
    call computeStressTensor(Ufilk,Vfilk,Wfilk,Wave,Sm11,Sm12,Sm13,Sm22,Sm23,Sm33,res) ! -> Smij = Sij|^
    T11%values = H11%values*Sm11%values + H22%values*Sm22%values + H33%values*Sm33%values + &
        2.0_WP *(H12%values*Sm12%values + H13%values*Sm13%values + H23%values*Sm23%values)    ! T11 -> (S|ik- S|jk)^ Sij|^
    call computeFilter(Wave,2.0_WP*deltaLoc,temp,T22,filter)                                  ! T22 -> (S|ik- S|jk Sij|)^
    SikmoinsSjkSijSijtoutchap=computeFieldAvg(T11,spec_rank)
    SikmoinsSjkSijchapSijchap=computeFieldAvg(T22,spec_rank)
    T33%values = deltaLoc**2 * ( T22%values - T11%values )  !-> delta**2 ( (Sik-| Sjk| Sij|)^ - (Sik-| Sjk|)^ Sij|^)
    !Calcul de TijSij|^
    call computeFilter(Wave,2.0_WP*deltaLoc,VecXk,S11,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,VecYk,S22,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,VecZk,S33,filter)
    call computeT_ijVecA(VecX,VecY,VecZ,S11,S22,S33,Wave,T11  ,T12  ,T13  ,T22  ,T23  ,T33  ,filter,res    ,2.0_WP*deltaLoc)
    Sm11%values = T11%values*Sm11%values + T22%values*Sm22%values + T33%values*Sm33%values + &
        & 2.0_WP*(T12%values*Sm12%values + T13%values*Sm13%values + T23%values*Sm23%values)
    avgtijsijbar = computeFieldAvg(Sm11,spec_rank)
    !Calcul de eps_nu
    call computeViscousDissResolFiltered(Sm11,VecX,VecXk,VecYk,VecZk,wave,nu,filter,2.0_WP*deltaLoc,res)
    avg_epsNu = computeFieldAvg(Sm11,spec_rank)
    coefLoc = ( avgtijsijbar + avg_epsNu ) / ( deltaLoc**2 * (SikmoinsSjkSijchapSijchap - SikmoinsSjkSijSijtoutchap))
    if (spec_rank .eq. 0) write (6,'(a,1x,g15.8)') '[INFO LES RGM SikmoinsSkj DynYouMoin] (S-S)^S^     ',deltaLoc**2*SikmoinsSjkSijchapSijchap
    if (spec_rank .eq. 0) write (6,'(a,1x,g15.8)') '[INFO LES RGM SikmoinsSkj DynYouMoin] (S-SS)^      ',deltaLoc**2*SikmoinsSjkSijSijtoutchap
    if (spec_rank .eq. 0) write (6,'(a,1x,g15.8)') '[INFO LES RGM SikmoinsSkj DynYouMoin] eps_nu       ',avg_epsNu
    if (spec_rank .eq. 0) write (6,'(a,1x,g15.8)') '[INFO LES RGM SikmoinsSkj DynYouMoin] TijSij|^     ',avgtijsijbar
    if (spec_rank .eq. 0) write (6,'(a,1x,g15.8)') '[INFO LES RGM SikmoinsSkj DynYouMoin] denom        ',deltaLoc**2*(SikmoinsSjkSijchapSijchap-SikmoinsSjkSijSijtoutchap)
    if (spec_rank .eq. 0) write (6,'(a,1x,g15.8)') '[INFO LES RGM SikmoinsSkj DynYouMoin] Dynanmic coef',coefLoc
    call RGMSikmoinsSkj_tij(T11,T12,T13,T22,T23,T33,VecX,VecXk,VecYk,VecZk,Wave,spec_rank,deltaLoc,coefLoc)

    call ftran(T11,Ufilk,res)
    call ftran(T12,Vfilk,res)
    call ftran(T13,Wfilk,res)
    call computeDivergenceK(Wave,Ufilk,Vfilk,Wfilk,T1)
    call ftran(T23,Ufilk,res)
    call ftran(T22,Wfilk,res)
    call computeDivergenceK(Wave,VfilK,Wfilk,Ufilk,T2)
    call ftran(T13,Vfilk,res)
    call ftran(T33,Wfilk,res)
    call computeDivergenceK(Wave,Vfilk,Ufilk,Wfilk,T3)

    if ( pdffield ) then
      call btran(T1,Sm11,res)
      call btran(T2,Sm22,res)
      call btran(T3,Sm33,res)
      Sm11%name="sgs_RGMSikSjk_dynYouMoin_T1"
      SM22%name="sgs_RGMSikSjk_dynYouMoin_T2"
      SM33%name="sgs_RGMSikSjk_dynYouMoin_T3"
      if (.not. computeFieldPDF(ite,Sm11,300,spec_rank,10000)) return
      if (.not. computeFieldPDF(ite,Sm22,300,spec_rank,10000)) return
      if (.not. computeFieldPDF(ite,Sm33,300,spec_rank,10000)) return
      call computeStressTensor(VecXk,VecYk,VecZk,Wave,S11,S12,S13,S22,S23,S33,res)
      T11%values =  ( T11%values*S11%values + T22%values*S22%values + T33%values*S33%values + &
          & 2.0_WP* ( T12%values*S12%values + T13%values*S13%values + T23%values*S23%values ))
      T11%name="diss_RGMSikSjk_dynYouMoin"
      if (.not. computeFieldPDF(ite,T11,300,spec_rank,10000)) return
      coefLoc = computeFieldAvg(T11,spec_rank)
      if (spec_rank .eq. 0) write (6,'(a,1x,g15.8)') '[INFO LES] RGM SikmoinsSkj dissipation over the box ',coefLoc
    endif

    ! ===== Free memory and delete temporary fields =====
    call deleteDataLayout(T11)
    call deleteDataLayout(T12)
    call deleteDataLayout(T13)
    call deleteDataLayout(T23)
    call deleteDataLayout(T22)
    call deleteDataLayout(T33)
    call deleteDataLayout(H11)
    call deleteDataLayout(H12)
    call deleteDataLayout(H13)
    call deleteDataLayout(H23)
    call deleteDataLayout(H22)
    call deleteDataLayout(H33)
    call deleteDataLayout(S11)
    call deleteDataLayout(S12)
    call deleteDataLayout(S13)
    call deleteDataLayout(S23)
    call deleteDataLayout(S22)
    call deleteDataLayout(S33)
    call deleteDataLayout(temp)
    call deleteDataLayout(Sm11)
    call deleteDataLayout(Sm12)
    call deleteDataLayout(Sm13)
    call deleteDataLayout(Sm23)
    call deleteDataLayout(Sm22)
    call deleteDataLayout(Sm33)
    call deleteDataLayout(Ufilk)
    call deleteDataLayout(Vfilk)
    call deleteDataLayout(Wfilk)

    res = .true.

end subroutine RGMSikmoinsSkjDynYouMoin 


!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> Modèle du gradient en prenant en compte l'incompressibilité. Calcul de 3 termes par composantes du
!> tenseur sous maille au lieu de 6. Calcul de d/dxj(tau_ij)=Ti
!> d/dxj(tau_ij) = C delta^2 (duj/dxk d^2ui/dxjdxk) 
!> @param [in] VecX velocity component along first direction in spectral space
!> @param [in] VecY velocity component along second direction in spectral space
!> @param [in] VecZ velocity component along third direction in spectral space
!> @param [in] FieldForParam field to get informations about grid (physical space) 
!> @param [in] Wave wave numbers
!> @param [in] spec_rank rank of processor
!> @param [inout] T1 = eps_1jk eps_jmn C delta^2 duk/dxl d^2un/dxmdxl
!> @param [inout] T2 = eps_2jk eps_jmn C delta^2 duk/dxl d^2un/dxmdxl
!> @param [inout] T3 = eps_3jk eps_jmn C delta^2 duk/dxl d^2un/dxmdxl
!------------------------------------------------------------------------------
subroutine RGMSikmoinsSkjDyn(T1,T2,T3,VecX,VecY,VecZ,Wave,spec_rank,filter,delta,ite)

  
   TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: VecX,VecY,VecZ
   TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT) :: T1,T2,T3
   TYPE(WaveNumbers), INTENT(IN) :: Wave
   INTEGER,INTENT(IN)            :: spec_rank
   INTEGER,INTENT(IN)            :: filter 
   REAL(WP),INTENT(IN),OPTIONAL  :: delta
   INTEGER,INTENT(IN),OPTIONAL   :: ite

    ! Local field
    type(REAL_DATA_LAYOUT) :: Sm11,Sm12,Sm13,Sm22,Sm23,Sm33
    type(REAL_DATA_LAYOUT) :: S11,S12,S13,S22,S23,S33
    type(REAL_DATA_LAYOUT) :: H11,H12,H13,H22,H23,H33
    type(REAL_DATA_LAYOUT) :: T11,T12,T13,T22,T23,T33
    type(COMPLEX_DATA_LAYOUT)::Ufilk,Vfilk,Wfilk 
    integer                  :: nbcpus
    real(WP) :: coefLoc,deltaLoc
    logical :: res,pdffield

    ! ===== Initialisation =====
    res = .false.
    pdffield=.false.
    nbcpus = getnbcpus()
    if (present(ite)) then
      if ( ite .ne. -1 ) pdffield = .true.
    endif
    if (present(delta)) then
      deltaLoc=delta
    else
      deltaLoc=computeDelta(VecX)
    endif 
    if (.not. initDataLayout("UFilk", UFilk,(VecX%nx/2)+1,VecX%ny,VecX%nz, &
       & VecX%Lx,VecX%Ly,VecX%Lz,nbcpus,spec_rank,alongZ)) then 
      write(6,'(a,i0,a)')'[ERROR] Not able to create Exk'
      return
    endif
    if (.not. copystructonly(VecX,T11) ) return
    if (.not. copystructonly(VecX,T12) ) return
    if (.not. copystructonly(VecX,T13) ) return
    if (.not. copystructonly(VecX,T22) ) return
    if (.not. copystructonly(VecX,T23) ) return
    if (.not. copystructonly(VecX,T33) ) return
    if (.not. copystructonly(VecX,S11) ) return
    if (.not. copystructonly(VecX,S12) ) return
    if (.not. copystructonly(VecX,S13) ) return
    if (.not. copystructonly(VecX,S22) ) return
    if (.not. copystructonly(VecX,S23) ) return
    if (.not. copystructonly(VecX,S33) ) return
    if (.not. copystructonly(VecX,H11) ) return
    if (.not. copystructonly(VecX,H12) ) return
    if (.not. copystructonly(VecX,H13) ) return
    if (.not. copystructonly(VecX,H22) ) return
    if (.not. copystructonly(VecX,H23) ) return
    if (.not. copystructonly(VecX,H33) ) return
    if (.not. copystructonly(VecX,Sm11) ) return
    if (.not. copystructonly(VecX,Sm12) ) return
    if (.not. copystructonly(VecX,Sm13) ) return
    if (.not. copystructonly(VecX,Sm22 ) )return
    if (.not. copystructonly(VecX,Sm23) ) return
    if (.not. copystructonly(VecX,Sm33) ) return
    if (.not. copystructonly(Ufilk,Vfilk) ) return
    if (.not. copystructonly(Ufilk,Wfilk) ) return
   

    !Calcul du coef dynamique
    !Cd = < (Lij Hij>/<Hij Hij>
    !avec
    ! Lij=(u|_i u|_j)^ - u|^_i u|^_j
    ! Hij = delta|^**2/12 S|^moinsik S|^jk - delta|**2/12 (S|moinsik S|kj)^

    !Calcul de Hij
    call ftran(VecX,Ufilk,res)
    call ftran(VecY,Vfilk,res)
    call ftran(VecZ,Wfilk,res)
    call computeStressTensor(Ufilk,Vfilk,Wfilk,Wave,S11,S12,S13,S22,S23,S33,res)
    if (.not. computeSplitingSij(S11,S12,S13,S22,S23,S33,Sm11,Sm12,Sm13,Sm22,Sm23,Sm33,'moins')) return
    T11%values = Sm11%values*S11%values+Sm12%values*S12%values+Sm13%values*S13%values
    T12%values = Sm11%values*S12%values+Sm12%values*S22%values+Sm13%values*S23%values
    T13%values = Sm11%values*S13%values+Sm12%values*S23%values+Sm13%values*S33%values
    T22%values = Sm12%values*S12%values+Sm22%values*S22%values+Sm23%values*S23%values
    T23%values = Sm12%values*S13%values+Sm22%values*S23%values+Sm23%values*S33%values
    T33%values = Sm13%values*S13%values+Sm23%values*S23%values+Sm33%values*S33%values
    call computeFilter(Wave,2.0_WP*deltaLoc,T11,H11,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,T12,H12,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,T13,H13,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,T22,H22,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,T23,H23,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,T33,H33,filter)

    call computeFilter(Wave,2.0_WP*deltaLoc,VecX,Ufilk,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,VecY,Vfilk,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,VecZ,Wfilk,filter)
    call computeStressTensor(Ufilk,Vfilk,Wfilk,Wave,S11,S12,S13,S22,S23,S33,res)
    if (.not. computeSplitingSij(S11,S12,S13,S22,S23,S33,Sm11,Sm12,Sm13,Sm22,Sm23,Sm33,'moins')) return
    T11%values = Sm11%values*S11%values+Sm12%values*S12%values+Sm13%values*S13%values
    T12%values = Sm11%values*S12%values+Sm12%values*S22%values+Sm13%values*S23%values
    T13%values = Sm11%values*S13%values+Sm12%values*S23%values+Sm13%values*S33%values
    T22%values = Sm12%values*S12%values+Sm22%values*S22%values+Sm23%values*S23%values
    T23%values = Sm12%values*S13%values+Sm22%values*S23%values+Sm23%values*S33%values
    T33%values = Sm13%values*S13%values+Sm23%values*S23%values+Sm33%values*S33%values
    H11%values = (2.0_WP*deltaLoc)**2 * T11%values - deltaLoc**2* H11%values
    H12%values = (2.0_WP*deltaLoc)**2 * T12%values - deltaLoc**2* H12%values
    H13%values = (2.0_WP*deltaLoc)**2 * T13%values - deltaLoc**2* H13%values
    H22%values = (2.0_WP*deltaLoc)**2 * T22%values - deltaLoc**2* H22%values
    H23%values = (2.0_WP*deltaLoc)**2 * T23%values - deltaLoc**2* H23%values
    H33%values = (2.0_WP*deltaLoc)**2 * T33%values - deltaLoc**2* H33%values

    !Calcul de Lij
    call computeFilter(Wave,2.0_WP*deltaLoc,VecX,Sm11,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,VecY,Sm22,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,VecZ,Sm33,filter)
!        computeT_ijVecA(in1 ,in2 ,in3 ,in1f,in2f,in3f,wave,out11,out12,out13,out22,out23,out33,filter,success,delta)
    call computeT_ijVecA(VecX,VecY,VecZ,Sm11,Sm22,Sm33,Wave,T11  ,T12  ,T13  ,T22  ,T23  ,T33  ,filter,res    ,2.0_WP*deltaLoc)

    Sm11%values = T11%values*H11%values + T22%values*H22%values + T33%values*H33%values + &
                & 2.0_WP*(T12%values*H12%values + T13%values*H13%values+T23%values*H23%values)
    Sm22%values = H11%values**2+H22%values**2+H33%values**2+2.0_WP*(H12%values**2+H13%values**2+H23%values**2) 
    coefLoc = computeFieldAvg(Sm11,spec_rank) / computeFieldAvg(Sm22,spec_rank) 
    if (spec_rank .eq. 0) write (6,'(a,1x,g15.8)') '[INFO LES RGM SikmoinsSkjDyn] Dynanmic coef',coefLoc

    call ftran(VecX,Ufilk,res)
    call ftran(VecY,Vfilk,res)
    call ftran(VecZ,Wfilk,res)
                                         !out11,out12,out13,out22,out23,out33
    call computeStressTensor(Ufilk,Vfilk,Wfilk,Wave,S11,S12,S13,S22,S23,S33,res)
    !if (.not. computeSplitingSij_analyt_field(S11in,S12in,S13in,S22in,S23in,S33in,S11out,S12out,S13out,S22out,S23out,S33out,'moins')) return
    if (.not. computeSplitingSij(S11,S12,S13,S22,S23,S33,Sm11,Sm12,Sm13,Sm22,Sm23,Sm33,'moins')) return
    T11%values = coefLoc*deltaLoc**2*(Sm11%values*S11%values+Sm12%values*S12%values+Sm13%values*S13%values)
    T12%values = coefLoc*deltaLoc**2*(Sm11%values*S12%values+Sm12%values*S22%values+Sm13%values*S23%values)
    T13%values = coefLoc*deltaLoc**2*(Sm11%values*S13%values+Sm12%values*S23%values+Sm13%values*S33%values)
    T22%values = coefLoc*deltaLoc**2*(Sm12%values*S12%values+Sm22%values*S22%values+Sm23%values*S23%values)
    T23%values = coefLoc*deltaLoc**2*(Sm12%values*S13%values+Sm22%values*S23%values+Sm23%values*S33%values)
    T33%values = coefLoc*deltaLoc**2*(Sm13%values*S13%values+Sm23%values*S23%values+Sm33%values*S33%values)

!    call ComputeDivSymTijVec1Vec2(T1,T2,T3,T12,T13,T23,T11,T22,T33,wave,res)
    call ftran(T11,Ufilk,res)
    call ftran(T12,Vfilk,res)
    call ftran(T13,Wfilk,res)
    call computeDivergenceK(Wave,Ufilk,Vfilk,Wfilk,T1)
    call ftran(T23,Ufilk,res)
    call ftran(T22,Wfilk,res)
    call computeDivergenceK(Wave,VfilK,Wfilk,Ufilk,T2)
    call ftran(T13,Vfilk,res)
    call ftran(T33,Wfilk,res)
    call computeDivergenceK(Wave,Vfilk,Ufilk,Wfilk,T3)

    if ( pdffield ) then
      call btran(T1,T11,res)
      call btran(T2,T22,res)
      call btran(T3,T33,res)
      T11%name="sgs_RGMSikOjkOikSjk_dyn_T1"
      T22%name="sgs_RGMSikOjkOikSjk_dyn_T2"
      T33%name="sgs_RGMSikOjkOikSjk_dyn_T3"
      if (.not. computeFieldPDF(ite,T11,300,spec_rank,10000)) return
      if (.not. computeFieldPDF(ite,T22,300,spec_rank,10000)) return
      if (.not. computeFieldPDF(ite,T33,300,spec_rank,10000)) return
    endif

    ! ===== Free memory and delete temporary fields =====
    call deleteDataLayout(T11)
    call deleteDataLayout(T12)
    call deleteDataLayout(T13)
    call deleteDataLayout(T23)
    call deleteDataLayout(T22)
    call deleteDataLayout(T33)
    call deleteDataLayout(H11)
    call deleteDataLayout(H12)
    call deleteDataLayout(H13)
    call deleteDataLayout(H23)
    call deleteDataLayout(H22)
    call deleteDataLayout(H33)
    call deleteDataLayout(S11)
    call deleteDataLayout(S12)
    call deleteDataLayout(S13)
    call deleteDataLayout(S23)
    call deleteDataLayout(S22)
    call deleteDataLayout(S33)
    call deleteDataLayout(Sm11)
    call deleteDataLayout(Sm12)
    call deleteDataLayout(Sm13)
    call deleteDataLayout(Sm23)
    call deleteDataLayout(Sm22)
    call deleteDataLayout(Sm33)
    call deleteDataLayout(Ufilk)
    call deleteDataLayout(Vfilk)
    call deleteDataLayout(Wfilk)

    res = .true.

end subroutine RGMSikmoinsSkjDyn 


!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> @param [in] VecX  velocity component along first direction in spectral space
!> @param [in] VecXk velocity component along first direction in spectral space
!> @param [in] VecYk velocity component along second direction in spectral space
!> @param [in] VecZk velocity component along third direction in spectral space
!> @param [in] 
!> @param [in] Wave wave numbers
!> @param [in] spec_rank rank of processor
!> @param [inout] 
!> @param [inout] 
!> @param [inout] 
!------------------------------------------------------------------------------
subroutine RGMSikmoinsSkj_div(T1,T2,T3,VecX,VecXk,VecYk,VecZk,Wave,spec_rank,delta,coef,ite)

    !I/O data
    type(real_data_layout),intent(in)       :: VecX
    type(complex_data_layout),intent(in)    :: VecXk,VecYk,VecZk
    type(complex_data_layout),intent(inout) :: T1,T2,T3
    type(WaveNumbers),intent(in)            :: Wave
    integer,intent(in)                      :: spec_rank
    real(WP),intent(in),optional            :: coef,delta
    integer,intent(in),optional             :: ite

    !Local data
    type(real_data_layout)  :: T11,T12,T13,T22,T23,T33
    type(real_data_layout)  :: F11,F12,F13,F22,F23,F33
    integer                 :: nbcpus
    logical                 :: pdffield,res
    real(WP)                :: coefLoc,deltaLoc

    pdffield= .false.
    nbcpus = getnbcpus()
    if (present(ite)) then
      if ( ite .ne. -1 ) pdffield = .true.
    endif
    if (present(coef)) then
      coefLoc = coef
    else
      coefLoc = 1.0_WP/12.0_WP
    endif
    if (present(delta)) then
      deltaLoc = delta
    else
      deltaLoc = computeDelta(VecX)
    endif

    if (.not. copystructonly(VecX,T11) ) return
    if (.not. copystructonly(VecX,T12) ) return
    if (.not. copystructonly(VecX,T13) ) return
    if (.not. copystructonly(VecX,T22) ) return
    if (.not. copystructonly(VecX,T23) ) return
    if (.not. copystructonly(VecX,T33) ) return
    if (spec_rank .eq. 0) write (6,'(a)') '[INFO LES] sgs model for velocity: RGMSikmoinsSkj'

    call RGMSikmoinsSkj_tij(T11,T12,T13,T22,T23,T33,VecX,VecXk,VecYk,VecZk,Wave,spec_rank,deltaLoc,coefLoc)
    call ComputeDivSymTijVec1Vec2(T1,T2,T3,T12,T13,T23,T11,T22,T33,Wave,res)
    if (.not. res) write(6,'(a)') '[ERROR] in RGMSikmoinsSkj_div'
    if ( pdffield ) then
      if (.not. copystructonly(VecX,F11) ) return
      if (.not. copystructonly(VecX,F12) ) return
      if (.not. copystructonly(VecX,F13) ) return
      if (.not. copystructonly(VecX,F22) ) return
      if (.not. copystructonly(VecX,F23) ) return
      if (.not. copystructonly(VecX,F33) ) return
      call computeStressTensor(VecXk,VecYk,VecZk,Wave,F11,F12,F13,F22,F23,F33,res)
      T11%values = T11%values*F11%values + T22%values*F22%values + T33%values*F33%values + &
         2.0_WP * (T12%values*F12%values + T13%values*F13%values + T23%values*F23%values )
      T11%name="diss_tijSij"
      if (.not. compute_spectrum(T11,ite,spec_rank)) return
      if (.not. computeFieldPDF(ite,T11,300,spec_rank,10000)) return
      call btran(VecYk,F12,res)
      call btran(VecZk,F13,res)
      F11%values = VecX%values *T11%values
      F22%values = VecX%values *T12%values
      F33%values = VecX%values *T13%values
      if (.not. computeFieldDivergence_scaField(Wave,F11,F22,F33,T11,nbcpus,spec_rank)) return
      F11%values = F12%values *T12%values
      F22%values = F12%values *T22%values
      F33%values = F12%values *T23%values
      if (.not. computeFieldDivergence_scaField(Wave,F11,F22,F33,T22,nbcpus,spec_rank)) return
      F11%values = F13%values *T13%values
      F22%values = F13%values *T23%values
      F33%values = F13%values *T33%values
      if (.not. computeFieldDivergence_scaField(Wave,F11,F22,F33,T33,nbcpus,spec_rank)) return
      T11%values = -T11%values - T22%values - T33%values
      T11%name="diss_-duitij_dxj"
      if (.not. compute_spectrum(T11,ite,spec_rank))return
      if (.not. computeFieldPDF(ite,T11,300,spec_rank,10000)) return
      call btran(T1,T11,res)
      call btran(T2,T22,res)
      call btran(T3,T33,res)
!      T11%name="sgs_RGMSikmoinsSjk_T1"
!      T22%name="sgs_RGMSikmoinsSjk_T2"
!      T33%name="sgs_RGMSikmoinsSjk_T3"
!      if (.not. computeFieldPDF(ite,T11,300,spec_rank,10000)) return
!      if (.not. computeFieldPDF(ite,T22,300,spec_rank,10000)) return
!      if (.not. computeFieldPDF(ite,T33,300,spec_rank,10000)) return
      F13%name="diss_-uidtij_dxj"
      F13%values = - (T11%values*VecX%values + T22%values*F12%values + T33%values*F13%values)
      if (.not. compute_spectrum(F13,ite,spec_rank)) return
      if (.not. computeFieldPDF(ite,F13,300,spec_rank,10000)) return
      if (spec_rank .eq. 0) write (6,'(a,1x,g15.8)') '[INFO LES] RGM SikmoinsSkj dissipation over the box ',coefLoc
      call deleteDataLayout(F11)
      call deleteDataLayout(F12)
      call deleteDataLayout(F13)
      call deleteDataLayout(F23)
      call deleteDataLayout(F22)
      call deleteDataLayout(F33)
    endif

    call deleteDataLayout(T11)
    call deleteDataLayout(T12)
    call deleteDataLayout(T13)
    call deleteDataLayout(T23)
    call deleteDataLayout(T22)
    call deleteDataLayout(T33)

end subroutine RGMSikmoinsSkj_div

!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> @param [in] VecX  velocity component along first direction in spectral space
!> @param [in] VecXk velocity component along first direction in spectral space
!> @param [in] VecYk velocity component along second direction in spectral space
!> @param [in] VecZk velocity component along third direction in spectral space
!> @param [in] 
!> @param [in] Wave wave numbers
!> @param [in] spec_rank rank of processor
!> @param [inout] 
!> @param [inout] 
!> @param [inout] 
!------------------------------------------------------------------------------
subroutine RGMSikmoinsSkj_tij(T11,T12,T13,T22,T23,T33,VecX,VecXk,VecYk,VecZk,Wave,spec_rank,delta,coef)

  
    TYPE(REAL_DATA_LAYOUT), INTENT(IN)     :: VecX
    TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN)  :: VecXk,VecYk,VecZk
    type(REAL_DATA_LAYOUT),intent(inout)   :: T11,T12,T13,T22,T23,T33
    TYPE(WaveNumbers), INTENT(IN)          :: Wave
    INTEGER,INTENT(IN)                     :: spec_rank
    REAL(WP),INTENT(IN)                    :: delta,coef

    ! Local field
    type(REAL_DATA_LAYOUT) :: Sm11,Sm12,Sm13,Sm22,Sm23,Sm33
    type(REAL_DATA_LAYOUT) :: S11,S12,S13,S22,S23,S33
    integer                :: i,j,k
    logical                :: res,positive

    ! ===== Initialisation =====
    res = .false.
    if (.not. copystructonly(VecX,S11) ) return
    if (.not. copystructonly(VecX,S12) ) return
    if (.not. copystructonly(VecX,S13) ) return
    if (.not. copystructonly(VecX,S22) ) return
    if (.not. copystructonly(VecX,S23) ) return
    if (.not. copystructonly(VecX,S33) ) return
    if (.not. copystructonly(VecX,Sm11) ) return
    if (.not. copystructonly(VecX,Sm12) ) return
    if (.not. copystructonly(VecX,Sm13) ) return
    if (.not. copystructonly(VecX,Sm22 ) )return
    if (.not. copystructonly(VecX,Sm23) ) return
    if (.not. copystructonly(VecX,Sm33) ) return
    
    call computeStressTensor(VecXk,VecYk,VecZk,Wave,S11,S12,S13,S22,S23,S33,res)
    if (.not. computeSplitingSij(S11,S12,S13,S22,S23,S33,Sm11,Sm12,Sm13,Sm22,Sm23,Sm33,'moins')) return
    T11%values = coef*delta**2*(Sm11%values*S11%values+Sm12%values*S12%values+Sm13%values*S13%values)
    T12%values = coef*delta**2*(Sm11%values*S12%values+Sm12%values*S22%values+Sm13%values*S23%values)
    T13%values = coef*delta**2*(Sm11%values*S13%values+Sm12%values*S23%values+Sm13%values*S33%values)
    T22%values = coef*delta**2*(Sm12%values*S12%values+Sm22%values*S22%values+Sm23%values*S23%values)
    T23%values = coef*delta**2*(Sm12%values*S13%values+Sm22%values*S23%values+Sm23%values*S33%values)
    T33%values = coef*delta**2*(Sm13%values*S13%values+Sm23%values*S23%values+Sm33%values*S33%values)

    Sm11%values = T11%values * S11%values + T22%values * S22%values + T33%values * S33%values + &
           & 2.0_WP * (T12%values*S12%values + S13%values *T13%values + T23%values * S23%values)
    positive = .false.
    do k=VecX%zmin,VecX%zmax
      do j=VecX%ymin,VecX%ymax
        do i=VecX%xmin,VecX%xmax
          if ( Sm11%values(i,j,k) .gt. -1e-20 ) then
            positive = .true.
          endif
        enddo
      enddo
    enddo
    if (positive) then
      if (spec_rank .eq. 0 ) print *,'ERROR dissipation positive !'
    endif

    ! ===== Free memory and delete temporary fields =====
    call deleteDataLayout(S11)
    call deleteDataLayout(S12)
    call deleteDataLayout(S13)
    call deleteDataLayout(S23)
    call deleteDataLayout(S22)
    call deleteDataLayout(S33)
    call deleteDataLayout(Sm11)
    call deleteDataLayout(Sm12)
    call deleteDataLayout(Sm13)
    call deleteDataLayout(Sm23)
    call deleteDataLayout(Sm22)
    call deleteDataLayout(Sm33)

    res = .true.

end subroutine RGMSikmoinsSkj_tij 


!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> Modèle du gradient en prenant en compte l'incompressibilité. Calcul de 3 termes par composantes du
!> tenseur sous maille au lieu de 6. Calcul de d/dxj(tau_ij)=Ti
!> d/dxj(tau_ij) = C delta^2 (duj/dxk d^2ui/dxjdxk) 
!> @param [in] VecX velocity component along first direction in spectral space
!> @param [in] VecY velocity component along second direction in spectral space
!> @param [in] VecZ velocity component along third direction in spectral space
!> @param [in] FieldForParam field to get informations about grid (physical space) 
!> @param [in] Wave wave numbers
!> @param [in] spec_rank rank of processor
!> @param [inout] T1 = eps_1jk eps_jmn C delta^2 duk/dxl d^2un/dxmdxl
!> @param [inout] T2 = eps_2jk eps_jmn C delta^2 duk/dxl d^2un/dxmdxl
!> @param [inout] T3 = eps_3jk eps_jmn C delta^2 duk/dxl d^2un/dxmdxl
!------------------------------------------------------------------------------
subroutine RGMSikmoinsSkjOikOjkDyn2Fab(T1,T2,T3,VecX,VecY,VecZ,Wave,spec_rank,filter,delta,coef,ite)

  
   TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: VecX,VecY,VecZ
   TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT) :: T1,T2,T3
   TYPE(WaveNumbers), INTENT(IN) :: Wave
   INTEGER,INTENT(IN)            :: spec_rank,filter
   REAL(WP),INTENT(IN),OPTIONAL  :: delta,coef
   INTEGER,INTENT(IN),OPTIONAL   :: ite

    ! Local field
    type(REAL_DATA_LAYOUT) :: Sm11,Sm12,Sm13,Sm22,Sm23,Sm33
    type(REAL_DATA_LAYOUT) :: S11,S12,S13,S22,S23,S33
    type(REAL_DATA_LAYOUT) :: T11,T12,T13,T22,T23,T33
    type(REAL_DATA_LAYOUT) :: H11,H12,H13,H22,H23,H33
    type(REAL_DATA_LAYOUT) :: O23,O31,O12
    type(COMPLEX_DATA_LAYOUT)::Ufilk,Vfilk,Wfilk 
    integer                  :: nbcpus
    real(WP) :: coefLoc,deltaLoc
    logical :: res,pdffield

    ! ===== Initialisation =====
    res = .false.
    pdffield = .false.
    nbcpus = getnbcpus()
    if (present(ite)) then
      if ( ite .ne. -1 ) pdffield = .true.
    endif
    if (present(coef)) then
      coefLoc=coef
    else
      coefLoc=1.0_WP/12.0_WP
    endif
    if (present(delta)) then
      deltaLoc=delta
    else
      deltaLoc=computeDelta(VecX)
    endif 
    if (.not. initDataLayout("UFilk", UFilk,(VecX%nx/2)+1,VecX%ny,VecX%nz, &
       & VecX%Lx,VecX%Ly,VecX%Lz,nbcpus,spec_rank,alongZ)) then 
      write(6,'(a,i0,a)')'[ERROR] Not able to create Exk'
      return
    endif
    if (.not. copystructonly(VecX,T11) ) return
    if (.not. copystructonly(VecX,T12) ) return
    if (.not. copystructonly(VecX,T13) ) return
    if (.not. copystructonly(VecX,T22) ) return
    if (.not. copystructonly(VecX,T23) ) return
    if (.not. copystructonly(VecX,T33) ) return
    if (.not. copystructonly(VecX,H11) ) return
    if (.not. copystructonly(VecX,H12) ) return
    if (.not. copystructonly(VecX,H13) ) return
    if (.not. copystructonly(VecX,H22) ) return
    if (.not. copystructonly(VecX,H23) ) return
    if (.not. copystructonly(VecX,H33) ) return
    if (.not. copystructonly(VecX,S11) ) return
    if (.not. copystructonly(VecX,S12) ) return
    if (.not. copystructonly(VecX,S13) ) return
    if (.not. copystructonly(VecX,S22) ) return
    if (.not. copystructonly(VecX,S23) ) return
    if (.not. copystructonly(VecX,S33) ) return
    if (.not. copystructonly(VecX,Sm11) ) return
    if (.not. copystructonly(VecX,Sm12) ) return
    if (.not. copystructonly(VecX,Sm13) ) return
    if (.not. copystructonly(VecX,Sm22) )return
    if (.not. copystructonly(VecX,Sm23) ) return
    if (.not. copystructonly(VecX,Sm33) ) return
    if (.not. copystructonly(VecX,O23) )return
    if (.not. copystructonly(VecX,O31) ) return
    if (.not. copystructonly(VecX,O12) ) return
    if (.not. copystructonly(Ufilk,Vfilk) ) return
    if (.not. copystructonly(Ufilk,Wfilk) ) return
    

    !Cd = < (Lij-Hij)Tij>/<Tij Tij>
    !avec
    ! Lij=(u|_i u|_j)^ - u|^_i u|^_j
    ! Hij = delta|^**2/12 S|^moinsik S|^jk
    ! Tij = delta|^**2 O|^ik O|^jk
    
    !Calcul de Hij
    call computeFilter(Wave,2.0_WP*deltaLoc,VecX,Ufilk,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,VecY,Vfilk,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,VecZ,Wfilk,filter)
    call computeStressTensor(Ufilk,Vfilk,Wfilk,Wave,S11,S12,S13,S22,S23,S33,res)
    if (.not. computeSplitingSij(S11,S12,S13,S22,S23,S33,Sm11,Sm12,Sm13,Sm22,Sm23,Sm33,'moins')) return
    H11%values = (2.0_WP*deltaLoc)**2*(Sm11%values*S11%values+Sm12%values*S12%values+Sm13%values*S13%values)
    H12%values = (2.0_WP*deltaLoc)**2*(Sm11%values*S12%values+Sm12%values*S22%values+Sm13%values*S23%values)
    H13%values = (2.0_WP*deltaLoc)**2*(Sm11%values*S13%values+Sm12%values*S23%values+Sm13%values*S33%values)
    H22%values = (2.0_WP*deltaLoc)**2*(Sm12%values*S12%values+Sm22%values*S22%values+Sm23%values*S23%values)
    H23%values = (2.0_WP*deltaLoc)**2*(Sm12%values*S13%values+Sm22%values*S23%values+Sm23%values*S33%values)
    H33%values = (2.0_WP*deltaLoc)**2*(Sm13%values*S13%values+Sm23%values*S23%values+Sm33%values*S33%values)

    !Calcul de Mij
    ! Tij = delta|^**2 O|^ik O|^jk
    call computeVorticity(Ufilk,Vfilk,Wfilk,Wave,O23,O31,O12,res)
    T11%values =   (2.0_WP*deltaLoc)**2/12.0_WP * (O12%values*O12%values + O31%values*O31%values)
    T12%values = - (2.0_WP*deltaLoc)**2/12.0_WP * O31%values*O23%values
    T13%values = - (2.0_WP*deltaLoc)**2/12.0_WP * O12%values*O23%values
    T22%values =   (2.0_WP*deltaLoc)**2/12.0_WP * (O12%values*O12%values + O23%values*O23%values)
    T23%values = - (2.0_WP*deltaLoc)**2/12.0_WP * O12%values*O31%values
    T33%values =   (2.0_WP*deltaLoc)**2/12.0_WP * (O31%values*O31%values + O23%values*O23%values)

    !Calcul de Lij
    call computeFilter(Wave,2.0_WP*deltaLoc,VecX,Sm11,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,VecY,Sm22,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,VecZ,Sm33,filter)
    call computeT_ijVecA(VecX,VecY,VecZ,Sm11,Sm22,Sm33,Wave,S11,S12,S13,S22,S23,S33,filter,res,2.0_WP*deltaLoc)
    Sm11%values = (S11%values-T11%values)*H11%values + (S22%values-T22%values)*H22%values + (S33%values-T33%values)*H33%values + &
                & 2.0_WP*((S12%values-T12%values)*H12%values + (S13%values-T13%values)*H13%values +(S23%values-T23%values)*H23%values)
    Sm22%values = H11%values**2+H22%values**2+H33%values**2+2.0_WP*(H12%values**2+H13%values**2+H23%values**2) 
    coefLoc = computeFieldAvg(Sm11,spec_rank) / computeFieldAvg(Sm22,spec_rank) 
    if (spec_rank .eq. 0) write (6,'(a,1x,g15.8)') '[INFO LES RGMSikmoinsSkj_OikOjk Dyn 2 Fab] Dynanmic coef',coefLoc

    call ftran(VecX,Ufilk,res)
    call ftran(VecY,Vfilk,res)
    call ftran(VecZ,Wfilk,res)
                                         !out11,out12,out13,out22,out23,out33
    call computeStressTensor(Ufilk,Vfilk,Wfilk,Wave,S11,S12,S13,S22,S23,S33,res)
    if (.not. computeSplitingSij(S11,S12,S13,S22,S23,S33,Sm11,Sm12,Sm13,Sm22,Sm23,Sm33,'moins')) return
    call computeVorticity(Ufilk,Vfilk,Wfilk,Wave,O23,O31,O12,res)
    T11%values = coefLoc*deltaLoc**2*(Sm11%values*S11%values+Sm12%values*S12%values+Sm13%values*S13%values + O12%values*O12%values + O31%values*O31%values )
    T12%values = coefLoc*deltaLoc**2*(Sm11%values*S12%values+Sm12%values*S22%values+Sm13%values*S23%values - O31%values*O23%values )
    T13%values = coefLoc*deltaLoc**2*(Sm11%values*S13%values+Sm12%values*S23%values+Sm13%values*S33%values - O12%values*O23%values )
    T22%values = coefLoc*deltaLoc**2*(Sm12%values*S12%values+Sm22%values*S22%values+Sm23%values*S23%values + O12%values*O12%values + O23%values*O23%values )
    T23%values = coefLoc*deltaLoc**2*(Sm12%values*S13%values+Sm22%values*S23%values+Sm23%values*S33%values - O12%values*O31%values )
    T33%values = coefLoc*deltaLoc**2*(Sm13%values*S13%values+Sm23%values*S23%values+Sm33%values*S33%values + O31%values*O31%values + O23%values*O23%values )
    call ftran(T11,Ufilk,res)
    call ftran(T12,Vfilk,res)
    call ftran(T13,Wfilk,res)
    call computeDivergenceK(Wave,Ufilk,Vfilk,Wfilk,T1)
    call ftran(T22,Ufilk,res)
    call ftran(T23,Wfilk,res)
    call computeDivergenceK(Wave,UfilK,Vfilk,Wfilk,T2)
    call ftran(T13,Vfilk,res)
    call ftran(T33,Ufilk,res)
    call computeDivergenceK(Wave,Vfilk,Ufilk,Wfilk,T3)

    if ( pdffield ) then
      call btran(T1,T11,res)
      call btran(T2,T22,res)
      call btran(T3,T33,res)
      T11%name="sgs_RGMSikmoinsSkjOikOjkDynFab_T1"
      T22%name="sgs_RGMSikmoinsSkjOikOjkDynFab_T2"
      T33%name="sgs_RGMSikmoinsSkjOikOjkDynFab_T3"
      if (.not. computeFieldPDF(ite,T11,300,spec_rank,10000)) return
      if (.not. computeFieldPDF(ite,T22,300,spec_rank,10000)) return
      if (.not. computeFieldPDF(ite,T33,300,spec_rank,10000)) return
    endif

    ! ===== Free memory and delete temporary fields =====
    call deleteDataLayout(T11)
    call deleteDataLayout(T12)
    call deleteDataLayout(T13)
    call deleteDataLayout(T23)
    call deleteDataLayout(T22)
    call deleteDataLayout(T33)
    call deleteDataLayout(S11)
    call deleteDataLayout(S12)
    call deleteDataLayout(S13)
    call deleteDataLayout(S23)
    call deleteDataLayout(S22)
    call deleteDataLayout(S33)
    call deleteDataLayout(H11)
    call deleteDataLayout(H12)
    call deleteDataLayout(H13)
    call deleteDataLayout(H23)
    call deleteDataLayout(H22)
    call deleteDataLayout(H33)
    call deleteDataLayout(Sm11)
    call deleteDataLayout(Sm12)
    call deleteDataLayout(Sm13)
    call deleteDataLayout(Sm23)
    call deleteDataLayout(Sm22)
    call deleteDataLayout(Sm33)
    call deleteDataLayout(O23)
    call deleteDataLayout(O31)
    call deleteDataLayout(O12)
    call deleteDataLayout(Ufilk)
    call deleteDataLayout(Vfilk)
    call deleteDataLayout(Wfilk)

    res = .true.

end subroutine RGMSikmoinsSkjOikOjkDyn2Fab


!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> Modèle du gradient en prenant en compte l'incompressibilité. Calcul de 3 termes par composantes du
!> tenseur sous maille au lieu de 6. Calcul de d/dxj(tau_ij)=Ti
!> d/dxj(tau_ij) = C delta^2 (duj/dxk d^2ui/dxjdxk) 
!> @param [in] VecX velocity component along first direction in spectral space
!> @param [in] VecY velocity component along second direction in spectral space
!> @param [in] VecZ velocity component along third direction in spectral space
!> @param [in] FieldForParam field to get informations about grid (physical space) 
!> @param [in] Wave wave numbers
!> @param [in] spec_rank rank of processor
!> @param [inout] T1 = eps_1jk eps_jmn C delta^2 duk/dxl d^2un/dxmdxl
!> @param [inout] T2 = eps_2jk eps_jmn C delta^2 duk/dxl d^2un/dxmdxl
!> @param [inout] T3 = eps_3jk eps_jmn C delta^2 duk/dxl d^2un/dxmdxl
!------------------------------------------------------------------------------
subroutine RGMSikmoinsSkjOikOjkDynFab(T1,T2,T3,VecX,VecY,VecZ,Wave,spec_rank,filter,delta,coef,ite)

  
   TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: VecX,VecY,VecZ
   TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT) :: T1,T2,T3
   TYPE(WaveNumbers), INTENT(IN) :: Wave
   INTEGER,INTENT(IN)            :: spec_rank,filter
   REAL(WP),INTENT(IN),OPTIONAL  :: delta,coef
   INTEGER,INTENT(IN),OPTIONAL   :: ite

    ! Local field
    type(REAL_DATA_LAYOUT) :: Sm11,Sm12,Sm13,Sm22,Sm23,Sm33
    type(REAL_DATA_LAYOUT) :: S11,S12,S13,S22,S23,S33
    type(REAL_DATA_LAYOUT) :: T11,T12,T13,T22,T23,T33
    type(REAL_DATA_LAYOUT) :: H11,H12,H13,H22,H23,H33
    type(REAL_DATA_LAYOUT) :: O23,O31,O12
    type(COMPLEX_DATA_LAYOUT)::Ufilk,Vfilk,Wfilk 
    integer                  :: nbcpus
    real(WP) :: coefLoc,deltaLoc
    logical :: res,pdffield

    ! ===== Initialisation =====
    res = .false.
    pdffield = .false.
    nbcpus = getnbcpus()
    if (present(ite)) then
      if ( ite .ne. -1 ) pdffield = .true.
    endif
    if (present(coef)) then
      coefLoc=coef
    else
      coefLoc=1.0_WP/12.0_WP
    endif
    if (present(delta)) then
      deltaLoc=delta
    else
      deltaLoc=computeDelta(VecX)
    endif 
    if (.not. initDataLayout("UFilk", UFilk,(VecX%nx/2)+1,VecX%ny,VecX%nz, &
       & VecX%Lx,VecX%Ly,VecX%Lz,nbcpus,spec_rank,alongZ)) then 
      write(6,'(a,i0,a)')'[ERROR] Not able to create Exk'
      return
    endif
    if (.not. copystructonly(VecX,T11) ) return
    if (.not. copystructonly(VecX,T12) ) return
    if (.not. copystructonly(VecX,T13) ) return
    if (.not. copystructonly(VecX,T22) ) return
    if (.not. copystructonly(VecX,T23) ) return
    if (.not. copystructonly(VecX,T33) ) return
    if (.not. copystructonly(VecX,H11) ) return
    if (.not. copystructonly(VecX,H12) ) return
    if (.not. copystructonly(VecX,H13) ) return
    if (.not. copystructonly(VecX,H22) ) return
    if (.not. copystructonly(VecX,H23) ) return
    if (.not. copystructonly(VecX,H33) ) return
    if (.not. copystructonly(VecX,S11) ) return
    if (.not. copystructonly(VecX,S12) ) return
    if (.not. copystructonly(VecX,S13) ) return
    if (.not. copystructonly(VecX,S22) ) return
    if (.not. copystructonly(VecX,S23) ) return
    if (.not. copystructonly(VecX,S33) ) return
    if (.not. copystructonly(VecX,Sm11) ) return
    if (.not. copystructonly(VecX,Sm12) ) return
    if (.not. copystructonly(VecX,Sm13) ) return
    if (.not. copystructonly(VecX,Sm22) )return
    if (.not. copystructonly(VecX,Sm23) ) return
    if (.not. copystructonly(VecX,Sm33) ) return
    if (.not. copystructonly(VecX,O23) )return
    if (.not. copystructonly(VecX,O31) ) return
    if (.not. copystructonly(VecX,O12) ) return
    if (.not. copystructonly(Ufilk,Vfilk) ) return
    if (.not. copystructonly(Ufilk,Wfilk) ) return
    

    !Cd = < (Lij-Hij)Tij>/<Tij Tij>
    !avec
    ! Lij=(u|_i u|_j)^ - u|^_i u|^_j
    ! Hij = delta|^**2/12 S|^moinsik S|^jk
    ! Tij = delta|^**2 O|^ik O|^jk
    
    !Calcul de Hij
    call computeFilter(Wave,2.0_WP*deltaLoc,VecX,Ufilk,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,VecY,Vfilk,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,VecZ,Wfilk,filter)
    call computeStressTensor(Ufilk,Vfilk,Wfilk,Wave,S11,S12,S13,S22,S23,S33,res)
    if (.not. computeSplitingSij(S11,S12,S13,S22,S23,S33,Sm11,Sm12,Sm13,Sm22,Sm23,Sm33,'moins')) return
    H11%values = (2.0_WP*deltaLoc)**2/12.0_WP *(Sm11%values*S11%values+Sm12%values*S12%values+Sm13%values*S13%values)
    H12%values = (2.0_WP*deltaLoc)**2/12.0_WP *(Sm11%values*S12%values+Sm12%values*S22%values+Sm13%values*S23%values)
    H13%values = (2.0_WP*deltaLoc)**2/12.0_WP *(Sm11%values*S13%values+Sm12%values*S23%values+Sm13%values*S33%values)
    H22%values = (2.0_WP*deltaLoc)**2/12.0_WP *(Sm12%values*S12%values+Sm22%values*S22%values+Sm23%values*S23%values)
    H23%values = (2.0_WP*deltaLoc)**2/12.0_WP *(Sm12%values*S13%values+Sm22%values*S23%values+Sm23%values*S33%values)
    H33%values = (2.0_WP*deltaLoc)**2/12.0_WP *(Sm13%values*S13%values+Sm23%values*S23%values+Sm33%values*S33%values)

    !Calcul de Mij
    ! Tij = delta|^**2 O|^ik O|^jk
    call computeVorticity(Ufilk,Vfilk,Wfilk,Wave,O23,O31,O12,res)
    T11%values =   (2.0_WP*deltaLoc)**2 * (O12%values*O12%values + O31%values*O31%values)
    T12%values = - (2.0_WP*deltaLoc)**2 * O31%values*O23%values
    T13%values = - (2.0_WP*deltaLoc)**2 * O12%values*O23%values
    T22%values =   (2.0_WP*deltaLoc)**2 * (O12%values*O12%values + O23%values*O23%values)
    T23%values = - (2.0_WP*deltaLoc)**2 * O12%values*O31%values
    T33%values =   (2.0_WP*deltaLoc)**2 * (O31%values*O31%values + O23%values*O23%values)

    !Calcul de Lij
    call computeFilter(Wave,2.0_WP*deltaLoc,VecX,Sm11,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,VecY,Sm22,filter)
    call computeFilter(Wave,2.0_WP*deltaLoc,VecZ,Sm33,filter)
    call computeT_ijVecA(VecX,VecY,VecZ,Sm11,Sm22,Sm33,Wave,S11,S12,S13,S22,S23,S33,filter,res,2.0_WP*deltaLoc)
    Sm11%values = (S11%values-H11%values)*T11%values + (S22%values-H22%values)*T22%values + (S33%values-H33%values)*T33%values + &
                & 2.0_WP*((S12%values-H12%values)*T12%values + (S13%values-H13%values)*T13%values +(S23%values-H23%values)*T23%values)
    Sm22%values = T11%values**2+T22%values**2+T33%values**2+2.0_WP*(T12%values**2+T13%values**2+T23%values**2) 
    coefLoc = computeFieldAvg(Sm11,spec_rank) / computeFieldAvg(Sm22,spec_rank) 
    if (spec_rank .eq. 0) write (6,'(a,1x,g15.8)') '[INFO LES RGMSikmoinsSkj_OikOjk Dyn Fab] Dynanmic coef',coefLoc

    call ftran(VecX,Ufilk,res)
    call ftran(VecY,Vfilk,res)
    call ftran(VecZ,Wfilk,res)
                                         !out11,out12,out13,out22,out23,out33
    call computeStressTensor(Ufilk,Vfilk,Wfilk,Wave,S11,S12,S13,S22,S23,S33,res)
    if (.not. computeSplitingSij(S11,S12,S13,S22,S23,S33,Sm11,Sm12,Sm13,Sm22,Sm23,Sm33,'moins')) return
    call computeVorticity(Ufilk,Vfilk,Wfilk,Wave,O23,O31,O12,res)
    T11%values = coefLoc*deltaLoc**2*(Sm11%values*S11%values+Sm12%values*S12%values+Sm13%values*S13%values + O12%values*O12%values + O31%values*O31%values )
    T12%values = coefLoc*deltaLoc**2*(Sm11%values*S12%values+Sm12%values*S22%values+Sm13%values*S23%values - O31%values*O23%values )
    T13%values = coefLoc*deltaLoc**2*(Sm11%values*S13%values+Sm12%values*S23%values+Sm13%values*S33%values - O12%values*O23%values )
    T22%values = coefLoc*deltaLoc**2*(Sm12%values*S12%values+Sm22%values*S22%values+Sm23%values*S23%values + O12%values*O12%values + O23%values*O23%values )
    T23%values = coefLoc*deltaLoc**2*(Sm12%values*S13%values+Sm22%values*S23%values+Sm23%values*S33%values - O12%values*O31%values )
    T33%values = coefLoc*deltaLoc**2*(Sm13%values*S13%values+Sm23%values*S23%values+Sm33%values*S33%values + O31%values*O31%values + O23%values*O23%values )
    call ftran(T11,Ufilk,res)
    call ftran(T12,Vfilk,res)
    call ftran(T13,Wfilk,res)
    call computeDivergenceK(Wave,Ufilk,Vfilk,Wfilk,T1)
    call ftran(T22,Ufilk,res)
    call ftran(T23,Wfilk,res)
    call computeDivergenceK(Wave,UfilK,Vfilk,Wfilk,T2)
    call ftran(T13,Vfilk,res)
    call ftran(T33,Ufilk,res)
    call computeDivergenceK(Wave,Vfilk,Ufilk,Wfilk,T3)

    if ( pdffield ) then
      call btran(T1,T11,res)
      call btran(T2,T22,res)
      call btran(T3,T33,res)
      T11%name="sgs_RGMSikmoinsSkjOikOjkDynFab_T1"
      T22%name="sgs_RGMSikmoinsSkjOikOjkDynFab_T2"
      T33%name="sgs_RGMSikmoinsSkjOikOjkDynFab_T3"
      if (.not. computeFieldPDF(ite,T11,300,spec_rank,10000)) return
      if (.not. computeFieldPDF(ite,T22,300,spec_rank,10000)) return
      if (.not. computeFieldPDF(ite,T33,300,spec_rank,10000)) return
    endif

    ! ===== Free memory and delete temporary fields =====
    call deleteDataLayout(T11)
    call deleteDataLayout(T12)
    call deleteDataLayout(T13)
    call deleteDataLayout(T23)
    call deleteDataLayout(T22)
    call deleteDataLayout(T33)
    call deleteDataLayout(S11)
    call deleteDataLayout(S12)
    call deleteDataLayout(S13)
    call deleteDataLayout(S23)
    call deleteDataLayout(S22)
    call deleteDataLayout(S33)
    call deleteDataLayout(H11)
    call deleteDataLayout(H12)
    call deleteDataLayout(H13)
    call deleteDataLayout(H23)
    call deleteDataLayout(H22)
    call deleteDataLayout(H33)
    call deleteDataLayout(Sm11)
    call deleteDataLayout(Sm12)
    call deleteDataLayout(Sm13)
    call deleteDataLayout(Sm23)
    call deleteDataLayout(Sm22)
    call deleteDataLayout(Sm33)
    call deleteDataLayout(O23)
    call deleteDataLayout(O31)
    call deleteDataLayout(O12)
    call deleteDataLayout(Ufilk)
    call deleteDataLayout(Vfilk)
    call deleteDataLayout(Wfilk)

    res = .true.

end subroutine RGMSikmoinsSkjOikOjkDynFab

!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> @param [in] VecX  velocity component along first direction in spectral space
!> @param [in] VecXk velocity component along first direction in spectral space
!> @param [in] VecYk velocity component along second direction in spectral space
!> @param [in] VecZk velocity component along third direction in spectral space
!> @param [in] 
!> @param [in] Wave wave numbers
!> @param [in] spec_rank rank of processor
!> @param [inout] 
!> @param [inout] 
!> @param [inout] 
!------------------------------------------------------------------------------
subroutine RGMSikmoinsSkjOikOjk_div(T1,T2,T3,VecX,VecXk,VecYk,VecZk,Wave,spec_rank,delta,coef1,coef2,ite)

    !I/O data
    type(real_data_layout),intent(in)       :: VecX
    type(complex_data_layout),intent(in)    :: VecXk,VecYk,VecZk
    type(complex_data_layout),intent(inout) :: T1,T2,T3
    type(WaveNumbers),intent(in)            :: Wave
    integer,intent(in)                      :: spec_rank
    real(WP),intent(in),optional            :: coef1,coef2,delta
    integer,intent(in),optional             :: ite

    !Local data
    type(real_data_layout)  :: T11,T12,T13,T22,T23,T33
    type(real_data_layout)  :: F11,F12,F13,F22,F23,F33
    logical                 :: pdffield,res
    real(WP)                :: coefLoc1,coefLoc2,deltaLoc

    pdffield= .false.
    if (present(ite)) then
      if ( ite .ne. -1 ) pdffield = .true.
    endif
    if (present(coef1)) then
      coefLoc1 = coef1
    else
      coefLoc1 = 1.0_WP/12.0_WP
    endif
    if (present(coef2)) then
      coefLoc2 = coef2
    else
      coefLoc2 = 1.0_WP/12.0_WP
    endif
    if (present(delta)) then
      deltaLoc = delta
    else
      deltaLoc = computeDelta(VecX)
    endif

    if (.not. copystructonly(VecX,T11) ) return
    if (.not. copystructonly(VecX,T12) ) return
    if (.not. copystructonly(VecX,T13) ) return
    if (.not. copystructonly(VecX,T22) ) return
    if (.not. copystructonly(VecX,T23) ) return
    if (.not. copystructonly(VecX,T33) ) return
    if (spec_rank .eq. 0) write (6,'(a)') '[INFO LES] sgs model for velocity: RGMSikmoinsSkj'

    call RGMSikmoinsSkjOikOjk_tij(T11,T12,T13,T22,T23,T33,VecX,VecXk,VecYk,VecZk,Wave,spec_rank,deltaLoc,coefLoc1,coefLoc2)
    call ComputeDivSymTijVec1Vec2(T1,T2,T3,T12,T13,T23,T11,T22,T33,Wave,res)
    if (.not. res) write(6,'(a)') '[ERROR] in RGMSikmoinsSkjOikOjk_div'
    if ( pdffield ) then
      if (.not. copystructonly(VecX,F11) ) return
      if (.not. copystructonly(VecX,F12) ) return
      if (.not. copystructonly(VecX,F13) ) return
      if (.not. copystructonly(VecX,F22) ) return
      if (.not. copystructonly(VecX,F23) ) return
      if (.not. copystructonly(VecX,F33) ) return
      call computeStressTensor(VecXk,VecYk,VecZk,Wave,F11,F12,F13,F22,F23,F33,res)
      T11%values = T11%values*F11%values + T22%values*F22%values + T33%values*F33%values + &
         2.0_WP * (T12%values*F12%values + T13%values*F13%values + T23%values*F23%values )
      T11%name="diss_tijSij"
      coefLoc1 = computeFieldAvg(T11,spec_rank)
      if (spec_rank .eq. 0) write (6,'(a,1x,g15.8)') '[INFO LES] RGM SikmoinsSkj+OikOjk dissipation over the box ',coefLoc1
      call btran(T1,T11,res)
      call btran(T2,T22,res)
      call btran(T3,T33,res)
      T11%name="sgs_S1kmoinsSkj+O1kOjk"
      T22%name="sgs_S2kmoinsSkj+O2kOjk"
      T33%name="sgs_S3kmoinsSkj+O3kOjk"
      if (.not. computeFieldPDF(ite,T11,300,spec_rank,10000)) return
      if (.not. computeFieldPDF(ite,T22,300,spec_rank,10000)) return
      if (.not. computeFieldPDF(ite,T33,300,spec_rank,10000)) return
      call deleteDataLayout(F11)
      call deleteDataLayout(F12)
      call deleteDataLayout(F13)
      call deleteDataLayout(F23)
      call deleteDataLayout(F22)
      call deleteDataLayout(F33)
    endif

    call deleteDataLayout(T11)
    call deleteDataLayout(T12)
    call deleteDataLayout(T13)
    call deleteDataLayout(T23)
    call deleteDataLayout(T22)
    call deleteDataLayout(T33)

end subroutine RGMSikmoinsSkjOikOjk_div



!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> @param [in] VecX velocity component along first direction in spectral space
!> @param [in] VecXk velocity component along first direction in spectral space
!> @param [in] VecYk velocity component along second direction in spectral space
!> @param [in] VecZk velocity component along third direction in spectral space
!> @param [in] FieldForParam field to get informations about grid (physical space) 
!> @param [in] Wave wave numbers
!> @param [in] spec_rank rank of processor
!> @param [in] delta is the lenght of filter
!> @param [in] coef is the coefficient in front of the model (optional) 
!> @param [inout] T11 =  delta**2 ( coef1*S1k-S1k + coef2*O1kO1k )
!> @param [inout] T22 =  delta**2 ( coef1*S2k-S2k + coef2*O2kO2k )
!> @param [inout] T33 =  delta**2 ( coef1*S3k-S3k + coef2*O3kO3k )
!> @param [inout] T12 =  delta**2 ( coef1*S1k-S2k + coef2*O1kO2k )
!> @param [inout] T13 =  delta**2 ( coef1*S2k-S3k + coef2*O2kO3k )
!> @param [inout] T23 =  delta**2 ( coef1*S3k-S3k + coef2*O3kO3k )
!------------------------------------------------------------------------------
subroutine RGMSikmoinsSkjOikOjk_tij(T11,T12,T13,T22,T23,T33,VecX,VecXk,VecYk,VecZk,Wave,spec_rank,delta,coef1,coef2)

  
   TYPE(REAL_DATA_LAYOUT), INTENT(IN)    :: VecX
   TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: VecXk,VecYk,VecZk
   TYPE(REAL_DATA_LAYOUT),INTENT(INOUT)  :: T11,T12,T13,T22,T23,T33
   TYPE(WaveNumbers), INTENT(IN)         :: Wave
   INTEGER,INTENT(IN)                    :: spec_rank
   REAL(WP),INTENT(IN)                   :: delta,coef1,coef2

    ! Local field
    type(REAL_DATA_LAYOUT) :: S11,S12,S13,S22,S23,S33
    logical                :: res

    ! ===== Initialisation =====
    res = .false.
    if (.not. copystructonly(VecX,S11) ) return
    if (.not. copystructonly(VecX,S12) ) return
    if (.not. copystructonly(VecX,S13) ) return
    if (.not. copystructonly(VecX,S22) ) return
    if (.not. copystructonly(VecX,S23) ) return
    if (.not. copystructonly(VecX,S33) ) return
    call RGMSikmoinsSkj_tij(T11,T12,T13,T22,T23,T33,VecX,VecXk,VecYk,VecZk,Wave,spec_rank,1.0_WP,1.0_WP)
    call RGMOikOjk_tij(S11,S12,S13,S22,S23,S33,VecX,VecXk,VecYk,VecZk,Wave,1.0_WP,1.0_WP)
    T11%values = delta**2*(coef1*T11%values + coef2*S11%values )
    T12%values = delta**2*(coef1*T12%values + coef2*S12%values )
    T13%values = delta**2*(coef1*T13%values + coef2*S23%values )
    T22%values = delta**2*(coef1*T22%values + coef2*S22%values )
    T23%values = delta**2*(coef1*T23%values + coef2*S23%values )
    T33%values = delta**2*(coef1*T33%values + coef2*S33%values )

    ! ===== Free memory and delete temporary fields =====
    call deleteDataLayout(S11)
    call deleteDataLayout(S12)
    call deleteDataLayout(S13)
    call deleteDataLayout(S23)
    call deleteDataLayout(S22)
    call deleteDataLayout(S33)

    res = .true.

end subroutine RGMSikmoinsSkjOikOjk_tij



!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> @param [in] VecX  velocity component along first direction in spectral space
!> @param [in] VecXk velocity component along first direction in spectral space
!> @param [in] VecYk velocity component along second direction in spectral space
!> @param [in] VecZk velocity component along third direction in spectral space
!> @param [in] 
!> @param [in] Wave wave numbers
!> @param [in] spec_rank rank of processor
!> @param [inout] 
!> @param [inout] 
!> @param [inout] 
!------------------------------------------------------------------------------
subroutine EVSGSMVreman_div(T1,T2,T3,VecX,VecXk,VecYk,VecZk,Wave,spec_rank,delta,coef,ite)

    !I/O data
    type(real_data_layout),intent(in)       :: VecX
    type(complex_data_layout),intent(in)    :: VecXk,VecYk,VecZk
    type(complex_data_layout),intent(inout) :: T1,T2,T3
    type(WaveNumbers),intent(in)            :: Wave
    integer,intent(in)                      :: spec_rank
    real(WP),intent(in),optional            :: coef,delta
    integer,intent(in),optional             :: ite

    !Local data
    type(real_data_layout)  :: T11,T12,T13,T22,T23,T33
    type(real_data_layout)  :: F11,F12,F13,F22,F23,F33
    logical                 :: pdffield,res
    real(WP)                :: coefLoc,deltaLoc

    pdffield= .false.
    if (present(ite)) then
      if ( ite .ne. -1 ) pdffield = .true.
    endif
    if (present(coef)) then
      coefLoc = coef
    else
      coefLoc = -5.0_WP*0.18_WP**2
    endif
    if (present(delta)) then
      deltaLoc = delta
    else
      deltaLoc = computeDelta(VecX)
    endif

    if (.not. copystructonly(VecX,T11) ) return
    if (.not. copystructonly(VecX,T12) ) return
    if (.not. copystructonly(VecX,T13) ) return
    if (.not. copystructonly(VecX,T22) ) return
    if (.not. copystructonly(VecX,T23) ) return
    if (.not. copystructonly(VecX,T33) ) return
    if (spec_rank .eq. 0) write (6,'(a)') '[INFO LES] sgs model for velocity: RGMSikmoinsSkj'

    call EVSGSMVreman_tij(T11,T12,T13,T22,T23,T33,VecX,VecXk,VecYk,VecZk,Wave,spec_rank,deltaLoc,coefLoc)
    call ComputeDivSymTijVec1Vec2(T1,T2,T3,T12,T13,T23,T11,T22,T33,Wave,res)
    if (.not. res) write(6,'(a)') '[ERROR] in RGMSikmoinsSkj_div'
    if ( pdffield ) then
      if (.not. copystructonly(VecX,F11) ) return
      if (.not. copystructonly(VecX,F12) ) return
      if (.not. copystructonly(VecX,F13) ) return
      if (.not. copystructonly(VecX,F22) ) return
      if (.not. copystructonly(VecX,F23) ) return
      if (.not. copystructonly(VecX,F33) ) return
      call computeStressTensor(VecXk,VecYk,VecZk,Wave,F11,F12,F13,F22,F23,F33,res)
      T11%values = T11%values*F11%values + T22%values*F22%values + T33%values*F33%values + &
         2.0_WP * (T12%values*F12%values + T13%values*F13%values + T23%values*F23%values )
      T11%name="diss_tijSij"
      coefLoc = computeFieldAvg(T11,spec_rank)
      if (spec_rank .eq. 0) write (6,'(a,1x,g15.8)') '[INFO LES] RGM OikOjk dissipation over the box ',coefLoc
      call btran(T1,T11,res)
      call btran(T2,T22,res)
      call btran(T3,T33,res)
      T11%name="sgs_vreman_T1"
      T22%name="sgs_vreman_T2"
      T33%name="sgs_vreman_T3"
      if (.not. computeFieldPDF(ite,T11,300,spec_rank,10000)) return
      if (.not. computeFieldPDF(ite,T22,300,spec_rank,10000)) return
      if (.not. computeFieldPDF(ite,T33,300,spec_rank,10000)) return
      call deleteDataLayout(F11)
      call deleteDataLayout(F12)
      call deleteDataLayout(F13)
      call deleteDataLayout(F23)
      call deleteDataLayout(F22)
      call deleteDataLayout(F33)
    endif

    call deleteDataLayout(T11)
    call deleteDataLayout(T12)
    call deleteDataLayout(T13)
    call deleteDataLayout(T23)
    call deleteDataLayout(T22)
    call deleteDataLayout(T33)

end subroutine EVSGSMVreman_div  


!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!> Vreman's model according to PoF 2004
!> @details
!> @param [in] VecX velocity component along first direction in spectral space
!> @param [in] VecXk velocity component along first direction in spectral space
!> @param [in] VecYk velocity component along second direction in spectral space
!> @param [in] VecZk velocity component along third direction in spectral space
!> @param [in] Wave wave numbers
!> @param [in] spec_rank rank of processor
!> @param [in] delta is the lenght of filter
!> @param [in] coef is the coefficient in front of the model (optional) 
!> @param [inout] T11 = coef * vu_sgs * S11 
!> @param [inout] T22 = coef * vu_sgs * S22  
!> @param [inout] T33 = coef * vu_sgs * S33  
!> @param [inout] T12 = coef * vu_sgs * S12  
!> @param [inout] T13 = coef * vu_sgs * S13  
!> @param [inout] T23 = coef * vu_sgs * S23  
!------------------------------------------------------------------------------
subroutine EVSGSMVreman_tij(T11,T12,T13,T22,T23,T33,VecX,VecXk,VecYk,VecZk,Wave,spec_rank,delta,coef)
  
    !I/O field
    TYPE(REAL_DATA_LAYOUT), INTENT(IN)    :: VecX
    TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: VecXk,VecYk,VecZk
    TYPE(REAL_DATA_LAYOUT),INTENT(INOUT)  :: T11,T12,T13,T22,T23,T33
    TYPE(WaveNumbers), INTENT(IN)         :: Wave
    INTEGER,INTENT(IN)                    :: spec_rank
    REAL(WP),INTENT(IN)                   :: delta,coef

    ! Local field
    type(REAL_DATA_LAYOUT) :: dudx,dudy,dudz
    type(REAL_DATA_LAYOUT) :: dvdx,dvdy,dvdz
    type(REAL_DATA_LAYOUT) :: dwdx,dwdy,dwdz
    type(REAL_DATA_LAYOUT) :: nu_sgs,beta
    TYPE(COMPLEX_DATA_LAYOUT) :: t1k,t2k,t3k
    logical                :: res

    ! ===== Initialisation =====
    res = .false.
    if (.not. copystructonly(VecXk,t1k) ) return
    if (.not. copystructonly(VecXk,t2k) ) return
    if (.not. copystructonly(VecXk,t3k) ) return
    if (.not. copystructonly(VecX,dudx) ) return
    if (.not. copystructonly(VecX,dudy) ) return
    if (.not. copystructonly(VecX,dudz) ) return
    if (.not. copystructonly(VecX,dvdx) ) return
    if (.not. copystructonly(VecX,dvdy) ) return
    if (.not. copystructonly(VecX,dvdz) ) return
    if (.not. copystructonly(VecX,dwdx) ) return
    if (.not. copystructonly(VecX,dwdy) ) return
    if (.not. copystructonly(VecX,dwdz) ) return
    if (.not. copystructonly(VecX,nu_sgs) ) return
    if (.not. copystructonly(VecX,beta) ) return
    call computeGradientKmain(Wave,VecXk,t1k,t2k,t3k)
    call btran(t1k,dudx,res)
    call btran(t2k,dudy,res)
    call btran(t3k,dudz,res)
    call computeGradientKmain(Wave,VecYk,t1k,t2k,t3k)
    call btran(t1k,dvdx,res)
    call btran(t2k,dvdy,res)
    call btran(t3k,dvdz,res)
    call computeGradientKmain(Wave,VecZk,t1k,t2k,t3k)
    call btran(t1k,dwdx,res)
    call btran(t2k,dwdy,res)
    call btran(t3k,dwdz,res)
    !Beta = B11B22 - B12^2 + B11B33 -B13^2 + B22B33 - B23^2
    !with Beta_ij=Delta^2*duidxk*dujdxk
    beta%values = (dudx%values**2+dudy%values**2+dudz%values**2)*(dvdx%values**2+dvdy%values**2+dvdz%values**2) +& !B11*B11
         & (dwdx%values**2+dwdy%values**2+dwdz%values**2)*(dvdx%values**2+dvdy%values**2+dvdz%values**2) +& !B22*B22
         & (dwdx%values**2+dwdy%values**2+dwdz%values**2)*(dudx%values**2+dudy%values**2+dudz%values**2) -& !B33*B33
         & (dudx%values * dvdx%values + dudy%values * dvdy%values + dudz%values * dvdz%values )**2 - &    !B12^2
         & (dvdx%values * dwdx%values + dvdy%values * dwdy%values + dvdz%values * dwdz%values )**2 - &    !B23^2
         & (dudx%values * dwdx%values + dudy%values * dwdy%values + dudz%values * dwdz%values )**2        !B13^2
    !alfa_ij = duidxj duidxj
    nu_sgs%values= dudx%values**2+dudy%values**2+dudz%values**2+&
                 &dvdx%values**2+dvdy%values**2+dvdz%values**2+&
                 &dwdx%values**2+dwdy%values**2+dwdz%values**2
    nu_sgs%values=-2.0_WP*coef*delta**2*sqrt(beta%values/nu_sgs%values)
    T11%values = coef*delta**2*nu_sgs%values*dudx%values
    T22%values = coef*delta**2*nu_sgs%values*dvdy%values
    T33%values = coef*delta**2*nu_sgs%values*dwdz%values
    T12%values = coef*delta**2*nu_sgs%values*0.5_WP*(dudy%values+dvdx%values)
    T13%values = coef*delta**2*nu_sgs%values*0.5_WP*(dudz%values+dwdx%values)
    T23%values = coef*delta**2*nu_sgs%values*0.5_WP*(dvdz%values+dwdy%values)

    ! ===== Free memory and delete temporary fields =====
    call deleteDataLayout(dudx)
    call deleteDataLayout(dudy)
    call deleteDataLayout(dudz)
    call deleteDataLayout(dvdx)
    call deleteDataLayout(dvdy)
    call deleteDataLayout(dvdz)
    call deleteDataLayout(dwdx)
    call deleteDataLayout(dwdy)
    call deleteDataLayout(dwdz)
    call deleteDataLayout(t1k)
    call deleteDataLayout(t2k)
    call deleteDataLayout(t3k)
    call deleteDataLayout(nu_sgs)
    call deleteDataLayout(beta)

    res = .true.

end subroutine EVSGSMVreman_tij 

!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!> Nicoud's model named Wall Adapting local Eddy viscosity (WALE) according to Flow Turbulence Combustion (1999) 
!> @details
!> @param [in] VecX velocity component along first direction in spectral space
!> @param [in] VecXk velocity component along first direction in spectral space
!> @param [in] VecYk velocity component along second direction in spectral space
!> @param [in] VecZk velocity component along third direction in spectral space
!> @param [in] Wave wave numbers
!> @param [in] spec_rank rank of processor
!> @param [in] delta is the lenght of filter
!> @param [in] coef is the coefficient in front of the model (optional) 
!> @param [inout] T11 = coef * vu_sgs * S11 
!> @param [inout] T22 = coef * vu_sgs * S22  
!> @param [inout] T33 = coef * vu_sgs * S33  
!> @param [inout] T12 = coef * vu_sgs * S12  
!> @param [inout] T13 = coef * vu_sgs * S13  
!> @param [inout] T23 = coef * vu_sgs * S23  
!------------------------------------------------------------------------------
subroutine WALE_tij(T11,T12,T13,T22,T23,T33,VecX,VecXk,VecYk,VecZk,Wave,spec_rank,delta,coef)
  
    !I/O field
    TYPE(REAL_DATA_LAYOUT), INTENT(IN)    :: VecX
    TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: VecXk,VecYk,VecZk
    TYPE(REAL_DATA_LAYOUT),INTENT(INOUT)  :: T11,T12,T13,T22,T23,T33
    TYPE(WaveNumbers), INTENT(IN)         :: Wave
    INTEGER,INTENT(IN)                    :: spec_rank
    REAL(WP),INTENT(IN)                   :: delta,coef

    ! Local field
    type(REAL_DATA_LAYOUT) :: dudx,dudy,dudz
    type(REAL_DATA_LAYOUT) :: dvdx,dvdy,dvdz
    type(REAL_DATA_LAYOUT) :: dwdx,dwdy,dwdz
    TYPE(COMPLEX_DATA_LAYOUT) :: t1k,t2k,t3k
    real(WP)               :: gsquarre11,gsquarre22,gsquarre33,gsquarre12,gsquarre13,gsquarre23
    real(WP)               :: gsquarre31,gsquarre32,gsquarre21,Sd11,Sd22,Sd33,Sd12,Sd13,Sd23,SdSd,SS
    real(WP)               :: nu_sgs,traceG2
    integer                :: i,j,k
    logical                :: res

    ! ===== Initialisation =====
    res = .false.
    if (.not. copystructonly(VecXk,t1k) ) return
    if (.not. copystructonly(VecXk,t2k) ) return
    if (.not. copystructonly(VecXk,t3k) ) return
    if (.not. copystructonly(VecX,dudx) ) return
    if (.not. copystructonly(VecX,dudy) ) return
    if (.not. copystructonly(VecX,dudz) ) return
    if (.not. copystructonly(VecX,dvdx) ) return
    if (.not. copystructonly(VecX,dvdy) ) return
    if (.not. copystructonly(VecX,dvdz) ) return
    if (.not. copystructonly(VecX,dwdx) ) return
    if (.not. copystructonly(VecX,dwdy) ) return
    if (.not. copystructonly(VecX,dwdz) ) return
    call computeGradientKmain(Wave,VecXk,t1k,t2k,t3k)
    call btran(t1k,dudx,res)
    call btran(t2k,dudy,res)
    call btran(t3k,dudz,res)
    call computeGradientKmain(Wave,VecYk,t1k,t2k,t3k)
    call btran(t1k,dvdx,res)
    call btran(t2k,dvdy,res)
    call btran(t3k,dvdz,res)
    call computeGradientKmain(Wave,VecZk,t1k,t2k,t3k)
    call btran(t1k,dwdx,res)
    call btran(t2k,dwdy,res)
    call btran(t3k,dwdz,res)
    do k= VecX%zmin , VecX%zmax
      do j= VecX%ymin , VecX%ymax
        do i= VecX%xmin , VecX%xmax
          traceG2 = 1.0_WP/3.0_WP*(dudx%values(i,j,k)**2         +dvdy%values(i,j,k)**2         +dwdz%values(i,j,k)**2 +&
                          &2.0_WP*(dudy%values(i,j,k)*dvdx%values(i,j,k)+dudz%values(i,j,k)*dwdx%values(i,j,k)+dvdz%values(i,j,k)*dwdy%values(i,j,k)))
          gsquarre11 = dudx%values(i,j,k)**2+dudy%values(i,j,k)*dvdx%values(i,j,k)+dudz%values(i,j,k)*dwdx%values(i,j,k)
          gsquarre22 = dvdx%values(i,j,k)*dudy%values(i,j,k)+dvdy%values(i,j,k)**2+dvdz%values(i,j,k)*dwdy%values(i,j,k)
          gsquarre33 = dwdx%values(i,j,k)*dudz%values(i,j,k)+dwdy%values(i,j,k)*dvdz%values(i,j,k)+dwdz%values(i,j,k)**2
          gsquarre12 = dudx%values(i,j,k)*dudy%values(i,j,k)+dudy%values(i,j,k)*dvdy%values(i,j,k)+dudz%values(i,j,k)*dwdy%values(i,j,k)
          gsquarre13 = dudx%values(i,j,k)*dudz%values(i,j,k)+dudy%values(i,j,k)*dvdz%values(i,j,k)+dudz%values(i,j,k)*dwdz%values(i,j,k) 
          gsquarre23 = dvdx%values(i,j,k)*dudz%values(i,j,k)+dvdy%values(i,j,k)*dvdz%values(i,j,k)+dvdz%values(i,j,k)*dwdz%values(i,j,k)
          gsquarre31 = dwdx%values(i,j,k)*dudx%values(i,j,k)+dwdy%values(i,j,k)*dvdx%values(i,j,k)+dwdz%values(i,j,k)*dwdx%values(i,j,k)
          gsquarre32 = dwdx%values(i,j,k)*dudy%values(i,j,k)+dwdy%values(i,j,k)*dvdy%values(i,j,k)+dwdz%values(i,j,k)*dwdy%values(i,j,k)
          gsquarre21 = dvdx%values(i,j,k)*dudx%values(i,j,k)+dvdy%values(i,j,k)*dvdx%values(i,j,k)+dvdz%values(i,j,k)*dwdx%values(i,j,k)
          Sd11 = gsquarre11 - traceG2
          Sd22 = gsquarre22 - traceG2
          Sd33 = gsquarre33 - traceG2
          Sd12 = 0.5_WP * ( gsquarre12 + gsquarre21 ) - traceG2
          Sd13 = 0.5_WP * ( gsquarre13 + gsquarre31 ) - traceG2
          Sd23 = 0.5_WP * ( gsquarre23 + gsquarre32 ) - traceG2
          SdSd = Sd11**2 + Sd22**2 + Sd33**2 + 2.0_WP*(Sd12**2 + Sd13**2 + Sd23**2)
          SS   = dudx%values(i,j,k)**2 + dvdy%values(i,j,k)**2 + dwdz%values(i,j,k)**2 +&
          2.0_WP*((0.5_WP*(dudy%values(i,j,k)+dvdx%values(i,j,k)))**2 +&
               & (0.5_WP*(dudz%values(i,j,k)+dwdx%values(i,j,k)))**2 +&
               & (0.5_WP*(dvdz%values(i,j,k)+dwdy%values(i,j,k)))**2)
          nu_sgs = SdSd**1.5_WP /( SS**2.5_WP + SdSd**0.8_WP ) 
          T11%values(i,j,k) = nu_sgs*dudx%values(i,j,k)
          T22%values(i,j,k) = nu_sgs*dvdy%values(i,j,k)
          T33%values(i,j,k) = nu_sgs*dwdz%values(i,j,k)
          T12%values(i,j,k) = nu_sgs*0.5_WP*(dudy%values(i,j,k)+dvdx%values(i,j,k))
          T13%values(i,j,k) = nu_sgs*0.5_WP*(dudz%values(i,j,k)+dwdx%values(i,j,k))
          T23%values(i,j,k) = nu_sgs*0.5_WP*(dvdz%values(i,j,k)+dwdy%values(i,j,k))
        enddo
      enddo
    enddo

    ! ===== Free memory and delete temporary fields =====
    call deleteDataLayout(dudx)
    call deleteDataLayout(dudy)
    call deleteDataLayout(dudz)
    call deleteDataLayout(dvdx)
    call deleteDataLayout(dvdy)
    call deleteDataLayout(dvdz)
    call deleteDataLayout(dwdx)
    call deleteDataLayout(dwdy)
    call deleteDataLayout(dwdz)
    call deleteDataLayout(t1k)
    call deleteDataLayout(t2k)
    call deleteDataLayout(t3k)

    res = .true.

end subroutine WALE_tij 

!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> Modèle du gradient en prenant en compte l'incompressibilité. Calcul de 3 termes par composantes du
!> tenseur sous maille au lieu de 6. Calcul de d/dxj(tau_ij)=Ti
!> d/dxj(tau_ij) = C delta^2 (duj/dxk d^2ui/dxjdxk) 
!> @param [in] VecX velocity component along first direction in spectral space
!> @param [in] VecY velocity component along second direction in spectral space
!> @param [in] VecZ velocity component along third direction in spectral space
!> @param [in] FieldForParam field to get informations about grid (physical space) 
!> @param [in] Wave wave numbers
!> @param [in] spec_rank rank of processor
!> @param [inout] T1 = eps_1jk eps_jmn C delta^2 duk/dxl d^2un/dxmdxl
!> @param [inout] T2 = eps_2jk eps_jmn C delta^2 duk/dxl d^2un/dxmdxl
!> @param [inout] T3 = eps_3jk eps_jmn C delta^2 duk/dxl d^2un/dxmdxl
!------------------------------------------------------------------------------
subroutine graddOkdxlOmegajl(T1,T2,T3,VecX,VecY,VecZ,Wave,spec_rank,delta,coef,ite)

  use mpilayout_tools , only : getnbcpus
  
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: VecX,VecY,VecZ
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT) :: T1,T2,T3
  TYPE(WaveNumbers), INTENT(IN) :: Wave
  INTEGER,INTENT(IN)            :: spec_rank
  REAL(WP),INTENT(IN),OPTIONAL  :: delta,coef
  INTEGER,INTENT(IN),OPTIONAL   :: ite

  TYPE(REAL_DATA_LAYOUT),DIMENSION(:),ALLOCATABLE :: doxd,doyd,dozd 
  TYPE(REAL_DATA_LAYOUT)                          :: O12,O23,O31
  TYPE(REAL_DATA_LAYOUT),DIMENSION(:),ALLOCATABLE :: TPhyBuf 
  REAL(WP) :: deltaLoc,coefLoc
  INTEGER :: nbcpus 
  logical :: res,pdffield=.false.

  if (present(ite)) then
    if ( ite .ne. -1 ) pdffield = .true.
  endif
  pdffield = .false.
  nbcpus=getnbcpus()
  if (.not. copystructonly(VecX,O23) )return
  if (.not. copystructonly(VecX,O31) ) return
  if (.not. copystructonly(VecX,O12) ) return
  IF (.NOT. samelayout(VecX,VecY,VecZ)) print *,'[ERROR]'
  IF (.NOT. initWorkArray(VecX,3,doxd) .OR. &
      .NOT. initWorkArray(VecX,3,doyd) .OR. &
      .NOT. initWorkArray(VecX,3,dozd) .OR. &
      .NOT. initWorkArray(VecX,3,TPhyBuf)) THEN 
    print *,'[ERROR]'
  ENDIF
  call ftran(VecX,T1,res)
  call ftran(VecY,T2,res)
  call ftran(VecZ,T3,res)
  call computeVorticity(T1,T2,T3,Wave,O23,O31,O12,res)
  if (.not. computeFieldGradient(Wave,O23,doxd(1),doxd(2),doxd(3),nbcpus,spec_rank)) return
  if (.not. computeFieldGradient(Wave,O31,doyd(1),doyd(2),doyd(3),nbcpus,spec_rank)) return
  if (.not. computeFieldGradient(Wave,O12,dozd(1),dozd(2),dozd(3),nbcpus,spec_rank)) return
  if (present(coef)) then
    coefLoc=coef
  else
    coefLoc=1.0_WP/12.0_WP
  endif
  if (present(delta)) then
    deltaLoc=delta
  else
    deltaLoc=computeDelta(VecX)
  endif 
  TPhyBuf(1)%values = coefLoc*deltaLoc**2*( dozd(1)%values*O12%values+dozd(3)%values*O23%values -&
                                   & doyd(1)%values*O31%values-doyd(2)%values*O23%values )
  TPhyBuf(2)%values = coefLoc*deltaLoc**2*( doxd(1)%values*O31%values+doxd(2)%values*O23%values -&
                                   & dozd(2)%values*O12%values-dozd(3)%values*O31%values )
  TPhyBuf(3)%values = coefLoc*deltaLoc**2*( doyd(2)%values*O12%values+doyd(3)%values*O31%values -&
                                   & doxd(1)%values*O12%values-doxd(3)%values*O23%values )
  if ( pdffield ) then
    TPhyBuf(1)%name="sgs_graddOkdxlSjlmoins_T1"
    TPhyBuf(2)%name="sgs_graddOkdxlSjlmoins_T2"
    TPhyBuf(3)%name="sgs_graddOkdxlSjlmoins_T3"
    if (.not. computeFieldPDF(ite,TPhyBuf(1),300,spec_rank,10000)) return
    if (.not. computeFieldPDF(ite,TPhyBuf(2),300,spec_rank,10000)) return
    if (.not. computeFieldPDF(ite,TPhyBuf(3),300,spec_rank,10000)) return
  endif
  CALL ftran(TPhyBuf(1),T1,res)
  if (.not. res) print *,'[ERROR]' 
  CALL ftran(TPhyBuf(2),T2,res)
  if (.not. res) print *,'[ERROR]' 
  CALL ftran(TPhyBuf(3),T3,res)
  if (.not. res) print *,'[ERROR]' 

  call deletedatalayout(O12)
  call deletedatalayout(O23)
  call deletedatalayout(O31)
  IF (.NOT. deleteWorkArray(doxd) .OR. &
     &.NOT. deleteWorkArray(doyd) .OR. &
     &.NOT. deleteWorkArray(dozd) .OR. &
     &.NOT. deleteWorkArray(TPhyBuf)) THEN
    print *,'[ERROR]'  
  ENDIF

end subroutine graddOkdxlOmegajl



!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> Modèle du gradient en prenant en compte l'incompressibilité. Calcul de 3 termes par composantes du
!> tenseur sous maille au lieu de 6. Calcul de d/dxj(tau_ij)=Ti
!> d/dxj(tau_ij) = C delta^2 (duj/dxk d^2ui/dxjdxk) 
!> @param [in] VecX velocity component along first direction in spectral space
!> @param [in] VecY velocity component along second direction in spectral space
!> @param [in] VecZ velocity component along third direction in spectral space
!> @param [in] FieldForParam field to get informations about grid (physical space) 
!> @param [in] Wave wave numbers
!> @param [in] spec_rank rank of processor
!> @param [inout] T1 = eps_1jk eps_jmn C delta^2 duk/dxl d^2un/dxmdxl
!> @param [inout] T2 = eps_2jk eps_jmn C delta^2 duk/dxl d^2un/dxmdxl
!> @param [inout] T3 = eps_3jk eps_jmn C delta^2 duk/dxl d^2un/dxmdxl
!------------------------------------------------------------------------------
subroutine graddOkdxlSjl(T1,T2,T3,VecX,VecY,VecZ,Wave,spec_rank,signe,delta,coef,ite)

  use mpilayout_tools , only : getnbcpus
  
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: VecX,VecY,VecZ
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT) :: T1,T2,T3
  TYPE(WaveNumbers), INTENT(IN) :: Wave
  INTEGER,INTENT(IN)            :: spec_rank
  CHARACTER(len=*),INTENT(IN)   :: signe
  REAL(WP),INTENT(IN),OPTIONAL  :: delta,coef
  INTEGER,INTENT(IN),OPTIONAL   :: ite

  TYPE(REAL_DATA_LAYOUT) :: Sm11,Sm12,Sm13,Sm22,Sm23,Sm33
  TYPE(REAL_DATA_LAYOUT),DIMENSION(:),ALLOCATABLE :: doxd,doyd,dozd 
  TYPE(REAL_DATA_LAYOUT),DIMENSION(:),ALLOCATABLE :: TPhyBuf 
  REAL(WP) :: deltaLoc,coefLoc
  INTEGER :: nbcpus 
  logical :: res,pdffield

  if (present(ite)) then
    if ( ite .ne. -1 ) pdffield = .true.
  endif
  pdffield = .false.
  nbcpus=getnbcpus()
  IF (.NOT. samelayout(VecX,VecY,VecZ)) print *,'[ERROR]'
  IF (.NOT. initWorkArray(VecX,3,doxd) .OR. &
      .NOT. initWorkArray(VecX,3,doyd) .OR. &
      .NOT. initWorkArray(VecX,3,dozd) .OR. &
      .NOT. initWorkArray(VecX,3,TPhyBuf) .OR. &
     &.NOT. copyStructOnly(VecX,Sm11) .OR. &
     &.NOT. copyStructOnly(VecX,Sm22) .OR. &
     &.NOT. copyStructOnly(VecX,Sm33) .OR. &
     &.NOT. copyStructOnly(VecX,Sm12) .OR. &
     &.NOT. copyStructOnly(VecX,Sm13) .OR. &
     &.NOT. copyStructOnly(VecX,Sm23)) THEN
    print *,'[ERROR]'
  ENDIF
  call ftran(VecX,T1,res)
  call ftran(VecY,T2,res)
  call ftran(VecZ,T3,res)
  call computeStressTensor(T1,T2,T3,Wave,doxd(1),doxd(2),doxd(3),doyd(1),doyd(2),doyd(3),res)
  if (.not. computeSplitingSij(doxd(1),doxd(2),doxd(3),doyd(1),doyd(2),doyd(3),Sm11,Sm12,Sm13,Sm22,Sm23,Sm33,signe)) return
  call computeVorticity(T1,T2,T3,Wave,TPhyBuf(1),TPhyBuf(2),TPhyBuf(3),res)
  if (.not. computeFieldGradient(Wave,TPhyBuf(1),doxd(1),doxd(2),doxd(3),nbcpus,spec_rank)) return
  if (.not. computeFieldGradient(Wave,TPhyBuf(2),doyd(1),doyd(2),doyd(3),nbcpus,spec_rank)) return
  if (.not. computeFieldGradient(Wave,TPhyBuf(3),dozd(1),dozd(2),dozd(3),nbcpus,spec_rank)) return
  if (present(coef)) then
    coefLoc=coef
  else
    coefLoc=1.0_WP/12.0_WP
  endif
  if (present(delta)) then
    deltaLoc=delta
  else
    deltaLoc=computeDelta(VecX)
  endif 
  TPhyBuf(1)%values = coefLoc*deltaLoc**2*( dozd(1)%values*Sm12%values+dozd(2)%values*Sm22%values+dozd(3)%values*Sm23%values -&
                                   & doyd(1)%values*Sm13%values-doyd(2)%values*Sm23%values-doyd(3)%values*Sm33%values )
  TPhyBuf(2)%values = coefLoc*deltaLoc**2*( doxd(1)%values*Sm13%values+doxd(2)%values*Sm23%values+doxd(3)%values*Sm33%values -&
                                   & dozd(1)%values*Sm11%values-dozd(2)%values*Sm12%values-dozd(3)%values*Sm13%values )
  TPhyBuf(3)%values = coefLoc*deltaLoc**2*( doyd(1)%values*Sm11%values+doyd(2)%values*Sm12%values+doyd(3)%values*Sm13%values -&
                                   & doxd(1)%values*Sm12%values-doxd(2)%values*Sm22%values-doxd(3)%values*Sm23%values )
  if ( pdffield ) then
    TPhyBuf(1)%name="sgs_graddOkdxlSjlmoins_T1"
    TPhyBuf(2)%name="sgs_graddOkdxlSjlmoins_T2"
    TPhyBuf(3)%name="sgs_graddOkdxlSjlmoins_T3"
    if (.not. computeFieldPDF(ite,TPhyBuf(1),300,spec_rank,10000)) return
    if (.not. computeFieldPDF(ite,TPhyBuf(2),300,spec_rank,10000)) return
    if (.not. computeFieldPDF(ite,TPhyBuf(3),300,spec_rank,10000)) return
  endif
  CALL ftran(TPhyBuf(1),T1,res)
  if (.not. res) print *,'[ERROR]' 
  CALL ftran(TPhyBuf(2),T2,res)
  if (.not. res) print *,'[ERROR]' 
  CALL ftran(TPhyBuf(3),T3,res)
  if (.not. res) print *,'[ERROR]' 

  if (.NOT. deleteWorkArray(doxd) .OR. &
     &.NOT. deleteWorkArray(doyd) .OR. &
     &.NOT. deleteWorkArray(dozd) .OR. &
     &.NOT. deleteWorkArray(TPhyBuf)) THEN
    print *,'[ERROR]'  
  ENDIF

end subroutine graddOkdxlSjl

!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> Modèle du gradient en prenant en compte l'incompressibilité. Calcul de 3 termes par composantes du
!> tenseur sous maille au lieu de 6. Calcul de d/dxj(tau_ij)=Ti
!> d/dxj(tau_ij) = C delta^2 (duj/dxk d^2ui/dxjdxk) 
!> @param [in] VecX velocity component along first direction in spectral space
!> @param [in] VecY velocity component along second direction in spectral space
!> @param [in] VecZ velocity component along third direction in spectral space
!> @param [in] FieldForParam field to get informations about grid (physical space) 
!> @param [in] Wave wave numbers
!> @param [in] spec_rank rank of processor
!> @param [inout] T1 = eps_1jk eps_jmn C delta^2 duk/dxl d^2un/dxmdxl
!> @param [inout] T2 = eps_2jk eps_jmn C delta^2 duk/dxl d^2un/dxmdxl
!> @param [inout] T3 = eps_3jk eps_jmn C delta^2 duk/dxl d^2un/dxmdxl
!------------------------------------------------------------------------------
subroutine graddOkdxldUjdxl(T1,T2,T3,VecX,VecY,VecZ,Wave,spec_rank,delta,coef,ite)

  use mpilayout_tools , only : getnbcpus
  
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: VecX,VecY,VecZ
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT) :: T1,T2,T3
  TYPE(WaveNumbers), INTENT(IN) :: Wave
  INTEGER,INTENT(IN)            :: spec_rank
  REAL(WP),INTENT(IN),OPTIONAL  :: delta,coef
  INTEGER,INTENT(IN),OPTIONAL   :: ite

  TYPE(REAL_DATA_LAYOUT),DIMENSION(:),ALLOCATABLE :: dud,dvd,dwd,doxd,doyd,dozd 
  TYPE(REAL_DATA_LAYOUT),DIMENSION(:),ALLOCATABLE :: TPhyBuf 
  REAL(WP) :: deltaLoc,coefLoc
  INTEGER :: nbcpus 
  logical :: res,pdffield

  if (present(ite)) then
    if ( ite .ne. -1 ) pdffield = .true.
  endif
  pdffield = .false.
  nbcpus=getnbcpus()
  IF (.NOT. samelayout(VecX,VecY,VecZ)) print *,'[ERROR]'
  IF (.NOT. initWorkArray(VecX,3,dud) .OR. &
      .NOT. initWorkArray(VecX,3,dvd) .OR. &
      .NOT. initWorkArray(VecX,3,dwd) .OR. &
      .NOT. initWorkArray(VecX,3,doxd) .OR. &
      .NOT. initWorkArray(VecX,3,doyd) .OR. &
      .NOT. initWorkArray(VecX,3,dozd) .OR. &
      .NOT. initWorkArray(VecX,3,TPhyBuf)) THEN
    print *,'[ERROR]'
  ENDIF
  call ftran(VecX,T1,res)
  call ftran(VecY,T2,res)
  call ftran(VecZ,T3,res)
  call computeVorticity(T1,T2,T3,Wave,dud(1),dud(2),dud(3),res)
  if (.not. computeFieldGradient(Wave,dud(1),doxd(1),doxd(2),doxd(3),nbcpus,spec_rank)) return
  if (.not. computeFieldGradient(Wave,dud(2),doyd(1),doyd(2),doyd(3),nbcpus,spec_rank)) return
  if (.not. computeFieldGradient(Wave,dud(3),dozd(1),dozd(2),dozd(3),nbcpus,spec_rank)) return
  if (.not. computeFieldGradient(Wave,VecX,dud(1),dud(2),dud(3),nbcpus,spec_rank)) return
  if (.not. computeFieldGradient(Wave,VecY,dvd(1),dvd(2),dvd(3),nbcpus,spec_rank)) return
  if (.not. computeFieldGradient(Wave,VecZ,dwd(1),dwd(2),dwd(3),nbcpus,spec_rank)) return
  if (present(coef)) then
    coefLoc=coef
  else
    coefLoc=1.0_WP/12.0_WP
  endif
  if (present(delta)) then
    deltaLoc=delta
  else
    deltaLoc=computeDelta(VecX)
  endif 
  TPhyBuf(1)%values = coefLoc*deltaLoc**2*( dozd(1)%values*dvd(1)%values+dozd(2)%values*dvd(2)%values+dozd(3)%values*dvd(3)%values -&
                                   & doyd(1)%values*dwd(1)%values-doyd(2)%values*dwd(2)%values-doyd(3)%values*dwd(3)%values )
  TPhyBuf(2)%values = coefLoc*deltaLoc**2*( doxd(1)%values*dwd(1)%values+doxd(2)%values*dwd(2)%values+doxd(3)%values*dwd(3)%values -&
                                   & dozd(1)%values*dud(1)%values-dozd(2)%values*dud(2)%values-dozd(3)%values*dud(3)%values )
  TPhyBuf(3)%values = coefLoc*deltaLoc**2*( doyd(1)%values*dvd(1)%values+doyd(2)%values*dvd(2)%values+doyd(3)%values*dvd(3)%values -&
                                   & doxd(1)%values*dvd(1)%values-doxd(2)%values*dvd(2)%values-doxd(3)%values*dvd(3)%values )
  if ( pdffield ) then
    TPhyBuf(1)%name="sgs_graddOkdxldUjdxl_T1"
    TPhyBuf(2)%name="sgs_graddOkdxldUjdxl_T2"
    TPhyBuf(3)%name="sgs_graddOkdxldUjdxl_T3"
    if (.not. computeFieldPDF(ite,TPhyBuf(1),300,spec_rank,10000)) return
    if (.not. computeFieldPDF(ite,TPhyBuf(2),300,spec_rank,10000)) return
    if (.not. computeFieldPDF(ite,TPhyBuf(3),300,spec_rank,10000)) return
  endif
  CALL ftran(TPhyBuf(1),T1,res)
  if (.not. res) print *,'[ERROR]' 
  CALL ftran(TPhyBuf(2),T2,res)
  if (.not. res) print *,'[ERROR]' 
  CALL ftran(TPhyBuf(3),T3,res)
  if (.not. res) print *,'[ERROR]' 

  IF (.NOT. deleteWorkArray(dud) .OR. &
     &.NOT. deleteWorkArray(dvd) .OR. &
     &.NOT. deleteWorkArray(dwd) .OR. &
     &.NOT. deleteWorkArray(doxd) .OR. &
     &.NOT. deleteWorkArray(doyd) .OR. &
     &.NOT. deleteWorkArray(dozd) .OR. &
     &.NOT. deleteWorkArray(TPhyBuf)) THEN
    print *,'[ERROR]'  
  ENDIF

end subroutine graddOkdxldUjdxl


!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> Modèle du gradient en prenant en compte l'incompressibilité. Calcul de 3 termes par composantes du
!> tenseur sous maille au lieu de 6. Calcul de d/dxj(tau_ij)=Ti
!> d/dxj(tau_ij) = C delta^2 (duj/dxk d^2ui/dxjdxk) 
!> @param [in] VecX velocity component along first direction in spectral space
!> @param [in] VecY velocity component along second direction in spectral space
!> @param [in] VecZ velocity component along third direction in spectral space
!> @param [in] FieldForParam field to get informations about grid (physical space) 
!> @param [in] Wave wave numbers
!> @param [in] spec_rank rank of processor
!> @param [inout] T1 = eps_1jk eps_jmn C delta^2 duk/dxl d^2un/dxmdxl
!> @param [inout] T2 = eps_2jk eps_jmn C delta^2 duk/dxl d^2un/dxmdxl
!> @param [inout] T3 = eps_3jk eps_jmn C delta^2 duk/dxl d^2un/dxmdxl
!------------------------------------------------------------------------------
subroutine gradientQDMvort(T1,T2,T3,VecX,VecY,VecZ,Wave,spec_rank,delta,coef,ite)

  
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: VecX,VecY,VecZ
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT) :: T1,T2,T3
  TYPE(WaveNumbers), INTENT(IN) :: Wave
  INTEGER,INTENT(IN)            :: spec_rank
  REAL(WP),INTENT(IN),OPTIONAL  :: delta,coef
  INTEGER,INTENT(IN),OPTIONAL   :: ite

  TYPE(REAL_DATA_LAYOUT),DIMENSION(:),ALLOCATABLE :: dukdxl
  TYPE(REAL_DATA_LAYOUT),DIMENSION(:),ALLOCATABLE :: dundxmdxl
  TYPE(REAL_DATA_LAYOUT),DIMENSION(:),ALLOCATABLE :: UBuf
  TYPE(REAL_DATA_LAYOUT),DIMENSION(:),ALLOCATABLE :: TPhyBuf 
  REAL(WP) :: deltaLoc,coefLoc
  INTEGER :: i,j,k,l,m,n
  logical :: res,pdffield

  pdffield=.false.
  if (present(ite)) then
    if ( ite .ne. -1 ) pdffield = .true.
  endif  

  IF (.NOT. samelayout(VecX,VecY,VecZ)) print *,'[ERROR]'
  IF (.NOT. initWorkArray(VecX,3,dukdxl) .OR. &
      .NOT. initWorkArray(VecX,3,dundxmdxl) .OR. &
      .NOT. initWorkArray(VecX,3,UBuf) .OR. &
      .NOT. initWorkArray(VecX,3,TPhyBuf)) THEN
    print *,'[ERROR]'
  ENDIF

  Ubuf(1)%values=VecX%values
  Ubuf(2)%values=VecY%values
  Ubuf(3)%values=VecZ%values
  if (present(coef)) then
    coefLoc=coef
  else
    coefLoc=1.0_WP/12.0_WP
  endif
  if (present(delta)) then
    deltaLoc=delta
  else
    deltaLoc=computeDelta(VecX)
  endif 
  do i=1,3
    TPhyBuf(i)%values=0.0_WP
    do j=1,3
      do k=1,3
        do l=1,3
          do m=1,3
            do n=1,3
              if (.not.( j == m .or. j == n .or. m == n .or. i == j .or. i == k .or. j == k) ) then
                if (.not. computeDerivationField(Ubuf(k),(/l/),dukdxl(k),spec_rank,Wave)) return
                if (.not. computeDerivationField(Ubuf(n),(/m,l/),dundxmdxl(n),spec_rank,Wave)) return 
                TPhyBuf(i)%values = levicivita(i,j,k)*levicivita(j,m,n)*dukdxl(k)%values*dundxmdxl(n)%values + TPhyBuf(i)%values
              endif 
            enddo
          enddo
        enddo
      enddo
    enddo
    TPhyBuf(i)%values=deltaLoc**2*coefLoc*TPhyBuf(i)%values
  enddo
  if ( pdffield ) then
    TPhyBuf(1)%name="sgs_gradientQDMVort_T1"
    TPhyBuf(2)%name="sgs_gradientQDMVort_T2"
    TPhyBuf(3)%name="sgs_gradientQDMVort_T3"
    if (.not. computeFieldPDF(ite,TPhyBuf(1),300,spec_rank,10000)) return
    if (.not. computeFieldPDF(ite,TPhyBuf(2),300,spec_rank,10000)) return
    if (.not. computeFieldPDF(ite,TPhyBuf(3),300,spec_rank,10000)) return
  endif
  CALL ftran(TPhyBuf(1),T1,res)
  if (.not. res) print *,'[ERROR]' 
  CALL ftran(TPhyBuf(2),T2,res)
  if (.not. res) print *,'[ERROR]' 
  CALL ftran(TPhyBuf(3),T3,res)
  if (.not. res) print *,'[ERROR]' 

  IF (.NOT. deleteWorkArray(dukdxl) .OR. &
     &.NOT. deleteWorkArray(dundxmdxl) .OR. &
     &.NOT. deleteWorkArray(UBuf) .OR. &
     &.NOT. deleteWorkArray(TPhyBuf)) THEN
    print *,'[ERROR]'  
  ENDIF

end subroutine gradientQDMVort


!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> Modèle du gradient en prenant en compte l'incompressibilité. Calcul de 3 termes par composantes du
!> tenseur sous maille au lieu de 6. Calcul de d/dxj(tau_ij)=Ti
!> d/dxj(tau_ij) = C delta^2 (duj/dxk d^2ui/dxjdxk) 
!> Therefore, the computational time is increased because of the function used....
!> @param [in] VecX velocity component along first direction in spectral space
!> @param [in] VecY velocity component along second direction in spectral space
!> @param [in] VecZ velocity component along third direction in spectral space
!> @param [in] FieldForParam field to get informations about grid (physical space) 
!> @param [in] Wave wave numbers
!> @param [in] spec_rank rank of processor
!> @param [inout] T1 = C delta^2 (duj/dxk d^2u1/dxjdxk) in spectral space
!> @param [inout] T2 = C delta^2 (duj/dxk d^2u2/dxjdxk) in spectral space
!> @param [inout] T3 = C delta^2 (duj/dxk d^2u3/dxjdxk) in spectral space
!------------------------------------------------------------------------------
subroutine gradientVelIncomp(T1,T2,T3,VecX,VecY,VecZ,Wave,spec_rank,delta,coef,ite)

  
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: VecX,VecY,VecZ
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT) :: T1,T2,T3
  TYPE(WaveNumbers), INTENT(IN) :: Wave
  INTEGER,INTENT(IN)            :: spec_rank
  REAL(WP),INTENT(IN),OPTIONAL  :: delta,coef
  INTEGER,INTENT(IN),OPTIONAL   :: ite

  TYPE(REAL_DATA_LAYOUT),DIMENSION(:),ALLOCATABLE :: dujdxk
  TYPE(REAL_DATA_LAYOUT),DIMENSION(:),ALLOCATABLE :: duidxkdxj
  TYPE(REAL_DATA_LAYOUT),DIMENSION(:),ALLOCATABLE :: UBuf
  TYPE(REAL_DATA_LAYOUT),DIMENSION(:),ALLOCATABLE :: TPhyBuf 
  REAL(WP) :: deltaLoc,coefLoc
  INTEGER :: i,j,k
  logical :: res,pdffield

  pdffield=.false.
  if (present(ite)) then
    if ( ite .ne. -1 ) pdffield = .true.
  endif

  IF (.NOT. samelayout(VecX,VecY,VecZ)) print *,'[ERROR]'
  IF (.NOT. initWorkArray(VecX,3,dujdxk) .OR. &
      .NOT. initWorkArray(VecX,3,duidxkdxj) .OR. &
      .NOT. initWorkArray(VecX,3,UBuf) .OR. &
      .NOT. initWorkArray(VecX,3,TPhyBuf)) THEN
    print *,'[ERROR]'
  ENDIF
  !TODO remplacer Ubuf par un pointeur qui pointe successivement sur VecX, VecY et VecZ
  !TODO remplacer TPhyBuf par un pointeur qui pointe successivement sur T1,T2 et T3 

  Ubuf(1)%values=VecX%values
  Ubuf(2)%values=VecY%values
  Ubuf(3)%values=VecZ%values
  if (present(coef)) then
    coefLoc=coef
  else
    coefLoc=1.0_WP/12.0_WP
  endif
  if (present(delta)) then
    deltaLoc=delta
  else
    deltaLoc=computeDelta(VecX)
  endif 
  do i=1,3
    TPhyBuf(i)%values=0.0_WP
    do j=1,3
      do k=1,3
        if (.not. computeDerivationField(Ubuf(j),(/k/),dujdxk(k),spec_rank,Wave)) then
          print *,'[ERROR]' 
        endif
        if (.not. computeDerivationField(Ubuf(i),(/j,k/),duidxkdxj(k),spec_rank,Wave)) then
          print *,'[ERROR]'
        endif
        TPhyBuf(i)%values = TPhyBuf(i)%values + dujdxk(k)%values * duidxkdxj(k)%values
      enddo
    enddo
    TPhyBuf(i)%values=deltaLoc**2*coefLoc*TPhyBuf(i)%values
  enddo
  if ( pdffield ) then
    TPhyBuf(1)%name="sgs_gradientVelIncomp_T1"
    TPhyBuf(2)%name="sgs_gradientVelIncomp_T2"
    TPhyBuf(3)%name="sgs_gradientVelIncomp_T3"
    if (.not. computeFieldPDF(ite,TPhyBuf(1),300,spec_rank,10000)) return
    if (.not. computeFieldPDF(ite,TPhyBuf(2),300,spec_rank,10000)) return
    if (.not. computeFieldPDF(ite,TPhyBuf(3),300,spec_rank,10000)) return
  endif
  CALL ftran(TPhyBuf(1),T1,res)
  if (.not. res) print *,'[ERROR]' 
  CALL ftran(TPhyBuf(2),T2,res)
  if (.not. res) print *,'[ERROR]' 
  CALL ftran(TPhyBuf(3),T3,res)
  if (.not. res) print *,'[ERROR]' 

  IF (.NOT. deleteWorkArray(dujdxk) .OR. &
     &.NOT. deleteWorkArray(duidxkdxj) .OR. &
     &.NOT. deleteWorkArray(UBuf) .OR. &
     &.NOT. deleteWorkArray(TPhyBuf)) THEN
    print *,'[ERROR]'  
  ENDIF

end subroutine gradientVelIncomp


!------------------------------------------------------------------------------
!> @author 
!> Mouloud Kessar, LEGI
!>
!> @details
!> @param [in] VecX velocity component along first direction
!> @param [in] VecY velocity component along second direction
!> @param [in] VecZ velocity component along third direction
!> @param [in] Wave wave numbers
!> @param [in] spec_rank rank of processor
!> @param [out]
!------------------------------------------------------------------------------
subroutine gradientVel(t11,t12,t13,t22,t23,t33,VecX,VecY,VecZ,Wave,delta,spec_rank,nbcpus)

implicit none

  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: VecX,VecY,VecZ
  TYPE(REAL_DATA_LAYOUT) :: dfiudx, dfiudy, dfiudz
  TYPE(REAL_DATA_LAYOUT) :: dfjudx, dfjudy, dfjudz
  INTEGER,INTENT(IN)     :: spec_rank,nbcpus
  TYPE(REAL_DATA_LAYOUT), INTENT(INOUT) :: t11,t12,t13,t22,t23,t33
  REAL(WP) :: delta
  logical :: res
  TYPE(WaveNumbers), INTENT(IN) :: Wave

  IF ((.NOT. copyStructOnly(VecX,dfiudx)) .OR. &
     &(.NOT. copyStructOnly(VecX,dfiudy)) .OR. &
     &(.NOT. copyStructOnly(VecX,dfiudz)) .OR. &
     &(.NOT. copyStructOnly(VecX,dfjudx)) .OR. &
     &(.NOT. copyStructOnly(VecX,dfjudy)) .OR. &
     &(.NOT. copyStructOnly(VecX,dfjudz))    ) RETURN 

    res = computeFieldGradient(wave, VecX, dfiudx, dfiudy, dfiudz,nbcpus,spec_rank)
    res = computeFieldGradient(wave, VecY, dfjudx, dfjudy, dfjudz,nbcpus,spec_rank)

    t11%values = (delta*delta/12)*( dfiudx%values*dfiudx%values + dfiudy%values*dfiudy%values + dfiudz%values*dfiudz%values )

    t12%values = (delta*delta/12)*( dfiudx%values*dfjudx%values + dfiudy%values*dfjudy%values + dfiudz%values*dfjudz%values )

    res = computeFieldGradient(wave, VecZ, dfiudx, dfiudy, dfiudz,nbcpus,spec_rank)

    t22%values = (delta*delta/12)*( dfjudx%values*dfjudx%values + dfjudy%values*dfjudy%values + dfjudz%values*dfjudz%values )

    t23%values = (delta*delta/12)*( dfiudx%values*dfjudx%values + dfiudy%values*dfjudy%values + dfiudz%values*dfjudz%values )

    t33%values = (delta*delta/12)*( dfiudx%values*dfiudx%values + dfiudy%values*dfiudy%values + dfiudz%values*dfiudz%values )

    res = computeFieldGradient(wave, VecX, dfjudx, dfjudy, dfjudz,nbcpus,spec_rank)

    t13%values = (delta*delta/12)*( dfiudx%values*dfjudx%values + dfiudy%values*dfjudy%values + dfiudz%values*dfjudz%values )


    CALL deleteDataLayout(dfiudx)
    CALL deleteDataLayout(dfiudy)
    CALL deleteDataLayout(dfiudz)
    CALL deleteDataLayout(dfjudx)
    CALL deleteDataLayout(dfjudy)
    CALL deleteDataLayout(dfjudz)

end subroutine gradientVel

!------------------------------------------------------------------------------
!> @author 
!> Mouloud Kessar, LEGI
!>
!> @details
!> @param [in] VecX velocity component along first direction
!> @param [in] VecY velocity component along second direction
!> @param [in] VecZ velocity component along third direction
!> @param [in] Wave wave numbers
!> @param [in] spec_rank rank of processor
!> @param [out]
!------------------------------------------------------------------------------

subroutine SimilVelDIV(T1,T2,T3,VecX,VecY,VecZ,VecXk,VecYk,VecZk,Wave,delta,filter,spec_rank,nbcpus)

implicit none

  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: VecX,VecY,VecZ
  INTEGER :: spec_rank,nbcpus,filter
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT) :: T1,T2,T3
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: VecXk,VecYk,VecZk

  REAL(WP) :: delta
  logical :: res
  TYPE(WaveNumbers), INTENT(IN) :: Wave
  TYPE(REAL_DATA_LAYOUT) :: t11,t12,t13,t22,t23,t33


  IF ((.NOT. copyStructOnly(VecX,t11)) .OR. &
     &(.NOT. copyStructOnly(VecX,t12)) .OR. &
     &(.NOT. copyStructOnly(VecX,t13)) .OR. &
     &(.NOT. copyStructOnly(VecX,t22)) .OR. &
     &(.NOT. copyStructOnly(VecX,t23)) .OR. &
     &(.NOT. copyStructOnly(VecX,t33))    ) RETURN 

  call SimilVel(t11,t12,t13,t22,t23,t33,VecX,VecY,VecZ,VecXk,VecYk,VecZk,Wave,delta,filter,spec_rank,nbcpus)

  call ComputeDivSymTijVec1Vec2(T1,T2,T3,T12,T13,T23,T11,T22,T33,wave,res)


    CALL deleteDataLayout(t11)
    CALL deleteDataLayout(t12)
    CALL deleteDataLayout(t13)
    CALL deleteDataLayout(t22)
    CALL deleteDataLayout(t23)
    CALL deleteDataLayout(t33)

end subroutine SimilVelDIV


!------------------------------------------------------------------------------
!> @author 
!> Mouloud Kessar, LEGI
!>
!> @details
!> @param [in] VecX velocity component along first direction
!> @param [in] VecY velocity component along second direction
!> @param [in] VecZ velocity component along third direction
!> @param [in] Wave wave numbers
!> @param [in] spec_rank rank of processor
!> @param [out]
!------------------------------------------------------------------------------
subroutine SimilVel(t11,t12,t13,t22,t23,t33,VecX,VecY,VecZ,VecXk,VecYk,VecZk,Wave,delta,filter,spec_rank,nbcpus)

implicit none

  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: VecX,VecY,VecZ
  TYPE(REAL_DATA_LAYOUT) :: VecXf,VecYf,VecZf
  TYPE(REAL_DATA_LAYOUT) :: Wtab1
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: VecXk,VecYk,VecZk
  TYPE(COMPLEX_DATA_LAYOUT) :: tabk1,tabk2
  INTEGER :: spec_rank,nbcpus,filter
  TYPE(REAL_DATA_LAYOUT), INTENT(INOUT) :: t11,t12,t13,t22,t23,t33
  REAL(WP) :: delta
  REAL(WP) :: delta2,Ct
  logical :: res
  TYPE(WaveNumbers), INTENT(IN) :: Wave
  delta2=2*delta
  Ct=0.91_WP
  IF ((.NOT. copyStructOnly(VecX,Wtab1)) .OR. &
     &(.NOT. copyStructOnly(VecXk,tabk1)) .OR. &
     &(.NOT. copyStructOnly(VecX,VecXf)) .OR. &
     &(.NOT. copyStructOnly(VecX,VecYf)) .OR. &
     &(.NOT. copyStructOnly(VecX,VecZf)) .OR. &
     &(.NOT. copyStructOnly(VecXk,tabk2))    ) RETURN 

    CALL computeFilter(wave,delta2,VecXk,tabk2,filter)
    call btran(tabk2,VecXf,res)

    CALL computeFilter(wave,delta2,VecYk,tabk2,filter)
    call btran(tabk2,VecYf,res)

    CALL computeFilter(wave,delta2,VecZk,tabk2,filter)
    call btran(tabk2,VecZf,res)

    !calcul de T11

    wtab1%values=VecX%values*VecX%values
    call ftran(wtab1,tabk1,res)
    CALL computeFilter(wave,delta2,tabk1,tabk2,filter)
    call btran(tabk2,Wtab1,res)
    
    t11%values = Ct*(Wtab1%values -VecXf%values*VecXf%values) 

    !calcul de T12

    wtab1%values=VecX%values*VecY%values
    call ftran(wtab1,tabk1,res)
    CALL computeFilter(wave,delta2,tabk1,tabk2,filter)
    call btran(tabk2,Wtab1,res)
    
    t12%values = Ct*(Wtab1%values -VecXf%values*VecYf%values) 

   !calcul de T13

    wtab1%values=VecX%values*VecZ%values
    call ftran(wtab1,tabk1,res)
    CALL computeFilter(wave,delta2,tabk1,tabk2,filter)
    call btran(tabk2,Wtab1,res)
    
    t13%values = Ct*(Wtab1%values -VecXf%values*VecZf%values) 

    !calcul de T22

    wtab1%values=VecY%values*VecY%values
    call ftran(wtab1,tabk1,res)
    CALL computeFilter(wave,delta2,tabk1,tabk2,filter)
    call btran(tabk2,Wtab1,res)
    
    t22%values = Ct*(Wtab1%values -VecYf%values*VecYf%values) 

    !calcul de T23

    wtab1%values=VecY%values*VecZ%values
    call ftran(wtab1,tabk1,res)
    CALL computeFilter(wave,delta2,tabk1,tabk2,filter)
    call btran(tabk2,Wtab1,res)
    
    t23%values = Ct*(Wtab1%values -VecYf%values*VecZf%values) 
    
    !calcul de T33

    wtab1%values=VecZ%values*VecZ%values
    call ftran(wtab1,tabk1,res)
    CALL computeFilter(wave,delta2,tabk1,tabk2,filter)
    call btran(tabk2,Wtab1,res)
    
    t33%values = Ct*(Wtab1%values -VecZf%values*VecZf%values) 
    
    CALL deleteDataLayout(tabk1)
    CALL deleteDataLayout(tabk2)
    CALL deleteDataLayout(Wtab1)
    CALL deleteDataLayout(VecXf)
    CALL deleteDataLayout(VecYf)
    CALL deleteDataLayout(VecZf)


end subroutine

end module velmodels
!> @}
