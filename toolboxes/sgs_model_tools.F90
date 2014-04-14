!> @addtogroup toolbox 
!! {@

module subgrid_tools

 use toolbox
 use stat_tools 
 use differential_tools 
 use transforms_tools
 use wavenumber_tools
 use datalayout   
 use filtering_tools 
 use physical_values_tools

 implicit none
 
 public

 interface computeLilly 
   module procedure computeLilly_2d 
   module procedure computeLilly_3d
 end interface computeLilly

 interface computeSplitingSij
   module procedure computeSplitingSij_analyt_sca
! these two following function cannot be set at the same time
! as they have the same interface. But you can choose to use
! Lapack or not in order to compute eigen values and eigen vectors.
! Of course be sur to use -DLAPACK flag during compilation
   module procedure computeSplitingSij_analyt_field
!   module procedure computeSplitingSij_lapack_field
 end interface computeSplitingSij


 contains
 
!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> Compute Lilly's term <LiMi>/<MiMi> 
!> @param [in] M1 is the first component of Mi in physical space
!> @param [in] M2 is the second component of Mi in physical space
!> @param [in] M3 is the third component of Mi in physical space
!> @param [in] L1 is the first component of Li in physical space
!> @param [in] L2 is the second component of Li in physical space
!> @param [in] L3 is the third component of Li in physical space
!> @param [out] Lilly = <LiMi>/<MiMi> with <.> the spatial average 
!> @param [in] spec_rank is the number of process in cpu pool 
!------------------------------------------------------------------------------
 subroutine computeLilly_2d(M1,M2,M3,L1,L2,L3,Lilly,spec_rank,dir2DAvg)

  !I/O
  type(real_data_layout),intent(in)   :: M1,M2,M3,L1,L2,L3
  integer,intent(in)                  :: spec_rank
  real(WP),intent(inout),dimension(:) :: Lilly
  integer,intent(in)                  :: dir2DAvg

  !Local Data
  type(real_data_layout)              :: temp
  real(WP),dimension(:),allocatable   :: tempArray
  integer                             :: ires,i
  integer                             :: nbpoints=0
 

  select case (dir2DAvg)
    case(1)
      nbpoints=M1%nx
    case(2)
      nbpoints=M1%ny
    case(3)
      nbpoints=M1%nz
    case default 
      write(6,'(a,i0)')'[ERROR] computeLilly : dir2DAvg unknown'
  end select
 
  if ( .not. sameLayout( M1,M2,M3,L1) .or. &
     & .not. sameLayout( L2,L3,L1 )) then 
     write(6,'(a,i0)')'[ERROR] computeLilly : not same data layout !',spec_rank
  endif
  if ( .not. copyStructOnly( M1,temp )) then
     write(6,'(a,i0)')'[ERROR] computeLilly : not enought memory !',spec_rank
  endif
  if ( size(Lilly,1) .ne. nbpoints ) then
     write(6,'(a,i0)')'[ERROR] size of array equal to ',M1%nz 
  endif 
  allocate(tempArray(nbpoints),stat=ires)
  if ( ires .ne. 0) then
     write(6,'(a,i0)')'[ERROR] computeLilly : not enought memory !',spec_rank
  endif
  temp%values = L1%values * M1%values + L2%values * M2%values + L3%values * M3%values
  if (.not. computeFieldAvg2D( temp , Lilly, spec_rank , dir2DAvg )) then
    write(6,'(a,i0)')'[ERROR] computeLilly_jet: computeFieldAvg2d failed !',spec_rank
  endif
  temp%values = M1%values * M1%values + M2%values * M2%values + M3%values * M3%values
  if (.not. computeFieldAvg2D( temp , tempArray, spec_rank , dir2DAvg )) then
    write(6,'(a,i0)')'[ERROR] computeLilly_jet: computeFieldAvg2d failed !',spec_rank
  endif
  do i = 1 , nbpoints
    Lilly(i) = Lilly(i) / tempArray(i)
  enddo

  call deleteDataLayout(temp)
  deallocate(tempArray)

 end subroutine computeLilly_2d

!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> Compute Lilly's term <LiMi>/<MiMi> 
!> @param [in] M1 is the first component of Mi in physical space
!> @param [in] M2 is the second component of Mi in physical space
!> @param [in] M3 is the third component of Mi in physical space
!> @param [in] L1 is the first component of Li in physical space
!> @param [in] L2 is the second component of Li in physical space
!> @param [in] L3 is the third component of Li in physical space
!> @param [out] Lilly = <LiMi>/<MiMi> with <.> the spatial average 
!> @param [in] spec_rank is the number of process in cpu pool 
!------------------------------------------------------------------------------
 subroutine computeLilly_3d(M1,M2,M3,L1,L2,L3,Lilly,spec_rank)

  !I/O
  type(real_data_layout),intent(in)   :: M1,M2,M3,L1,L2,L3
  integer,intent(in)                  :: spec_rank
  real(WP),intent(inout)              :: Lilly

  !Local Data
  type(real_data_layout)              :: temp
  
  if ( .not. sameLayout( M1,M2,M3,L1) .or. &
     & .not. sameLayout( L2,L3,L1 )) then 
     write(6,'(a,i0)')'[ERROR] computeLilly : not same data layout !',spec_rank
  endif

  if ( .not. copyStructOnly( M1,temp ) ) then
     write(6,'(a,i0)')'[ERROR] computeLilly : not enought memory !',spec_rank
  endif

  temp%values = L1%values * M1%values + L2%values * M2%values + L3%values * M3%values
  Lilly = computeFieldAvg( temp , spec_rank )
  temp%values = M1%values * M1%values + M2%values * M2%values + M3%values * M3%values
  Lilly = Lilly / computeFieldAvg( temp , spec_rank )

  call deleteDataLayout(temp)


 end subroutine computeLilly_3d


!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> Compute Germano identity for Scalar
!> @param [in] wave is the waves number of quantities 
!> @param [in] U is the first component of velocity in physical space and filtered 
!> @param [in] V is the second component of velocity in physical space and filtered 
!> @param [in] W is the third component of velocity in physical space and filtered 
!> @param [in] Scal is the scalar in physical space and filtered 
!> @param [in] Uk is the first component of velocity in Fourier space and filtered 
!> @param [in] Vk is the second component of velocity in Fourier space and filtered 
!> @param [in] Wk is the third component of velocity in Fourier space and filtered 
!> @param [in] Scalk is the scalar in Fourier space and filtered 
!> @param [in] Filter is the type test filter applied see "computeFilter" in sub_grid_tools.F90 
!> @param [in] delta is the size of a posteriori or a priori test 
!> @param [inout] L1 is (Ubar*Zbar)Tilde - UbarTilde * ZbarTilde (tilde is the test filter) 
!> @param [inout] L2 is (Vbar*Zbar)Tilde - VbarTilde * ZbarTilde (tilde is the test filter) 
!> @param [inout] L3 is (Wbar*Zbar)Tilde - WbarTilde * ZbarTilde (tilde is the test filter) 
!------------------------------------------------------------------------------
 subroutine computeGermanoSca(wave,U,V,W,Scal,Uk,Vk,Wk,Scalk,Filter,delta,L1,L2,L3)

  !I/O
  type(complex_data_layout),intent(in)    :: Uk,Vk,Wk,Scalk
  type(real_data_layout),intent(in)       :: U,V,W,Scal
  type(waveNumbers),intent(in)            :: wave
  integer,intent(in)                      :: Filter
  real(WP),intent(in)                     :: delta
  type(real_data_layout),intent(inout)    :: L1,L2,L3

  !Local Data
  logical                   :: res
  type(complex_data_layout) :: TildeK,tempK
  type(real_data_layout)    :: temp2,temp1,ScalTilde
  real(WP)                  :: sizeTestFilter 
  integer                   :: spec_rank

  sizeTestFilter = 2.0_WP

  if ( .not. sameLayout(U,V,W,Scal) .or. &
     & .not. sameLayout(Uk,Vk,Wk,Scalk) .or. &
     & .not. sameLayout(L1,L2,U) .or. &
     & .not. sameLayout(L3,U) ) then 
      write(6,'(a)')'[ERROR] in computeLSca : data layout are not the same !'
  endif
  if ( .not. copyStructOnly(Uk,TildeK) .or. &
     & .not. copyStructOnly(Uk,tempK) .or. &
     & .not. copyStructOnly(U,temp1) .or. &
     & .not. copyStructOnly(U,temp2) .or. &
     & .not. copyStructOnly(Scal,ScalTilde)) then
      write(6,'(a)')'[ERROR] in computeLSca : not enought memory !'
  endif
  spec_rank=getnbmycpu() 

  ! Sca|~
  call computeFilter(wave,sizeTestFilter*delta,Scalk,TildeK,filter)
  call btran(TildeK,ScalTilde,res)
  if ( .not. res ) write(6,'(a)')'[ERROR] in computeGermanoSca : btran failed !' 

  ! U| x Sca|
  temp2%values = U%values * Scal%values
  call ftran(temp2,tempK,res)
  if ( .not. res ) write(6,'(a)')'[ERROR] in computeGermanoSca : ftran failed !' 
  ! (U| x Sca|)~
  call computeFilter(wave,sizeTestFilter*delta,tempK,TildeK,filter)
  call btran(TildeK,temp2,res)
  if ( .not. res ) write(6,'(a)')'[ERROR] in computeGermanoSca : btran failed !' 
  ! U|~ 
  call computeFilter(wave,sizeTestFilter*delta,Uk,TildeK,filter)
  call btran(TildeK,temp1,res)
  if ( .not. res ) write(6,'(a)')'[ERROR] in computeGermanoSca : btran failed !' 
  ! (U| x Sca|)~ - U|~ x Sca|~
  L1%values = temp2%values - temp1%values * ScalTilde%values

  ! V| x Sca|
  temp2%values = V%values * Scal%values
  call ftran(temp2,tempK,res)
  if ( .not. res ) write(6,'(a)')'[ERROR] in computeGermanoSca : ftran failed !' 
  ! (V| x Sca|)~
  call computeFilter(wave,sizeTestFilter*delta,tempK,TildeK,filter)
  call btran(TildeK,temp2,res)
  if ( .not. res ) write(6,'(a)')'[ERROR] in computeGermanoSca : btran failed !' 
  ! V|~
  call computeFilter(wave,sizeTestFilter*delta,Vk,TildeK,filter)
  call btran(TildeK,temp1,res)
  if ( .not. res ) write(6,'(a)')'[ERROR] in computeGermanoSca : btran failed !' 
  ! (V| x Sca|)~ - V|~ x Sca|~
  L2%values = temp2%values - temp1%values * ScalTilde%values

  ! W| x Sca|
  temp2%values = W%values * Scal%values
  call ftran(temp2,tempK,res)
  if ( .not. res ) write(6,'(a)')'[ERROR] in computeGermanoSca : ftran failed !' 
  ! (W| x Sca|)~
  call computeFilter(wave,sizeTestFilter*delta,tempK,TildeK,filter)
  call btran(TildeK,temp2,res)
  if ( .not. res ) write(6,'(a)')'[ERROR] in computeGermanoSca : btran failed !' 
  ! W|~
  call computeFilter(wave,sizeTestFilter*delta,Wk,TildeK,filter)
  call btran(TildeK,temp1,res)
  if ( .not. res ) write(6,'(a)')'[ERROR] in computeGermanoSca : btran failed !' 
  ! ( W| x Sca|)~ - W|~ x Sca|~
  L3%values = temp2%values - temp1%values * ScalTilde%values

  call deleteDataLayout(TildeK)
  call deleteDataLayout(tempK)
  call deleteDataLayout(temp1)
  call deleteDataLayout(temp2)
  call deleteDataLayout(ScalTilde)

 end subroutine computeGermanoSca




!------------------------------------------------------------------------------
!> @author 
!> Guillaume Balarac and Antoine Vollant, LEGI
!> 
!> @details
!> Compute a recurrent term |S|dZdxi 
!> @param [in]    U is an example of realdatalayout to create work realdatalayout 
!> @param [in]    Uk is the first component of velocity in Fourier space
!> @param [in]    Vk is the second component of velocity in Fourier space
!> @param [in]    Wk is the third component of velocity in Fourier space
!> @param [in]    ScalK is the scalar in Fourier space
!> @param [inout] outSdZ1 is the fisrt component of |S|dZdxi
!> @param [inout] outSdZ2 is the second component of |S|dZdxi
!> @param [inout] outSdZ3 is the third component of |S|dZdxi
!> @param [inout] res is a logical value equal to .true. if the routine is ok .false. otherwize 
!------------------------------------------------------------------------------
 subroutine computeSdZ(U,Uk,Vk,Wk,ScalK,ScalWN,outSdZ1,outSdZ2,outSdZ3,res)

  !I/O
  type(complex_data_layout),intent(in)    :: Uk,Vk,Wk,ScalK 
  type(real_data_layout),intent(in)       :: U
  type(real_data_layout),intent(inout)    :: outSdZ1,outSdZ2,outSdZ3
  type(WaveNumbers), intent(in)           :: ScalWN
  logical, intent(inout)                  :: res

  !Local data
  type(real_data_layout),dimension(:),allocatable     :: temp 
  type(complex_data_layout),dimension(:),allocatable  :: tempK  
  
  res = .true.

  if (.not. initWorkArray(U,4,temp) ) then
    res = .false.
    write(6,'(a)')'[ERROR] in computeSdZ for scalar : not enought memory'
  endif
  if (.not. initWorkArray(Uk,3,tempK) ) then
    res = .false.
    write(6,'(a,i0)')'[ERROR] in computeSdZ for scalar : not enought memory'
  endif

  call computeStressTensor(Uk,VK,WK,ScalWN,temp(1),temp(2),temp(3),outSdZ1,outSdZ2,outSdZ3,res)
  !    computeStressTensor(Uk,VK,WK,ScalWN,S11    ,S12    ,S13    ,S22    ,S23    ,S33    ,res)
#ifdef DEBUGLES
  if (.not. res ) write(6,'(a,i0)')'[ERROR] in computeSdZ for scalar : compute StressTensor failed'
#endif
  call computeTensorNorme1(temp(1),temp(2),temp(3),outSdZ1,outSdZ2,outSdZ3,temp(4),res)
  !    computeTensorNorme1(S11    ,S12    ,S13    ,S22    ,S23    ,S33    ,norm   ,res)
#ifdef DEBUGLES
  if (.not. res ) write(6,'(a,i0)')'[ERROR] in computeSdZ for scalar : compute StressTensor failed'
#endif
  call computeGradientK(ScalWN,ScalK,tempK(1),tempK(2),tempK(3))
  call btran(tempK(1),outSdZ1,res)
#ifdef DEBUGLES
  if (.not. res ) write(6,'(a,i0)')'[ERROR] in computeSdZ for scalar : compute btran 1 failed'
#endif
  call btran(tempK(2),outSdZ2,res)
#ifdef DEBUGLES
  if (.not. res ) write(6,'(a,i0)')'[ERROR] in computeSdZ for scalar : compute btran 2 failed'
#endif
  call btran(tempK(3),outSdZ3,res)
#ifdef DEBUGLES
  if (.not. res ) write(6,'(a,i0)')'[ERROR] in computeSdZ for scalar : compute btran 3 failed'
#endif

  temp(1)%values = temp(4)%values * outSdZ1%values
  outSdZ1%values = temp(1)%values
  temp(1)%values = temp(4)%values * outSdZ2%values
  outSdZ2%values = temp(1)%values
  temp(1)%values = temp(4)%values * outSdZ3%values
  outSdZ3%values = temp(1)%values 

  if (.not. deleteWorkArray(temp) .or. &
     &.not. deleteWorkArray(tempK)) then
    res = .false.
    write(6,'(a,i0)')'[ERROR] in computeSdZ for scalar : memory leak is occuring...'
  endif

 end subroutine computeSdZ

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
!> @param [in] filterSize is the value of the filter for a priori test
!> @param [in] spec_rank is the number of process in the cpu pool
!> @param [inout] CoeffDynScalar is dynamic coefficient which is computed by the routine. Average is done all over the domain 
!> @param [inout] CoeffDynArray is dynamic coefficient which is computed by the routine. Average is done in 2D slice 
!> @return      success is a logical equal to .true. if everything is Ok and .false. otherewize.
!------------------------------------------------------------------------------
 function coefDynClarkFabreSca(U,V,W,Uk,Vk,Wk,ScalIn,ScalKIn,ScalWN,filter,filterSize,spec_rank,&
                              &CoeffDynScalar,CoeffDynArray,dir2DAvg) result(success)

  !I/O
  type(complex_data_layout),intent(in)          :: Uk,Vk,Wk,ScalKIn 
  type(real_data_layout),intent(in)             :: U,V,W,ScalIn
  type(WaveNumbers), intent(in)                 :: ScalWN
  integer, intent(in)                           :: spec_rank,filter
  integer, intent(in)                           :: filterSize
  integer, intent(in),optional                  :: dir2DAvg 
  real(WP), intent(inout),dimension(:),optional :: CoeffDynArray
  real(WP), intent(inout),optional              :: CoeffDynScalar
  logical                                       :: success

  !Local data
  real(WP)                                           :: sizeTestFilter
  type(real_data_layout),dimension(:),allocatable    :: L
  type(real_data_layout),dimension(:),allocatable    :: K 
  type(real_data_layout),dimension(:),allocatable    :: N
  type(complex_data_layout),dimension(:),allocatable :: ScalKhat
  type(complex_data_layout),dimension(:),allocatable :: UKhat
  real(WP)                                           :: delta
  logical                                            :: res 

  success = .false.

  !Init des tableaux de travail
  if( .not. initWorkArray(Uk,3,Ukhat) .or. &
  & .not. initWorkArray(U,3,L) .or. &
  & .not. initWorkArray(U,3,K) .or. &
  & .not. initWorkArray(U,3,N) .or. &
  & .not. initWorkArray(Uk,1,ScalKhat)) then
    write(6,'(a)') '[ERROR] not enought memory in coefDynClarkSca'
    return
  endif
  !valeur de la taille de filtre
  delta = real(filterSize*computeDelta(ScalIn),WP)
  !coef multiplicateur entre le filtre et le filtre test
  sizeTestFilter  = 2.0_WP
  !filtre test pour le calcul des coefs
  !Filtrage de tous les termes par le filtre 'filtre' à la longueur de coupure nd*delta
  call computeFilter(ScalWN,sizeTestFilter*delta,Uk,UKhat(1),filter)
  call computeFilter(ScalWN,sizeTestFilter*delta,Vk,UKhat(2),filter)
  call computeFilter(ScalWN,sizeTestFilter*delta,Wk,UKhat(3),filter)
  call computeFilter(ScalWN,sizeTestFilter*delta,ScalKIn,ScalKhat(1),filter)
  !calcul des Li
  call computeGermanoSca(ScalWN,U,V,W,ScalIn,Uk,Vk,Wk,ScalkIn,Filter,delta,L(1),L(2),L(3))
  !calcul des Ki
  call compute_Hi_Ki( K(1),K(2),K(3),ScalWN,delta,sizeTestFilter,filter,Uk,.false.,res,&
                         & Vk=Vk,Wk=Wk,Scalk=Scalkin,&
                         & Ukhat=Ukhat(1),Vkhat=Ukhat(2),Wkhat=Ukhat(3),Scalkhat=ScalKhat(1) )
  !calcul des Ni   
  call compute_Mi_Ni ( N(1),N(2),N(3),ScalWN,delta,2.0_WP,filter,Uk,.false.,res,&
                     & Vk=Vk,Wk=Wk,Scalk=Scalkin,&
                     & Ukhat=UKhat(1),Vkhat=UKhat(2),Wkhat=UKhat(3),Scalkhat=ScalKhat(1) )

  L(1)%values = L(1)%values-K(1)%values
  L(2)%values = L(2)%values-K(2)%values
  L(3)%values = L(3)%values-K(3)%values
  
  !Moyenne isotrope
  if (present(CoeffDynScalar)) then
    call computeLilly(N(1),N(2),N(3),L(1),L(2),L(3),CoeffDynScalar,spec_rank)
    if (spec_rank .eq. 0) write (6,'(a,1x,f9.6)') '[INFO] CoefDyn ClarkFabreModel=',CoeffDynScalar
!  Moyenne anisotrope
  elseif (present(CoeffDynArray) .and. present(dir2DAvg)) then
    call computeLilly(N(1),N(2),N(3),L(1),L(2),L(3),CoeffDynArray,spec_rank,dir2DAvg)
  else
    if (spec_rank .eq. 0) write (6,'(a)') '[ERROR] something is missing (CoeffDynScalar or &
                                          & CoeffDynArray or dir2DAvg) coef not computed'
  endif
  !Désallocation des tableaux de travail
  if ( .not. deleteWorkArray(Ukhat) .or. &
  & .not. deleteWorkArray(L) .or. &
  & .not. deleteWorkArray(K) .or. &
  & .not. deleteWorkArray(N) .or. &
  & .not. deleteWorkArray(ScalKhat)) then
    write(6,'(a)') '[ERROR] memory leakage in coefDynClarkSca'
    return
  endif

  success = .true.

 end function coefDynClarkFabreSca


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
!> @param [in] filterSize is the value of the filter for a priori test
!> @param [in] spec_rank is the number of process in the cpu pool
!> @param [inout] CoeffDynScalar is dynamic coefficient which is computed by the routine. Average is done all over the domain 
!> @param [inout] CoeffDynArray is dynamic coefficient which is computed by the routine. Average is done in 2D slice 
!> @return      success is a logical equal to .true. if everything is Ok and .false. otherewize.
!------------------------------------------------------------------------------
 function coefDynClarkSca(U,V,W,Uk,Vk,Wk,ScalIn,ScalKIn,ScalWN,filter,filterSize,spec_rank,CoeffDynScalar,&
                         &CoeffDynArray,dir2DAvg) result(success)

  !I/O
  type(complex_data_layout),intent(in)          :: Uk,Vk,Wk,ScalKIn 
  type(real_data_layout),intent(in)             :: U,V,W,ScalIn
  type(WaveNumbers), intent(in)                 :: ScalWN
  integer, intent(in)                           :: spec_rank,filter
  real(WP), intent(inout),dimension(:),optional :: CoeffDynArray
  real(WP), intent(inout),optional              :: CoeffDynScalar
  integer, intent(in),optional                  :: dir2DAvg 
  integer, intent(in)                           :: filterSize
  logical                                       :: success

  !Local data
  real(WP)                                           :: sizeTestFilter
  type(real_data_layout),dimension(:),allocatable    :: L
  type(real_data_layout),dimension(:),allocatable    :: H 
  type(real_data_layout),dimension(:),allocatable    :: M
  type(complex_data_layout),dimension(:),allocatable :: UKhat
  type(complex_data_layout),dimension(:),allocatable :: ScalKhat
  real(WP)                                           :: delta
  logical                                            :: res
 
  success = .false.

  !Init des tableaux de travail
  if (.not. initWorkArray(Uk,3,Ukhat) .or. &
  & .not. initWorkArray(U,3,L) .or. &
  & .not. initWorkArray(U,3,H) .or. &
  & .not. initWorkArray(U,3,M) .or. &
  & .not. initWorkArray(Uk,1,ScalKhat)) then
    write(6,'(a)') '[ERROR] not enought memory in coefDynClarkSca.'
    return
  endif

  !valeur de la taille de filtre
  delta = real(filterSize*computeDelta(ScalIn),WP)
  !coef multiplicateur entre le filtre et le filtre test
  sizeTestFilter  = 2.0_WP

  !filtre test pour le calcul des coefs
  !Filtrage de tous les termes par le filtre 'filtre' à la longueur de coupure nd*delta
  call computeFilter(ScalWN,sizeTestFilter*delta,Uk,UKhat(1),filter)
  call computeFilter(ScalWN,sizeTestFilter*delta,Vk,UKhat(2),filter)
  call computeFilter(ScalWN,sizeTestFilter*delta,Wk,UKhat(3),filter)
  call computeFilter(ScalWN,sizeTestFilter*delta,ScalKIn,ScalKhat(1),filter)

  !calcul des Li
  call computeGermanoSca(ScalWN,U,V,W,ScalIn,Uk,Vk,Wk,ScalkIn,filter,delta,L(1),L(2),L(3))
  !calcul des Mi   
  call compute_Mi_Ni ( M(1),M(2),M(3),ScalWN,delta,2.0_WP,filter,Uk,.true.,res,&
                     & Vk=Vk,Wk=Wk,Scalk=ScalkIn,&
                     & Ukhat=UKhat(1),Vkhat=UKhat(2),Wkhat=UKhat(3),Scalkhat=ScalKhat(1) )
  !calcul des Hi
  call compute_Hi_Ki( H(1),H(2),H(3),ScalWN,delta,sizeTestFilter,filter,Uk,.true.,res,&
                         & Vk=Vk,Wk=Wk,Scalk=Scalkin,&
                         & Ukhat=Ukhat(1),Vkhat=Ukhat(2),Wkhat=Ukhat(3),Scalkhat=ScalKhat(1) )

  L(1)%values = L(1)%values-H(1)%values
  L(2)%values = L(2)%values-H(2)%values
  L(3)%values = L(3)%values-H(3)%values
  !Moyenne isotrope
  if (present(CoeffDynScalar)) then
    call computeLilly(M(1),M(2),M(3),L(1),L(2),L(3),CoeffDynScalar,spec_rank)
    if (spec_rank .eq. 0) write (6,'(a,1x,f9.6)') '[INFO] CoefDyn ClarkModel=',CoeffDynScalar
  !Moyenne anisotrope
  elseif (present(CoeffDynArray) .and. present(dir2DAvg)) then
    call computeLilly(M(1),M(2),M(3),L(1),L(2),L(3),CoeffDynArray,spec_rank,dir2DAvg)
  else
    if (spec_rank .eq. 0) write (6,'(a)') '[ERROR] something is missing (CoeffDynScalar or &
                                          & CoeffDynArray or dir2DAvg) coef not computed'
  endif

  !Désallocation des tableaux de travail
  if( .not. deleteWorkArray(Ukhat) .or. &
  & .not. deleteWorkArray(H) .or. &
  & .not. deleteWorkArray(L) .or. &
  & .not. deleteWorkArray(M) .or. &
  & .not. deleteWorkArray(ScalKhat)) then
    write(6,'(a)') '[ERROR] memory leakage in coefDynClarkSca.'
    return
  endif

  success = .true.

 end function coefDynClarkSca



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
 function coefDynSmagSca(U,V,W,Uk,Vk,Wk,ScalIn,ScalKIn,ScalWN,filter,filterSize,spec_rank,&
                        &CoeffDynScalar,CoeffDynArray,dir2DAvg) result(success)

  !I/O
  type(complex_data_layout),intent(in)                     :: Uk,Vk,Wk,ScalKIn 
  type(real_data_layout),intent(in)                        :: U,V,W,ScalIn
  type(WaveNumbers), intent(in)                            :: ScalWN
  integer, intent(in)                                      :: spec_rank,filter
  integer,intent(in)                                       :: filterSize 
  real(WP),intent(inout),optional                          :: coeffDynScalar 
  real(WP),intent(inout),dimension(:),allocatable,optional :: coeffDynArray 
  integer, intent(in),optional                             :: dir2DAvg 
  logical                                                  :: success

  !Local data
  real(WP)                                           :: delta
  real(WP)                                           :: sizeTestFilter
  type(real_data_layout),dimension(:),allocatable    :: L
  type(real_data_layout),dimension(:),allocatable    :: M
  type(complex_data_layout),dimension(:),allocatable :: UKhat
  type(complex_data_layout),dimension(:),allocatable :: ScalKhat
  logical                                            :: res
 
  success = .false.

  !Init des tableaux de travail
  if (.not. initWorkArray(Uk,3,UKhat) .or. &
     &.not. initWorkArray(U,3,L).or. &
     &.not. initWorkArray(U,3,M) .or. &
     &.not. initWorkArray(Uk,1,ScalKhat)) then
    write(6,'(a)') '[ERROR] not enought memory in coefDynClarkSca'
    return
  endif

  sizeTestFilter  = 2.0_WP
  delta = real(filterSize*computeDelta(ScalIn),WP)
  !Filtering 
  call computeFilter(ScalWN,sizeTestFilter*delta,Uk,UKhat(1),filter) 
  call computeFilter(ScalWN,sizeTestFilter*delta,Vk,UKhat(2),filter) 
  call computeFilter(ScalWN,sizeTestFilter*delta,Wk,UKhat(3),filter) 
  call computeFilter(ScalWN,sizeTestFilter*delta,ScalKIn,ScalKhat(1),filter) 
  !compute Germano's term
  !calcul des Li
  call computeGermanoSca(ScalWN,U,V,W,ScalIn,Uk,Vk,Wk,ScalkIn,Filter,delta,L(1),L(2),L(3))
  !Compute Mi
  call compute_Mi_Ni ( M(1),M(2),M(3),ScalWN,delta,2.0_WP,filter,Uk,.true.,res,&
                     & Vk=Vk,Wk=Wk,Scalk=ScalkIn,&
                     & Ukhat=UKhat(1),Vkhat=UKhat(2),Wkhat=UKhat(3),Scalkhat=ScalKhat(1) )

  !Moyenne isotrope
  if (present(CoeffDynScalar)) then
    call computeLilly(M(1),M(2),M(3),L(1),L(2),L(3),CoeffDynScalar,spec_rank)
    if (spec_rank .eq. 0) write (6,'(a,1x,f9.6)') '[INFO] CoefDyn SmagModel=',CoeffDynScalar
  !Moyenne anisotrope
  elseif (present(CoeffDynArray) .and. present(dir2DAvg)) then
    call computeLilly(M(1),M(2),M(3),L(1),L(2),L(3),CoeffDynArray,spec_rank,dir2DAvg)
  else
    if (spec_rank .eq. 0) write (6,'(a)') '[ERROR] something is missing (CoeffDynScalar or &
                                          & CoeffDynArray or dir2DAvg) coef not computed'
  endif

  !Free memory
  if (.not. deleteWorkArray(UKhat) .or. &
  & .not. deleteWorkArray(L) .or. &
  & .not. deleteWorkArray(M) .or. &
  & .not. deleteWorkArray(ScalKhat)) then
    write (6,'(a)') '[ERROR] memory leackage in coefDynClarkSca'
    return
  endif

  success = .true.

 end function coefDynSmagSca


!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> M_i = delt^hat |S^bar^hat| d z^bar^hat /dxi - delt^bar ( |S^bar| d z^bar /dxi ) ^ hat
!> N_i = delt^hat |S^bar^hat| d z^bar^hat /dxi if rightRHS is set to .false. 
!> cf Y. Fabre and G. Balarac
!> @param [inout] M1 is the first component
!> @param [inout] M2 is the second component
!> @param [inout] M3 is the third component
!> @param [in] ScalWN is the waves number of Scalar 
!> @param [in] delta is the size of main filter 
!> @param [in] sizeTestFilter is coefficient which set the size of test filter from delta (usually 2) 
!> @param [in] filter is the type of filter 
!> @param [in] Uk is the first component of velocity in Fourier space and filtered
!> @param [in] rightRHS if you don't want to compute the second term 
!> @param [out] res is a logical value which is equal to .true. if everything is fine and .false. otherwize
!> @param [in] U is the first component of velocity in physical space and filtered (optional)
!> @param [in] V is the second component of velocity in physical space and filtered (optional)
!> @param [in] W is the third component of velocity in physical space and filtered (optional)
!> @param [in] Scal is the scalar in physical space and filtered (optional)
!> @param [in] Vk is the second component of velocity in Fourier space and filtered (optional)
!> @param [in] Wk is the third component of velocity in Fourier space and filtered (optional)
!> @param [in] Scalk is the scalar in Fourier space and filtered (optional)
!> @param [in] Ukhat the first component of velocity in Fourier space filtered at test filter (optional)
!> @param [in] Vkhat the second component of velocity in Fourier space filtered at test filter (optional)
!> @param [in] Wkhat the third component of velocity in Fourier space filtered at test filter (optional)
!> @param [in] Scalkhat the scalar in Fourier space filtered at test filter (optional)
!------------------------------------------------------------------------------
 subroutine compute_Mi_Ni( M1,M2,M3,ScalWN,delta,sizeTestFilter,filter,Uk,rightRHS,res,&
                         & U,V,W,Scal,&
                         & Vk,Wk,Scalk,&
                         & Ukhat,Vkhat,Wkhat,Scalkhat )

  !I/O data
  TYPE(real_data_layout),intent(inout)          :: M1,M2,M3
  TYPE(complex_data_layout),intent(in)          :: Uk 
  TYPE(wavenumbers),intent(in)                  :: ScalWN
  REAL(WP), intent(in)                          :: delta,sizeTestFilter
  INTEGER,intent(in)                            :: filter
  LOGICAL,intent(in)                            :: rightRHS
  TYPE(real_data_layout),intent(in),optional    :: U,V,W,Scal 
  TYPE(complex_data_layout),intent(in),optional :: Vk,Wk,Scalk,Ukhat,Vkhat,Wkhat,Scalkhat

  !Local data
  TYPE(complex_data_layout)                     :: UKhat_comp,VKhat_comp,WKhat_comp,ScalKhat_comp
  TYPE(complex_data_layout)                     :: UK_comp,VK_comp,WK_comp,ScalK_comp
  TYPE(complex_data_layout),dimension(:),allocatable ::  SdZK
  TYPE(complex_data_layout),dimension(:),allocatable :: SdZKhat
  TYPE(real_data_layout),dimension(:),allocatable ::  SdZ
  TYPE(real_data_layout),dimension(:),allocatable :: ShatdZhat 
  TYPE(real_data_layout),dimension(:),allocatable :: SdZhat
  LOGICAL                            :: res,success
  LOGICAL                            :: k_comp=.false. 
  LOGICAL                            :: khat_comp=.false. 

  if (.not. initWorkArray(M1,3,ShatdZhat) ) then
    write(6,'(a)') '[ERROR] not enought memory in computeMi'
    return
  endif

  if (.not.( present(Vk) .and. present(Wk) .and. present(Scalk)) .and. &
     &.not.( present(U) .and. present(V) .and. present(W) .and. present(Scal)))then
    write (6,'(a)') '[ERROR] not able to compute, you need to give at least (Uk,Vk,Wk,Scalk) &
                      & or (U,V,W,Scal)'
  endif
  
  if ( present(Ukhat) .and. present(Vkhat) .and. present (Wkhat) .and. present(Scalkhat) ) then
    call computeSdZ(M1,UKhat,VKhat,WKhat,ScalKhat,&
                   &ScalWN,ShatdZhat(1),ShatdZhat(2),ShatdZhat(3),res)
  elseif (present(Vk) .and. present(Wk) .and. present(Scalk)) then
    if (.not. copyStructOnly(UK,UKhat_comp) .or. &
       &.not. copyStructOnly(UK,VKhat_comp) .or. &
       &.not. copyStructOnly(UK,Wkhat_comp) .or. &
       &.not. copyStructOnly(UK,Scalkhat_comp) )then 
        success = .false.
       write(6,'(a)')'[ERROR] in computeMi: not enought memory !'
    endif 
    call computeFilter(ScalWN,sizeTestFilter*delta,Uk,UKhat_comp,filter) 
    call computeFilter(ScalWN,sizeTestFilter*delta,Vk,VKhat_comp,filter) 
    call computeFilter(ScalWN,sizeTestFilter*delta,Wk,WKhat_comp,filter) 
    call computeFilter(ScalWN,sizeTestFilter*delta,Scalk,Scalkhat_comp,filter) 
    call computeSdZ(M1,UKhat_comp,VKhat_comp,WKhat_comp,ScalKhat_comp,&
                   &ScalWN,ShatdZhat(1),ShatdZhat(2),ShatdZhat(3),res)
    khat_comp=.true.
  elseif (present(U) .and. present(V) .and. present(W) .and. present(Scal)) then
    if (.not. copyStructOnly(UK,UK_comp) .or. &
       &.not. copyStructOnly(UK,VK_comp) .or. &
       &.not. copyStructOnly(UK,Wk_comp) .or. &
       &.not. copyStructOnly(UK,Scalk_comp) )then 
        success = .false.
       write(6,'(a)')'[ERROR] in computeMi: not enought memory !'
    endif 
    call ftran(U,Uk_comp,res)
    call ftran(V,Vk_comp,res)
    call ftran(W,Wk_comp,res)
    call ftran(Scal,Scalk_comp,res)
    if (.not. copyStructOnly(UK,UKhat_comp) .or. &
       &.not. copyStructOnly(UK,VKhat_comp) .or. &
       &.not. copyStructOnly(UK,Wkhat_comp) .or. &
       &.not. copyStructOnly(UK,Scalkhat_comp) )then 
        success = .false.
       write(6,'(a)')'[ERROR] in computeMi: not enought memory !'
    endif 
    call computeFilter(ScalWN,sizeTestFilter*delta,Uk,UKhat_comp,filter) 
    call computeFilter(ScalWN,sizeTestFilter*delta,Vk,VKhat_comp,filter) 
    call computeFilter(ScalWN,sizeTestFilter*delta,Wk,WKhat_comp,filter) 
    call computeFilter(ScalWN,sizeTestFilter*delta,Scalk,Scalkhat_comp,filter) 
    call computeSdZ(M1,UKhat_comp,VKhat_comp,WKhat_comp,ScalKhat_comp,&
                   &ScalWN,ShatdZhat(1),ShatdZhat(2),ShatdZhat(3),res)
    k_comp=.true.
    khat_comp=.true.
  endif
  M1%values = delta**2.0 * sizeTestFilter**2.0 * ShatdZhat(1)%values 
  M2%values = delta**2.0 * sizeTestFilter**2.0 * ShatdZhat(2)%values 
  M3%values = delta**2.0 * sizeTestFilter**2.0 * ShatdZhat(3)%values 

  if ( rightRHS ) then
    !Compute ShatdZhat_i
    if ( .not. initWorkArray(M1,3,SdZ) .or. &
       &.not. initWorkArray(Uk,3,SdZK) .or. &
       &.not. initWorkArray(M1,3,SdZhat) .or. &
       &.not. initWorkArray(Uk,3,SdZKhat)) then
      write(6,'(a)') '[ERROR] not enought memory in computeMi'
      return
    endif
    if ( present(Vk) .and. present(Wk) .and. present(Scalk)) then
      call computeSdZ(M1,Uk,Vk,Wk,ScalK,ScalWN,SdZ(1),SdZ(2),SdZ(3),res)
    elseif ( present(U) .and. present(V) .and. present(W) .and. present(Scal)) then
      call computeSdZ(M1,Uk_comp,Vk_comp,Wk_comp,ScalK_comp,ScalWN,SdZ(1),SdZ(2),SdZ(3),res)
    endif
    call ftran(SdZ(1),SdZK(1),res)
    call ftran(SdZ(2),SdZK(2),res)
    call ftran(SdZ(3),SdZK(3),res)
    call computeFilter(ScalWN,sizeTestFilter*delta,SdZK(1),SdZKhat(1),filter) 
    call computeFilter(ScalWN,sizeTestFilter*delta,SdZK(2),SdZKhat(2),filter) 
    call computeFilter(ScalWN,sizeTestFilter*delta,SdZK(3),SdZKhat(3),filter) 
    call btran(SdZKhat(1),SdZhat(1),res)
    call btran(SdZKhat(2),SdZhat(2),res)
    call btran(SdZKhat(3),SdZhat(3),res)
    M1%values = M1%values - delta**2 * SdZhat(1)%values 
    M2%values = M2%values - delta**2 * SdZhat(2)%values
    M3%values = M3%values - delta**2 * SdZhat(3)%values 
  endif

  if (.not. deleteWorkArray(ShatdZhat)) then
    write (6,'(a)') '[ERROR]'
  endif
  if (k_comp) then
    call deletedatalayout(Uk_comp)
    call deletedatalayout(Vk_comp)
    call deletedatalayout(Wk_comp)
    call deletedatalayout(Scalk_comp)
  endif
  if (khat_comp) then
    call deletedatalayout(Ukhat_comp)
    call deletedatalayout(Vkhat_comp)
    call deletedatalayout(Wkhat_comp)
    call deletedatalayout(Scalkhat_comp)
  endif
  if ( rightRHS ) then
    if (.not. deleteWorkArray(SdZ) .or. &
       &.not. deleteWorkArray(SdZK) .or. &
       &.not. deleteWorkArray(SdZhat) .or. &
       &.not. deleteWorkArray(SdZKhat)) then
      write (6,'(a)') '[ERROR]'
    endif
  endif

 end subroutine compute_Mi_Ni



!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> compute those terms:
!> H_i = delt^hat/12 d(U_i^bar^hat)/dxj d(z^bar^hat)/dxj - delt^bar/12 ( d(U_i^bar)/dxj d(z^bar)/dxj )^hat
!> K_i = delt^hat/12 d(U_i^bar^hat)/dxj d(z^bar^hat)/dxj if rightRHS set to .false. 
!> cf Y. Fabre and G. Balarac
!> @param [inout] H1 is the first component
!> @param [inout] H2 is the second component
!> @param [inout] H3 is the third component
!> @param [in] ScalWN is the waves number of Scalar 
!> @param [in] delta is the size of main filter 
!> @param [in] sizeTestFilter is coefficient which set the size of test filter from delta (usually 2) 
!> @param [in] filter is the type of filter 
!> @param [in] Uk is the first component of velocity in Fourier space and filtered
!> @param [in] rightRHS if you don't want to compute the second term 
!> @param [out] res is a logical value which is equal to .true. if everything is fine and .false. otherwize
!> @param [in] U is the first component of velocity in physical space and filtered (optional)
!> @param [in] V is the second component of velocity in physical space and filtered (optional)
!> @param [in] W is the third component of velocity in physical space and filtered (optional)
!> @param [in] Scal is the scalar in physical space and filtered (optional)
!> @param [in] Vk is the second component of velocity in Fourier space and filtered (optional)
!> @param [in] Wk is the third component of velocity in Fourier space and filtered (optional)
!> @param [in] Scalk is the scalar in Fourier space and filtered (optional)
!> @param [in] Ukhat the first component of velocity in Fourier space filtered at test filter (optional)
!> @param [in] Vkhat the second component of velocity in Fourier space filtered at test filter (optional)
!> @param [in] Wkhat the third component of velocity in Fourier space filtered at test filter (optional)
!> @param [in] Scalkhat the scalar in Fourier space filtered at test filter (optional)
! ------------------------------------------------------------------------------

 subroutine compute_Hi_Ki( H1,H2,H3,ScalWN,delta,sizeTestFilter,filter,Uk,rightRHS,res,&
                         & U,V,W,Scal,&
                         & Vk,Wk,Scalk,&
                         & Ukhat,Vkhat,Wkhat,Scalkhat )

  !I/O data
  TYPE(real_data_layout),intent(inout)          :: H1,H2,H3
  TYPE(complex_data_layout),intent(in)          :: Uk 
  TYPE(wavenumbers),intent(in)                  :: ScalWN
  REAL(WP), intent(in)                          :: delta,sizeTestFilter
  INTEGER,intent(in)                            :: filter
  LOGICAL,intent(in)                            :: rightRHS
  TYPE(real_data_layout),intent(in),optional    :: U,V,W,Scal 
  TYPE(complex_data_layout),intent(in),optional :: Vk,Wk,Scalk,Ukhat,Vkhat,Wkhat,Scalkhat

  !Local data
  TYPE(complex_data_layout)                     :: UKhat_comp,VKhat_comp,WKhat_comp,ScalKhat_comp
  TYPE(complex_data_layout)                     :: UK_comp,VK_comp,WK_comp,ScalK_comp
  TYPE(complex_data_layout),dimension(:),allocatable ::  H1K
  TYPE(complex_data_layout),dimension(:),allocatable ::  H2K
  TYPE(complex_data_layout),dimension(:),allocatable :: tempk 
  TYPE(real_data_layout),dimension(:),allocatable ::  temp
  TYPE(real_data_layout),dimension(:),allocatable :: H 
  LOGICAL                            :: res,success
  LOGICAL                            :: k_comp=.false. 
  LOGICAL                            :: khat_comp=.false. 

  if (.not. initWorkArray(H1,12,temp) .or. & 
     &.not. initWorkArray(Vk,12,tempK) ) then
    write(6,'(a)') '[ERROR] not enought memory in computeMi'
    return
  endif
  
 
  if ( present(Ukhat) .and. present(Vkhat) .and. present (Wkhat) .and. present(Scalkhat) ) then
    call computeGradientK(ScalWN,UKhat,tempK(1),tempK(2),tempK(3))
    call computeGradientK(ScalWN,VKhat,tempK(4),tempK(5),tempK(6))
    call computeGradientK(ScalWN,WKhat,tempK(7),tempK(8),tempK(9))
    call computeGradientK(ScalWN,ScalKhat,tempK(10),tempK(11),tempK(12))
  elseif ( present(Vk) .and. present(Wk) .and. present(Scalk)) then
    if (.not. copyStructOnly(UK,Ukhat_comp) .or. &
       &.not. copyStructOnly(UK,Vkhat_comp) .or. &
       &.not. copyStructOnly(UK,Wkhat_comp) .or. &
       &.not. copyStructOnly(UK,Scalkhat_comp) )then 
        success = .false.
       write(6,'(a)')'[ERROR] in computeMi: not enought memory !'
    endif 
    call computeFilter(ScalWN,sizeTestFilter*delta,Uk,Ukhat_comp,filter) 
    call computeFilter(ScalWN,sizeTestFilter*delta,Vk,Vkhat_comp,filter) 
    call computeFilter(ScalWN,sizeTestFilter*delta,Wk,Wkhat_comp,filter) 
    call computeFilter(ScalWN,sizeTestFilter*delta,Scalk,Scalkhat_comp,filter) 
    call computeGradientK(ScalWN,UKhat_comp,tempK(1),tempK(2),tempK(3))
    call computeGradientK(ScalWN,VKhat_comp,tempK(4),tempK(5),tempK(6))
    call computeGradientK(ScalWN,WKhat_comp,tempK(7),tempK(8),tempK(9))
    call computeGradientK(ScalWN,ScalKhat_comp,tempK(10),tempK(11),tempK(12))
    khat_comp = .true.
  elseif (present(U) .and. present(V) .and. present (W) .and. present(Scal) ) then
    if (.not. copyStructOnly(UK,Ukhat_comp) .or. &
       &.not. copyStructOnly(UK,UK_comp) .or. &
       &.not. copyStructOnly(UK,VK_comp) .or. &
       &.not. copyStructOnly(UK,WK_comp) .or. &
       &.not. copyStructOnly(UK,ScalK_comp) .or. &
       &.not. copyStructOnly(UK,Vkhat_comp) .or. &
       &.not. copyStructOnly(UK,Wkhat_comp) .or. &
       &.not. copyStructOnly(UK,Scalkhat_comp) )then 
        success = .false.
       write(6,'(a)')'[ERROR] in computeMi: not enought memory !'
    endif 
    call ftran(U,Uk_comp,res)
    call ftran(V,Vk_comp,res)
    call ftran(W,Wk_comp,res)
    call ftran(Scal,Scalk_comp,res)
    call computeFilter(ScalWN,sizeTestFilter*delta,Uk_comp,Ukhat_comp,filter) 
    call computeFilter(ScalWN,sizeTestFilter*delta,Vk_comp,Vkhat_comp,filter) 
    call computeFilter(ScalWN,sizeTestFilter*delta,Wk_comp,Wkhat_comp,filter) 
    call computeFilter(ScalWN,sizeTestFilter*delta,Scalk_comp,Scalkhat_comp,filter) 
    call computeGradientK(ScalWN,UKhat_comp,tempK(1),tempK(2),tempK(3))
    call computeGradientK(ScalWN,VKhat_comp,tempK(4),tempK(5),tempK(6))
    call computeGradientK(ScalWN,WKhat_comp,tempK(7),tempK(8),tempK(9))
    call computeGradientK(ScalWN,ScalKhat_comp,tempK(10),tempK(11),tempK(12))
    k_comp=.true.
    khat_comp=.true.
  endif
  call btran(tempK(1),temp(1),res) 
  call btran(tempK(2),temp(2),res) 
  call btran(tempK(3),temp(3),res) 
  call btran(tempK(4),temp(4),res) 
  call btran(tempK(5),temp(5),res) 
  call btran(tempK(6),temp(6),res) 
  call btran(tempK(7),temp(7),res) 
  call btran(tempK(8),temp(8),res) 
  call btran(tempK(9),temp(9),res) 
  call btran(tempK(10),temp(10),res) 
  call btran(tempK(11),temp(11),res) 
  call btran(tempK(12),temp(12),res) 
  H1%values = (sizeTestFilter*delta)**2/12.0_WP   * &
                 &(temp(1)%values * temp(10)%values + &
                 & temp(2)%values * temp(11)%values + &
                 & temp(3)%values * temp(12)%values)
  H2%values = (sizeTestFilter*delta)**2/12.0_WP   * & 
                 &(temp(4)%values * temp(10)%values + &
                 & temp(5)%values * temp(11)%values + &
                 & temp(6)%values * temp(12)%values)
  H3%values = (sizeTestFilter*delta)**2/12.0_WP   * & 
                 &(temp(7)%values * temp(10)%values + &
                 & temp(8)%values * temp(11)%values + &
                 & temp(9)%values * temp(12)%values)
  if (rightRHS) then
    if ( .not. initWorkArray(H1,3,H) .or. &
       &.not. initWorkArray(Uk,3,H1K) .or. &
       &.not. initWorkArray(Uk,3,H2K)) then
      write(6,'(a)') '[ERROR] not enought memory in computeMi'
      return
    endif
    if ( present(Vk) .and. present(Wk) .and. present(Scalk) ) then
      call computeGradientK(ScalWN,Uk,tempK(1),tempK(2),tempK(3))
      call computeGradientK(ScalWN,Vk,tempK(4),tempK(5),tempK(6))
      call computeGradientK(ScalWN,Wk,tempK(7),tempK(8),tempK(9))
      call computeGradientK(ScalWN,ScalK,tempK(10),tempK(11),tempK(12))
    elseif (present(U) .and. present(V) .and. present (W) .and. present(Scal)) then
      if (k_comp) then
        call computeGradientK(ScalWN,Uk_comp,tempK(1),tempK(2),tempK(3))
        call computeGradientK(ScalWN,Vk_comp,tempK(4),tempK(5),tempK(6))
        call computeGradientK(ScalWN,Wk_comp,tempK(7),tempK(8),tempK(9))
        call computeGradientK(ScalWN,ScalK_comp,tempK(10),tempK(11),tempK(12))
      else
        call ftran(U,Uk_comp,res)
        call ftran(V,Vk_comp,res)
        call ftran(W,Wk_comp,res)
        call ftran(Scal,Scalk_comp,res)
        call computeGradientK(ScalWN,Uk_comp,tempK(1),tempK(2),tempK(3))
        call computeGradientK(ScalWN,Vk_comp,tempK(4),tempK(5),tempK(6))
        call computeGradientK(ScalWN,Wk_comp,tempK(7),tempK(8),tempK(9))
        call computeGradientK(ScalWN,ScalK_comp,tempK(10),tempK(11),tempK(12))
      endif
    endif
    call btran(tempK(1),temp(1),res) 
    call btran(tempK(2),temp(2),res) 
    call btran(tempK(3),temp(3),res) 
    call btran(tempK(4),temp(4),res) 
    call btran(tempK(5),temp(5),res) 
    call btran(tempK(6),temp(6),res) 
    call btran(tempK(7),temp(7),res) 
    call btran(tempK(8),temp(8),res) 
    call btran(tempK(9),temp(9),res) 
    call btran(tempK(10),temp(10),res) 
    call btran(tempK(11),temp(11),res) 
    call btran(tempK(12),temp(12),res) 
    H(1)%values = (temp(1)%values * temp(10)%values + &
                   &temp(2)%values * temp(11)%values + &
                   &temp(3)%values * temp(12)%values)
    H(2)%values = (temp(4)%values * temp(10)%values + &
                   &temp(5)%values * temp(11)%values + &
                   &temp(6)%values * temp(12)%values)
    H(3)%values = (temp(7)%values * temp(10)%values + &
                   &temp(8)%values * temp(11)%values + &
                   &temp(9)%values * temp(12)%values)
    call ftran(H(1),H1K(1),res)
    call ftran(H(2),H1K(2),res)
    call ftran(H(3),H1K(3),res)
    call computeFilter(ScalWN,sizeTestFilter*delta,H1K(1),H2K(1),filter)
    call computeFilter(ScalWN,sizeTestFilter*delta,H1K(2),H2K(2),filter)
    call computeFilter(ScalWN,sizeTestFilter*delta,H1K(3),H2K(3),filter)
    call btran(H2K(1),H(1),res) 
    call btran(H2K(2),H(2),res) 
    call btran(H2K(3),H(3),res)  
  
    H1%values = H1%values - delta**2/12.0_WP*H(1)%values
    H2%values = H2%values - delta**2/12.0_WP*H(2)%values
    H3%values = H3%values - delta**2/12.0_WP*H(3)%values
  endif

  if (.not. deleteWorkArray(temp) .or. &
     &.not. deleteWorkArray(tempk)) then
    write (6,'(a)') '[ERROR] in computeMi : memory leackage'
  endif
  if (k_comp) then
    call deletedatalayout(Uk_comp)
    call deletedatalayout(Vk_comp)
    call deletedatalayout(Wk_comp)
    call deletedatalayout(Scalk_comp)
  endif
  if (khat_comp) then
    call deletedatalayout(Ukhat_comp)
    call deletedatalayout(Vkhat_comp)
    call deletedatalayout(Wkhat_comp)
    call deletedatalayout(Scalkhat_comp)
  endif
  if ( rightRHS ) then
    if (.not. deleteWorkArray(H) .or. &
       &.not. deleteWorkArray(H1K) .or. &
       &.not. deleteWorkArray(H2K)) then
      write (6,'(a)') '[ERROR]'
    endif
  endif


 end subroutine compute_Hi_Ki


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
!> @param [in] filterSize is the value of the filter for a priori test
!> @param [in] spec_rank is the number of process in the cpu pool
!> @param [inout] CoeffDynScalar is dynamic coefficient which is computed by the routine. Average is done all over the domain 
!> @param [inout] CoeffDynArray is dynamic coefficient which is computed by the routine. Average is done in 2D slice 
!> @return      success is a logical equal to .true. if everything is Ok and .false. otherewize.
!------------------------------------------------------------------------------
 function coefDynClarkScaTEST(U,V,W,Uk,Vk,Wk,ScalIn,ScalKIn,ScalWN,filter,filterSize,spec_rank,CoeffDynScalar,CoeffDynArray) result(success)

  !I/O
  type(complex_data_layout),intent(in)          :: Uk,Vk,Wk,ScalKIn 
  type(real_data_layout),intent(in)             :: U,V,W,ScalIn
  type(WaveNumbers), intent(in)                 :: ScalWN
  integer, intent(in)                           :: spec_rank,filter
  real(WP), intent(inout),dimension(:),optional :: CoeffDynArray
  real(WP), intent(inout),optional              :: CoeffDynScalar
  integer, intent(in)                           :: filterSize
  logical                                       :: success

  !Local data
  real(WP)                                           :: sizeTestFilter
  type(real_data_layout),dimension(:),allocatable    :: Uhat
  type(real_data_layout),dimension(:),allocatable    :: L
  type(real_data_layout),dimension(:),allocatable    :: UZhat
  type(real_data_layout),dimension(:),allocatable    :: UZ
  type(real_data_layout),dimension(:),allocatable    :: H 
  type(real_data_layout),dimension(:),allocatable    :: H_test 
  type(real_data_layout),dimension(:),allocatable    :: H1
  type(real_data_layout),dimension(:),allocatable    :: M
  type(real_data_layout),dimension(:),allocatable    :: M_test
  type(real_data_layout),dimension(:),allocatable    :: Scalhat
  type(real_data_layout),dimension(:),allocatable    :: temp
  type(real_data_layout),dimension(:),allocatable    :: numer
  type(real_data_layout),dimension(:),allocatable    :: denom
  type(real_data_layout),dimension(:),allocatable    :: SdZ
  type(real_data_layout),dimension(:),allocatable    :: SdZhat
  type(real_data_layout),dimension(:),allocatable    :: ShatdZhat 
  type(complex_data_layout),dimension(:),allocatable :: SdZK
  type(complex_data_layout),dimension(:),allocatable :: SdZKhat
  type(complex_data_layout),dimension(:),allocatable :: tempK
  type(complex_data_layout),dimension(:),allocatable :: H1K
  type(complex_data_layout),dimension(:),allocatable :: H2K
  type(complex_data_layout),dimension(:),allocatable :: UZK
  type(complex_data_layout),dimension(:),allocatable :: UZKhat
  type(complex_data_layout),dimension(:),allocatable :: UKhat
  type(complex_data_layout),dimension(:),allocatable :: ScalKhat
  real(WP)                                           :: delta
  real(WP),dimension(:),allocatable                  :: num1D
  real(WP),dimension(:),allocatable                  :: den1D
  integer                                            :: ires
  logical                                            :: res

  logical :: erreur 
  integer :: i,j,k
  real(WP) :: Lilly

 
  success = .false.

  !Init des tableaux de travail
  if (.not. initWorkArray(U,3,Uhat) .or. &
  & .not. initWorkArray(U,3,L) .or. &
  & .not. initWorkArray(Uk,12,tempK) .or. & 
  & .not. initWorkArray(U,12,temp) .or. &
  & .not. initWorkArray(U,3,H) .or. &
  & .not. initWorkArray(U,3,H_test) .or. &
  & .not. initWorkArray(U,3,H1) .or. &
  & .not. initWorkArray(U,3,UZ) .or. &
  & .not. initWorkArray(U,3,ShatdZhat) .or. &
  & .not. initWorkArray(U,3,SdZhat) .or. &
  & .not. initWorkArray(U,3,SdZ) .or. &
  & .not. initWorkArray(Uk,3,SdZK) .or. &
  & .not. initWorkArray(Uk,3,SdZKhat) .or. &
  & .not. initWorkArray(U,3,UZhat) .or. &
  & .not. initWorkArray(U,3,M) .or. &
  & .not. initWorkArray(U,3,M_test) .or. &
  & .not. initWorkArray(U,1,Scalhat) .or. &
  & .not. initWorkArray(U,1,numer) .or. &
  & .not. initWorkArray(U,1,denom) .or. &
  & .not. initWorkArray(Uk,3,H1K) .or. &
  & .not. initWorkArray(Uk,3,H2K) .or. &
  & .not. initWorkArray(Uk,3,UZK) .or. &
  & .not. initWorkArray(Uk,3,UZKhat) .or. &
  & .not. initWorkArray(Uk,3,UKhat) .or. &
  & .not. initWorkArray(Uk,1,ScalKhat)) then
    write(6,'(a)') '[ERROR] not enought memory in coefDynClarkSca.'
    return
  endif

  !valeur de la taille de filtre
  delta = real(filterSize*computeDelta(ScalIn),WP)
  !coef multiplicateur entre le filtre et le filtre test
  sizeTestFilter  = 2.0_WP

  !filtre test pour le calcul des coefs
  !Filtrage de tous les termes par le filtre 'filtre' à la longueur de coupure nd*delta
  call computeFilter(ScalWN,sizeTestFilter*delta,Uk,UKhat(1),filter)
  call computeFilter(ScalWN,sizeTestFilter*delta,Vk,UKhat(2),filter)
  call computeFilter(ScalWN,sizeTestFilter*delta,Wk,UKhat(3),filter)
  call computeFilter(ScalWN,sizeTestFilter*delta,ScalKIn,ScalKhat(1),filter)
  !Retour espace physique des quantitées filtrées
  call btran(UKhat(1),Uhat(1),res) 
  call btran(UKhat(2),Uhat(2),res) 
  call btran(UKhat(3),Uhat(3),res)  
  call btran(ScalKhat(1),Scalhat(1),res)  !Zhatbar

  !calcul des Li
  UZ(1)%values = U%values * ScalIn%values
  UZ(2)%values = V%values * ScalIn%values
  UZ(3)%values = W%values * ScalIn%values
  call ftran(UZ(1),UZK(1),res) 
  call ftran(UZ(2),UZK(2),res) 
  call ftran(UZ(3),UZK(3),res)  
  call computeFilter(ScalWN,sizeTestFilter*delta,UZK(1),UZKhat(1),filter)
  call computeFilter(ScalWN,sizeTestFilter*delta,UZK(2),UZKhat(2),filter)
  call computeFilter(ScalWN,sizeTestFilter*delta,UZK(3),UZKhat(3),filter)
  call btran(UZKhat(1),UZhat(1),res) 
  call btran(UZKhat(2),UZhat(2),res) 
  call btran(UZKhat(3),UZhat(3),res)  
  L(1)%values=UZhat(1)%values - Scalhat(1)%values * Uhat(1)%values 
  L(2)%values=UZhat(2)%values - Scalhat(1)%values * Uhat(2)%values
  L(3)%values=UZhat(3)%values - Scalhat(1)%values * Uhat(3)%values

  !calcul des Mi   
  !Compute ShatdZhat_i
  call computeSdZ(U,UKhat(1),UKhat(2),UKhat(3),ScalKhat(1),&
                 &ScalWN,ShatdZhat(1),ShatdZhat(2),ShatdZhat(3),res)
  !Compute (SdZ)hat_i
  call computeSdZ(U,Uk,Vk,Wk,ScalKIn,ScalWN,SdZ(1),SdZ(2),SdZ(3),res)
  call ftran(SdZ(1),SdZK(1),res)
  call ftran(SdZ(2),SdZK(2),res)
  call ftran(SdZ(3),SdZK(3),res)
  call computeFilter(ScalWN,sizeTestFilter*delta,SdZK(1),SdZKhat(1),filter) 
  call computeFilter(ScalWN,sizeTestFilter*delta,SdZK(2),SdZKhat(2),filter) 
  call computeFilter(ScalWN,sizeTestFilter*delta,SdZK(3),SdZKhat(3),filter) 
  call btran(SdZKhat(1),SdZhat(1),res)
  call btran(SdZKhat(2),SdZhat(2),res)
  call btran(SdZKhat(3),SdZhat(3),res)
  M(1)%values = delta**2.0 *( sizeTestFilter**2.0 * ShatdZhat(1)%values - SdZhat(1)%values )
  M(2)%values = delta**2.0 *( sizeTestFilter**2.0 * ShatdZhat(2)%values - SdZhat(2)%values )
  M(3)%values = delta**2.0 *( sizeTestFilter**2.0 * ShatdZhat(3)%values - SdZhat(3)%values )

!!!!!!!!!!!!!!!!!!!!!!!!!!!TESTS!!!!!!!!!!!!!!!!!!!!!!!!  
  call compute_Mi_Ni ( M_test(1),M_test(2),M_test(3),ScalWN,delta,2.0_WP,filter,Uk,.true.,res,&
                     & Vk=Vk,Wk=Wk,Scalk=Scalkin,&
                     & Ukhat=UKhat(1),Vkhat=UKhat(2),Wkhat=UKhat(3),Scalkhat=ScalKhat(1) )
  erreur = .false.
  do k = U%zmin,U%zmax
  do j = U%ymin,U%ymax
  do i = U%xmin,U%xmax
    if ( .not. equalToVal (M_test(1)%values(i,j,k) ,M(1)%values(i,j,k)) .or. &
       & .not. equalToVal (M_test(2)%values(i,j,k) ,M(2)%values(i,j,k)) .or. &
       & .not. equalToVal (M_test(3)%values(i,j,k) ,M(3)%values(i,j,k)) ) then
      if (spec_rank .eq. 0 .and. .not. erreur ) print *,'problèmes dans M_test'
      erreur = .true.
    endif
  enddo
  enddo
  enddo

  call compute_Mi_Ni ( M_test(1),M_test(2),M_test(3),ScalWN,delta,2.0_WP,filter,Uk,.true.,res,&
                     & Vk=Vk,Wk=Wk,Scalk=Scalkin)
                     
  erreur = .false.
  do k = U%zmin,U%zmax
  do j = U%ymin,U%ymax
  do i = U%xmin,U%xmax
    if ( .not. equalToVal (M_test(1)%values(i,j,k) ,M(1)%values(i,j,k)) .or. &
       & .not. equalToVal (M_test(2)%values(i,j,k) ,M(2)%values(i,j,k)) .or. &
       & .not. equalToVal (M_test(3)%values(i,j,k) ,M(3)%values(i,j,k))) then
      if (spec_rank .eq. 0 .and. .not. erreur ) print *,'problèmes dans M_test'
      erreur = .true.
    endif
  enddo
  enddo
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!TESTS!!!!!!!!!!!!!!!!!!!!!!!!  

  !calcul des Hi
  call computeGradientK(ScalWN,UKhat(1),tempK(1),tempK(2),tempK(3))
  call computeGradientK(ScalWN,UKhat(2),tempK(4),tempK(5),tempK(6))
  call computeGradientK(ScalWN,UKhat(3),tempK(7),tempK(8),tempK(9))
  call computeGradientK(ScalWN,ScalKhat(1),tempK(10),tempK(11),tempK(12))
  call btran(tempK(1),temp(1),res) 
  call btran(tempK(2),temp(2),res) 
  call btran(tempK(3),temp(3),res) 
  call btran(tempK(4),temp(4),res) 
  call btran(tempK(5),temp(5),res) 
  call btran(tempK(6),temp(6),res) 
  call btran(tempK(7),temp(7),res) 
  call btran(tempK(8),temp(8),res) 
  call btran(tempK(9),temp(9),res) 
  call btran(tempK(10),temp(10),res) 
  call btran(tempK(11),temp(11),res) 
  call btran(tempK(12),temp(12),res) 
  H(1)%values =  (temp(1)%values * temp(10)%values + &
                 &temp(2)%values * temp(11)%values + &
                 &temp(3)%values * temp(12)%values)
  H(2)%values =  (temp(4)%values * temp(10)%values + &
                 &temp(5)%values * temp(11)%values + &
                 &temp(6)%values * temp(12)%values)
  H(3)%values =  (temp(7)%values * temp(10)%values + &
                 &temp(8)%values * temp(11)%values + &
                 &temp(9)%values * temp(12)%values)
  call computeGradientK(ScalWN,Uk,tempK(1),tempK(2),tempK(3))
  call computeGradientK(ScalWN,Vk,tempK(4),tempK(5),tempK(6))
  call computeGradientK(ScalWN,Wk,tempK(7),tempK(8),tempK(9))
  call computeGradientK(ScalWN,ScalKIn,tempK(10),tempK(11),tempK(12))
  call btran(tempK(1),temp(1),res) 
  call btran(tempK(2),temp(2),res) 
  call btran(tempK(3),temp(3),res) 
  call btran(tempK(4),temp(4),res) 
  call btran(tempK(5),temp(5),res) 
  call btran(tempK(6),temp(6),res) 
  call btran(tempK(7),temp(7),res) 
  call btran(tempK(8),temp(8),res) 
  call btran(tempK(9),temp(9),res) 
  call btran(tempK(10),temp(10),res) 
  call btran(tempK(11),temp(11),res) 
  call btran(tempK(12),temp(12),res) 
  H1(1)%values =  (temp(1)%values * temp(10)%values + &
                  &temp(2)%values * temp(11)%values + &
                  &temp(3)%values * temp(12)%values)
  H1(2)%values =  (temp(4)%values * temp(10)%values + &
                  &temp(5)%values * temp(11)%values + &
                  &temp(6)%values * temp(12)%values)
  H1(3)%values =  (temp(7)%values * temp(10)%values + &
                  &temp(8)%values * temp(11)%values + &
                  &temp(9)%values * temp(12)%values)
  call ftran(H1(1),H1K(1),res)
  call ftran(H1(2),H1K(2),res)
  call ftran(H1(3),H1K(3),res)
  call computeFilter(ScalWN,sizeTestFilter*delta,H1K(1),H2K(1),filter)
  call computeFilter(ScalWN,sizeTestFilter*delta,H1K(2),H2K(2),filter)
  call computeFilter(ScalWN,sizeTestFilter*delta,H1K(3),H2K(3),filter)
  call btran(H2K(1),H1(1),res) 
  call btran(H2K(2),H1(2),res) 
  call btran(H2K(3),H1(3),res)  
  H(1)%values = (sizeTestFilter*delta)**2/12.0_WP*H(1)%values - (delta)**2/12.0_WP*H1(1)%values
  H(2)%values = (sizeTestFilter*delta)**2/12.0_WP*H(2)%values - (delta)**2/12.0_WP*H1(2)%values
  H(3)%values = (sizeTestFilter*delta)**2/12.0_WP*H(3)%values - (delta)**2/12.0_WP*H1(3)%values

!!!!!!!!!!!!!!!!!!!!!!!!!!!TESTS!!!!!!!!!!!!!!!!!!!!!!!!  
  call compute_Hi_Ki( H_test(1),H_test(2),H_test(3),ScalWN,delta,sizeTestFilter,filter,Uk,.true.,res,&
                         & Vk=Vk,Wk=Wk,Scalk=Scalkin,&
                         & Ukhat=Ukhat(1),Vkhat=Ukhat(2),Wkhat=Ukhat(3),Scalkhat=ScalKhat(1) )
  erreur = .false.
  do k = U%zmin,U%zmax
  do j = U%ymin,U%ymax
  do i = U%xmin,U%xmax
    if ( .not. equalToVal (H_test(1)%values(i,j,k) ,H(1)%values(i,j,k)) .or. &
       & .not. equalToVal (H_test(2)%values(i,j,k) ,H(2)%values(i,j,k)) .or. &
       & .not. equalToVal (H_test(3)%values(i,j,k) ,H(3)%values(i,j,k)) ) then
      if (spec_rank .eq. 0 .and. .not. erreur ) print *,'problèmes dans M_test'
      erreur = .true.
    endif
  enddo
  enddo
  enddo

  call compute_Hi_Ki ( H_test(1),H_test(2),H_test(3),ScalWN,delta,2.0_WP,filter,Uk,.true.,res,&
                     & Vk=Vk,Wk=Wk,Scalk=Scalkin)
                     
  erreur = .false.
  do k = U%zmin,U%zmax
  do j = U%ymin,U%ymax
  do i = U%xmin,U%xmax
    if ( .not. equalToVal (H_test(1)%values(i,j,k) ,H(1)%values(i,j,k)) .or. &
       & .not. equalToVal (H_test(2)%values(i,j,k) ,H(2)%values(i,j,k)) .or. &
       & .not. equalToVal (H_test(3)%values(i,j,k) ,H(3)%values(i,j,k))) then
      if (spec_rank .eq. 0 .and. .not. erreur ) print *,'problèmes dans H_test'
      erreur = .true.
    endif
  enddo
  enddo
  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!TESTS!!!!!!!!!!!!!!!!!!!!!!!!  

  !Somme numérateur et dénominateur
  numer(1)%values = (L(1)%values-H(1)%values)*M(1)%values + (L(2)%values-H(2)%values)*M(2)%values +&
                  & (L(3)%values-H(3)%values)*M(3)%values
  denom(1)%values = M(1)%values*M(1)%values + M(2)%values*M(2)%values + M(3)%values*M(3)%values

  L(1)%values = L(1)%values - H(1)%values
  L(2)%values = L(2)%values - H(2)%values
  L(3)%values = L(3)%values - H(3)%values

  !Moyenne isotrope
  if (present(CoeffDynScalar)) then
    CoeffDynScalar = computeFieldAvg(numer(1),spec_rank) / computeFieldAvg(denom(1),spec_rank)
    if (spec_rank .eq. 0) write (6,'(a,1x,f9.6)') '[INFO] CoefDyn ClarkModel TEST=',CoeffDynScalar
    call computeLilly(M(1),M(2),M(3),L(1),L(2),L(3),Lilly,spec_rank)
    if ( .not. equalToVal (Lilly ,CoeffDynScalar )) then
      if (spec_rank .eq. 0) write (6,'(a)') '[ERROR] in coefDyn ClarkModel'
    endif
  endif
!  Moyenne anisotrope
  if (present(CoeffDynArray)) then
    if (.not.( (lbound(CoeffDynArray,1) .eq. 1 ) .and. (ubound(CoeffDynArray,1) .eq. U%nz) ) )then
       write (6,'(a)') '[ERROR] CoeffDynArray not declared between good bounds'
       return
    endif
    allocate(num1D(U%nz),stat=ires)
    if ( ires .ne. 0 ) then
      write (6,'(a)') '[ERROR] Not able to allocate num1D in coefDynClarkSca'
      return
    endif
    allocate(den1D(U%nz),stat=ires)
    if ( ires .ne. 0 ) then
      write (6,'(a)') '[ERROR] Not able to allocate den1D in coefDynClarkSca'
      return
    endif
    if (.not. computeFieldAvg2D(numer(1),num1D,spec_rank)) then
      write(*,'(a)')'[ERROR] Could not compute AVG for temp(4)'
      return
    endif
    if (.not. computeFieldAvg2D(denom(1),den1D,spec_rank)) then
      write(*,'(a)')'[ERROR] Could not compute AVG for temp(1)'
      return
    endif
    CoeffDynArray = num1D / den1D
    call computeLilly(M(1),M(2),M(3),L(1),L(2),L(3),CoeffDynArray,spec_rank,3)
    do i = 1,U%nz
      if ( .not. equalToVal (CoeffDynArray(i) , num1D(i) / den1D(i) )) then
        if (spec_rank .eq. 0 ) print *,'[ERROR] dans le calcul du coef dyn clark'
      endif
    enddo
    deallocate(num1D)
    deallocate(den1D)
  endif

  !Désallocation des tableaux de travail
  if( .not. deleteWorkArray(Uhat) .or. &
  & .not. deleteWorkArray(L) .or. &
  & .not. deleteWorkArray(H) .or. &
  & .not. deleteWorkArray(H_test) .or. &
  & .not. deleteWorkArray(ShatdZhat) .or. &
  & .not. deleteWorkArray(SdZhat) .or. &
  & .not. deleteWorkArray(SdZK) .or. &
  & .not. deleteWorkArray(SdZ) .or. &
  & .not. deleteWorkArray(SdZKhat) .or. &
  & .not. deleteWorkArray(UZ) .or. &
  & .not. deleteWorkArray(UZhat) .or. &
  & .not. deleteWorkArray(H1) .or. &
  & .not. deleteWorkArray(M) .or. &
  & .not. deleteWorkArray(M_test) .or. &
  & .not. deleteWorkArray(Scalhat) .or. &
  & .not. deleteWorkArray(H1K) .or. &
  & .not. deleteWorkArray(UZK) .or. &
  & .not. deleteWorkArray(UZKhat) .or. &
  & .not. deleteWorkArray(H2K) .or. &
  & .not. deleteWorkArray(UKhat) .or. &
  & .not. deleteWorkArray(ScalKhat) .or. &
  & .not. deleteWorkArray(temp) .or. &
  & .not. deleteWorkArray(numer) .or. &
  & .not. deleteWorkArray(denom) .or. &
  & .not. deleteWorkArray(tempK)) then
    write(6,'(a)') '[ERROR] memory leakage in coefDynClarkSca.'
    return
  endif

  success = .true.

 end function coefDynClarkScaTEST

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
 function coefDynSmagScaTEST(U,V,W,Uk,Vk,Wk,ScalIn,ScalKIn,ScalWN,filter,filterSize,spec_rank,CoeffDynScalar,CoeffDynArray) result(success)

  !I/O
  type(complex_data_layout),intent(in)                     :: Uk,Vk,Wk,ScalKIn 
  type(real_data_layout),intent(in)                        :: U,V,W,ScalIn
  type(WaveNumbers), intent(in)                            :: ScalWN
  integer, intent(in)                                      :: spec_rank,filter
  integer,intent(in)                                      :: filterSize 
  real(WP),intent(inout),optional                          :: coeffDynScalar 
  real(WP),intent(inout),dimension(:),allocatable,optional :: coeffDynArray 
  logical                                                  :: success

  !Local data
  real(WP)                                           :: delta
  real(WP)                                           :: sizeTestFilter
  type(real_data_layout),dimension(:),allocatable    :: Uhat
  type(real_data_layout),dimension(:),allocatable    :: L
  type(real_data_layout),dimension(:),allocatable    :: M
  type(real_data_layout),dimension(:),allocatable    :: M_test
  type(real_data_layout),dimension(:),allocatable    :: Scal
  type(real_data_layout),dimension(:),allocatable    :: Scalhat
  type(real_data_layout),dimension(:),allocatable    :: numer
  type(real_data_layout),dimension(:),allocatable    :: denom 
  type(real_data_layout),dimension(:),allocatable    :: SdZ
  type(real_data_layout),dimension(:),allocatable    :: ShatdZhat
  type(real_data_layout),dimension(:),allocatable    :: SdZhat
  type(complex_data_layout),dimension(:),allocatable :: SdZK
  type(complex_data_layout),dimension(:),allocatable :: UKhat
  type(complex_data_layout),dimension(:),allocatable :: ScalK
  type(complex_data_layout),dimension(:),allocatable :: ScalKhat
  type(complex_data_layout),dimension(:),allocatable :: SdZKhat
  real(WP),dimension(:),allocatable                  :: num1D
  real(WP),dimension(:),allocatable                  :: den1D
  integer                                            :: ires,i,j,k
  logical                                            :: res,erreur
  real(WP)  :: Lilly 
 
  success = .false.

  !Init des tableaux de travail
  if (.not. initWorkArray(U,3,Uhat) .or. &
     &.not. initWorkArray(Uk,3,UKhat) .or. &
     &.not. initWorkArray(U,3,L).or. &
     &.not. initWorkArray(U,3,M) .or. &
     &.not. initWorkArray(U,3,M_test) .or. &
     &.not. initWorkArray(U,1,Scal) .or. &
     &.not. initWorkArray(U,1,Scalhat) .or. &
     &.not. initWorkArray(Uk,1,ScalK) .or. &
     &.not. initWorkArray(Uk,1,ScalKhat) .or. &
     &.not. initWorkArray(U,3,SdZ) .or. &
     &.not. initWorkArray(Uk,3,SdZK) .or. &
     &.not. initWorkArray(U,3,ShatdZhat) .or. &
     &.not. initWorkArray(U,3,SdZhat) .or. &
     &.not. initWorkArray(Uk,3,SdZKhat) .or. &
     &.not. initWorkArray(U,1,numer) .or. &
     &.not. initWorkArray(U,1,denom)) then
    write(6,'(a)') '[ERROR] not enought memory in coefDynClarkSca'
    return
  endif

  sizeTestFilter  = 2.0_WP
  delta = real(filterSize*computeDelta(ScalIn),WP)
  !Filtering 
  call computeFilter(ScalWN,sizeTestFilter*delta,Uk,UKhat(1),filter) 
  call computeFilter(ScalWN,sizeTestFilter*delta,Vk,UKhat(2),filter) 
  call computeFilter(ScalWN,sizeTestFilter*delta,Wk,UKhat(3),filter) 
  call computeFilter(ScalWN,sizeTestFilter*delta,ScalKIn,ScalKhat(1),filter) 
  !compute Germano's term
  !calcul des Li
  call computeGermanoSca(ScalWN,U,V,W,ScalIn,Uk,Vk,Wk,ScalkIn,Filter,delta,L(1),L(2),L(3))
  !Compute Mi
  !Compute ShatdZhat_i
  call computeSdZ(U,UKhat(1),UKhat(2),UKhat(3),ScalKhat(1),ScalWN,ShatdZhat(1),ShatdZhat(2),ShatdZhat(3),res)
  !Compute (SdZ)hat_i
  call computeSdZ(U,Uk,Vk,Wk,ScalKIn,ScalWN,SdZ(1),SdZ(2),SdZ(3),res)
  call ftran(SdZ(1),SdZK(1),res)
  call ftran(SdZ(2),SdZK(2),res)
  call ftran(SdZ(3),SdZK(3),res)
  call computeFilter(ScalWN,sizeTestFilter*delta,SdZK(1),SdZKhat(1),filter) 
  call computeFilter(ScalWN,sizeTestFilter*delta,SdZK(2),SdZKhat(2),filter) 
  call computeFilter(ScalWN,sizeTestFilter*delta,SdZK(3),SdZKhat(3),filter) 
  call btran(SdZKhat(1),SdZhat(1),res)
  call btran(SdZKhat(2),SdZhat(2),res)
  call btran(SdZKhat(3),SdZhat(3),res)
  M(1)%values = (delta*sizeTestFilter)**2.0 * ShatdZhat(1)%values - delta**2.0*SdZhat(1)%values 
  M(2)%values = (delta*sizeTestFilter)**2.0 * ShatdZhat(2)%values - delta**2.0*SdZhat(2)%values 
  M(3)%values = (delta*sizeTestFilter)**2.0 * ShatdZhat(3)%values - delta**2.0*SdZhat(3)%values 

  call compute_Mi_Ni ( M(1),M(2),M(3),ScalWN,delta,2.0_WP,filter,Uk,.true.,res,&
                     & Vk=Vk,Wk=Wk,Scalk=ScalkIn,&
                     & Ukhat=UKhat(1),Vkhat=UKhat(2),Wkhat=UKhat(3),Scalkhat=ScalKhat(1) )
!!!!!!!!!!!!!!!!!!!!!!!!!!!TESTS!!!!!!!!!!!!!!!!!!!!!!!!  
  call compute_Mi_Ni ( M_test(1),M_test(2),M_test(3),ScalWN,delta,2.0_WP,filter,Uk,.true.,res,&
                     & Vk=Vk,Wk=Wk,Scalk=Scalkin,&
                     & Ukhat=UKhat(1),Vkhat=UKhat(2),Wkhat=UKhat(3),Scalkhat=ScalKhat(1) )
  erreur = .false.
  do k = U%zmin,U%zmax
  do j = U%ymin,U%ymax
  do i = U%xmin,U%xmax
    if ( .not. equalToVal (M_test(1)%values(i,j,k) ,M(1)%values(i,j,k)) .or. &
       & .not. equalToVal (M_test(2)%values(i,j,k) ,M(2)%values(i,j,k)) .or. &
       & .not. equalToVal (M_test(3)%values(i,j,k) ,M(3)%values(i,j,k)) ) then
      if (spec_rank .eq. 0 .and. .not. erreur ) print *,'[ERROR] problèmes dans M_test'
      erreur = .true.
    endif
  enddo
  enddo
  enddo

  call compute_Mi_Ni ( M_test(1),M_test(2),M_test(3),ScalWN,delta,2.0_WP,filter,Uk,.true.,res,&
                     & Vk=Vk,Wk=Wk,Scalk=Scalkin)
                     
  erreur = .false.
  do k = U%zmin,U%zmax
  do j = U%ymin,U%ymax
  do i = U%xmin,U%xmax
    if ( .not. equalToVal (M_test(1)%values(i,j,k) ,M(1)%values(i,j,k)) .or. &
       & .not. equalToVal (M_test(2)%values(i,j,k) ,M(2)%values(i,j,k)) .or. &
       & .not. equalToVal (M_test(3)%values(i,j,k) ,M(3)%values(i,j,k))) then
      if (spec_rank .eq. 0 .and. .not. erreur ) print *,'[ERROR] problèmes dans M_test'
      erreur = .true.
    endif
  enddo
  enddo
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!TESTS!!!!!!!!!!!!!!!!!!!!!!!!  

  !Somme numérateur et dénominateur
  numer(1)%values = L(1)%values * M(1)%values + L(2)%values * M(2)%values + L(3)%values * M(3)%values
  denom(1)%values = M(1)%values * M(1)%values + M(2)%values * M(2)%values + M(3)%values * M(3)%values
  !Moyenne isotrope
  if (present(CoeffDynScalar)) then
    CoeffDynScalar = computeFieldAvg(numer(1),spec_rank) / computeFieldAvg(denom(1),spec_rank)
    if (spec_rank .eq. 0) write (6,'(a,1x,f9.6)') '[INFO] CoefDyn SmagModel TEST=',CoeffDynScalar
    call computeLilly(M(1),M(2),M(3),L(1),L(2),L(3),Lilly,spec_rank)
    if ( .not. equalToVal (Lilly ,CoeffDynScalar )) then
      if (spec_rank .eq. 0) write (6,'(a)') '[ERROR] in coefDyn SmagModel'
    endif
  endif
  !Moyenne anisotrope
  if (present(CoeffDynArray)) then
    if (.not.( (lbound(CoeffDynArray,1) .eq. 1 ) .and. (ubound(CoeffDynArray,1) .eq. U%nz) ) )then
       write (6,'(a)') '[ERROR] CoeffDynArray not declared between good bounds'
       return
    endif
    allocate(num1D(U%nz),stat=ires)
    if ( ires .ne. 0 ) then
      write (6,'(a)') '[ERROR] Not able to allocate num1D in coefDynClarkSca'
      return
    endif
    allocate(den1D(U%nz),stat=ires)
    if ( ires .ne. 0 ) then
      write (6,'(a)') '[ERROR] Not able to allocate den1D in coefDynClarkSca'
      return
    endif
    if (.not. computeFieldAvg2D(numer(1),num1D,spec_rank)) then
      write(*,'(a)')'[ERROR] Could not compute AVG for temp(4)'
      return
    endif
    if (.not. computeFieldAvg2D(denom(1),den1D,spec_rank)) then
      write(*,'(a)')'[ERROR] Could not compute AVG for temp(1)'
      return
    endif
    do i = 1,U%nz
      CoeffDynArray(i) = num1D(i) / den1D(i)
    enddo
    call computeLilly(M(1),M(2),M(3),L(1),L(2),L(3),CoeffDynArray,spec_rank,3)
    do i = 1,U%nz
      if ( .not. equalToVal (CoeffDynArray(i) , num1D(i) / den1D(i) )) then
        if (spec_rank .eq. 0 ) print *,'[ERROR] dans le calcul du coef dyn clark'
      endif
    enddo
    deallocate(num1D)
    deallocate(den1D)
  endif
  !Free memory
  if ( .not. deleteWorkArray(Uhat) .or. &
  & .not. deleteWorkArray(UKhat) .or. &
  & .not. deleteWorkArray(L) .or. &
  & .not. deleteWorkArray(M) .or. &
  & .not. deleteWorkArray(M_test) .or. &
  & .not. deleteWorkArray(Scal) .or. &
  & .not. deleteWorkArray(Scalhat) .or. &
  & .not. deleteWorkArray(ScalK) .or. &
  & .not. deleteWorkArray(ScalKhat) .or. &
  & .not. deleteWorkArray(SdZ) .or. &
  & .not. deleteWorkArray(SdZK) .or. &
  & .not. deleteWorkArray(ShatdZhat) .or. &
  & .not. deleteWorkArray(SdZhat) .or. &
  & .not. deleteWorkArray(SdZKhat) .or. &
  & .not. deleteWorkArray(numer) .or. &
  & .not. deleteWorkArray(denom) ) then
    write (6,'(a)') '[ERROR] memory leackage in coefDynClarkSca'
    return
  endif

  success = .true.


 end function coefDynSmagScaTEST


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
!> @param [in] filterSize is the value of the filter for a priori test
!> @param [in] spec_rank is the number of process in the cpu pool
!> @param [inout] CoeffDynScalar is dynamic coefficient which is computed by the routine. Average is done all over the domain 
!> @param [inout] CoeffDynArray is dynamic coefficient which is computed by the routine. Average is done in 2D slice 
!> @return      success is a logical equal to .true. if everything is Ok and .false. otherewize.
!------------------------------------------------------------------------------
 function coefDynClarkFabreScaTEST(U,V,W,Uk,Vk,Wk,ScalIn,ScalKIn,ScalWN,filter,filterSize,spec_rank,&
                              &CoeffDynScalar,CoeffDynArray) result(success)

  !I/O
  type(complex_data_layout),intent(in)          :: Uk,Vk,Wk,ScalKIn 
  type(real_data_layout),intent(in)             :: U,V,W,ScalIn
  type(WaveNumbers), intent(in)                 :: ScalWN
  integer, intent(in)                           :: spec_rank,filter
  integer, intent(in)                           :: filterSize
  real(WP), intent(inout),dimension(:),optional :: CoeffDynArray
  real(WP), intent(inout),optional              :: CoeffDynScalar
  logical                                       :: success

  !Local data
  real(WP)                                           :: sizeTestFilter
  type(real_data_layout),dimension(:),allocatable    :: Uhat
  type(real_data_layout),dimension(:),allocatable    :: L
  type(real_data_layout),dimension(:),allocatable    :: L_test
  type(real_data_layout),dimension(:),allocatable    :: UZhat
  type(real_data_layout),dimension(:),allocatable    :: UZ
  type(real_data_layout),dimension(:),allocatable    :: K 
  type(real_data_layout),dimension(:),allocatable    :: N
  type(real_data_layout),dimension(:),allocatable    :: K_test 
  type(real_data_layout),dimension(:),allocatable    :: N_test 
  type(real_data_layout),dimension(:),allocatable    :: Scalhat
  type(real_data_layout),dimension(:),allocatable    :: temp
  type(real_data_layout),dimension(:),allocatable    :: numer
  type(real_data_layout),dimension(:),allocatable    :: denom
  type(complex_data_layout),dimension(:),allocatable :: tempK
  type(complex_data_layout),dimension(:),allocatable :: UZK
  type(complex_data_layout),dimension(:),allocatable :: UZKhat
  type(complex_data_layout),dimension(:),allocatable :: UKhat
  type(complex_data_layout),dimension(:),allocatable :: ScalKhat
  real(WP)                                           :: delta,Lilly
  real(WP),dimension(:),allocatable                  :: num1D
  real(WP),dimension(:),allocatable                  :: den1D
  integer                                            :: ires,i,j,o
  logical                                            :: res,erreur 

  success = .false.

  !Init des tableaux de travail
  if( .not. initWorkArray(U,3,Uhat) .or. &
  & .not. initWorkArray(U,3,L) .or. &
  & .not. initWorkArray(U,3,L_test) .or. &
  & .not. initWorkArray(U,12,temp) .or. &
  & .not. initWorkArray(Uk,12,tempK) .or. &
  & .not. initWorkArray(U,3,K) .or. &
  & .not. initWorkArray(U,3,N) .or. &
  & .not. initWorkArray(U,3,K_test) .or. &
  & .not. initWorkArray(U,3,N_test) .or. &
  & .not. initWorkArray(U,3,UZ) .or. &
  & .not. initWorkArray(U,3,UZhat) .or. &
  & .not. initWorkArray(U,1,Scalhat) .or. &
  & .not. initWorkArray(U,1,numer) .or. &
  & .not. initWorkArray(U,1,denom) .or. &
  & .not. initWorkArray(Uk,3,UZK) .or. &
  & .not. initWorkArray(Uk,3,UZKhat) .or. &
  & .not. initWorkArray(Uk,3,UKhat) .or. &
  & .not. initWorkArray(Uk,1,ScalKhat)) then
    write(6,'(a)') '[ERROR] not enought memory in coefDynClarkSca'
    return
  endif
  !valeur de la taille de filtre
  delta = real(filterSize*computeDelta(ScalIn),WP)
  !coef multiplicateur entre le filtre et le filtre test
  sizeTestFilter  = 2.0_WP
  !filtre test pour le calcul des coefs
  !Filtrage de tous les termes par le filtre 'filtre' à la longueur de coupure nd*delta
  call computeFilter(ScalWN,sizeTestFilter*delta,Uk,UKhat(1),filter)
  call computeFilter(ScalWN,sizeTestFilter*delta,Vk,UKhat(2),filter)
  call computeFilter(ScalWN,sizeTestFilter*delta,Wk,UKhat(3),filter)
  call computeFilter(ScalWN,sizeTestFilter*delta,ScalKIn,ScalKhat(1),filter)
  !Retour espace physique des quantitées filtrées
  call btran(UKhat(1),Uhat(1),res) 
  call btran(UKhat(2),Uhat(2),res) 
  call btran(UKhat(3),Uhat(3),res)  
  call btran(ScalKhat(1),Scalhat(1),res)  !Zhatbar
  !calcul des Li
  UZ(1)%values = U%values * ScalIn%values
  UZ(2)%values = V%values * ScalIn%values
  UZ(3)%values = W%values * ScalIn%values
  call ftran(UZ(1),UZK(1),res) 
  call ftran(UZ(2),UZK(2),res) 
  call ftran(UZ(3),UZK(3),res)  
  call computeFilter(ScalWN,sizeTestFilter*delta,UZK(1),UZKhat(1),filter)
  call computeFilter(ScalWN,sizeTestFilter*delta,UZK(2),UZKhat(2),filter)
  call computeFilter(ScalWN,sizeTestFilter*delta,UZK(3),UZKhat(3),filter)
  call btran(UZKhat(1),UZhat(1),res) 
  call btran(UZKhat(2),UZhat(2),res) 
  call btran(UZKhat(3),UZhat(3),res)  
  L(1)%values=UZhat(1)%values - Scalhat(1)%values * Uhat(1)%values 
  L(2)%values=UZhat(2)%values - Scalhat(1)%values * Uhat(2)%values
  L(3)%values=UZhat(3)%values - Scalhat(1)%values * Uhat(3)%values

  call computeGermanoSca(ScalWN,U,V,W,ScalIn,Uk,Vk,Wk,ScalkIn,Filter,delta,L_test(1),L_test(2),L_test(3))

  erreur = .false.
  do o = U%zmin,U%zmax
  do j = U%ymin,U%ymax
  do i = U%xmin,U%xmax
    if ( .not. equalToVal (L_test(1)%values(i,j,o) ,L(1)%values(i,j,o)) .or. &
       & .not. equalToVal (L_test(2)%values(i,j,o) ,L(2)%values(i,j,o)) .or. &
       & .not. equalToVal (L_test(3)%values(i,j,o) ,L(3)%values(i,j,o)) ) then
      if (spec_rank .eq. 0 .and. .not. erreur ) print *,'[ERROR] problèmes dans M_test'
      erreur = .true.
    endif
  enddo
  enddo
  enddo

  !calcul des Ki
  call computeGradientK(ScalWN,UKhat(1),tempK(1),tempK(2),tempK(3))
  call computeGradientK(ScalWN,UKhat(2),tempK(4),tempK(5),tempK(6))
  call computeGradientK(ScalWN,UKhat(3),tempK(7),tempK(8),tempK(9))
  call computeGradientK(ScalWN,ScalKhat(1),tempK(10),tempK(11),tempK(12))
  call btran(tempK(1),temp(1),res) 
  call btran(tempK(2),temp(2),res) 
  call btran(tempK(3),temp(3),res) 
  call btran(tempK(4),temp(4),res) 
  call btran(tempK(5),temp(5),res) 
  call btran(tempK(6),temp(6),res) 
  call btran(tempK(7),temp(7),res) 
  call btran(tempK(8),temp(8),res) 
  call btran(tempK(9),temp(9),res) 
  call btran(tempK(10),temp(10),res) 
  call btran(tempK(11),temp(11),res) 
  call btran(tempK(12),temp(12),res) 
  K(1)%values = (sizeTestFilter*delta)**2/12.0_WP * (temp(1)%values * temp(10)%values + &
                                                    &temp(2)%values * temp(11)%values + &
                                                    &temp(3)%values * temp(12)%values)
  K(2)%values = (sizeTestFilter*delta)**2/12.0_WP * (temp(4)%values * temp(10)%values + &
                                                    &temp(5)%values * temp(11)%values + &
                                                    &temp(6)%values * temp(12)%values)
  K(3)%values = (sizeTestFilter*delta)**2/12.0_WP * (temp(7)%values * temp(10)%values + &
                                                    &temp(8)%values * temp(11)%values + &
                                                    &temp(9)%values * temp(12)%values)
!!!!!!!!!!!!!!!!!!!!!!!!!!!TESTS!!!!!!!!!!!!!!!!!!!!!!!!  
  call compute_Hi_Ki( K_test(1),K_test(2),K_test(3),ScalWN,delta,sizeTestFilter,filter,Uk,.false.,res,&
                         & Vk=Vk,Wk=Wk,Scalk=Scalkin,&
                         & Ukhat=Ukhat(1),Vkhat=Ukhat(2),Wkhat=Ukhat(3),Scalkhat=ScalKhat(1) )
  erreur = .false.
  do o = U%zmin,U%zmax
  do j = U%ymin,U%ymax
  do i = U%xmin,U%xmax
    if ( .not. equalToVal (K_test(1)%values(i,j,o) ,K(1)%values(i,j,o)) .or. &
       & .not. equalToVal (K_test(2)%values(i,j,o) ,K(2)%values(i,j,o)) .or. &
       & .not. equalToVal (K_test(3)%values(i,j,o) ,K(3)%values(i,j,o)) ) then
      if (spec_rank .eq. 0 .and. .not. erreur ) print *,'problèmes dans M_test'
      erreur = .true.
    endif
  enddo
  enddo
  enddo

  call compute_Hi_Ki ( K_test(1),K_test(2),K_test(3),ScalWN,delta,2.0_WP,filter,Uk,.false.,res,&
                     & Vk=Vk,Wk=Wk,Scalk=Scalkin)
                     
  erreur = .false.
  do o = U%zmin,U%zmax
  do j = U%ymin,U%ymax
  do i = U%xmin,U%xmax
    if ( .not. equalToVal (K_test(1)%values(i,j,o) ,K(1)%values(i,j,o)) .or. &
       & .not. equalToVal (K_test(2)%values(i,j,o) ,K(2)%values(i,j,o)) .or. &
       & .not. equalToVal (K_test(3)%values(i,j,o) ,K(3)%values(i,j,o))) then
      if (spec_rank .eq. 0 .and. .not. erreur ) print *,'problèmes dans H_test'
      erreur = .true.
    endif
  enddo
  enddo
  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!TESTS!!!!!!!!!!!!!!!!!!!!!!!!  
  !calcul des Ni   
  call computeSdZ(Uhat(1),UKhat(1),UKhat(2),UKhat(3),ScalKhat(1),&
                 &ScalWN,N(1),N(2),N(3),res)
  N(1)%values = (sizeTestFilter*delta)**2*N(1)%values
  N(2)%values = (sizeTestFilter*delta)**2*N(2)%values
  N(3)%values = (sizeTestFilter*delta)**2*N(3)%values
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!TESTS!!!!!!!!!!!!!!!!!!!!!!!!  
  call compute_Mi_Ni ( N_test(1),N_test(2),N_test(3),ScalWN,delta,2.0_WP,filter,Uk,.false.,res,&
                     & Vk=Vk,Wk=Wk,Scalk=Scalkin,&
                     & Ukhat=UKhat(1),Vkhat=UKhat(2),Wkhat=UKhat(3),Scalkhat=ScalKhat(1) )
  erreur = .false.
  do o = U%zmin,U%zmax
  do j = U%ymin,U%ymax
  do i = U%xmin,U%xmax
    if ( .not. equalToVal (N_test(1)%values(i,j,o) ,N(1)%values(i,j,o)) .or. &
       & .not. equalToVal (N_test(2)%values(i,j,o) ,N(2)%values(i,j,o)) .or. &
       & .not. equalToVal (N_test(3)%values(i,j,o) ,N(3)%values(i,j,o)) ) then
      if (spec_rank .eq. 0 .and. .not. erreur ) print *,'[ERROR] problemes dans N_test'
      erreur = .true.
    endif
  enddo
  enddo
  enddo

  call compute_Mi_Ni ( N_test(1),N_test(2),N_test(3),ScalWN,delta,2.0_WP,filter,Uk,.true.,res,&
                     & Vk=Vk,Wk=Wk,Scalk=Scalkin)
                     
  erreur = .false.
  do o = U%zmin,U%zmax
  do j = U%ymin,U%ymax
  do i = U%xmin,U%xmax
    if ( .not. equalToVal (N_test(1)%values(i,j,o) ,N(1)%values(i,j,o)) .or. &
       & .not. equalToVal (N_test(2)%values(i,j,o) ,N(2)%values(i,j,o)) .or. &
       & .not. equalToVal (N_test(3)%values(i,j,o) ,N(3)%values(i,j,o))) then
      if (spec_rank .eq. 0 .and. .not. erreur ) print *,'[ERROR] problemes dans N_test'
      erreur = .true.
    endif
  enddo
  enddo
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!TESTS!!!!!!!!!!!!!!!!!!!!!!!!  

  L(1)%values = L(1)%values-K(1)%values
  L(2)%values = L(2)%values-K(2)%values
  L(3)%values = L(3)%values-K(3)%values

  !Somme numérateur et dénominateur
  numer(1)%values = L(1)%values * N(1)%values + L(2)%values * N(2)%values + L(3)%values * N(3)%values
  denom(1)%values = N(1)%values * N(1)%values + N(2)%values * N(2)%values + N(3)%values * N(3)%values
  
  !Moyenne isotrope
  if (present(CoeffDynScalar)) then
    call computeLilly(N(1),N(2),N(3),L(1),L(2),L(3),CoeffDynScalar,spec_rank)
    if (spec_rank .eq. 0) write (6,'(a,1x,f9.6)') '[INFO] CoefDyn ClarkFabreModel TEST=',CoeffDynScalar
    call computeLilly(N(1),N(2),N(3),L(1),L(2),L(3),Lilly,spec_rank)
    if ( .not. equalToVal (Lilly ,CoeffDynScalar )) then
      if (spec_rank .eq. 0) write (6,'(a)') '[ERROR] in coefDyn ClarkFabreModel'
    endif
  endif
  !Moyenne anisotrope
  if (present(CoeffDynArray)) then
    allocate(num1D(U%nz),stat=ires)
    if ( ires .ne. 0 ) then
      write (6,'(a)') '[ERROR] Not able to allocate num1D in coefDynClarkFabreSca'
      return
    endif
    allocate(den1D(U%nz),stat=ires)
    if ( ires .ne. 0 ) then
      write (6,'(a)') '[ERROR] Not able to allocate den1D in coefDynClarkFabreSca'
      return
    endif
    if (.not. computeFieldAvg2D(numer(1),num1D,spec_rank)) then
      write(*,'(a)')'[ERROR] Could not compute AVG for temp(4)'
      return
    endif
    if (.not. computeFieldAvg2D(denom(1),den1D,spec_rank)) then
      write(*,'(a)')'[ERROR] Could not compute AVG for temp(1)'
      return
    endif
    CoeffDynArray = num1D / den1D
    call computeLilly(N(1),N(2),N(3),L(1),L(2),L(3),CoeffDynArray,spec_rank,3)
    deallocate(num1D)
    deallocate(den1D)
  endif
  !Désallocation des tableaux de travail
  if ( .not. deleteWorkArray(Uhat) .or. &
  & .not. deleteWorkArray(L) .or. &
  & .not. deleteWorkArray(K) .or. &
  & .not. deleteWorkArray(N) .or. &
  & .not. deleteWorkArray(L_test) .or. &
  & .not. deleteWorkArray(K_test) .or. &
  & .not. deleteWorkArray(N_test) .or. &
  & .not. deleteWorkArray(UZ) .or. &
  & .not. deleteWorkArray(UZhat) .or. &
  & .not. deleteWorkArray(Scalhat) .or. &
  & .not. deleteWorkArray(UZK) .or. &
  & .not. deleteWorkArray(UZKhat) .or. &
  & .not. deleteWorkArray(UKhat) .or. &
  & .not. deleteWorkArray(ScalKhat) .or. &
  & .not. deleteWorkArray(temp) .or. &
  & .not. deleteWorkArray(numer) .or. &
  & .not. deleteWorkArray(denom) .or. &
  & .not. deleteWorkArray(tempK)) then
    write(6,'(a)') '[ERROR] memory leakage in coefDynClarkSca'
    return
  endif

  success = .true.

 end function coefDynClarkFabreScaTEST



!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> A partir d'une matrice symétrique réelle, calcul de sa partie positive ou négative tq:
!> Sij^+ = \sum_k=1^3 max(lam(i)) VP_i(k) VP_j(k)^t
!> Sij^- = \sum_k=1^3 min(lam(i)) VP_i(k) VP_j(k)^t
!> avec lam(i) les valeurs propre de Sij et VP(i) les vecteurs propres unitaires
!> @param [in] S11in composante (1,1) de la matrice symétrique à décomposer
!> @param [in] S12in composante (1,2) de la matrice symétrique à décomposer
!> @param [in] S13in composante (1,3) de la matrice symétrique à décomposer
!> @param [in] S22in composante (2,2) de la matrice symétrique à décomposer
!> @param [in] S23in composante (2,3) de la matrice symétrique à décomposer
!> @param [in] S33in composante (3,3) de la matrice symétrique à décomposer
!> @param [inout] S11out composante (1,1) de la matrice symétrique décomposée
!> @param [inout] S12out composante (1,2) de la matrice symétrique décomposée
!> @param [inout] S13out composante (1,3) de la matrice symétrique décomposée
!> @param [inout] S22out composante (2,2) de la matrice symétrique décomposée
!> @param [inout] S23out composante (2,3) de la matrice symétrique décomposée
!> @param [inout] S33out composante (3,3) de la matrice symétrique décomposée
!> @param [in] signe peut prendre les valeurs "moins" ou "plus" pour le calcul de Sij^+ ou Sij^- 
!> @param [inout] lam1 première vp (optionnel)
!> @param [inout] lam2 deuxième vp (optionnel)
!> @param [inout] lam3 troisième vp (optionnel)
!> @param [inout] VPout contient les trois vecteurs propres
!>
!> @return      success is a logical equal to .true. if everything is Ok and .false. otherewize.
!------------------------------------------------------------------------------

 function computeSplitingSij_analyt_sca(S11in,S12in,S13in,S22in,S23in,S33in,&
                          &S11out,S12out,S13out,S22out,S23out,S33out,&
                          &signe,lam1,lam2,lam3,VPout) result(success)

  use stat_tools , only : equalToVal
  use mpilayout_tools , only : getnbmycpu

  !I/O data
  real(WP),intent(in)                             :: S11in,S12in,S13in,S22in,S23in,S33in
  real(WP),intent(inout)                          :: S11out,S12out,S13out,S22out,S23out,S33out
  real(WP),intent(inout),optional                 :: lam1,lam2,lam3 
  real(WP),dimension(3,3),intent(inout),optional  :: VPout
  character(len=*),intent(in)                     :: signe
  logical                                         :: success

  !Local data
  integer                 :: i,j,k,spec_rank
  real(WP)                :: u,v,w,q,p
  real(WP)                :: temp
  real(WP),dimension(3)   :: lam
  real(WP),dimension(3,3) :: VP

  !Init
  success = .false.
  VP = 0.0_WP
  lam = 0.0_WP 
  spec_rank = getnbmycpu()

  ! Calcul des coefficients du polynome caracteristique
  ! det (lam x Id - Sij) = 0
  ! lam^3 - tr(Sij) lam^2 - 0.5 (tr(Sij^2)-tr(Sij)^2)lam - det(Sij)
  ! donc:
  ! a= 1
  ! b= -tr(Sij) = div(u) = 0
  ! b= - (S11 + S22 + S33)
  ! c= -0.5 (tr(Sij^2)-tr(Sij)^2)
  ! c= S11*S22 + S22*S33 + S11*S33 - S12^2 - S23^2 - S13^2
  ! d= -det(Sij)
  ! d= S13^2*S22 + S12^2*S33 + S23^2*S11 - S11*S22*S33 - 2*S12*S23*S13
  p = S11in*S33in + S22in*S33in + S11in*S22in - S13in**2 - S12in**2 - S23in**2
  q = S13in**2*S22in + S12in**2*S33in + S23in**2*S11in - S11in*S22in*S33in - 2.0_WP*S12in*S23in*S13in
  !Calcul des racines du polynome
  !Formule de Cardan avec p = c et q = d
  if ( .not. equalToVal (0.0_WP , S11in + S22in + S33in)) then
    write (6,'(a,1x,g15.8)') '[ERROR] in computeSplitingSij: divU != 0 :',S11in + S22in + S33in
    return
  endif
  !Calcul des racines du polynôme caractéristique
  do i = 1,3
    u = -q/2.0_WP*sqrt(27.0_WP/ ( -p**3 ))
    if ( u .lt. -1.0_WP .or. u .gt. 1.0_WP ) then
      write (6,'(a,1x,g15.8)') '[ERROR] u out of bounds',u
      return
    endif
    u = 1.0_WP/3.0_WP * acos( u )
    v = 2.0_WP * real(i,WP) * acos(-1.0_WP) / 3.0_WP 
    w = 2.0_WP*sqrt(-p/3.0_WP)
    lam(i) = w * cos ( u + v ) 
  enddo
  !Calcul des vecteurs propres
  i=1
  do while ( (S13in*(S22in-lam(i))-S22in*S12in .lt. 1e-12) .and. &
           & (S23in*(S11in-lam(i))-S12in*S13in .lt. 1e-12) .and. &
           & (S12in**2-(S22in-lam(i))*(S11in-lam(i)) .lt. 1e-12) .and. (i .le. 3) )
    i=i+1
  end do
  VP(1,i) = S13in*S23in - (S33in-lam(i))*S12in
  VP(1,i) = VP(1,i)/(S12in*S13in-(S11in-lam(i))*S23in)
  VP(2,i) = (S33in-lam(i))*S12in - S13in*S23in 
  VP(2,i) = VP(2,i)/(S13in*(S22in-lam(i))-S12in*S23in)
  VP(3,i) = 1.0_WP/sqrt(1+VP(1,i)**2+VP(2,i)**2)
  VP(1,i) = VP(1,i) * VP(3,i)
  VP(2,i) = VP(2,i) * VP(3,i)
  j = i + 1
  if ( i .eq. 3) j=1
  VP(1,j) = S13in*S23in - (S33in-lam(j))*S12in
  VP(1,j) = VP(1,j)/(S12in*S13in-(S11in-lam(j))*S23in)
  VP(2,j) = (S33in-lam(j))*S12in - S13in*S23in 
  VP(2,j) = VP(2,j)/(S13in*(S22in-lam(j))-S12in*S23in)
  VP(3,j) = 1.0_WP/sqrt(1+VP(1,j)**2+VP(2,j)**2)
  VP(1,j) = VP(1,j) * VP(3,j)
  VP(2,j) = VP(2,j) * VP(3,j)

  k = j + 1
  if ( j .eq. 3) k=1
  VP(1,k) = (VP(2,i)*VP(3,j)-VP(3,i)*VP(2,j))
  VP(2,k) = (VP(3,i)*VP(1,j)-VP(1,i)*VP(3,j))
  VP(3,k) = (VP(1,i)*VP(2,j)-VP(2,i)*VP(1,j))
  !selection des vp >0 ou <0
  ! D^(+||-)
  do i=1,3
    if (trim(adjustl(signe)) .eq. 'plus') then
      temp = max(lam(i),0.0_WP)
      lam(i) = temp
    elseif( trim(adjustl(signe)) .eq. 'moins') then
      temp = min(lam(i),0.0_WP)
      lam(i) = temp
    else
      write (6,'(a)') '[ERROR] in computeSplitingSij : "signe" unknown!' 
      return
    endif
  enddo
  !Calcul de Sij^+ ou Sij^- 
  ! S = VP D^(+||-) VP^t = Sij = Sum_k lam(k) e(i) e(j)
  S11out=lam(1)*VP(1,1)**2 + lam(2)*VP(1,2)**2 + lam(3)*VP(1,3)**2
  S12out=lam(1)*VP(1,1)*VP(2,1) + lam(2)*VP(1,2)*VP(2,2) + lam(3)*VP(1,3)*VP(2,3)
  S13out=lam(1)*VP(1,1)*VP(3,1) + lam(2)*VP(1,2)*VP(3,2) + lam(3)*VP(1,3)*VP(3,3)
  S22out=lam(1)*VP(2,1)**2 + lam(2)*VP(2,2)**2 + lam(3)*VP(2,3)**2
  S23out=lam(1)*VP(2,1)*VP(3,1) + lam(2)*VP(2,2)*VP(3,2) + lam(3)*VP(2,3)*VP(3,3)
  S33out=lam(1)*VP(3,1)**2 + lam(2)*VP(3,2)**2 + lam(3)*VP(3,3)**2
  !Sortie optionnelle
  if ( present(lam1) .and. present(lam2) .and. present(lam3)) then
    lam1 = lam(1)
    lam2 = lam(2)
    lam3 = lam(3)
  endif
  if ( present(VPout)) then
    VPout = VP
  endif

  success = .true.

 end function computeSplitingSij_analyt_sca

!------------------------------------------------------------------------------
!> @author 
!>
!> @details
!> 
!> @param [in] Uin,Vin,Win,Zin
!> @param [in] sizeFilt
!> @param [in] filter
!> @param [inout] res
!> @param [in,out] Ufilout,Vfilout,Wfilout,Zfilout
!> @param [in,out] Ukfilout,Vkfilout,Wkfilout,Zkfilout
!> @param [out] success .TRUE. if all is succcessfull, .FALSE. otherwize.
!------------------------------------------------------------------------------

 function computeSplitingSij_analyt_field(S11in,S12in,S13in,S22in,S23in,S33in,&
                          &S11out,S12out,S13out,S22out,S23out,S33out,&
                          &signe) result(success)

  use mpilayout_tools , only : getnbmycpu

  !I/O data
  type(real_data_layout),intent(in)               :: S11in,S12in,S13in,S22in,S23in,S33in
  type(real_data_layout),intent(inout)            :: S11out,S12out,S13out,S22out,S23out,S33out
  character(len=*),intent(in)                     :: signe
  logical                                         :: success

  !Local data
  integer                  :: i,j,k

  success =.false.
  do k = S11in%zmin, S11in%zmax
    do j = S11in%ymin, S11in%ymax
      do i = S11in%xmin, S11in%xmax
        if (.not. computeSplitingSij_analyt_sca(S11in%values(i,j,k),S12in%values(i,j,k),S13in%values(i,j,k),&
                                             &S22in%values(i,j,k),S23in%values(i,j,k),S33in%values(i,j,k),&
                                             &S11out%values(i,j,k),S12out%values(i,j,k),S13out%values(i,j,k),&
                                             &S22out%values(i,j,k),S23out%values(i,j,k),S33out%values(i,j,k),signe)) then
          write (6,'(3(a,1x,i0,1x),a,1x,i0)') '[ERROR] in computeSplitingSij_analyt_field at this point i=',i,'j=',j,'k=',k,'on processus',getnbmycpu() 
          return
        endif
      enddo
    enddo
  enddo
  success =.true.

 end function computeSplitingSij_analyt_field

!------------------------------------------------------------------------------
!> @author 
!>
!> @details
!> 
!> @param [in] Uin,Vin,Win,Zin
!> @param [in] sizeFilt
!> @param [in] filter
!> @param [inout] res
!> @param [in,out] Ufilout,Vfilout,Wfilout,Zfilout
!> @param [in,out] Ukfilout,Vkfilout,Wkfilout,Zkfilout
!> @param [out] success .TRUE. if all is succcessfull, .FALSE. otherwize.
!------------------------------------------------------------------------------
!Test LAPACK
 function computeSplitingSij_lapack_field(S11in,S12in,S13in,S22in,S23in,S33in,&
                          &S11out,S12out,S13out,S22out,S23out,S33out,&
                          &signe) result(success)

  use stat_tools , only : equalToVal
  use mpilayout_tools , only : getnbmycpu

  !I/O data
  type(real_data_layout),intent(in)               :: S11in,S12in,S13in,S22in,S23in,S33in
  type(real_data_layout),intent(inout)            :: S11out,S12out,S13out,S22out,S23out,S33out
  character(len=*),intent(in)                     :: signe
  logical                                         :: success

  !Local data
  integer                  :: i,j,k
  real(WP)                 :: temp
  integer                  :: a,INFO
  real(WP), dimension(3,3) :: tab,VecP,VL
  real(WP), dimension(3)   :: valR,valI
  real(WP)                 :: tabwork(102)

  success =.false.
#ifdef LAPACK
    if ( trim(adjustl(signe)) == 'moins' ) then
      do k = S11in%zmin, S11in%zmax
         do j = S11in%ymin, S11in%ymax
            do i = S11in%xmin, S11in%xmax
               tab(1,1) = S11in%values(i,j,k)
               tab(1,2) = S12in%values(i,j,k)
               tab(1,3) = S13in%values(i,j,k)
               tab(2,1) = S12in%values(i,j,k)
               tab(2,2) = S22in%values(i,j,k)
               tab(2,3) = S23in%values(i,j,k)
               tab(3,1) = S13in%values(i,j,k)
               tab(3,2) = S23in%values(i,j,k)
               tab(3,3) = S33in%values(i,j,k)
               call DGEEV('N','V',3,tab,3,valR,valI,VL,3,VecP,3,tabwork,102,INFO)
               do a=1,3
                 temp=min(0.0_WP,valR(a))
                 valR(a)=temp
               enddo
               S11out%values(i,j,k)=ValR(1)*VecP(1,1)**2   + ValR(2)*VecP(1,2)**2   + ValR(3)*VecP(1,3)**2
               S12out%values(i,j,k)=ValR(1)*VecP(1,1)*VecP(2,1) + ValR(2)*VecP(1,2)*VecP(2,2) + ValR(3)*VecP(1,3)*VecP(2,3)
               S13out%values(i,j,k)=ValR(1)*VecP(1,1)*VecP(3,1) + ValR(2)*VecP(1,2)*VecP(3,2) + ValR(3)*VecP(1,3)*VecP(3,3)
               S22out%values(i,j,k)=ValR(1)*VecP(2,1)**2   + ValR(2)*VecP(2,2)**2   + ValR(3)*VecP(2,3)**2
               S23out%values(i,j,k)=ValR(1)*VecP(2,1)*VecP(3,1) + ValR(2)*VecP(2,2)*VecP(3,2) + ValR(3)*VecP(2,3)*VecP(3,3)
               S33out%values(i,j,k)=ValR(1)*VecP(3,1)**2   + ValR(2)*VecP(3,2)**2   + ValR(3)*VecP(3,3)**2
            end do
         end do
      end do
    elseif (trim(adjustl(signe)) == 'plus' ) then
      do k = S11in%zmin, S11in%zmax
         do j = S11in%ymin, S11in%ymax
            do i = S11in%xmin, S11in%xmax
               tab(1,1) = S11in%values(i,j,k)
               tab(1,2) = S12in%values(i,j,k)
               tab(1,3) = S13in%values(i,j,k)
               tab(2,1) = S12in%values(i,j,k)
               tab(2,2) = S22in%values(i,j,k)
               tab(2,3) = S23in%values(i,j,k)
               tab(3,1) = S13in%values(i,j,k)
               tab(3,2) = S23in%values(i,j,k)
               tab(3,3) = S33in%values(i,j,k)
               call DGEEV('N','V',3,tab,3,valR,valI,VL,3,VecP,3,tabwork,102,INFO)
               do a=1,3
                 temp=max(0.0_WP,valR(a))
                 valR(a)=temp
               enddo
               S11out%values(i,j,k)=ValR(1)*VecP(1,1)**2   + ValR(2)*VecP(1,2)**2   + ValR(3)*VecP(1,3)**2
               S12out%values(i,j,k)=ValR(1)*VecP(1,1)*VecP(2,1) + ValR(2)*VecP(1,2)*VecP(2,2) + ValR(3)*VecP(1,3)*VecP(2,3)
               S13out%values(i,j,k)=ValR(1)*VecP(1,1)*VecP(3,1) + ValR(2)*VecP(1,2)*VecP(3,2) + ValR(3)*VecP(1,3)*VecP(3,3)
               S22out%values(i,j,k)=ValR(1)*VecP(2,1)**2   + ValR(2)*VecP(2,2)**2   + ValR(3)*VecP(2,3)**2
               S23out%values(i,j,k)=ValR(1)*VecP(2,1)*VecP(3,1) + ValR(2)*VecP(2,2)*VecP(3,2) + ValR(3)*VecP(2,3)*VecP(3,3)
               S33out%values(i,j,k)=ValR(1)*VecP(3,1)**2   + ValR(2)*VecP(3,2)**2   + ValR(3)*VecP(3,3)**2
            end do
         end do
      end do
    else
      return
    endif
#else
  print *,'[ERROR] You use LAPACK but you did not use the flag -DLAPACK at the compilation'
#endif

  success = .true.

 end function computeSplitingSij_lapack_field



!!!! tools FOr MHD
!------------------------------------------------------------------------------
!> @author 
!> Mouloud Kessar 
!>
!> @details
!> 
!> @param [in] 
!> @param [inout] 
!------------------------------------------------------------------------------

subroutine ComputeDivTijVec1Vec2(div1,div2,div3,T12,T13,T23,wave,res)

 !I/O
  type(complex_data_layout),intent(inout) :: div1,div2,div3
  type(WaveNumbers), intent(in)           :: wave
  logical, intent(inout)                  :: res
  type(real_data_layout), intent(in)      :: T23,T12,T13
  !Local data
  type(complex_data_layout)          :: T23k,T12K,T13K
  type(complex_data_layout)          :: Tab1k,Tab2K,Tab3K
    
  if (.not. copyStructOnly(div1,T12K) .or. &
     &.not. copyStructOnly(div1,T13K) .or. &
     &.not. copyStructOnly(div1,T23K) .or. &
     &.not. copyStructOnly(div1,Tab1K) .or. &
     &.not. copyStructOnly(div1,Tab2K) .or. &
     &.not. copyStructOnly(div1,Tab3K)  ) then
       res = .false.
       write(6,'(a,i0)')'[ERROR] in ComputeDivTijVec1Vec2 : not enought memory !'
  endif

  call ftran(T12,T12K,res)
  call ftran(T13,T13K,res)
  call ftran(T23,T23K,res)

  Tab1k%values = 0.0_WP  
  call computeDivergenceK(wave,tab1k,T12K,T13K,div1)

  Tab2k%values = -T12K%values
  call computeDivergenceK(wave,tab2k,tab1k,T23K,div2)

  Tab2k%values = -T13K%values  
  Tab3K%values = -T23K%values  
  call computeDivergenceK(wave,Tab2k,Tab3K,tab1k,div3)

  call deleteDataLayout(T12K)
  call deleteDataLayout(T13K)
  call deleteDataLayout(T23K)
  call deleteDataLayout(Tab1K)
  call deleteDataLayout(Tab2K)
  call deleteDataLayout(Tab3K)

end subroutine ComputeDivTijVec1Vec2

!------------------------------------------------------------------------------
!> @author 
!>
!> @details
!> 
!> @param [in] Uin,Vin,Win,Zin
!> @param [in] sizeFilt
!> @param [in] filter
!> @param [inout] res
!> @param [in,out] Ufilout,Vfilout,Wfilout,Zfilout
!> @param [in,out] Ukfilout,Vkfilout,Wkfilout,Zkfilout
!> @param [out] success .TRUE. if all is succcessfull, .FALSE. otherwize.
!------------------------------------------------------------------------------

subroutine computeLijMHD(Uin,Vin,Win,Uinf,Vinf,Winf,Bxin,Byin,Bzin,Bxinf,Byinf,Bzinf&
                     &,Bxk,L12,L13,L23,delta,wave,filter)
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Uin,Vin,Win
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Bxin,Byin,Bzin
!   TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: t12,t13,t23
  
  TYPE(REAL_DATA_LAYOUT), INTENT(INOUT) :: Uinf,Vinf,Winf
  TYPE(REAL_DATA_LAYOUT), INTENT(INOUT) :: Bxinf,Byinf,Bzinf
  TYPE(REAL_DATA_LAYOUT), INTENT(INOUT) :: L12,L13,L23 
  TYPE(COMPLEX_DATA_LAYOUT) :: Bxk

  TYPE(REAL_DATA_LAYOUT) :: tab,tabf
  TYPE(COMPLEX_DATA_LAYOUT) :: tabk,tabkf
  TYPE(WAVENUMBERS), intent(in) :: wave

  REAL(WP), INTENT(IN) :: delta
  logical :: res
  integer :: filter


  IF ( (.NOT. copyStructOnly(Bxk,tabk)) .OR. &
      &(.NOT. copyStructOnly(Bxk,tabkf)) .OR. &
      &(.NOT. copyStructOnly(Bxin,tab)) .OR. &
      &(.NOT. copyStructOnly(Bxin,tabf)) ) RETURN 

  tab%values = Bxin%values*Vin%values - Byin%values*Uin%values 
     CALL ftran(tab,tabk,res)
     CALL computeFilter(wave,2.0*delta,tabk,tabkf,Filter) 
     CALL btran(tabkf,tabf,res)
  L12%values = tabf%values - (Bxinf%values*Vinf%values - Byinf%values*Uinf%values)
  
  tab%values = Bxin%values*Win%values - Bzin%values*Uin%values 
     CALL ftran(tab,tabk,res)
     CALL computeFilter(wave,2.0*delta,tabk,tabkf,Filter) 
     CALL btran(tabkf,tabf,res)
  L13%values = tabf%values - (Bxinf%values*Winf%values - Bzinf%values*Uinf%values)

  tab%values = Byin%values*Win%values - Bzin%values*Vin%values 
     CALL ftran(tab,tabk,res)
     CALL computeFilter(wave,2.0*delta,tabk,tabkf,Filter) 
     CALL btran(tabkf,tabf,res)
  L23%values = tabf%values - (Byinf%values*Winf%values - Bzinf%values*Vinf%values)

  call deleteDataLayout(tab)
  call deleteDataLayout(tabf)
  call deleteDataLayout(tabk)
  call deleteDataLayout(tabkf)

end subroutine computeLijMHD

 subroutine computeTensorJdens1(in12,in13,in23,norm,res)

   !I/O
   type(real_data_layout),intent(in)     :: in23,in12,in13
   type(real_data_layout),intent(inout)  :: norm 
   logical, intent(inout)                :: res  
   
   res = .false.

   norm%values = sqrt(4.0_WP *(  in12%values**2 + &
               & in13%values**2 + in23%values**2) )

   res = .true.

 end subroutine computeTensorJdens1

 subroutine computeJijTensor(in1K,in2K,in3K,wave,out12,out13,out23, &
            & success)

   implicit none

   !I/O
   type(complex_data_layout),intent(in)  :: in1K,in2K,in3K
   type(real_data_layout),intent(inout)  :: out12,out13,out23
   type(WaveNumbers),intent(in)          :: wave
   logical, intent(inout)                :: success  

   !Local data
   type(complex_data_layout)    :: dxin1K,dxin2K,dxin3K
   type(complex_data_layout)    :: dyin1K,dyin2K,dyin3K
   type(complex_data_layout)    :: dzin1K,dzin2K,dzin3K
   type(real_data_layout)       :: dyin1,dzin1,dxin2
   type(real_data_layout)       :: dzin2,dxin3,dyin3
   logical                      :: res  


   success = .true.

   if (.not. copyStructOnly(in1K,dxin1K) .or. &
      &.not. copyStructOnly(in1K,dyin1K) .or. &
      &.not. copyStructOnly(in1K,dzin1K) .or. &
      &.not. copyStructOnly(in1K,dxin2K) .or. &
      &.not. copyStructOnly(in1K,dyin2K) .or. &
      &.not. copyStructOnly(in1K,dzin2K) .or. &
      &.not. copyStructOnly(in1K,dxin3K) .or. &
      &.not. copyStructOnly(in1K,dyin3K) .or. &
      &.not. copyStructOnly(in1K,dzin3K) )then 
       success = .false.
       write(6,'(a)')'[ERROR] in computeStressTensor : not enought memory !'
   endif 

   if (.not. copyStructOnly(out12,dyin1) .or. &
      &.not. copyStructOnly(out12,dzin1) .or. &
      &.not. copyStructOnly(out12,dxin2) .or. &
      &.not. copyStructOnly(out12,dzin2) .or. &
      &.not. copyStructOnly(out12,dxin3) .or. &
      &.not. copyStructOnly(out12,dyin3) )then
       success = .false.
       write(6,'(a)')'[ERROR] in computeStressTensor : not enought memory !'
   endif 

   call computeGradientK(wave,in1K,dxin1K,dyin1K,dzin1K)
   call computeGradientK(wave,in2K,dxin2K,dyin2K,dzin2K)
   call computeGradientK(wave,in3K,dxin3K,dyin3K,dzin3K) 

   call btran(dyin1K,dyin1,res)
   if ( .not. res) success = .false.

   call btran(dzin1K,dzin1,res)
   if ( .not. res) success = .false.

   call btran(dxin2K,dxin2,res)
   if ( .not. res) success = .false.

   call btran(dzin2K,dzin2,res)
   if ( .not. res) success = .false.

   call btran(dxin3K,dxin3,res)
   if ( .not. res) success = .false.

   call btran(dyin3K,dyin3,res)
   if ( .not. res) success = .false.

   out12%values = 0.5_WP*(dyin1%values - dxin2%values)
   out13%values = 0.5_WP*(dzin1%values - dxin3%values)
   out23%values = 0.5_WP*(dzin2%values - dyin3%values)

   call deleteDataLayout(dxin1K)
   call deleteDataLayout(dyin1K)
   call deleteDataLayout(dzin1K)
   call deleteDataLayout(dxin2K)
   call deleteDataLayout(dyin2K)
   call deleteDataLayout(dzin2K)
   call deleteDataLayout(dxin3K)
   call deleteDataLayout(dyin3K)
   call deleteDataLayout(dzin3K)
   call deleteDataLayout(dyin1)
   call deleteDataLayout(dzin1)
   call deleteDataLayout(dxin2)
   call deleteDataLayout(dzin2)
   call deleteDataLayout(dxin3)
   call deleteDataLayout(dyin3)

 end subroutine computeJijTensor


!------------------------------------------------------------------------------
!> @author 
!> Mouloud KESSAR, LEGI
!>
!> @details
!> Compute the 6 components of stress tensor 
!> @param [in] wave the WaveNumbers associated to inK variable
!> @param [in] delta size of filter 
!> @param [in] in1 first component in physical space 
!> @param [in] in2 second component in physical space
!> @param [in] in3 third component in physical space 
!> @param [in] in1f first component filtered in physical space 
!> @param [in] in2f second component filtered in physical space
!> @param [in] in3f third component filtered in physical space 
!> @param [in,out] out11 calculated term in physical space 
!> @param [in,out] out22 calculated term in physical space 
!> @param [in,out] out33 calculated term in physical space 
!> @param [in,out] out12 calculated term in physical space 
!> @param [in,out] out13 calculated term in physical space 
!> @param [in,out] out23 calculated term in physical space 
!> @param [out] success .TRUE. if all is succcessfull, .FALSE. otherwize.
!------------------------------------------------------------------------------
 subroutine computeT_ijVecA(in1,in2,in3,in1f,in2f,in3f,wave,out11,out12,out13, &
            & out22,out23,out33,filter,success,delta)

   implicit none

   !I/O
   type(real_data_layout),intent(inout)     :: out11,out12,out13
   type(real_data_layout),intent(inout)     :: out22,out23,out33
   type(real_data_layout),intent(in)     :: in1,in2,in3
   type(real_data_layout),intent(in)     :: in1f,in2f,in3f
   type(WaveNumbers),intent(in)          :: wave
   logical, intent(inout)                :: success  
   integer, intent(in) :: filter
   REAL(WP) :: delta
   !Local data
   type(real_data_layout)       :: Wtab


   success = .false.

   if (.not. copyStructOnly(out11,Wtab)   )then
       success = .false.
       write(6,'(a)')'[ERROR] in computeT_ijVecA : not enought memory !'
   endif 

   !calcul de out11
   Wtab%values = in1%values**2
   call computeFilter(Wave, delta, Wtab, out11, filter)
   out11%values= out11%values - in1f%values**2

   !calcul de out22
   Wtab%values = in2%values**2
   call computeFilter(Wave, delta, Wtab, out22 , filter)
   out22%values= out22%values - in2f%values**2
   
   !calcul de out33
   Wtab%values = in3%values**2
   call computeFilter(Wave, delta, Wtab, out33 , filter)
   out33%values= out33%values - in3f%values**2
   
   !calcul de out12
   Wtab%values = in1%values*in2%values
   call computeFilter(Wave, delta, Wtab, out12 , filter)
   out12%values= out12%values - in1f%values*in2f%values
   
   !calcul de out13
   Wtab%values = in1%values*in3%values
   call computeFilter(Wave, delta, Wtab, out13 , filter)
   out13%values= out13%values - in1f%values*in3f%values
   
   !calcul de out23
   Wtab%values = in2%values*in3%values
   call computeFilter(Wave, delta, Wtab, out23 , filter)
   out23%values= out23%values - in2f%values*in3f%values

   call deleteDataLayout(Wtab)
   success = .true.

 end subroutine computeT_ijVecA


!------------------------------------------------------------------------------
!> @author 
!> Mouloud KESSAR, LEGI
!>
!> @details
!> Compute the 6 components of stress tensor 
!> @param [in] in1K first component in Fourier space
!> @param [in] in2K second component in Fourier space
!> @param [in] in3K third component in Fourier space
!> @param [in] wave the WaveNumbers associated to inK variable
!> @param [in,out] out11 calculated term in physical space 
!> @param [in,out] out22 calculated term in physical space 
!> @param [in,out] out33 calculated term in physical space 
!> @param [in,out] out12 calculated term in physical space 
!> @param [in,out] out13 calculated term in physical space 
!> @param [in,out] out23 calculated term in physical space 
!> @param [out] success .TRUE. if all is succcessfull, .FALSE. otherwize.
!------------------------------------------------------------------------------
 subroutine computeT_ijVecAVec_B(Ain1K,Ain1,Ain2,Ain3,Ain1f,Ain2f,Ain3f,&
                                &Bin1,Bin2,Bin3,Bin1f,Bin2f,Bin3f,&
                                &wave,out12,out13,out23,filter,success,delta)

   implicit none

   !I/O
   type(complex_data_layout),intent(in)  :: Ain1K
   type(real_data_layout),intent(inout)     :: out12,out13
   type(real_data_layout),intent(inout)     :: out23
   type(real_data_layout),intent(in)     :: Ain1,Ain2,Ain3
   type(real_data_layout),intent(in)     :: Ain1f,Ain2f,Ain3f
   type(real_data_layout),intent(in)     :: Bin1,Bin2,Bin3
   type(real_data_layout),intent(in)     :: Bin1f,Bin2f,Bin3f
   type(WaveNumbers),intent(in)          :: wave
   logical, intent(inout)                :: success  
   integer, intent(in) :: filter
   REAL(WP) :: delta
   !Local data
   type(complex_data_layout)    :: Wtabk,Wtabkf
   type(real_data_layout)       :: Wtab
   logical                      :: res  


   success = .true.

   if (.not. copyStructOnly(Ain1K,WtabK).OR.&
      &.not. copyStructOnly(Ain1K,Wtabkf)  )then 
       success = .false.
       write(6,'(a)')'[ERROR] in computeT_ijVecAVecB : not enought memory !'
   endif 

   if (.not. copyStructOnly(out12,Wtab)   )then
       success = .false.
       write(6,'(a)')'[ERROR] in computeT_ijVecAVecB : not enought memory !'
   endif 



   !calcul de out12
   Wtab%values = Ain1%values*Bin2%values - Ain2%values*Bin1%values
   call ftran(Wtab,Wtabk,res)
   call computeFilter(Wave, delta, Wtabk, Wtabkf, filter)
   call btran(Wtabkf,out12,res)
   out12%values= out12%values - Ain1f%values*Bin2f%values + Ain2f%values*Bin1f%values
   
   if ( .not. res) success = .false.

   !calcul de out13
   Wtab%values = Ain1%values*Bin3%values - Ain3%values*Bin1%values
   call ftran(Wtab,Wtabk,res)
   call computeFilter(Wave, delta, Wtabk, Wtabkf, filter)
   call btran(Wtabkf,out13,res)
   out13%values= out13%values - Ain1f%values*Bin3f%values + Ain3f%values*Bin1f%values
   
   if ( .not. res) success = .false.

   !calcul de out23

   Wtab%values = Ain2%values*Bin3%values - Ain3%values*Bin2%values
   call ftran(Wtab,Wtabk,res)
   call computeFilter(Wave, delta, Wtabk, Wtabkf, filter)
   call btran(Wtabkf,out23,res)
   out23%values= out23%values - Ain2f%values*Bin3f%values + Ain3f%values*Bin2f%values
 
   if ( .not. res) success = .false.

   call deleteDataLayout(Wtab)
   call deleteDataLayout(Wtabk)
   call deleteDataLayout(Wtabkf)

 end subroutine computeT_ijVecAVec_B

!------------------------------------------------------------------------------
!> @author 
!>
!> @details
!> 
!> @param [in] Uin,Vin,Win,Zin
!> @param [in] sizeFilt
!> @param [in] filter
!> @param [inout] res
!> @param [in,out] Ufilout,Vfilout,Wfilout,Zfilout
!> @param [in,out] Ukfilout,Vkfilout,Wkfilout,Zkfilout
!> @param [out] success .TRUE. if all is succcessfull, .FALSE. otherwize.
!------------------------------------------------------------------------------

subroutine ComputeDivSymTijVec1Vec2(div1,div2,div3,T12,T13,T23,T11,T22,T33,wave,res)

 !I/O
  type(complex_data_layout),intent(inout) :: div1,div2,div3
  type(WaveNumbers), intent(in)           :: wave
  logical, intent(inout)                  :: res
  type(real_data_layout), intent(in)      :: T23,T12,T13,T11,T22,T33
  !Local data
  type(complex_data_layout)          :: T23k,T12K,T13K
  type(complex_data_layout)          :: Tab1k

  res = .false.    
  if (.not. copyStructOnly(div1,T12K) .or. &
     &.not. copyStructOnly(div1,T13K) .or. &
     &.not. copyStructOnly(div1,T23K) .or. &
     &.not. copyStructOnly(div1,Tab1K)) then
       write(6,'(a,i0)')'[ERROR] in ComputeDivTijVec1Vec2 : not enought memory !'
       return
  endif

  call ftran(T12,T12K,res)
  call ftran(T13,T13K,res)
  call ftran(T23,T23K,res)
  call ftran(T11,Tab1k,res)
  call computeDivergenceK(wave,tab1k,T12K,T13K,div1)

  call ftran(T22,Tab1k,res)
  call computeDivergenceK(wave,T12K,tab1k,T23K,div2)

  call ftran(T33,Tab1k,res)
  call computeDivergenceK(wave,T13K,T23K,tab1k,div3)

  call deleteDataLayout(T12K)
  call deleteDataLayout(T13K)
  call deleteDataLayout(T23K)
  call deleteDataLayout(Tab1K)

  res = .true.

end subroutine ComputeDivSymTijVec1Vec2

!------------------------------------------------------------------------------
!> @author 
!>
!> @details
!> 
!> @param [in] Uin,Vin,Win,Zin
!> @param [in] sizeFilt
!> @param [in] filter
!> @param [inout] res
!> @param [in,out] Ufilout,Vfilout,Wfilout,Zfilout
!> @param [in,out] Ukfilout,Vkfilout,Wkfilout,Zkfilout
!> @param [out] success .TRUE. if all is succcessfull, .FALSE. otherwize.
!------------------------------------------------------------------------------

 subroutine interfVelSca(Uin,Vin,Win,Zin,sizeFilt,filter,res,&
                        &Ufilout,Vfilout,Wfilout,Zfilout,&
                        &Ukfilout,Vkfilout,Wkfilout,Zkfilout,WN)

   use datalayout
   use mpilayout_tools , only : getnbmycpu
   use transforms_tools, only : ftran,btran
   use wavenumber_tools

   !I/O Data
   type(real_data_layout),intent(in)                :: Uin,Vin,Win,Zin
   integer,intent(in)                               :: sizeFilt,filter
   type(real_data_layout),intent(inout),optional    :: Ufilout,Vfilout,Wfilout,Zfilout
   type(complex_data_layout),intent(inout),optional :: Ukfilout,Vkfilout,Wkfilout,Zkfilout
   type(wavenumbers),intent(inout),optional         :: WN 
   logical,intent(inout)                            :: res

   !Local data
   integer                                          :: me
   type(wavenumbers)                                :: wave
   

   res=.true.
   me = getnbmycpu()
   !Calcul des nombres d'onde
   if (.not. initWN(wave,me,Uin%nx,Uin%ny,Uin%nz) .or. &
      &.not.  computeWN(wave,me,Uin%Lx,Uin%Ly,Uin%Lz,Uin%nx,Uin%ny,Uin%nz) ) res = .false.
   if (.not. samelayout(Uin,Vin,Win,Zin) ) res = .false.

   if (present(Ufilout)) then
     call computeFilter(wave,sizeFilt*computedelta(Uin),Uin,Ufilout,filter)
   endif
   if (present(Vfilout)) then
     call computeFilter(wave,sizeFilt*computedelta(Vin),Vin,Vfilout,filter)
   endif
   if (present(Wfilout)) then
     call computeFilter(wave,sizeFilt*computedelta(Win),Win,Wfilout,filter)
   endif
   if (present(Zfilout)) then
     call computeFilter(wave,sizeFilt*computedelta(Zin),Zin,Zfilout,filter)
   endif
   if (present(Ukfilout)) then
     call computeFilter(wave,sizeFilt*computedelta(Uin),Uin,Ukfilout,filter)
   endif
   if (present(Vkfilout)) then
     call computeFilter(wave,sizeFilt*computedelta(Vin),Vin,Vkfilout,filter)
   endif
   if (present(Wkfilout)) then
     call computeFilter(wave,sizeFilt*computedelta(Win),Win,Wkfilout,filter)
   endif
   if (present(Zkfilout)) then
     call computeFilter(wave,sizeFilt*computedelta(Zin),Zin,Zkfilout,filter)
   endif
   if (present(WN)) then
     if (.not. initWN(WN,me,Uin%nx,Uin%ny,Uin%nz) .or. &
        &.not.  computeWN(WN,me,Uin%Lx,Uin%Ly,Uin%Lz,Uin%nx,Uin%ny,Uin%nz) ) res = .false.
   endif
   if(.not. deleteWN(wave)) res = .false.

 end subroutine interfVelSca

end module subgrid_tools
!> @}
