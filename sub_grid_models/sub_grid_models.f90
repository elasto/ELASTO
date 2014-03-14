!> @addtogroup sub_gird_models
!! @{

!------------------------------------------------------------------------------
!
! MODULE: subgridmodels 
!
!> @author
!> Antoine Vollant
!
! DESCRIPTION: 
!> The aim of this module is to select the model for the right field
!> 
!------------------------------------------------------------------------------
 

 
module subgridmodels

 use mhd_model
 use datalayout
 use precision_tools 
 use param
 use scamodels
 use velmodels
 use wavenumber_tools
 use conditional_mean_tools 
 use opt_estim_models
 use Eqdm_lorentz_models
 implicit none

 !Number of velocity model to use
 integer, private                                        :: numVel
 !Number of filter type for velocity 
 integer, private                                        :: numVelFilter
 !Type of averaging for dynamic proc for vel
 integer, private                                        :: type_avg
 !Number of scalar model to use 
 integer, dimension(:), pointer, private                 :: numSca
 !Number of filter type for scalars 
 integer, dimension(:), pointer, private                 :: numScaFilter
 !Coef for model for scalars 
 real(WP), dimension(:), pointer, private                :: coef_sca
 !Type of averaging for dynamic proc for scalars 
 integer, dimension(:), pointer, private                :: type_avg_sca
 !Number of filter type for scalars 
 character(len=str_long), dimension(:), pointer, private:: nameEOFile
 !Number of induction model to use
 integer, private                                        :: numInduc
 !Number of filter type for induction 
 integer, private                                        :: numInducFilter
 !Type of averaging for dynamic proc for induction
 integer, private                                        :: type_avg_induc
  integer, private,save :: iterationBmodel
 integer,private                                         :: ite_loc_vel

real(WP),private,save :: C_mu=0.046_WP
real(WP),private,save :: C_lambda=(5/7*0.046_WP)


 contains

!---------------------------------------------------------------------------
!> @details 
!> Initialize the context of LES.
!! @author = Antoine Vollant, LEGI
!!    @param[in]    nbscal is the number of scalars
!!    @param[in]    U is a component of velocity in physical space
!!    @param[in]    Scal is the array which contains all scalars in physical space
!!    @param[in]    Uk is a component of velocity in spectral space
!!    @param[in]    Scalk is the array which contains all scalars in spectral space
!!    @param[in]    spec_rank is mpi processus rank
!!    @return       res is logical value equal to TRUE if the initialization is
!!                  Ok, false in the other case
!---------------------------------------------------------------------------

 function init_LES(nbscal,U,Scal,Uk,ScalK,spec_rank) result (res)
 
  use opt_estim_models , only : read_models_param 
 
  !I/O
  integer, intent(in)                                :: nbscal
  integer, intent(in)                                :: spec_rank 
  type(complex_data_layout),intent(in)               :: Uk
  type(complex_data_layout),dimension(:),intent(in)  :: ScalK 
  type(real_data_layout),intent(in)                  :: U
  type(real_data_layout),dimension(:),intent(in)     :: Scal 
  logical                                            :: res

  !Local data
  integer                                            :: i,ires
  character(str_long)                                :: request
  character(str_long)                                :: nameles 
  character(str_long)                                :: namelessc

  res = .false.

  ! make the link between the name of velocity model and the number of model
  if (parser_is_defined('LES model for velocity')) then
    call parser_read('LES model for velocity',nameles)
    if ( nameles .eq. 'smagorinsky') then
      numVel = 1
    elseif ( nameles .eq. 'dynsmagorinsky') then
      if (.not. parseFilter(numVelFilter,spec_rank)) return
      write(request,'(a,i0)') 'Type of averaging for velocity'
      call parser_read(request,type_avg)
      numVel = 2
    elseif ( nameles .eq. 'sf') then
      numVel = 3
    elseif ( nameles .eq. 'fsf') then
      numVel = 4
    elseif ( nameles .eq. 'gradient') then
      numVel = 5
    elseif ( nameles .eq. 'Similarity') then
      numVel = 6
    elseif ( nameles .eq. 'gradientIncomp') then
      numVel = 7
    elseif ( nameles .eq. 'RGMVel') then
      numVel = 8
    elseif ( nameles .eq. 'gradientvort') then
      numVel = 9
    elseif ( nameles .eq. 'rgmdivdynv1') then
      if (.not. parseFilter(numVelFilter,spec_rank)) return
      numVel = 10
    elseif ( nameles .eq. 'rgmdivdynv2') then
      if (.not. parseFilter(numVelFilter,spec_rank)) return
      numVel = 11
    elseif ( nameles .eq. 'rgmoikojk') then
      numVel = 12
    elseif ( nameles .eq. 'rgmsikmoinskjoikokj') then
      numVel = 13
    elseif ( nameles .eq. 'rgmsikmoinsskj') then
      numVel = 14
    elseif ( nameles .eq. 'rgmsikmoinsskjsikojkoiksjk') then
      numVel = 15
    elseif ( nameles .eq. 'rgmsikmoinsskjsikojkoiksjkdyn') then
      if (.not. parseFilter(numVelFilter,spec_rank)) return
      numVel = 16
    elseif ( nameles .eq. 'rgmsikmoinsskjdyn') then
      if (.not. parseFilter(numVelFilter,spec_rank)) return
      numVel = 17
    elseif ( nameles .eq. 'rgmsikmoinsskjsikojkoiksjkdynInv') then
      if (.not. parseFilter(numVelFilter,spec_rank)) return
      numVel = 18
    elseif ( nameles .eq. 'rgmsikmoinsskjdynfab') then
      if (.not. parseFilter(numVelFilter,spec_rank)) return
      numVel = 19
    elseif ( nameles .eq. 'rgmsikmoinsskjsikojkoiksjkdynInvFab') then
      if (.not. parseFilter(numVelFilter,spec_rank)) return
      numVel = 20
    elseif ( nameles .eq. 'rgmsikmoinsskjsikojkoiksjkdynFab') then
      if (.not. parseFilter(numVelFilter,spec_rank)) return
      numVel = 21
    elseif ( nameles .eq. 'rgmsikmoinsskjdynYouMoin') then
      if (.not. parseFilter(numVelFilter,spec_rank)) return
      numVel = 22
    elseif ( nameles .eq. 'Hamba') then
!       if (.not. parseFilter(numVelFilter,spec_rank)) return
      numVel = 23
      C_mu=0.046_WP
      C_lambda=5/7*C_mu
    elseif ( nameles .eq. 'Carrati') then
      numVel = 24

    else
      numVel = 0
    endif
  else
    numVel = 0
  endif
  ite_loc_vel = 0
  !if (spec_rank.EQ.0) print*, 'numvel:', numvel
  ! make the link between the name of scalar model and the number of model
  if (nbscal .gt. 0) then
    call init_opt_estim_models
    allocate(numSca(nbscal),stat=ires)
    if (ires .ne. 0) return
    allocate(numScaFilter(nbscal),stat=ires)
    if (ires .ne. 0) return
    allocate(coef_sca(nbscal),stat=ires)
    if (ires .ne. 0) return
    allocate(type_avg_sca(nbscal),stat=ires)
    if (ires .ne. 0) return
    allocate(nameEOFile(nbscal),stat=ires)
    if (ires .ne. 0) return
    do i=1,nbscal
      !read options for scalars' model
      write(request,'(a,i0)') 'LES model for scalar ',i
      call parser_read(request,namelessc)
      if (spec_rank .eq. 0 .and. namelessc .ne. "0") then 
        write (6,'(a,1x,i2,a,1x,a20)') '[INFO] model asked for',i,'th scalar',namelessc
      endif
      if ( namelessc .eq. 'smagorinsky') then
        numSca(i) = 1
      elseif ( (namelessc .eq. 'dynsmagorinsky_3d')) then
        if (.not. parseFilter(numScaFilter(i),spec_rank,sca=i)) return
        numSca(i) = 2
      elseif ( namelessc .eq. 'gradient') then
        write(request,'(a,i0)') 'Model coefficient for scalar ',i
        call parser_read(request,coef_sca(i))
        numSca(i) = 3 
      elseif ( namelessc .eq. 'RGM') then
        write(request,'(a,i0)') 'Model coefficient for scalar ',i
        call parser_read(request,coef_sca(i))
        numSca(i) = 5 
        numScaFilter(i) = 0
      elseif ( namelessc .eq. 'dynRGM') then
        if (.not. parseFilter(numScaFilter(i),spec_rank,sca=i)) return
        write(request,'(a,i0)') 'Type of averaging for scalar ',i
        call parser_read(request,type_avg_sca(i))
        numSca(i) = 6 
      elseif ( namelessc .eq. 'RGMO') then
        write(request,'(a,i0)') 'Model coefficient for scalar ',i
        call parser_read(request,coef_sca(i))
        numSca(i) = 7 
        numScaFilter(i) = 0
      elseif ( namelessc .eq. 'dynWANG') then
        if (.not. parseFilter(numScaFilter(i),spec_rank,sca=i)) return
        numSca(i) = 8 
      elseif ((namelessc .eq. 'dynsmagorinsky_2d_z')) then
        if (.not. parseFilter(numScaFilter(i),spec_rank,sca=i)) return
        numSca(i) = 9 
      elseif ( namelessc .eq. 'dynclark_3d') then
        if (.not. parseFilter(numScaFilter(i),spec_rank,sca=i)) return
        numSca(i) = 10 
      elseif ( namelessc .eq. 'dynclarkfabre_3d') then
        if (.not. parseFilter(numScaFilter(i),spec_rank,sca=i)) return
        numSca(i) = 11
      elseif ( namelessc .eq. 'dynclark_2d_z') then
        if (.not. parseFilter(numScaFilter(i),spec_rank,sca=i)) return
        numSca(i) = 12 
      elseif ( namelessc .eq. 'dynclarkfabre_2d_z') then
        if (.not. parseFilter(numScaFilter(i),spec_rank,sca=i)) return
        numSca(i) = 13
      elseif ( namelessc .eq. 'RGM_analyt') then
        write(request,'(a,i0)') 'Model coefficient for scalar ',i
        call parser_read(request,coef_sca(i))
        numSca(i) = 14
      elseif ( namelessc .eq. 'smagdyn_verif') then
        if (.not. parseFilter(numScaFilter(i),spec_rank,sca=i)) return
        write(request,'(a,i0)') 'Type of averaging for scalar ',i
        call parser_read(request,type_avg_sca(i))
        numSca(i) = 15
      elseif ( namelessc .eq. 'dynsmagorinsky') then
        if (.not. parseFilter(numScaFilter(i),spec_rank,sca=i)) return
        write(request,'(a,i0)') 'Type of averaging for scalar ',i
        call parser_read(request,type_avg_sca(i))
        numSca(i) = 16
      elseif ( namelessc .eq. 'dynclark') then
        if (.not. parseFilter(numScaFilter(i),spec_rank,sca=i)) return
        write(request,'(a,i0)') 'Type of averaging for scalar ',i
        call parser_read(request,type_avg_sca(i))
        numSca(i) = 17
      elseif ( namelessc .eq. 'dynclarkfabre') then
        if (.not. parseFilter(numScaFilter(i),spec_rank,sca=i)) return
        write(request,'(a,i0)') 'Type of averaging for scalar ',i
        call parser_read(request,type_avg_sca(i))
        numSca(i) = 18
      elseif ( namelessc .eq. 'scalesimilarity') then
        if (.not. parseFilter(numScaFilter(i),spec_rank,sca=i)) return
        numSca(i) = 19
      elseif ( namelessc .eq. 'optimalestimator') then
        if (.not. read_models_param(i) ) return
        numSca(i) = 20
      elseif ( namelessc .eq. "0" ) then
        if (spec_rank .eq. 0) write (6,'(a,1x,i2,a,1x,a20)') '[INFO] NO MODEL  asked for',i,'th scalar',namelessc
      else
        if (spec_rank .eq. 0 ) then 
          write (6,'(a,1x,i2,a,1x,a20)') '[INFO] the model asked for',i,'th scalar is unknown',namelessc
          return
        endif
        numSca(i) = 0
      endif
    enddo
  endif
  iterationBmodel=0

 ! make the link between the name of velocity model and the number of model
  if (parser_is_defined('LES model for Induction')) then
    call parser_read('LES model for Induction',nameles)

    if ( nameles .eq. 'smagorinsky') then
      numInduc = 1
    elseif ( nameles .eq. 'gradient') then
      numInduc = 2
    elseif ( nameles .eq. 'similarity') then
      numInduc = 3
    elseif ( nameles .eq. 'dynsmagorinsky') then
      numInduc = 4
    elseif ( nameles .eq. 'Hamba') then
      numInduc = 5
    elseif ( nameles .eq. 'Carrati') then
      numInduc = 6

    else
      numInduc = 0
    endif
  else
    numInduc = 0
  endif
 
  res = .true.

 end function init_LES


!---------------------------------------------------------------------------
!> @details 
!> This subroutine gives an interface for the different velocity models.
!! @author = Antoine Vollant, LEGI
!>    @param [in]    U  longitudinal velocity in physical space
!>    @param [in]    V  velocity in physical space
!>    @param [in]    W  velocity in physical space
!>    @param [in]    Uk longitudinal velocity in spectral space
!>    @param [in]    Vk velocity in spectral space
!>    @param [in]    Wk velocity in spectral space
!>    @param [inout] nlx non linear term for velocity 
!>    @param [inout] nly non linear term for velocity 
!>    @param [inout] nlz non linear term for velocity 
!>    @param [in]    VelWN wave numbers for velocity 
!>    @param [in]    spec_rank is the number of processus in cpu pool 
!!    @return       res is logical value equal to TRUE if the initialization is
!!                  Ok, false in the other case
!---------------------------------------------------------------------------

 subroutine computeVelocityModel(U,V,W,Uk,Vk,Wk,Bxk,Byk,Bzk,nlx,nly,nlz,VelWN,spec_rank,res,lorentz)

  !I/O
  type(complex_data_layout),intent(in)    :: Uk,Vk,Wk 
  type(complex_data_layout),intent(in)    :: Bxk,Byk,Bzk
  type(complex_data_layout),intent(inout) :: nlx,nly,nlz
  type(real_data_layout),intent(in)       :: U,V,W
  type(WaveNumbers), intent(in)           :: VelWN
  integer, intent(in)                     :: spec_rank 
  logical, intent(inout)                  :: res
  logical,intent(in),optional             :: lorentz
  !Local data
  type(complex_data_layout) :: T1vel,T2vel,T3vel 
!  type(real_data_layout)    :: temp0,temp1,temp2
  REAL(WP) :: deltx,maxDiv,moyDiv,minDiv
  integer :: nbcpus
  logical :: inc_ite_loc_vel=.false.
  !allocation of memory
  if (.not. copyStructOnly(Uk,T1vel) .or. &
     &.not. copyStructOnly(Uk,T2vel) .or. &
!     &.not. copyStructOnly(U,temp0) .or. &
     &.not. copyStructOnly(Uk,T3vel) ) then
       res = .false.
       write(6,'(a,i0)')'[ERROR] in LES for vel : not enought memory !',spec_rank
  endif 
  call parser_read('Get pdf of sgs term for qdm',inc_ite_loc_vel)
  deltx = computeDelta(U) 
  nbcpus = getnbcpus()

  if ( .not. inc_ite_loc_vel ) ite_loc_vel = -1
!  if (.not. computeFieldDivergence(VelWN,U,V,W,temp0,nbcpus,spec_rank) ) return
!  maxDiv=doglobalmax(spec_rank,maxval(temp0%values)) 
!  moyDiv=computeFieldAvg(temp0,spec_rank)
!  minDiv=doglobalmin(spec_rank,minval(temp0%values)) 
!  if (spec_rank.eq.0) print *,'Divergence min-> ',minDiv,' moy-> ',moyDiv,' max-> ',maxDiv

  select case (numVel)
  case(1)
    call smagorinskyVel(T1vel,T2vel,T3vel,U,V,W,Uk,Vk,Wk,VelWN,numVel,spec_rank,res,ite=ite_loc_vel)
  case(2)
    call smagorinskyVel(T1vel,T2vel,T3vel,U,V,W,Uk,Vk,Wk,VelWN,numVel,spec_rank,res,numVelFilter=numVelFilter,type_avg=type_avg,ite=ite_loc_vel)
  case(3)
    call SfVel(T1vel,T2vel,T3vel,U,V,W,Uk,Vk,Wk,VelWN,spec_rank,numVel,res)
  case(4)
    call SfVel(T1vel,T2vel,T3vel,U,V,W,Uk,Vk,Wk,VelWN,spec_rank,numVel,res)
  case(5)
    call gradientVelDIV(T1vel,T2vel,T3vel,U,V,W,VelWN,deltx,spec_rank,nbcpus)
  case(6)
    call SimilVelDIV(T1vel,T2vel,T3vel,U,V,W,Uk,Vk,Wk,VelWN,deltx,numVelFilter,spec_rank,nbcpus)
  case(7)
    call gradientVelIncomp(T1vel,T2vel,T3vel,U,V,W,VelWN,spec_rank,ite=ite_loc_vel) 
  case(8)
    call RGMVel(T1vel,T2vel,T3vel,U,U,V,W,VelWN,spec_rank,signe='moins',ite=ite_loc_vel) 
  case(9)
    call gradientQDMvort(T1vel,T2vel,T3vel,U,V,W,VelWN,spec_rank,ite=ite_loc_vel)
  case(10)
    call RGMDivDynV1(T1vel,T2vel,T3vel,U,U,V,W,VelWN,spec_rank,numVelFilter,ite=ite_loc_vel)
  case(11)
    call RGMDivDynV2(T1vel,T2vel,T3vel,U,U,V,W,VelWN,spec_rank,numVelFilter,ite=ite_loc_vel)
  case(12)
    call RGMOikOjk(T1vel,T2vel,T3vel,U,Uk,Vk,Wk,VelWN,spec_rank,ite=ite_loc_vel)
  case(13)
    call RGMSikmoinsSkjOikOjk(T1vel,T2vel,T3vel,U,Uk,Vk,Wk,VelWN,spec_rank,ite=ite_loc_vel)
  case(14)
    call RGMSikmoinsSkj(T1vel,T2vel,T3vel,U,Uk,Vk,Wk,VelWN,spec_rank,ite=ite_loc_vel)
  case(15)
    call RGMSikmoinsSkjSikOjkOikSjk(T1vel,T2vel,T3vel,U,Uk,Vk,Wk,VelWN,spec_rank,ite=ite_loc_vel)
  case(16)
    call RGMSikmoinsSkjSikOjkOikSjkdyn(T1vel,T2vel,T3vel,U,V,W,VelWN,spec_rank,numVelFilter,ite=ite_loc_vel)
!    if (.not. success) res = .false.
  case(17)
    call RGMSikmoinsSkjDyn(T1vel,T2vel,T3vel,U,V,W,VelWN,spec_rank,numVelFilter,ite=ite_loc_vel)
  case(18)
    call RGMSikmoinsSkjSikOjkOikSjkdyn2(T1vel,T2vel,T3vel,U,V,W,VelWN,spec_rank,numVelFilter,ite=ite_loc_vel)
  case(19)
    call RGMSikmoinsSkjDynFab(T1vel,T2vel,T3vel,U,Uk,Vk,Wk,VelWN,spec_rank,numVelFilter,ite=ite_loc_vel)
  case(20)
    call RGMSikmoinsSkjSikOjkOikSjkDyn2Fab(T1vel,T2vel,T3vel,U,V,W,VelWN,spec_rank,numVelFilter,ite=ite_loc_vel)
  case(21)
    call RGMSikmoinsSkjSikOjkOikSjkDynFab(T1vel,T2vel,T3vel,U,V,W,VelWN,spec_rank,numVelFilter,ite=ite_loc_vel)
  case(22)
    call RGMSikmoinsSkjDynYouMoin(T1vel,T2vel,T3vel,U,V,W,Uk,Vk,Wk,VelWN,spec_rank,numVelFilter,ite=ite_loc_vel)
  case(23)
    call DivHambaYoshizawaNSequation(U,Uk,Vk,Wk,Bxk,Byk,Bzk,C_Lambda,C_mu,VelWN,spec_rank,res, &
                          &  deltx,T1vel,T2vel,T3vel)
  case(24)
    call DivMullerCarratiNSequation(U,Uk,Vk,Wk,Bxk,Byk,Bzk,VelWN,spec_rank,res, &
                          &  deltx,T1vel,T2vel,T3vel)
  end select

  if (present(lorentz)) then
     nlx%values = nlx%values - T1vel%values
     nly%values = nly%values - T2vel%values
     nlz%values = nlz%values - T3vel%values
  else
     nlx%values = nlx%values + T1vel%values
     nly%values = nly%values + T2vel%values
     nlz%values = nlz%values + T3vel%values
  endif

!  call deleteDataLayout(temp0)
  call deleteDataLayout(T1vel)
  call deleteDataLayout(T2vel)
  call deleteDataLayout(T3vel)

  ite_loc_vel = ite_loc_vel + 1

 end subroutine computeVelocityModel

!---------------------------------------------------------------------------
!> @details 
!> This subroutine gives an interface for the different scalar models.
!! @author = Antoine Vollant, LEGI
!!    @param[in]     sca is the number of scalars
!>    @param [in]    U  longitudinal velocity in physical space
!>    @param [in]    V  velocity in physical space
!>    @param [in]    W  velocity in physical space
!>    @param [in]    Uk longitudinal velocity in spectral space
!>    @param [in]    Vk velocity in spectral space
!>    @param [in]    Wk velocity in spectral space
!>    @param [in]    Scal is the sca^th scalar in phycical space
!>    @param [in]    ScalK is the sca^th scalar in spectral space
!>    @param [inout] ScalNL is the non linear term for scalar 
!>    @param [in]    ScalWN are waves number for scalar
!>    @param [in]    spec_rank is the number of processus in cpu pool 
!!    @return       res is logical value equal to TRUE if the initialization is
!!                  Ok, false in the other case
!---------------------------------------------------------------------------
 subroutine computeScalarModel( sca,U,V,W,Uk,Vk,Wk,Scal,Scalk,ScalNL,ScalWN,spec_rank,res )

  !I/O
  integer, intent(in)                      :: sca  !number of scalar considered
  integer, intent(in)                      :: spec_rank
  type(complex_data_layout),intent(in)     :: Uk,Vk,Wk,Scalk
  type(real_data_layout),intent(in)        :: U,V,W,Scal
  type(WaveNumbers),intent(in)             :: ScalWN
  type(complex_data_layout),intent(inout)  :: ScalNL
  logical,intent(inout)                    :: res 
 
  !Local data
  type(complex_data_layout)                :: T1sca
  real(WP), dimension(:),allocatable       :: DynCoef1D

  if (.not. copyStructOnly(Scalk,T1sca) ) then
    write(6,'(a,i0)')'[ERROR] in initialization for scalar LES : not enought memory !',spec_rank
  endif

  select case (numSca(sca))
  case(1)
    call smagorinskySca(T1Sca,U,Uk,Vk,Wk,ScalK,ScalWN,res,spec_rank)
  case(2)
    call DynSmagorinskySca_3d(T1Sca,U,V,W,Scal,Uk,Vk,Wk,ScalK,ScalWN,res,numScaFilter(sca),spec_rank)
  case(3)
    call gradientSca(T1sca,U,Uk,Vk,Wk,ScalK,Scal,ScalWN,res,1,coef_sca(sca))
  case(5)
    call RGMSca_lapack(T1sca,U,Uk,Vk,Wk,ScalK,ScalWN,res,spec_rank,coef_sca(sca))
  case(6)
    call DynRGMSca(T1sca,U,V,W,Uk,Vk,Wk,Scal,ScalK,ScalWN,numScaFilter(sca),type_avg_sca(sca),spec_rank,res)
  case(7)
    call RGMOSca(T1sca,U,Uk,Vk,Wk,ScalK,ScalWN,res,spec_rank,coef_sca(sca))    !!! TO FINISH TRY DYNAMIC PROCEDURE ALSO
  case(8)
    call DynWangSca(T1sca,U,V,W,Uk,Vk,Wk,Scal,ScalK,ScalWN,numScaFilter(sca),spec_rank,res)
  case(9)
    call DynSmagorinskySca_2d_z(T1sca,U,V,W,Scal,Uk,Vk,Wk,ScalK,ScalWN,res,numScaFilter(sca),spec_rank)
  case(10)
    call DynClarkSca_3d(T1sca,U,V,W,Uk,Vk,Wk,Scal,ScalK,ScalWN,numScafilter(sca),spec_rank,res)
  case(11)
    call DynClarkFabreSca_3d(T1Sca,U,V,W,Uk,Vk,Wk,Scal,ScalK,ScalWN,numScafilter(sca),spec_rank,res)
  case(12)
    call DynClarkSca_2d_z(T1sca,U,V,W,Uk,Vk,Wk,Scal,ScalK,ScalWN,numScafilter(sca),spec_rank,res)
  case(13)
    call DynClarkFabreSca_2d_z(T1Sca,U,V,W,Uk,Vk,Wk,Scal,ScalK,ScalWN,numScafilter(sca),spec_rank,res)
  case(14)
    call RGMSca_analyt(T1sca,U,Uk,Vk,Wk,ScalK,ScalWN,res,spec_rank,coef_sca(sca))
!    call RGMScaForTest(T1sca,U,Uk,Vk,Wk,ScalK,ScalWN,res,spec_rank,coef_sca(sca))
  case(15)
    allocate(Dyncoef1D(U%zmin:U%zmax))
    call coefDynSca2(U,V,W,Uk,Vk,Wk,Scal,Scalk,ScalWN,numScaFilter(sca),Dyncoef1D,U%zmin,U%zmax,type_avg_sca(sca),spec_rank,res)
    call smagorinskySca2(T1Sca,U,Uk,Vk,Wk,Scalk,ScalWN,res,spec_rank,Dyncoef1D,U%zmin,U%zmax)
    deallocate(dyncoef1D)
  case(16)
    call DynSmagorinskySca(T1Sca,U,V,W,Scal,Uk,Vk,Wk,ScalK,ScalWN,res,numScaFilter(sca),spec_rank,type_avg_sca(sca))
  case(17)
    call DynClarkSca(T1sca,U,V,W,Uk,Vk,Wk,Scal,ScalK,ScalWN,numScafilter(sca),spec_rank,res,type_avg_sca(sca))
  case(18)
    call DynClarkFabreSca(T1Sca,U,V,W,Uk,Vk,Wk,Scal,ScalK,ScalWN,numScafilter(sca),spec_rank,res,type_avg_sca(sca))
  case(19)
    call ScaleSimilaritySca(T1Sca,U,V,W,Scal,Uk,Vk,Wk,ScalK,ScalWN,res,numScafilter(sca),1,coef_sca(sca))
  case(20)
    call optimalEstimatorModelsSca(T1Sca,Scal,Uk,Vk,Wk,ScalK,ScalWN,sca,res)
  end select

  !Compute the NL term
  ScalNL%values = ScalNL%values - T1sca%values

  call deleteDataLayout(T1Sca)

 end subroutine computeScalarModel

!---------------------------------------------------------------------------
!> @details 
!> This subroutine gives an interface for the different velocity models.
!! @author = Mouloud KESSAR, LEGI
!>    @param [in]    Bx magnetic field in physical space (x)
!>    @param [in]    By magnetic field in physical space (y)
!>    @param [in]    Bz magnetic field in physical space (z)
!>    @param [in]    Bxk magnetic field in spectral space
!>    @param [in]    Byk magnetic field in spectral space
!>    @param [in]    Bzk magnetic field in spectral space
!>    @param [inout] nlx non linear term for velocity 
!>    @param [inout] nly non linear term for velocity 
!>    @param [inout] nlz non linear term for velocity 
!>    @param [in]    BfieldWN wave numbers for velocity 
!>    @param [in]    spec_rank is the number of processus in cpu pool 
!!    @return       res is logical value equal to TRUE if the initialization is
!!                  Ok, false in the other case
!---------------------------------------------------------------------------

 subroutine computeLorentzVelModel(Bx,By,Bz,Bxk,Byk,Bzk,nlx,nly,nlz,BfieldWN,spec_rank,res)

  !I/O
  type(complex_data_layout),intent(in)    :: Bxk,Byk,Bzk 
  type(complex_data_layout),intent(inout) :: nlx,nly,nlz
  type(real_data_layout),intent(in)       :: Bx,By,Bz
  type(WaveNumbers), intent(in)           :: BfieldWN
  integer, intent(in)                     :: spec_rank 
  logical, intent(inout)                  :: res
  logical :: lorentz
  lorentz =.true.
  if (.not.Lorentz) then
   call computeVelocityModel(Bx,By,Bz,Bxk,Byk,Bzk,Bxk,Byk,Bzk,nlx,nly,nlz,BfieldWN,spec_rank,res,lorentz)
  endif
 end subroutine computeLorentzVelModel

!---------------------------------------------------------------------------
!> @details 
!> This subroutine gives an interface for the different velocity models.
!! @author = Mouloud KESSAR, LEGI
!>    @param [in]    Bx magnetic field in physical space (x)
!>    @param [in]    By magnetic field in physical space (y)
!>    @param [in]    Bz magnetic field in physical space (z)
!>    @param [in]    Bxk magnetic field in spectral space
!>    @param [in]    Byk magnetic field in spectral space
!>    @param [in]    Bzk magnetic field in spectral space
!>    @param [in]    U  longitudinal velocity in physical space
!>    @param [in]    V  velocity in physical space
!>    @param [in]    W  velocity in physical space
!>    @param [in]    Uk longitudinal velocity in spectral space
!>    @param [in]    Vk velocity in spectral space
!>    @param [in]    Wk velocity in spectral space
!>    @param [inout] nlx non linear term for velocity 
!>    @param [inout] nly non linear term for velocity 
!>    @param [inout] nlz non linear term for velocity 
!>    @param [in]    BfieldWN wave numbers for B field 
!>    @param [in]    spec_rank is the number of processus in cpu pool 
!!    @return       res is logical value equal to TRUE if the initialization is
!!                  Ok, false in the other case
!---------------------------------------------------------------------------

 subroutine computeBfieldModel(Bx,By,Bz,Bxk,Byk,Bzk,U,V,W,Uk,Vk,Wk,nlx,nly,nlz,BfieldWN,spec_rank,res)

  !I/O
  type(complex_data_layout),intent(in)    :: Bxk,Byk,Bzk 
  type(complex_data_layout),intent(inout) :: nlx,nly,nlz
  type(real_data_layout),intent(in)       :: Bx,By,Bz
  type(WaveNumbers), intent(in)           :: BfieldWN
  type(complex_data_layout),intent(in)    :: Uk,Vk,Wk 
  type(real_data_layout),intent(in)       :: U,V,W
  integer, intent(in)                     :: spec_rank
  logical, intent(inout)                  :: res
  !Local data
  type(complex_data_layout) :: T1Bfield,T2Bfield,T3Bfield 
  real(WP)                  :: deltx!,meanEmTBbmodel
  type(real_data_layout) :: T12,T13,T23 
!   type(real_data_layout) :: J12,J13,J23,EmTBbmodel
  integer                :: nbcpus

  !allocation of memory
  if (.not. copyStructOnly(Bxk,T1Bfield) .or. &
     &.not. copyStructOnly(Bxk,T2Bfield) .or. &
     &.not. copyStructOnly(Bxk,T3Bfield) ) then
       res = .false.
       write(6,'(a,i0)')'[ERROR] in LES for Bfield : not enought memory !',spec_rank
  endif 

  if (.not. copyStructOnly(U,T12) .or. &
     &.not. copyStructOnly(U,T13) .or. &
     &.not. copyStructOnly(U,T23) )then
       res = .false.
       write(6,'(a,i0)')'[ERROR] in LES for Bfield : not enought memory !',spec_rank
   endif

  deltx = computeDelta(U)

  select case (numInduc)
  case(1)
    call smagorinskyMHDTij(T12,T13,T23,Bxk,Byk,Bzk,BfieldWN,spec_rank,res, deltx)

    call ComputeDivTijVec1Vec2(T1Bfield,T2Bfield,T3Bfield,T12,T13,T23,BfieldWN,res)

  case(2)
    call InducGradient(U,V,W,Bx,By,Bz,BfieldWN,deltx,T12,T13,T23,numVelFilter,spec_rank )

    call ComputeDivTijVec1Vec2(T1Bfield,T2Bfield,T3Bfield,T12,T13,T23,BfieldWN,res)

  case(3)  
    call SimilInduc(T12,T13,T23,U,V,W,Uk,Vk,Wk,Bx,By,Bz,Bxk,Byk,Bzk,BfieldWN,deltx,numVelFilter,spec_rank)

    call ComputeDivTijVec1Vec2(T1Bfield,T2Bfield,T3Bfield,T12,T13,T23,BfieldWN,res)

  case(4)  
    call DynsmagorinskyMHDTij(T12,T13,T23,U,V,W,Uk,Vk,Wk,Bxk,Byk,Bzk,Bx,By,Bz,BfieldWN,spec_rank,res, &
                          &  numVelFilter,deltx)
    call ComputeDivTijVec1Vec2(T1Bfield,T2Bfield,T3Bfield,T12,T13,T23,BfieldWN,res)

  case(5)  
    call HambaYoshizawaInduc(Bx,Uk,Vk,Wk,Bxk,Byk,Bzk,C_Lambda,C_mu,BfieldWN,spec_rank,res, &
                          &  deltx,T1Bfield,T2Bfield,T3Bfield)
  case(6)
    call MullerCarratiInduc(Uk,Vk,Wk,Bxk,Byk,Bzk,BfieldWN,spec_rank,res, &
                          &  deltx,T23,T12,T13)
    call ComputeDivTijVec1Vec2(T1Bfield,T2Bfield,T3Bfield,T12,T13,T23,BfieldWN,res)

  end select


!   IF ((.NOT. copyStructOnly(Bx,J12)) .OR. &
!      &(.NOT. copyStructOnly(Bx,J13)) .OR. &
!      &(.NOT. copyStructOnly(Bx,J23)) .OR. &
!      &(.NOT. copyStructOnly(Bx,EmTBbmodel))    ) RETURN 
! 
!    nbcpus = getnbcpus()
!      call computeJijTensor(BxK,ByK,BzK,BfieldWN,J12,J13,J23, res)
! !!!!! facteur 4, car on a un deux dans l'expression et le tenseur est antisymétrique
!     EmTBbmodel%values = 4.0_WP*( J12%values*T12%values &
!                                &+J13%values*T13%values &
!                                &+J23%values*T23%values )
!     meanEmTBbmodel  = computeFieldAvg(EmTBbmodel, spec_rank,nbcpus)
!     if (spec_rank .EQ.0) then
!        if (iterationBmodel .EQ.0) then
! 
!        open(10,file='meanEmTBbmodelSmagtime.out',form='formatted')
!          write(10,*)  iterationBmodel,meanEmTBbmodel
!        close(10)
!        else
!        open(10,file='meanEmTBbmodelSmagtime.out',form='formatted',position='append')
!          write(10,*)  iterationBmodel,meanEmTBbmodel
!        close(10)
!        endif
!     endif
!computeFieldPDF(iterationBmodel,EmTBbmode,100,spec_rank,n_iter)
!     iterationBmodel = iterationBmodel + 1
! 
!     CALL deleteDataLayout(J12)
!     CALL deleteDataLayout(J13)
!     CALL deleteDataLayout(J23)
!     CALL deleteDataLayout(EmTBbmodel)

  CALL deleteDataLayout(T12)
  CALL deleteDataLayout(T13)
  CALL deleteDataLayout(T23)

  nlx%values = nlx%values + T1Bfield%values
  nly%values = nly%values + T2Bfield%values
  nlz%values = nlz%values + T3Bfield%values

  call deleteDataLayout(T1Bfield)
  call deleteDataLayout(T2Bfield)
  call deleteDataLayout(T3Bfield)

 end subroutine computeBfieldModel

!---------------------------------------------------------------------------
!> @details 
!> This subroutine close the subgridmodels context 
!! @author = Antoine Vollant, LEGI
!---------------------------------------------------------------------------
 subroutine delete_LES()

  if ( associated(numSca) ) deallocate(numSca)
  if ( associated(numScaFilter) ) deallocate(numScaFilter)
  if ( associated(coef_sca) ) deallocate(coef_sca)
  if ( associated(type_avg_sca) ) deallocate(type_avg_sca)
  if ( associated(nameEOFile) ) deallocate(nameEOFile)
  call del_opt_estim_models()

 end subroutine delete_LES


end module subgridmodels
!> @}
