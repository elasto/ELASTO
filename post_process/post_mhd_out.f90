!> @addtogroup post_process
!! @{

!------------------------------------------------------------------------------

!
! MODULE: post_mhd
!
!
! DESCRIPTION: 
!>  This module provide a library of post-process routines, specifics to MHD simulations. This library
!! provide all that you need to build your own post-process.
!
!> @details 
!! 
!> @author
!! Mouloud KESSAR, LEGI
!
!------------------------------------------------------------------------------

module post_mhd_out

    use data
    use datalayout
    use toolbox 
    use stat_tools 
    use differential_tools 
    USE wavenumber_tools 
    use transforms_tools
    use subgrid_tools
    use differential_tools 
    use filtering_tools 
    use mhd_model
    use velmodels
    use post_mhd_ener_lib
    implicit none
    private

     ! ==== Public procedures ==== 
!     public          :: computealltheterms_out
    public          :: computeallthetermsMHDout
    ! ======= Mean terms ========
    real(WP), PROTECTED, SAVE :: meanadvecQDMsumT, meandiffvisQDMsumT, meandissvisQDMsumT, meanlorrQDMsumT, sommesumT ,sommesumTmhd
    real(WP), PROTECTED, SAVE :: meanadvecMHDsumT, meandiffvisMHDsumT, meandissvisMHDsumT, meanlorrMHDsumT, meanforMHDsumT 
    real(WP), PROTECTED, SAVE :: meanTotEsumT,meanEksumT,meanEmsumT,meanTotE,meanEk,meanEm
    real(WP), PROTECTED, SAVE :: eta,mu,Prm,dt,timesimu

    TYPE(COMPLEX_DATA_LAYOUT), PROTECTED, SAVE :: Valout3k,Valout2k,Valout1k,ValInk
!     real(WP), PROTECTED, ALLOCATABLE, dimension(:,:), SAVE :: meanDissGSEQDM, meanAdvecGSEQDMsumT, &
!                                                             & meanDifGSEQDM, meanTransEQDMGSsumT, meanTransEqdmSGSsumT
!     real(WP), PROTECTED, ALLOCATABLE, dimension(:,:), SAVE :: meanDiffUSGSEqdm, meanDissUSGSEqdmsumT, &
!                                                             & meanDiffBSGSEqdm, meanDissBSGSEqdmsumT, meanDissSGSEQDMsumT
    real(WP), PROTECTED, ALLOCATABLE, dimension(:,:), SAVE :: meanDifGSEMHD, meanAdvecGSEMHD, &
                                                            & meanDissGSEMHD, meanDiffUBSGSEMHD,&
                                                            & meanProdEMHD, meanDissUBSGSEMHD, meanDissSGSEMHD,meanTransGSEmEk
    real(WP), PROTECTED, ALLOCATABLE, dimension(:,:), SAVE :: meanEkGS, meanEmGS

    real(WP), PROTECTED, ALLOCATABLE, dimension(:,:), SAVE ::meanEmTBbmodelSmag,meanEqdmTUuumodelSmag, meanEqdmTUubmodelSmag
    real(WP), PROTECTED, ALLOCATABLE, dimension(:,:), SAVE ::meanEmTBbmodelSimil,meanEqdmTUuumodelSimil, meanEqdmTUubmodelSimil
    real(WP), PROTECTED, ALLOCATABLE, dimension(:,:), SAVE ::meanEmTBbmodelGrad,meanEqdmTUuumodelGrad, meanEqdmTUubmodelGrad

!     real(WP), PROTECTED, ALLOCATABLE, dimension(:,:), SAVE ::meanEmTBbmodelSmag,meanEqdmTUuumodelSmag, meanEqdmTUubmodelSmag
!     real(WP), PROTECTED, ALLOCATABLE, dimension(:,:), SAVE ::meanEmTBbmodelSimil,meanEqdmTUuumodelSimil, meanEqdmTUubmodelSimil
    real(WP), PROTECTED, ALLOCATABLE, dimension(:,:), SAVE :: meanbdivTubexact, meanbdivTubSmag
    
    INTEGER, PROTECTED, SAVE,PRIVATE :: mhdndeltaout,filterType
 contains

! ===== post mhd initialisation =====
!------------------------------------------------------------------------------
!> Compute and save spectrum of all field (scalars and velocity).
!! @author Mouloud KESSAR
!!    @param[in]    muIn       = viscosity
!!    @param[in]    PrmIn      = magnetic Prandl number
!!    @param[in]    UkIn       = spectral U field, for initialisation of post process
!! @details
!!        This function initialise the post_mhd module 
!------------------------------------------------------------------------------

! subroutine post_mhd_init_out(muIn,PrmIn,Ukin)
!    implicit none
!    REAL(WP), intent(in)                     :: muIn,PrmIn
!    TYPE(COMPLEX_DATA_LAYOUT), intent(in)    :: Ukin
! 
! 
!     
! end subroutine post_mhd_init_out

subroutine post_MHD_init(Uin,Vin,Win,Bin,&
                         & nbcpus,spec_rank)

    implicit none
    TYPE(REAL_DATA_LAYOUT), dimension(:), INTENT(IN) :: Bin(:)

    TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Uin,Vin,Win
    INTEGER :: nbcpus,spec_rank,i,ierr
    logical :: res
   IF(.NOT.wavesNumber_init(VelWN,spec_rank,Uin%Lx,Uin%Ly,Uin%Lz,Uin%nx,Uin%ny,Uin%nz)) RETURN

    IF (.NOT. initDataLayout("U_Spectral", &
            & Uk,(Uin%nx/2)+1,Uin%ny,Uin%nz,Uin%Lx,Uin%Ly,Uin%Lz,nbcpus,spec_rank,alongZ)) THEN
        WRITE(6,'(a,i0,a)')'[ERROR] initSolver on process ',spec_rank,&
            & ': not enought memory for U velocities in fourier space.'
        RETURN
    ENDIF
    IF (.NOT. initDataLayout("V_Spectral", &
            & Vk,(Vin%nx/2)+1,Vin%ny,Vin%nz,Vin%Lx,Vin%Ly,Vin%Lz,nbcpus,spec_rank,alongZ)) THEN
        WRITE(6,'(a,i0,a)')'[ERROR] initSolver on process ',spec_rank,&
            & ': not enought memory for V velocities in fourier space.'
        RETURN
    ENDIF
    IF (.NOT. initDataLayout("W_Spectral", &
            & Wk,(Win%nx/2)+1,Win%ny,Win%nz,Win%Lx,Win%Ly,Win%Lz,nbcpus,spec_rank,alongZ)) THEN
        WRITE(6,'(a,i0,a)')'[ERROR] initSolver on process ',spec_rank,&
            & ': not enought memory for W velocities in fourier space.'
        RETURN
    ENDIF
   CALL ftran(Uin,Uk,res)
   IF (.NOT.res) RETURN
   CALL ftran(Vin,Vk,res)
   IF (.NOT.res) RETURN
   CALL ftran(Win,Wk,res)
   IF (.NOT.res) RETURN

       allocate(Bk(3),stat=ierr)
       DO i=1,3
          IF (.NOT. initDataLayout("Bfield_Spectral", &
             & Bk(i),(Bin(1)%nx/2)+1,Bin(1)%ny,Bin(1)%nz,Bin(1)%Lx,Bin(1)%Ly,Bin(1)%Lz,nbcpus,spec_rank,alongZ)) THEN
             WRITE(6,'(a,i0,a)')'[ERROR] initSolver on process ',spec_rank,&
                               & ': not enought memory for Bfields in fourier space.'
            RETURN
          ENDIF
       ENDDO

   CALL ftran(Bin(1),Bk(1),res)
   IF (.NOT.res) RETURN
   CALL ftran(Bin(2),Bk(2),res)
   IF (.NOT.res) RETURN
   CALL ftran(Bin(3),Bk(3),res)
   IF (.NOT.res) RETURN

end subroutine post_MHD_init

subroutine computeallthetermsMHDout(Uin,Vin,Win,Bin,&
                         & nbcpus,spec_rank)

    implicit none

    TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Uin,Vin,Win
    INTEGER :: nbcpus,spec_rank
    TYPE(REAL_DATA_LAYOUT), dimension(:), INTENT(IN) :: Bin(:)
! 
!   CALL parser_read('MHD filter size',hdndelta)
!   CALL parser_read('MHD filter type',filterType)
!   CALL parser_read('MHD Ener Post',enerpostHD)
!   CALL parser_read('HD Hel Post',HelpostHD)
  CALL parser_read('Viscosity',mu)
  CALL parser_read('Prandl Number',Prm)

!   muout=muIn
!   Prm=PrmIn
  eta=mu/Prm
  CALL parser_read('MHD filter size',mhdndeltaout)
  CALL parser_read('MHD filter type',filterType)

  allocate( meanProdEMHD(2,mhdndeltaout)      )
  allocate( meanDifGSEMHD(2,mhdndeltaout)     )
  allocate( meanDissGSEMHD(2,mhdndeltaout)    )
  allocate( meanDissSGSEMHD(2,mhdndeltaout)   )
  allocate( meanAdvecGSEMHD(2,mhdndeltaout)   )
  allocate( meanDiffUBSGSEMHD(2,mhdndeltaout) )
  allocate( meanDissUBSGSEMHD(2,mhdndeltaout) )
  allocate( meanTransGSEmEk(2,mhdndeltaout) )

  allocate( meanEkGS(2,mhdndeltaout) )
  allocate( meanEmGS(2,mhdndeltaout) )

  allocate(meanEmTBbmodelSmag(2,mhdndeltaout))
  allocate(meanEqdmTUuumodelSmag(2,mhdndeltaout))
  allocate(meanEqdmTUubmodelSmag(2,mhdndeltaout))

  allocate(meanEmTBbmodelSimil(2,mhdndeltaout))
  allocate(meanEqdmTUuumodelSimil(2,mhdndeltaout))
  allocate(meanEqdmTUubmodelSimil(2,mhdndeltaout))

  allocate(meanEmTBbmodelGrad(2,mhdndeltaout))
  allocate(meanEqdmTUuumodelGrad(2,mhdndeltaout))
  allocate(meanEqdmTUubmodelGrad(2,mhdndeltaout))

  allocate(meanbdivTubexact(2,mhdndeltaout))
  allocate(meanbdivTubSmag(2,mhdndeltaout))

  call post_MHD_init(Uin,Vin,Win,Bin,&
                         & nbcpus,spec_rank)

  call computealltheterms_out(Uin,Vin,Win,Uk,Vk,Wk,Bin,Bk,&
                         & VelWN,BfieldWN,nbcpus,spec_rank)

!  if (spec_rank.EQ.0) print*, "mu", mu

end subroutine computeallthetermsMHDout

! ===== Compute all the terms =====

!------------------------------------------------------------------------------
!> Compute and save mean values of all quantities in Ennergie équations .
!! @author Mouloud KESSAR
!!    @param[in]    spec_rank     = mpi rank inside spectral communicator
!!    @param[in]    Uin,Vin,Win   = Velocity fields ins physical space
!!    @param[in]    Bin           = magnetic Field
!!    @param[in]    WaveVel       = velocity wavenumber
!!    @param[in]    WaveB         = Bfield wavenumber
!!    @param[in]    nbcpus        = numbers of cpus
!!    @param[in]    simtime       = time in the simulation
!!    @param[in]    tstep         = step time in the simulation
!!    
!! @details
!!        This function compute the mean values of all quantities in Ennergie équations
!!    stores it, and print this vlaues at each iteration
!------------------------------------------------------------------------------
subroutine computealltheterms_out(Uin,Vin,Win,Ukin,Vkin,Wkin,Bin,Bkin,&
                         & WaveVel,WaveB,nbcpus,spec_rank)


    implicit none

    TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Uin,Vin,Win
    TYPE(REAL_DATA_LAYOUT), dimension(:), INTENT(IN) :: Bin(:)
    TYPE(WaveNumbers),INTENT(IN) :: WaveVel,WaveB
    TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Ukin,Vkin,Wkin
    TYPE(COMPLEX_DATA_LAYOUT),dimension(:), INTENT(IN) :: Bkin(:)

    INTEGER :: i,nbcpus,spec_rank
!     REAL(WP) :: tstep,simtime
!     dt=tstep
!     timesimu=simtime
 

    call computeEQDM_EQMHD_out(Uin,Vin,Win,Bin,Ukin,Vkin,Wkin,Bkin,WaveVel,nbcpus,spec_rank)
    call computeallGSandSGS_out(Uin,Vin,Win,Bin,Ukin,Vkin,Wkin,Bkin,WaveVel,nbcpus,spec_rank)

    if(spec_rank.EQ.0) then

!      if (timesimu .LT. 0.00000001)then

       open(10,file='mean_Ek_GS.out',form='formatted')
       do i=1,mhdndeltaout
       write(10,*) timesimu, meanEkGS(1,i),meanEkGS(2,i)
       enddo
       close(10)

       open(10,file='mean_Em_GS.out',form='formatted')
        do i=1,mhdndeltaout
           write(10,*) timesimu, meanEmGS(1,i),meanEmGS(2,i)
        enddo
       close(10)

       open(10,file='meanDissUBSGSEMHDsumT.out',form='formatted')
        do i=1,mhdndeltaout
         write(10,*)  meanDissUBSGSEMHD(1,i), meanDissUBSGSEMHD(2,i)
        enddo
       close(10)

       open(10,file='meanDiffUBSGSEMHDsumT.out',form='formatted')
        do i=1,mhdndeltaout
         write(10,*)  meanDiffUBSGSEMHD(1,i), meanDiffUBSGSEMHD(2,i)
        enddo
       close(10)
      
       open(10,file='meanAdvecGSEMHDsumT.out',form='formatted')
        do i=1,mhdndeltaout
         write(10,*)  meanAdvecGSEMHD(1,i), meanAdvecGSEMHD(2,i)
        enddo
       close(10)

       open(10,file='meanDifGSEMHD.out',form='formatted')
        do i=1,mhdndeltaout
         write(10,*)  meanDifGSEMHD(1,i), meanDifGSEMHD(2,i)
        enddo
       close(10)

       open(10,file='meanProdEMHD.out',form='formatted')
        do i=1,mhdndeltaout
         write(10,*)  meanProdEMHD(1,i), meanProdEMHD(2,i)
        enddo
       close(10)

       open(10,file='meanDissGSEMHD.out',form='formatted')
        do i=1,mhdndeltaout
         write(10,*)  meanDissGSEMHD(1,i), meanDissGSEMHD(2,i)
        enddo
       close(10)

   endif
    call aprioriallMHDEm(Uin,Vin,Win,Bin(1),Bin(2),Bin(3),Ukin,Vkin,Wkin,&
                    & Bkin(1),Bkin(2),Bkin(3),WaveB,spec_rank,nbcpus)


end subroutine computealltheterms_out




subroutine computeEQDM_EQMHD_out(Uin,Vin,Win,Bin,Ukin,Vkin,Wkin,Bkin,Wave,nbcpus,spec_rank)

  implicit none
  TYPE(REAL_DATA_LAYOUT), INTENT(IN)          :: Uin,Vin,Win
  TYPE(REAL_DATA_LAYOUT), dimension(:), INTENT(IN) :: Bin(:)
  TYPE(COMPLEX_DATA_LAYOUT) :: Ukin,Vkin,Wkin
  TYPE(COMPLEX_DATA_LAYOUT), dimension(:), INTENT(IN) :: Bkin(:)
  REAL(WP) :: etain

  INTEGER :: nbcpus, spec_rank, ndelta
  TYPE(WAVENUMBERS) :: Wave
  REAL(WP) :: meanadvecMHD, meandiffvisMHD, meandissvisMHD, meanlorrMHD, meanforMHD ,somme
  logical :: res
  
  call eqmhd_lib(Uin,Vin,Win,Bin(1),Bin(2),Bin(3),Ukin,Vkin,Wkin,&
               & Bkin(1),Bkin(2),Bkin(3),Wave,nbcpus,spec_rank,&
               & meanadvecMHD, meandiffvisMHD, meandissvisMHD, meanlorrMHD, meanforMHD ,somme,etain,meanEm)

    if(spec_rank.EQ.0) then


!        open(10,file='mean_values_QDM.out',form='formatted')
!        write(10,*) timesimu, meanadvecQDMsumT, meandiffvisQDMsumT, meandissvisQDMsumT, meanlorrQDMsumT,&
!                  & sommesumT,meanEk,meanEm
!        close(10)

       open(10,file='mean_values_MHD.out',form='formatted')
       write(10,*) meanadvecMHD, meandiffvisMHD, meandissvisMHD, meanlorrMHD, meanforMHD ,somme
       close(10)
    endif

end subroutine computeEQDM_EQMHD_out


!======= compute all the grid and subgrid terms =========
!------------------------------------------------------------------------------
!> Compute and save mean values of all quantities in GS and SGS énergie equations
!! @author Mouloud KESSAR
!!    @param[in]    spec_rank     = mpi rank inside spectral communicator
!!    @param[in]    Uin,Vin,Win   = Velocity fields ins physical space
!!    @param[in]    Bin           = magnetic Field
!!    @param[in]    Wave          = wavenumber
!!    @param[in]    nbcpus        = numbers of cpus
!! @details
!!    
!------------------------------------------------------------------------------


subroutine computeallGSandSGS_out(Uin,Vin,Win,Bin,Ukin,Vkin,Wkin,Bkin,Wave,nbcpus,spec_rank)

  implicit none
  TYPE(REAL_DATA_LAYOUT), INTENT(IN)          :: Uin,Vin,Win
  TYPE(REAL_DATA_LAYOUT), dimension(:), INTENT(IN) :: Bin(:)
  TYPE(REAL_DATA_LAYOUT) :: Bfx,Bfy,Bfz,Uf,Vf,Wf
  TYPE(COMPLEX_DATA_LAYOUT) :: Bfxk,Bfyk,Bfzk,Ukf,Vkf,Wkf
  TYPE(COMPLEX_DATA_LAYOUT) :: Ukin,Vkin,Wkin
  TYPE(COMPLEX_DATA_LAYOUT), dimension(:), INTENT(IN) :: Bkin(:)


  INTEGER :: nbcpus, spec_rank, ndelta
  TYPE(WAVENUMBERS) :: Wave
  REAL(WP) :: delta,sommeGS
  logical :: res
  
  IF ((.NOT. copyStructOnly(Uin,Uf)) .OR. &
     &(.NOT. copyStructOnly(Uin,Vf)) .OR. &
     &(.NOT. copyStructOnly(Uin,Wf)) .OR. &
     &(.NOT. copyStructOnly(Bin(1),Bfx)) .OR. &
     &(.NOT. copyStructOnly(Bin(1),Bfy)) .OR. &
     &(.NOT. copyStructOnly(Bin(1),Bfz)) ) RETURN

  IF ((.NOT. copyStructOnly(Ukin,Ukf)) .OR. &
     &(.NOT. copyStructOnly(Ukin,Vkf)) .OR. &
     &(.NOT. copyStructOnly(Ukin,Wkf)) .OR. &
     &(.NOT. copyStructOnly(Bkin(1),Bfxk)) .OR. &
     &(.NOT. copyStructOnly(Bkin(1),Bfyk)) .OR. &
     &(.NOT. copyStructOnly(Bkin(1),Bfzk)) ) RETURN

    DO ndelta=1,mhdndeltaout
    delta = Uin%Lx*(2*ndelta)/Uin%nx

!     meanTransEQDMGS(1,ndelta)  = 2*ndelta
!     meanTransEqdmSGS(1,ndelta) = 2*ndelta
!     meanDifGSEQDM(1,ndelta)    = 2*ndelta
!     meanAdvecGSEQDM(1,ndelta)  = 2*ndelta
!     meanDissGSEQDMsumT(1,ndelta)   = 2*ndelta
!     meanDiffUSGSEqdmsumT(1,ndelta) = 2*ndelta
!     meanDiffBSGSEqdmsumT(1,ndelta) = 2*ndelta
!     meanDissUSGSEqdmsumT(1,ndelta) = 2*ndelta
!     meanDissBSGSEqdmsumT(1,ndelta) = 2*ndelta
!     meanDissSGSEQDMsumT(1,ndelta)  = 2*ndelta

    meanDifGSEMHD(1,ndelta)     = 2*ndelta
    meanAdvecGSEMHD(1,ndelta)   = 2*ndelta
    meanDissGSEMHD(1,ndelta)    = 2*ndelta
    meanDissSGSEMHD(1,ndelta)   = 2*ndelta
    meanDiffUBSGSEMHD(1,ndelta) = 2*ndelta
    meanDissUBSGSEMHD(1,ndelta) = 2*ndelta
    meanProdEMHD(1,ndelta)      = 2*ndelta
    meanTransGSEmEk(1,ndelta)   = 2*ndelta
!     meanEkGSsumT(1,ndelta) = 2*ndelta
    meanEmGS(1,ndelta) = 2*ndelta

   call computeFilter(Wave, delta, Ukin, Ukf, filterType)
   CALL btran(Ukf,Uf,res)
   IF (.NOT.res) RETURN

   call computeFilter(Wave, delta, Vkin, Vkf, filterType)
   CALL btran(Vkf,Vf,res)
   IF (.NOT.res) RETURN

   call computeFilter(Wave, delta, Wkin, Wkf, filterType)
   CALL btran(Wkf,Wf,res)
   IF (.NOT.res) RETURN


   call computeFilter(Wave, delta, Bin(1), Bfxk, filterType)
   CALL btran(Bfxk,Bfx,res)
   IF (.NOT.res) RETURN

   call computeFilter(Wave, delta, Bin(2), Bfyk, filterType)
   CALL btran(Bfyk,Bfy,res)
   IF (.NOT.res) RETURN

   call computeFilter(Wave, delta, Bin(3), Bfzk, filterType)
   CALL btran(Bfzk,Bfz,res)
   IF (.NOT.res) RETURN


    CALL EmhdGS_lib( Uf,Vf,Wf,Bfx,Bfy,Bfz,Ukin,Vkin,Wkin,Uin,Vin,Win,Bin(1),Bin(2),Bin(3),Bfxk,Bfyk,Bfzk,&
                   & Wave,delta,ndelta,filterType,nbcpus,spec_rank,eta,&
                   & meanDiffUBSGSEMHD(2,ndelta),meanDissUBSGSEMHD(2,ndelta),meanadvecGSEMHD(2,ndelta),&
                   & meanDifGSEMHD(2,ndelta), meanDissGSEMHD(2,ndelta), meanTransGSEmEk(2,ndelta), meanProdEMHD(2,ndelta) ,sommeGS,meanEmGS(2,ndelta))
    ENDDO


    CALL deleteDataLayout(Uf)
    CALL deleteDataLayout(Vf)
    CALL deleteDataLayout(Wf)
    CALL deleteDataLayout(Bfx)
    CALL deleteDataLayout(Bfy)
    CALL deleteDataLayout(Bfz)

    CALL deleteDataLayout(Ukf)
    CALL deleteDataLayout(Vkf)
    CALL deleteDataLayout(Wkf)
    CALL deleteDataLayout(Bfxk)
    CALL deleteDataLayout(Bfyk)
    CALL deleteDataLayout(Bfzk)

end subroutine computeallGSandSGS_out


subroutine aprioriallMHDEm(Uin,Vin,Win,Bxin,Byin,Bzin,Ukin,Vkin,Wkin,&
                    & Bxkin,Bykin,Bzkin,Wave,spec_rank,nbcpus)

  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Uin,Vin,Win
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Bxin,Byin,Bzin 
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Ukin,Vkin,Wkin
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Bxkin,Bykin,Bzkin
  TYPE(WaveNumbers), INTENT(IN) :: Wave

  TYPE(REAL_DATA_LAYOUT) :: Uf,Vf,Wf
  TYPE(REAL_DATA_LAYOUT) :: Bxf,Byf,Bzf 
  TYPE(COMPLEX_DATA_LAYOUT) :: Ukf,Vkf,Wkf
  TYPE(COMPLEX_DATA_LAYOUT) :: Bxkf,Bykf,Bzkf

  INTEGER :: spec_rank,nbcpus,i
  logical :: res
  REAL(WP) :: delta
  Integer :: ndelta

  IF ((.NOT. copyStructOnly(Bxkin,Ukf)) .OR. &
     &(.NOT. copyStructOnly(Bxkin,Vkf)) .OR. &
     &(.NOT. copyStructOnly(Bxkin,Wkf)) .OR. &
     &(.NOT. copyStructOnly(Bxkin,Bxkf)) .OR. &
     &(.NOT. copyStructOnly(Bxkin,Bykf)) .OR. &
     &(.NOT. copyStructOnly(Bxkin,Bzkf)) .OR. &
     &(.NOT. copyStructOnly(Bxin,Bxf)) .OR. &
     &(.NOT. copyStructOnly(Bxin,Byf)) .OR. &
     &(.NOT. copyStructOnly(Bxin,Bzf)) .OR. &
     &(.NOT. copyStructOnly(Bxin,Uf)) .OR. &
     &(.NOT. copyStructOnly(Bxin,Vf)) .OR. &
     &(.NOT. copyStructOnly(Bxin,Wf))    ) RETURN 

    DO ndelta=1,mhdndeltaout
    delta = Uin%Lx*(2*ndelta)/Uin%nx

    meanEmTBbmodelSmag(1,ndelta)    = 2*ndelta
    meanEqdmTUuumodelSmag(1,ndelta) = 2*ndelta
    meanEqdmTUubmodelSmag(1,ndelta) = 2*ndelta

    meanEmTBbmodelGrad(1,ndelta)    = 2*ndelta
    meanEqdmTUuumodelGrad(1,ndelta) = 2*ndelta
    meanEqdmTUubmodelGrad(1,ndelta) = 2*ndelta

    meanEmTBbmodelSimil(1,ndelta)    = 2*ndelta
    meanEqdmTUuumodelSimil(1,ndelta) = 2*ndelta
    meanEqdmTUubmodelSimil(1,ndelta) = 2*ndelta

    meanbdivTubexact(1,ndelta) = 2*ndelta
    meanbdivTubSmag(1,ndelta) = 2*ndelta
    
     CALL computeFilter(wave,delta,Bxkin,Bxkf,filterType) 
     CALL btran(Bxkf,Bxf,res)
     CALL computeFilter(wave,delta,Bykin,Bykf,filterType) 
     CALL btran(Bykf,Byf,res)
     CALL computeFilter(wave,delta,Bzkin,Bzkf,filterType) 
     CALL btran(Bzkf,Bzf,res)

     CALL computeFilter(wave,delta,Ukin,Ukf,filterType) 
     CALL btran(Ukf,Uf,res)
     CALL computeFilter(wave,delta,Vkin,Vkf,filterType) 
     CALL btran(Vkf,Vf,res)
     CALL computeFilter(wave,delta,Wkin,Wkf,filterType) 
     CALL btran(Wkf,Wf,res)

     call  aprioriMHDEm(Uf,Vf,Wf,Bxf,Byf,Bzf,Ukf,Vkf,Wkf,&
               & Bxkf,Bykf,Bzkf,Wave,spec_rank,delta,filterType,meanEmTBbmodelSmag(2,ndelta),nbcpus)

     call  aprioriEQDmMHD(Uf,Vf,Wf,Bxf,Byf,Bzf,Ukf,Vkf,Wkf,&
                    & Bxkf,Bykf,Bzkf,Wave,spec_rank,filterType,delta,meanEqdmTUuumodelSmag(2,ndelta),&
                    & meanEqdmTUubmodelSmag(2,ndelta),meanEqdmTUuumodelSimil(2,ndelta),meanEqdmTUubmodelSimil(2,ndelta),&
                    & meanEqdmTUuumodelGrad(2,ndelta), meanEqdmTUubmodelGrad(2,ndelta),nbcpus)

     call aprioBDivtub(Uin,Vin,Win,Bxin,Byin,Bzin,Ukin,Vkin,Wkin,&
                     & Bxkin,Bykin,Bzkin,Uf,Vf,Wf,Bxf,Byf,Bzf,& 
                     & Wave,spec_rank,delta,filterType,nbcpus,meanbdivTubexact(2,ndelta), meanbdivTubSmag(2,ndelta))
 
    enddo

    if(spec_rank.EQ.0) then

       open(10,file='meanEmTBbmodelSmag.out',form='formatted')
        do i=1,mhdndeltaout
         write(10,*)  meanEmTBbmodelSmag(1,i), meanEmTBbmodelSmag(2,i)
        enddo
       close(10)

       open(10,file='meanEqdmTUubSmag.out',form='formatted')
        do i=1,mhdndeltaout
         write(10,*)  meanEqdmTUubmodelSmag(1,i), meanEqdmTUubmodelSmag(2,i)
        enddo
       close(10)

       open(10,file='meanEqdmTUuuSmag.out',form='formatted')
        do i=1,mhdndeltaout
         write(10,*)  meanEqdmTUuumodelSmag(1,i), meanEqdmTUuumodelSmag(2,i)
        enddo
       close(10)

       open(10,file=' meanbdivTubSmag.out',form='formatted')
        do i=1,mhdndeltaout
         write(10,*)   meanbdivTubSmag(1,i),  meanbdivTubSmag(2,i)
        enddo
       close(10)

       open(10,file='meanbdivTubexact.out',form='formatted')
        do i=1,mhdndeltaout
         write(10,*)  meanbdivTubexact(1,i), meanbdivTubexact(2,i)
        enddo
       close(10)
    endif

    CALL deleteDataLayout(Ukf)
    CALL deleteDataLayout(Vkf)
    CALL deleteDataLayout(Wkf)
    CALL deleteDataLayout(Uf)
    CALL deleteDataLayout(Vf)
    CALL deleteDataLayout(Wf)

    CALL deleteDataLayout(Bxkf)
    CALL deleteDataLayout(Bykf)
    CALL deleteDataLayout(Bzkf)
    CALL deleteDataLayout(Bxf)
    CALL deleteDataLayout(Byf)
    CALL deleteDataLayout(Bzf)


end subroutine aprioriallMHDEm

end module post_mhd_out
!> @}
