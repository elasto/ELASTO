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

module post_mhd

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
    use Eqdm_lorentz_models
    implicit none
    private

     ! ==== Public procedures ==== 
    public          :: computealltheterms
    public          :: post_mhd_init
    ! ======= Mean terms ========
    real(WP), PROTECTED, SAVE :: meanadvecQDMsumT, meandiffvisQDMsumT, meandissvisQDMsumT, meanlorrQDMsumT, sommesumT ,sommesumTmhd
    real(WP), PROTECTED, SAVE :: meanadvecMHDsumT, meandiffvisMHDsumT, meandissvisMHDsumT, meanlorrMHDsumT, meanforMHDsumT 
    real(WP), PROTECTED, SAVE :: meanTotEsumT,meanEksumT,meanEmsumT,meanTotE,meanEk,meanEm
    real(WP), PROTECTED, SAVE :: eta,mu,Prm,dt,timesimu
    TYPE(COMPLEX_DATA_LAYOUT), PROTECTED, SAVE :: Valout3k,Valout2k,Valout1k,ValInk
    real(WP), PROTECTED, ALLOCATABLE, dimension(:,:), SAVE :: meanDissGSEQDMsumT, meanAdvecGSEQDMsumT, &
                                                            & meanDifGSEQDMsumT, meanTransEQDMGSsumT, meanTransEqdmSGSsumT
    real(WP), PROTECTED, ALLOCATABLE, dimension(:,:), SAVE :: meanDiffUSGSEqdmsumT, meanDissUSGSEqdmsumT, &
                                                            & meanDiffBSGSEqdmsumT, meanDissBSGSEqdmsumT, meanDissSGSEQDMsumT
    real(WP), PROTECTED, ALLOCATABLE, dimension(:,:), SAVE :: meanDifGSEMHDsumT, meanAdvecGSEMHDsumT, &
                                                            & meanDissGSEMHDsumT, meanDiffUBSGSEMHDsumT,&
                                                            & meanProdEMHDsumT, meanDissUBSGSEMHDsumT, meanDissSGSEMHDsumT
    real(WP), PROTECTED, ALLOCATABLE, dimension(:,:), SAVE :: meanEkGSsumT, meanEmGSsumT

    real(WP), PROTECTED, ALLOCATABLE, dimension(:,:), SAVE ::meanEmTBbmodelSmag,meanEqdmTUuumodelSmag, meanEqdmTUubmodelSmag

    real(WP), PROTECTED, ALLOCATABLE, dimension(:,:), SAVE ::meanEmTBbmodelSimil,meanEqdmTUuumodelSimil, meanEqdmTUubmodelSimil
    real(WP), PROTECTED, ALLOCATABLE, dimension(:,:), SAVE ::meanEmTBbmodelGrad,meanEqdmTUuumodelGrad, meanEqdmTUubmodelGrad

    real(WP), PROTECTED, ALLOCATABLE, dimension(:,:), SAVE ::meanEqdmTUumodelHamba!,meanEmTBbmodelHamba
    real(WP), PROTECTED, ALLOCATABLE, dimension(:,:), SAVE ::meanEmhdTBbmodelCarrati,meanEqdmTUumodelCarrati


    real(WP),private,save :: C_mu=0.046_WP
    real(WP),private,save :: C_lambda=(5/7*0.046_WP)

    INTEGER, PROTECTED, SAVE,PRIVATE :: mhdndelta,filterType
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

subroutine post_mhd_init(muIn,PrmIn,Ukin)
   implicit none
   REAL(WP), intent(in)                     :: muIn,PrmIn
   TYPE(COMPLEX_DATA_LAYOUT), intent(in)    :: Ukin
   IF ((.NOT. copyStructOnly(UkIn,valOut1k)) .OR. &

      &(.NOT. copyStructOnly(UkIn,valOut2k)) .OR. &
      &(.NOT. copyStructOnly(UkIn,valOut3k)) .OR. &
      &(.NOT. copyStructOnly(UkIn,valInk)) ) RETURN

  mu=muIn
  Prm=PrmIn
  eta=mu/Prm
  CALL parser_read('MHD filter size',mhdndelta)
  CALL parser_read('MHD filter type',filterType)
  allocate( meanTransEQDMGSsumT(2,mhdndelta) )  
  allocate( meanTransEqdmSGSsumT(2,mhdndelta) ) 
  allocate( meanDifGSEQDMsumT(2,mhdndelta) )    
  allocate( meanAdvecGSEQDMsumT(2,mhdndelta) )  
  allocate( meanDissGSEQDMsumT(2,mhdndelta) )  
  allocate( meanDiffUSGSEqdmsumT(2,mhdndelta) )
  allocate( meanDiffBSGSEqdmsumT(2,mhdndelta) )
  allocate( meanDissUSGSEqdmsumT(2,mhdndelta) )
  allocate( meanDissBSGSEqdmsumT(2,mhdndelta) )
  allocate( meanDissSGSEQDMsumT(2,mhdndelta) )

  allocate( meanProdEMHDsumT(2,mhdndelta)      )
  allocate( meanDifGSEMHDsumT(2,mhdndelta)     )
  allocate( meanDissGSEMHDsumT(2,mhdndelta)    )
  allocate( meanDissSGSEMHDsumT(2,mhdndelta)   )
  allocate( meanAdvecGSEMHDsumT(2,mhdndelta)   )
  allocate( meanDiffUBSGSEMHDsumT(2,mhdndelta) )
  allocate( meanDissUBSGSEMHDsumT(2,mhdndelta) )

  allocate( meanEkGSsumT(2,mhdndelta) )
  allocate( meanEmGSsumT(2,mhdndelta) )

  allocate(meanEmTBbmodelSmag(2,mhdndelta))
  allocate(meanEqdmTUuumodelSmag(2,mhdndelta))
  allocate(meanEqdmTUubmodelSmag(2,mhdndelta))

  allocate(meanEmTBbmodelGrad(2,mhdndelta))
  allocate(meanEqdmTUuumodelGrad(2,mhdndelta))
  allocate(meanEqdmTUubmodelGrad(2,mhdndelta))

  allocate(meanEmTBbmodelSimil(2,mhdndelta))
  allocate(meanEqdmTUuumodelSimil(2,mhdndelta))
  allocate(meanEqdmTUubmodelSimil(2,mhdndelta))

  allocate(meanEqdmTUumodelCarrati(2,mhdndelta))
  allocate(meanEmhdTBbmodelCarrati(2,mhdndelta))

  allocate(meanEqdmTUumodelHamba(2,mhdndelta))
!   allocate(meanEmTBbmodelHamba(2,mhdndelta))

end subroutine    



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
subroutine computealltheterms(Uin,Vin,Win,Ukin,Vkin,Wkin,Bin,Bkin,&
                         & WaveVel,WaveB,nbcpus,spec_rank,simtime,tstep)


    implicit none

    TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Uin,Vin,Win
    TYPE(REAL_DATA_LAYOUT), dimension(:), INTENT(IN) :: Bin(:)
    TYPE(WaveNumbers),INTENT(IN) :: WaveVel,WaveB
    TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Ukin,Vkin,Wkin
    TYPE(COMPLEX_DATA_LAYOUT),dimension(:), INTENT(IN) :: Bkin(:)

    INTEGER :: i,nbcpus,spec_rank
    REAL(WP) :: tstep,simtime
    dt=tstep
    timesimu=simtime
 

    call eqdm(Uin,Vin,Win,Bin,WaveVel,nbcpus,spec_rank)
    call eqmhd(Uin,Vin,Win,Bin,WaveB,nbcpus,spec_rank)
    call computeallGSandSGS(Uin,Vin,Win,Bin,WaveVel,nbcpus,spec_rank)
    call aprioriallMHDEm(Uin,Vin,Win,Bin(1),Bin(2),Bin(3),Ukin,Vkin,Wkin,&
                    & Bkin(1),Bkin(2),Bkin(3),WaveB,spec_rank,nbcpus)
    if(spec_rank.EQ.0) then

     if (timesimu .LT. 0.00000001)then
       open(10,file='mean_values_QDM.out',form='formatted')
       write(10,*) timesimu, meanadvecQDMsumT, meandiffvisQDMsumT, meandissvisQDMsumT, meanlorrQDMsumT,&
                 & sommesumT,meanEk,meanEm
       close(10)

       open(10,file='mean_values_MHD.out',form='formatted')
       write(10,*) timesimu,meanadvecMHDsumT, meandiffvisMHDsumT, meandissvisMHDsumT, meanlorrMHDsumT,&
                 & meanforMHDsumT  ,sommesumTmhd
       close(10)

       open(10,file='mean_Ek_GS.out',form='formatted')
       do i=1,mhdndelta
       write(10,*) timesimu, meanEkGSsumT(1,i),meanEkGSsumT(2,i)
       enddo
       close(10)

       open(10,file='mean_Em_GS.out',form='formatted')
        do i=1,mhdndelta
           write(10,*) timesimu, meanEmGSsumT(1,i),meanEmGSsumT(2,i)
        enddo
       close(10)

       open(10,file='meanDissGSEQDM.out',form='formatted')
        do i=1,mhdndelta
          write(10,*) timesimu, meanDissGSEQDMsumT(1,i),meanDissGSEQDMsumT(2,i)
        enddo
       close(10)

       open(10,file='meanAdvecGSEQDM.out',form='formatted')
        do i=1,mhdndelta
         write(10,*) timesimu, meanAdvecGSEQDMsumT(1,i), meanAdvecGSEQDMsumT(2,i)
        enddo
       close(10)

       open(10,file='meanDifGSEQDM.out',form='formatted')
        do i=1,mhdndelta
         write(10,*) timesimu, meanDifGSEQDMsumT(1,i), meanDifGSEQDMsumT(2,i)
        enddo
       close(10)     
  
       open(10,file='meanTransEQDMGS.out',form='formatted')
        do i=1,mhdndelta
          write(10,*) timesimu, meanTransEQDMGSsumT(1,i), meanTransEQDMGSsumT(2,i)
        enddo
       close(10)

       open(10,file='meanTransEqdmSGS.out',form='formatted')
        do i=1,mhdndelta
         write(10,*) timesimu, meanTransEqdmSGSsumT(1,i), meanTransEqdmSGSsumT(2,i)
        enddo
       close(10)

       open(10,file='meanDiffUSGSEqdm.out',form='formatted')
        do i=1,mhdndelta                                       
         write(10,*) timesimu, meanDiffUSGSEqdmsumT(1,i), meanDiffUSGSEqdmsumT(2,i)
        enddo                                                
       close(10)

       open(10,file='meanDiffBSGSEqdm.out',form='formatted')
        do i=1,mhdndelta
         write(10,*) timesimu, meanDiffBSGSEqdmsumT(1,i), meanDiffBSGSEqdmsumT(2,i)
        enddo
       close(10)

       open(10,file='meanDissBSGSEqdm.out',form='formatted')
        do i=1,mhdndelta
         write(10,*) timesimu, meanDissBSGSEqdmsumT(1,i), meanDissBSGSEqdmsumT(2,i)
        enddo
       close(10)
       open(10,file='meanDissUSGSEqdm.out',form='formatted')
        do i=1,mhdndelta
         write(10,*) timesimu, meanDissUSGSEqdmsumT(1,i), meanDissUSGSEqdmsumT(2,i)
        enddo
       close(10)

       open(10,file='meanDissUBSGSEMHDsumT.out',form='formatted')
        do i=1,mhdndelta
         write(10,*)  timesimu, meanDissUBSGSEMHDsumT(1,i), meanDissUBSGSEMHDsumT(2,i)
        enddo
       close(10)

       open(10,file='meanDiffUBSGSEMHDsumT.out',form='formatted')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanDiffUBSGSEMHDsumT(1,i), meanDiffUBSGSEMHDsumT(2,i)
        enddo
       close(10)
      
       open(10,file='meanAdvecGSEMHDsumT.out',form='formatted')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanAdvecGSEMHDsumT(1,i), meanAdvecGSEMHDsumT(2,i)
        enddo
       close(10)

       open(10,file='meanDifGSEMHDsumT.out',form='formatted')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanDifGSEMHDsumT(1,i), meanDifGSEMHDsumT(2,i)
        enddo
       close(10)

       open(10,file='meanProdEMHDsumT.out',form='formatted')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanProdEMHDsumT(1,i), meanProdEMHDsumT(2,i)
        enddo
       close(10)

       open(10,file='meanDissGSEMHDsumT.out',form='formatted')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanDissGSEMHDsumT(1,i), meanDissGSEMHDsumT(2,i)
        enddo
       close(10)

       open(10,file='meanDissSGSEQDMsumT.out',form='formatted')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanDissSGSEQDMsumT(1,i), meanDissSGSEQDMsumT(2,i)
        enddo
       close(10)

       open(10,file='meanDissSGSEMHDsumT.out',form='formatted')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanDissSGSEMHDsumT(1,i), meanDissSGSEMHDsumT(2,i)
        enddo
       close(10)

       open(10,file='meanEmTBbmodelSmag.out',form='formatted')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanEmTBbmodelSmag(1,i), meanEmTBbmodelSmag(2,i)
        enddo
       close(10)

       open(10,file='meanEqdmTUubSmag.out',form='formatted')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanEqdmTUubmodelSmag(1,i), meanEqdmTUubmodelSmag(2,i)
        enddo
       close(10)

       open(10,file='meanEqdmTUuuSmag.out',form='formatted')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanEqdmTUuumodelSmag(1,i), meanEqdmTUuumodelSmag(2,i)
        enddo
       close(10)

       open(10,file='meanEmTBbmodelGrad.out',form='formatted')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanEmTBbmodelGrad(1,i), meanEmTBbmodelGrad(2,i)
        enddo
       close(10)

       open(10,file='meanEqdmTUubGrad.out',form='formatted')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanEqdmTUubmodelGrad(1,i), meanEqdmTUubmodelGrad(2,i)
        enddo
       close(10)

       open(10,file='meanEqdmTUuuGrad.out',form='formatted')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanEqdmTUuumodelGrad(1,i), meanEqdmTUuumodelGrad(2,i)
        enddo
       close(10)

       open(10,file='meanEmTBbmodelSimil.out',form='formatted')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanEmTBbmodelSimil(1,i), meanEmTBbmodelSimil(2,i)
        enddo
       close(10)

       open(10,file='meanEqdmTUubSimil.out',form='formatted')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanEqdmTUubmodelSimil(1,i), meanEqdmTUubmodelSimil(2,i)
        enddo
       close(10)

       open(10,file='meanEqdmTUuuSimil.out',form='formatted')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanEqdmTUuumodelSimil(1,i), meanEqdmTUuumodelSimil(2,i)
        enddo
       close(10)

       open(10,file='meanEqdmTUuHamba.out',form='formatted')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanEqdmTUumodelHamba(1,i), meanEqdmTUumodelHamba(2,i)
        enddo
       close(10)

       open(10,file='meanEqdmTUuCarrati.out',form='formatted')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanEqdmTUumodelCarrati(1,i), meanEqdmTUumodelCarrati(2,i)
        enddo
       close(10)


       open(10,file='meanEmhdTBbCarrati.out',form='formatted')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanEmhdTBbmodelCarrati(1,i), meanEmhdTBbmodelCarrati(2,i)
        enddo
       close(10)

     else

      open(10,file='mean_values_QDM.out',form='formatted',position='append')
      write(10,*) timesimu, meanadvecQDMsumT, meandiffvisQDMsumT, meandissvisQDMsumT, meanlorrQDMsumT, sommesumT
      close(10)

      open(10,file='mean_values_MHD.out',form='formatted',position='append')
      write(10,*) timesimu,meanadvecMHDsumT, meandiffvisMHDsumT, meandissvisMHDsumT, meanlorrMHDsumT,&
                 & meanforMHDsumT  ,sommesumTmhd
      close(10)

       open(10,file='meanDissGSEQDM.out',form='formatted',position='append')
        do i=1,mhdndelta
          write(10,*) timesimu, meanDissGSEQDMsumT(1,i),meanDissGSEQDMsumT(2,i)
        enddo
       close(10)

       open(10,file='meanAdvecGSEQDM.out',form='formatted',position='append')
        do i=1,mhdndelta
         write(10,*) timesimu, meanAdvecGSEQDMsumT(1,i), meanAdvecGSEQDMsumT(2,i)
        enddo
       close(10)

       open(10,file='meanDifGSEQDM.out',form='formatted',position='append')
        do i=1,mhdndelta
         write(10,*) timesimu, meanDifGSEQDMsumT(1,i), meanDifGSEQDMsumT(2,i)
        enddo
       close(10)

       open(10,file='meanTransEQDMGS.out',form='formatted',position='append')
        do i=1,mhdndelta
          write(10,*) timesimu, meanTransEQDMGSsumT(1,i), meanTransEQDMGSsumT(2,i)
        enddo
       close(10)

       open(10,file='meanTransEqdmSGS.out',form='formatted',position='append')
        do i=1,mhdndelta
         write(10,*)  timesimu, meanTransEqdmSGSsumT(1,i), meanTransEqdmSGSsumT(2,i)
        enddo
       close(10)

       open(10,file='meanDiffUSGSEqdm.out',form='formatted',position='append')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanDiffUSGSEqdmsumT(1,i), meanDiffUSGSEqdmsumT(2,i)
        enddo
       close(10)
      
       open(10,file='meanDiffBSGSEqdm.out',form='formatted',position='append')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanDiffBSGSEqdmsumT(1,i), meanDiffBSGSEqdmsumT(2,i)
        enddo
       close(10)

       open(10,file='meanDissBSGSEqdm.out',form='formatted',position='append')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanDissBSGSEqdmsumT(1,i), meanDissBSGSEqdmsumT(2,i)
        enddo
       close(10)

       open(10,file='meanDissUSGSEqdm.out',form='formatted',position='append')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanDissUSGSEqdmsumT(1,i), meanDissUSGSEqdmsumT(2,i)
        enddo
       close(10)


       open(10,file='meanDissUBSGSEMHDsumT.out',form='formatted',position='append')
        do i=1,mhdndelta
         write(10,*)  timesimu, meanDissUBSGSEMHDsumT(1,i), meanDissUBSGSEMHDsumT(2,i)
        enddo
       close(10)

       open(10,file='meanDiffUBSGSEMHDsumT.out',form='formatted',position='append')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanDiffUBSGSEMHDsumT(1,i), meanDiffUBSGSEMHDsumT(2,i)
        enddo
       close(10)
      
       open(10,file='meanAdvecGSEMHDsumT.out',form='formatted',position='append')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanAdvecGSEMHDsumT(1,i), meanAdvecGSEMHDsumT(2,i)
        enddo
       close(10)

       open(10,file='meanDifGSEMHDsumT.out',form='formatted',position='append')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanDifGSEMHDsumT(1,i), meanDifGSEMHDsumT(2,i)
        enddo
       close(10)
       open(10,file='meanProdEMHDsumT.out',form='formatted',position='append')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanProdEMHDsumT(1,i), meanProdEMHDsumT(2,i)
        enddo
       close(10)
       open(10,file='meanDissGSEMHDsumT.out',form='formatted',position='append')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanDissGSEMHDsumT(1,i), meanDissGSEMHDsumT(2,i)
        enddo
       close(10)
       endif

       open(10,file='meanDissSGSEMHDsumT.out',form='formatted',position='append')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanDissSGSEMHDsumT(1,i), meanDissSGSEMHDsumT(2,i)
        enddo
       close(10)

       open(10,file='meanDissSGSEQDMsumT.out',form='formatted',position='append')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanDissSGSEQDMsumT(1,i), meanDissSGSEQDMsumT(2,i)
        enddo
       close(10)

       open(10,file='meanEmTBbmodelSmag.out',form='formatted',position='append')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanEmTBbmodelSmag(1,i), meanEmTBbmodelSmag(2,i)
        enddo
       close(10)

       open(10,file='meanEqdmTUubSmag.out',form='formatted',position='append')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanEqdmTUubmodelSmag(1,i), meanEqdmTUubmodelSmag(2,i)
        enddo
       close(10)

       open(10,file='meanEqdmTUuuSmag.out',form='formatted',position='append')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanEqdmTUuumodelSmag(1,i), meanEqdmTUuumodelSmag(2,i)
        enddo
       close(10)

       open(10,file='meanEmTBbmodelGrad.out',form='formatted',position='append')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanEmTBbmodelGrad(1,i), meanEmTBbmodelGrad(2,i)
        enddo
       close(10)

       open(10,file='meanEqdmTUubGrad.out',form='formatted',position='append')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanEqdmTUubmodelGrad(1,i), meanEqdmTUubmodelGrad(2,i)
        enddo
       close(10)

       open(10,file='meanEqdmTUuuGrad.out',form='formatted',position='append')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanEqdmTUuumodelGrad(1,i), meanEqdmTUuumodelGrad(2,i)
        enddo
       close(10)

       open(10,file='meanEmTBbmodelSimil.out',form='formatted',position='append')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanEmTBbmodelSimil(1,i), meanEmTBbmodelSimil(2,i)
        enddo
       close(10)

       open(10,file='meanEqdmTUubSimil.out',form='formatted',position='append')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanEqdmTUubmodelSimil(1,i), meanEqdmTUubmodelSimil(2,i)
        enddo
       close(10)

       open(10,file='meanEqdmTUuuSimil.out',form='formatted',position='append')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanEqdmTUuumodelSimil(1,i), meanEqdmTUuumodelSimil(2,i)
        enddo
       close(10)


       open(10,file='meanEqdmTUuHamba.out',form='formatted',position='append')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanEqdmTUumodelHamba(1,i), meanEqdmTUumodelHamba(2,i)
        enddo
       close(10)

       open(10,file='meanEqdmTUuCarrati.out',form='formatted',position='append')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanEqdmTUumodelCarrati(1,i), meanEqdmTUumodelCarrati(2,i)
        enddo
       close(10)


       open(10,file='meanEqdmTBbCarrati.out',form='formatted',position='append')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanEmhdTBbmodelCarrati(1,i), meanEmhdTBbmodelCarrati(2,i)
        enddo
       close(10)

    endif



end subroutine computealltheterms

!======== EQDM==========
!------------------------------------------------------------------------------
!> Compute and save mean values of all quantities in kinetic energie equation.
!! @author Mouloud KESSAR
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in]    Uin,Vin,Win     = Velocity fields ins physical space
!!    @param[in]    Bin             = magnetic Field
!!    @param[in]    Wave            = wavenumber
!! @details
!!    
!------------------------------------------------------------------------------

subroutine eqdm(Uin,Vin,Win,Bin,Wave,nbcpus,spec_rank)


   implicit none

   integer, intent(in)    :: nbcpus,spec_rank
   TYPE(REAL_DATA_LAYOUT), INTENT(IN)          :: Uin,Vin,Win
   TYPE(REAL_DATA_LAYOUT),dimension(:), INTENT(IN) :: Bin(:)
   TYPE(WaveNumbers), INTENT(IN) :: Wave
   TYPE(REAL_DATA_LAYOUT) :: tab1, tab2, tab3
   TYPE(REAL_DATA_LAYOUT) :: dfxdx, dfxdy, dfxdz
   TYPE(REAL_DATA_LAYOUT) :: dfydx, dfydy, dfydz
   TYPE(REAL_DATA_LAYOUT) :: dfzdx, dfzdy, dfzdz
   TYPE(REAL_DATA_LAYOUT) :: advecQDM,diffvisQDM,dissvisQDM,lorrQDM,EnQDM,EnMHD
   REAL(WP) :: meanadvecQDM, meandiffvisQDM, meandissvisQDM, meanlorrQDM, somme 


     IF ((.NOT. copyStructOnly(Uin,advecQDM)) .OR. &
        &(.NOT. copyStructOnly(Uin,diffvisQDM)) .OR. &
        &(.NOT. copyStructOnly(Uin,dissvisQDM)) .OR. &
        &(.NOT. copyStructOnly(Uin,lorrQDM)) .OR. &
        &(.NOT. copyStructOnly(Uin,EnQDM)) .OR. &
        &(.NOT. copyStructOnly(Uin,EnMHD)) .OR. &
        &(.NOT. copyStructOnly(Uin,dfxdx)) .OR. &
        &(.NOT. copyStructOnly(Uin,dfxdy)) .OR. &
        &(.NOT. copyStructOnly(Uin,dfxdz)) .OR. &
        &(.NOT. copyStructOnly(Uin,dfydx)) .OR. &
        &(.NOT. copyStructOnly(Uin,dfydy)) .OR. &
        &(.NOT. copyStructOnly(Uin,dfydz)) .OR. &
        &(.NOT. copyStructOnly(Uin,dfzdx)) .OR. &
        &(.NOT. copyStructOnly(Uin,dfzdy)) .OR. &
        &(.NOT. copyStructOnly(Uin,dfzdz)) .OR. &
        &(.NOT. copyStructOnly(Uin,tab1)) .OR. &
        &(.NOT. copyStructOnly(Uin,tab2)) .OR. &
        &(.NOT. copyStructOnly(Uin,tab3)) ) RETURN


  
     !!Energie!!
    EnQDM%values = 0.5_WP*( Uin%values*Uin%values + Vin%values*Vin%values + Win%values*Win%values)
    EnMHD%values = 0.5_WP*( Bin(1)%values*Bin(1)%values + Bin(2)%values*Bin(2)%values + Bin(3)%values*Bin(3)%values)

    !!!Advection
    call computePhysicalGradient(Wave,EnQDM,dfxdx,dfxdy,dfxdz)
    advecQDM%values = -Uin%values*dfxdx%values -Vin%values*dfxdy%values -Win%values*dfxdz%values
 
  
    !!!Diffusion visqueuse
    tab1%values = dfxdx%values
    tab2%values = dfxdy%values
    tab3%values = dfxdz%values 
    call computePhysicalGradient(Wave,tab1,dfxdx,dfxdy,dfxdz)
    call computePhysicalGradient(Wave,tab2,dfydx,dfydy,dfydz)
    call computePhysicalGradient(Wave,tab3,dfzdx,dfzdy,dfzdz)
    diffvisQDM%values = mu*( dfxdx%values + dfydy%values + dfzdz%values )
  
    !!!Dissipation visqueuse
    call computePhysicalGradient(Wave,Uin,dfxdx,dfxdy,dfxdz)
    call computePhysicalGradient(Wave,Vin,dfydx,dfydy,dfydz)
    call computePhysicalGradient(Wave,Win,dfzdx,dfzdy,dfzdz)

    dissvisQDM%values = -2*mu*(dfxdx%values*dfxdx%values + dfxdy%values*dfxdy%values + dfxdz%values*dfxdz%values + &
                             & dfydx%values*dfydx%values + dfydy%values*dfydy%values + dfydz%values*dfydz%values + &
                             & dfzdx%values*dfzdx%values + dfzdy%values*dfzdy%values + dfzdz%values*dfzdz%values )
  
    !!!Force de lorrentz
    tab1%values = Bin(1)%values*Bin(1)%values
    tab2%values = Bin(1)%values*Bin(2)%values
    tab3%values = Bin(1)%values*Bin(3)%values
    call computePhysicalGradient(Wave,tab1,dfxdx,dfxdy,dfxdz)
    call computePhysicalGradient(Wave,tab2,dfydx,dfydy,dfydz)
    call computePhysicalGradient(Wave,tab3,dfzdx,dfzdy,dfzdz)

    lorrQDM%values = 2.*Uin%values*( dfxdx%values + dfydy%values + dfzdz%values )

    tab1%values = Bin(2)%values*Bin(1)%values
    tab2%values = Bin(2)%values*Bin(2)%values
    tab3%values = Bin(2)%values*Bin(3)%values
    call computePhysicalGradient(Wave,tab1,dfxdx,dfxdy,dfxdz)
    call computePhysicalGradient(Wave,tab2,dfydx,dfydy,dfydz)
    call computePhysicalGradient(Wave,tab3,dfzdx,dfzdy,dfzdz)

    lorrQDM%values = lorrQDM%values + 2.*Vin%values*(dfxdx%values + dfydy%values + dfzdz%values )

    tab1%values = Bin(3)%values*Bin(1)%values
    tab2%values = Bin(3)%values*Bin(2)%values
    tab3%values = Bin(3)%values*Bin(3)%values
    call computePhysicalGradient(Wave,tab1,dfxdx,dfxdy,dfxdz)
    call computePhysicalGradient(Wave,tab2,dfydx,dfydy,dfydz)
    call computePhysicalGradient(Wave,tab3,dfzdx,dfzdy,dfzdz)

    lorrQDM%values = lorrQDM%values + 2.*Win%values*( dfxdx%values + dfydy%values + dfzdz%values )
      
  
!     if (rank .eq. 0) print*,'Mean calculation'
    meanadvecQDM   = computeFieldAvg(advecQDM,spec_rank,nbcpus) 
    meandiffvisQDM = computeFieldAvg(diffvisQDM,spec_rank,nbcpus) 
    meandissvisQDM = computeFieldAvg(dissvisQDM,spec_rank,nbcpus) 
    meanlorrQDM    = computeFieldAvg(lorrQDM,spec_rank,nbcpus) 
    meanEk = computeFieldAvg(EnQDM,spec_rank,nbcpus) 
    meanEm = computeFieldAvg(EnMHD,spec_rank,nbcpus) 

    somme = meanadvecQDM + meandiffvisQDM + meandissvisQDM + meanlorrQDM
    meanTotE = meanEm + meanEk
    if(timesimu.EQ.0.0)  then
     meanadvecQDMsumT = meanadvecQDM
     meandiffvisQDMsumT = meandiffvisQDM
     meandissvisQDMsumT = meandissvisQDM
     meanlorrQDMsumT = meanlorrQDM
     sommesumT = somme
     meanTotEsumT = meanTotE
     meanEksumT = meanEk
     meanEmsumT = meanEm
    else
     meanadvecQDMsumT = ((timesimu-dt)*meanadvecQDMsumT + dt*meanadvecQDM)/timesimu
     meandiffvisQDMsumT = ((timesimu-dt)* meandiffvisQDMsumT + dt*meandiffvisQDM)/timesimu
     meandissvisQDMsumT = ((timesimu-dt)*meandissvisQDMsumT + dt*meandissvisQDM)/timesimu
     meanlorrQDMsumT = ((timesimu-dt)*meanlorrQDMsumT + dt*meanlorrQDM)/timesimu
     sommesumT = ((timesimu-dt)*sommesumT + dt*somme)/timesimu
     meanTotEsumT = ((timesimu-dt)*meanTotEsumT + dt*meanTotE)/timesimu
     meanEksumT = ((timesimu-dt)*meanEksumT + dt*meanEk)/timesimu
     meanEmsumT = ((timesimu-dt)*meanEmsumT + dt*meanEm)/timesimu

    endif

  CALL deleteDataLayout(advecQDM)
  CALL deleteDataLayout(diffvisQDM)
  CALL deleteDataLayout(dissvisQDM)
  CALL deleteDataLayout(lorrQDM)
  CALL deleteDataLayout(EnQDM)
  CALL deleteDataLayout(dfxdx)
  CALL deleteDataLayout(dfxdy)
  CALL deleteDataLayout(dfxdz)
  CALL deleteDataLayout(dfydx)
  CALL deleteDataLayout(dfydy)
  CALL deleteDataLayout(dfydz)
  CALL deleteDataLayout(dfzdx)
  CALL deleteDataLayout(dfzdy)
  CALL deleteDataLayout(dfzdz)
  CALL deleteDataLayout(tab1)
  CALL deleteDataLayout(tab2)
  CALL deleteDataLayout(tab3)

end subroutine eqdm
 
!======== EQMHD==========
!------------------------------------------------------------------------------
!> Compute and save mean values of all quantities in magnetic energie equation.
!! @author Mouloud KESSAR
!!    @param[in]    spec_rank     = mpi rank inside spectral communicator
!!    @param[in]    Uin,Vin,Win   = Velocity fields ins physical space
!!    @param[in]    Bin           = magnetic Field
!!    @param[in]    Wave          = wavenumber
!!    @param[in]    nbcpus        = numbers of cpus
!! @details
!!    
!------------------------------------------------------------------------------
subroutine eqmhd(Uin,Vin,Win,Bin,Wave,nbcpus,spec_rank)


 implicit none
 integer, intent(in) :: nbcpus,spec_rank
 TYPE(REAL_DATA_LAYOUT), INTENT(IN)          :: Uin,Vin,Win
 TYPE(REAL_DATA_LAYOUT), dimension(:), INTENT(IN) :: Bin(:)
 TYPE(WaveNumbers), INTENT(IN) :: Wave
 TYPE(REAL_DATA_LAYOUT) :: tab1, tab2, tab3
 TYPE(REAL_DATA_LAYOUT) :: dfxdx, dfxdy, dfxdz
 TYPE(REAL_DATA_LAYOUT) :: dfydx, dfydy, dfydz
 TYPE(REAL_DATA_LAYOUT) :: dfzdx, dfzdy, dfzdz
 TYPE(REAL_DATA_LAYOUT) :: advecMHD,diffvisMHD,dissvisMHD,lorrMHD,forMHD,EnMHD
 REAL(WP) :: meanadvecMHD, meandiffvisMHD, meandissvisMHD, meanlorrMHD, meanforMHD ,somme
  

     IF ((.NOT. copyStructOnly(Uin,advecMHD)) .OR. &
        &(.NOT. copyStructOnly(Uin,diffvisMHD)) .OR. &
        &(.NOT. copyStructOnly(Uin,dissvisMHD)) .OR. &
        &(.NOT. copyStructOnly(Uin,lorrMHD)) .OR. &
        &(.NOT. copyStructOnly(Uin,forMHD)) .OR. &
        &(.NOT. copyStructOnly(Uin,EnMHD)) .OR. &
        &(.NOT. copyStructOnly(Uin,dfxdx)) .OR. &
        &(.NOT. copyStructOnly(Uin,dfxdy)) .OR. &
        &(.NOT. copyStructOnly(Uin,dfxdz)) .OR. &
        &(.NOT. copyStructOnly(Uin,dfydx)) .OR. &
        &(.NOT. copyStructOnly(Uin,dfydy)) .OR. &
        &(.NOT. copyStructOnly(Uin,dfydz)) .OR. &
        &(.NOT. copyStructOnly(Uin,dfzdx)) .OR. &
        &(.NOT. copyStructOnly(Uin,dfzdy)) .OR. &
        &(.NOT. copyStructOnly(Uin,dfzdz)) .OR. &
        &(.NOT. copyStructOnly(Uin,tab1)) .OR. &
        &(.NOT. copyStructOnly(Uin,tab2)) .OR. &
        &(.NOT. copyStructOnly(Uin,tab3)) ) RETURN
 

  !!!Advection
  EnMHD%values = 0.5_WP*(Bin(1)%values*Bin(1)%values + Bin(2)%values*Bin(2)%values + Bin(3)%values*Bin(3)%values)
  call computePhysicalGradient(Wave,EnMHD,dfxdx,dfxdy,dfxdz)
  advecMHD%values = -Uin%values*dfxdx%values - Vin%values*dfxdy%values - Win%values*dfxdz%values
  
  
  !!!Diffusion visqueuse
  tab1%values = dfxdx%values
  tab2%values = dfxdy%values
  tab3%values = dfxdz%values
  call computePhysicalGradient(Wave,tab1,dfxdx,dfxdy,dfxdz)
  call computePhysicalGradient(Wave,tab2,dfydx,dfydy,dfydz)
  call computePhysicalGradient(Wave,tab3,dfzdx,dfzdy,dfzdz)
  diffvisMHD%values = eta*( dfxdx%values + dfydy%values + dfzdz%values )
  
  
  !!!Dissipation visqueuse
  call computePhysicalGradient(Wave,Bin(1),dfxdx,dfxdy,dfxdz)
  call computePhysicalGradient(Wave,Bin(2),dfydx,dfydy,dfydz)
  call computePhysicalGradient(Wave,Bin(3),dfzdx,dfzdy,dfzdz)
  dissvisMHD%values = -2.*eta*( dfxdx%values*dfxdx%values + dfxdy%values*dfxdy%values + dfxdz%values*dfxdz%values + &
                              & dfydx%values*dfydx%values + dfydy%values*dfydy%values + dfydz%values*dfydz%values + &
                              & dfzdx%values*dfzdx%values + dfzdy%values*dfzdy%values + dfzdz%values*dfzdz%values )

  !!Force de lorrentz
  tab1%values = Bin(1)%values*Bin(1)%values
  tab2%values = Bin(1)%values*Bin(2)%values
  tab3%values = Bin(1)%values*Bin(3)%values


  call computePhysicalGradient(Wave,tab1,dfxdx,dfxdy,dfxdz)
  call computePhysicalGradient(Wave,tab2,dfydx,dfydy,dfydz)
  call computePhysicalGradient(Wave,tab3,dfzdx,dfzdy,dfzdz)

  lorrMHD%values = -2.*Uin%values*( dfxdx%values + dfydy%values + dfzdz%values )
  tab1%values = Bin(2)%values*Bin(1)%values
  tab2%values = Bin(2)%values*Bin(2)%values
  tab3%values = Bin(2)%values*Bin(3)%values

  call computePhysicalGradient(Wave,tab1,dfxdx,dfxdy,dfxdz)
  call computePhysicalGradient(Wave,tab2,dfydx,dfydy,dfydz)
  call computePhysicalGradient(Wave,tab3,dfzdx,dfzdy,dfzdz)

  lorrMHD%values = lorrMHD%values - 2.*Vin%values*( dfxdx%values + dfydy%values + dfzdz%values )

  tab1%values = Bin(3)%values*Bin(1)%values
  tab2%values = Bin(3)%values*Bin(2)%values
  tab3%values = Bin(3)%values*Bin(3)%values

  call computePhysicalGradient(Wave,tab1,dfxdx,dfxdy,dfxdz)
  call computePhysicalGradient(Wave,tab2,dfydx,dfydy,dfydz)
  call computePhysicalGradient(Wave,tab3,dfzdx,dfzdy,dfzdz)

  lorrMHD%values = lorrMHD%values - 2.*Win%values*( dfxdx%values + dfydy%values + dfzdz%values )
 
  
  !!!Force
  tab1%values = Uin%values*Bin(1)%values*Bin(1)%values
  tab2%values = Uin%values*Bin(1)%values*Bin(2)%values
  tab3%values = Uin%values*Bin(1)%values*Bin(3)%values

  call computePhysicalGradient(Wave,tab1,dfxdx,dfxdy,dfxdz)
  call computePhysicalGradient(Wave,tab2,dfydx,dfydy,dfydz)
  call computePhysicalGradient(Wave,tab3,dfzdx,dfzdy,dfzdz)

  forMHD%values = 2.*( dfxdx%values + dfydy%values + dfzdz%values )

  tab1%values = Vin%values*Bin(2)%values*Bin(1)%values
  tab2%values = Vin%values*Bin(2)%values*Bin(2)%values
  tab3%values = Vin%values*Bin(2)%values*Bin(3)%values

  call computePhysicalGradient(Wave,tab1,dfxdx,dfxdy,dfxdz)
  call computePhysicalGradient(Wave,tab2,dfydx,dfydy,dfydz)
  call computePhysicalGradient(Wave,tab3,dfzdx,dfzdy,dfzdz)

  forMHD%values = forMHD%values + 2.*( dfxdx%values + dfydy%values + dfzdz%values )

  tab1%values = Win%values*Bin(3)%values*Bin(1)%values
  tab2%values = Win%values*Bin(3)%values*Bin(2)%values
  tab3%values = Win%values*Bin(3)%values*Bin(3)%values
  call computePhysicalGradient(Wave,tab1,dfxdx,dfxdy,dfxdz)
  call computePhysicalGradient(Wave,tab2,dfydx,dfydy,dfydz)
  call computePhysicalGradient(Wave,tab3,dfzdx,dfzdy,dfzdz)
  forMHD%values = forMHD%values + 2.*(dfxdx%values + dfydy%values + dfzdz%values )
 
 
 
!   if (rank .eq. 0) print*,'Mean calculation'

    meanadvecMHD   = computeFieldAvg(advecMHD,spec_rank,nbcpus) 
    meandiffvisMHD = computeFieldAvg(diffvisMHD,spec_rank,nbcpus) 
    meandissvisMHD = computeFieldAvg(dissvisMHD,spec_rank,nbcpus) 
    meanlorrMHD    = computeFieldAvg(lorrMHD,spec_rank,nbcpus)
    meanforMHD     = computeFieldAvg(forMHD,spec_rank,nbcpus)

    somme = meanadvecMHD + meandiffvisMHD + meandissvisMHD + meanlorrMHD + meanforMHD

    if(timesimu.EQ.0.0)  then
     meanadvecMHDsumT = meanadvecMHD
     meandiffvisMHDsumT = meandiffvisMHD
     meandissvisMHDsumT = meandissvisMHD
     meanlorrMHDsumT = meanlorrMHD
     meanforMHDsumT = meanforMHD
     sommesumTMHD = somme
    else
     meanadvecMHDsumT = ((timesimu-dt)*meanadvecMHDsumT + dt*meanadvecMHD)/timesimu
     meandiffvisMHDsumT = ((timesimu-dt)* meandiffvisMHDsumT + dt*meandiffvisMHD)/timesimu
     meandissvisMHDsumT = ((timesimu-dt)*meandissvisMHDsumT + dt*meandissvisMHD)/timesimu
     meanlorrMHDsumT = ((timesimu-dt)*meanlorrMHDsumT + dt*meanlorrMHD)/timesimu
     meanforMHDsumT = ((timesimu-dt)*meanforMHDsumT + dt*meanforMHD)/timesimu
     sommesumTMHD = ((timesimu-dt)*sommesumTMHD + dt*somme)/timesimu
    endif

  CALL deleteDataLayout(advecMHD)
  CALL deleteDataLayout(diffvisMHD)
  CALL deleteDataLayout(dissvisMHD)
  CALL deleteDataLayout(lorrMHD)
  CALL deleteDataLayout(forMHD)
  CALL deleteDataLayout(EnMHD)
  CALL deleteDataLayout(dfxdx)
  CALL deleteDataLayout(dfxdy)
  CALL deleteDataLayout(dfxdz)
  CALL deleteDataLayout(dfydx)
  CALL deleteDataLayout(dfydy)
  CALL deleteDataLayout(dfydz)
  CALL deleteDataLayout(dfzdx)
  CALL deleteDataLayout(dfzdy)
  CALL deleteDataLayout(dfzdz) 
  CALL deleteDataLayout(tab1)
  CALL deleteDataLayout(tab2)
  CALL deleteDataLayout(tab3)

end subroutine eqmhd



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


subroutine computeallGSandSGS(Uin,Vin,Win,Bin,Wave,nbcpus,spec_rank)

  implicit none
  TYPE(REAL_DATA_LAYOUT), INTENT(IN)          :: Uin,Vin,Win
  TYPE(REAL_DATA_LAYOUT), dimension(:), INTENT(IN) :: Bin(:)
  TYPE(REAL_DATA_LAYOUT) :: Bfx,Bfy,Bfz,Uf,Vf,Wf

  INTEGER :: nbcpus, spec_rank, ndelta
  TYPE(WAVENUMBERS) :: Wave
  REAL(WP) :: delta

  
  IF ((.NOT. copyStructOnly(Uin,Uf)) .OR. &
     &(.NOT. copyStructOnly(Uin,Vf)) .OR. &
     &(.NOT. copyStructOnly(Uin,Wf)) .OR. &
     &(.NOT. copyStructOnly(Bin(1),Bfx)) .OR. &
     &(.NOT. copyStructOnly(Bin(1),Bfy)) .OR. &
     &(.NOT. copyStructOnly(Bin(1),Bfz)) ) RETURN

    DO ndelta=1,mhdndelta
    delta = Uin%Lx*(2*ndelta)/Uin%nx

    meanTransEQDMGSsumT(1,ndelta)  = 2*ndelta
    meanTransEqdmSGSsumT(1,ndelta) = 2*ndelta
    meanDifGSEQDMsumT(1,ndelta)    = 2*ndelta
    meanAdvecGSEQDMsumT(1,ndelta)  = 2*ndelta
    meanDissGSEQDMsumT(1,ndelta)   = 2*ndelta
    meanDiffUSGSEqdmsumT(1,ndelta) = 2*ndelta
    meanDiffBSGSEqdmsumT(1,ndelta) = 2*ndelta
    meanDissUSGSEqdmsumT(1,ndelta) = 2*ndelta
    meanDissBSGSEqdmsumT(1,ndelta) = 2*ndelta
    meanDissSGSEQDMsumT(1,ndelta)  = 2*ndelta

    meanDifGSEMHDsumT(1,ndelta)     = 2*ndelta
    meanAdvecGSEMHDsumT(1,ndelta)   = 2*ndelta
    meanDissGSEMHDsumT(1,ndelta)    = 2*ndelta
    meanDissSGSEMHDsumT(1,ndelta)   = 2*ndelta
    meanDiffUBSGSEMHDsumT(1,ndelta) = 2*ndelta
    meanDissUBSGSEMHDsumT(1,ndelta) = 2*ndelta
    meanProdEMHDsumT(1,ndelta)      = 2*ndelta

    meanEkGSsumT(1,ndelta) = 2*ndelta
    meanEmGSsumT(1,ndelta) = 2*ndelta

    CALL computePhysicalFilter(Wave,Uin,Uf,delta,filterType)
    CALL computePhysicalFilter(Wave,Vin,Vf,delta,filterType)
    CALL computePhysicalFilter(Wave,Win,Wf,delta,filterType)
    CALL computePhysicalFilter(Wave,Bin(1),Bfx,delta,filterType)
    CALL computePhysicalFilter(Wave,Bin(2),Bfy,delta,filterType)
    CALL computePhysicalFilter(Wave,Bin(3),Bfz,delta,filterType)

    CALL EQDMGS(Uf,Vf,Wf,Bfx,Bfy,Bfz,Uin,Vin,Win,Bin,Wave,delta,ndelta,filterType,nbcpus,spec_rank)
    CALL EmhdGS(Uf,Vf,Wf,Bfx,Bfy,Bfz,Uin,Vin,Win,Bin,Wave,delta,ndelta,filterType,nbcpus,spec_rank)
    ENDDO

    CALL deleteDataLayout(Uf)
    CALL deleteDataLayout(Vf)
    CALL deleteDataLayout(Wf)
    CALL deleteDataLayout(Bfx)
    CALL deleteDataLayout(Bfy)
    CALL deleteDataLayout(Bfz)

end subroutine computeallGSandSGS
!======= EQDMGS =========
!------------------------------------------------------------------------------
!> Compute and save mean values of all quantities in 
!! @author Mouloud KESSAR
!!    @param[in]    spec_rank     = mpi rank inside spectral communicator
!!    @param[in]    spec_rank     = mpi rank inside spectral communicator
!!    @param[in]    Uin,Vin,Win   = Velocity fields ins physical space
!!    @param[in]    Bin           = magnetic Field
!!    @param[in]    Wave          = wavenumber
!!    @param[in]    nbcpus        = numbers of cpus
!!    @param[in]    spec_rank     = mpi rank inside spectral communicator
!!    @param[in]    Ufin,Vfin,Wfin   = filtered Velocity fields ins physical space
!!    @param[in]    Bfxin, Bfyin, Bfzin    = filtered magnetic Field
!!    @param[in]    delta     = filtersize
!!    @param[in]    ndelta    = filter number
!!    @param[in]    filter    = filter type
!! @details
!!    
!------------------------------------------------------------------------------


subroutine EQDMGS(Ufin,Vfin,Wfin,Bfxin,Bfyin,Bfzin,Uin,Vin,Win,Bin,Wave,delta,ndelta,filter,nbcpus,spec_rank)

  IMPLICIT NONE
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Bfxin,Bfyin,Bfzin,Ufin,Vfin,Wfin
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Uin,Vin,Win
  TYPE(REAL_DATA_LAYOUT), dimension(:), INTENT(IN) :: Bin(:)
  TYPE(REAL_DATA_LAYOUT) :: DissGSEQDM, AdvecGSEQDM, DifGSEQDM, TransEQDMGS, TransEqdmSGS,EkGS,EmGS
  TYPE(REAL_DATA_LAYOUT) :: DiffUSGSEqdm, DissUSGSEqdm, DiffBSGSEqdm, DissBSGSEqdm,DissSGSEQDM
  REAL(WP) :: delta
  REAL(WP) :: meanDissGSEQDM, meanAdvecGSEQDM, meanDifGSEQDM, meanTransEQDMGS, meanTransEqdmSGS
  REAL(WP) :: meanDiffUSGSEqdm, meanDissUSGSEqdm, meanDiffBSGSEqdm, meanDissBSGSEqdm, meanDissSGSEQDM 
  REAL(WP) :: meanEkGS,meanEmGS
  INTEGER, INTENT(IN) :: nbcpus, spec_rank, filter,ndelta
  TYPE(REAL_DATA_LAYOUT) :: dfxdx, dfxdy, dfxdz
  TYPE(REAL_DATA_LAYOUT) :: dfydx, dfydy, dfydz
  TYPE(REAL_DATA_LAYOUT) :: dfzdx, dfzdy, dfzdz
  TYPE(REAL_DATA_LAYOUT) :: tab1, tab2, tab3
  type(waveNumbers), intent(in)        :: wave

     IF ((.NOT. copyStructOnly(Uin,dfxdx)) .OR. &
        &(.NOT. copyStructOnly(Uin,dfxdy)) .OR. &
        &(.NOT. copyStructOnly(Uin,dfxdz)) .OR. &
        &(.NOT. copyStructOnly(Uin,dfydx)) .OR. &
        &(.NOT. copyStructOnly(Uin,dfydy)) .OR. &
        &(.NOT. copyStructOnly(Uin,dfydz)) .OR. &
        &(.NOT. copyStructOnly(Uin,dfzdx)) .OR. &
        &(.NOT. copyStructOnly(Uin,dfzdy)) .OR. &
        &(.NOT. copyStructOnly(Uin,dfzdz)) .OR. &
        &(.NOT. copyStructOnly(Uin,tab1)) .OR. &
        &(.NOT. copyStructOnly(Uin,tab2)) .OR. &
        &(.NOT. copyStructOnly(Uin,tab3)) ) RETURN

     IF ( (.NOT. copyStructOnly(Uin,DissGSEQDM ))  .OR. &
         &(.NOT. copyStructOnly(Uin,TransEQDMGS))  .OR. &
         &(.NOT. copyStructOnly(Uin,DissSGSEQDM )) .OR. &
         &(.NOT. copyStructOnly(Uin,TransEqdmSGS)) .OR. &
         &(.NOT. copyStructOnly(Uin,AdvecGSEQDM))  .OR. &
         &(.NOT. copyStructOnly(Uin,EkGS))  .OR. &
         &(.NOT. copyStructOnly(Uin,EmGS))  .OR. &
         &(.NOT. copyStructOnly(Uin,DifGSEQDM)) ) RETURN

  EkGS%values = 0.5*( Ufin%values**2 + Vfin%values**2 + Wfin%values**2  )

  EmGS%values = 0.5*( Bfxin%values**2 + Bfyin%values**2 + Bfxin%values**2 )

  CALL computePhysicalGradient(wave,Ufin,dfxdx,dfxdy,dfxdz)
  CALL computePhysicalGradient(wave,Vfin,dfydx,dfydy,dfydz)
  CALL computePhysicalGradient(wave,Wfin,dfzdx,dfzdy,dfzdz)

  DissGSEQDM%values = -2.0*mu*( dfxdx%values**2 + dfxdy%values**2 + dfxdz%values**2 &
                            & + dfydx%values**2 + dfydy%values**2 + dfydz%values**2 &
                            & + dfzdx%values**2 + dfzdy%values**2 + dfzdz%values**2 &
                            &   )

 
  CALL computePhysicalGradient(wave,Uin,dfxdx,dfxdy,dfxdz)
  CALL computePhysicalGradient(wave,Vin,dfydx,dfydy,dfydz)
  CALL computePhysicalGradient(wave,Win,dfzdx,dfzdy,dfzdz)


  tab1%values  = -2.0*mu*( dfxdx%values**2 + dfxdy%values**2 + dfxdz%values**2 &
                       & + dfydx%values**2 + dfydy%values**2 + dfydz%values**2 &
                       & + dfzdx%values**2 + dfzdy%values**2 + dfzdz%values**2 &
                       &   )

  CALL computePhysicalFilter(wave,tab1,DissSGSEQDM,delta,filter)

  DissSGSEQDM%values= DissSGSEQDM%values - DissGSEQDM%values 


  CALL computePhysicalGradient(wave,Bfxin,dfxdx,dfxdy,dfxdz)
  CALL computePhysicalGradient(wave,Bfyin,dfydx,dfydy,dfydz)
  CALL computePhysicalGradient(wave,Bfzin,dfzdx,dfzdy,dfzdz)


  TransEQDMGS%values = 2.0*Ufin%values*( Bfxin%values*dfxdx%values + Bfyin%values*dfxdy%values + Bfzin%values*dfxdz%values) &
                   & + 2.0*Vfin%values*( Bfxin%values*dfydx%values + Bfyin%values*dfydy%values + Bfzin%values*dfydz%values) &
                   & + 2.0*Wfin%values*( Bfxin%values*dfzdx%values + Bfyin%values*dfzdy%values + Bfzin%values*dfzdz%values) 
  

  CALL computePhysicalGradient(wave,Bin(1),dfxdx,dfxdy,dfxdz)
  CALL computePhysicalGradient(wave,Bin(2),dfydx,dfydy,dfydz)
  CALL computePhysicalGradient(wave,Bin(3),dfzdx,dfzdy,dfzdz)

  tab1%values=2*( Uin%values*Bin(1)%values*dfxdx%values &
              & + Uin%values*Bin(2)%values*dfxdy%values &
              & + Uin%values*Bin(3)%values*dfxdz%values &
              & + Vin%values*Bin(1)%values*dfydx%values &
              & + Vin%values*Bin(2)%values*dfydy%values &
              & + Vin%values*Bin(3)%values*dfydz%values &
              & + Win%values*Bin(1)%values*dfzdx%values &
              & + Win%values*Bin(2)%values*dfzdy%values &
              & + Win%values*Bin(3)%values*dfzdz%values &
              & )

  CALL computePhysicalFilter(wave,tab1,tab2,delta,filter)

  TransEqdmSGS%values= tab2%values - TransEqdmGS%values

  tab1%values = Ufin%values*Ufin%values + Vfin%values*Vfin%values + Wfin%values*Wfin%values
  call computePhysicalGradient(wave,tab1,dfxdx,dfxdy,dfxdz)
  AdvecGSEQDM%values = Ufin%values*dfxdx%values + Vfin%values*dfxdy%values + Wfin%values*dfxdz%values

  tab1%values = dfxdx%values
  tab2%values = dfxdy%values
  tab3%values = dfxdz%values
  call computePhysicalGradient(wave,tab1,dfxdx,dfxdy,dfxdz)
  call computePhysicalGradient(wave,tab2,dfydx,dfydy,dfydz)
  call computePhysicalGradient(wave,tab3,dfzdx,dfzdy,dfzdz)

  DifGSEQDM%values = mu*( dfxdx%values + dfydy%values + dfzdz%values )

  CALL deleteDataLayout(dfxdx)
  CALL deleteDataLayout(dfxdy)
  CALL deleteDataLayout(dfxdz)
  CALL deleteDataLayout(dfydx)
  CALL deleteDataLayout(dfydy)
  CALL deleteDataLayout(dfydz)
  CALL deleteDataLayout(dfzdx)
  CALL deleteDataLayout(dfzdy)
  CALL deleteDataLayout(dfzdz)
  CALL deleteDataLayout(tab1)
  CALL deleteDataLayout(tab2)
  CALL deleteDataLayout(tab3)

  IF ((.NOT. copyStructOnly(Uin,DiffUSGSEqdm)) .OR. &
     &(.NOT. copyStructOnly(Uin,DissUSGSEqdm)) .OR. &
     &(.NOT. copyStructOnly(Uin,DiffBSGSEqdm)) .OR. &
     &(.NOT. copyStructOnly(Uin,DissBSGSEqdm)) ) RETURN



  DiffUSGSEqdm%values = 0.0_WP
  DissUSGSEqdm%values = 0.0_WP

  !Tau_xx^u
  CALL computetauijA(Ufin,Uin,Uin,Uin,Uin,Ufin,Ufin,1,1,wave,DiffUSGSEqdm,DissUSGSEqdm,delta,filter,-1.0_WP)

  !Tau_xy^u
  CALL computetauijA(Ufin,Uin,Vin,Uin,Vin,Ufin,Vfin,1,2,wave,DiffUSGSEqdm,DissUSGSEqdm,delta,filter,-1.0_WP)

  !Tau_xz^u
  CALL computetauijA(Ufin,Uin,Win,Uin,Win,Ufin,Wfin,1,3,wave,DiffUSGSEqdm,DissUSGSEqdm,delta,filter,-1.0_WP)

  !Tau_yx^u
  CALL computetauijA(Vfin,Vin,Uin,Vin,Uin,Vfin,Ufin,2,1,wave,DiffUSGSEqdm,DissUSGSEqdm,delta,filter,-1.0_WP)

  !Tau_yy^u
  CALL computetauijA(Vfin,Vin,Vin,Vin,Vin,Vfin,Vfin,2,2,wave,DiffUSGSEqdm,DissUSGSEqdm,delta,filter,-1.0_WP)

  !Tau_yz^u
  CALL computetauijA(Vfin,Vin,Win,Vin,Win,Vfin,Wfin,2,3,wave,DiffUSGSEqdm,DissUSGSEqdm,delta,filter,-1.0_WP)

  !Tau_zx^u
  CALL computetauijA(Wfin,Win,Uin,Win,Uin,Wfin,Ufin,3,1,wave,DiffUSGSEqdm,DissUSGSEqdm,delta,filter,-1.0_WP)

  !Tau_zy^u
  CALL computetauijA(Wfin,Win,Vin,Win,Vin,Wfin,Vfin,3,2,wave,DiffUSGSEqdm,DissUSGSEqdm,delta,filter,-1.0_WP)

  !Tau_zz^u
  CALL computetauijA(Wfin,Win,Win,Win,Win,Wfin,Wfin,3,3,wave,DiffUSGSEqdm,DissUSGSEqdm,delta,filter,-1.0_WP)

  DiffBSGSEqdm%values = 0.0_WP
  DissBSGSEqdm%values = 0.0_WP

  !Tau_xx^b
  CALL computetauijA(Ufin,Uin,Uin,Bin(1),Bin(1),Bfxin,Bfxin,1,1,wave,DiffBSGSEqdm,DissBSGSEqdm,delta,filter,1.0_WP)

  !Tau_xy^b
  CALL computetauijA(Ufin,Uin,Vin,Bin(1),Bin(2),Bfxin,Bfyin,1,2,wave,DiffBSGSEqdm,DissBSGSEqdm,delta,filter,1.0_WP)

  !Tau_xz^b
  CALL computetauijA(Ufin,Uin,Win,Bin(1),Bin(3),Bfxin,Bfzin,1,3,wave,DiffBSGSEqdm,DissBSGSEqdm,delta,filter,1.0_WP)

  !Tau_yx^b
  CALL computetauijA(Vfin,Vin,Uin,Bin(2),Bin(1),Bfyin,Bfxin,2,1,wave,DiffBSGSEqdm,DissBSGSEqdm,delta,filter,1.0_WP)

  !Tau_yy^b
  CALL computetauijA(Vfin,Vin,Vin,Bin(2),Bin(2),Bfyin,Bfyin,2,2,wave,DiffBSGSEqdm,DissBSGSEqdm,delta,filter,1.0_WP)

  !Tau_yz^b
  CALL computetauijA(Vfin,Vin,Win,Bin(2),Bin(3),Bfyin,Bfzin,2,3,wave,DiffBSGSEqdm,DissBSGSEqdm,delta,filter,1.0_WP)

  !Tau_zx^b
  CALL computetauijA(Wfin,Win,Uin,Bin(3),Bin(1),Bfzin,Bfxin,3,1,wave,DiffBSGSEqdm,DissBSGSEqdm,delta,filter,1.0_WP)

  !Tau_zy^b
  CALL computetauijA(Wfin,Win,Vin,Bin(3),Bin(2),Bfzin,Bfyin,3,2,wave,DiffBSGSEqdm,DissBSGSEqdm,delta,filter,1.0_WP)

  !Tau_zz^b
  CALL computetauijA(Wfin,Win,Win,Bin(3),Bin(3),Bfzin,Bfzin,3,3,wave,DiffBSGSEqdm,DissBSGSEqdm,delta,filter,1.0_WP)


  meanTransEQDMGS  = computeFieldAvg(TransEQDMGS, spec_rank,nbcpus)
  meanTransEqdmSGS = computeFieldAvg(TransEqdmSGS, spec_rank,nbcpus)
  meanDifGSEQDM    = computeFieldAvg(DifGSEQDM,  spec_rank,nbcpus)
  meanAdvecGSEQDM  = computeFieldAvg(AdvecGSEQDM, spec_rank,nbcpus)
  meanDissGSEQDM   = computeFieldAvg(DissGSEQDM, spec_rank,nbcpus)
  meanDiffUSGSEqdm = computeFieldAvg(DiffUSGSEqdm, spec_rank,nbcpus)
  meanDiffBSGSEqdm = computeFieldAvg(DiffBSGSEqdm, spec_rank,nbcpus)
  meanDissUSGSEqdm = computeFieldAvg(DissUSGSEqdm, spec_rank,nbcpus)
  meanDissBSGSEqdm = computeFieldAvg(DissBSGSEqdm, spec_rank,nbcpus)
  meanDissSGSEQDM   = computeFieldAvg(DissSGSEQDM, spec_rank,nbcpus)
  meanEkGS  = computeFieldAvg(EkGS, spec_rank,nbcpus)
  meanEmGS  = computeFieldAvg(EmGS, spec_rank,nbcpus)

  if(timesimu.EQ.0.0)  then

     meanTransEQDMGSsumT(2,ndelta)  = meanTransEQDMGS
     meanTransEqdmSGSsumT(2,ndelta) = meanTransEqdmSGS
     meanDifGSEQDMsumT(2,ndelta)    = meanDifGSEQDM
     meanAdvecGSEQDMsumT(2,ndelta)  = meanAdvecGSEQDM
     meanDissGSEQDMsumT(2,ndelta)   = meanDissGSEQDM
     meanDiffUSGSEqdmsumT(2,ndelta) = meanDiffUSGSEqdm
     meanDiffBSGSEqdmsumT(2,ndelta) = meanDiffBSGSEqdm
     meanDissUSGSEqdmsumT(2,ndelta) = meanDissUSGSEqdm
     meanDissBSGSEqdmsumT(2,ndelta) = meanDissBSGSEqdm
     meanDissSGSEQDMsumT(2,ndelta)   = meanDissSGSEQDM
     meanEkGSsumT(2,ndelta)   = meanEkGS
     meanEmGSsumT(2,ndelta)   = meanEmGS
  else

     meanTransEQDMGSsumT(2,ndelta)  = ( (timesimu-dt)*meanTransEQDMGSsumT(2,ndelta) + dt*meanTransEQDMGS )/timesimu
     meanTransEqdmSGSsumT(2,ndelta) = ( (timesimu-dt)*meanTransEqdmSGSsumT(2,ndelta) + dt*meanTransEqdmSGS )/timesimu
     meanDifGSEQDMsumT(2,ndelta)    = ( (timesimu-dt)*meanDifGSEQDMsumT(2,ndelta) + dt*meanDifGSEQDM )/timesimu
     meanAdvecGSEQDMsumT(2,ndelta)  = ( (timesimu-dt)*meanAdvecGSEQDMsumT(2,ndelta) + dt*meanAdvecGSEQDM )/timesimu
     meanDissGSEQDMsumT(2,ndelta)   = ( (timesimu-dt)*meanDissGSEQDMsumT(2,ndelta) + dt*meanDissGSEQDM )/timesimu
     meanDiffUSGSEqdmsumT(2,ndelta) = ( (timesimu-dt)*meanDiffUSGSEqdmsumT(2,ndelta) + dt*meanDiffUSGSEqdm )/timesimu
     meanDiffBSGSEqdmsumT(2,ndelta) = ( (timesimu-dt)*meanDiffBSGSEqdmsumT(2,ndelta) + dt*meanDiffBSGSEqdm )/timesimu
     meanDissUSGSEqdmsumT(2,ndelta) = ( (timesimu-dt)*meanDissUSGSEqdmsumT(2,ndelta) + dt*meanDissUSGSEqdm )/timesimu
     meanDissBSGSEqdmsumT(2,ndelta) = ( (timesimu-dt)*meanDissBSGSEqdmsumT(2,ndelta) + dt*meanDissBSGSEqdm )/timesimu
     meanDissSGSEQDMsumT(2,ndelta)  = ( (timesimu-dt)*meanDissSGSEQDMsumT(2,ndelta) + dt*meanDissSGSEQDM )/timesimu
     meanEkGSsumT(2,ndelta)  = ( (timesimu-dt)*meanEkGSsumT(2,ndelta) + dt*meanEkGS )/timesimu
     meanEmGSsumT(2,ndelta)  = ( (timesimu-dt)*meanEmGSsumT(2,ndelta) + dt*meanEmGS )/timesimu

  endif

  CALL deleteDataLayout(TransEQDMGS)
  CALL deleteDataLayout(TransEqdmSGS)
  CALL deleteDataLayout(DifGSEQDM)
  CALL deleteDataLayout(AdvecGSEQDM)
  CALL deleteDataLayout(DissGSEQDM)
  CALL deleteDataLayout(DiffUSGSEqdm)
  CALL deleteDataLayout(DiffBSGSEqdm)
  CALL deleteDataLayout(DissUSGSEqdm)
  CALL deleteDataLayout(DissBSGSEqdm)
  CALL deleteDataLayout(EkGS)
  CALL deleteDataLayout(EmGS)

end subroutine EQDMGS

!======= EmhdGS =========
!------------------------------------------------------------------------------
!> Compute and save mean values of all quantities in
!! @author Mouloud KESSAR
!!    @param[in]    spec_rank     = mpi rank inside spectral communicator
!!    @param[in]    spec_rank     = mpi rank inside spectral communicator
!!    @param[in]    Uin,Vin,Win   = Velocity fields ins physical space
!!    @param[in]    Bin           = magnetic Field
!!    @param[in]    Wave          = wavenumber
!!    @param[in]    nbcpus        = numbers of cpus
!!    @param[in]    spec_rank     = mpi rank inside spectral communicator
!!    @param[in]    Ufin,Vfin,Wfin   = filtered Velocity fields ins physical space
!!    @param[in]    Bfxin, Bfyin, Bfzin    = filtered magnetic Field
!!    @param[in]    delta     = filtersize
!!    @param[in]    ndelta    = filter number
!!    @param[in]    filter    = filter type
!! @details
!!
!------------------------------------------------------------------------------

subroutine EmhdGS(Ufin,Vfin,Wfin,Bfxin,Bfyin,Bfzin,Uin,Vin,Win,Bin,Wave,delta,ndelta,filter,nbcpus,spec_rank)

  IMPLICIT NONE
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Bfxin,Bfyin,Bfzin,Ufin,Vfin,Wfin
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Uin,Vin,Win
  TYPE(REAL_DATA_LAYOUT), dimension(:), INTENT(IN) :: Bin(:)
  TYPE(REAL_DATA_LAYOUT) :: DissGSEMHD, AdvecGSEMHD, DifGSEMHD
  TYPE(REAL_DATA_LAYOUT) :: DiffUBSGSEMHD, DissUBSGSEMHD,ProdEMhd, DissSGSEMHD
  REAL(WP) :: delta
  REAL(WP) :: meanDissGSEMHD, meanAdvecGSEMHD, meanDifGSEMHD,meanDiffUBSGSEMHD
  REAL(WP) :: meanDissUBSGSEMHD,meanProdEMHD , meanDissSGSEMHD 
  INTEGER, INTENT(IN) :: nbcpus, spec_rank, filter,ndelta
  TYPE(REAL_DATA_LAYOUT) :: dfxdx, dfxdy, dfxdz
  TYPE(REAL_DATA_LAYOUT) :: dfydx, dfydy, dfydz
  TYPE(REAL_DATA_LAYOUT) :: dfzdx, dfzdy, dfzdz
  TYPE(REAL_DATA_LAYOUT) :: tab1, tab2, tab3
  type(waveNumbers), intent(in)        :: wave

     IF ((.NOT. copyStructOnly(Uin,dfxdx)) .OR. &
        &(.NOT. copyStructOnly(Uin,dfxdy)) .OR. &
        &(.NOT. copyStructOnly(Uin,dfxdz)) .OR. &
        &(.NOT. copyStructOnly(Uin,dfydx)) .OR. &
        &(.NOT. copyStructOnly(Uin,dfydy)) .OR. &
        &(.NOT. copyStructOnly(Uin,dfydz)) .OR. &
        &(.NOT. copyStructOnly(Uin,dfzdx)) .OR. &
        &(.NOT. copyStructOnly(Uin,dfzdy)) .OR. &
        &(.NOT. copyStructOnly(Uin,dfzdz)) .OR. &
        &(.NOT. copyStructOnly(Uin,tab1)) .OR. &
        &(.NOT. copyStructOnly(Uin,tab2)) .OR. &
        &(.NOT. copyStructOnly(Uin,tab3)) ) RETURN

     IF ( (.NOT. copyStructOnly(Uin,DissGSEMHD )) .OR. &
         &(.NOT. copyStructOnly(Uin,AdvecGSEMHD)) .OR. &
         &(.NOT. copyStructOnly(Uin,DissSGSEMHD)) .OR. &
         &(.NOT. copyStructOnly(Uin,ProdEMHD)) .OR. &
         &(.NOT. copyStructOnly(Uin,DifGSEMHD)) ) RETURN

  CALL computePhysicalGradient(wave,Bfxin,dfxdx,dfxdy,dfxdz)
  CALL computePhysicalGradient(wave,Bfyin,dfydx,dfydy,dfydz)
  CALL computePhysicalGradient(wave,Bfzin,dfzdx,dfzdy,dfzdz)

  DissGSEMHD%values = -2.0*eta*( dfxdx%values**2 + dfxdy%values**2 + dfxdz%values**2 &
                             & + dfydx%values**2 + dfydy%values**2 + dfydz%values**2 & 
                             & + dfzdx%values**2 + dfzdy%values**2 + dfzdz%values**2 &
                             & )



  tab1%values = Bfxin%values*Bfxin%values + Bfyin%values*Bfyin%values + Bfzin%values*Bfzin%values
  call computePhysicalGradient(wave,tab1,dfxdx,dfxdy,dfxdz)
  AdvecGSEMHD%values = Ufin%values*dfxdx%values + Vfin%values*dfxdy%values + Wfin%values*dfxdz%values

  tab1%values = dfxdx%values 
  tab2%values = dfxdy%values
  tab3%values = dfxdz%values
  call computePhysicalGradient(wave,tab1,dfxdx,dfxdy,dfxdz)
  call computePhysicalGradient(wave,tab2,dfydx,dfydy,dfydz)
  call computePhysicalGradient(wave,tab3,dfzdx,dfzdy,dfzdz)

  DifGSEMHD%values = eta*( dfxdx%values + dfydy%values + dfzdz%values )

  !!!Dissipation visqueuse
  call computePhysicalGradient(Wave,Bin(1),dfxdx,dfxdy,dfxdz)
  call computePhysicalGradient(Wave,Bin(2),dfydx,dfydy,dfydz)
  call computePhysicalGradient(Wave,Bin(3),dfzdx,dfzdy,dfzdz)

  tab1%values = -2.*eta*( dfxdx%values*dfxdx%values + dfxdy%values*dfxdy%values + dfxdz%values*dfxdz%values + &
                        & dfydx%values*dfydx%values + dfydy%values*dfydy%values + dfydz%values*dfydz%values + &
                        & dfzdx%values*dfzdx%values + dfzdy%values*dfzdy%values + dfzdz%values*dfzdz%values )

  CALL computePhysicalFilter(wave,tab1,tab2,delta,filter) 

  DissSGSEMHD%values = tab2%values - DissGSEMHD%values 


  tab1%values = Ufin%values*Bfxin%values + Vfin%values*Bfyin%values + Wfin%values*Bfzin%values 
  call computePhysicalGradient(wave,tab1,dfxdx,dfxdy,dfxdz)

  ProdEMHD%values = 2.0*( dfxdx%values*Bfxin%values + dfxdy%values*Bfyin%values + dfxdz%values*Bfzin%values )
  
  CALL deleteDataLayout(dfxdx)
  CALL deleteDataLayout(dfxdy)
  CALL deleteDataLayout(dfxdz)
  CALL deleteDataLayout(dfydx)
  CALL deleteDataLayout(dfydy)
  CALL deleteDataLayout(dfydz)
  CALL deleteDataLayout(dfzdx)
  CALL deleteDataLayout(dfzdy)
  CALL deleteDataLayout(dfzdz)
  CALL deleteDataLayout(tab1)
  CALL deleteDataLayout(tab2)
  CALL deleteDataLayout(tab3)
  
  IF ((.NOT. copyStructOnly(Uin,DiffUBSGSEMHD)) .OR. &
     &(.NOT. copyStructOnly(Uin,DissUBSGSEMHD)) ) RETURN



  DiffUBSGSEMHD%values = 0.0_WP
  DissUBSGSEMHD%values = 0.0_WP

  !Tau_xx^ub
  CALL computetauijAB(Uin,Uin,Ufin,Ufin,Bin(1),Bin(1),Bfxin,Bfxin,1,1,wave,&
                     & DiffUBSGSEMHD,DissUBSGSEMHD,delta,filter)

  !Tau_xy^ub
  CALL computetauijAB(Uin,Vin,Ufin,Vfin,Bin(1),Bin(2),Bfxin,Bfyin,1,2,wave,&
                    & DiffUBSGSEMHD,DissUBSGSEMHD,delta,filter)

  !Tau_xz^ub
  CALL computetauijAB(Uin,Win,Ufin,Wfin,Bin(1),Bin(3),Bfxin,Bfzin,1,3,wave,&
                     &DiffUBSGSEMHD,DissUBSGSEMHD,delta,filter)

  !Tau_yx^ub
  CALL computetauijAB(Vin,Uin,Vfin,Ufin,Bin(2),Bin(1),Bfyin,Bfxin,2,1,wave,&
                    &DiffUBSGSEMHD,DissUBSGSEMHD,delta,filter)

  !Tau_yy^ub
  CALL computetauijAB(Vin,Vin,Vfin,Vfin,Bin(2),Bin(2),Bfyin,Bfyin,2,2,wave,&
                     &DiffUBSGSEMHD,DissUBSGSEMHD,delta,filter)

  !Tau_yz^ub
  CALL computetauijAB(Vin,Win,Vfin,Wfin,Bin(2),Bin(3),Bfyin,Bfzin,2,3,wave,&
                     &DiffUBSGSEMHD,DissUBSGSEMHD,delta,filter)

  !Tau_zx^ub
  CALL computetauijAB(Win,Uin,Wfin,Ufin,Bin(3),Bin(1),Bfzin,Bfxin,3,1,wave,&
                     & DiffUBSGSEMHD,DissUBSGSEMHD,delta,filter)

  !Tau_zy^ub
  CALL computetauijAB(Win,Vin,Wfin,Vfin,Bin(3),Bin(2),Bfzin,Bfyin,3,2,wave,&
                    & DiffUBSGSEMHD,DissUBSGSEMHD,delta,filter)

  !Tau_zz^ub
  CALL computetauijAB(Win,Win,Wfin,Wfin,Bin(3),Bin(3),Bfzin,Bfzin,3,3,wave,&
                    & DiffUBSGSEMHD,DissUBSGSEMHD,delta,filter)



  meanDifGSEMHD     = computeFieldAvg(DifGSEMHD,  spec_rank,nbcpus)
  meanAdvecGSEMHD   = computeFieldAvg(AdvecGSEMHD, spec_rank,nbcpus)
  meanDissGSEMHD    = computeFieldAvg(DissGSEMHD, spec_rank,nbcpus)
  meanDiffUBSGSEMHD = computeFieldAvg(DiffUBSGSEMHD, spec_rank,nbcpus)
  meanDissUBSGSEMHD = computeFieldAvg(DissUBSGSEMHD, spec_rank,nbcpus)
  meanProdEMHD      = computeFieldAvg(ProdEMHD, spec_rank,nbcpus)
  meanDissSGSEMHD    = computeFieldAvg(DissSGSEMHD, spec_rank,nbcpus)

    if(timesimu.EQ.0.0)  then
     meanDifGSEMHDsumT(2,ndelta)     = meanDifGSEMHD
     meanAdvecGSEMHDsumT(2,ndelta)   = meanAdvecGSEMHD
     meanDissGSEMHDsumT(2,ndelta)    = meanDissGSEMHD
     meanDiffUBSGSEMHDsumT(2,ndelta) = meanDiffUBSGSEMHD
     meanDissUBSGSEMHDsumT(2,ndelta) = meanDissUBSGSEMHD
     meanProdEMHDsumT(2,ndelta)      = meanProdEMHD
     meanDissSGSEMHDsumT(2,ndelta)    = meanDissSGSEMHD
    else
     meanDifGSEMHDsumT(2,ndelta)     = ( (timesimu-dt)*meanDifGSEMHDsumT(2,ndelta) + dt*meanDifGSEMHD )/timesimu
     meanAdvecGSEMHDsumT(2,ndelta)   = ( (timesimu-dt)*meanAdvecGSEMHDsumT(2,ndelta) + dt*meanAdvecGSEMHD )/timesimu
     meanDissGSEMHDsumT(2,ndelta)    = ( (timesimu-dt)*meanDissGSEMHDsumT(2,ndelta) + dt*meanDissGSEMHD )/timesimu
     meanDiffUBSGSEMHDsumT(2,ndelta) = ( (timesimu-dt)*meanDiffUBSGSEMHDsumT(2,ndelta) + dt*meanDiffUBSGSEMHD )/timesimu
     meanDissUBSGSEMHDsumT(2,ndelta) = ( (timesimu-dt)*meanDissUBSGSEMHDsumT(2,ndelta) + dt*meanDissUBSGSEMHD )/timesimu
     meanProdEMHDsumT(2,ndelta)      = ( (timesimu-dt)*meanProdEMHDsumT(2,ndelta) + dt*meanProdEMHD )/timesimu
     meanDissSGSEMHDsumT(2,ndelta)    = ( (timesimu-dt)*meanDissSGSEMHDsumT(2,ndelta)+  dt*meanDissSGSEMHD)/timesimu
    endif

  CALL deleteDataLayout(DifGSEMHD) 
  CALL deleteDataLayout(AdvecGSEMHD) 
  CALL deleteDataLayout(DissGSEMHD)
  CALL deleteDataLayout(ProdEMHD)
  CALL deleteDataLayout(DiffUBSGSEMHD) 
  CALL deleteDataLayout(DissUBSGSEMHD)
  CALL deleteDataLayout(DissSGSEMHD) 

end subroutine EmhdGS

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

  INTEGER :: spec_rank,nbcpus
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

    DO ndelta=1,mhdndelta
    delta = Uin%Lx*(2*ndelta)/Uin%nx
    !!!! Smagorinsky
    meanEmTBbmodelSmag(1,ndelta)    = 2*ndelta
    meanEqdmTUuumodelSmag(1,ndelta) = 2*ndelta
    meanEqdmTUubmodelSmag(1,ndelta) = 2*ndelta
    !!!! Simil
    meanEmTBbmodelSimil(1,ndelta) = 2*ndelta
    meanEqdmTUuumodelSimil(1,ndelta) = 2*ndelta
    meanEqdmTUubmodelSimil(1,ndelta) = 2*ndelta
    !!!!Gradient
    meanEmTBbmodelGrad(1,ndelta)  = 2*ndelta
    meanEqdmTUuumodelGrad(1,ndelta) = 2*ndelta
    meanEqdmTUubmodelGrad(1,ndelta) = 2*ndelta

    meanEqdmTUumodelCarrati(1,ndelta) = 2*ndelta
    meanEmhdTBbmodelCarrati(1,ndelta) = 2*ndelta

    meanEqdmTUumodelHamba(1,ndelta) = 2*ndelta

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
                    & Bxkf,Bykf,Bzkf,Wave,spec_rank,delta,&
                    & meanEmTBbmodelSmag(2,ndelta),meanEmTBbmodelGrad(2,ndelta),meanEmTBbmodelSimil(2,ndelta),meanEmhdTBbmodelCarrati(2,ndelta),nbcpus)

    call  aprioriEQDmMHD(Uf,Vf,Wf,Bxf,Byf,Bzf,Ukf,Vkf,Wkf,&
                    & Bxkf,Bykf,Bzkf,Wave,spec_rank,delta,meanEqdmTUuumodelSmag(2,ndelta),&
                    & meanEqdmTUubmodelSmag(2,ndelta),meanEqdmTUuumodelSimil(2,ndelta),&
                    & meanEqdmTUubmodelSimil(2,ndelta),meanEqdmTUuumodelGrad(2,ndelta),&
                    & meanEqdmTUubmodelGrad(2,ndelta),meanEqdmTUumodelCarrati(2,ndelta),meanEqdmTUumodelHamba(2,ndelta),nbcpus)

    enddo

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


end subroutine

subroutine aprioriMHDEm(Uin,Vin,Win,Bxin,Byin,Bzin,Ukin,Vkin,Wkin,&
                    & Bxkin,Bykin,Bzkin,Wave,spec_rank,delta,meanEmTBbSmag,meanEmTBbGradient,meanEmTBbSimilarity,meanEmTBbCarrati,nbcpus)

  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Uin,Vin,Win
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Bxin,Byin,Bzin 
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Ukin,Vkin,Wkin
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Bxkin,Bykin,Bzkin
  TYPE(WaveNumbers), INTENT(IN) :: Wave

  INTEGER :: spec_rank,nbcpus
  TYPE(REAL_DATA_LAYOUT) :: T12,T13,T23
  TYPE(REAL_DATA_LAYOUT) :: J12,J13,J23
  TYPE(REAL_DATA_LAYOUT) :: EmTBbmodel
  logical :: res
  REAL(WP) :: delta,meanEmTBbSmag,meanEmTBbGradient,meanEmTBbSimilarity,meanEmTBbCarrati

  IF ((.NOT. copyStructOnly(Bxin,J12)) .OR. &
     &(.NOT. copyStructOnly(Bxin,J13)) .OR. &
     &(.NOT. copyStructOnly(Bxin,J23)) .OR. &
     &(.NOT. copyStructOnly(Bxin,EmTBbmodel)) .OR. &
     &(.NOT. copyStructOnly(Bxin,t12)) .OR. &
     &(.NOT. copyStructOnly(Bxin,t13)) .OR. &
     &(.NOT. copyStructOnly(Bxin,t23))    ) RETURN 


     call computeJijTensor(BxKin,ByKin,BzKin,wave,J12,J13,J23, res)


     call smagorinskyMHDTij(T12,T13,T23,Bxkin,Bykin,Bzkin,&
                           &wave,spec_rank,res, delta)

!!!!! facteur 4, car on a un deux dans l'expression et le tenseur est antisymétrique
     EmTBbmodel%values = 4.0_WP*(J12%values*T12%values &
                               &+J13%values*T13%values &
                               &+J23%values*T23%values )
    meanEmTBbSmag  = computeFieldAvg(EmTBbmodel, spec_rank,nbcpus)
    !!!!!Gradient
    call InducGradient(Uin,Vin,Win,Bxin,Byin,Bzin,&
                      &wave,delta,T12,T13,T23,filterType,spec_rank )
     EmTBbmodel%values = 4.0_WP*(J12%values*T12%values &
                               &+J13%values*T13%values &
                               &+J23%values*T23%values )
    meanEmTBbGradient  = computeFieldAvg(EmTBbmodel, spec_rank,nbcpus)
    !!! Similarity
    call SimilInduc(T12,T13,T23,Uin,Vin,Win,Ukin,Vkin,Wkin,Bxin,Byin,Bzin,Bxkin,Bykin,Bzkin,&
                    &wave,delta,filterType,spec_rank)
     EmTBbmodel%values = 4.0_WP*(J12%values*T12%values &
                               &+J13%values*T13%values &
                               &+J23%values*T23%values )
    meanEmTBbSimilarity  = computeFieldAvg(EmTBbmodel, spec_rank,nbcpus)

    call MullerCarratiInduc(Ukin,Vkin,Wkin,Bxkin,Bykin,Bzkin,wave,spec_rank,res, &
                          &  delta,T23,T12,T13)
     EmTBbmodel%values = 4.0_WP*(J12%values*T12%values &
                               &+J13%values*T13%values &
                               &+J23%values*T23%values )
    meanEmTBbCarrati  = computeFieldAvg(EmTBbmodel, spec_rank,nbcpus)

    CALL deleteDataLayout(J12)
    CALL deleteDataLayout(J13)
    CALL deleteDataLayout(J23)
    CALL deleteDataLayout(T12)
    CALL deleteDataLayout(T13)
    CALL deleteDataLayout(T23)
    CALL deleteDataLayout(EmTBbmodel)



end subroutine

subroutine aprioriEQDmMHD(Uin,Vin,Win,Bxin,Byin,Bzin,Ukin,Vkin,Wkin,&
                    & Bxkin,Bykin,Bzkin,Wave,spec_rank,delta,meanEqdmTUuuSmag,&
                    & meanEqdmTUubSmag,meanEqdmTUuuSimil,meanEqdmTUubSimil,meanEqdmTUuuGrad, meanEqdmTUubGrad,&
                    & meanEqdmTUuCarrati,meanEqdmTUuHamba,nbcpus)

  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Uin,Vin,Win
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Bxin,Byin,Bzin 
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Ukin,Vkin,Wkin
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Bxkin,Bykin,Bzkin
  TYPE(WaveNumbers), INTENT(IN) :: Wave
  
  TYPE(COMPLEX_DATA_LAYOUT) :: div1,div2,div3
  INTEGER :: spec_rank,nbcpus,numVel,type_avg
  TYPE(REAL_DATA_LAYOUT) :: T12,T13,T23
  TYPE(REAL_DATA_LAYOUT) :: T11,T22,T33
  TYPE(REAL_DATA_LAYOUT) :: S12,S13,S23
  TYPE(REAL_DATA_LAYOUT) :: S11,S22,S33
  TYPE(REAL_DATA_LAYOUT) :: Wtab
  logical :: res
  REAL(WP) :: delta,meanEqdmTUuuSmag, meanEqdmTUubSmag,meanEqdmTUuuSimil,meanEqdmTUubSimil
  REAL(WP) :: meanEqdmTUuuGrad, meanEqdmTUubGrad,meanEqdmTUuHamba,meanEqdmTUuCarrati

  IF ((.NOT. copyStructOnly(Bxin,S12)) .OR. &
     &(.NOT. copyStructOnly(Bxin,S13)) .OR. &
     &(.NOT. copyStructOnly(Bxin,S23)) .OR. &
     &(.NOT. copyStructOnly(Bxin,S11)) .OR. &
     &(.NOT. copyStructOnly(Bxin,S22)) .OR. &
     &(.NOT. copyStructOnly(Bxin,S33)) .OR. &
     &(.NOT. copyStructOnly(Bxin,Wtab)) .OR. &
     &(.NOT. copyStructOnly(Bxin,T12)) .OR. &
     &(.NOT. copyStructOnly(Bxin,T13)) .OR. &
     &(.NOT. copyStructOnly(Bxin,T23))  .OR. &
     &(.NOT. copyStructOnly(Bxin,T11)) .OR. &
     &(.NOT. copyStructOnly(Bxin,T22)) .OR. &
     &(.NOT. copyStructOnly(Bxin,T33))   ) RETURN 

  IF ((.NOT. copyStructOnly(Bxkin,div1)) .OR. &
     &(.NOT. copyStructOnly(Bxkin,div2)) .OR. &
     &(.NOT. copyStructOnly(Bxkin,div3))   ) RETURN 

     numVel=1
     type_avg=0
     call computeStressTensor(Ukin,VKin,WKin,wave,S11,S12,S13,S22,S23,S33,res)

     call smagorinskyVel(div1,div2,div3,Uin,Vin,Win,Ukin,Vkin,Wkin,wave,numVel,spec_rank,res, &
                          &  filterType,type_avg,T11,T12,T13,T22,T23,T33,delta)


     Wtab%values = 4.0_WP*(   S12%values*T12%values &
                          & + S13%values*T13%values &
                          & + S23%values*T23%values) &
               & + 2.0_WP*(   S11%values*T11%values &
                          & + S22%values*T22%values &
                          & + S23%values*T33%values )

    meanEqdmTUuuSmag  = computeFieldAvg(Wtab, spec_rank,nbcpus)


     call smagorinskyVel(div1,div2,div3,Bxin,Byin,Bzin,Bxkin,Bykin,Bzkin,wave,numVel,spec_rank,res, &
                          &  filterType,type_avg,T11,T12,T13,T22,T23,T33,delta)


     Wtab%values = - 4.0_WP*(  S12%values*T12%values &
                           & + S13%values*T13%values &
                           & + S23%values*T23%values ) &
                 & - 2.0_WP*(  S11%values*T11%values &
                           & + S22%values*T22%values &
                           & + S23%values*T33%values )

    meanEqdmTUubSmag = computeFieldAvg(Wtab, spec_rank,nbcpus)

    !!!!!!similarity

    call SimilVel(t11,t12,t13,t22,t23,t33,Uin,Vin,Win,Ukin,Vkin,Wkin,Wave,delta,filterType,spec_rank,nbcpus)

     Wtab%values = 4.0_WP*(   S12%values*T12%values &
                          & + S13%values*T13%values &
                          & + S23%values*T23%values) &
               & + 2.0_WP*(   S11%values*T11%values &
                          & + S22%values*T22%values &
                          & + S23%values*T33%values )

    meanEqdmTUuuSimil  = computeFieldAvg(Wtab, spec_rank,nbcpus)

    call SimilVel(t11,t12,t13,t22,t23,t33,Bxin,Byin,Bzin,Bxkin,Bykin,Bzkin,Wave,delta,filterType,spec_rank,nbcpus)

     Wtab%values = - 4.0_WP*(   S12%values*T12%values &
                            & + S13%values*T13%values &
                            & + S23%values*T23%values) &
                 & - 2.0_WP*(   S11%values*T11%values &
                            & + S22%values*T22%values &
                            & + S23%values*T33%values )

    meanEqdmTUubSimil  = computeFieldAvg(Wtab, spec_rank,nbcpus)

    !!!!!!Gradient

!     call SimilVel(t11,t12,t13,t22,t23,t33,Uin,Vin,Win,Ukin,Vkin,Wkin,Wave,delta,filterType,spec_rank,nbcpus)
    call gradientVel(t11,t12,t13,t22,t23,t33,Uin,Vin,Win,Wave,delta,spec_rank,nbcpus)
     Wtab%values = 4.0_WP*(   S12%values*T12%values &
                          & + S13%values*T13%values &
                          & + S23%values*T23%values) &
               & + 2.0_WP*(   S11%values*T11%values &
                          & + S22%values*T22%values &
                          & + S23%values*T33%values )

    meanEqdmTUuuGrad  = computeFieldAvg(Wtab, spec_rank,nbcpus)

!     call SimilVel(t11,t12,t13,t22,t23,t33,Bxin,Byin,Bzin,Bxkin,Bykin,Bzkin,Wave,delta,filterType,spec_rank,nbcpus)
    call gradientVel(t11,t12,t13,t22,t23,t33,Bxin,Byin,Bzin,Wave,delta,spec_rank,nbcpus)
     Wtab%values = - 4.0_WP*(   S12%values*T12%values &
                            & + S13%values*T13%values &
                            & + S23%values*T23%values) &
                 & - 2.0_WP*(   S11%values*T11%values &
                            & + S22%values*T22%values &
                            & + S23%values*T33%values )

    meanEqdmTUubGrad  = computeFieldAvg(Wtab, spec_rank,nbcpus)
    call MullerCarratiNSequation(Ukin,Vkin,Wkin,Bxkin,Bykin,Bzkin,wave,spec_rank,res, &
                          &  delta,T23,T12,T13,T11,T22,T33)

         Wtab%values = 4.0_WP*(   S12%values*T12%values &
                          & + S13%values*T13%values &
                          & + S23%values*T23%values) &
               & + 2.0_WP*(   S11%values*T11%values &
                          & + S22%values*T22%values &
                          & + S23%values*T33%values )

    meanEqdmTUuCarrati = computeFieldAvg(Wtab, spec_rank,nbcpus)

    call HambaYoshizawaNSequation(Ukin,Vkin,Wkin,Bxkin,Bykin,Bzkin,C_Lambda,C_mu,wave,spec_rank,res, &
                          &  delta,T11,T22,T33,T12,T13,T23)
    Wtab%values = 4.0_WP*(    S12%values*T12%values &
                          & + S13%values*T13%values &
                          & + S23%values*T23%values) &
               & + 2.0_WP*(   S11%values*T11%values &
                          & + S22%values*T22%values &
                          & + S23%values*T33%values )

    meanEqdmTUuHamba = computeFieldAvg(Wtab, spec_rank,nbcpus)


    CALL deleteDataLayout(S11)
    CALL deleteDataLayout(S22)
    CALL deleteDataLayout(S33)
    CALL deleteDataLayout(S12)
    CALL deleteDataLayout(S13)
    CALL deleteDataLayout(S23)

    CALL deleteDataLayout(T12)
    CALL deleteDataLayout(T13)
    CALL deleteDataLayout(T23)
    CALL deleteDataLayout(T11)
    CALL deleteDataLayout(T22)
    CALL deleteDataLayout(T33)

    CALL deleteDataLayout(Wtab)

    CALL deleteDataLayout(div3)
    CALL deleteDataLayout(div2)
    CALL deleteDataLayout(div1)

end subroutine

!======== computePhysicalGradient==========
!--------------------------------------------------------------------------
!> @details 
!> Compute Gradient in physical space
!! @author  Mouloud KESSAR, LEGI
!! @param[in]    wave waves number for field valIn
!! @param[in]    valIn component of scalar field in physical space
!! @param[out]   valOut1 first component of gradient in physical space
!! @param[out]   valOut2 second component of gradient in physical space
!! @param[out]   valOut3 third component of gradient in physical space
!---------------------------------------------------------------------------
 subroutine computePhysicalGradient(wave,valIn,valOut1,valOut2,valOut3)

   implicit none
   !IO
   type(REAL_data_layout), intent(in)    :: valIn
   type(REAL_data_layout), intent(inout) :: valOut1,valOut2,valOut3
   type(waveNumbers), intent(in)        :: wave

! XXX nbcpus n'est pas utilisé ! A enlever ... le flag -Wall permet de détecter
! les variables inutilisée (il y en a pas mal).
   logical :: res

   CALL ftran(valIn,valInk,res)
   IF (.NOT.res) RETURN

   call computeGradientKmain(wave,valInk,valOut1k,valOut2k,valOut3k)

   CALL btran(valOut1k,valOut1,res)
   IF (.NOT.res) RETURN
   CALL btran(valOut2k,valOut2,res)
   IF (.NOT.res) RETURN
   CALL btran(valOut3k,valOut3,res)
   IF (.NOT.res) RETURN

 end subroutine computePhysicalGradient

!======== computePhysicalFilter==========
!---------------------------------------------------------------------------
!> @details 
!> Compute filter field from physical space
!! @author  Mouloud KESSAR, LEGI
!! @param[in]    wave waves number for field valIn
!! @param[in]    valIn component of scalar field in physical space
!! @param[out]   valOut1 first component of gradient in physical space
!! @param[out]   valOut2 second component of gradient in physical space
!! @param[out]   valOut3 third component of gradient in physical space
!---------------------------------------------------------------------------
 subroutine computePhysicalFilter(wave,valIn,valOut,delta,filter)

   implicit none

   type(REAL_data_layout), intent(in)    :: valIn
   type(REAL_data_layout), intent(inout) :: valOut
   type(waveNumbers), intent(in)        :: wave
   REAL(WP) :: delta
   !Internal
   integer :: filter
   logical :: res

   CALL ftran(valIn,valInk,res)
   IF (.NOT.res) RETURN

   call computeFilter(Wave, delta, valInk, valOut1k, filter)

   CALL btran(valOut1k,valOut,res)
   IF (.NOT.res) RETURN


 end subroutine computePhysicalFilter

!========= computetau_ijA =========

 subroutine computetauijA(UifIn,UiIn,UjIn,AiIn,AjIn,AfiIn,AfjIn,Iin,Jin,wave,Diff,Diss,delta,filter,SgnIn)

  implicit none
  TYPE(REAL_data_layout), intent(in)     :: AiIn,AjIn,AfiIn,AfjIn,UifIn,UiIn,UjIn
  TYPE(REAL_data_layout), intent(inout)  :: Diff,Diss
  TYPE(REAL_data_layout) :: dfidx,dfidy,dfidz,dfjdx,dfjdy,dfjdz
  TYPE(REAL_data_layout) :: tabA,tabB
  INTEGER, INTENT(IN) :: filter,Iin,Jin
  type(waveNumbers), intent(in)        :: wave
  REAL(WP), intent(in) :: delta,SgnIn

     IF ((.NOT. copyStructOnly(AiIn,dfidx)) .OR. &
        &(.NOT. copyStructOnly(AiIn,dfidy)) .OR. &
        &(.NOT. copyStructOnly(AiIn,dfidz)) .OR. &
        &(.NOT. copyStructOnly(AiIn,dfjdx)) .OR. &
        &(.NOT. copyStructOnly(AiIn,dfjdy)) .OR. &
        &(.NOT. copyStructOnly(AiIn,dfjdz)) .OR. &
        &(.NOT. copyStructOnly(AiIn,tabA))  .OR. &
        &(.NOT. copyStructOnly(AiIn,tabB)) ) RETURN

  tabA%values = AiIn%values*AjIn%values
  CALL computePhysicalFilter(wave,tabA,tabB,delta,filter)
  tabA%values = tabB%values - AfiIn%values*AfjIn%values !T_ij^A
  
  tabB%values = tabA%values*UifIn%values
  CALL computePhysicalGradient(wave,tabB,dfidx,dfidy,dfidz)
  if(Jin.EQ.1) Diff%values = Diff%values + SgnIn*2.0*dfidx%values
  if(Jin.EQ.2) Diff%values = Diff%values + SgnIn*2.0*dfidy%values
  if(Jin.EQ.3) Diff%values = Diff%values + SgnIn*2.0*dfidz%values

  CALL computePhysicalGradient(wave,UiIn,dfidx,dfidy,dfidz)
  CALL computePhysicalGradient(wave,UjIn,dfjdx,dfjdy,dfjdz)
  if(Jin.EQ.1) then
    if (Iin.EQ.1) Diss%values = Diss%values - SgnIn*(dfidx%values + dfjdx%values)*tabA%values
    if (Iin.EQ.2) Diss%values = Diss%values - SgnIn*(dfidx%values + dfjdy%values)*tabA%values
    if (Iin.EQ.3) Diss%values = Diss%values - SgnIn*(dfidx%values + dfjdz%values)*tabA%values
  endif  
  if(Jin.EQ.2) then
    if (Iin.EQ.1) Diss%values = Diss%values - SgnIn*(dfidy%values + dfjdx%values)*tabA%values
    if (Iin.EQ.2) Diss%values = Diss%values - SgnIn*(dfidy%values + dfjdy%values)*tabA%values
    if (Iin.EQ.3) Diss%values = Diss%values - SgnIn*(dfidy%values + dfjdz%values)*tabA%values
  endif
  if(Jin.EQ.3) then
    if (Iin.EQ.1) Diss%values = Diss%values - SgnIn*(dfidz%values + dfjdx%values)*tabA%values
    if (Iin.EQ.2) Diss%values = Diss%values - SgnIn*(dfidz%values + dfjdy%values)*tabA%values
    if (Iin.EQ.3) Diss%values = Diss%values - SgnIn*(dfidz%values + dfjdz%values)*tabA%values
  endif


  CALL deleteDataLayout(dfidx)
  CALL deleteDataLayout(dfidy)
  CALL deleteDataLayout(dfidz)
  CALL deleteDataLayout(dfjdx)
  CALL deleteDataLayout(dfjdy)
  CALL deleteDataLayout(dfjdz)
  CALL deleteDataLayout(tabA)
  CALL deleteDataLayout(tabB)
 
 end subroutine computetauijA

!========= computetau_ijAB =========
 subroutine computetauijAB(AiIn,AjIn,AfiIn,AfjIn,BiIn,BjIn,BfiIn,BfjIn,Iin,Jin,wave,Diff,Diss,delta,filter)

  implicit none
  TYPE(REAL_data_layout), intent(in)     :: AiIn,AjIn,AfiIn,AfjIn
  TYPE(REAL_data_layout), intent(in)     :: BiIn,BjIn,BfiIn,BfjIn
  TYPE(REAL_data_layout), intent(inout)  :: Diff,Diss
  TYPE(REAL_data_layout) :: dfidx,dfidy,dfidz,dfjdx,dfjdy,dfjdz
  TYPE(REAL_data_layout) :: tabA,tabB,tabC
  INTEGER, INTENT(IN) :: filter,Iin,Jin
  type(waveNumbers), intent(in)        :: wave
  REAL(WP) :: delta

     IF ((.NOT. copyStructOnly(AiIn,dfidx)) .OR. &
        &(.NOT. copyStructOnly(AiIn,dfidy)) .OR. &
        &(.NOT. copyStructOnly(AiIn,dfidz)) .OR. &
        &(.NOT. copyStructOnly(AiIn,dfjdx)) .OR. &
        &(.NOT. copyStructOnly(AiIn,dfjdy)) .OR. &
        &(.NOT. copyStructOnly(AiIn,dfjdz)) .OR. &
        &(.NOT. copyStructOnly(AiIn,tabA))  .OR. &
        &(.NOT. copyStructOnly(AiIn,tabB))  .OR. &
        &(.NOT. copyStructOnly(AiIn,tabC)) ) RETURN

  tabA%values = BiIn%values*AjIn%values
  CALL computePhysicalFilter(wave,tabA,tabB,delta,filter)
  tabC%values = tabB%values - BfiIn%values*AfjIn%values

  tabA%values = AiIn%values*BjIn%values
  CALL computePhysicalFilter(wave,tabA,tabB,delta,filter)
  tabC%values = tabC%values - tabB%values + AfiIn%values*BfjIn%values !Tau_ij^UB
  
  tabB%values = tabC%values*BfiIn%values
  CALL computePhysicalGradient(wave,tabB,dfidx,dfidy,dfidz)
  if(Jin.EQ.1) Diff%values = Diff%values - 2.0*dfidx%values
  if(Jin.EQ.2) Diff%values = Diff%values - 2.0*dfidy%values
  if(Jin.EQ.3) Diff%values = Diff%values - 2.0*dfidz%values

  CALL computePhysicalGradient(wave,BfiIn,dfidx,dfidy,dfidz)
  CALL computePhysicalGradient(wave,BfjIn,dfjdx,dfjdy,dfjdz)
  if(Jin.EQ.1) then
    if (Iin.EQ.1) Diss%values = Diss%values + (dfidx%values - dfjdx%values)*tabC%values
    if (Iin.EQ.2) Diss%values = Diss%values + (dfidx%values - dfjdy%values)*tabC%values
    if (Iin.EQ.3) Diss%values = Diss%values + (dfidx%values - dfjdz%values)*tabC%values
  endif  
  if(Jin.EQ.2) then
    if (Iin.EQ.1) Diss%values = Diss%values + (dfidy%values - dfjdx%values)*tabC%values
    if (Iin.EQ.2) Diss%values = Diss%values + (dfidy%values - dfjdy%values)*tabC%values
    if (Iin.EQ.3) Diss%values = Diss%values + (dfidy%values - dfjdz%values)*tabC%values
  endif
  if(Jin.EQ.3) then
    if (Iin.EQ.1) Diss%values = Diss%values + (dfidz%values - dfjdx%values)*tabC%values
    if (Iin.EQ.2) Diss%values = Diss%values + (dfidz%values - dfjdy%values)*tabC%values
    if (Iin.EQ.3) Diss%values = Diss%values + (dfidz%values - dfjdz%values)*tabC%values
  endif


  CALL deleteDataLayout(dfidx)
  CALL deleteDataLayout(dfidy)
  CALL deleteDataLayout(dfidz)
  CALL deleteDataLayout(dfjdx)
  CALL deleteDataLayout(dfjdy)
  CALL deleteDataLayout(dfjdz)
  CALL deleteDataLayout(tabA)
  CALL deleteDataLayout(tabB)
 
  CALL deleteDataLayout(tabC)
 
 end subroutine computetauijAB

end module post_mhd
!> @}
