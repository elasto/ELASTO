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

module post_hd

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
    use post_hd_ener_lib
    implicit none
    private

     ! ==== Public procedures ==== 
    public          :: computeallthetermsHD
    public          :: post_hd_init
    ! ======= Mean terms ========
    real(WP), PROTECTED, SAVE :: meanadvecQDMsumT, meandiffvisQDMsumT, meandissvisQDMsumT, sommesumT 
    real(WP), PROTECTED, SAVE :: meanEksumT,meanEk
    real(WP), PROTECTED, SAVE :: eta,mu,Prm,dt,timesimu
    TYPE(COMPLEX_DATA_LAYOUT), PROTECTED, SAVE :: Valout3k,Valout2k,Valout1k,ValInk
    real(WP), PROTECTED, ALLOCATABLE, dimension(:,:), SAVE :: meanDissGSEQDMsumT, meanAdvecGSEQDMsumT, &
                                                            & meanDifGSEQDMsumT
    real(WP), PROTECTED, ALLOCATABLE, dimension(:,:), SAVE :: meanDiffUSGSEqdmsumT, meanDissUSGSEqdmsumT, &
                                                            & meanDissSGSEQDMsumT
    real(WP), PROTECTED, ALLOCATABLE, dimension(:,:), SAVE :: meanEkGSsumT

    real(WP), PROTECTED, ALLOCATABLE, dimension(:,:), SAVE :: meanEqdmTUuumodelSmag


    INTEGER, PROTECTED, SAVE :: mhdndelta,filterType
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

subroutine post_hd_init(muIn,Ukin)
   implicit none
   REAL(WP), intent(in)                     :: muIn
   TYPE(COMPLEX_DATA_LAYOUT), intent(in)    :: Ukin
   IF ((.NOT. copyStructOnly(UkIn,valOut1k)) .OR. &

      &(.NOT. copyStructOnly(UkIn,valOut2k)) .OR. &
      &(.NOT. copyStructOnly(UkIn,valOut3k)) .OR. &
      &(.NOT. copyStructOnly(UkIn,valInk)) ) RETURN

  mu=muIn
  CALL parser_read('HD filter size',mhdndelta)
  CALL parser_read('HD filter type',filterType)
  allocate( meanDifGSEQDMsumT(2,mhdndelta) )    
  allocate( meanAdvecGSEQDMsumT(2,mhdndelta) )  
  allocate( meanDissGSEQDMsumT(2,mhdndelta) )  
  allocate( meanDiffUSGSEqdmsumT(2,mhdndelta) )
  allocate( meanDissUSGSEqdmsumT(2,mhdndelta) )
  allocate( meanDissSGSEQDMsumT(2,mhdndelta) )

  allocate( meanEkGSsumT(2,mhdndelta) )


  allocate(meanEqdmTUuumodelSmag(2,mhdndelta))



end subroutine  post_hd_init



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
subroutine computeallthetermsHD(Uin,Vin,Win,Ukin,Vkin,Wkin,&
                         & Wave,nbcpus,spec_rank,simtime,tstep)


    implicit none

    TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Uin,Vin,Win
    TYPE(WaveNumbers),INTENT(IN) :: Wave
    TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Ukin,Vkin,Wkin
    real(WP):: meanadvecQDM, meandiffvisQDM, meandissvisQDM, somme
    INTEGER :: i,nbcpus,spec_rank
    REAL(WP) :: tstep,simtime
    dt=tstep
    timesimu=simtime
 

    call eqdm(Uin,Vin,Win,Ukin,Vkin,Wkin,Wave,nbcpus,spec_rank,mu,&
             &meanadvecQDM, meandiffvisQDM, meandissvisQDM, somme,meanEk)
    if(timesimu.EQ.0.0)  then
     meanadvecQDMsumT = meanadvecQDM
     meandiffvisQDMsumT = meandiffvisQDM
     meandissvisQDMsumT = meandissvisQDM
     sommesumT = somme
     meanEksumT = meanEk
    else
     meanadvecQDMsumT = ((timesimu-dt)*meanadvecQDMsumT + dt*meanadvecQDM)/timesimu
     meandiffvisQDMsumT = ((timesimu-dt)* meandiffvisQDMsumT + dt*meandiffvisQDM)/timesimu
     meandissvisQDMsumT = ((timesimu-dt)*meandissvisQDMsumT + dt*meandissvisQDM)/timesimu
     sommesumT = ((timesimu-dt)*sommesumT + dt*somme)/timesimu
     meanEksumT = ((timesimu-dt)*meanEksumT + dt*meanEk)/timesimu

    endif

    if(spec_rank.EQ.0) then

     if (timesimu .LT. 0.00000001)then
       open(10,file='mean_values_QDM.out',form='formatted')
       write(10,*) timesimu, meanadvecQDMsumT, meandiffvisQDMsumT, meandissvisQDMsumT,&
                 & sommesumT,meanEk
       close(10)
 
      else

      open(10,file='mean_values_QDM.out',form='formatted',position='append')
      write(10,*) timesimu, meanadvecQDMsumT, meandiffvisQDMsumT, meandissvisQDMsumT, sommesumT,meanEk
      close(10)


     endif
    endif

    call computeallGSandSGSHD(Uin,Vin,Win,Ukin,Vkin,Wkin,Wave,nbcpus,spec_rank)
!     call aprioriallEQDHD(Uin,Vin,Win,Ukin,Vkin,Wkin,&
!                     & Wave,spec_rank,nbcpus)
    if(spec_rank.EQ.0) then

     if (timesimu .LT. 0.00000001)then

       open(10,file='mean_Ek_GS.out',form='formatted')
       do i=1,mhdndelta
       write(10,*) timesimu, meanEkGSsumT(1,i),meanEkGSsumT(2,i)
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
  
       open(10,file='meanDiffUSGSEqdm.out',form='formatted')
        do i=1,mhdndelta                                       
         write(10,*) timesimu, meanDiffUSGSEqdmsumT(1,i), meanDiffUSGSEqdmsumT(2,i)
        enddo                                                
       close(10)

       open(10,file='meanDissUSGSEqdm.out',form='formatted')
        do i=1,mhdndelta
         write(10,*) timesimu, meanDissUSGSEqdmsumT(1,i), meanDissUSGSEqdmsumT(2,i)
        enddo
       close(10)

      
!   
!        open(10,file='meanEqdmTUuuSmag.out',form='formatted')
!         do i=1,mhdndelta
!          write(10,*)  timesimu,meanEqdmTUuumodelSmag(1,i), meanEqdmTUuumodelSmag(2,i)
!         enddo
!        close(10)

       open(10,file='meanDissSGSEQDMsumT.out',form='formatted')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanDissSGSEQDMsumT(1,i), meanDissSGSEQDMsumT(2,i)
        enddo
       close(10)

     else

       open(10,file='mean_Ek_GS.out',form='formatted',position='append')
       do i=1,mhdndelta
       write(10,*) timesimu, meanEkGSsumT(1,i),meanEkGSsumT(2,i)
       enddo
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


       open(10,file='meanDiffUSGSEqdm.out',form='formatted',position='append')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanDiffUSGSEqdmsumT(1,i), meanDiffUSGSEqdmsumT(2,i)
        enddo
       close(10)

       open(10,file='meanDissUSGSEqdm.out',form='formatted',position='append')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanDissUSGSEqdmsumT(1,i), meanDissUSGSEqdmsumT(2,i)
        enddo
       close(10)

 
       open(10,file='meanDissSGSEQDMsumT.out',form='formatted',position='append')
        do i=1,mhdndelta
         write(10,*)  timesimu,meanDissSGSEQDMsumT(1,i), meanDissSGSEQDMsumT(2,i)
        enddo
       close(10)
! 
!        open(10,file='meanEqdmTUuuSmag.out',form='formatted',position='append')
!         do i=1,mhdndelta
!          write(10,*)  timesimu,meanEqdmTUuumodelSmag(1,i), meanEqdmTUuumodelSmag(2,i)
!         enddo
!        close(10)

    endif
endif


end subroutine computeallthetermsHD


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


subroutine computeallGSandSGSHD(Uin,Vin,Win,Ukin,Vkin,Wkin,Wave,nbcpus,spec_rank)

  implicit none
  TYPE(REAL_DATA_LAYOUT), INTENT(IN)          :: Uin,Vin,Win
  TYPE(REAL_DATA_LAYOUT) :: Uf,Vf,Wf
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Ukin,Vkin,Wkin
  TYPE(COMPLEX_DATA_LAYOUT) :: Ukf,Vkf,Wkf

  INTEGER :: nbcpus, spec_rank, ndelta
  TYPE(WAVENUMBERS) :: Wave
  REAL(WP) :: delta
  REAL(WP) :: meanDissGSEQDM, meanAdvecGSEQDM, meanDifGSEQDM
  REAL(WP) :: meanDiffUSGSEqdm, meanDissUSGSEqdm, meanDissSGSEQDM, meanEkGS
  logical :: res
  IF ((.NOT. copyStructOnly(Uin,Uf)) .OR. &
     &(.NOT. copyStructOnly(Uin,Vf)) .OR. &
     &(.NOT. copyStructOnly(Uin,Wf))  ) RETURN

  IF ((.NOT. copyStructOnly(Ukin,Ukf)) .OR. &
     &(.NOT. copyStructOnly(Ukin,Vkf)) .OR. &
     &(.NOT. copyStructOnly(Ukin,Wkf))  ) RETURN

    DO ndelta=1,mhdndelta
    delta = Uin%Lx*(2*ndelta)/Uin%nx

    meanDifGSEQDMsumT(1,ndelta)    = 2*ndelta
    meanAdvecGSEQDMsumT(1,ndelta)  = 2*ndelta
    meanDissGSEQDMsumT(1,ndelta)   = 2*ndelta
    meanDiffUSGSEqdmsumT(1,ndelta) = 2*ndelta
    meanDissUSGSEqdmsumT(1,ndelta) = 2*ndelta
    meanDissSGSEQDMsumT(1,ndelta)  = 2*ndelta

    meanEkGSsumT(1,ndelta) = 2*ndelta

    CALL computePhysicalFilter(Wave,Uin,Uf,delta,filterType)
    CALL computePhysicalFilter(Wave,Vin,Vf,delta,filterType)
    CALL computePhysicalFilter(Wave,Win,Wf,delta,filterType)

   CALL ftran(Uf,Ukf,res)
   IF (.NOT.res) RETURN

   CALL ftran(Vf,Vkf,res)
   IF (.NOT.res) RETURN

   CALL ftran(Wf,Wkf,res)
   IF (.NOT.res) RETURN

    CALL EQDMGSHD(Uf,Vf,Wf,Uin,Vin,Win,Ukf,Vkf,Wkf,Ukin,Vkin,Wkin,Wave,delta,ndelta,filterType,nbcpus,spec_rank,mu,&
                   & meanDissGSEQDM, meanAdvecGSEQDM, meanDifGSEQDM,&
                   & meanDiffUSGSEqdm, meanDissUSGSEqdm, meanDissSGSEQDM, meanEkGS)

    if(timesimu.EQ.0.0)  then

     meanDifGSEQDMsumT(2,ndelta)    = meanDifGSEQDM
     meanAdvecGSEQDMsumT(2,ndelta)  = meanAdvecGSEQDM
     meanDissGSEQDMsumT(2,ndelta)   = meanDissGSEQDM
     meanDiffUSGSEqdmsumT(2,ndelta) = meanDiffUSGSEqdm
     meanDissUSGSEqdmsumT(2,ndelta) = meanDissUSGSEqdm
     meanDissSGSEQDMsumT(2,ndelta)   = meanDissSGSEQDM
     meanEkGSsumT(2,ndelta)   = meanEkGS
    else

     meanDifGSEQDMsumT(2,ndelta)    = ( (timesimu-dt)*meanDifGSEQDMsumT(2,ndelta) + dt*meanDifGSEQDM )/timesimu
     meanAdvecGSEQDMsumT(2,ndelta)  = ( (timesimu-dt)*meanAdvecGSEQDMsumT(2,ndelta) + dt*meanAdvecGSEQDM )/timesimu
     meanDissGSEQDMsumT(2,ndelta)   = ( (timesimu-dt)*meanDissGSEQDMsumT(2,ndelta) + dt*meanDissGSEQDM )/timesimu
     meanDiffUSGSEqdmsumT(2,ndelta) = ( (timesimu-dt)*meanDiffUSGSEqdmsumT(2,ndelta) + dt*meanDiffUSGSEqdm )/timesimu
     meanDissUSGSEqdmsumT(2,ndelta) = ( (timesimu-dt)*meanDissUSGSEqdmsumT(2,ndelta) + dt*meanDissUSGSEqdm )/timesimu
     meanDissSGSEQDMsumT(2,ndelta)  = ( (timesimu-dt)*meanDissSGSEQDMsumT(2,ndelta) + dt*meanDissSGSEQDM )/timesimu
     meanEkGSsumT(2,ndelta)  = ( (timesimu-dt)*meanEkGSsumT(2,ndelta) + dt*meanEkGS )/timesimu

    endif

    ENDDO

    CALL deleteDataLayout(Uf)
    CALL deleteDataLayout(Vf)
    CALL deleteDataLayout(Wf)

    CALL deleteDataLayout(Ukf)
    CALL deleteDataLayout(Vkf)
    CALL deleteDataLayout(Wkf)
end subroutine computeallGSandSGSHD
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



subroutine aprioriallEQDHD(Uin,Vin,Win,Ukin,Vkin,Wkin,Wave,spec_rank,nbcpus)

  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Uin,Vin,Win
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Ukin,Vkin,Wkin
  TYPE(WaveNumbers), INTENT(IN) :: Wave

  TYPE(REAL_DATA_LAYOUT) :: Uf,Vf,Wf
  TYPE(COMPLEX_DATA_LAYOUT) :: Ukf,Vkf,Wkf

  INTEGER :: spec_rank,nbcpus,type_avg
  logical :: res
  REAL(WP) :: delta
  Integer :: ndelta

  IF ((.NOT. copyStructOnly(Ukin,Ukf)) .OR. &
     &(.NOT. copyStructOnly(Ukin,Vkf)) .OR. &
     &(.NOT. copyStructOnly(Ukin,Wkf)) .OR. &
     &(.NOT. copyStructOnly(Uin,Uf)) .OR. &
     &(.NOT. copyStructOnly(Uin,Vf)) .OR. &
     &(.NOT. copyStructOnly(Uin,Wf))    ) RETURN 

    DO ndelta=1,mhdndelta
    delta = Uin%Lx*(2*ndelta)/Uin%nx
type_avg=0
    meanEqdmTUuumodelSmag(1,ndelta) = 2*ndelta

     CALL computeFilter(wave,delta,Ukin,Ukf,filterType) 
     CALL btran(Ukf,Uf,res)
     CALL computeFilter(wave,delta,Vkin,Vkf,filterType) 
     CALL btran(Vkf,Vf,res)
     CALL computeFilter(wave,delta,Wkin,Wkf,filterType) 
     CALL btran(Wkf,Wf,res)

    call  aprioriEQDHD(Uf,Vf,Wf,Ukf,Vkf,Wkf,&
                    & Wave,spec_rank,delta,meanEqdmTUuumodelSmag(2,ndelta),nbcpus,filterType,type_avg)

    enddo

    CALL deleteDataLayout(Ukf)
    CALL deleteDataLayout(Vkf)
    CALL deleteDataLayout(Wkf)
    CALL deleteDataLayout(Uf)
    CALL deleteDataLayout(Vf)
    CALL deleteDataLayout(Wf)


end subroutine aprioriallEQDHD


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



end module post_hd
!> @}
