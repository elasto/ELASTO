!> @addtogroup post_process
!! @{

!------------------------------------------------------------------------------

!
! MODULE: post_hd_out
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

module post_hd_out

    use datalayout
    use toolbox 
    use stat_tools 
    use differential_tools 
    USE wavenumber_tools 
    use transforms_tools
    use subgrid_tools
    use differential_tools 
    use filtering_tools 
    use velmodels
    use post_hd_ener_lib
    use post_hd_kinHel_lib
    use data
    implicit none
    private

     ! ==== Public procedures ==== 
    public          :: computeallthetermsHDout
    public          :: post_HD_init
    ! ======= Mean terms ========
!     real(WP), PROTECTED, SAVE :: meanadvecQDMsumT, meandiffvisQDMsumT, meandissvisQDMsumT, sommesumT 
!     real(WP), PROTECTED, SAVE :: meanEksumT,meanEk

!     real(WP), PROTECTED, ALLOCATABLE, dimension(:,:), SAVE :: meanEqdmTUuumodelSmag


    INTEGER, PROTECTED, SAVE :: hdndelta,filterType

contains

subroutine post_HD_init(Uin,Vin,Win,&
                         & nbcpus,spec_rank)

    implicit none

    TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Uin,Vin,Win
    INTEGER :: nbcpus,spec_rank
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

end subroutine post_HD_init

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
subroutine computeallthetermsHDout(Uin,Vin,Win,&
                         & nbcpus,spec_rank)

    implicit none

    TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Uin,Vin,Win
    INTEGER :: nbcpus,spec_rank
    real(WP):: mu
    logical :: enerpostHD,HelpostHD

  CALL parser_read('HD filter size',hdndelta)
  CALL parser_read('HD filter type',filterType)
  CALL parser_read('HD Ener Post',enerpostHD)
  CALL parser_read('HD Hel Post',HelpostHD)
  CALL parser_read('Viscosity',mu)


  call post_HD_init(Uin,Vin,Win,&
                         & nbcpus,spec_rank)
  if (enerpostHD) then
   call computeallEnertermsHDout(Uin,Vin,Win,Uk,Vk,Wk,&
                         & VelWN,nbcpus,spec_rank,mu)
  endif

  if (HelpostHD) then
   call computeallHkintermsHDout(Uin,Vin,Win,Uk,Vk,Wk,&
                         & VelWN,nbcpus,spec_rank,mu)
  endif
 if (spec_rank.EQ.0) print*, "mu", mu

end subroutine computeallthetermsHDout

! ===== Compute all the Energy related terms =====

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
subroutine computeallEnertermsHDout(Uin,Vin,Win,Ukin,Vkin,Wkin,&
                         & Wave,nbcpus,spec_rank,muin)

    implicit none

    TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Uin,Vin,Win
    TYPE(WaveNumbers),INTENT(IN) :: Wave
    TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Ukin,Vkin,Wkin
    real(WP):: meanadvecQDM, meandiffvisQDM, meandissvisQDM, somme,meanEk
    INTEGER :: nbcpus,spec_rank
    real(WP) ::muin

    call eqdm(Uin,Vin,Win,Ukin,Vkin,Wkin,Wave,nbcpus,spec_rank,muin,&
             &meanadvecQDM, meandiffvisQDM, meandissvisQDM, somme,meanEk)

    if(spec_rank.EQ.0) then

       open(10,file='mean_values_QDM.out',form='formatted')
       write(10,*) meanadvecQDM, meandiffvisQDM, meandissvisQDM,&
                 & somme,meanEk
       close(10)
 
    endif

    call computeGSandSGSEnerHD(Uin,Vin,Win,Ukin,Vkin,Wkin,Wave,nbcpus,spec_rank,muin)
!     call aprioriallEQDHD(Uin,Vin,Win,Ukin,Vkin,Wkin,&
!                     & Wave,spec_rank,nbcpus)
 


end subroutine computeallEnertermsHDout


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


subroutine computeGSandSGSEnerHD(Uin,Vin,Win,Ukin,Vkin,Wkin,Wave,nbcpus,spec_rank,muin)

  implicit none
  TYPE(REAL_DATA_LAYOUT), INTENT(IN)          :: Uin,Vin,Win
  TYPE(REAL_DATA_LAYOUT) :: Uf,Vf,Wf
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Ukin,Vkin,Wkin
  TYPE(COMPLEX_DATA_LAYOUT) :: Ukf,Vkf,Wkf

  INTEGER :: nbcpus, spec_rank, ndelta,i
  TYPE(WAVENUMBERS) :: Wave
  REAL(WP) :: delta
  REAL(WP) :: muin
  logical :: res
    real(WP), ALLOCATABLE, dimension(:,:) :: meanDissGSEQDMsumT, meanAdvecGSEQDMsumT, &
                                                            & meanDifGSEQDMsumT
    real(WP), ALLOCATABLE, dimension(:,:) :: meanDiffUSGSEqdmsumT, meanDissUSGSEqdmsumT, &
                                                            & meanDissSGSEQDMsumT
    real(WP), ALLOCATABLE, dimension(:,:) :: meanEkGSsumT

  allocate( meanDifGSEQDMsumT(2,hdndelta) )    
  allocate( meanAdvecGSEQDMsumT(2,hdndelta) )  
  allocate( meanDissGSEQDMsumT(2,hdndelta) )  
  allocate( meanDiffUSGSEqdmsumT(2,hdndelta) )
  allocate( meanDissUSGSEqdmsumT(2,hdndelta) )
  allocate( meanDissSGSEQDMsumT(2,hdndelta) )

  allocate( meanEkGSsumT(2,hdndelta) )

  IF ((.NOT. copyStructOnly(Uin,Uf)) .OR. &
     &(.NOT. copyStructOnly(Uin,Vf)) .OR. &
     &(.NOT. copyStructOnly(Uin,Wf))  ) RETURN

  IF ((.NOT. copyStructOnly(Ukin,Ukf)) .OR. &
     &(.NOT. copyStructOnly(Ukin,Vkf)) .OR. &
     &(.NOT. copyStructOnly(Ukin,Wkf))  ) RETURN

    DO ndelta=1,hdndelta
    delta = Uin%Lx*(2*ndelta)/Uin%nx

    meanDifGSEQDMsumT(1,ndelta)    = 2*ndelta
    meanAdvecGSEQDMsumT(1,ndelta)  = 2*ndelta
    meanDissGSEQDMsumT(1,ndelta)   = 2*ndelta
    meanDiffUSGSEqdmsumT(1,ndelta) = 2*ndelta
    meanDissUSGSEqdmsumT(1,ndelta) = 2*ndelta
    meanDissSGSEQDMsumT(1,ndelta)  = 2*ndelta

    meanEkGSsumT(1,ndelta) = 2*ndelta

   call computeFilter(Wave, delta, Ukin, Ukf, filterType)
   CALL btran(Ukf,Uf,res)
   IF (.NOT.res) RETURN

   call computeFilter(Wave, delta, Vkin, Vkf, filterType)
   CALL btran(Vkf,Vf,res)
   IF (.NOT.res) RETURN

   call computeFilter(Wave, delta, Wkin, Wkf, filterType)
   CALL btran(Wkf,Wf,res)
   IF (.NOT.res) RETURN

    CALL EQDMGSHD(Uf,Vf,Wf,Uin,Vin,Win,Ukf,Vkf,Wkf,Ukin,Vkin,Wkin,Wave,delta,ndelta,filterType,nbcpus,spec_rank,muin,&
                   & meanDissGSEQDMsumT(2,ndelta), meanAdvecGSEQDMsumT(2,ndelta),meanDifGSEQDMsumT(2,ndelta) ,&
                   & meanDiffUSGSEqdmsumT(2,ndelta), meanDissUSGSEqdmsumT(2,ndelta), meanDissSGSEQDMsumT(2,ndelta), meanEkGSsumT(2,ndelta))

    ENDDO

    CALL deleteDataLayout(Uf)
    CALL deleteDataLayout(Vf)
    CALL deleteDataLayout(Wf)

    CALL deleteDataLayout(Ukf)
    CALL deleteDataLayout(Vkf)
    CALL deleteDataLayout(Wkf)

   if(spec_rank.EQ.0) then


       open(10,file='mean_Ek_GS.out',form='formatted')
       do i=1,hdndelta
       write(10,*) meanEkGSsumT(1,i),meanEkGSsumT(2,i)
       enddo
       close(10)

       open(10,file='meanDissGSEQDM.out',form='formatted')
        do i=1,hdndelta
          write(10,*) meanDissGSEQDMsumT(1,i),meanDissGSEQDMsumT(2,i)
        enddo
       close(10)

       open(10,file='meanAdvecGSEQDM.out',form='formatted')
        do i=1,hdndelta
         write(10,*) meanAdvecGSEQDMsumT(1,i), meanAdvecGSEQDMsumT(2,i)
        enddo
       close(10)

       open(10,file='meanDifGSEQDM.out',form='formatted')
        do i=1,hdndelta
         write(10,*) meanDifGSEQDMsumT(1,i), meanDifGSEQDMsumT(2,i)
        enddo
       close(10)     
  
       open(10,file='meanDiffUSGSEqdm.out',form='formatted')
        do i=1,hdndelta                                       
         write(10,*) meanDiffUSGSEqdmsumT(1,i), meanDiffUSGSEqdmsumT(2,i)
        enddo                                                
       close(10)

       open(10,file='meanDissUSGSEqdm.out',form='formatted')
        do i=1,hdndelta
         write(10,*) meanDissUSGSEqdmsumT(1,i), meanDissUSGSEqdmsumT(2,i)
        enddo
       close(10)
      
       open(10,file='meanDissSGSEQDMsumT.out',form='formatted')
        do i=1,hdndelta
         write(10,*)  meanDissSGSEQDMsumT(1,i), meanDissSGSEQDMsumT(2,i)
        enddo
       close(10)


endif
  deallocate( meanDifGSEQDMsumT )    
  deallocate( meanAdvecGSEQDMsumT )  
  deallocate( meanDissGSEQDMsumT )  
  deallocate( meanDiffUSGSEqdmsumT )
  deallocate( meanDissUSGSEqdmsumT )
  deallocate( meanDissSGSEQDMsumT )

  deallocate( meanEkGSsumT )

end subroutine computeGSandSGSEnerHD


! ===== Compute all the kinetic helicity related terms =====

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
subroutine computeallHkintermsHDout(Uin,Vin,Win,Ukin,Vkin,Wkin,&
                         & Wave,nbcpus,spec_rank,muin)

    implicit none

    TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Uin,Vin,Win
    TYPE(WaveNumbers),INTENT(IN) :: Wave
    TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Ukin,Vkin,Wkin
    TYPE(REAL_DATA_LAYOUT) :: Vortx,Vorty,Vortz
    TYPE(COMPLEX_DATA_LAYOUT) :: Vortxk,Vortyk,Vortzk
    real(WP) :: meanadvecHkin, meandiffvisHkin, meandissvisHkin, meanHkin
    INTEGER  :: nbcpus,spec_rank
    real(WP) :: muin
    logical :: res
     IF ( (.NOT. copyStructOnly(Ukin,Vortxk ))) RETURN
     IF ( (.NOT. copyStructOnly(Ukin,Vortyk ))) RETURN
     IF ( (.NOT. copyStructOnly(Ukin,Vortzk ))) RETURN

     IF ( (.NOT. copyStructOnly(Uin,Vortx ))) RETURN
     IF ( (.NOT. copyStructOnly(Uin,Vorty ))) RETURN
     IF ( (.NOT. copyStructOnly(Uin,Vortz ))) RETURN



    call computeVorticity(Ukin,Vkin,Wkin,wave,Vortx,Vorty,Vortz,res)

           CALL ftran(Vortx,Vortxk,res)
           IF (.NOT.res) RETURN
           CALL ftran(Vorty,Vortyk,res)
           IF (.NOT.res) RETURN
           CALL ftran(Vortz,Vortzk,res)
           IF (.NOT.res) RETURN

    call EqHkin(Uin,Vin,Win,Ukin,VKin,WKin,Wave,Vortx,Vorty,Vortz,Vortxk,Vortyk,Vortzk,&
                 &nbcpus,spec_rank,muin,meanadvecHkin, meandiffvisHkin, meandissvisHkin, meanHkin)

    if(spec_rank.EQ.0) then

       open(10,file='mean_values_EqHkin.out',form='formatted')
       write(10,*) meanadvecHkin, meandiffvisHkin, meandissvisHkin, meanHkin
       close(10)
 
    endif

  call computeGSandSGSHkinHD(Uin,Vin,Win,Ukin,Vkin,Wkin,Vortxk,Vortyk,Vortzk,&
                                &Vortx,Vorty,Vortz,Wave,nbcpus,spec_rank,muin)
 
  CALL deleteDataLayout(Vortxk)
  CALL deleteDataLayout(Vortyk)
  CALL deleteDataLayout(Vortzk)

  CALL deleteDataLayout(Vortx)
  CALL deleteDataLayout(Vorty)
  CALL deleteDataLayout(Vortz)

end subroutine computeallHkintermsHDout


subroutine computeGSandSGSHkinHD(Uin,Vin,Win,Ukin,Vkin,Wkin,Vortxk,Vortyk,Vortzk,&
                                &Vortx,Vorty,Vortz,Wave,nbcpus,spec_rank,muin)

  implicit none
  TYPE(REAL_DATA_LAYOUT), INTENT(IN)          :: Uin,Vin,Win
  TYPE(REAL_DATA_LAYOUT) :: Uf,Vf,Wf
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Ukin,Vkin,Wkin
  TYPE(COMPLEX_DATA_LAYOUT) :: Ukf,Vkf,Wkf

  TYPE(REAL_DATA_LAYOUT) :: Vortx,Vorty,Vortz
  TYPE(COMPLEX_DATA_LAYOUT) :: Vortxk,Vortyk,Vortzk

  TYPE(REAL_DATA_LAYOUT) :: Vortxf,Vortyf,Vortzf
  TYPE(COMPLEX_DATA_LAYOUT) :: Vortxkf,Vortykf,Vortzkf

  INTEGER :: nbcpus, spec_rank, ndelta,i
  TYPE(WAVENUMBERS) :: Wave
  REAL(WP) :: delta
  REAL(WP) :: muin
  logical :: res
    real(WP), ALLOCATABLE, dimension(:,:) :: meanDissGSEqHkinSumT, meanAdvecGSEqHkinSumT, meanDifGSEqHkinSumT,&
                & meanDissTUWUijSGSHkinSumT, meanDiffTUWSGSHkinSumT,meanDiffTUSGSHkinSumT,&
                & meanDissTUSwSGSHkinSumT, meanHkinGSSumT

  allocate( meanDissGSEqHkinSumT(2,hdndelta) )    
  allocate( meanAdvecGSEqHkinSumT(2,hdndelta) )  
  allocate( meanDifGSEqHkinSumT(2,hdndelta) )  
  allocate( meanDissTUWUijSGSHkinSumT(2,hdndelta) )
  allocate( meanDiffTUWSGSHkinSumT(2,hdndelta) )
  allocate( meanDiffTUSGSHkinSumT(2,hdndelta) )
  allocate( meanDissTUSwSGSHkinSumT(2,hdndelta) )
  allocate( meanHkinGSSumT(2,hdndelta) )

  IF ((.NOT. copyStructOnly(Uin,Uf)) .OR. &
     &(.NOT. copyStructOnly(Uin,Vf)) .OR. &
     &(.NOT. copyStructOnly(Uin,Wf))  ) RETURN

  IF ((.NOT. copyStructOnly(Ukin,Ukf)) .OR. &
     &(.NOT. copyStructOnly(Ukin,Vkf)) .OR. &
     &(.NOT. copyStructOnly(Ukin,Wkf))  ) RETURN

     IF ( (.NOT. copyStructOnly(Ukin,Vortxkf ))) RETURN
     IF ( (.NOT. copyStructOnly(Ukin,Vortykf ))) RETURN
     IF ( (.NOT. copyStructOnly(Ukin,Vortzkf ))) RETURN

     IF ( (.NOT. copyStructOnly(Uin,Vortxf ))) RETURN
     IF ( (.NOT. copyStructOnly(Uin,Vortyf ))) RETURN
     IF ( (.NOT. copyStructOnly(Uin,Vortzf ))) RETURN


    DO ndelta=1,hdndelta
    delta = Uin%Lx*(2*ndelta)/Uin%nx


    meanDissGSEqHkinSumT(1,ndelta)      = 2*ndelta
    meanAdvecGSEqHkinSumT(1,ndelta)     = 2*ndelta 
    meanDifGSEqHkinSumT(1,ndelta)       = 2*ndelta 
    meanDissTUWUijSGSHkinSumT(1,ndelta) = 2*ndelta
    meanDiffTUWSGSHkinSumT(1,ndelta)    = 2*ndelta
    meanDiffTUSGSHkinSumT(1,ndelta)     = 2*ndelta
    meanDissTUSwSGSHkinSumT(1,ndelta)   = 2*ndelta
    meanHkinGSSumT(1,ndelta)            = 2*ndelta

   call computeFilter(Wave, delta, Ukin, Ukf, filterType)
   CALL btran(Ukf,Uf,res)
   IF (.NOT.res) RETURN

   call computeFilter(Wave, delta, Vkin, Vkf, filterType)
   CALL btran(Vkf,Vf,res)
   IF (.NOT.res) RETURN

   call computeFilter(Wave, delta, Wkin, Wkf, filterType)
   CALL btran(Wkf,Wf,res)
   IF (.NOT.res) RETURN

   call computeFilter(Wave, delta, Vortxk, Vortxkf, filterType)
   CALL btran(Vortxkf,Vortxf,res)
   IF (.NOT.res) RETURN

   call computeFilter(Wave, delta, Vortyk, Vortykf, filterType)
   CALL btran(Vortykf,Vortyf,res)
   IF (.NOT.res) RETURN

   call computeFilter(Wave, delta, Vortzk, Vortzkf, filterType)
   CALL btran(Vortzkf,Vortzf,res)
   IF (.NOT.res) RETURN

    CALL EqHkinGS(Uf,Vf,Wf,Uin,Vin,Win,Ukf,Vkf,Wkf,Ukin,VKin,WKin,&
                & Vortx,Vorty,Vortz,Vortxf,Vortyf,Vortzf,&
                & Vortxk,Vortyk,Vortzk,Vortxkf,Vortykf,Vortzkf,&
                & Wave,delta,filterType,nbcpus,spec_rank,muin,&
                & meanDissGSEqHkinSumT(2,ndelta), meanAdvecGSEqHkinSumT(2,ndelta), meanDifGSEqHkinSumT(2,ndelta),&
                & meanDissTUWUijSGSHkinSumT(2,ndelta), meanDiffTUWSGSHkinSumT(2,ndelta),meanDiffTUSGSHkinSumT(2,ndelta),&
                & meanDissTUSwSGSHkinSumT(2,ndelta), meanHkinGSSumT(2,ndelta))

    ENDDO

    CALL deleteDataLayout(Uf)
    CALL deleteDataLayout(Vf)
    CALL deleteDataLayout(Wf)

    CALL deleteDataLayout(Ukf)
    CALL deleteDataLayout(Vkf)
    CALL deleteDataLayout(Wkf)

  CALL deleteDataLayout(Vortxkf)
  CALL deleteDataLayout(Vortykf)
  CALL deleteDataLayout(Vortzkf)

  CALL deleteDataLayout(Vortxf)
  CALL deleteDataLayout(Vortyf)
  CALL deleteDataLayout(Vortzf)

   if(spec_rank.EQ.0) then

       open(10,file='mean_Hk_GS.out',form='formatted')
       do i=1,hdndelta
       write(10,*) meanHkinGSSumT(1,i),meanHkinGSSumT(2,i)
       enddo
       close(10)

       open(10,file='meanDissGSEqHkin.out',form='formatted')
       do i=1,hdndelta
       write(10,*) meanDissGSEqHkinSumT(1,i),meanDissGSEqHkinSumT(2,i)
       enddo
       close(10)

       open(10,file='meanAdvecGSEqHkin.out',form='formatted')
       do i=1,hdndelta
       write(10,*) meanAdvecGSEqHkinSumT(1,i),meanAdvecGSEqHkinSumT(2,i)
       enddo
       close(10)


       open(10,file='meanDifGSEqHkin.out',form='formatted')
       do i=1,hdndelta
       write(10,*) meanDifGSEqHkinSumT(1,i),meanDifGSEqHkinSumT(2,i)
       enddo
       close(10)

       open(10,file='meanDissTUWUijSGSHkin.out',form='formatted')
       do i=1,hdndelta
       write(10,*) meanDissTUWUijSGSHkinSumT(1,i),meanDissTUWUijSGSHkinSumT(2,i)
       enddo
       close(10)

       open(10,file='meanDiffTUWSGSHkin.out',form='formatted')
       do i=1,hdndelta
       write(10,*) meanDiffTUWSGSHkinSumT(1,i),meanDiffTUWSGSHkinSumT(2,i)
       enddo
       close(10)

       open(10,file='meanDiffTUSGSHkin.out',form='formatted')
       do i=1,hdndelta
       write(10,*) meanDiffTUSGSHkinSumT(1,i),meanDiffTUSGSHkinSumT(2,i)
       enddo
       close(10)

       open(10,file='meanDissTUSwSGSHkin.out',form='formatted')
       do i=1,hdndelta
       write(10,*) meanDissTUSwSGSHkinSumT(1,i),meanDissTUSwSGSHkinSumT(2,i)
       enddo
       close(10)


  endif

  deallocate( meanDissGSEqHkinSumT )    
  deallocate( meanAdvecGSEqHkinSumT )  
  deallocate( meanDifGSEqHkinSumT )  
  deallocate( meanDissTUWUijSGSHkinSumT )
  deallocate( meanDiffTUWSGSHkinSumT )
  deallocate( meanDiffTUSGSHkinSumT )
  deallocate( meanDissTUSwSGSHkinSumT )
  deallocate( meanHkinGSSumT )

end subroutine computeGSandSGSHkinHD

end module post_hd_out
!> @}
