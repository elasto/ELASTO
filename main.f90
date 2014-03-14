!===============================================================================
PROGRAM main
!==============================================================================
  USE options
  USE param
  USE mpilayout_tools
  USE communicators_tools
  USE random_tools
  USE data
  USE datalayout
  USE precision_tools
  USE avs
  USE mpi
  USE solver
  USE post_lib
  USE post_in_simu
  USE post_out_simu
  USE stat_tools
 !USE subgridmodels

  USE wavenumber_tools
  USE io_interface    

  IMPLICIT NONE
  
  INTEGER :: nbcpus=1
  !> the number of MPI processes
  INTEGER :: ncpus1=1
  !> the number of MPI processes in the first distributed array dimension
  INTEGER :: ncpus2=1
  !> the number of MPI processes in the second distributed array dimension
  INTEGER :: spec_rank
  !> my id in the pool of MPI processes inside the "spectral" communicator
  INTEGER :: world_rank
  !> my id in the pool of MPI processes inside the default "MPI_COMM_WOLRD" communicator
  INTEGER :: ierr
  !> Dimension of the mpi-topology used (if 0 the none)
  INTEGER :: topo_dim
  !> an integer to store error codes returned by subroutines
  LOGICAL :: res
  !> To determine postprocessing to perform (none => computation)
  CHARACTER(LEN=32):: postprocess
  !> Reached time (by the simulation)
  REAL(WP)    :: current_time=0.0_WP
  ! Time step
  REAL(WP)    :: dt, med

  !LUCA
  REAL(WP)    :: ps,T,gammadot_min,corrXYw
  INTEGER     :: Nw,freq_sigma_corr

  !VAR DEFINED FOR TESTING
  TYPE(REAL_DATA_LAYOUT)      :: sigma, sigma_old,mult
  TYPE(REAL_DATA_LAYOUT)      :: Esigma,Lsigma,n_activ,convol
  TYPE(COMPLEX_DATA_LAYOUT)   :: G,sigmak
  LOGICAL                     :: success

  !> number of solver iterations requested
  INTEGER::n_iter
  !> perform post-process ipost timestep
  INTEGER:: ipost
  INTEGER:: ite, sca
  REAL(WP) :: coeff, final_time
  REAL:: time1, time2,times=0.0,timesPost=0.0,tx

  ! LUCA
  INTEGER:: i,j,k

  !== INITIALIZE MPI ENVIRONMENT ==
  CALL MPI_INIT(ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,world_rank,ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nbcpus,ierr)

  ! Gets the command line options.
  IF (.NOT. get_options(world_rank)) THEN
    CALL MPI_FINALIZE(ierr)
    STOP 'FAILED'
  ENDIF

  CALL computeProcsMatrix(nbcpus,ncpus1,ncpus2,world_rank,res)
  IF (.NOT. res) THEN
    CALL MPI_FINALIZE(ierr)
    STOP 'computeProcsMatrix FAILED'
  ELSE
    IF (world_rank.EQ.0) WRITE(6,'(2(a,i3),a)') &
      & '[INFO] Cpus distribution is [',ncpus1,'] x [',ncpus2,']'
  ENDIF

  !== READS THE INPUT FILE CALLED "input" AND INITIALIZE PARSER ==
  IF (.NOT. param_init(GetInputFile())) THEN
    CALL MPI_FINALIZE(ierr)
    STOP 'param_init FAILED'
  ENDIF
  IF (world_rank.EQ.0) CALL parser_show()

  !== INITIALIZE PARALLEL CONTEXT ==
  !ATTENTION TO THE TOPOLOGY
  ! Read input to initialize mpi topology:
    CALL parser_read('Default topology', topo_dim)
!    topo_dim = 3

  ! Initialize the mpi topology and the communicators
  CALL initCommunicators(ncpus1, ncpus2, topo_dim, spec_rank)
  ! Send new rank to mpilayout_tools
  CALL init_spec_rank(spec_rank)

  !== INITIALIZE THE COMPUTATION PART OF THE APPLICATION ==
  ! Check if we have asked for the debug option with GetDebugRandom()
  !------------------------------------------------
  CALL random_init(spec_rank,GetDebugRandom())

  !== INITIALIZE THE DATA STRUCTRE
  IF(.NOT.data_init_luca(nbcpus,spec_rank,mult,sigma,sigma_old,Esigma&
          & ,Lsigma,G,sigmak,n_activ,convol)) & 
&     STOP 'FAILED to initialize the data structure'


  !------------------------------------------------
  ! Perform post out simu OR run a simulation
  !------------------------------------------------
  CALL parser_read('Post-process out simulation',postprocess)
  IF ( postprocess.eq.'yes' ) THEN
  !------------------------------------------------
  ! Post out simu part
  !------------------------------------------------
    IF(spec_rank.EQ.0) WRITE(6,'(a)')'[INFO] Post-process out of simulation'
    CALL cpu_time(time1)

!    IF (.NOT. post_out_simu_init_luca(spec_rank,sigma)) THEN
!      CALL MPI_FINALIZE(ierr)
!      STOP 'Error in initialization of post-processing'
!    ENDIF

!    IF (.NOT. post_out(nbcpus ,spec_rank, nbscal, scalArray, nbscal_part,scal_partArray, U, V ,W,B,imhd)) THEN
!      CALL MPI_FINALIZE(ierr)
!      STOP 'Error in post-processing'
!    ENDIF

    CALL cpu_time(time2)
    IF(spec_rank.EQ.0) WRITE(6,'(a,x,g16.6,a)')'[INFO] Postprocess requires ',time2-time1,' s'

!    IF (.NOT. post_out_simu_del(nbscal,nbscal_part) ) THEN
!      CALL MPI_FINALIZE(ierr)
!      STOP 'Error in delete of post-processing'
!    ENDIF

  ELSE ! ( postprocess.eq.'no' )

  !=== COMPUTATIONAL PART ===
  
    !=== INITIALIZATION OF THE SOLVER ===
    IF(spec_rank.EQ.0) PRINT*,'-----------------------------------------------------------------'

    !== INITIALIZE GAMMA DOT ==
    CALL parser_read('Gamma dot',gammadot_min)
    
    !== ALL THE FIELDS ARE INITIALIZED ==
    IF (.NOT. initSolver_luca(convol,n_activ,sigmak,G,sigma,sigma_old,Esigma,Lsigma&
              & ,gammadot_min,current_time,nbcpus,spec_rank)) THEN
        CALL MPI_FINALIZE(ierr)
        STOP 'FAILED to init solver luca'
    ENDIF

    IF(spec_rank.EQ.0) PRINT*,'-----------------------------------------------------------------'
   
     !== INITIALIZE TIME STEP dt AND FINAL TIME ==
    CALL parser_read('Simulation iterations',n_iter)
    CALL parser_read('Fixed time step',dt)
    IF (parser_is_defined('Final time')) then
      CALL parser_Read('Final time', final_time)
      coeff = 1._WP
    ELSE
      coeff = -1._WP
      final_time = 1._WP
    END IF

    !== READS TEMPERATURE ELASTIC AND PLASTIC TAU ==
    CALL parser_read('Temperature',T)
    CALL parser_read('Sigma correlation XY iter',freq_sigma_corr)
    CALL parser_read('Waiting time iter',Nw)

    !== INITIALIZE POST-PROCESS AND GETS THE OUTPUTS FIELDS FREQUENCY 
    IF(.NOT.post_in_simu_init_luca(sigma,spec_rank,ipost)) THEN
      CALL MPI_FINALIZE(ierr)
      STOP 'unable to init post-process'
    ENDIF

!    IF(.not. post_simu_luca(0,current_time, 0.0_WP,nbcpus ,spec_rank,mult,sigma,sigma_old&
!       & ,n_activ,n_iter,corrXYw)) write(*,'(a)') &
!      & '[Warning] Unable to perform post-process and to show info'
!    IF(.not. post_simu_deltaT_luca(current_time,nbcpus,spec_rank,sigma,n_iter)) write(*,'(a)') &
!      & '[Warning] Unable to perform post-process and to show info -delta T'

    ! ===== INITIALIZE RANDON NUM GEN AND ITER =====
    ite = 1

    ! ===== TIME LOOP =====
    DO WHILE ((ite<=n_iter).and.(coeff*current_time<final_time))

      !== SOLVE EQUATION ==
      CALL cpu_time(time1)
      IF (.NOT. solve_Euler_luca(convol,n_activ,mult,sigma,sigma_old,Esigma,Lsigma,gammadot_min,sigmak&
          &,G,dt,T, nbcpus,spec_rank,ite,Nw,freq_sigma_corr)) THEN
        CALL MPI_FINALIZE(ierr)
        STOP 'Solver Error'
      ENDIF

      CALL cpu_time(time2)
      IF(spec_rank.EQ.0) WRITE(6,'(a,x,i0,x,a,f8.4,a, f12.4)') &
        &'[PROGRESSION] Iteration', ite, 'requires ',time2-time1,' s, and t = ', current_time
      
      times=times+(time2-time1)

        !== OUTPUT INSIDE THE TIME LOOP ==
        CALL cpu_time(time1)
   
        !=== COMPUTE SIGMA WAITING AVERAGE ===
        IF(Nw.EQ.ite) THEN
          corrXYw = computeFieldAvg(sigma_old,spec_rank)
          print*, 'PORCO DIO:',corrXYw
        ENDIF
 
        !=== PERFORMS POST-PROCESS ===
        IF(.not. post_simu_luca(ite,current_time,dt,nbcpus,spec_rank,mult,sigma,sigma_old,&
           &n_activ,n_iter,corrXYw)) WRITE(*,'(a)') &
          & '[Warning] Unable to perform post-process and to show info'

        !== PERFEROMANCE OUTPUT ==
        CALL cpu_time(time2)
        IF(spec_rank.EQ.0) WRITE(6,'(a,x,i0,x,a,f8.4,a)') &
          & '[INFO] PostSimu at Iteration', ite,'requires ',time2-time1,' s'
        timesPost=timesPost+(time2-time1)

      IF(.not. post_simu_deltaT_luca(current_time,nbcpus,spec_rank,sigma,n_iter)) write(*,'(a)') &
        & '[Warning] Unable to perform post-process and to show info'

      ite = ite + 1

      !LUCA: the time advection step put outside the fucntion  solver_step_luca 
      current_time=current_time + dt
    ENDDO 


    !== TIME PROFILING FOR THE OUTPUT PART == 
    tx=timesPost/(timesPost+times)*100.0_WP
    IF(spec_rank.EQ.0) WRITE(6,'(a,g16.6,a)')'[INFO] Total run requires  ',times,'s'
    IF(spec_rank.EQ.0) WRITE(6,'(a,g16.6,a,f5.2,a)')'[INFO] Total post-processing requires  ',timesPost,'s soit ',tx,'%'
    IF(spec_rank.EQ.0) WRITE(6,'(a,f16.6,a)')'[INFO] Average timestep is ',times/max(1,(ite-1)),'s'

  END IF ! ( postprocess.eq.'no' )

  !== DEALLOCATE MEMORY ==
  CALL parser_reset()
  CALL deleteDataLayout(sigma)
  CALL deleteDataLayout(sigma_old)
  CALL deleteDataLayout(G)
  CALL deleteDataLayout(sigmak)
  CALL deleteDataLayout(n_activ)
  CALL deleteDataLayout(convol)
  CALL deleteDataLayout(mult)

  !=== TEMINATE ALL MPI PROCESS ===
  CALL MPI_FINALIZE(ierr)
  STOP 'SUCCESS'

END PROGRAM main


!************************ TEST PART **********************************
  !SIGMA INITIALIZATION
!  IF (.NOT. initSolver_luca(sigma,sigma_old,current_time,nbcpus,spec_rank)) THEN
!        CALL MPI_FINALIZE(ierr)
!        STOP 'FAILED to init solver luca'
!    ENDIF

!! AVERAGE TESTED and it works by pay attention read todo file
!  med  = computeFieldAvg(sigma,spec_rank)

!!TESTED
!  print*,"sigma average:", med
    !FFT TEST
!!  CALL ftran(sigma, G,success)

    !TESTED
    !MPI PRINTS THE G INITIAL CONDITION
!  IF(.NOT.cmplx_WriteToFile(spec_rank,trim(G%name),G))  THEN
!     WRITE(*,*) "does not work!!!"
!  ENDIF

!    CALL btran(G,sigma_old,res)

    !MPI PRINTS THE SIGMA INITIAL CONDITION
!  IF(.NOT.real_WriteToFile(spec_rank,trim(sigma_old%name),sigma_old))  THEN
!     WRITE(*,*) "does not work!!!"
!  ENDIF


  !TEST RANDOM_NUMBER
!!  IF(spec_rank.eq.0) THEN
!!  OPEN(400,file='random.dat',form='formatted')
!!   DO k=1,sigma%nz
!!    DO j=1,sigma%ny
!!     DO i=1,((sigma%nx)/2)+1
!!       CALL random_number(psr)
!!       write(400,'(g15.8,1x)') psr
!!     ENDDO
!!    ENDDO
!!  ENDDO
!! ENDIF

  !TESTED
  !!CALL convolution(risu,G,sigmak,res)

    !MPI PRINTS THE CONVOLUTION
  !!IF(.NOT.cmplx_WriteToFile(spec_rank,trim(risu%name),risu))  THEN
  !!   WRITE(*,*) "does not work!!!"
  !!ENDIF

    
  !HDF5 PRINT OF SIGMA INITIAL CONDITION
!!   IF(.not.write_datalayout_scalar(sigma,ite,spec_rank,"hdf5") ) THEN
!!       WRITE (6,'(a)') 'FAILED to write output for velocity field'
!!       RETURN
! !  ENDIF 

  !HDF5 PRINT OF SIGMA INITIAL CONDITION
!!   IF(.not.write_datalayout_scalar(sigma_old,ite,spec_rank) ) THEN
!!       WRITE (6,'(a)') 'FAILED to write output for velocity field'
!!       RETURN
!!   ENDIF

  !PRINT OF G PROPAGATOR INITIAL CONDITION
!! IF(spec_rank.eq.0) THEN
!!  OPEN(400,file='G-propagator.dat',form='formatted')
!!   DO k=1,sigma%nz
!!    DO j=1,sigma%ny
!!     DO i=1,((sigma%nx)/2)+1
!!      IF(G_full(i,j,k).lt.1.E-20) THEN
!!        G_full(i,j,k) = 0.d0
!!       write(400,'(g15.8,1x)') G_full(i,j,k)
!!      ELSE
!!       write(400,'(g15.8,1x)') G_full(i,j,k)
!!      ENDIF
!!    ENDDO
!!    ENDDO
!!   ENDDO

!!  CLOSE(400)
!! ENDIF
