!------------------------------------------------------------------------------
!
! MODULE: solver
!
!! DESCRIPTION:
!> This module implements the NS solver for scalars and velocities.
!> @author
!> Guillaume Balarac et Patrick BEGOU, LEGI
!>
!> @details
!! This solver implements the numerical schemes to solve diffusion, advection and NS equations.
!> @author
!> Guillaume Balarac et Patrick BEGOU, LEGI
!------------------------------------------------------------------------------
MODULE solver

    use differential_tools
    USE precision_tools
    USE datalayout
    !USE real_4D_datalayout
    USE param
    ! MY PLASTICFLOW INIT       
    USE plasticflow
    USE wavenumber_tools
    USE parallel_tools
    USE data
    USE rediscretization_tools
    USE forcing
    USE forcing_scalar
!    USE subgridmodels

    IMPLICIT NONE

    !> Storage for L.E.S. models for scalars
    INTEGER, DIMENSION(:), POINTER, PRIVATE    :: ilessc
    !> Storage for the L.E.S. model for velocities
    INTEGER, PRIVATE    :: iles
    INTEGER, PRIVATE    :: time_integration
    !> Storage of setup about forcing term on scalar i- scalars advected with pseudo-spectral solver
    integer, dimension(:), allocatable, private     :: sc_spec_forc
    !> Viscosity, Default to 1E-02
    REAL(WP), PRIVATE :: mu

! XXX TODO - Delete when new RK is validate --
    ! ===== Storage of the spectral values of field at the RK intermediate step =====
    !> Storage for the Spectral values of U at step i-1 for RK2
    TYPE(COMPLEX_DATA_LAYOUT), PRIVATE, SAVE :: UkOld,sigmakOld
    !> Storage for the Spectral values of V at step i-1 for RK2
    TYPE(COMPLEX_DATA_LAYOUT), PRIVATE, SAVE :: VkOld
    !> Storage for the Spectral values of W at step i-1 for RK2
    TYPE(COMPLEX_DATA_LAYOUT), PRIVATE, SAVE :: WkOld
    !> Storage for the spectral values of each scalar at step i-1 for RK2
    TYPE(COMPLEX_DATA_LAYOUT),ALLOCATABLE, DIMENSION(:), PRIVATE :: ScalArraykOld

    !> Storage of setup about forcing term on scalar i- scalars advected with particles method
    integer, dimension(:), allocatable, private     :: sc_part_forc
! XXX Todo
!    !> Storage for the spectral values of each scalar - scalar advected with particles method
!    TYPE(COMPLEX_VECTOR_DATA_LAYOUT), PROTECTED :: Scal_partVectorK
!    !> Storage for Schmidt numbers for scalars - scalar advected with particles method
!    REAL(WP), DIMENSION(:), POINTER, PRIVATE :: schmidt_partVector
!    !> Storage of setup about forcing term on scalar i- scalars advected with particles method
!    integer, dimension(:), allocatable, private     :: sc_partVector_forc
! XXX End Todo

    ! MHD variables
    !> Magnetic Prandl number, Default to 0E+00
    REAL(WP), PRIVATE :: Prm
    !> Storage for the Spectral values of B at step i-1 for RK2
    TYPE(COMPLEX_DATA_LAYOUT), ALLOCATABLE, DIMENSION(:), PRIVATE :: BkOld
    ! External imposed magnetic field
    !> External imposed Bx
    REAL(WP), PRIVATE :: Bxext
    !> External imposed By
    REAL(WP), PRIVATE :: Byext
    !> External imposed Bz
    REAL(WP), PRIVATE :: Bzext

    !> G (propagator) and sigmak complex layout for the Fourier transform
    !TYPE(COMPLEX_DATA_LAYOUT)     :: sigmak,G!,risu

    !> Velocity used to advect scalar field solved with particle/spectral method
    !! (as velocity is computed separatly from the velocity - a "trapeze" method is used)
    TYPE(REAL_DATA_LAYOUT)     , PRIVATE, SAVE    :: U_mean
    TYPE(REAL_DATA_LAYOUT)     , PRIVATE, SAVE    :: V_mean
    TYPE(REAL_DATA_LAYOUT)     , PRIVATE, SAVE    :: W_mean
    TYPE(COMPLEX_DATA_LAYOUT), PRIVATE, SAVE    :: Uk_mean
    TYPE(COMPLEX_DATA_LAYOUT), PRIVATE, SAVE    :: Vk_mean
    TYPE(COMPLEX_DATA_LAYOUT), PRIVATE, SAVE    :: Wk_mean

    !> To choose the right time integration
    !LUCA: I made public the solve_Euler_luca 
    private ::solve_spec_Euler, solve_spec_RK2, solve_spec_RK4, solve_spec_AB2, solve_spec_AB3
    !procedure(solve_spec_Euler), pointer, private :: solve_spec_basic => solve_spec_RK2
    procedure(solve_spec_Euler), pointer, private :: solve_spec_basic => NULL()

    !> Storage for intermadiate point, save and reconstruction of non-linear
    !term in RK and Adam-Smith time schemes
    TYPE(COMPLEX_DATA_LAYOUT), ALLOCATABLE, DIMENSION(:),PRIVATE, SAVE :: Uk_sav,Vk_sav,Wk_sav
    TYPE(COMPLEX_DATA_LAYOUT), ALLOCATABLE, DIMENSION(:,:),PRIVATE, SAVE :: Bk_sav, Scalk_sav
    TYPE(COMPLEX_DATA_LAYOUT), DIMENSION(3),PRIVATE, SAVE :: Fk_sav
    !> Number of field to save
    INTEGER :: save_size, dt_save
    LOGICAL :: forcing_save
    !> Save of previous time step for AB schemes
    real(WP), allocatable, dimension(:), private, save :: dt_prev

    !> total simulation time
    REAL(WP), PROTECTED :: sim_time
    !>timestep
    REAL(WP), PROTECTED ::step_time
    !>timestep
    REAL(WP), PROTECTED ::dtn0,dtn22

    !> number of scalars fields requested
    INTEGER,    PRIVATE :: nbscal=0
    !> number of scalars fields requested
    INTEGER,    PRIVATE :: nbscal_part=0
    !> forcing requested
    INTEGER,    PRIVATE :: iforce=0
    !> forcing amplitude for scalar forcing
    REAL(WP), DIMENSION(:), POINTER, PRIVATE :: amplitude_sf
    !> lower bound of wavenumber shell of scalar forcing
    REAL(WP), DIMENSION(:), POINTER, PRIVATE :: kmin_sf
    !> upper bound of wavenumber shell of scalar forcing
    REAL(WP), DIMENSION(:), POINTER, PRIVATE :: kmax_sf
    !> forcing amplitude for scalar forcing
    REAL(WP), DIMENSION(:), POINTER, PRIVATE :: amplitude_part_sf
    !> lower bound of wavenumber shell of scalar forcing
    REAL(WP), DIMENSION(:), POINTER, PRIVATE :: kmin_part_sf
    !> upper bound of wavenumber shell of scalar forcing
    REAL(WP), DIMENSION(:), POINTER, PRIVATE :: kmax_part_sf

    !> coefficient for stability condition in the particle solver
    REAL(WP), PRIVATE :: part_stab_coeff=0

    type(COMPLEX_DATA_LAYOUT), allocatable, dimension(:), public:: divTub

    PRIVATE wavesNumber_init
    PRIVATE deleteWavesNumber
    PRIVATE firstTerm
    PRIVATE secondTerm
    PRIVATE non_linear_vel
    PRIVATE non_linear_scal
    PRIVATE non_linear_mhd
    PRIVATE mhd_velocity_coupling
    PRIVATE get_timestep

!    PRIVATE part_spec_step
!    PRIVATE part_spec_step_fullPart
    !procedure(part_spec_step), pointer, private :: part_solver_step => part_spec_step
!    procedure(part_spec_step), pointer, private :: part_solver_step => NULL()

!    PRIVATE solve_velo_ScalSpec_basic
!    PRIVATE solve_velo_ScalSpec_subcycle

    LOGICAL,    PRIVATE :: initialized=.FALSE.

    COMPLEX(WP), PRIVATE, PARAMETER :: ii=(0.0_WP,1.0_WP)

    !> Solve Navier-Stokes equation for velocity and advection-diffusion equation
    !! for scalar (solved with full pseudo-spectral solver).

!  interface solve_velo_ScalSpec
!      module procedure solve_velo_ScalSpec_basic
!      module procedure solve_velo_ScalSpec_subcycle
!  end interface solve_velo_ScalSpec

  !> Compute first non-linear term for scalar or vector field
  interface firstTerm
      module procedure firstTerm_scal
      module procedure firstTerm_vect
  end interface firstTerm

  !> Update fields (velocity, B, scalarS) for one time step with the input "source term".
  interface spectral_update_fields
    module procedure spectral_update_fields_basic
    module procedure spectral_update_fields_two_dt
  end interface spectral_update_fields

CONTAINS

!------------------------------------------------------------------------------
!> @author
!> Luca Marradi
!
!> @details
!> This function initialize the solver module
!> @param [in,out]  G the propagator
!> @param [in,out]  sigma longitudinal velocity
!> @param [in]      nbcpus the total amount of processes (or cpus)
!> @param [in]      spec_rank my MPI processe rank in the range [0:nbcpus-1]
!> @param [out]     current_time = simulation time
!> @return .TRUE. if initialization is successfull.
!------------------------------------------------------------------------------
FUNCTION initSolver_luca (convol,n_activ,sigmak,G,sigma,sigma_old,Esigma,Lsigma,gammadot_min,current_time, nbcpus,spec_rank)

    USE cart_topology
    USE differential_tools

    IMPLICIT NONE

    !=== INPUT/OUTPUT ===
    INTEGER,                INTENT(IN)              :: nbcpus,spec_rank
    TYPE(REAL_DATA_LAYOUT), INTENT(INOUT)           :: sigma,sigma_old
    TYPE(REAL_DATA_LAYOUT), INTENT(INOUT)           :: Esigma,Lsigma,n_activ,convol
    TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT)        :: G,sigmak
    TYPE(WaveNumbers)                               :: kWN
    REAL(WP), INTENT(IN)                            :: gammadot_min
    REAL(WP), INTENT(OUT)                           :: current_time

    !=== LOCAL DATA ===
    LOGICAL initSolver_luca
    REAL(WP) :: delt
    INTEGER                 :: i,j,k,ierr,start
    CHARACTER(len=str_long) :: request
    CHARACTER(len=4)        :: advec_solver ! choose advection solver (p_ON = particle method at order N, N = 2 or 4)
    CHARACTER(len=4)        :: advec_interp! choose advection solver (p_ON = particle method at order N, N = 2 or 4)
    INTEGER                 :: ishowdiv
    CHARACTER(LEN=str_long) :: sim_name
    CHARACTER(LEN=str_long) :: les_name
    CHARACTER(LEN=str_long) :: TimeMethod_name
    CHARACTER(LEN=str_long) :: null_firstWN
    LOGICAL                 :: res          ! result return by the function
    integer                 :: disc_init    ! to init scalars to a disc shape.
    ! For particles method and assoicated scalar
    LOGICAL                 :: advanced     ! to  know if group size is set manually or not.
    INTEGER                 :: gp_size      ! particle solver gather line by group. It allow to choose group size
    INTEGER,DIMENSION(3)    :: gp_size_array ! particle solver gather line by group. It allow to choose group size

    initSolver_luca=.FALSE.

    ! Initial time
    if (parser_is_defined('Initial time')) then
      call parser_read('Initial time', sim_time)
    else
      sim_time = 0.0_WP
    end if
    current_time = sim_time

    !=== INITIALIZES THE WAVENUMBER FOR SIGMA ===
    IF(.NOT.wavesNumber_init(kWN,spec_rank,sigma%Lx,sigma%Ly,sigma%Lz,sigma%nx,sigma%ny,sigma%nz)) THEN
      WRITE(6,'(a,i0,a)')'[ERROR] init solver on process ',spec_rank,': wave-number for velocity field not initialize'
      RETURN
    ENDIF

    !=== GET INTEGRATION SCHEME ===
    CALL parser_read('Integration Method',TimeMethod_name)
    forcing_save = .false.
    time_integration = 1
    if (trim(TimeMethod_name)=='RK4') then
        solve_spec_basic => solve_spec_RK4
        save_size = 2
        dt_save = 0
    else if (trim(TimeMethod_name)=='RK2') then
        solve_spec_basic => solve_spec_RK2
        save_size = 1
        dt_save = 0
    else if (trim(TimeMethod_name)=='AB2') then
        solve_spec_basic => solve_spec_AB2_init
        save_size = 1
        dt_save = 1
    else if (trim(TimeMethod_name)=='AB3') then
        solve_spec_basic => solve_spec_AB3_init
        save_size = 2
        dt_save = 2
    else if (trim(TimeMethod_name)=='EULER') then
        solve_spec_basic => solve_spec_Euler
        save_size = 0
        dt_save = 0
    else
        solve_spec_basic => solve_spec_RK2
        time_integration = 1
        save_size = 1
        dt_save = 0
    end if

    ! Ask for the name of particular initial condition
    CALL parser_read('Simulation name', sim_name)

    !INITIAL CONDITION
    IF (trim(sim_name).EQ.'plasticflow') THEN
        IF(.NOT.initplasticflow(convol,n_activ,sigmak,G,sigma,sigma_old,Esigma,Lsigma &
                & ,gammadot_min,kWN,nbcpus,spec_rank)) THEN 
           RETURN
        ENDIF
    ELSEIF(trim(sim_name).EQ.'Zero') THEN
        sigma%values = 0.0_WP
    ELSE
        IF(spec_rank.EQ.0) WRITE(6,'(a)') '[ERROR] initSolver: unknown Simulation Name['// &
            & TRIM(sim_name)//']'
        RETURN
    ENDIF
    
   ! I should define delt just once in the code
   !delt=1./2.0_WP*acos(-1.0_WP)
 
 
   ! INITIALIZATION OF THE PROPAGATOR G
!   DO k=G%zmin,G%zmax
!    DO j=G%ymin,G%ymax
!     DO i=G%xmin,G%xmax
!        delt = (kWN%ky(j)**2_WP + kWN%kx(i)**2_WP)**2_WP
!        IF (delt.LT.1.E-20) THEN
!          delt = 10000_WP
!        ENDIF
!        delt=1_WP/delt
!        !TESTED
!        !!G%values(i,j,k) = spec_rank!delt * exp(-0.5*( (0.5)**2.) ) 
!        !TESTED FOR NX=2500, NY=10, NZ=10
!        !!G%values(i,j,k) = delt * exp(-0.5*( (kWN%kx(i) - kWN%kx(500))**2.)/100. )
!        !TESTED
!        !!G%values(i,j,k) =  kWN%kx(i)**2_WP !* delt
!        G%values(i,j,k) = -4_WP * (kWN%kx(i)**2_WP * kWN%ky(j)**2_WP) * delt
!        !IF(spec_rank.EQ.0) print*, "G(",i,",",j,",",k,"):",REALPART(G%values(i,j,k))
!     ENDDO
!    ENDDO
!   ENDDO

!TESTED
!! WAVENUMBER VALUES
!! IF(spec_rank.eq.0) THEN
!!  OPEN(400,file='G-propagator.dat',form='formatted')
!!   DO k=1,sigma%nz
!!    DO j=1,sigma%ny
!!     DO i=1,((sigma%nx)/2)+1
!!       write(400,'(g15.8,1x)') kWN%kx(i), kWN%ky(j), kWN%kz(k)
!!     ENDDO
!!    ENDDO
!!   ENDDO
!!  CLOSE(400)
!!ENDIF

!!    call btran(G,sigma,res)

!!   DO k=1,sigma%nz
!!    DO j=1,sigma%ny
!!     DO i=1,(sigma%nx)/2 + 1
        !!G_full(i,j,k) = delt * exp(-0.5*( (kWN%kx(i) - kWN%kx(500))**2.)/100.) 
        !G_full(i,j,k) =  REALPART(G%values(i,j,k))
!!     ENDDO
!!    ENDDO
!!   ENDDO

   ! Now check if we want to overide initialization with existing datafiles
   !! IF (.NOT. data_read(sigma, spec_rank)) then
   !!     WRITE(6,'(a,i0)')'[ERROR] init solver - failed to read data on process',spec_rank
   !!     RETURN
   !! end if

!   CALL ftran(sigma,sigmaK,res)
!    IF (.NOT. res) THEN
!        WRITE(6,'(a)')'[ERROR] initSolver_luca: cannot initialize sigma in fourier space'
!        RETURN
!    ENDIF

    !!if(.not. init_save(sigmak,save_size)) then
    !!  write(6,'(a,i0)') '[ERROR] iinit_solver: init_save failed'
    !!  return
    !!end if


    initialized=.TRUE.
    initSolver_luca=.TRUE.
    RETURN

END FUNCTION initSolver_luca

!------------------------------------------------------------------------------
!> @author
!> Patrick BEGOU, LEGI
!
!>
!> @details
!> This function free the memory used by solver module
!------------------------------------------------------------------------------
!SUBROUTINE deleteSolver (imhd, me)
!
!    use differential_tools
!
!    logical, INTENT(IN)             :: imhd
!    INTEGER, INTENT(IN), OPTIONAL   ::me
!    INTEGER :: i
!
!    IF (.NOT. initialized) RETURN
!
!    IF (PRESENT(me)) THEN
!         IF(me .EQ. 0) WRITE(6,'(a,g15.8)') '[PROGRESSION] Simulation time reached: ',sim_time
!    ENDIF
!
!    if(.not. delete_save(save_size, imhd)) then
!      write(6,'(a,i0)') '[ERROR] delete_solver: delete_save failed'
!      return
!    end if
!
!    ! free velocities wave numbers
!    CALL deleteWavesNumber(VelWN)
!
!    ! free scalars wave numbers
!    IF (ASSOCIATED(ScalWN)) THEN
!        DO I=1, SIZE(ScalWN)
!             CALL deleteWavesNumber(ScalWN(i))
!        ENDDO
!        DEALLOCATE(ScalWN)
!    ENDIF
!
!    IF (ASSOCIATED(Scal_partWN)) THEN
!        DO I=1, SIZE(Scal_partWN)
!             CALL deleteWavesNumber(Scal_partWN(i))
!        ENDDO
!        DEALLOCATE(Scal_partWN)
!    ENDIF
!
!    IF (ALLOCATED(ScalArrayK)) THEN
!        DO I=1, SIZE(ScalArrayK)
!             CALL deleteDataLayout(ScalArrayK(i))
!        ENDDO
!        DEALLOCATE(ScalArrayK)
!    ENDIF
!
!    IF (ALLOCATED(Scal_partArrayK)) THEN
!        DO I=1, SIZE(Scal_partArrayK)
!             CALL deleteDataLayout(Scal_partArrayK(i))
!        ENDDO
!        DEALLOCATE(Scal_partArrayK)
!    ENDIF
!
!!!!    CALL delete_LES(SIZE(ScalArrayK))
!    CALL delete_LES()
!
!    IF (ASSOCIATED(schmidt)) DEALLOCATE(schmidt)
!    IF (ASSOCIATED(ilessc)) DEALLOCATE(ilessc)
!    IF (ASSOCIATED(schmidt_part)) DEALLOCATE(schmidt_part)
!    IF (ALLOCATED(res_t)) DEALLOCATE(res_t)
!    CALL deleteDataLayout(Uk)
!    CALL deleteDataLayout(Vk)
!    CALL deleteDataLayout(Wk)
!
!    if(imhd) then
!        CALL deleteDataLayout(Bk(1))
!        CALL deleteDataLayout(Bk(2))
!        CALL deleteDataLayout(Bk(3))
!        deallocate(Bk)
!    end if
!
!
!    call diff_tool_deleteBuff(imhd)
!    !CALL deleteDataLayout(phybuff)
!
!    CALL deleteDataLayout(U_mean)
!    CALL deleteDataLayout(V_mean)
!    CALL deleteDataLayout(W_mean)
!    CALL deleteDataLayout(Uk_mean)
!    CALL deleteDataLayout(Vk_mean)
!    CALL deleteDataLayout(Wk_mean)
!
!END SUBROUTINE deleteSolver

!------------------------------------------------------------------------------
!> @author
!> Luca MARRADI, LIPhy 2014
!
!> This function execute one time step
!> @param [in,out]  sigma stress (physical space)
!> @param [in,out]  L_sigma probability
!> @param [in]      G the propagator 
!> @param [in]      nbcpus the total amount of processes (or cpus)
!> @param [in]      spec_rank my MPI processe rank in the range [0:nbcpus-1]
!> @param [out]     time    = reached time
!> @return .TRUE. if the time-step is successfull and U,V,W,ScalArray after this
!new timestep.
!> @details
!! This function solve equation for one time step. It updates the field in both
!! real and spectral space. Be aware that, for optimisation reason, scalar
!! advected with particle method are up to date in Fourrier space but not in real
!! space (to avoid a inverse FFT transform, as computations start and end in
!! Fourrier space). This does not affect solver precision, but only the
!output.
!------------------------------------------------------------------------------
 FUNCTION solver_step_luca(convol,n_activ,sigma,sigma_old,Esigma,Lsigma,G,sigmak,nbcpus,spec_rank,t_step,T) result(res)
   
   USE avs
   USE mpi
 
   !=== INPUT/OUTPUT DATA
   INTEGER,                INTENT(IN)              :: nbcpus,spec_rank
   TYPE(REAL_DATA_LAYOUT), INTENT(INOUT)           :: sigma,sigma_old,n_activ,convol
   TYPE(REAL_DATA_LAYOUT), INTENT(INOUT)           :: Esigma,Lsigma
   TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT)        :: G,sigmak
   REAL(WP), INTENT(IN)                            :: t_step,T
   LOGICAL res
   
  res=.FALSE.
  IF (.NOT.initialized) THEN
      IF(spec_rank.EQ.0) WRITE(6,'(a,f14.6,a)') '[ERROR] solver_step called before initialization!'
      RETURN
  ENDIF

 END FUNCTION solver_step_luca

!------------------------------------------------------------------------------
!> @author
!> Luca MARRADI, LIPhy 2014
!
!> Evolve the equations by EULER methods
!> @param [in,out]  sigma (physical space)
!> @param [in,out]  sigma_old
!> @param [in,out]  n_activ
!> @param [in,out]  mult
!> @param [in]      G the propagator 
!> @param [in]      nbcpus the total amount of processes (or cpus)
!> @param [in]      spec_rank my MPI process rank in the range [0:nbcpus-1]
!> @param [out]     time    = reached time
!> @param [in]      T temperature
!> @return .TRUE. if the time-step is successfull 
!> @details
!! This function solve equation for one time step. It update the field in 
!! real space by using an Euler scheme. 
!------------------------------------------------------------------------------
 FUNCTION solve_Euler_luca(convol,n_activ,mult,sigma,sigma_old,Esigma,Lsigma,gammadot_min,sigmak,G,t_step,T,nbcpus,spec_rank,ite,Nw,freq_sigma_corr) RESULT(res)

 
   !USE avs
   USE mpi
   USE plasticflow
   USE transforms_tools
   USE random_tools
   
   IMPLICIT NONE

   !== INPUT/OUTPUT DATA ==
   INTEGER,                  INTENT(IN)       :: nbcpus,spec_rank,ite,Nw,freq_sigma_corr
   TYPE(REAL_DATA_LAYOUT),   INTENT(INOUT)    :: sigma,sigma_old,n_activ
   TYPE(REAL_DATA_LAYOUT),   INTENT(INOUT)    :: mult,convol
   TYPE(REAL_DATA_LAYOUT),   INTENT(INOUT)    :: Esigma,Lsigma
   TYPE(COMPLEX_DATA_LAYOUT),INTENT(INOUT)    :: sigmak,G
   REAL(WP), INTENT(IN)                       :: t_step,T,gammadot_min
   
   !== LOCAL DATA ==
   INTEGER                                    :: i,j,k,ierr
   INTEGER                                    :: seed,freq,l
!!!   REAL(WP), DIMENSION(:,:,:), ALLOCATABLE    :: L_sigma
   REAL(WP)                                   :: psr
   LOGICAL                                    :: success
   LOGICAL                                    :: res

   res=.FALSE.
   freq = 10
   seed = spec_rank
   !== ALLOCATE MEMORY ==
!!   ALLOCATE(L_sigma(sigma%xmin:sigma%xmax,sigma%ymin:sigma%ymax,sigma%zmin:sigma%zmax))

   !== CHECK FOR ELASTICITY PLASTICITY EVENTS ==
   IF(.NOT.probabilityL(sigma,Lsigma,T,t_step,spec_rank)) THEN
      CALL MPI_FINALIZE(ierr)
      STOP 'FAILED to compute probability from elastic to plastic'
   ENDIF   
   IF(.NOT.probabilityE(sigma,Esigma,T,t_step,spec_rank)) THEN
      CALL MPI_FINALIZE(ierr)
      STOP 'FAILED to compute probability from plastic to elastic'
   ENDIF

 DO k=sigma%zmin,sigma%zmax
  DO j=sigma%ymin,sigma%ymax
   DO i=sigma%xmin,sigma%xmax
     !CALL random_number(psr)
     psr = ran3(seed)
     IF ( (n_activ%values(i,j,k).EQ.0_WP) .AND. (psr.LT.(Lsigma%values(i,j,k))) &
           .AND. (sigma%values(i,j,k) .GT. 1_WP ))  THEN
       n_activ%values(i,j,k) = 1_WP
     ELSE
     !CALL random_number(psr)
     psr = ran3(seed)
      if ( (n_activ%values(i,j,k).EQ.1_WP) .AND. (psr.LT.Esigma%values(i,j,k)) )  then
           n_activ%values(i,j,k) = 0_WP
      endif
     ENDIF
     convol%values(i,j,k) = n_activ%values(i,j,k)*sigma%values(i,j,k)
    ENDDO
  ENDDO
 ENDDO

   !== FROM REAL TO FOURIER SPACE ==
   CALL ftran(convol,sigmak,success)
  
   !== CONVOLUTION ==
   CALL convolution(sigmak,sigmak,G,success)

   !== FROM FOURIER TO REAL SPACE ==
   CALL btran(sigmak,convol,res)

   !== UPDATE SIGMA ==
   DO k=sigma%zmin,sigma%zmax
     DO j=sigma%ymin,sigma%ymax
      DO i=sigma%xmin,sigma%xmax
       sigma%values(i,j,k) = sigma%values(i,j,k) + t_step*(convol%values(i,j,k) + gammadot_min)
      ENDDO
     ENDDO
    ENDDO

  !=== COMPUTE SIGMA AT WAITING TIME === 
  IF(Nw.EQ.ite) THEN
     do k=sigma%zmin,sigma%zmax
      do j=sigma%ymin,sigma%ymax
       do i=sigma%xmin,sigma%xmax
        sigma_old%values(i,j,k) = sigma%values(i,j,k)
       enddo
      enddo
     enddo
  ENDIF

  IF (freq_sigma_corr>0 .AND. (mod((ite-Nw), freq_sigma_corr)== 0).AND.(ite>Nw)) THEN
    do k=sigma%zmin,sigma%zmax
     do j=sigma%ymin,sigma%ymax
      do i=sigma%xmin,sigma%xmax
        mult%values(i,j,k)= sigma_old%values(i,j,k) * sigma%values(i,j,k)
      enddo
     enddo
    enddo
  ENDIF



   !DEALLOCATE MEMORY
!!   DEALLOCATE(L_sigma)

    !=== TEST PRINT OUT PROPAGATOR ===
    !TESTED: 
!   IF(.NOT.real_WriteToFile(spec_rank,trim(convol%name),convol))  THEN
!     WRITE(*,*) "does not work!!!"
!   ENDIF

   res=.TRUE.

END FUNCTION solve_Euler_luca 

!> Explicit Euler Solve Navier-Stokes equation for velocity and advection-diffusion equation
!! for scalar (solved with full pseudo-spectral solver).
!! @author
!! Mouloud Kessar and Jean-Baptiste Lagaert, LEGI
!!     @param[in,out]  U            = longitudinal velocity (physical space)
!!     @param[in,out]  V            = vertical velocity  (physical space)
!!     @param[in,out]  W            = spanwise velocity  (physical space)
!!     @param[in,out]  Uk           = longitudinal velocity (spectral space)
!!     @param[in,out]  Vk           = vertical velocity  (spectral space)
!!     @param[in,out]  Wk           = spanwise velocity  (spectral space)
!!     @param[in,out]  nlx          = nonlinear terms for U (fourier space), initialy set by firstTerm
!!     @param[in,out]  nly          = nonlinear terms for V (fourier space), initialy set by firstTerm
!!     @param[in,out]  nlz          = nonlinear terms for W (fourier space), initialy set by firstTerm
!!     @param[in]      spec_rank    = my MPI processe rank in the range [0:nbcpus-1]
!!     @param[in]      imhd         = 1 if MHD is active, 0 else
!!     @param[out]     time         = reached time
!!     @param[in, out] B            = magnetic field
!!     @param[out]     res          = true if no error occurs
!! @details
!!        Velocity follows and Navier-Stokes equation and the scalar fields follow a
!!    diffusion-advection equation. This subroutine solve it for a time intervall of
!!    [T ; T+dt]. This function update these fields (by solving the equations)
!!    both spectral and real space.
!!        This function provide a first order intergrator: explicit Euler Scheme.
!!    If the velocity field in real space is already on the good mesh
!!    (considering the scalar resolution), then solver used it. Else, it compute
!!    this vector on the right resolution directly from the Fourrier space
!!    (and thus perform Fourrier interpolation).
subroutine solve_spec_Euler(U, V, W, Uk, Vk, Wk, nlx, nly, nlz, &
                      & B, nlB, ScalS, ScalP, ScalNL,         &
                      & dt, imhd, Srank, res)

  use differential_tools

  ! Input/Output
  INTEGER,                INTENT(IN)                      :: Srank
  TYPE(REAL_DATA_LAYOUT), INTENT(INOUT)                   :: U,V,W
  TYPE(REAL_DATA_LAYOUT)                                  :: sigma
  type(COMPLEX_DATA_LAYOUT)                               :: sigmak
  type(COMPLEX_DATA_LAYOUT), intent(inout)                :: Uk,Vk,Wk
  TYPE(COMPLEX_DATA_LAYOUT),INTENT(INOUT)                 :: nlx, nly, nlz
  TYPE(REAL_DATA_LAYOUT), POINTER, INTENT(INOUT)          :: B(:)
  TYPE(COMPLEX_DATA_LAYOUT), DIMENSION(:), INTENT(INOUT)  :: nlB
  TYPE(REAL_DATA_LAYOUT), POINTER, INTENT(INOUT)          :: ScalS(:)
  TYPE(REAL_DATA_LAYOUT), POINTER, INTENT(IN)             :: ScalP(:)
  TYPE(COMPLEX_DATA_LAYOUT), DIMENSION(:), INTENT(INOUT)  :: ScalNL
  LOGICAL, INTENT(IN)                                     :: imhd
  real(WP)                                                :: dt
  LOGICAL                                                 :: res

  ! Initialisation
  res = .false.

  ! == Compute non-linear term, forcing and projection (div(U,V,W)=0) ==
  call spectral_source_term(U,V,W,Uk,Vk,Wk, &
                          & nlx, nly, nlz,  &
                          & imhd, B, nlB,   &
                          & ScalS,ScalNL,   &
                          & ScalP, dt,      &
                          & Srank, res)
  if(.not.res) then
    write(6,'(a)') '[ERROR] solve_spec_Euler : fail to compute source_term'
    return
  end if

  ! == compute diffusion and update complex fields ==
  call spectral_update_fields(Uk,Vk,Wk, Uk,     &
                     & Vk, Wk, nlx, nly, nlz,   &
                     & imhd,Bk,Bk,nlB,          &
                     & ScalArrayk, ScalArrayk, ScalNL,    &
                     & dt, res)
  if(.not.res) then
    write(6,'(a)') '[ERROR] solve_spec_Euler : fail to update complex fields'
    return
  end if

  ! == Update real fields ==
  call spectral_update_real_fields(U,V,W,Uk,Vk,Wk, &
                            & imhd, B, Bk,    &
                            & ScalS,ScalArrayk,   &
                            & res)
  if(.not.res) then
    write(6,'(a)') '[ERROR] solve_spec_Euler : fail to update real fields'
    return
  end if

  ! == Update time ==
  sim_time = sim_time + dt

end subroutine solve_spec_Euler

!-------------------------------------------------------------------------------
!> Explicit Euler solves plasticity equations for sigma stress tensor
!! for scalar (solved with full pseudo-spectral solver).
!! @author 
!! Luca MARRADI LIPhy, 2014
!!  @param[in,out] sigma       = stress tensor
!!  @param[in]     G           = propagator
!!  @param[in]     mu          = elasticity tensor
!!  @param[in]     sigmak      = stress tensor (spectral space)
!!  @param[in]      spec_rank  = my MPI processe rank in the range [0:nbcpus-1]
!!  @param[out]     time         = reached time
!!   @param[out]     res          = true if no error occurs
!! @details
!! This subroutine solve it for a time intervall of [T ; T+dt]. This function 
!! update these fields (by solving the equations) both spectral and real space.
!! This function provide a first order intergrator: explicit Euler Scheme.

!!SUBROUTINE solve_spec_Euler_luca(sigma,sigma_old,sigmak,G_prop,step,T,Srank,res)

!!  USE transforms_tools
!!  USE plasticflow

  ! Input/Output
!!  INTEGER                                         :: i,j,k
!!  INTEGER,                   INTENT(IN)           :: Srank
!!  TYPE(REAL_DATA_LAYOUT),    INTENT(INOUT)        :: sigma,sigma_old
!!  TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT)           :: G_prop,sigmak
!!  REAL(WP),DIMENSION(:,:,:), ALLOCATABLE          :: Lsigma
!!  REAL(WP)                                        :: step,T,psr
!!  LOGICAL                                         :: res

!! ALLOCATE(Lsigma(sigma%xmin:sigma%xmax,sigma%ymin:sigma%ymax,sigma%zmin:sigma%zmax)) 

 ! INITIALIZATION
!!  res = .FALSE.
 
 !== STOCHASTIC PART ==
!! CALL probability(sigma,Lsigma,T,step)  

 !== TEST OF PLASTICITY EVENT ==
!!  psr=REAL(Srank)
!!  DO k=sigma%zmin,sigma%zmax
!!   DO j=sigma%ymin,sigma%ymax
!!    DO i=sigma%xmin,sigma%xmax
!!     CALL random_number(psr)
!!    IF(Lsigma(i,j,k).LT.psr) sigma%values(i,j,k) = 0_WP
!!    ENDDO
!!   ENDDO
!!  ENDDO

 !== FROM REAL TO FOURIER SPACE ==
!!  CALL ftran(sigma,sigmak,res)

 !== CONVOLUTION ==
!!  CALL convolution(sigmak,sigmak,G_prop,res)

 !== FROM FOURIER TO REAL SPACE ==
!!  CALL btran(sigmak,sigma,res)

 !== UPDATE FIELD SIGMA ==
!! DO k = sigma%zmin,sigma%zmax
!!   DO j = sigma%ymin,sigma%ymax
!!     DO i = sigma%xmin,sigma%xmax
!!      sigma%values(i,j,k) = sigma_old%values(i,j,k) + step*sigma%values(i,j,k) !+ mu%values(i,j,k) )
!!     ENDDO
!!   ENDDO
!! ENDDO

 !== STORE SIGMA IN SIGMA_OLD FOR NEXT STEP ==
!! DO k = sigma%zmin,sigma%zmax
!!   DO j = sigma%ymin,sigma%ymax
!!     DO i = sigma%xmin,sigma%xmax
!!      sigma_old%values(i,j,k) = sigma%values(i,j,k) 
!!     ENDDO
!!   ENDDO
!! ENDDO

  !DEALLOCATE MEMORY
!! DEALLOCATE(Lsigma)

!! res=.TRUE.
!!END SUBROUTINE solve_spec_Euler_luca
!-------------------------------------------------------------------------------


!> Runge-Kutta 2 for spectral sovler (Navier-Stokes equation for velocity and advection-diffusion equation
!! for scalar).
!! @author
!! Patrick Begou, Jean-Baptiste Lagaert, Guillaume Balarac, LEGI
!! Mouloud Kessar and Jean-Baptiste Lagaert, LEGI
!!     @param[in,out]  U            = longitudinal velocity (physical space)
!!     @param[in,out]  V            = vertical velocity  (physical space)
!!     @param[in,out]  W            = spanwise velocity  (physical space)
!!     @param[in,out]  Uk           = longitudinal velocity (spectral space)
!!     @param[in,out]  Vk           = vertical velocity  (spectral space)
!!     @param[in,out]  Wk           = spanwise velocity  (spectral space)
!!     @param[in,out]  nlx          = nonlinear terms for U (fourier space), initialy set by firstTerm
!!     @param[in,out]  nly          = nonlinear terms for V (fourier space), initialy set by firstTerm
!!     @param[in,out]  nlz          = nonlinear terms for W (fourier space), initialy set by firstTerm
!!     @param[in]      spec_rank    = my MPI processe rank in the range [0:nbcpus-1]
!!     @param[in]      imhd         = 1 if MHD is active, 0 else
!!     @param[out]     time         = reached time
!!     @param[in, out] B            = magnetic field
!!     @param[out]     res          = true if no error occurs
!! @details
!!        Velocity follows and Navier-Stokes equation and the scalar fields follow a
!!    diffusion-advection equation. This subroutine solve it for a time intervall of
!!    [T ; T+dt]. This function update these fields (by solving the equations)
!!    both spectral and real space.
!!        This is a second order in timer integrator.
!!    If the velocity field in real space is already on the good mesh
!!    (considering the scalar resolution), then solver used it. Else, it compute
!!    this vector on the right resolution directly from the Fourrier space
!!    (and thus perform Fourrier interpolation).
subroutine solve_spec_RK2(U, V, W, Uk, Vk, Wk, nlx, nly, nlz, &
                      & B, nlB, ScalS, ScalP, ScalNL,         &
                      & dt, imhd, Srank, res)

  use differential_tools

  ! Input/Output
  INTEGER,                INTENT(IN)                      :: Srank
  TYPE(REAL_DATA_LAYOUT), INTENT(INOUT)                   :: U,V,W
  type(COMPLEX_DATA_LAYOUT), intent(inout)                :: Uk,Vk,Wk
  TYPE(REAL_DATA_LAYOUT)                                  :: sigma
  type(COMPLEX_DATA_LAYOUT)                               :: sigmak
  TYPE(COMPLEX_DATA_LAYOUT),INTENT(INOUT)                 :: nlx, nly, nlz
  TYPE(REAL_DATA_LAYOUT), POINTER, INTENT(INOUT)          :: B(:)
  TYPE(COMPLEX_DATA_LAYOUT), DIMENSION(:), INTENT(INOUT)  :: nlB
  TYPE(REAL_DATA_LAYOUT), POINTER, INTENT(INOUT)          :: ScalS(:)
  TYPE(REAL_DATA_LAYOUT), POINTER, INTENT(IN)             :: ScalP(:)
  TYPE(COMPLEX_DATA_LAYOUT), DIMENSION(:), INTENT(INOUT)  :: ScalNL
  LOGICAL, INTENT(IN)                                     :: imhd
  real(WP)                                                :: dt
  LOGICAL                                                 :: res

  ! Local
  !> Loop indices
  INTEGER  :: i
  REAL(WP) :: step

  ! == Initialisation ==
  res = .false.
  step = dt/2._WP

  ! ##### First Step : predict solution at time t+(dt/2) #####
  ! == Compute non-linear term, forcing and projection (div(U,V,W)=0) ==
  call spectral_source_term(U,V,W,Uk,Vk,Wk, &
                          & nlx, nly, nlz,  &
                          & imhd, B, nlB,   &
                          & ScalS,ScalNL,   &
                          & ScalP, step,    &
                          & Srank, res)
  if(.not.res) then
    write(6,'(a)') '[ERROR] solve_RK2 -step 1: fail to compute source_term'
    return
  end if

  ! == Save fields
  if(.not. (permute(Uk,Uk_sav(1)) .and. &
          & permute(Vk,Vk_sav(1)) .and. &
          & permute(Wk,Wk_sav(1)))) then
    write(6,'(a)') '[ERROR] solve_RK2 -step 1: fail to save velocity'
    return
  end if
  do i=1, nbscal
    if(.not. permute(ScalArrayk(i), Scalk_sav(i,1))) then
        write(6,'(a,i0)') '[ERROR] solve_RK2 -step 1: fail to save scalarS ', i
        return
    end if
  end do
  if(imhd) then
    do i=1, 3
      if(.not. permute(Bk(i),Bk_sav(i,1))) then
        write(6,'(a,i0)') '[ERROR] solve_RK2 -step 1: fail to save Bk ', i
        return
    end if
    end do
  end if

  ! == compute diffusion and update complex fields ==
  call spectral_update_fields(Uk,Vk,Wk, Uk_sav(1), &
                     & Vk_sav(1),Wk_sav(1),        &
                     & nlx, nly, nlz,   &
                     & imhd,Bk,Bk_sav(:,1),nlB,    &
                     & ScalArrayk,Scalk_sav(:,1),ScalNL,&
                     & step, res)
  if(.not.res) then
    write(6,'(a)') '[ERROR] solve_RK2 -step 1: fail to update complex fields'
    return
  end if

  ! == Update real fields ==
  call spectral_update_real_fields(U,V,W,Uk,Vk,Wk, &
                            & imhd, B, Bk,    &
                            & ScalS,ScalArrayk,   &
                            & res)
  if(.not.res) then
    write(6,'(a)') '[ERROR] solve_RK2 -step 1: fail to update real fields'
    return
  end if

  ! ##### Second Step : compute solution at time t+dt with source term from time t+dt/2 #####
  ! == Compute non-linear term, forcing and projection (div(U,V,W)=0) ==
  call spectral_source_term(U,V,W,Uk,Vk,Wk, &
                          & nlx, nly, nlz,  &
                          & imhd, B, nlB,   &
                          & ScalS,ScalNL,   &
                          & ScalP, dt,      &
                          & Srank, res)
  if(.not.res) then
    write(6,'(a)') '[ERROR] solve_RK2 -step 2: fail to compute source_term'
    return
  end if

  ! == compute diffusion and update complex fields ==
  call spectral_update_fields(Uk,Vk,Wk, Uk_sav(1), &
                     & Vk_sav(1),Wk_sav(1),        &
                     & nlx, nly, nlz,   &
                     & imhd,Bk,Bk_sav(:,1),nlB,    &
                     & ScalArrayk,Scalk_sav(:,1),ScalNL,&
                     & dt, step, res)
  if(.not.res) then
    write(6,'(a)') '[ERROR] solve_RK2 -step 2: fail to update complex fields'
    return
  end if

  ! == Update real fields ==
  call spectral_update_real_fields(U,V,W,Uk,Vk,Wk, &
                            & imhd, B, Bk,    &
                            & ScalS,ScalArrayk,   &
                            & res)
  if(.not.res) then
    write(6,'(a)') '[ERROR] solve_RK2 -step 2: fail to update real fields'
    return
  end if

  ! == Update time ==
  sim_time = sim_time + dt

end subroutine solve_spec_RK2


!> Runge-Kutta 4 for spectral sovler (Navier-Stokes equation for velocity and advection-diffusion equation
!! for scalar).
!! @author
!! Mouloud Kessar and Jean-Baptiste Lagaert, LEGI
!!     @param[in,out]  U            = longitudinal velocity (physical space)
!!     @param[in,out]  V            = vertical velocity  (physical space)
!!     @param[in,out]  W            = spanwise velocity  (physical space)
!!     @param[in,out]  Uk           = longitudinal velocity (spectral space)
!!     @param[in,out]  Vk           = vertical velocity  (spectral space)
!!     @param[in,out]  Wk           = spanwise velocity  (spectral space)
!!     @param[in,out]  nlx          = nonlinear terms for U (fourier space), initialy set by firstTerm
!!     @param[in,out]  nly          = nonlinear terms for V (fourier space), initialy set by firstTerm
!!     @param[in,out]  nlz          = nonlinear terms for W (fourier space), initialy set by firstTerm
!!     @param[in]      spec_rank    = my MPI processe rank in the range [0:nbcpus-1]
!!     @param[in]      imhd         = 1 if MHD is active, 0 else
!!     @param[out]     time         = reached time
!!     @param[in, out] B            = magnetic field
!!     @param[out]     res          = true if no error occurs
!! @details
!!        Velocity follows and Navier-Stokes equation and the scalar fields follow a
!!    diffusion-advection equation. This subroutine solve it for a time intervall of
!!    [T ; T+dt]. This function update these fields (by solving the equations)
!!    both spectral and real space.
!!        This is a fourth order in time RK scheme. See Appendix E in
!!    "Shell Models of MHD turubulence", F. Plunian, R. Stepanov and P. Frick in
!!    Physics Reports.
!!    If the velocity field in real space is already on the good mesh
!!    (considering the scalar resolution), then solver used it. Else, it compute
!!    this vector on the right resolution directly from the Fourrier space
!!    (and thus perform Fourrier interpolation).
subroutine solve_spec_RK4(U, V, W, Uk, Vk, Wk, nlx, nly, nlz, &
                      & B, nlB, ScalS, ScalP, ScalNL,         &
                      & dt, imhd, Srank, res)

  use differential_tools

  ! Input/Output
  INTEGER,                INTENT(IN)                      :: Srank
  TYPE(REAL_DATA_LAYOUT), INTENT(INOUT)                   :: U,V,W
  type(COMPLEX_DATA_LAYOUT), intent(inout)                :: Uk,Vk,Wk
  TYPE(REAL_DATA_LAYOUT)                                  :: sigma
  type(COMPLEX_DATA_LAYOUT)                               :: sigmak
  TYPE(COMPLEX_DATA_LAYOUT),INTENT(INOUT)                 :: nlx, nly, nlz
  TYPE(REAL_DATA_LAYOUT), POINTER, INTENT(INOUT)          :: B(:)
  TYPE(COMPLEX_DATA_LAYOUT), DIMENSION(:), INTENT(INOUT)  :: nlB
  TYPE(REAL_DATA_LAYOUT), POINTER, INTENT(INOUT)          :: ScalS(:)
  TYPE(REAL_DATA_LAYOUT), POINTER, INTENT(IN)             :: ScalP(:)
  TYPE(COMPLEX_DATA_LAYOUT), DIMENSION(:), INTENT(INOUT)  :: ScalNL
  LOGICAL, INTENT(IN)                                     :: imhd
  real(WP)                                                :: dt
  LOGICAL                                                 :: res
  ! Local
  !> Loop indices
  INTEGER  :: i,j,k,sca, force_sav
  !> To compute exponential decay (diffusion term)
  REAL(WP) :: kk, kkz, kkyz, ekk, ekk_bis, step, coef
  complex(WP) :: Axk, Ayk, Azk

  ! == Initialisation ==
  res = .false.
  step = dt/2._WP
  coef = - step*mu

  ! ##### First Step : predict solution at time t+(dt/2) #####
  ! == Compute non-linear term, forcing and projection (div(U,V,W)=0) ==
  call spectral_source_term(U,V,W,Uk,Vk,Wk, &
                          & nlx, nly, nlz,  &
                          & imhd, B, nlB,   &
                          & ScalS,ScalNL,   &
                          & ScalP, step,    &
                          & Srank, res, forcing_save)
  if(.not.res) then
    write(6,'(a)') '[ERROR] solve_RK4 -step 1: fail to compute source_term'
    return
  end if

  ! == Save fields
  if(.not. (permute(Uk,Uk_sav(1)) .and. &
          & permute(Vk,Vk_sav(1)) .and. &
          & permute(Wk,Wk_sav(1)))) then
    write(6,'(a)') '[ERROR] solve_RK4 -step 1: fail to save velocity'
    return
  end if
  do i=1, nbscal
    if(.not. permute(ScalArrayk(i), Scalk_sav(1,i))) then
        write(6,'(a,i0)') '[ERROR] solve_RK4 -step 1: fail to save scalarS ', i
        return
    end if
  end do
  if(imhd) then
    do i=1, 3
      if(.not. permute(Bk(i),Bk_sav(i,1))) then
        write(6,'(a,i0)') '[ERROR] solve_RK4 -step 1: fail to save Bk ', i
        return
    end if
    end do
  end if

  ! == compute diffusion and update complex fields ==
  call spectral_update_fields(Uk,Vk,Wk, Uk_sav(1), &
                     & Vk_sav(1),Wk_sav(1),        &
                     & nlx, nly, nlz,   &
                     & imhd,Bk,Bk_sav(:,1),nlB,    &
                     & ScalArrayk,Scalk_sav(:,1),ScalNL,&
                     & step, res)
  if(.not.res) then
    write(6,'(a)') '[ERROR] solve_RK4 -step 1: fail to update complex fields'
    return
  end if

  ! == Update final non-linear term ==
  ! We actually compute here (final NL)/(dt*exp(-mu*kk*dt/2))
  ! For velocty
  do k = Uk%zmin,Uk%zmax
    kkz = VelWN%kz(k)*VelWN%kz(k)
    do j = Uk%ymin,Uk%ymax
      kkyz = kkz + VelWN%ky(j)*VelWN%ky(j)
      do i = Uk%xmin,Uk%xmax
        kk = kkyz + VelWN%kx(i)*VelWN%kx(i)
        ekk=exp(coef*kk)
        !Uk_sav(2)%values(i,j,k) = ((sqrt(2._WP))*Fk_sav(1)%values(i,j,k) + nlx%values(i,j,k))*ekk
        !Vk_sav(2)%values(i,j,k) = ((sqrt(2._WP))*Fk_sav(2)%values(i,j,k) + nly%values(i,j,k))*ekk
        !Wk_sav(2)%values(i,j,k) = ((sqrt(2._WP))*Fk_sav(3)%values(i,j,k) + nlz%values(i,j,k))*ekk
        Uk_sav(2)%values(i,j,k) = nlx%values(i,j,k)*ekk
        Vk_sav(2)%values(i,j,k) = nly%values(i,j,k)*ekk
        Wk_sav(2)%values(i,j,k) = nlz%values(i,j,k)*ekk
      end do
    end do
  end do
  ! For scalars
  do sca=1,nbscal
    do k = ScalArrayK(sca)%zmin,ScalArrayK(sca)%zmax
      kkz = ScalWN(sca)%kz(k)*ScalWN(sca)%kz(k)
      do j = ScalArrayK(sca)%ymin,ScalArrayK(sca)%ymax
        kkyz = kkz + ScalWN(sca)%ky(j)*ScalWN(sca)%ky(j)
        do i = ScalArrayK(sca)%xmin,ScalArrayK(sca)%xmax
          kk = ScalWN(sca)%kx(i)*ScalWN(sca)%kx(i) + kkyz
          ekk=exp(kk*coef/schmidt(sca))
          Scalk_sav(sca,2)%values(i,j,k) = ScalNL(sca)%values(i,j,k)*ekk
        end do
      end do
    end do
  end do
  ! For MHD if needed
  if (imhd) then
    do k = Bk(1)%zmin,Bk(1)%zmax
      kkz = BfieldWN%kz(k)*BfieldWN%kz(k)
      do j = Bk(1)%ymin,Bk(1)%ymax
        kkyz = kkz + BfieldWN%ky(j)*BfieldWN%ky(j)
        do i = Bk(1)%xmin,Bk(1)%xmax
          kk = BfieldWN%kx(i)*BfieldWN%kx(i) + kkyz
          ekk=exp(kk*coef/Prm)
          Bk_sav(1,2)%values(i,j,k) = nlB(1)%values(i,j,k)*ekk
          Bk_sav(2,2)%values(i,j,k) = nlB(2)%values(i,j,k)*ekk
          Bk_sav(3,2)%values(i,j,k) = nlB(3)%values(i,j,k)*ekk
        end do
      end do
    end do
  end if ! mhd

  ! == Update real fields ==
  call spectral_update_real_fields(U,V,W,Uk,Vk,Wk, &
                            & imhd, B, Bk,    &
                            & ScalS,ScalArrayk,   &
                            & res)
  if(.not.res) then
    write(6,'(a)') '[ERROR] solve_RK4 -step 1: fail to update real fields'
    return
  end if

  ! ##### Second Step : correction at time t+dt/2 #####
  ! == Compute non-linear term, forcing and projection (div(U,V,W)=0) ==
  call spectral_source_term(U,V,W,Uk,Vk,Wk, &
                          & nlx, nly, nlz,  &
                          & imhd, B, nlB,   &
                          & ScalS,ScalNL,   &
                          & ScalP, step,      &
                          & Srank, res)
  if(.not.res) then
    write(6,'(a)') '[ERROR] solve_RK4 -step 2: fail to compute source_term'
    return
  end if

  ! == compute diffusion and update complex fields ==
  call spectral_update_fields(Uk,Vk,Wk, Uk_sav(1), &
                     & Vk_sav(1),Wk_sav(1),        &
                     & nlx, nly, nlz,   &
                     & imhd,Bk,Bk_sav(:,1),nlB,    &
                     & ScalArrayk,Scalk_sav(:,1),ScalNL,&
                     & dt, step, res)
  if(.not.res) then
    write(6,'(a)') '[ERROR] solve_RK4 -step 2: fail to update complex fields'
    return
  end if

  ! == Update real fields ==
  call spectral_update_real_fields(U,V,W,Uk,Vk,Wk, &
                            & imhd, B, Bk,    &
                            & ScalS,ScalArrayk,   &
                            & res)
  if(.not.res) then
    write(6,'(a)') '[ERROR] solve_RK4 -step 2: fail to update real fields'
    return
  end if

  ! == Update final non-linear term ==
  ! We actually compute here (final NL)/(dt*exp(-mu*kk*dt/2))
  ! For velocity
  Uk_sav(2)%values = Uk_sav(2)%values + 2._WP*nlx%values
  Vk_sav(2)%values = Vk_sav(2)%values + 2._WP*nly%values
  Wk_sav(2)%values = Wk_sav(2)%values + 2._WP*nlz%values
  ! For mhd
  if(imhd) then
    Bk_sav(1,2)%values = Bk_sav(1,2)%values + 2._WP*nlB(1)%values
    Bk_sav(2,2)%values = Bk_sav(2,2)%values + 2._WP*nlB(2)%values
    Bk_sav(3,2)%values = Bk_sav(3,2)%values + 2._WP*nlB(3)%values
  end if
  ! For scalar
  do sca=1,nbscal
    Scalk_sav(sca,2)%values = Scalk_sav(sca,2)%values + 2._WP*ScalNL(sca)%values
  end do

  ! ##### Third Step : prediction at time t+dt starting from time t+dt/2 #####
  ! == Compute non-linear term, forcing and projection (div(U,V,W)=0) ==
  call spectral_source_term(U,V,W,Uk,Vk,Wk, &
                          & nlx, nly, nlz,  &
                          & imhd, B, nlB,   &
                          & ScalS,ScalNL,   &
                          & ScalP, step,    &
                          & Srank, res)
  if(.not.res) then
    write(6,'(a)') '[ERROR] solve_RK4 -step 3: fail to compute source_term'
    return
  end if

  ! == compute diffusion and update complex fields ==
  call spectral_update_fields(Uk,Vk,Wk, Uk_sav(1), &
                     & Vk_sav(1),Wk_sav(1),        &
                     & nlx, nly, nlz,   &
                     & imhd,Bk,Bk_sav(:,1),nlB,    &
                     & ScalArrayk,Scalk_sav(:,1),ScalNL,&
                     & dt, step, res)
  if(.not.res) then
    write(6,'(a)') '[ERROR] solve_RK4 -step 3: fail to update complex fields'
    return
  end if

  ! == Update real fields ==
  call spectral_update_real_fields(U,V,W,Uk,Vk,Wk, &
                            & imhd, B, Bk,    &
                            & ScalS,ScalArrayk,   &
                            & res)
  if(.not.res) then
    write(6,'(a)') '[ERROR] solve_RK4 -step 3: fail to update real fields'
    return
  end if

  ! == Update final non-linear term ==
  ! We actually compute here (final NL)/(dt*exp(-mu*kk*dt/2))
  ! For velocity
  Uk_sav(2)%values = Uk_sav(2)%values + 2._WP*nlx%values
  Vk_sav(2)%values = Vk_sav(2)%values + 2._WP*nly%values
  Wk_sav(2)%values = Wk_sav(2)%values + 2._WP*nlz%values
  ! For mhd
  if(imhd) then
    Bk_sav(1,2)%values = Bk_sav(1,2)%values + 2._WP*nlB(1)%values
    Bk_sav(2,2)%values = Bk_sav(2,2)%values + 2._WP*nlB(2)%values
    Bk_sav(3,2)%values = Bk_sav(3,2)%values + 2._WP*nlB(3)%values
  end if
  ! For scalar
  do sca=1,nbscal
    Scalk_sav(sca,2)%values = Scalk_sav(sca,2)%values + 2._WP*ScalNL(sca)%values
  end do

  ! ##### Fourth Step : prediction at time t+dt using ponderation of previous results #####
  ! == Compute non-linear term, forcing and projection (div(U,V,W)=0) ==
  force_sav = iforce
  iforce = 0 ! no more forcing for this step
  call spectral_source_term(U,V,W,Uk,Vk,Wk, &
                          & nlx, nly, nlz,  &
                          & imhd, B, nlB,   &
                          & ScalS,ScalNL,   &
                          & ScalP, step,    &
                          & Srank, res)
  if(.not.res) then
    write(6,'(a)') '[ERROR] solve_RK4 -step 4: fail to compute source_term'
    return
  end if
  iforce = force_sav

  ! == Update spectral fields ==
  coef = -mu*dt
  ! For velocty
  do k = Uk%zmin,Uk%zmax
    kkz = VelWN%kz(k)*VelWN%kz(k)
    do j = Uk%ymin,Uk%ymax
      kkyz = kkz + VelWN%ky(j)*VelWN%ky(j)
      do i = Uk%xmin,Uk%xmax
        kk = kkyz + VelWN%kx(i)*VelWN%kx(i)
        ekk=exp(coef*kk)
        ekk_bis = exp(coef*kk/2._WP)
        Uk%values(i,j,k) = (Uk_sav(1)%values(i,j,k)*ekk) + (dt/6._WP)*(nlx%values(i,j,k)+ekk_bis*Uk_sav(2)%values(i,j,k))
        Vk%values(i,j,k) = (Vk_sav(1)%values(i,j,k)*ekk) + (dt/6._WP)*(nly%values(i,j,k)+ekk_bis*Vk_sav(2)%values(i,j,k))
        Wk%values(i,j,k) = (Wk_sav(1)%values(i,j,k)*ekk) + (dt/6._WP)*(nlz%values(i,j,k)+ekk_bis*Wk_sav(2)%values(i,j,k))
      end do
    end do
  end do
  ! For scalars
  do sca=1,nbscal
    do k = ScalArrayK(sca)%zmin,ScalArrayK(sca)%zmax
      kkz = ScalWN(sca)%kz(k)*ScalWN(sca)%kz(k)
      do j = ScalArrayK(sca)%ymin,ScalArrayK(sca)%ymax
        kkyz = kkz + ScalWN(sca)%ky(j)*ScalWN(sca)%ky(j)
        do i = ScalArrayK(sca)%xmin,ScalArrayK(sca)%xmax
          kk = ScalWN(sca)%kx(i)*ScalWN(sca)%kx(i) + kkyz
          ekk=exp(kk*coef/schmidt(sca))
          ekk_bis = exp(coef*kk/(2._WP*schmidt(sca)))
          ScalArrayk(sca)%values(i,j,k) = (Scalk_sav(sca,1)%values(i,j,k)*ekk) + &
                & (dt/6._WP)*(nlx%values(i,j,k)+ekk_bis* Scalk_sav(sca,2)%values(i,j,k))
        end do
      end do
    end do
  end do
  ! For MHD if needed
  if (imhd) then
    do k = Bk(1)%zmin,Bk(1)%zmax
      kkz = BfieldWN%kz(k)*BfieldWN%kz(k)
      do j = Bk(1)%ymin,Bk(1)%ymax
        kkyz = kkz + BfieldWN%ky(j)*BfieldWN%ky(j)
        do i = Bk(1)%xmin,Bk(1)%xmax
          kk = BfieldWN%kx(i)*BfieldWN%kx(i) + kkyz
          ekk=exp(kk*coef/Prm)
          ekk_bis=exp(kk*coef/(2._WP*Prm))
          nlB(1)%values(i,j,k) = nlB(1)%values(i,j,k) + ekk_bis*Bk_sav(1,2)%values(i,j,k)
          nlB(2)%values(i,j,k) = nlB(2)%values(i,j,k) + ekk_bis*Bk_sav(2,2)%values(i,j,k)
          nlB(3)%values(i,j,k) = nlB(3)%values(i,j,k) + ekk_bis*Bk_sav(3,2)%values(i,j,k)
          Axk = ii*BfieldWN%ky(j)*nlB(3)%values(i,j,k) - ii*BfieldWN%kz(k)*nlB(2)%values(i,j,k)
          Ayk = ii*BfieldWN%kz(k)*nlB(1)%values(i,j,k) - ii*BfieldWN%kx(i)*nlB(3)%values(i,j,k)
          Azk = ii*BfieldWN%kx(i)*nlB(2)%values(i,j,k) - ii*BfieldWN%ky(j)*nlB(1)%values(i,j,k)
          Bk(1)%values(i,j,k) = (Bk_sav(1,1)%values(i,j,k)*ekk) + (dt/6._WP)*Axk
          Bk(2)%values(i,j,k) = (Bk_sav(2,1)%values(i,j,k)*ekk) + (dt/6._WP)*Ayk
          Bk(3)%values(i,j,k) = (Bk_sav(3,1)%values(i,j,k)*ekk) + (dt/6._WP)*Azk
        end do
      end do
    end do
  end if ! mhd


  ! == Update real fields ==
  call spectral_update_real_fields(U,V,W,Uk,Vk,Wk, &
                            & imhd, B, Bk,    &
                            & ScalS,ScalArrayk,   &
                            & res)
  if(.not.res) then
    write(6,'(a)') '[ERROR] solve_RK4 -step 4: fail to update real fields'
    return
  end if

  ! == Update time ==
  sim_time = sim_time + dt

end subroutine solve_spec_RK4


!> Adams-Barshforth for spectral sovler (Navier-Stokes equation for velocity and advection-diffusion equation
!! for scalar).
!! @author
!! Mouloud Kessar and Jean-Baptiste Lagaert, LEGI
!!     @param[in,out]  U            = longitudinal velocity (physical space)
!!     @param[in,out]  V            = vertical velocity  (physical space)
!!     @param[in,out]  W            = spanwise velocity  (physical space)
!!     @param[in,out]  Uk           = longitudinal velocity (spectral space)
!!     @param[in,out]  Vk           = vertical velocity  (spectral space)
!!     @param[in,out]  Wk           = spanwise velocity  (spectral space)
!!     @param[in,out]  nlx          = nonlinear terms for U (fourier space), initialy set by firstTerm
!!     @param[in,out]  nly          = nonlinear terms for V (fourier space), initialy set by firstTerm
!!     @param[in,out]  nlz          = nonlinear terms for W (fourier space), initialy set by firstTerm
!!     @param[in]      spec_rank    = my MPI processe rank in the range [0:nbcpus-1]
!!     @param[in]      imhd         = 1 if MHD is active, 0 else
!!     @param[out]     time         = reached time
!!     @param[in, out] B            = magnetic field
!!     @param[out]     res          = true if no error occurs
!! @details
!!        Velocity follows and Navier-Stokes equation and the scalar fields follow a
!!    diffusion-advection equation. This subroutine solve it for a time intervall of
!!    [T ; T+dt]. This function update these fields (by solving the equations)
!!    both spectral and real space.
!!        This is a second order multpile step Adams-Bashforth. See Section 3.3.1 in
!!    "Simulation of confined MHD flows usinf a pseudo-spectral method with
!!    volume penalization"', J.A. Morales, M. Leroy, W.J.T. Bos and K. Schneider
!!    to appear in JCP (available in HAL, ref hal-00719737)
!!    If the velocity field in real space is already on the good mesh
!!    (considering the scalar resolution), then solver used it. Else, it compute
!!    this vector on the right resolution directly from the Fourrier space
!!    (and thus perform Fourrier interpolation).
!!        For an alvelius forcing, at the statiscal steady state, the effective
!!    forcing power is, by assuming dt constant and using the notation from [AB]
!!    P_eff = P*(beta_10 + beta_11/2) = (3/2-1/2*1/2)P = 1.25*P
subroutine solve_spec_AB2(U,V,W,Uk,Vk,Wk,nlx,nly,nlz, &
                      & B,nlB,ScalS,ScalP,ScalNL,   &
                      & dt, imhd, Srank, res)

  use differential_tools

  ! Input/Output
  INTEGER,                INTENT(IN)                      :: Srank
  TYPE(REAL_DATA_LAYOUT), INTENT(INOUT)                   :: U,V,W
  type(COMPLEX_DATA_LAYOUT), intent(inout)                :: Uk,Vk,Wk
  TYPE(REAL_DATA_LAYOUT)                                  :: sigma
  type(COMPLEX_DATA_LAYOUT)                               :: sigmak
  TYPE(COMPLEX_DATA_LAYOUT),INTENT(INOUT)                 :: nlx, nly, nlz
  TYPE(REAL_DATA_LAYOUT), POINTER, INTENT(INOUT)          :: B(:)
  TYPE(COMPLEX_DATA_LAYOUT), DIMENSION(:), INTENT(INOUT)  :: nlB
  TYPE(REAL_DATA_LAYOUT), POINTER, INTENT(INOUT)          :: ScalS(:)
  TYPE(REAL_DATA_LAYOUT), POINTER, INTENT(IN)             :: ScalP(:)
  TYPE(COMPLEX_DATA_LAYOUT), DIMENSION(:), INTENT(INOUT)  :: ScalNL
  LOGICAL, INTENT(IN)                                     :: imhd
  real(WP)                                                :: dt
  LOGICAL                                                 :: res
  ! Local
  !> Loop indices
  INTEGER  :: i
  !> Ponderation of non-linear term computed at different time-steps
  real(WP)  :: coef1, coef2

  ! == Initialisation ==
  res = .false.

  ! == Compute non-linear term, forcing and projection (div(U,V,W)=0) ==
  call spectral_source_term(U,V,W,Uk,Vk,Wk, &
                          & nlx, nly, nlz,  &
                          & imhd, B, nlB,   &
                          & ScalS,ScalNL,   &
                          & ScalP, dt,      &
                          & Srank, res)
  if(.not.res) then
    write(6,'(a)') '[ERROR] solve_AB2: fail to compute source_term'
    return
  end if

  ! == Save Non Linear terms for next time-step ==
  if(.not. (permute(nlx,Uk_sav(1)) .and. &
          & permute(nly,Vk_sav(1)) .and. &
          & permute(nlz,Wk_sav(1)))) then
    write(6,'(a)') '[ERROR] solve_AB2: fail to save non-linear terms for velocity'
    return
  end if
  do i=1, nbscal
    if(.not. permute(ScalNL(i), Scalk_sav(i,1))) then
        write(6,'(a,i0)') '[ERROR] solve_AB2: fail to save scalarS_nl ', i
        return
    end if
  end do
  if(imhd) then
    do i=1, 3
      if(.not. permute(nlB(i),Bk_sav(i,1))) then
        write(6,'(a,i0)') '[ERROR] solve_AB2: fail to save nlB ', i
        return
    end if
    end do
  end if

  ! == Update NL terms from ones computed in previous step ==
  ! Be careful that due to previous permutation:
  !  1 = nl     = non-linear terms computed at previous time step
  !  2 = nl_sav = non-linear terms computed at current time step
  ! Compute coefficients (divided by dt)
  coef1= -dt/(2._WP*dt_prev(1))
  coef2 = (dt+2._WP*dt_prev(1))/(2._WP*dt_prev(1))
  call compute_AB_NL(nlx, nly, nlz, Uk_sav(1), Vk_sav(1), Wk_sav(1),&
              & imhd,nlB,Bk_sav(1,1),Bk_sav(2,1),Bk_sav(3,1),       &
              & ScalNL, Scalk_sav(:,1), coef1, coef2,               &
              & dt_prev(1)/2._WP, res)
  if(.not.res) then
    write(6,'(a)') '[ERROR] solve_AB2: fail to update non-linear term'
    return
  end if

  ! == compute diffusion and update complex fields ==
  call spectral_update_fields(Uk,Vk,Wk, Uk, &
                     & Vk,Wk,nlx, nly, nlz, &
                     & imhd,Bk,Bk,nlB,      &
                     & ScalArrayk,ScalArrayk,ScalNL,&
                     & dt, res)
  if(.not.res) then
    write(6,'(a)') '[ERROR] solve_AB2: fail to update complex fields'
    return
  end if

  ! == Update real fields ==
  call spectral_update_real_fields(U,V,W,Uk,Vk,Wk, &
                            & imhd, B, Bk,    &
                            & ScalS,ScalArrayk,   &
                            & res)
  if(.not.res) then
    write(6,'(a)') '[ERROR] solve_AB2: fail to update real fields'
    return
  end if

  ! == Update time ==
  sim_time = sim_time + dt
  dt_prev(1) = dt

end subroutine solve_spec_AB2


!> Adams-Barshforth initialisation for spectral sovler (Navier-Stokes equation for velocity and advection-diffusion equation
!! for scalar): this multi-step method cannot be used for the first step.
!! @author
!! Mouloud Kessar and Jean-Baptiste Lagaert, LEGI
!!     @param[in,out]  U            = longitudinal velocity (physical space)
!!     @param[in,out]  V            = vertical velocity  (physical space)
!!     @param[in,out]  W            = spanwise velocity  (physical space)
!!     @param[in,out]  Uk           = longitudinal velocity (spectral space)
!!     @param[in,out]  Vk           = vertical velocity  (spectral space)
!!     @param[in,out]  Wk           = spanwise velocity  (spectral space)
!!     @param[in,out]  nlx          = nonlinear terms for U (fourier space), initialy set by firstTerm
!!     @param[in,out]  nly          = nonlinear terms for V (fourier space), initialy set by firstTerm
!!     @param[in,out]  nlz          = nonlinear terms for W (fourier space), initialy set by firstTerm
!!     @param[in]      spec_rank    = my MPI processe rank in the range [0:nbcpus-1]
!!     @param[in]      imhd         = 1 if MHD is active, 0 else
!!     @param[out]     time         = reached time
!!     @param[in, out] B            = magnetic field
!!     @param[out]     res          = true if no error occurs
!! @details
!!        Velocity follows and Navier-Stokes equation and the scalar fields follow a
!!    diffusion-advection equation. This subroutine solve it for a time intervall of
!!    [T ; T+dt]. This function update these fields (by solving the equations)
!!    both spectral and real space.
!!        This called an Euler or RK method for the first time step.
!!    If the velocity field in real space is already on the good mesh
!!    (considering the scalar resolution), then solver used it. Else, it compute
!!    this vector on the right resolution directly from the Fourrier space
!!    (and thus perform Fourrier interpolation).
subroutine solve_spec_AB2_init(U,V,W,Uk,Vk,Wk,nlx,nly,nlz,  &
                      & B, nlB, ScalS, ScalP, ScalNL,     &
                      & dt, imhd, Srank, res)

  use differential_tools

  ! Input/Output
  INTEGER,                INTENT(IN)                      :: Srank
  TYPE(REAL_DATA_LAYOUT), INTENT(INOUT)                   :: U,V,W
  type(COMPLEX_DATA_LAYOUT), intent(inout)                :: Uk,Vk,Wk
  TYPE(REAL_DATA_LAYOUT)                                  :: sigma
  type(COMPLEX_DATA_LAYOUT)                               :: sigmak
  TYPE(COMPLEX_DATA_LAYOUT),INTENT(INOUT)                 :: nlx, nly, nlz
  TYPE(REAL_DATA_LAYOUT), POINTER, INTENT(INOUT)          :: B(:)
  TYPE(COMPLEX_DATA_LAYOUT), DIMENSION(:), INTENT(INOUT)  :: nlB
  TYPE(REAL_DATA_LAYOUT), POINTER, INTENT(INOUT)          :: ScalS(:)
  TYPE(REAL_DATA_LAYOUT), POINTER, INTENT(IN)             :: ScalP(:)
  TYPE(COMPLEX_DATA_LAYOUT), DIMENSION(:), INTENT(INOUT)  :: ScalNL
  LOGICAL, INTENT(IN)                                     :: imhd
  real(WP)                                                :: dt
  LOGICAL                                                 :: res
  ! Local
  !> Loop indices
  INTEGER  :: i

  ! == Initialisation ==
  res = .false.

  ! == Use Euler scheme for this first step ==
  call solve_spec_Euler(U,V,W,Uk,Vk,Wk,nlx,nly,nlz, &
                      & B,nlB,ScalS,ScalP,ScalNL,   &
                      & dt, imhd, Srank, res)
  if(.not. res) then
    write(6,'(a)') '[ERROR] solve_AB2_init: fail to compute first time step with Euler scheme.'
    return
  end if

  ! == Save Non Linear terms for next time-step ==
  if(.not. (permute(nlx,Uk_sav(1)) .and. &
          & permute(nly,Vk_sav(1)) .and. &
          & permute(nlz,Wk_sav(1)))) then
    write(6,'(a)') '[ERROR] solve_AB2_init: fail to save non-linear terms for velocity'
    return
  end if
  do i=1, nbscal
    if(.not. permute(ScalNL(i), Scalk_sav(i,1))) then
        write(6,'(a,i0)') '[ERROR] solve_AB2_init: fail to save scalarS_nl ', i
        return
    end if
  end do
  if(imhd) then
    do i=1, 3
      if(.not. permute(nlB(i),Bk_sav(i,1))) then
        write(6,'(a,i0)') '[ERROR] solve_AB2_init: fail to save nlB ', i
        return
    end if
    end do
  end if

  ! == Save time step ==
  dt_prev(1) = dt

  ! == Use AB2 at the next step ==
  solve_spec_basic => solve_spec_AB2

end subroutine solve_spec_AB2_init


!> Adams-Barshforth (order 3) for spectral sovler (Navier-Stokes equation for velocity and advection-diffusion equation
!! for scalar).
!! @author
!! Mouloud Kessar and Jean-Baptiste Lagaert, LEGI
!!     @param[in,out]  U            = longitudinal velocity (physical space)
!!     @param[in,out]  V            = vertical velocity  (physical space)
!!     @param[in,out]  W            = spanwise velocity  (physical space)
!!     @param[in,out]  Uk           = longitudinal velocity (spectral space)
!!     @param[in,out]  Vk           = vertical velocity  (spectral space)
!!     @param[in,out]  Wk           = spanwise velocity  (spectral space)
!!     @param[in,out]  nlx          = nonlinear terms for U (fourier space), initialy set by firstTerm
!!     @param[in,out]  nly          = nonlinear terms for V (fourier space), initialy set by firstTerm
!!     @param[in,out]  nlz          = nonlinear terms for W (fourier space), initialy set by firstTerm
!!     @param[in]      spec_rank    = my MPI processe rank in the range [0:nbcpus-1]
!!     @param[in]      imhd         = 1 if MHD is active, 0 else
!!     @param[out]     time         = reached time
!!     @param[in, out] B            = magnetic field
!!     @param[out]     res          = true if no error occurs
!! @details
!!        Velocity follows and Navier-Stokes equation and the scalar fields follow a
!!    diffusion-advection equation. This subroutine solve it for a time intervall of
!!    [T ; T+dt]. This function update these fields (by solving the equations)
!!    both spectral and real space.
!!        This is a third order multpile step Adams-Bashforth. See Section 3.3.1 in
!!    [AB] : "Simulation of confined MHD flows usinf a pseudo-spectral method with
!!    volume penalization"', J.A. Morales, M. Leroy, W.J.T. Bos and K. Schneider
!!    to appear in JCP (available in HAL, ref hal-00719737)
!!    If the velocity field in real space is already on the good mesh
!!    (considering the scalar resolution), then solver used it. Else, it compute
!!    this vector on the right resolution directly from the Fourrier space
!!    (and thus perform Fourrier interpolation).
!!        For an alvelius forcing, at the statiscal steady state, the effective
!!    forcing power is, by assuming dt constant and using the notation from [AB]
!!    P_eff = P*(beta_20 + beta_21/2 + beta_22/3) = 1.38*P
subroutine solve_spec_AB3(U,V,W,Uk,Vk,Wk,nlx,nly,nlz, &
                      & B,nlB,ScalS,ScalP,ScalNL,   &
                      & dt, imhd, Srank, res)

  use differential_tools

  ! Input/Output
  INTEGER,                INTENT(IN)                      :: Srank
  TYPE(REAL_DATA_LAYOUT), INTENT(INOUT)                   :: U,V,W
  type(COMPLEX_DATA_LAYOUT), intent(inout)                :: Uk,Vk,Wk
  TYPE(REAL_DATA_LAYOUT)                                  :: sigma
  type(COMPLEX_DATA_LAYOUT)                               :: sigmak
  TYPE(COMPLEX_DATA_LAYOUT),INTENT(INOUT)                 :: nlx, nly, nlz
  TYPE(REAL_DATA_LAYOUT), POINTER, INTENT(INOUT)          :: B(:)
  TYPE(COMPLEX_DATA_LAYOUT), DIMENSION(:), INTENT(INOUT)  :: nlB
  TYPE(REAL_DATA_LAYOUT), POINTER, INTENT(INOUT)          :: ScalS(:)
  TYPE(REAL_DATA_LAYOUT), POINTER, INTENT(IN)             :: ScalP(:)
  TYPE(COMPLEX_DATA_LAYOUT), DIMENSION(:), INTENT(INOUT)  :: ScalNL
  LOGICAL, INTENT(IN)                                     :: imhd
  real(WP)                                                :: dt
  LOGICAL                                                 :: res
  ! Local
  !> Loop indices
  INTEGER  :: i
  !> Ponderation of non-linear term computed at different time-steps
  real(WP)  :: coef1, coef2

  ! == Initialisation ==
  res = .false.

  ! == Compute non-linear term, forcing and projection (div(U,V,W)=0) ==
  call spectral_source_term(U,V,W,Uk,Vk,Wk, &
                          & nlx, nly, nlz,  &
                          & imhd, B, nlB,   &
                          & ScalS,ScalNL,   &
                          & ScalP, dt,      &
                          & Srank, res)
  if(.not.res) then
    write(6,'(a)') '[ERROR] solve_AB3: fail to compute source_term'
    return
  end if

  ! == Save Non Linear terms for next time-step ==
  ! After permutation, we want:
  ! sav(1) -> non-linear at current time step t(n)            , to save
  ! sav(2) -> non-linear at previous time step t(n-1)         , to save
  ! NL     -> non-linear at previous previous time step t(n-2), could be erased after use
  if(.not. (permute(nlx,Uk_sav(1),Uk_sav(2)) .and. &
          & permute(nly,Vk_sav(1),Vk_sav(2)) .and. &
          & permute(nlz,Wk_sav(1),Wk_sav(2)))) then
    write(6,'(a)') '[ERROR] solve_AB3: fail to save non-linear terms for velocity'
    return
  end if
  do i=1, nbscal
    if(.not. permute(ScalNL(i), Scalk_sav(i,1), Scalk_sav(i,2))) then
        write(6,'(a,i0)') '[ERROR] solve_AB3: fail to save scalarS_nl ', i
        return
    end if
  end do
  if(imhd) then
    do i=1, 3
      if(.not. permute(nlB(i),Bk_sav(i,1),Bk_sav(i,2))) then
        write(6,'(a,i0)') '[ERROR] solve_AB3: fail to save nlB ', i
        return
    end if
    end do
  end if

  ! == Update NL terms from ones computed in previous step ==
  ! Compute combination of time step (n-1) and (n-2):
  ! 1 = nl      = time (n-2)
  ! 2 = sav(2)  = time (n-1)
  coef1 = dt*((2._WP*dt)+(3._WP*dt_prev(1))) &
          & /(6._WP*dt_prev(2)*(dt_prev(1)+dt_prev(2)))
  coef2 = -dt*((2._WP*dt)+(3._WP*dt_prev(1))+(3._WP*dt_prev(2))) &
          & /(6._WP*dt_prev(2)*dt_prev(1))
  call compute_AB_NL(nlx, nly, nlz, Uk_sav(2), Vk_sav(2), Wk_sav(2),&
              & imhd,nlB,Bk_sav(1,2),Bk_sav(2,2),Bk_sav(3,2),       &
              & ScalNL, Scalk_sav(:,2), coef1, coef2,               &
              & dt_prev(2)/2._WP, res)
  if(.not.res) then
    write(6,'(a)') '[ERROR] solve_AB3: fail to update non-linear term with times (n-2) and (n-1)'
    return
  end if
  ! Compute combination of time step (n-1) and (n-2):
  ! 1 = NL      = combinaison of times (n-1) and (n-2)
  ! 2 = sav(1)  = current time (n)
  coef1 = 1._WP
  coef2 = ((dt*((2._WP*dt)+(6._WP*dt_prev(1))+(3._WP*dt_prev(2))))  &
          &   + (6._WP*dt_prev(1)*(dt_prev(1)+dt_prev(2))))         &
          & /(6._WP*dt_prev(2)*(dt_prev(1)+dt_prev(2)))
  call compute_AB_NL(nlx, nly, nlz, Uk_sav(1), Vk_sav(1), Wk_sav(1),&
              & imhd,nlB,Bk_sav(1,1),Bk_sav(2,1),Bk_sav(3,1),       &
              & ScalNL, Scalk_sav(:,1), coef1, coef2,               &
              & dt_prev(1)/2._WP, res)
  if(.not.res) then
    write(6,'(a)') '[ERROR] solve_AB3: fail to update non-linear term'
    return
  end if

  ! == compute diffusion and update complex fields ==
  call spectral_update_fields(Uk,Vk,Wk, Uk, &
                     & Vk,Wk,nlx, nly, nlz, &
                     & imhd,Bk,Bk,nlB,      &
                     & ScalArrayk,ScalArrayk,ScalNL,&
                     & dt, res)
  if(.not.res) then
    write(6,'(a)') '[ERROR] solve_AB3: fail to update complex fields'
    return
  end if

  ! == Update real fields ==
  call spectral_update_real_fields(U,V,W,Uk,Vk,Wk, &
                            & imhd, B, Bk,    &
                            & ScalS,ScalArrayk,   &
                            & res)
  if(.not.res) then
    write(6,'(a)') '[ERROR] solve_AB3: fail to update real fields'
    return
  end if

  ! == Update time ==
  sim_time = sim_time + dt
  dt_prev(2) = dt_prev(1)
  dt_prev(1) = dt

end subroutine solve_spec_AB3


!> Adams-Barshforth initialisation for spectral sovler (Navier-Stokes equation for velocity and advection-diffusion equation
!! for scalar): this multi-step method cannot be used for the first step.
!! @author
!! Mouloud Kessar and Jean-Baptiste Lagaert, LEGI
!!     @param[in,out]  U            = longitudinal velocity (physical space)
!!     @param[in,out]  V            = vertical velocity  (physical space)
!!     @param[in,out]  W            = spanwise velocity  (physical space)
!!     @param[in,out]  Uk           = longitudinal velocity (spectral space)
!!     @param[in,out]  Vk           = vertical velocity  (spectral space)
!!     @param[in,out]  Wk           = spanwise velocity  (spectral space)
!!     @param[in,out]  nlx          = nonlinear terms for U (fourier space), initialy set by firstTerm
!!     @param[in,out]  nly          = nonlinear terms for V (fourier space), initialy set by firstTerm
!!     @param[in,out]  nlz          = nonlinear terms for W (fourier space), initialy set by firstTerm
!!     @param[in]      spec_rank    = my MPI processe rank in the range [0:nbcpus-1]
!!     @param[in]      imhd         = 1 if MHD is active, 0 else
!!     @param[out]     time         = reached time
!!     @param[in, out] B            = magnetic field
!!     @param[out]     res          = true if no error occurs
!! @details
!!        Velocity follows and Navier-Stokes equation and the scalar fields follow a
!!    diffusion-advection equation. This subroutine solve it for a time intervall of
!!    [T ; T+dt]. This function update these fields (by solving the equations)
!!    both spectral and real space.
!!        This called an Euler or RK method for the first time step.
!!    If the velocity field in real space is already on the good mesh
!!    (considering the scalar resolution), then solver used it. Else, it compute
!!    this vector on the right resolution directly from the Fourrier space
!!    (and thus perform Fourrier interpolation).
subroutine solve_spec_AB3_init(U,V,W,Uk,Vk,Wk,nlx,nly,nlz,  &
                      & B, nlB, ScalS, ScalP, ScalNL,     &
                      & dt, imhd, Srank, res)

  use differential_tools

  ! Input/Output
  INTEGER,                INTENT(IN)                      :: Srank
  TYPE(REAL_DATA_LAYOUT), INTENT(INOUT)                   :: U,V,W
  type(COMPLEX_DATA_LAYOUT), intent(inout)                :: Uk,Vk,Wk
  TYPE(REAL_DATA_LAYOUT)                                  :: sigma
  type(COMPLEX_DATA_LAYOUT)                               :: sigmak
  TYPE(COMPLEX_DATA_LAYOUT),INTENT(INOUT)                 :: nlx, nly, nlz
  TYPE(REAL_DATA_LAYOUT), POINTER, INTENT(INOUT)          :: B(:)
  TYPE(COMPLEX_DATA_LAYOUT), DIMENSION(:), INTENT(INOUT)  :: nlB
  TYPE(REAL_DATA_LAYOUT), POINTER, INTENT(INOUT)          :: ScalS(:)
  TYPE(REAL_DATA_LAYOUT), POINTER, INTENT(IN)             :: ScalP(:)
  TYPE(COMPLEX_DATA_LAYOUT), DIMENSION(:), INTENT(INOUT)  :: ScalNL
  LOGICAL, INTENT(IN)                                     :: imhd
  real(WP)                                                :: dt
  LOGICAL                                                 :: res
  ! Local
  !> Loop indices
  INTEGER  :: i

  ! == Initialisation ==
  res = .false.

  if (dt_prev(1)<0._WP)  then  ! First time step
    ! -- Use Euler scheme for this first step --
    call solve_spec_Euler(U,V,W,Uk,Vk,Wk,nlx,nly,nlz, &
                        & B,nlB,ScalS,ScalP,ScalNL,   &
                        & dt, imhd, Srank, res)
    if(.not. res) then
      write(6,'(a)') '[ERROR] solve_AB3_init: fail to compute first time step with Euler scheme.'
      return
    end if
    ! -- Save Non Linear terms for next time-step --
    if(.not. (permute(nlx,Uk_sav(1)) .and. &
            & permute(nly,Vk_sav(1)) .and. &
            & permute(nlz,Wk_sav(1)))) then
      write(6,'(a)') '[ERROR] solve_AB3_init: fail to save non-linear terms for velocity'
      return
    end if
    do i=1, nbscal
      if(.not. permute(ScalNL(i), Scalk_sav(i,1))) then
          write(6,'(a,i0)') '[ERROR] solve_AB3_init: fail to save scalarS_nl ', i
          return
      end if
    end do
    if(imhd) then
      do i=1, 3
        if(.not. permute(nlB(i),Bk_sav(i,1))) then
          write(6,'(a,i0)') '[ERROR] solve_AB3_init: fail to save nlB ', i
          return
      end if
      end do
    end if
    ! -- Save time step --
    dt_prev(1) = dt
  else  ! Second time step
    ! -- Save previous time step and previous NL terms --
    ! AB3 will replace sav(1) and dt_prev(1) by the "new"
    ! save at time(n=2)
    dt_prev(2) = dt_prev(1)
    Uk_sav(2)%values = Uk_sav(1)%values
    Vk_sav(2)%values = Vk_sav(1)%values
    Wk_sav(2)%values = Wk_sav(1)%values
    do i=1, nbscal
      Scalk_sav(i,2)%values = Scalk_sav(i,1)%values
    end do
    if(imhd) then
      Bk_sav(1,2)%values = Bk_sav(1,1)%values
      Bk_sav(2,2)%values = Bk_sav(2,1)%values
      Bk_sav(3,2)%values = Bk_sav(3,1)%values
    end if
    ! -- Use AB3 scheme for this second step --
    call solve_spec_AB3(U,V,W,Uk,Vk,Wk,nlx,nly,nlz, &
                        & B,nlB,ScalS,ScalP,ScalNL, &
                        & dt, imhd, Srank, res)
    if(.not. res) then
      write(6,'(a)') '[ERROR] solve_AB3_init: fail to compute first time step with Euler scheme.'
      return
    end if
    ! -- Init done => use AB3 at the next step --
    solve_spec_basic => solve_spec_AB3
  end if

end subroutine solve_spec_AB3_init


!------------------------------------------------------------------------------
!> This function calculate the time step (dt) value.
!> @author
!> Patrick BEGOU, LEGI
!
!> @details
!!     This function compute the right time step for each component : velocity (time
!! step is limited by the CFL), scalar compute with pseudo-spectral methods (CFL
!! again) and scalar solved by particle/spectral methods (dt is on limited by the
!! velocity gradient)
!!     @param[in]   U           = longitudinal velocity (physical space)
!!     @param[in]   V           = vertical velocity  (physical space)
!!     @param[in]   W           = spanwise velocity  (physical space)
!!     @param[in]   B           = magnetic fields  (vector in physical space)
!!     @param[in]   imhd        = true if magnetic fields is present
!!     @param[in]   ScalArray   = an array of scalars (physical space).
!!     @param[in]   ScalArray   = an array of scalars (physical space) - for particle method
!!     @param[out]  dt_velo     = maximal time step allowed to compute velocity
!!     @param[out]  dt_spec     = maximal time step allowed to compute scalar with
!!                                 full pseudo-spectral scheme and for MHD
!!     @param[out]  dt_spec     = maximal time step allowed to compute scalar with
!!                                 mixed particle/spectral scheme.
!!     @param[in]   spec_rank   = my MPI processe rank in the range [0:nbcpus-1]
!------------------------------------------------------------------------------
SUBROUTINE get_timestep(U,V,W,B,imhd,ScalArray,ScalArrayP,dt_velo, dt_spec, dt_part,spec_rank)
  IMPLICIT NONE

  TYPE(REAL_DATA_LAYOUT), INTENT(IN)              :: U,V,W
  TYPE(REAL_DATA_LAYOUT), INTENT(IN),DIMENSION(:) :: B
  TYPE(REAL_DATA_LAYOUT), POINTER, INTENT(IN)     :: ScalArray(:)
  TYPE(REAL_DATA_LAYOUT), POINTER, INTENT(IN)     :: ScalArrayP(:)
  LOGICAL, INTENT(IN)                             :: imhd
  INTEGER, INTENT(IN)                             :: spec_rank
  REAL(WP), INTENT(OUT)                           :: dt_velo, dt_spec, dt_part

  REAL(WP) :: dt_visq, dt_lor, dx,dt_induc

  real(WP) :: maxVel,minVel
  real(WP) :: maxDxU, maxDyV, maxDzW
  real(WP) :: maxGradVel
  !real(WP) :: maxB2,maxB,minB
  real(WP) :: minLorTemp,MinLor,minInducTemp,MinInduc

  integer  :: maxgridx, maxgridy, maxgridz
  !> Constant <= 1 to ensure not breaking the stability conditions
  real(WP) :: c_conv, c_vis, c_lor, c_part, dt_max
  !> Loop indice
  integer  :: i, j, k
  !> To check FFT result
  logical  :: res

  CALL parser_read('Fixed time step',dt_velo)
  if (dt_velo > 0.0_WP) then
    dt_spec = dt_velo
    dt_part = dt_velo
    return
  end if

  CALL parser_read('Maximum time step', dt_max)
  IF (parser_is_defined('Fixed dt - spec')) THEN
    CALL parser_read('Fixed dt - spec',dt_velo)
    dt_spec = dt_velo
  ELSE
    !get the absolute higher velocity
    maxVel=maxval(abs(U%values) + abs(V%values) + abs(W%values))
    maxVel=DoGlobalMax(spec_rank,maxVel)
    IF (maxvel.LE.0.0_wp) THEN
      IF (spec_rank.EQ.0) WRITE(6,'(a,f14.6)')'[ERROR] get_timestep: The maximum velocity value is wrong: ',maxVel
      STOP
      RETURN
    ENDIF

    ! Compute dt_velo
    !get the higher grid points numbernumber
    maxgridx=U%nx
    maxgridy=U%ny
    maxgridz=U%nz
    dx = min(U%Lx/maxgridx, U%Ly/maxgridy, U%Lz/maxgridz)
    dt_velo = dx / maxVel

    ! Compute dt_spec
    IF (ASSOCIATED(ScalArray)) THEN
      ! get the higher grid points numbernumber
      DO i=1, SIZE(ScalArray)
        maxgridx=max(maxgridx,ScalArray(i)%nx)
        maxgridy=max(maxgridy,ScalArray(i)%ny)
        maxgridz=max(maxgridz,ScalArray(i)%nz)
      ENDDO
      dx = min(U%Lx/maxgridx, U%Ly/maxgridy, U%Lz/maxgridz)
      dt_spec = dx / maxVel
    ELSE
      dt_spec = dt_velo
    ENDIF

    dt_visq = dx*dx / mu
    CALL parser_read('Convective time step number', c_conv)
    CALL parser_read('Viscous time step number', c_vis)

    IF (spec_rank.eq.0) THEN
      WRITE(6,'(2(a,g15.8))') '    [INFO-TIME-STEP] Convective time: ', &
      & min(dt_velo, dt_spec), ' Viscous time: ', dt_visq
    ENDIF
    !dt = min( c_conv*dt_conv, c_vis*dt_visq )
    ! Viscous time step is not a limit
    dt_velo = min( c_conv*dt_velo, dt_max)
    dt_spec = min( c_conv*dt_spec, dt_max)


    IF (ASSOCIATED(ScalArray)) THEN
      DO i=1, SIZE(ScalArray)
        !dt = min(dt,c_vis*schmidt(i)*dx*dx / mu)
        IF (spec_rank.eq.0) WRITE(6,'(a,i0,a,g15.8)') &
          & '    [INFO-TIME-STEP] Viscous time for scalar ',i,' : ',schmidt(i)*dx*dx / mu
      ENDDO
    ENDIF
! 
!     IF (imhd) THEN
!       !maxB2_U = maxval( ( B(1)%values**2 + &
!       !          &         B(2)%values**2 + &
!       !          &         B(3)%values**2 ) / &
!       !          & sqrt(   U%values**2    + &
!       !          &         V%values**2    + &
!       !          &         W%values**2    ) )
!       !maxB2_U = DoGlobalMax(spec_rank,maxB2_U)
!       !!!!! Lorentz time
!       !maxB2=maxval(abs(B(1)%values)**2 + abs(B(2)%values)**2 + abs(B(3)%values)**2 +&
!       !                &abs(B(1)%values)*abs(B(2)%values) + abs(B(1)%values)*abs(B(3)%values) +&
!       !                &abs(B(3)%values)*abs(B(2)%values))
!       !maxB2=DoGlobalMax(spec_rank,maxB2)
!       !minVel=minval((abs(U%values) + abs(V%values) + abs(W%values)))
!       !minVel=DoGlobalMax(spec_rank,minVel)
! 
!       minLorTemp=minval(abs((U%values/(B(1)%values)+V%values/(B(2)%values)+W%values/(B(3)%values))/B(1)%values))
!       minLor=DoGlobalMin(spec_rank,minLorTemp)
! 
!       minLorTemp=minval(abs((U%values/(B(1)%values)+V%values/(B(2)%values)+W%values/(B(3)%values))/B(2)%values))
!       minLorTemp=DoGlobalMin(spec_rank,minLorTemp)
! 
!       minLor=min(minLor,minLorTemp)
! 
!       minLorTemp=minval(abs((U%values/(B(1)%values)+V%values/(B(2)%values)+W%values/(B(3)%values))/B(3)%values))
!       minLorTemp=DoGlobalMin(spec_rank,minLorTemp)
! 
!       minLor=min(minLor,minLorTemp)
! 
!       IF (minLor.LE.0.0_wp) THEN
!         IF (spec_rank.EQ.0) WRITE(6,'(a,f14.6)')'[ERROR] get_timestep: The minimum velocity value is wrong: ',minVel
!         STOP
!         RETURN
!       ENDIF
! 
!       dt_lor = dx*minLor !dx/maxB2_U
!       CALL parser_read('Lorentz time step number', c_lor)
!       dt_spec = min( dt_spec, c_lor*dt_lor )
!       IF (spec_rank.eq.0) WRITE(6,'(2(a,g14.8))') '    [INFO-TIME-STEP] Lorentz force time:',dt_lor, '    dt_lor / dt ', dt_lor/dt_spec
! 
!       !!!!Induction time
!       !maxB=maxval(abs(B(1)%values) + abs(B(2)%values) + abs(B(3)%values))
!       !maxB=DoGlobalMax(spec_rank,maxB)
! 
!       !minB=maxval(abs(B(1)%values) + abs(B(2)%values) + abs(B(3)%values))
!       ! minB=DoGlobalMax(spec_rank,minB)
!       minInducTemp=minval(abs(((B(1)%values/U%values)+(B(2)%values/V%values)+(B(3)%values/W%values))/B(1)%values))
!       minInduc=DoGlobalMin(spec_rank,minInducTemp)
! 
!       minInducTemp=minval(abs(((B(1)%values/U%values)+(B(2)%values/V%values)+(B(3)%values/W%values))/B(2)%values))
!       minInducTemp=DoGlobalMin(spec_rank,minInducTemp)
! 
!       minInduc=min(minInduc,minInducTemp)
! 
!       minInducTemp=minval(abs((U%values/(B(1)%values)+V%values/(B(2)%values)+W%values/(B(3)%values))/B(3)%values))
!       minInducTemp=DoGlobalMin(spec_rank,minInducTemp)
!
!       minInduc=min(minInduc,minInducTemp)
!
!       dt_induc = dx*minInduc !dx/maxB2_U
!       CALL parser_read('Lorentz time step number', c_lor)
!       dt_spec = min( dt_spec, c_lor*dt_induc )
!       IF (spec_rank.eq.0) WRITE(6,'(2(a,g14.8))') '    [INFO-TIME-STEP] induction time:',dt_induc, '    dt_induc / dt ', c_lor*dt_induc/dt_spec
!     END IF ! mhd
  END IF ! fixed time step for spectral solver (velocity AND scalar)

  if(nbscal_part >0) then
    ! CFL for 1D explicit diffusion
    dx = min(U%Lx/ScalArrayP(1)%nx, U%Ly/ScalArrayP(1)%ny, U%Lz/ScalArrayP(1)%nz)
    if (spec_rank.eq.0) write(6,'(a,g15.8,a)') &
      & '    [INFO-TIME-STEP] Viscous time for scalarP (CFL 1D)  : ',minval(schmidt_part)*(dx**2) / (2._WP*mu),' (subcycle => dt/2< dt_diff)'
    dt_part = 0
    ! Stability condition for advection with particle methods.
    ! time step < 1/[2*bl_size*norm_inf{gradient(velocity)}]
    ! Due to the splitting, we have time step = dt/2, thus we get :
    ! dt <  1/[bl_size*norm_inf{gradient(veelocity)}]
    ! Compute max(dxU,dyV,dzW)
    ! Derivation in spectral space
    do k=Uk_mean%zmin,Uk_mean%zmax
      do j=Uk_mean%ymin,Uk_mean%ymax
        do i=Uk_mean%xmin,Uk_mean%xmax
          Uk_mean%values(i,j,k) = (0.0_WP,1.0_WP)*VelWN%kx(i)*Uk%values(i,j,k)
          Vk_mean%values(i,j,k) = (0.0_WP,1.0_WP)*VelWN%ky(j)*Vk%values(i,j,k)
          Wk_mean%values(i,j,k) = (0.0_WP,1.0_WP)*VelWN%kz(k)*Wk%values(i,j,k)
        enddo
      enddo
    enddo
    ! Go back to real space
    call btran(Uk_mean,U_mean,res)
    if (.not.res) return
    call btran(Vk_mean,V_mean,res)
    if (.not.res) return
    call btran(Wk_mean,W_mean,res)
    if (.not.res) return
    ! And take the max
    maxDxU = maxval(U_mean%values)
    maxDyV = maxval(V_mean%values)
    maxDzW = maxval(W_mean%values)
    maxGradVel = max(maxDxU,maxDyV,maxDzW)
    maxGradVel = DoGlobalMax(spec_rank, maxGradVel)
    if (maxGradVel>0) then
      call parser_read('Constant for particle time step', c_part)
      c_part = min(c_part,1.0_WP)
      dt_part = c_part*part_stab_coeff/maxGradVel
    else
      dt_part = dt_velo
    end if
      dt_part = min(dt_part, dt_max)
      dt_velo = min(dt_velo, dt_part)
      dt_spec = min(dt_spec, dt_part)
  else    ! nb_scal > 0
    dt_part = dt_velo
  end if  ! nb_scal > 0

END SUBROUTINE get_timestep


!-------------------------------------------------------------------------------
!> @author
!! Jean-Baptiste Lagaert
!
!> Compute all non-diffusive term for velocity, magnetic field and scalars
!! solved with pseudo-spectral method
!! @details
!! Scal_part is only used to compute residence time
!! @param[in]     U     = velocity along X
!! @param[in]     V     = velocity along Y
!! @param[in]     W     = velocity along Z
!! @param[in]     Uk    = spectral velocity along X
!! @param[in]     Vk    = spectral velocity along Y
!! @param[in]     Wk    = spectral velocity along Z
!! @param[in,out] nlx   = (spectral) non linear term for U
!! @param[in,out] nly   = (spectral) non linear term for V
!! @param[in,out] nlz   = (spectral) non linear term for W
!! @param[in]     mhd   = to active or not magnetic field
!! @param[in]     B     = magnetic field
!! @param[in,out] nlB   = (spectral) non linear term for B
!! @param[in]     ScalS = scalar (array) solved with pseudo-spectral method
!! @param[in,out] ScalNL= (spectral) non linear term for Scal
!! @param[in]     ScalP = scalar (array) solved with particles method - only
!!                        used for update residence time solved with spectral method
!! @param[in]     step        = time step
!! @param[in]     Srank = rank inside spectral MPI communicator
!! @param[in]     res   = true if succeed
!-------------------------------------------------------------------------------
SUBROUTINE spectral_source_term(U,V,W,Uk,Vk,Wk, &
                              & nlx, nly, nlz,  &
                              & mhd, B, nlB,    &
                              & ScalS,ScalNL,   &
                              & ScalP,          &
                              & step, Srank, res&
                              & , save_forc)

  ! Input/Output
  type(REAL_DATA_LAYOUT), intent(in)                      :: U,V,W
  type(COMPLEX_DATA_LAYOUT), intent(in)                   :: Uk, Vk, Wk
  type(COMPLEX_DATA_LAYOUT),intent(inout)                 :: nlx, nly, nlz
  logical, intent(in)                                     :: mhd
  type(REAL_DATA_LAYOUT), pointer, intent(in)             :: B(:)
  type(COMPLEX_DATA_LAYOUT), dimension(:), intent(inout)  :: nlB
  type(REAL_DATA_LAYOUT),intent(in),dimension(:)          :: ScalS, ScalP
  type(COMPLEX_DATA_LAYOUT), intent(inout), dimension(:)  :: ScalNL
  real(WP), intent(in)                                    :: step
  integer,                intent(in)                      :: Srank
  logical , intent(out)                                   :: res
  logical, intent(in), optional                           :: save_forc
  ! Local
  integer  :: sca  ! Loop on scalars

  ! == Init ==
  res = .false.

  ! == compute non linear term ==
  ! For velocities
  call non_linear_vel(U,V,W,nlx,nly,nlz,Srank,res)
  if (.not. res) then
    write(6,'(a)') '[ERROR] solver: non_linear_vel'
    return
  end if
  ! For scalar if we use a pseudo-spectral method
  call non_linear_scal(U,V,W,ScalS,ScalNL,ScalP,Srank,res)
  if (.not. res) then
    write(6,'(a)') '[ERROR] solver: non_linear_scal'
    return
  end if
  ! For MHD if needed
  if (mhd) then
    call non_linear_mhd(U,V,W,B,nlB,Srank,res)
    if (.not. res) then
      write(6,'(a)') '[ERROR] solver: non_linear_mhd'
      return
    end if
    IF ( iles .eq. 1 ) THEN
        DivTub(1)%values=0.0_WP
        DivTub(2)%values=0.0_WP
        DivTub(3)%values=0.0_WP
        CALL computeBfieldModel(B(1),B(2),B(3),Bk(1),Bk(2),Bk(3),U,V,W,Uk,Vk,Wk,DivTub(1),DivTub(2),DivTub(3),BfieldWN,Srank,res)
    ENDIF

    !Coupling with NS equation : Lorentz Force computation
    call mhd_velocity_coupling(B,nlx,nly,nlz,Srank,res)
    if (.not. res) then
      write(6,'(a)') '[ERROR] solver: mhd_velocity_coupling'
      return
    end if
  end if

  ! == Spectral forcing ==
  ! For velocity
  if (iforce==1) then
!   if (present(save_forc)) then
!     Fk_sav(1)%values = 0._WP
!     Fk_sav(2)%values = 0._WP
!     Fk_sav(3)%values = 0._WP
!     call force_compute(Uk,Vk,Wk,VelWN,Fk_sav(1),Fk_sav(2),Fk_sav(3),step,Srank,U,V,W)
!     Uk%values = Uk%values + Fk_sav(1)%values
!     Vk%values = Vk%values + Fk_sav(2)%values
!     Wk%values = Wk%values + Fk_sav(3)%values
!   else
      call force_compute(Uk,Vk,Wk,VelWN,nlx,nly,nlz,step,Srank,U,V,W)
!   end if
  end if
  ! For scalars
  do sca=1,nbscal
    if (amplitude_sf(sca).ne.0) call forcing_spec_scal(ScalWN(sca),ScalNL(sca), &
                                      & amplitude_sf(sca),kmin_sf(sca),kmax_sf(sca))
  end do

  ! == Projection step ==
  if(.not.NullifyDiv(nlx, nly, nlz, VelWN, Srank)) then
    write(6,'(a)') '[ERROR] solver_step: in projection'
    res = .false.
    return
  endif


END SUBROUTINE spectral_source_term

!------------------------------------------------------------------------------
!> @author
!> Guillaume Balarac and Mouloud KESSAR, LEGI
!
!>
!> @details
!> This function calculate the Lorentz force to add at the non linear velocity term.
!> The Lorentz force term is :
!> (jxB)_i = d/dxj( BiBj)    [B^2/2 delta_ij is in the pressure term]

!> @param [in] B Magnetic field  (physical space)
!> @param [in,out] nlx array for nonlinear terms for U (fourier space)
!> @param [in,out] nly array for nonlinear terms for V (fourier space)
!> @param [in,out] nlz array for nonlinear terms for W (fourier space)
!> @param [in] spec_rank my MPI processe rank in the range [0:nbcpus-1]
!> @param [out] res .TRUE. if the calculation is successfull.
!------------------------------------------------------------------------------
SUBROUTINE mhd_velocity_coupling(B,nlx,nly,nlz,spec_rank,res)

    use dealias_tools

    TYPE(REAL_DATA_LAYOUT), INTENT(IN), DIMENSION(:):: B
    TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT)        :: nlx, nly, nlz
    INTEGER, INTENT(IN)                             :: spec_rank
    LOGICAL, INTENT(OUT)                            :: res
    !
    TYPE(REAL_DATA_LAYOUT), DIMENSION(3)    :: vecReal
    TYPE(COMPLEX_DATA_LAYOUT)               :: bufComp
    TYPE(COMPLEX_DATA_LAYOUT), DIMENSION(3) :: vecComp
    TYPE(COMPLEX_DATA_LAYOUT)               :: nik


    IF (.NOT. samelayout(Uk,Vk,Wk,Bk(1))) THEN
        WRITE(6,'(a,i0)')'[ERROR] different meshing for velocity and magnetic fields not yet implemented',spec_rank
        RETURN
    ELSE
         IF ((.NOT. copyStructOnly(B(1),vecReal(1))) .OR. &
                &(.NOT. copyStructOnly(B(1),vecReal(2))) .OR. &
                &(.NOT. copyStructOnly(B(1),vecReal(3))) .OR. &
                &(.NOT. copyStructOnly(Bk(1),vecComp(1))) .OR. &
                &(.NOT. copyStructOnly(Bk(1),vecComp(2))) .OR. &
                &(.NOT. copyStructOnly(Bk(1),vecComp(3))) .OR. &
                &(.NOT. copyStructOnly(Bk(1),bufComp)) .OR. &
                &(.NOT. copyStructOnly(Bk(1),nik)) )RETURN

        ! x direction
        vecReal(1)%values = ( B(1)%values + Bxext) * ( B(1)%values + Bxext)
        vecReal(2)%values = ( B(1)%values + Bxext) * ( B(2)%values + Byext)
        vecReal(3)%values = ( B(1)%values + Bxext) * ( B(3)%values + Bzext)
        CALL ftran(vecReal(1),vecComp(1),res)
        CALL ftran(vecReal(2),vecComp(2),res)
        CALL ftran(vecReal(3),vecComp(3),res)
        ! for z direction computation
        bufComp%values = vecComp(3)%values
        CALL computeDivergenceK(BfieldWN,vecComp(1),vecComp(2),vecComp(3),nik)
        nlx%values = nlx%values + nik%values

        ! y direction - vecComp(2) in x direction is vecComp(1) in y direction
        vecReal(2)%values = ( B(2)%values + Byext) * ( B(2)%values + Byext)
        vecReal(3)%values = ( B(2)%values + Byext) * ( B(3)%values + Bzext)
        vecComp(1)%values = vecComp(2)%values
        CALL ftran(vecReal(2),vecComp(2),res)
        CALL ftran(vecReal(3),vecComp(3),res)
        CALL computeDivergenceK(BfieldWN,vecComp(1),vecComp(2),vecComp(3),nik)
        nly%values = nly%values + nik%values

        ! z direction - vecComp(3) in y direction is vecComp(2) in z direction AND vecComp(3)%values in x direction is vecComp(1) in z direction
        vecReal(3)%values = ( B(3)%values + Bzext) * ( B(3)%values + Bzext)
        vecComp(1)%values = bufComp%values
        vecComp(2)%values = vecComp(3)%values
        CALL ftran(vecReal(3),vecComp(3),res)
        CALL computeDivergenceK(BfieldWN,vecComp(1),vecComp(2),vecComp(3),nik)
        nlz%values = nlz%values + nik%values
       ! insert LES models for velocities
       IF ( iles .eq. 1 ) THEN
        CALL computeLorentzVelModel(B(1),B(2),B(3),Bk(1),Bk(2),Bk(3),nlx,nly,nlz,VelWN,spec_rank,res)
       ENDIF
        ! dealiasing
        SELECT CASE (get_dealiasMethod(spec_rank))
        !CASE(0)
        !    IF(spec_rank.EQ.0) PRINT*,'No dealiasing for velocities'
        !    CONTINUE
        CASE(1)
                ! IF(spec_rank.EQ.0) PRINT*,'Isotropic dealiasing for velocities'
                CALL dealiasIsotrope(nlx,VelWN)
                CALL dealiasIsotrope(nly,VelWN)
                CALL dealiasIsotrope(nlz,VelWN)
        CASE(2)
                ! IF(spec_rank.EQ.0) PRINT*,'Anisotropic dealiasing for velocities'
                CALL dealiasAnisotrope(nlx,VelWN)
                CALL dealiasAnisotrope(nly,VelWN)
                CALL dealiasAnisotrope(nlz,VelWN)
        END SELECT

    ENDIF

    CALL deleteDataLayout(nik)
    CALL deleteDataLayout(vecComp(1))
    CALL deleteDataLayout(vecComp(2))
    CALL deleteDataLayout(vecComp(3))
    CALL deleteDataLayout(bufComp)
    CALL deleteDataLayout(vecReal(1))
    CALL deleteDataLayout(vecReal(2))
    CALL deleteDataLayout(vecReal(3))

    res=.TRUE.


END SUBROUTINE mhd_velocity_coupling

!------------------------------------------------------------------------------
!> @author
!> Guillaume Balarac and Mouloud KESSAR, LEGI
!
!>
!> @details
!> This function calculate the non linear terms for mhd field.
!> The non linear term are :
!> nlbx = - (v.bz - w.by)
!> nlby = - (w.bx - u.bz)
!> nlbz = - (u.by - v.bx)
!> the opposite of these definitions are taken because the terms are used directly in the rhs

!> @param [in] U longitudinal velocity (physical space)
!> @param [in] V vertical velocity (physical space)
!> @param [in] W spanwise velocity (physical space)
!> @param [in] B Magnetic field (physical space)
!> @param [in,out] nlB array for nonlinear terms for B field (fourier space), initialy set by firstTerm
!> @param [in] spec_rank my MPI processe rank in the range [0:nbcpus-1]
!> @param [out] res .TRUE. if the calculation is successfull.
!------------------------------------------------------------------------------
SUBROUTINE non_linear_mhd(U,V,W,B,nlB,spec_rank,res)

    use dealias_tools

    TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: U,V,W
    TYPE(REAL_DATA_LAYOUT), INTENT(IN), DIMENSION(:) :: B
    TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT), DIMENSION(:) :: nlB
    INTEGER, INTENT(IN)    :: spec_rank
    LOGICAL, INTENT(OUT) :: res

    TYPE(REAL_DATA_LAYOUT)        :: nis


    IF (.NOT. samelayout(U,V,W,B(1))) THEN
        WRITE(6,'(a,i0)')'[ERROR] different meshing for velocity and magnetic fields not yet implemented',spec_rank
        RETURN
    ENDIF

    IF (.NOT. copyStructOnly(B(1),nis)) RETURN

    !! pour nlBx
    nis%values = V%values*(B(3)%values + Bzext) - W%values*(B(2)%values + Byext)
    CALL ftran(nis,nlB(1),res)

    !! pour nlBy
    nis%values = W%values*(B(1)%values + Bxext) - U%values*(B(3)%values + Bzext)
    CALL ftran(nis,nlB(2),res)

    !! pour nlBz
    nis%values = U%values*(B(2)%values + Byext) - V%values*(B(1)%values + Bxext)
    CALL ftran(nis,nlB(3),res)

    CALL deleteDataLayout(nis)

!        IF ( iles .eq. 1 ) THEN
!         CALL computeBfieldModel(B(1),B(2),B(3),Bk(1),Bk(2),Bk(3),U,V,W,Uk,Vk,Wk,nlB(1),nlB(2),nlB(3),BfieldWN,spec_rank,res)
!        ENDIF

    ! dealiasing as the same method used for velocity
    SELECT CASE (get_dealiasMethod(spec_rank))
    !CASE(0)
    !    IF(spec_rank.EQ.0) PRINT*,'No dealiasing for MHD ',sca
    !    CONTINUE
    CASE(1)
        !IF(spec_rank.EQ.0) PRINT*,'Isotrope dealiasing for MHD ',sca
        CALL dealiasIsotrope(nlB(1),BfieldWN)
        CALL dealiasIsotrope(nlB(2),BfieldWN)
        CALL dealiasIsotrope(nlB(3),BfieldWN)
    CASE(2)
        !IF(spec_rank.EQ.0) PRINT*,'Anisotrope dealiasing for MHD ',sca
        CALL dealiasAnisotrope(nlB(1),BfieldWN)
        CALL dealiasAnisotrope(nlB(2),BfieldWN)
        CALL dealiasAnisotrope(nlB(3),BfieldWN)
    END SELECT

    res=.TRUE.

END SUBROUTINE non_linear_mhd

!------------------------------------------------------------------------------
!> @author
!> Guillaume Balarac and Patrick BEGOU, LEGI
!
!>
!> @details
!> This function calculate the non linear terms for velocities.
!> @param [in] U longitudinal velocity (physical space)
!> @param [in] V vertical velocity (physical space)
!> @param [in] W spanwise velocity (physical space)
!> @param [in] spec_rank my MPI processe rank in the range [0:nbcpus-1]
!> @param [in,out] nlx array for nonlinear terms for U (fourier space), initialy set by firstTerm
!> @param [in,out] nly array for nonlinear terms for V (fourier space), initialy set by firstTerm
!> @param [in,out] nlz array for nonlinear terms for W (fourier space), initialy set by firstTerm
!> @param [out] res .TRUE. if the calculation is successfull.
!------------------------------------------------------------------------------
SUBROUTINE non_linear_vel(U,V,W,nlx,nly,nlz,spec_rank,res)

    use dealias_tools

    TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: U,V,W
    TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT) :: nlx, nly, nlz
    INTEGER, INTENT(IN)    :: spec_rank
    LOGICAL, INTENT(OUT) :: res

    TYPE(REAL_DATA_LAYOUT)    :: ddx, nis
    TYPE(COMPLEX_DATA_LAYOUT) :: ddxk


    res=.FALSE.

    ! Allocate working array ddx, ddxk, and nis
    ! We can use copyStructOnly to create the structure,
    ! should be faster than init_data_layout
    ! Warning, it do not care with datas and sets values(:,:,:) to 0.0
    IF (    .NOT. copyStructOnly(U,ddx)     .OR. &
        &   .NOT. copyStructOnly(U,nis)     .OR. &
        &   .NOT. copyStructOnly(Uk,ddxk)) RETURN

    ! Call firstTerm
    ! XXX Old version
!     CALL firstTerm(U,V,W,U,nlx,VelWN,ddx,ddxk,spec_rank,res)
!     IF(.NOT. res) RETURN
!     CALL firstTerm(U,V,W,V,nly,VelWN,ddx,ddxk,spec_rank,res)
!     IF(.NOT. res) RETURN
!     CALL firstTerm(U,V,W,W,nlz,VelWN,ddx,ddxk,spec_rank,res)
!     IF(.NOT. res) RETURN
  ! New version
  CALL firstTerm(U,V,W,nlx,nly,nlz,VelWN,ddx,ddxk,spec_rank,res)
  IF(.NOT. res) RETURN

    ! No more call second term: it is probably no more presice than only first
    ! term and it induces error on the average off velocity component.
    !!Call second term and sum nlx? for each Uk, Vk, Wk
    CALL secondTerm(U,V,W,Uk,ddxk,VelWN,ddx,nis,spec_rank,res)
    IF(.NOT. res) RETURN
    nlx%values=0.5*(nlx%values+ddxk%values)
    CALL secondTerm(U,V,W,Vk,ddxk,VelWN,ddx,nis,spec_rank,res)
    IF(.NOT. res) RETURN
    nly%values=0.5*(nly%values+ddxk%values)
    CALL secondTerm(U,V,W,Wk,ddxk,VelWN,ddx,nis,spec_rank,res)
    IF(.NOT. res) RETURN
    nlz%values=0.5*(nlz%values+ddxk%values)

    ! insert LES models for velocities
    IF ( iles .eq. 1 ) THEN
        CALL computeVelocityModel(U,V,W,Uk,Vk,Wk,Bk(1),Bk(2),Bk(3),nlx,nly,nlz,VelWN,spec_rank,res)
        IF (.NOT. res ) RETURN
    ENDIF

    ! dealiasing
    SELECT CASE (get_dealiasMethod(spec_rank))
    !CASE(0)
    !    IF(spec_rank.EQ.0) PRINT*,'No dealiasing for velocities'
    !    CONTINUE
    CASE(1)
        ! IF(spec_rank.EQ.0) PRINT*,'Isotropic dealiasing for velocities'
        CALL dealiasIsotrope(nlx,VelWN)
        CALL dealiasIsotrope(nly,VelWN)
        CALL dealiasIsotrope(nlz,VelWN)
    CASE(2)
        ! IF(spec_rank.EQ.0) PRINT*,'Anisotropic dealiasing for velocities'
        CALL dealiasAnisotrope(nlx,VelWN)
        CALL dealiasAnisotrope(nly,VelWN)
        CALL dealiasAnisotrope(nlz,VelWN)
    END SELECT

    nlx%values = -nlx%values
    nly%values = -nly%values
    nlz%values = -nlz%values

    !free module Working arrays
    CALL deleteDataLayout(ddx)
    CALL deleteDataLayout(nis)
    CALL deleteDataLayout(ddxk)

    res=.TRUE.
END SUBROUTINE non_linear_vel

!------------------------------------------------------------------------------
!> @author
!> Guillaume Balarac and Patrick BEGOU, LEGI
!>
!> DESCRIPTION
!> This function calculate the non linear terms for scalars.
!> @details
!> This function calculate the non linear terms for scalars. The prerequisit is that U,V,W and ScalArray
!> (physical space) are coherent with Uk, Vk, Wk and ScalArrayk (fourier space) wich are
!> internal to the solver module in the current time-step.
!>
!> @param [in] U longitudinal velocity (physical space)
!> @param [in] V vertical velocity (physical space)
!> @param [in] W spanwise velocity (physical space)
!> @param [in] ScalArray a pointer to an array of scalars (physical space).
!> @param [in] ScalArray a pointer to an array of scalars solved with spec/part
!! method (needed for some special source terms)(physical space).
!> @param [in,out] ScalArrayNL a pointer to an array for scalars non linear terms (fourier space).
!> @param [in] spec_rank my MPI processe rank in the range [0:nbcpus-1]
!> @param [out] res .TRUE. if the calculation is successfull.
!------------------------------------------------------------------------------

SUBROUTINE non_linear_scal(U,V,W,ScalArray,ScalArrayNL,Scal_partArray,spec_rank,res)

    use forcing_scalar
    use dealias_tools

    IMPLICIT NONE
    TYPE(REAL_DATA_LAYOUT), INTENT(IN)                      :: U,V,W
    TYPE(REAL_DATA_LAYOUT),INTENT(IN),DIMENSION(:)          :: ScalArray, Scal_partArray
    TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT), DIMENSION(:)  :: ScalArrayNL
    INTEGER, INTENT(IN)                                     :: spec_rank
    LOGICAL, INTENT(OUT)                                    :: res

    TYPE(REAL_DATA_LAYOUT)      :: ddx, nis, Uscal, Vscal, Wscal
    TYPE(COMPLEX_DATA_LAYOUT)   :: ddxk
    TYPE(COMPLEX_DATA_LAYOUT)   :: UscalK,VscalK,WscalK
    INTEGER :: sca, sca2

    res=.FALSE.

    DO sca=1,nbscal
        IF (.NOT. copyStructOnly(ScalArray(sca),ddx) .OR. &
        & .NOT. copyStructOnly(ScalArray(sca),nis) .OR. &
        & .NOT. copyStructOnly(ScalArrayk(sca),ddxk)) RETURN

        IF (.NOT. samelayout(U,V,W,ScalArray(sca))) THEN
            IF (.NOT. copyStructOnly(ScalArray(sca),Uscal) .OR. &
                & .NOT. copyStructOnly(ScalArray(sca),Vscal) .OR. &
                & .NOT. copyStructOnly(ScalArray(sca),Wscal)) RETURN
            IF ((ilessc(sca)==1).or.(sc_spec_forc(sca)==1)) THEN
                IF (.NOT. copyStructOnly(ScalArrayNL(sca),UscalK) .OR. &
                    & .NOT. copyStructOnly(ScalArrayNL(sca),VscalK) .OR. &
                    & .NOT. copyStructOnly(ScalArrayNL(sca),WscalK)) RETURN
                IF (.NOT. discretization (Uk,ScalArray(sca),Uscal,UscalK) .OR. &
                    & .NOT. discretization (Vk,ScalArray(sca),Vscal,VscalK) .OR. &
                    & .NOT. discretization (Wk,ScalArray(sca),Wscal,WscalK)) RETURN
            ELSE
                IF (.NOT. transdiscretization (Uk,ScalArray(sca),Uscal) .OR. &
                    & .NOT. transdiscretization (Vk,ScalArray(sca),Vscal) .OR. &
                    & .NOT. transdiscretization (Wk,ScalArray(sca),Wscal)) RETURN
            ENDIF

            CALL firstTerm(Uscal,Vscal,Wscal,ScalArray(sca),ScalArrayNL(sca),ScalWN(sca),ddx,ddxk,spec_rank,res)
            IF(.NOT. res) RETURN
            !ScalArrayNL(sca)%values= - ScalArrayNL(sca)%values

            CALL secondTerm(Uscal,Vscal,Wscal,ScalArrayk(sca),ddxk,ScalWN(sca),ddx,nis,spec_rank,res)
            IF(.NOT. res) RETURN
            ScalArrayNL(sca)%values= - 0.5_WP*(ScalArrayNL(sca)%values+ddxk%values)

            ! Forcing
            IF(sc_spec_forc(sca)==1) call forcing_scal_mean(ScalArrayNL(sca), VscalK, spec_rank)

            IF ( ilessc(sca)==1 ) THEN
                CALL computeScalarModel(sca,Uscal,Vscal,Wscal,UscalK,VscalK,WscalK,ScalArray(sca),ScalArrayk(sca),&
                    & ScalArrayNL(sca),ScalWN(sca),spec_rank,res)
                IF (.NOT. res ) RETURN
            ENDIF

            IF ((ilessc(sca)==1).or.(sc_spec_forc(sca)==1)) THEN
                CALL deleteDataLayout(UscalK)
                CALL deleteDataLayout(VscalK)
                CALL deleteDataLayout(WscalK)
            ENDIF

            ! For residence time
            ! -- Source term if it is a residence term --
            if((res_t(sca)>0).and.(res_t(sca)<=nbscal)) then
                ! The current scalar a residence time for a spectral scalar
                if (samelayout(ScalArray(sca),scalArray(res_t(sca)))) then
                    where ((ScalArray(res_t(sca))%values<1.-mix_threshold).and. &
                        & (ScalArray(res_t(sca))%values>mix_threshold))
                        ScalArrayNL(sca)%values = ScalArrayNL(sca)%values + 1.0_WP
                    end where
                else
                    ! Change discretisation
                    if (.not. transdiscretization (scalArrayK(res_t(sca)),Scalarray(sca),Uscal)) return
                    where ((Uscal%values<1.-mix_threshold).and. &
                        & (Uscal%values>mix_threshold))
                        ScalArrayNL(sca)%values = ScalArrayNL(sca)%values + 1.0_WP
                    end where
                end if
            else if (res_t(sca)>0) then
                ! The current scalar a residence time for a spectral scalar
                sca2 = res_t(sca) - nbscal
                if (samelayout(ScalArray(sca),scal_partArray(sca2))) then
                    where ((Scal_partArray(sca2)%values<1.-mix_threshold).and. &
                        & (Scal_partArray(sca2)%values>mix_threshold))
                        ScalArrayNL(sca)%values = ScalArrayNL(sca)%values + 1.0_WP
                    end where
                else
                    ! Change discretisation
                    if (.not. transdiscretization (scal_partArrayK(sca2),Scalarray(sca),Uscal)) return
                    where ((Uscal%values<1.-mix_threshold).and. &
                        & (Uscal%values>mix_threshold))
                        ScalArrayNL(sca)%values = ScalArrayNL(sca)%values + 1.0_WP
                    end where
                end if
            end if

            CALL deleteDataLayout(Uscal)
            CALL deleteDataLayout(Vscal)
            CALL deleteDataLayout(Wscal)



        ELSE
            CALL firstTerm(U,V,W,ScalArray(sca),ScalArrayNL(sca),ScalWN(sca),ddx,ddxk,spec_rank,res)
            IF(.NOT. res) RETURN
            !ScalArrayNL(sca)%values= - ScalArrayNL(sca)%values

            CALL secondTerm(U,V,W,ScalArrayk(sca),ddxk,ScalWN(sca),ddx,nis,spec_rank,res)
            IF(.NOT. res) RETURN
            ScalArrayNL(sca)%values= - 0.500000_WP*(ScalArrayNL(sca)%values+ddxk%values)

            ! Todo add forcing : ScalArrayNL is in complex space ..
            IF(sc_spec_forc(sca)==1) call forcing_scal_mean(ScalArrayNL(sca), Vk, spec_rank)

            IF ( ilessc(sca) .eq. 1 ) THEN
                CALL computeScalarModel(sca,U,V,W,Uk,Vk,Wk,ScalArray(sca),ScalArrayk(sca),&
                                 &ScalArrayNL(sca),ScalWN(sca),spec_rank,res)
                IF (.NOT. res ) RETURN
            ENDIF

            ! For residence time
            ! -- Source term if it is a residence term --
            if((res_t(sca)>0).and.(res_t(sca)<=nbscal)) then
                ! The current scalar a residence time for a spectral scalar
                if (samelayout(ScalArray(sca),scalArray(res_t(sca)))) then
                    where ((ScalArray(res_t(sca))%values<1.-mix_threshold).and. &
                        & (ScalArray(res_t(sca))%values>mix_threshold))
                        ScalArrayNL(sca)%values = ScalArrayNL(sca)%values + 1.0_WP
                    end where
                else
                    ! Change discretisation
                    if (.not. copystructonly(scalarray(sca),Uscal)) return
                    if (.not. transdiscretization (scalArrayK(res_t(sca)),Scalarray(sca),Uscal)) return
                    where ((Uscal%values<1.-mix_threshold).and. &
                        & (Uscal%values>mix_threshold))
                        ScalArrayNL(sca)%values = ScalArrayNL(sca)%values + 1.0_WP
                    end where
                    call deleteDataLayout(Uscal)
                end if
            else if (res_t(sca)>0) then
                ! The current scalar a residence time for a spectral scalar
                sca2 = res_t(sca) - nbscal
                if (samelayout(ScalArray(sca),scal_partArray(sca2))) then
                    where ((Scal_partArray(sca2)%values<1.-mix_threshold).and. &
                        & (Scal_partArray(sca2)%values>mix_threshold))
                        ScalArrayNL(sca)%values = ScalArrayNL(sca)%values + 1.0_WP
                    end where
                else
                    ! Change discretisation
                    if (.not. copystructonly(scalarray(sca),Uscal)) return
                    if (.not. transdiscretization (scal_partArrayK(sca2),Scalarray(sca),Uscal)) return
                    where ((Uscal%values<1.-mix_threshold).and. &
                        & (Uscal%values>mix_threshold))
                        ScalArrayNL(sca)%values = ScalArrayNL(sca)%values + 1.0_WP
                    end where
                    call deleteDataLayout(Uscal)
                end if
            end if

        ENDIF

        ! dealiasing
        SELECT CASE (get_dealiasMethod(spec_rank,sca))
        !CASE(0)
        !    !IF(spec_rank.EQ.0) PRINT*,'No dealiasing for scalar ',sca
        !    CONTINUE
        CASE(1)
            !IF(spec_rank.EQ.0) PRINT*,'Isotrope dealiasing for scalar ',sca
            CALL dealiasIsotrope(ScalArrayNL(sca),ScalWN(sca))
        CASE(2)
            !IF(spec_rank.EQ.0) PRINT*,'Anisotrope dealiasing for scalar ',sca
            CALL dealiasAnisotrope(ScalArrayNL(sca),ScalWN(sca))
        END SELECT

        !free module Working arrays
        CALL deleteDataLayout(ddx)
        CALL deleteDataLayout(nis)
        CALL deleteDataLayout(ddxk)

    ENDDO

    res=.TRUE.

END SUBROUTINE non_linear_scal

!------------------------------------------------------------------------------
!> @author
!> Guillaume Balarac and Patrick BEGOU, LEGI
!
!>
!> @details
!> Compute the first term in the skew symmetric form
!> this version only need 2 work arrays: ddx and ddxk instead of 6:
!> ddx is also used for ddy and ddz,
!> ddxk is also used for ddyk and ddzk.
!> @param [in] U longitudinal velocity (physical space)
!> @param [in] V vertical velocity (physical space)
!> @param [in] W spanwise velocity (physical space)
!> @param [in] cur U, V or W velocity or scalar in physical space
!> @param [in,out] nlu2 The calculated term (fourier space)
!> @param [in] wavenum the WaveNumbers associated to cur variable
!> @param [in,out] ddx work array (physical space)
!> @param [in,out] ddxk work array (fourier space)
!> @param [in] spec_rank my MPI processe rank in the range [0:nbcpus-1]
!> @param [out] res .TRUE. if all is succcessfull, .FALSE. elsewhere.
!------------------------------------------------------------------------------
SUBROUTINE firstTerm_scal(U,V,W,cur,nlu2,wavenum,ddx,ddxk,spec_rank,res)
  IMPLICIT NONE
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: U,V,W, cur
  TYPE(REAL_DATA_LAYOUT), INTENT(INOUT) :: ddx
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT) :: nlu2,ddxk
  TYPE(WaveNumbers), INTENT(IN) ::wavenum
  INTEGER, INTENT(IN)  :: spec_rank
  LOGICAL, INTENT(OUT) :: res

  INTEGER:: i,j,k
  COMPLEX(WP), PARAMETER :: ii=(0.0_WP,1.0_WP)


  res=.FALSE.

  IF (.NOT. (samelayout(U,V,W,ddx,cur) &
       &    .AND.samelayout(nlu2,ddxk))) THEN
    WRITE(6,'(a,i0)')'[ERROR]firstTerm: internal error, some array do not match on ',spec_rank
    RETURN
  ENDIF

  ddx%values = U%values*cur%values
  CALL ftran(ddx,ddxk,res)
  IF (.NOT.res) RETURN
  DO k = nlu2%zmin,nlu2%zmax
    DO j = nlu2%ymin,nlu2%ymax
      DO i = nlu2%xmin,nlu2%xmax
        nlu2%values(i,j,k) = wavenum%kx(i)*ddxk%values(i,j,k)
      ENDDO
    ENDDO
  ENDDO

  ddx%values = V%values*cur%values
  CALL ftran(ddx,ddxk,res)
  IF (.NOT.res) RETURN
  DO k = nlu2%zmin,nlu2%zmax
    DO j = nlu2%ymin,nlu2%ymax
      DO i = nlu2%xmin,nlu2%xmax
          nlu2%values(i,j,k) = nlu2%values(i,j,k) + wavenum%ky(j)*ddxk%values(i,j,k)
      ENDDO
    ENDDO
  ENDDO

  ddx%values = W%values*cur%values
  CALL ftran(ddx,ddxk,res)
  IF (.NOT.res) RETURN
  DO k = nlu2%zmin,nlu2%zmax
    DO j = nlu2%ymin,nlu2%ymax
      DO i = nlu2%xmin,nlu2%xmax
        nlu2%values(i,j,k) = nlu2%values(i,j,k) + wavenum%kz(k)*ddxk%values(i,j,k)
      ENDDO
    ENDDO
  ENDDO

  nlu2%values=ii*nlu2%values

  res=.TRUE.
END SUBROUTINE firstTerm_scal


!------------------------------------------------------------------------------
!> @author
!> Guillaume Balarac, JB Lagaert and Patrick BEGOU, LEGI
!
!>
!> @details
!! Compute the first term in the skew symmetric form for a vector
!! this version only need 2 work arrays: ddx and ddxk instead of 6 (they are re-used for other direction)
!! and avoid to compute twice the "non-diagonal terms like UV=VU as done if
!! "first_term_scal" is used for a vector.
!! @param [in] U longitudinal velocity (physical space)
!! @param [in] V vertical velocity (physical space)
!! @param [in] W spanwise velocity (physical space)
!! @param [in,out] nlu2 The calculated term (fourier space)
!! @param [in] wavenum the WaveNumbers associated to cur variable
!! @param [in,out] ddx work array (physical space)
!! @param [in,out] ddxk work array (fourier space)
!! @param [in] spec_rank my MPI processe rank in the range [0:nbcpus-1]
!! @param [out] res .TRUE. if all is succcessfull, .FALSE. elsewhere.
!------------------------------------------------------------------------------
SUBROUTINE firstTerm_vect(U,V,W,nlU,nlV,nlW,wavenum,ddx,ddxk,spec_rank,res)
  IMPLICIT NONE
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: U,V,W
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT) :: nlU, nlV, nlW
  TYPE(REAL_DATA_LAYOUT), INTENT(INOUT) :: ddx
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT) :: ddxk
  TYPE(WaveNumbers), INTENT(IN) ::wavenum
  INTEGER, INTENT(IN)  :: spec_rank
  LOGICAL, INTENT(OUT) :: res

  INTEGER:: i,j,k
  COMPLEX(WP), PARAMETER :: ii=(0.0_WP,1.0_WP)


  res=.FALSE.
! this version only need 2 work arrays: ddx and ddxk instead of 6
! ddx  is also used for ddy  and ddz
! ddxk is also used for ddyk and ddzk

  IF (.NOT. (samelayout(U,V,W,ddx) &
       &    .AND.samelayout(nlU,nlV,nlW,ddxk))) THEN
    WRITE(6,'(a,i0)')'[ERROR]firstTerm: internal error, some array do not match on ',spec_rank
    RETURN
  ENDIF

  ddx%values = U%values*U%values
  CALL ftran(ddx,ddxk,res)
  IF (.NOT.res) RETURN
  DO k = nlU%zmin,nlU%zmax
    DO j = nlU%ymin,nlU%ymax
      DO i = nlU%xmin,nlU%xmax
          nlU%values(i,j,k) = wavenum%kx(i)*ddxk%values(i,j,k)
       ENDDO
    ENDDO
  ENDDO

  ddx%values = V%values*U%values
  CALL ftran(ddx,ddxk,res)
  IF (.NOT.res) RETURN
  DO k = nlU%zmin,nlU%zmax
    DO j = nlU%ymin,nlU%ymax
      DO i = nlU%xmin,nlU%xmax
          nlU%values(i,j,k) = nlU%values(i,j,k) + wavenum%ky(j)*ddxk%values(i,j,k)
          nlV%values(i,j,k) =                     wavenum%kx(i)*ddxk%values(i,j,k)
       ENDDO
    ENDDO
  ENDDO

  ddx%values = W%values*U%values
  CALL ftran(ddx,ddxk,res)
  IF (.NOT.res) RETURN
  DO k = nlU%zmin,nlU%zmax
    DO j = nlU%ymin,nlU%ymax
      DO i = nlU%xmin,nlU%xmax
          nlU%values(i,j,k) = nlU%values(i,j,k) + wavenum%kz(k)*ddxk%values(i,j,k)
          nlW%values(i,j,k) =                     wavenum%kx(i)*ddxk%values(i,j,k)
       ENDDO
    ENDDO
  ENDDO

  ddx%values = V%values*V%values
  CALL ftran(ddx,ddxk,res)
  IF (.NOT.res) RETURN
  DO k = nlU%zmin,nlU%zmax
    DO j = nlU%ymin,nlU%ymax
      DO i = nlU%xmin,nlU%xmax
          nlV%values(i,j,k) = nlV%values(i,j,k) + wavenum%ky(j)*ddxk%values(i,j,k)
       ENDDO
    ENDDO
  ENDDO

  ddx%values = W%values*V%values
  CALL ftran(ddx,ddxk,res)
  IF (.NOT.res) RETURN
  DO k = nlU%zmin,nlU%zmax
    DO j = nlU%ymin,nlU%ymax
      DO i = nlU%xmin,nlU%xmax
          nlV%values(i,j,k) = nlV%values(i,j,k) + wavenum%kz(k)*ddxk%values(i,j,k)
          nlW%values(i,j,k) = nlW%values(i,j,k) + wavenum%ky(j)*ddxk%values(i,j,k)
       ENDDO
    ENDDO
  ENDDO

  ddx%values = W%values*W%values
  CALL ftran(ddx,ddxk,res)
  IF (.NOT.res) RETURN
  DO k = nlU%zmin,nlU%zmax
    DO j = nlU%ymin,nlU%ymax
      DO i = nlU%xmin,nlU%xmax
          nlW%values(i,j,k) = nlW%values(i,j,k) + wavenum%kz(k)*ddxk%values(i,j,k)
       ENDDO
    ENDDO
  ENDDO

  nlU%values=ii*nlU%values
  nlV%values=ii*nlV%values
  nlW%values=ii*nlW%values

  res=.TRUE.
END SUBROUTINE firstTerm_vect


!------------------------------------------------------------------------------
!> @author
!> Guillaume Balarac and Patrick BEGOU, LEGI
!
!>
!> @details
!> Compute the second term in the skew symmetric form
!> @param [in] U longitudinal velocity (physical space)
!> @param [in] V vertical velocity  (physical space)
!> @param [in] W spanwise velocity  (physical space)
!> @param [in] curk is Uk, Vk or Wk velocity or scalar in fourier space
!> @param [in,out] nlx The calculated term  (fourier space)
!> @param [in] wavenum the WaveNumbers associated to cur variable
!> @param [in,out] ddx work array (physical space)
!> @param [in,out] nis work array (physical space)
!> @param [in,out] ddxk work array (fourier space)
!> @param [in] spec_rank my MPI processe rank in the range [0:nbcpus-1]
!> @param [out] res .TRUE. if all is succcessfull, .FALSE. elsewhere.
!> @param [out] maxDeriv = maximum value of the gradient of cur (ie the norm inf of the gradient) - optional.
!------------------------------------------------------------------------------
SUBROUTINE secondTerm(U,V,W,curk,ddxk,wavenum,ddx,nis,spec_rank,res, maxDeriv)

  IMPLICIT NONE
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: U,V,W
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: curk
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT) :: ddxk
  TYPE(WaveNumbers), INTENT(IN) ::wavenum
  TYPE(REAL_DATA_LAYOUT), INTENT(INOUT) :: ddx,nis
  INTEGER, INTENT(IN)  :: spec_rank
  LOGICAL, INTENT(OUT) :: res
  REAL(WP), DIMENSION(4), INTENT(OUT), OPTIONAL :: maxDeriv

  INTEGER:: i,j,k
  COMPLEX(WP), PARAMETER :: ii=(0.0_WP,1.0_WP)

  res=.FALSE.
  if(present(maxDeriv)) maxDeriv(4) = 0

  IF (.NOT. (samelayout(U,V,W,ddx,nis) &
      & .AND.samelayout(curk,ddxk))) THEN
    WRITE(6,'(a,i0)')'[ERROR]secondTerm: internal error, some array do not match on ',spec_rank
    RETURN
  ENDIF

  DO k = ddxk%zmin,ddxk%zmax
    DO j =  ddxk%ymin,ddxk%ymax
      DO i =  ddxk%xmin,ddxk%xmax
        ddxk%values(i,j,k) = ii*wavenum%kx(i)*curk%values(i,j,k)
      ENDDO
    ENDDO
  ENDDO
  CALL btran(ddxk,ddx,res)
  IF (.NOT.res) RETURN

  IF(present(maxDeriv)) THEN
    IF(maxDeriv(1)==1) maxDeriv(4) = max(maxDeriv(4),maxval(abs(ddx%values)))
  END IF

  !general case
  nis%values = U%values*ddx%values

  DO k = ddxk%zmin,ddxk%zmax
    DO j =  ddxk%ymin,ddxk%ymax
      DO i =  ddxk%xmin,ddxk%xmax
        ddxk%values(i,j,k) = ii*wavenum%ky(j)*curk%values(i,j,k)
      ENDDO
    ENDDO
  ENDDO
  CALL btran(ddxk,ddx,res)
  IF (.NOT.res) RETURN

  IF(present(maxDeriv)) THEN
    IF(maxDeriv(2)==1) maxDeriv(4) = max(maxDeriv(4),maxval(abs(ddx%values)))
  END IF

  !general case
  nis%values = nis%values + V%values*ddx%values

  DO k = ddxk%zmin,ddxk%zmax
    DO j =  ddxk%ymin,ddxk%ymax
      DO i =  ddxk%xmin,ddxk%xmax
        ddxk%values(i,j,k) = ii*wavenum%kz(k)*curk%values(i,j,k)
      ENDDO
    ENDDO
  ENDDO
  CALL btran(ddxk,ddx,res)
  IF (.NOT.res) RETURN

  IF(present(maxDeriv)) THEN
    IF(maxDeriv(3)==1) maxDeriv(4) = max(maxDeriv(4),maxval(abs(ddx%values)))
  END IF

  !general case
  nis%values = nis%values + W%values*ddx%values

  CALL ftran(nis,ddxk,res)
  RETURN

END SUBROUTINE secondTerm


!-------------------------------------------------------------------------------
!> Update fields (velocity, B, scalarS) for one time step with the input "source term".
!! @author
!! Jean-Baptiste Lagaert, Patrick Begou, LEGI
!! @details
!!    Solve the diffusion equation for velocity scalar and B with a given source
!! term. For B, the actual source_term is rot(nlB)
!! @param[in,out] Uk          = new spectral velocity along X
!! @param[in,out] Vk          = new spectral velocity along Y
!! @param[in,out] Wk          = new spectral velocity along Z
!! @param[in]     Uk_old      = previous spectral velocity along X
!! @param[in]     Vk_old      = previous spectral velocity along Y
!! @param[in]     Wk_old      = previous spectral velocity along Z
!! @param[in]     nlx         = (spectral) non linear term for U
!! @param[in]     nly         = (spectral) non linear term for V
!! @param[in]     nlz         = (spectral) non linear term for W
!! @param[in]     mhd         = to active or not magnetic field
!! @param[in,out] Bk          = new magnetic field
!! @param[in]     Bk_old      = previous magnetic field
!! @param[in]     nlB         = (spectral) non linear term for B
!! @param[in,out] ScalSk      = scalar (array) solved with pseudo-spectral method
!! @param[in]     ScalSk_old  = scalar (array) solved with pseudo-spectral method
!! @param[in,out] ScalNL      = (spectral) non linear term for Scal
!!                              used for update residence time solved with spectral method
!! @param[in]     step        = time step
!! @param[in]     res         = true if succeed
!-------------------------------------------------------------------------------
subroutine spectral_update_fields_basic(Uk,Vk,Wk, Uk_old,     &
                       & Vk_old, Wk_old, nlx, nly, nlz, &
                       & mhd,Bk,Bk_old,nlB,             &
                       & ScalSk, ScalSk_old, ScalNL,    &
                       & step, res)

  ! Input/Output
  type(COMPLEX_DATA_LAYOUT), intent(in)                   :: Uk_old,Vk_old,Wk_old
  type(COMPLEX_DATA_LAYOUT), intent(inout)                :: Uk, Vk, Wk
  type(COMPLEX_DATA_LAYOUT),intent(in)                    :: nlx, nly, nlz
  logical, intent(in)                                     :: mhd
  type(COMPLEX_DATA_LAYOUT), intent(inout)                :: Bk(:)
  type(COMPLEX_DATA_LAYOUT), dimension(:), intent(in)     :: nlB, Bk_old
  type(COMPLEX_DATA_LAYOUT),intent(in),dimension(:)       :: ScalSk_old, ScalNL
  type(COMPLEX_DATA_LAYOUT), intent(inout), dimension(:)  :: ScalSk
  real(WP), intent(in)                                    :: step
  logical , intent(out)                                   :: res

  ! local
  real(WP)  :: coef   ! to compute exponential term for diffusion in Fourrier space.
  real(WP)  :: kk, ekk! wave number norm and exponential of (coeff*kk)
  real(WP)  :: kkz    ! = square of third componnent of wave numbers
  real(WP)  :: kkyz   ! = sum of square of third and second componnent of wave numbers ! kk = kkyz + (WN_x)**2
  integer   :: i,j,k  ! loop indices
  integer   :: sca    ! loop indice
  !> For MHD
  complex(WP) :: Axk,Ayk,Azk

  ! == Init ==
  res = .false.
  coef = -mu*step

  ! == For velocty ==
  do k = Uk%zmin,Uk%zmax
    kkz = VelWN%kz(k)*VelWN%kz(k)
    do j = Uk%ymin,Uk%ymax
      kkyz = kkz + VelWN%ky(j)*VelWN%ky(j)
      do i = Uk%xmin,Uk%xmax
        kk = kkyz + VelWN%kx(i)*VelWN%kx(i)
        ekk=exp(coef*kk)
        Uk%values(i,j,k) = (Uk_old%values(i,j,k) + step*nlx%values(i,j,k))*ekk
        Vk%values(i,j,k) = (Vk_old%values(i,j,k) + step*nly%values(i,j,k))*ekk
        Wk%values(i,j,k) = (Wk_old%values(i,j,k) + step*nlz%values(i,j,k))*ekk
      end do
    end do
  end do

  ! == For scalars ==
  do sca=1,nbscal
    do k = ScalArrayK(sca)%zmin,ScalArrayK(sca)%zmax
      kkz = ScalWN(sca)%kz(k)*ScalWN(sca)%kz(k)
      do j = ScalArrayK(sca)%ymin,ScalArrayK(sca)%ymax
        kkyz = kkz + ScalWN(sca)%ky(j)*ScalWN(sca)%ky(j)
        do i = ScalArrayK(sca)%xmin,ScalArrayK(sca)%xmax
          kk = ScalWN(sca)%kx(i)*ScalWN(sca)%kx(i) + kkyz
          ekk=exp(kk*coef/schmidt(sca))
          ScalSk(sca)%values(i,j,k) = (ScalSk_old(sca)%values(i,j,k) + step*ScalNL(sca)%values(i,j,k))*ekk
        end do
      end do
    end do
  end do

  ! == For MHD if needed ==
  if (mhd) then
    do k = Bk(1)%zmin,Bk(1)%zmax
      kkz = BfieldWN%kz(k)*BfieldWN%kz(k)
      do j = Bk(1)%ymin,Bk(1)%ymax
        kkyz = kkz + BfieldWN%ky(j)*BfieldWN%ky(j)
        do i = Bk(1)%xmin,Bk(1)%xmax
          kk = BfieldWN%kx(i)*BfieldWN%kx(i) + kkyz
          ekk=exp(kk*coef/Prm)
          Axk = ii*BfieldWN%ky(j)*nlB(3)%values(i,j,k) - ii*BfieldWN%kz(k)*nlB(2)%values(i,j,k)
          Ayk = ii*BfieldWN%kz(k)*nlB(1)%values(i,j,k) - ii*BfieldWN%kx(i)*nlB(3)%values(i,j,k)
          Azk = ii*BfieldWN%kx(i)*nlB(2)%values(i,j,k) - ii*BfieldWN%ky(j)*nlB(1)%values(i,j,k)
          Bk(1)%values(i,j,k) = (Bk_old(1)%values(i,j,k) + step*Axk)*ekk
          Bk(2)%values(i,j,k) = (Bk_old(2)%values(i,j,k) + step*Ayk)*ekk
          Bk(3)%values(i,j,k) = (Bk_old(3)%values(i,j,k) + step*Azk)*ekk
        end do
      end do
    end do
  end if ! mhd

    res = .true.

end subroutine spectral_update_fields_basic


!-------------------------------------------------------------------------------
!> Update fields (velocity, B, scalarS) for one time step with the input "source term".
!! Variante 1 for RK scheme, with two different time scale.
!! @author
!! Jean-Baptiste Lagaert, Patrick Begou, LEGI
!! @details
!!    Solve the diffusion equation for velocity scalar and B with a given source
!! term. For B, the actual source_term is rot(nlB)
!! @param[in,out] Uk          = new spectral velocity along X
!! @param[in,out] Vk          = new spectral velocity along Y
!! @param[in,out] Wk          = new spectral velocity along Z
!! @param[in]     Uk_old      = previous spectral velocity along X
!! @param[in]     Vk_old      = previous spectral velocity along Y
!! @param[in]     Wk_old      = previous spectral velocity along Z
!! @param[in]     nlx         = (spectral) non linear term for U
!! @param[in]     nly         = (spectral) non linear term for V
!! @param[in]     nlz         = (spectral) non linear term for W
!! @param[in]     mhd         = to active or not magnetic field
!! @param[in,out] Bk          = new magnetic field
!! @param[in]     Bk_old      = previous magnetic field
!! @param[in]     nlB         = (spectral) non linear term for B
!! @param[in,out] ScalSk      = scalar (array) solved with pseudo-spectral method
!! @param[in]     ScalSk_old  = scalar (array) solved with pseudo-spectral method
!! @param[in,out] ScalNL      = (spectral) non linear term for Scal
!!                              used for update residence time solved with spectral method
!! @param[in]     step1       = time step
!! @param[in]     step2       = second time scale
!! @param[in]     res         = true if succeed
!-------------------------------------------------------------------------------
subroutine spectral_update_fields_two_dt(Uk,Vk,Wk, Uk_old,     &
                       & Vk_old, Wk_old, nlx, nly, nlz, &
                       & mhd,Bk,Bk_old,nlB,             &
                       & ScalSk, ScalSk_old, ScalNL,    &
                       & step1, step2, res)

  ! Input/Output
  type(COMPLEX_DATA_LAYOUT), intent(in)                   :: Uk_old,Vk_old,Wk_old
  type(COMPLEX_DATA_LAYOUT), intent(inout)                :: Uk, Vk, Wk
  type(COMPLEX_DATA_LAYOUT),intent(in)                    :: nlx, nly, nlz
  logical, intent(in)                                     :: mhd
  type(COMPLEX_DATA_LAYOUT), intent(inout)                :: Bk(:)
  type(COMPLEX_DATA_LAYOUT), dimension(:), intent(in)     :: nlB, Bk_old
  type(COMPLEX_DATA_LAYOUT),intent(in),dimension(:)       :: ScalSk_old, ScalNL
  type(COMPLEX_DATA_LAYOUT), intent(inout), dimension(:)  :: ScalSk
  real(WP), intent(in)                                    :: step1, step2
  logical , intent(out)                                   :: res

  ! local
  real(WP)  :: coef   ! to compute exponential term for diffusion in Fourrier space.
  real(WP)  :: coef_bis! to compute exponential term for diffusion in Fourrier space.
  real(WP)  :: kk, ekk! wave number norm and exponential of (coeff*kk)
  real(WP)  :: ekk_bis! wave number exponential of (coeff_bis*kk)
  real(WP)  :: kkz    ! = square of third componnent of wave numbers
  real(WP)  :: kkyz   ! = sum of square of third and second componnent of wave numbers ! kk = kkyz + (WN_x)**2
  integer   :: i,j,k  ! loop indices
  integer   :: sca    ! loop indice
  !> For MHD
  complex(WP) :: Axk,Ayk,Azk

  ! == Init ==
  res = .false.
  coef = -mu*step1
  coef_bis = mu*step2

  ! == For velocty ==
  do k = Uk%zmin,Uk%zmax
    kkz = VelWN%kz(k)*VelWN%kz(k)
    do j = Uk%ymin,Uk%ymax
      kkyz = kkz + VelWN%ky(j)*VelWN%ky(j)
      do i = Uk%xmin,Uk%xmax
        kk = kkyz + VelWN%kx(i)*VelWN%kx(i)
        ekk=exp(coef*kk)
        ekk_bis=exp(coef_bis*kk)
        Uk%values(i,j,k) = (Uk_old%values(i,j,k) + step1*nlx%values(i,j,k)*ekk_bis)*ekk
        Vk%values(i,j,k) = (Vk_old%values(i,j,k) + step1*nly%values(i,j,k)*ekk_bis)*ekk
        Wk%values(i,j,k) = (Wk_old%values(i,j,k) + step1*nlz%values(i,j,k)*ekk_bis)*ekk
      end do
    end do
  end do

  ! == For scalars ==
  do sca=1,nbscal
    do k = ScalArrayK(sca)%zmin,ScalArrayK(sca)%zmax
      kkz = ScalWN(sca)%kz(k)*ScalWN(sca)%kz(k)
      do j = ScalArrayK(sca)%ymin,ScalArrayK(sca)%ymax
        kkyz = kkz + ScalWN(sca)%ky(j)*ScalWN(sca)%ky(j)
        do i = ScalArrayK(sca)%xmin,ScalArrayK(sca)%xmax
          kk = ScalWN(sca)%kx(i)*ScalWN(sca)%kx(i) + kkyz
          ekk=exp(kk*coef/schmidt(sca))
          ekk_bis=exp(kk*coef_bis/schmidt(sca))
          ScalSk(sca)%values(i,j,k) = (ScalSk_old(sca)%values(i,j,k) &
            & + step1*ScalNL(sca)%values(i,j,k)*ekk_bis)*ekk
        end do
      end do
    end do
  end do

  ! == For MHD if needed ==
  if (mhd) then
    do k = Bk(1)%zmin,Bk(1)%zmax
      kkz = BfieldWN%kz(k)*BfieldWN%kz(k)
      do j = Bk(1)%ymin,Bk(1)%ymax
        kkyz = kkz + BfieldWN%ky(j)*BfieldWN%ky(j)
        do i = Bk(1)%xmin,Bk(1)%xmax
          kk = BfieldWN%kx(i)*BfieldWN%kx(i) + kkyz
          ekk=exp(kk*coef/Prm)
          ekk_bis=exp(kk*coef_bis/Prm)
          Axk = ii*BfieldWN%ky(j)*nlB(3)%values(i,j,k) - ii*BfieldWN%kz(k)*nlB(2)%values(i,j,k)
          Ayk = ii*BfieldWN%kz(k)*nlB(1)%values(i,j,k) - ii*BfieldWN%kx(i)*nlB(3)%values(i,j,k)
          Azk = ii*BfieldWN%kx(i)*nlB(2)%values(i,j,k) - ii*BfieldWN%ky(j)*nlB(1)%values(i,j,k)
          Bk(1)%values(i,j,k) = (Bk_old(1)%values(i,j,k) + step1*Axk*ekk_bis)*ekk
          Bk(2)%values(i,j,k) = (Bk_old(2)%values(i,j,k) + step1*Ayk*ekk_bis)*ekk
          Bk(3)%values(i,j,k) = (Bk_old(3)%values(i,j,k) + step1*Azk*ekk_bis)*ekk
        end do
      end do
    end do
  end if ! mhd

  res = .true.

end subroutine spectral_update_fields_two_dt


!-------------------------------------------------------------------------------
!> @author
!! Jean-Baptiste Lagaert
!
!> Update real fields (solved with spectral method) from complex fields via FFT
!! @param[in,out] U     = velocity along X
!! @param[in,out] V     = velocity along Y
!! @param[in,out] W     = velocity along Z
!! @param[in]     Uk    = spectral velocity along X
!! @param[in]     Vk    = spectral velocity along Y
!! @param[in]     Wk    = spectral velocity along Z
!! @param[in,out] B     = magnetic field
!! @param[in]     Bk    = spectral magnetic field
!! @param[in,out] ScalS = scalar (array) solved with pseudo-spectral method
!! @param[in]     ScalSk= spectral scalar (array) solved with pseudo-spectral method
!! @param[in]     Srank = rank inside spectral MPI communicator
!! @param[in]     res   = true if succeed
!-------------------------------------------------------------------------------
subroutine spectral_update_real_fields(U,V,W,Uk,Vk,Wk, &
                              & mhd, B, Bk,     &
                              & ScalS,ScalSk,   &
                              & res)

  ! Input/Output
  type(REAL_DATA_LAYOUT), intent(inout)               :: U,V,W
  type(COMPLEX_DATA_LAYOUT), intent(in)               :: Uk, Vk, Wk
  logical, intent(in)                                 :: mhd
  type(REAL_DATA_LAYOUT), pointer, intent(inout)      :: B(:)
  type(COMPLEX_DATA_LAYOUT), dimension(:), intent(in) :: Bk
  type(REAL_DATA_LAYOUT),intent(inout),dimension(:)   :: ScalS
  type(COMPLEX_DATA_LAYOUT), intent(in), dimension(:) :: ScalSk
  logical , intent(out)                               :: res
  ! Local
  integer   :: sca ! for loop among the scalar fields

  res = .false.

  ! For velocities
  CALL btran(Uk,U,res)
  IF (.NOT.res) RETURN
  CALL btran(Vk,V,res)
  IF (.NOT.res) RETURN
  CALL btran(Wk,W,res)
  IF (.NOT.res) RETURN
  ! For scalaras
  DO sca=1,nbscal
    CALL btran(ScalSk(sca),ScalS(sca),res)
    IF (.NOT.res) RETURN
  ENDDO
  ! For MHD if needed
  IF (mhd) THEN
    CALL btran(Bk(1),B(1),res)
    IF (.NOT.res) RETURN
    CALL btran(Bk(2),B(2),res)
    IF (.NOT.res) RETURN
    CALL btran(Bk(3),B(3),res)
    IF (.NOT.res) RETURN
  ENDIF

  res = .true.

end subroutine spectral_update_real_fields


!-------------------------------------------------------------------------------
!> Update NL fields (velocity, B, scalarS) for AB scheme with one computed on
!! current step and the previous steps.
!! @author
!! Jean-Baptiste Lagaert, LEGI
!! @details
!! Get as input the two non-linear terms nl1 and nl2 and the right ponderation
!! and then give, as output, nl1 = coef1*nl1 + coef2*exp(-mu*step*kk)*nl2.
!! nl1 is an 'input and output' and nl2 is input (read-only)
!! @param[in,out] nlU1        = (spectral) non linear term 1 for U
!! @param[in,out] nlV1        = (spectral) non linear term 1 for V
!! @param[in,out] nlW1        = (spectral) non linear term 1 for W
!! @param[in]     nlU2        = (spectral) non linear term 2 for U, not modified
!! @param[in]     nlV2        = (spectral) non linear term 2 for V, not modified
!! @param[in]     nlW2        = (spectral) non linear term 2 for W, not modified
!! @param[in]     mhd         = to active or not magnetic field
!! @param[in,out] nlB1        = non linear magnetic field 1
!! @param[in]     nlBx2       = non linear magnetic field 2 along X
!! @param[in]     nlBy2       = non linear magnetic field 2 along Y
!! @param[in]     nlBz2       = non linear magnetic field 2 along Z
!! @param[in,out] ScalNL      = (spectral) non linear term 1 for Scal
!! @param[in]     ScalNL2     = (spectral) non linear term 2 for Scal
!! @param[in]     coef1       = weight for NL1
!! @param[in]     coef2       = weight for NL2
!! @param[in]     step        = time scale for exponential term with NL2
!! @param[out]    res         = true if succeed
!-------------------------------------------------------------------------------
subroutine compute_AB_NL(nlU1, nlV1, nlW1, nlU2, nlV2, nlW2,  &
              & mhd,nlB1,nlBx2,nlBy2,nlBz2, ScalNL1,ScalNL2, &
              & coef1, coef2, step, res)

  type(COMPLEX_DATA_LAYOUT),intent(inout)                 :: nlU1,nlV1,nlW1
  type(COMPLEX_DATA_LAYOUT),intent(in)                    :: nlU2,nlV2,nlW2
  logical, intent(in)                                     :: mhd
  type(COMPLEX_DATA_LAYOUT), dimension(:), intent(inout)  :: nlB1
  type(COMPLEX_DATA_LAYOUT), intent(in)                   :: nlBx2, nlBy2, nlBz2
  type(COMPLEX_DATA_LAYOUT),intent(inout),dimension(:)    :: ScalNL1
  type(COMPLEX_DATA_LAYOUT),intent(in),dimension(:)       :: ScalNL2
  real(WP), intent(in)                                    :: coef1, coef2, step
  logical , intent(out)                                   :: res
  ! local
  real(WP)  :: coef   ! to compute exponential term for diffusion in Fourrier space.
  real(WP)  :: kk, ekk! wave number norm and exponential of (coeff*kk)
  real(WP)  :: kkz    ! = square of third componnent of wave numbers
  real(WP)  :: kkyz   ! = sum of square of third and second componnent of wave numbers ! kk = kkyz + (WN_x)**2
  integer   :: i,j,k  ! loop indices
  integer   :: sca    ! loop indice

  ! == Init ==
  res = .false.
  coef = -mu*step

  ! == For velocty ==
  do k = nlU1%zmin,nlU1%zmax
    kkz = VelWN%kz(k)*VelWN%kz(k)
    do j = nlU1%ymin,nlU1%ymax
      kkyz = kkz + VelWN%ky(j)*VelWN%ky(j)
      do i = nlU1%xmin,nlU1%xmax
        kk = kkyz + VelWN%kx(i)*VelWN%kx(i)
        ekk=exp(coef*kk)
        nlU1%values(i,j,k) = (coef1*nlU1%values(i,j,k)) + (ekk*coef2*nlU2%values(i,j,k))
        nlV1%values(i,j,k) = (coef1*nlV1%values(i,j,k)) + (ekk*coef2*nlV2%values(i,j,k))
        nlW1%values(i,j,k) = (coef1*nlW1%values(i,j,k)) + (ekk*coef2*nlW2%values(i,j,k))
      end do
    end do
  end do

  ! == For scalars ==
  do sca=1,nbscal
    do k = ScalArrayK(sca)%zmin,ScalArrayK(sca)%zmax
      kkz = ScalWN(sca)%kz(k)*ScalWN(sca)%kz(k)
      do j = ScalArrayK(sca)%ymin,ScalArrayK(sca)%ymax
        kkyz = kkz + ScalWN(sca)%ky(j)*ScalWN(sca)%ky(j)
        do i = ScalArrayK(sca)%xmin,ScalArrayK(sca)%xmax
          kk = ScalWN(sca)%kx(i)*ScalWN(sca)%kx(i) + kkyz
          ekk=exp(kk*coef/schmidt(sca))
          ScalNL1(sca)%values(i,j,k) = (coef1*ScalNL1(sca)%values(i,j,k)) &
            & + (ekk*coef2*ScalNL2(sca)%values(i,j,k))
        end do
      end do
    end do
  end do

  ! == For MHD if needed ==
  if (mhd) then
    do k = Bk(1)%zmin,Bk(1)%zmax
      kkz = BfieldWN%kz(k)*BfieldWN%kz(k)
      do j = Bk(1)%ymin,Bk(1)%ymax
        kkyz = kkz + BfieldWN%ky(j)*BfieldWN%ky(j)
        do i = Bk(1)%xmin,Bk(1)%xmax
          kk = BfieldWN%kx(i)*BfieldWN%kx(i) + kkyz
          ekk=exp(kk*coef/Prm)
          nlB1(1)%values(i,j,k) = (coef1*nlB1(1)%values(i,j,k)) + (ekk*coef2*nlBx2%values(i,j,k))
          nlB1(2)%values(i,j,k) = (coef1*nlB1(2)%values(i,j,k)) + (ekk*coef2*nlBy2%values(i,j,k))
          nlB1(3)%values(i,j,k) = (coef1*nlB1(3)%values(i,j,k)) + (ekk*coef2*nlBz2%values(i,j,k))
        end do
      end do
    end do
  end if ! mhd

  res = .true.

end subroutine compute_AB_NL

!------------------------------------------------------------------------------
!> Init save array for time integration of order bigger than 1 (k_old, final
!! non linear term for RK4
!! @autor Jean-Baptiste Lagaert
!!     @param [in]     imhd         = 1 if MHD is active, 0 else
!!     @return         res          = true if no error occurs
!------------------------------------------------------------------------------
function init_save(Uk,Vk,Wk,Bk,ScalSk,size_save, imhd) result(res)

  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN)              :: Uk, Vk, Wk
  TYPE(COMPLEX_DATA_LAYOUT), DIMENSION(:), INTENT(INOUT)  :: Bk
  TYPE(COMPLEX_DATA_LAYOUT), DIMENSION(:), INTENT(INOUT)  :: ScalSk
  INTEGER,                INTENT(IN)              :: size_save
  LOGICAL, INTENT(IN)                             :: imhd
  LOGICAL                                         :: res
  ! Local
  integer :: i,j

  res = .false.

  allocate(Uk_sav(size_save))
  allocate(Vk_sav(size_save))
  allocate(Wk_sav(size_save))
  do i = 1, size_save
    if (.not.(copyStructOnly(Uk,Uk_sav(i)).and. copyStructOnly(Vk,Vk_sav(i)) &
        & .and. (copyStructOnly(Wk,Wk_sav(i))))) then
      write(6,'(a,i0)') '[ERROR] init_save not enought memory for Uk_sav, Vk_sav and Wk_sav'
      return
    end if
  end do

  if (forcing_save) then
    if (.not.(copyStructOnly(Uk,Fk_sav(1)).and. copyStructOnly(Vk,Fk_sav(2)) &
        & .and. (copyStructOnly(Wk,Fk_sav(3))))) then
      write(6,'(a,i0)') '[ERROR] init_save not enought memory for Fk_sav'
      return
    end if
  end if

  allocate(Bk_sav(3,size_save))
  if (imhd) then
    do i = 1, size_save
      if (.not.(copyStructOnly(Bk(1),Bk_sav(1,i)).and. copyStructOnly(Bk(1),Bk_sav(2,i)) &
          & .and. (copyStructOnly(Bk(1),Bk_sav(3,i))))) then
        write(6,'(a,i0)') '[ERROR] init_save not enought memory for Bk_sav'
        return
      end if
    end do
  end if

  if (nbscal>0) then
    allocate(Scalk_sav(nbscal,size_save))
    do i = 1, size_save
      do j =1,nbscal
        if (.not.(copyStructOnly(ScalSk(j),Scalk_sav(j,i)))) then
          write(6,'(a,i0)') '[ERROR] init_save not enought memory for Scalk_sav'
          return
        end if
      end do
    end do
  else
    allocate(Scalk_sav(1,size_save))
  end if

  if(dt_save>0) then
    allocate(dt_prev(dt_save))
    dt_prev = -1
  end if

  res = .true.

end function init_save


!------------------------------------------------------------------------------
!> Delete save array for time integration of order bigger than 1 (k_old, final
!! non linear term for RK4
!! @autor Jean-Baptiste Lagaert
!!     @param [in]     nbcpus       = the total amount of processes (or cpus)
!!     @param [in]     spec_rank    = my MPI processe rank in the range [0:nbcpus-1]
!!     @param [in]     imhd         = 1 if MHD is active, 0 else
!!     @return         res          = true if no error occurs
!------------------------------------------------------------------------------
function delete_save(save_size, imhd) result(res)

  INTEGER,                INTENT(IN)              :: save_size
  LOGICAL, INTENT(IN)                             :: imhd
  LOGICAL                                         :: res
  ! Local
  integer :: i,j

  res = .false.

  do i = 1, save_size
    call deleteDataLayout(Uk_sav(i))
    call deleteDataLayout(Vk_sav(i))
    call deleteDataLayout(Wk_sav(i))
  end do
  deallocate(Uk_sav)
  deallocate(Vk_sav)
  deallocate(Wk_sav)

  if (forcing_save) then
    call deleteDataLayout(Fk_sav(1))
    call deleteDataLayout(Fk_sav(2))
    call deleteDataLayout(Fk_sav(3))
  end if

  if (imhd) then
    do i = 1, save_size
      call deleteDataLayout(Bk_sav(1,i))
      call deleteDataLayout(Bk_sav(2,i))
      call deleteDataLayout(Bk_sav(3,i))
    end do
    deallocate(Bk_sav)
  end if

  if (nbscal>0) then
    do i = 1, save_size
      do j =1,nbscal
        call deleteDataLayout(Scalk_sav(j,i))
      end do
    end do
    deallocate(Scalk_sav)
  end if

  if(dt_save>0) deallocate(dt_prev)

  res = .true.

end function delete_save

!------------------------------------------------------------------------------
!> Init work array for non linear terms and for RK scheme.
!! @autor Patrick Begou
!!     @param [in]     nbcpus       = the total amount of processes (or cpus)
!!     @param [in]     spec_rank    = my MPI processe rank in the range [0:nbcpus-1]
!!     @param [in]     imhd         = 1 if MHD is active, 0 else
!!     @return         res          = true if no error occurs
!! @details
!!     Observe that the datalayout for "kOld" terms (ie terms use for RK scheme)
!! is initialized in the same time that these terms are affected (the "=" do the
!! datalayout initialisation).
!------------------------------------------------------------------------------
function init_workArray(U,B,ScalArray,imhd, nbcpus, spec_rank) result(res)

    TYPE(REAL_DATA_LAYOUT), INTENT(IN)              :: U
    TYPE(REAL_DATA_LAYOUT), POINTER, INTENT(INOUT)  :: B(:)
    TYPE(REAL_DATA_LAYOUT), POINTER, INTENT(INOUT)  :: ScalArray(:)
    INTEGER,                INTENT(IN)              :: nbcpus,spec_rank
    LOGICAL, INTENT(IN)                             :: imhd
    LOGICAL                                         :: res

    INTEGER  :: i,ierr

    ! Init
    res = .false.

    ! Allocate work arrays nl_x, nl_y and nl_z for velocities
    IF (.NOT.initDataLayout("workarray_Spectral",nl_x,(U%nx/2)+1,U%ny,U%nz,U%Lx,U%Ly,U%Lz,nbcpus,spec_rank,AlongZ)) THEN
        WRITE(6,'(a,i0)') '[ERROR] solver_step not enought memory for nlx, nly and nlz on ',spec_rank
        RETURN
    ENDIF
    IF (.NOT. copyStructOnly(nl_x,nl_y) .OR. .NOT. copyStructOnly(nl_x,nl_z)) RETURN

    ! Allocate work arrays for all the scalars
    IF (nbscal.GT.0) THEN
        IF(ALLOCATED(nl_ScalArray)) then
            DO i = 1, size(nl_ScalArray)
                CALL deleteDataLayout(nl_ScalArray(i))
            END DO
            DEALLOCATE(nl_ScalArray)
        END IF
        IF(ALLOCATED(ScalArraykOld)) then
            DO i = 1, size(ScalArraykOld)
                CALL deleteDataLayout(ScalArraykOld(i))
            END DO
            DEALLOCATE(ScalArraykOld)
        END IF
        ALLOCATE(nl_ScalArray(nbscal),ScalArraykOld(nbscal),stat=ierr)
       ! Save for RK2 scheme
        IF (ierr.NE.0) THEN
            WRITE(6,'(a,i0)') '[ERROR] solver_step not enought memory for ScalArrayNL and ScalArraykOld array on ',spec_rank
            RETURN
        ENDIF
        ! Init each component for the non linear terms
        DO i=1,nbscal
            IF (.NOT.initDataLayout("workarray_Spectral",nl_ScalArray(i), (ScalArray(i)%nx/2)+1, &
                        & ScalArray(i)%ny, &
                        & ScalArray(i)%nz,U%Lx,U%Ly,U%Lz,nbcpus,spec_rank,AlongZ)) THEN
                    WRITE(6,'(2(a,i0))') '[ERROR] solver_step not enought memory for ScalArrayNL(',i,') on ',spec_rank
                RETURN
            ENDIF
        ENDDO
    ENDIF

    ! Allocate work array for B field, nlB(1:3)
    IF (imhd) THEN
      IF(ALLOCATED(nl_B)) then
          DO i = 1, size(nl_B)
              CALL deleteDataLayout(nl_B(i))
          END DO
          DEALLOCATE(nl_B)
      END IF

      ALLOCATE(BkOld(3),stat=ierr)
      ALLOCATE(nl_B(3),stat=ierr)

      IF (.NOT.initDataLayout("workBarray_Spectral",nl_B(1),(B(1)%nx/2)+1,B(1)%ny,B(1)%nz,&
          & B(1)%Lx,B(1)%Ly,B(1)%Lz,nbcpus,spec_rank,AlongZ)) THEN
          WRITE(6,'(a,i0)') '[ERROR] solver_step not enought memory for nlBx, nlBy and nlBz on ',spec_rank
          RETURN
      ENDIF
      IF (.NOT. copyStructOnly(nl_B(1),nl_B(2)) .OR. .NOT. copyStructOnly(nl_B(1),nl_B(3))) RETURN
!     END IF

      ALLOCATE(divTub(3),stat=ierr)

      IF (.NOT.initDataLayout("workBarray_Spectral",divTub(1),(B(1)%nx/2)+1,B(1)%ny,B(1)%nz,&
          & B(1)%Lx,B(1)%Ly,B(1)%Lz,nbcpus,spec_rank,AlongZ)) THEN
          WRITE(6,'(a,i0)') '[ERROR] solver_step not enought memory for divTubx, divTuby and divTubz on ',spec_rank
          RETURN
      ENDIF
      IF (.NOT. copyStructOnly(divTub(1),divTub(2)) .OR. .NOT. copyStructOnly(divTub(1),divTub(3))) RETURN
    END IF


    res = .true.

end function init_workArray


!------------------------------------------------------------------------------
!> Delete work array for non linear terms and RK scheme
!! @autor Patrick Begou
!!     @param [in]     imhd         = 1 if MHD is active, 0 else
!!     @return         res          = true if no error occurs
!------------------------------------------------------------------------------
function delete_workArray(imhd) result(res)

    LOGICAL, INTENT(IN)                             :: imhd
    LOGICAL                                         :: res

    INTEGER  :: i

    ! Init
    res = .false.

    !free module Working arrays
    CALL DeleteDataLayout(nl_x)
    CALL DeleteDataLayout(nl_y)
    CALL DeleteDataLayout(nl_z)

    DO i=1,nbscal
        CALL DeleteDataLayout(nl_ScalArray(i))
        CALL DeleteDataLayout(ScalArraykOld(i))
    ENDDO
    IF(nbscal.GT.0) DEALLOCATE(nl_ScalArray,ScalArraykOld)

    IF (imhd) THEN
       CALL DeleteDataLayout(nl_B(1))
       CALL DeleteDataLayout(nl_B(2))
       CALL DeleteDataLayout(nl_B(3))
       CALL DeleteDataLayout(BkOld(1))
       CALL DeleteDataLayout(BkOld(2))
       CALL DeleteDataLayout(BkOld(3))
       CALL DeleteDataLayout(divTub(1))
       CALL DeleteDataLayout(divTub(2))
       CALL DeleteDataLayout(divTub(3))
       DEALLOCATE(divTub)
       DEALLOCATE(nl_B)
       DEALLOCATE(BkOld)
    ENDIF

    res = .true.

end function delete_workArray

END MODULE solver


