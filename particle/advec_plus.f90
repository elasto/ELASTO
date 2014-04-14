!USEFORTEST advec
!> @addtogroup part
!! @{
!------------------------------------------------------------------------------
!
! MODULE: advec
!
!
! DESCRIPTION:
!> The module advec provides all public interfaces to solve an advection equation
!! with a particle method.
!
!> @details
!!     This module contains the generic procedure to initialize and parametrise the
!! advection solver based on particles method. It also contains the subroutine
!! "advec_step" wich solves the equation for a given time step. It is the only one
!! module which is supposed to be included by a code using this library of
!! particle methods.
!
!> @author
!! Jean-Baptiste Lagaert, LEGI
!
!------------------------------------------------------------------------------
module advec_Plus

    use precision_tools
    use datalayout
    use advec_abstract_proc
    implicit none

    ! ===== Private variables =====
    !> numerical method use to advect the scalar
    character(len=str_short), private   :: type_part_solv
    !> dimensionnal splitting (eg classical, Strang or particle)
    character(len=str_short), private   :: dim_splitting
    !> Group size along current direction
    integer, private, dimension(2)  :: gsX, gsY, gsZ


    ! ===== Public procedures =====
    ! Scheme used to advec the scalar (order 2 or 4 ?)
!    public                              :: type_part_solver

    ! Advection methods
!    public                              :: advec_init           ! initialize the scalar solver
    public                              :: advec_step_Plus      ! advec the scalar field during a time step.
!    procedure(advec_step_Torder2), pointer, public    :: advec_step => null()
!    public                              :: advec_step_Torder1   ! advec the scalar field during a time step.
!    public                              :: advec_step_Torder2   ! advec the scalar field during a time step.
!
    ! Remeshing formula
    procedure(AC_remesh), pointer, private :: advec_remesh_plus => null()

contains

! ===== Public methods =====

!> Return the name of the particle method used for the advection
!!    @return type_part_solver      = numerical method used for advection
function type_part_solver()
    character(len=str_short)    :: type_part_solver

    type_part_solver = type_part_solv
end function


!> Initialise the particle advection methods
!!    @param[in]    order       = to choose the remeshing method (and thus the order)
!!    @param[out]   stab_coeff  = stability coefficient (condition stability is
!!                                  dt< stab_coeff/norm_inf(V))
!!    @param[in]    dim_split   = dimensionnal splitting (eg classical,
!!                                    Strang splitting or particle splitting)
!!    @param[in]    verbosity   = to display info about chosen remeshing formula (optional)
subroutine advec_init_Plus(order, stab_coeff, verbosity, dim_split)

    use advec_variables ! contains info about solver parameters and others.
    use cart_topology   ! Description of mesh and of mpi topology
    use advecX          ! solver for advection along X
    use advecY          ! solver for advection along Y
    use advecZ          ! solver for advection along Z
    use advec_common    ! some procedures common to advection along all directions

    ! Input/Output
    character(len=*), optional, intent(in)  ::  order, dim_split
    logical, optional, intent(in)           ::  verbosity
    real(WP), optional, intent(out)         ::  stab_coeff

    ! Use default solver if it is not chosen by the user.
    if(present(order)) then
        type_part_solv = order
    else
        type_part_solv = 'p_O2'
    end if

    ! Initialize the solver
    if (present(verbosity)) then
        call AC_solver_init(type_part_solv, verbosity)
    else
        call AC_solver_init(type_part_solv)
    end if

    ! Compute stability coefficient - as each dimension is solved in
    ! dt/2, stab_coef is 2 times bigger
    if (present(stab_coeff)) stab_coeff = 1.0/(dble(bl_size))

    ! Call the right remeshing formula
    select case(type_part_solv)
        case('p_O2')
            advec_remesh_plus => AC_remesh_lambda_group ! or Xremesh_O2
        case('p_O4')
            advec_remesh_plus => AC_remesh_lambda_group ! or Xremesh_O4
        case('p_L2')
            advec_remesh_plus => AC_remesh_limit_lambda_group ! lambda 2 limited
        case('p_M4')
            advec_remesh_plus => AC_remesh_Mprime_group ! Xremesh_Mprime6
        case('p_M6')
            advec_remesh_plus => AC_remesh_Mprime_group ! Xremesh_Mprime6
        case('p_L4')
            advec_remesh_plus => AC_remesh_Mprime_group ! Lambda 4,4
        case('p_L6')
            advec_remesh_plus => AC_remesh_Mprime_group ! Lambda 6,6
        case('p_M8')
            advec_remesh_plus => AC_remesh_Mprime_group ! Xremesh_Mprime6
        case('d_M4')
            advec_remesh_plus => AC_remesh_Mprime_group ! Xremesh_Mprime6
        case default
            advec_remesh_plus => AC_remesh_lambda_group ! or Xremesh_O2
    end select

    call AC_setup_init()
    call advecX_remesh_init()

    ! Save group size
    gsX =group_size(1,:)
    gsY =group_size(2,:)
    gsZ =group_size(3,:)

end subroutine advec_init_Plus


!> Solve advection equation - order 2 in time (order 2 dimensional splitting)
!!    @param[in]        dt          = time step
!!    @param[in]        Vx          = velocity along x (could be discretised on a bigger mesh then the scalar)
!!    @param[in]        Vy          = velocity along y
!!    @param[in]        Vz          = velocity along z
!!    @param[in,out]    ScalArray   = scalar fields to advect
subroutine advec_step_Plus(dt, Vx, Vy, Vz, ScalArray)

    use advecX          ! Method to advec along X
    use advecY          ! Method to advec along Y
    use advecZ          ! Method to advec along Z
    use advec, only : advec_setup_alongX, advec_setup_alongY, advec_setup_alongZ
    use advec_remeshing_Mprime, only : current_dir

    ! Input/Output
    real(WP), intent(in)                        :: dt
    real(WP), dimension(:,:,:), intent(in)      :: Vx, Vy, Vz
    type(REAL_DATA_LAYOUT), dimension(:), intent(inout) :: ScalArray

    current_dir = 1
    call advec_setup_alongX()
    call advec_plus_X_basic_no_com(dt/2.0, gsX, Vx, ScalArray)
    current_dir = 2
    call advec_setup_alongY()
    call advec_plus_1D_basic(dt/2.0, gsY, Vy, ScalArray)
    current_dir = 3
    call advec_setup_alongZ()
    call advec_plus_1D_basic(dt/2.0, gsZ, Vz, ScalArray)
    call advec_plus_1D_basic(dt/2.0, gsZ, Vz, ScalArray)
    current_dir = 2
    call advec_setup_alongY()
    call advec_plus_1D_basic(dt/2.0, gsY, Vy, ScalArray)
    current_dir = 1
    call advec_setup_alongX()
    call advec_plus_X_basic_no_com(dt/2.0, gsX, Vx, ScalArray)

end subroutine advec_step_Plus


!> Solve advection equation - order 2 in time (order 2 dimensional splitting)
!! Velocity is interpolated at particle position directly from its (coarser) grid with High Order formula.
!!    @param[in]        dt          = time step
!!    @param[in]        Vx          = velocity along x (could be discretised on a bigger mesh then the scalar)
!!    @param[in]        Vy          = velocity along y
!!    @param[in]        Vz          = velocity along z
!!    @param[in,out]    ScalArray   = scalar fields to advect
subroutine advec_step_Plus_Inter_HO(dt, Vx, Vy, Vz, ScalArray)

    use Interpolation_velo

    ! Input/Output
    real(WP), intent(in)                        :: dt
    real(WP), dimension(:,:,:), intent(in)      :: Vx, Vy, Vz
    type(REAL_DATA_LAYOUT), dimension(:), intent(inout) :: ScalArray
    ! Local
    real(WP), dimension(:,:,:), allocatable   :: Vx_c, Vy_c, Vz_c
    real(WP), dimension(:,:,:), allocatable   :: Vx_f, Vy_f, Vz_f
    integer                                   :: ierr                ! Error code.

    allocate(Vx_c(mesh_V%N_proc(1),mesh_sc%N_proc(2),mesh_sc%N_proc(3)),stat=ierr)
    if (ierr/=0) write(6,'(a,i0,a)') '[ERROR] on cart_rank ', cart_rank, ' - not enough memory for Vx_c'
    allocate(Vx_f(mesh_sc%N_proc(1),mesh_sc%N_proc(2),mesh_sc%N_proc(3)),stat=ierr)
    if (ierr/=0) write(6,'(a,i0,a)') '[ERROR] on cart_rank ', cart_rank, ' - not enough memory for Vx_f'
    allocate(Vy_c(mesh_V%N_proc(2),mesh_sc%N_proc(1),mesh_sc%N_proc(3)),stat=ierr)
    if (ierr/=0) write(6,'(a,i0,a)') '[ERROR] on cart_rank ', cart_rank, ' - not enough memory for Vy_c'
    allocate(Vy_f(mesh_sc%N_proc(2),mesh_sc%N_proc(1),mesh_sc%N_proc(3)),stat=ierr)
    if (ierr/=0) write(6,'(a,i0,a)') '[ERROR] on cart_rank ', cart_rank, ' - not enough memory for Vy_f'
    allocate(Vz_c(mesh_V%N_proc(3),mesh_sc%N_proc(1),mesh_sc%N_proc(2)),stat=ierr)
    if (ierr/=0) write(6,'(a,i0,a)') '[ERROR] on cart_rank ', cart_rank, ' - not enough memory for Vz_c'
    allocate(Vz_f(mesh_sc%N_proc(3),mesh_sc%N_proc(1),mesh_sc%N_proc(2)),stat=ierr)
    if (ierr/=0) write(6,'(a,i0,a)') '[ERROR] on cart_rank ', cart_rank, ' - not enough memory for Vz_f'

    call Interpol_2D_3D_vect(mesh_sc%dx, mesh_V%dx, Vx, Vy, Vz, Vx_c, Vx_f, Vy_c, Vy_f, Vz_c, Vz_f)

    call advec_step_Plus_Vcoarse(dt, Vx_c, Vy_c, Vz_c, Vx_f, Vy_f, Vz_f, ScalArray)

    deallocate(Vx_f)
    deallocate(Vy_f)
    deallocate(Vz_f)

    deallocate(Vx_c)
    deallocate(Vy_c)
    deallocate(Vz_c)

end subroutine advec_step_Plus_Inter_HO


!> Solve advection equation - order 2 in time (order 2 dimensional splitting)
!! Velocity is interpolated from velocity grid to scalar resolution with high
!! order formula and then the classical solver is used.
!!    @param[in]        dt          = time step
!!    @param[in]        Vx          = velocity along x (could be discretised on a bigger mesh then the scalar)
!!    @param[in]        Vy          = velocity along y
!!    @param[in]        Vz          = velocity along z
!!    @param[in,out]    ScalArray   = scalar fields to advect
subroutine advec_step_Plus_Inter_basic(dt, Vx, Vy, Vz, ScalArray)

    use Interpolation_velo

    ! Input/Output
    real(WP), intent(in)                        :: dt
    real(WP), dimension(:,:,:), intent(in)      :: Vx, Vy, Vz
    type(REAL_DATA_LAYOUT), dimension(:), intent(inout) :: ScalArray
    ! Local
    real(WP), dimension(:,:,:), allocatable   :: Vx_f, Vy_f, Vz_f
    integer                                   :: ierr                ! Error code.

    allocate(Vx_f(mesh_sc%N_proc(1),mesh_sc%N_proc(2),mesh_sc%N_proc(3)),stat=ierr)
    if (ierr/=0) write(6,'(a,i0,a)') '[ERROR] on cart_rank ', cart_rank, ' - not enough memory for Vx_f'
    allocate(Vy_f(mesh_sc%N_proc(1),mesh_sc%N_proc(2),mesh_sc%N_proc(3)),stat=ierr)
    if (ierr/=0) write(6,'(a,i0,a)') '[ERROR] on cart_rank ', cart_rank, ' - not enough memory for Vy_f'
    allocate(Vz_f(mesh_sc%N_proc(1),mesh_sc%N_proc(2),mesh_sc%N_proc(3)),stat=ierr)
    if (ierr/=0) write(6,'(a,i0,a)') '[ERROR] on cart_rank ', cart_rank, ' - not enough memory for Vz_f'

    call Interpol_3D(Vx, mesh_V%dx, Vx_f, mesh_sc%dx)
    call Interpol_3D(Vy, mesh_V%dx, Vy_f, mesh_sc%dx)
    call Interpol_3D(Vz, mesh_V%dx, Vz_f, mesh_sc%dx)
    if (cart_rank==0) write(*,'(a)') "        [INFO PARTICLES] Interpolation done"

    call advec_step_Plus(dt, Vx_f, Vy_f, Vz_f, ScalArray)

    deallocate(Vx_f)
    deallocate(Vy_f)
    deallocate(Vz_f)

end subroutine advec_step_Plus_Inter_Basic

!> Solve advection equation - order 2 in time (order 2 dimensional splitting) -
!! variant with velocity given on coarser grid.
!!    @param[in]        dt          = time step
!!    @param[in]        Vx_c        = velocity along x on coarser grid than the scalar one
!!    @param[in]        Vy_c        = velocity along y on coarser grid than the scalar one
!!    @param[in]        Vz_c        = velocity along z on coarser grid than the scalar one
!!    @param[in]        Vx_f        = velocity along x on the same grid than the scalar field
!!    @param[in]        Vy_f        = velocity along y on the same grid than the scalar field
!!    @param[in]        Vz_f        = velocity along z on the same grid than the scalar field
!!    @param[in,out]    ScalArray   = scalar fields to advect
subroutine advec_step_Plus_Vcoarse(dt, Vx_c, Vy_c, Vz_c, Vx_f, Vy_f, Vz_f, ScalArray)

    use advecX          ! Method to advec along X
    use advecY          ! Method to advec along Y
    use advecZ          ! Method to advec along Z
    use advec, only : advec_setup_alongX, advec_setup_alongY, advec_setup_alongZ
    use advec_remeshing_Mprime, only : current_dir

    ! Input/Output
    real(WP), intent(in)                        :: dt
    real(WP), dimension(:,:,:), intent(in)      :: Vx_c, Vy_c, Vz_c
    real(WP), dimension(:,:,:), intent(in)      :: Vx_f, Vy_f, Vz_f
    type(REAL_DATA_LAYOUT), dimension(:), intent(inout) :: ScalArray

    current_dir = 1
    call advec_setup_alongX()
    call advec_plus_X_Vcoarse_no_com(dt/2.0, gsX, Vx_c, Vx_f, ScalArray)
    current_dir = 2
    call advec_setup_alongY()
    call advec_plus_1D_Vcoarse(dt/2.0, gsY, Vy_c, Vy_f, ScalArray)
    current_dir = 3
    call advec_setup_alongZ()
    call advec_plus_1D_Vcoarse(dt/2.0, gsZ, Vz_c, Vz_f, ScalArray)
    call advec_plus_1D_Vcoarse(dt/2.0, gsZ, Vz_c, Vz_f, ScalArray)
    current_dir = 2
    call advec_setup_alongY()
    call advec_plus_1D_Vcoarse(dt/2.0, gsY, Vy_c, Vy_f, ScalArray)
    current_dir = 1
    call advec_setup_alongX()
    call advec_plus_X_Vcoarse_no_com(dt/2.0, gsX, Vx_c, Vx_f, ScalArray)

end subroutine advec_step_Plus_Vcoarse


!> ===== Private procedure =====
!> Scalar advection along one direction - variant for cases with no communication
!!    @param[in]        dt          = time step
!!    @param[in]        gs          = size of the work item along transverse direction
!!    @param[in]        V_comp      = velocity along X (could be discretised on a bigger mesh then the scalar)
!!    @param[in,out]    ScalArray   = array of scalar fields to advect
!> @details
!!   Work only for direction = X. Basic (and very simple) remeshing has just to
!! be add for other direction.
subroutine advec_plus_X_basic_no_com(dt, gs, V_comp, ScalArray)

    use advec, only : advec_init_velo, advec_remesh, line_dir, gp_dir1, gp_dir2
    use advecX          ! Procedure specific to advection along X
    use advec_common    ! Some procedures common to advection along all directions
    use advec_variables ! contains info about solver parameters and others.
    use cart_topology   ! Description of mesh and of mpi topology
    use advec_remeshing_Mprime, only : sc_remesh_ind

    ! Input/Output
    real(WP), intent(in)                          :: dt
    integer, dimension(2), intent(in)             :: gs
    real(WP), dimension(:,:,:), intent(in)        :: V_comp
    type(REAL_DATA_LAYOUT), dimension(:), intent(inout) :: ScalArray
    ! Other local variables
    integer                                                   :: j,k          ! indice of the currend mesh point
    integer                                                   :: sca          ! indice of the currend scalar indice
    integer, dimension(2)                                     :: ind_group    ! indice of the currend group of line (=(i,k) by default)
    real(WP),dimension(mesh_sc%N_proc(line_dir),gs(1),gs(2))  :: p_pos_adim   ! adimensionned particles position
    real(WP),dimension(mesh_sc%N_proc(line_dir),gs(1),gs(2))  :: p_V          ! particles velocity

    ind_group = 0

    do k = 1, mesh_sc%N_proc(gp_dir2), gs(2)
        ind_group(2) = ind_group(2) + 1
        ind_group(1) = 0
        do j = 1, mesh_sc%N_proc(gp_dir1), gs(1)
            ind_group(1) = ind_group(1) + 1

            ! ===== Init particles =====
            ! p_pos is used to store velocity at grid point
            call advec_init_velo(V_comp, j, k, gs, p_pos_adim)
            ! p_V = middle point position = position at middle point for RK2 scheme
            call AC_get_p_pos_adim(p_V, p_pos_adim, 0.5_WP*dt, mesh_sc%dx(line_dir), mesh_sc%N_proc(line_dir))

            ! ===== Advection =====
            ! -- Compute velocity (with a RK2 scheme): p_V = velocity at middle point position --
            ! Note that p_pos is used as velocity component storage
            call AC_interpol_lin_no_com(line_dir, gs, p_pos_adim, p_V)
            ! p_v = velocity at middle point position
            ! -- Push particles --
            call AC_get_p_pos_adim(p_pos_adim, p_V, dt, mesh_sc%dx(line_dir), mesh_sc%N_proc(line_dir))
            ! Now p_pos = particle position and p_V = particle velocity

            ! ===== Remeshing =====
            do sca = 1, size(ScalArray)
                sc_remesh_ind = sca
                call advecX_remesh_no_com(ind_group, gs, p_pos_adim, p_V, j, k, ScalArray(sca)%values, dt)
            end do

        end do
    end do

end subroutine advec_plus_X_basic_no_com


!> Scalar advection along one direction (this procedure call the right solver, depending on the simulation setup)
!!    @param[in]        dt          = time step
!!    @param[in]        gs          = size of the work item along transverse direction
!!    @param[in]        V_comp      = velocity component
!!    @param[in,out]    ScalArray   = array of scalar fields to advect
subroutine advec_plus_1D_basic(dt, gs, V_comp, ScalArray)

    use advec, only : advec_init_velo, advec_remesh, line_dir, gp_dir1, gp_dir2
    use advecX, only : advecX_init_group    ! procdure devoted to advection along Z
    use advecY, only : advecY_init_group    ! procdure devoted to advection along Z
    use advecZ, only : advecZ_init_group    ! procdure devoted to advection along Z
    use advec_variables ! contains info about solver parameters and others.
    use cart_topology   ! Description of mesh and of mpi topology
    use advec_common    ! some procedures common to advection along all line_dirs
    use advec_remeshing_Mprime, only : sc_remesh_ind

    ! Input/Output
    real(WP), intent(in)                                :: dt
    integer, dimension(2), intent(in)                   :: gs
    real(WP), dimension(:,:,:), intent(in)              :: V_comp
    type(REAL_DATA_LAYOUT), dimension(:), intent(inout) :: ScalArray
    ! Other local variables
    integer                                                   :: i,j          ! indice of the currend mesh point
    integer                                                   :: sca          ! indice of the currend scalar indice
    integer, dimension(2)                                     :: ind_group    ! indice of the currend group of line (=(i,k) by default)
    real(WP), dimension(mesh_sc%N_proc(line_dir),gs(1),gs(2)) :: p_pos_adim ! adimensionned particles position
    real(WP), dimension(mesh_sc%N_proc(line_dir),gs(1),gs(2)) :: p_V        ! particles velocity

    ind_group = 0

    do j = 1, mesh_sc%N_proc(gp_dir2), gs(2)
        ind_group(2) = ind_group(2) + 1
        ind_group(1) = 0
        do i = 1, mesh_sc%N_proc(gp_dir1), gs(1)
            ind_group(1) = ind_group(1) + 1

            ! ===== Init particles =====
            call advec_init_velo(V_comp, i, j, gs, p_pos_adim)
            ! p_pos is used to store velocity at grid point
            call AC_get_p_pos_adim(p_V, p_pos_adim, 0.5_WP*dt, mesh_sc%dx(line_dir), mesh_sc%N_proc(line_dir))
            ! p_V = middle point position = position at middle point for RK2 scheme

            ! ===== Advection =====
            ! -- Compute velocity (with a RK2 scheme) --
            ! Note that p_pos is used as velocity component storage
            call AC_interpol_lin(line_dir, gs, ind_group, p_pos_adim, p_V)
            ! p_v = velocity at middle point position
            ! -- Push particles --
            call AC_get_p_pos_adim(p_pos_adim, p_V, dt, mesh_sc%dx(line_dir), mesh_sc%N_proc(line_dir))
            ! Now p_pos = particle position and p_V = particle velocity

            ! ===== Remeshing =====
            do sca = 1, size(ScalArray)
                sc_remesh_ind = sca
                call advec_remesh_plus(line_dir, ind_group, gs, &
                    & p_pos_adim, p_V, i, j, ScalArray(sca)%values, dt)
            end do

        end do
    end do

end subroutine advec_plus_1D_basic


!> Scalar advection along one direction - variant for cases with no communication and velocity computed on a coarser grid.
!!    @param[in]        dt          = time step
!!    @param[in]        gs          = size of the work item along transverse direction
!!    @param[in]        V_comp      = velocity along X (could be discretised on a bigger mesh then the scalar)
!!    @param[in,out]    ScalArray   = array of scalar fields to advect
!> @details
!!   Work only for direction = X. Basic (and very simple) remeshing has just to
!! be add for other direction.
subroutine advec_plus_X_Vcoarse_no_com(dt, gs, V_coarse, V_fine, ScalArray)

    use advec, only : advec_init_velo, advec_remesh, line_dir, gp_dir1, gp_dir2
    use advecX          ! Procedure specific to advection along X
    use advec_common    ! Some procedures common to advection along all directions
    use advec_variables ! contains info about solver parameters and others.
    use cart_topology   ! Description of mesh and of mpi topology
    use advec_remeshing_Mprime, only : sc_remesh_ind

    ! Input/Output
    real(WP), intent(in)                          :: dt
    integer, dimension(2), intent(in)             :: gs
    real(WP), dimension(:,:,:), intent(in)        :: V_coarse
    real(WP), dimension(:,:,:), intent(in)        :: V_fine
    type(REAL_DATA_LAYOUT), dimension(:), intent(inout) :: ScalArray
    ! Other local variables
    integer                                                   :: j,k          ! indice of the currend mesh point
    integer                                                   :: sca          ! indice of the currend scalar indice
    integer, dimension(2)                                     :: ind_group    ! indice of the currend group of line (=(i,k) by default)
    real(WP),dimension(mesh_sc%N_proc(line_dir),gs(1),gs(2))  :: p_pos_adim   ! adimensionned particles position
    real(WP),dimension(mesh_sc%N_proc(line_dir)+1,gs(1),gs(2))  :: p_V          ! particles velocity

    ind_group = 0

    do k = 1, mesh_sc%N_proc(gp_dir2), gs(2)
        ind_group(2) = ind_group(2) + 1
        ind_group(1) = 0
        do j = 1, mesh_sc%N_proc(gp_dir1), gs(1)
            ind_group(1) = ind_group(1) + 1

            ! ===== Init particles =====
            call AC_get_p_pos_adim(p_V, V_fine, 0.5_WP*dt, &
                  & mesh_sc%dx(line_dir), mesh_V%dx(line_dir), mesh_sc%N_proc(line_dir),j, k)
            ! p_V = middle point position = position at middle point for RK2 scheme

            ! ===== Advection =====
            ! -- Compute velocity (with a RK2 scheme): p_V = velocity at middle point position --
            ! Note that p_pos is used as velocity component storage
            call AC_interpol_plus_no_com(line_dir, gs, j, k, V_coarse, p_V)
            ! p_v = velocity at middle point position
            ! -- Push particles --
            call AC_get_p_pos_adim(p_pos_adim, p_V, dt, mesh_sc%dx(line_dir), mesh_sc%N_proc(line_dir))
            ! Now p_pos = particle position and p_V = particle velocity

            ! ===== Remeshing =====
            do sca = 1, size(ScalArray)
                sc_remesh_ind = sca
                call advecX_remesh_no_com(ind_group, gs, p_pos_adim, p_V, j, k, ScalArray(sca)%values, dt)
            end do

        end do
    end do

end subroutine advec_plus_X_Vcoarse_no_com


!> Scalar advection along one direction (this procedure call the right solver, depending on the simulation setup)
!!    @param[in]        dt          = time step
!!    @param[in]        gs          = size of the work item along transverse direction
!!    @param[in]        V_comp      = velocity component
!!    @param[in,out]    ScalArray   = array of scalar fields to advect
subroutine advec_plus_1D_Vcoarse(dt, gs, V_coarse, V_fine, ScalArray)

    use advec, only : advec_init_velo, advec_remesh, line_dir, gp_dir1, gp_dir2
    use advecX, only : advecX_init_group    ! procdure devoted to advection along Z
    use advecY, only : advecY_init_group    ! procdure devoted to advection along Z
    use advecZ, only : advecZ_init_group    ! procdure devoted to advection along Z
    use advec_variables ! contains info about solver parameters and others.
    use cart_topology   ! Description of mesh and of mpi topology
    use advec_common    ! some procedures common to advection along all line_dirs
    use advec_remeshing_Mprime, only : sc_remesh_ind

    ! Input/Output
    real(WP), intent(in)                                :: dt
    integer, dimension(2), intent(in)                   :: gs
    real(WP), dimension(:,:,:), intent(in)              :: V_coarse
    real(WP), dimension(:,:,:), intent(in)              :: V_fine
    type(REAL_DATA_LAYOUT), dimension(:), intent(inout) :: ScalArray
    ! Other local variables
    integer                                                   :: i,j        ! indice of the currend mesh point
    integer                                                   :: sca        ! indice of the currend scalar indice
    integer, dimension(2)                                     :: ind_group  ! indice of the currend group of line (=(i,k) by default)
    real(WP), dimension(mesh_sc%N_proc(line_dir),gs(1),gs(2)) :: p_pos_adim ! adimensionned particles position
    real(WP), dimension(mesh_sc%N_proc(line_dir)+1,gs(1),gs(2)) :: p_V        ! particles velocity

    ind_group = 0

    do j = 1, mesh_sc%N_proc(gp_dir2), gs(2)
        ind_group(2) = ind_group(2) + 1
        ind_group(1) = 0
        do i = 1, mesh_sc%N_proc(gp_dir1), gs(1)
            ind_group(1) = ind_group(1) + 1

            ! ===== Init particles =====
            call AC_get_p_pos_adim(p_V, V_fine, 0.5_WP*dt, &
                  & mesh_sc%dx(line_dir), mesh_V%dx(line_dir), mesh_sc%N_proc(line_dir),i, j)
            ! p_V = middle point position = position at middle point for RK2 scheme

            ! ===== Advection =====
            ! -- Compute velocity (with a RK2 scheme) --
            call AC_interpol_plus(line_dir, gs, ind_group, i, j, V_coarse, p_V)
            ! p_v = velocity at middle point position
            ! -- Push particles --
            call AC_get_p_pos_adim(p_pos_adim, p_V, dt, mesh_sc%dx(line_dir), mesh_sc%N_proc(line_dir))
            ! Now p_pos = particle position and p_V = particle velocity

            ! ===== Remeshing =====
            do sca = 1, size(ScalArray)
                sc_remesh_ind = sca
                call advec_remesh_plus(line_dir, ind_group, gs, &
                    & p_pos_adim, p_V, i, j, ScalArray(sca)%values, dt)
            end do

        end do
    end do

end subroutine advec_plus_1D_Vcoarse

!> ===== Private procedure =====
end module advec_Plus
!> @}
