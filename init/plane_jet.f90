!> @addtogroup init
!! @{ 
!---------------------------------------------------------------------------------
!> Plane jet initialisation for velocity and scalars.
!! @autor = Guillaume Balarac and Jean-Baptiste Lagaert, LEGI
!
!
!
!---------------------------------------------------------------------------------
module planejet

  use precision_tools

  implicit none

  private :: computeJet
  public  :: initPlaneJet

contains

!---------------------------------------------------------------------------------
!> Initialisation for a plane-jet.
!! @autor = Jean-Baptiste Lagaert, Guillaume Balarac, LEGI
!
!> @param[in,out]   U               = longitudinal velocity storage to initialize
!> @param[in,out]   V               = vertical velocity  storage to initialize
!> @param[in,out]   W               = spanwise velocity  storage to initialize
!> @param[in]       VelWN           = the wavenumbers for U,V and W
!> @param[in,out]   B               = magnetic field (3D vector)
!> @param[in,out]   BfieldWN        = wave numbers associated to the magnetic field
!> @param[in,out]   ScalArray       = scalars solved with full pseudo-psectral method (to initialized)
!> @param[in]       ScalWN          = wavenumbers for the scalars
!> @param[in,out]   Scal_partArray  = scalars solved with full pseudo-psectral method (to initialized)
!> @param[in]       Scal_partWN     = wavenumbers for the scalars
!> @param[in]       prm             = Prandt number
!> @param[in]       nbcpus          = the total amount of processes (or cpus)
!> @param[in]       spec_rank       = my MPI processe rank in the range [0:nbcpus-1]
!> @return          ok              = logical set to .TRUE. if all initialization success.
!---------------------------------------------------------------------------------
function initPlaneJet(U,V,W,B,VelWN,ScalArray, ScalWN, Scal_partArray, &
        &  Scal_partWN, BfieldWN,nbcpus,spec_rank) result(ok)

    use stat_tools 
    use isotropicturbulence

    !I/O data
    type(REAL_DATA_LAYOUT), intent(inout)           :: U
    type(REAL_DATA_LAYOUT), intent(inout)           :: V
    type(REAL_DATA_LAYOUT), intent(inout)           :: W
    type(REAL_DATA_LAYOUT), pointer, intent(inout)  :: B(:)
    type(REAL_DATA_LAYOUT), pointer, intent(inout)  :: ScalArray(:)
    type(REAL_DATA_LAYOUT), pointer, intent(inout)  :: Scal_partArray(:)
    type(WaveNumbers), intent(in)                   :: VelWN
    type(WaveNumbers), dimension(:), intent(in)     :: ScalWN
    type(WaveNumbers), dimension(:), intent(in)     :: Scal_partWN
    type(WaveNumbers), intent(in)                   :: BfieldWN
    integer                                         :: nbcpus, spec_rank
    logical                                         :: ok

    !Local data
    real(WP)        :: zz, Lz, dz                       ! to manage jet position
    real(WP)        :: stiffness, mask, vel, u_prime, ke,b_prime! to deal with noise and tanh shape
    real(WP)        :: U_jet, U_coflow, Htheta, epsil, nozzle ! Jet parameters
    integer         :: i,k    ! Loop indices
    type(REAL_DATA_LAYOUT)  :: BufField                 ! to compute normalisation
    real(wp) :: coef

    ok = .false.

    if     (.not. parser_is_defined('H/theta')) then
       write(6,'(a)')'[warning] H/theta is not set: abort !!'
       return
    elseif (.not. parser_is_defined('Jet velocity')) then
       write(6,'(a)')'[warning] Jet velocity is not set: abort !!'
       return
    elseif (.not. parser_is_defined('Co-flow velocity'))then
       write(6,'(a)')'[warning] Co-flow velocity is not set: abort !!'
       return
    elseif (.not. parser_is_defined('noise level')) then
       write(6,'(a)')'[warning] noise level is not set: abort !!'
       return
    endif


    !init parameters from input file
    call parser_read('H/theta',Htheta)
    call parser_read('Jet velocity',U_jet)
    call parser_read('Co-flow velocity',U_coflow)
    call parser_read('noise level',epsil)
    call parser_read('noise wave number ke',ke)

    stiffness = Htheta/4_WP
    nozzle = 1      ! all dimension is adimensionned by the jet size.

    ! ===== Noise term for velocity =====
    ! Put noise - use thi init
    Lz = U%Lx
    ok = hit_vector_init(U, V, W, VelWN, Lz, nbcpus, spec_rank,ke,1._WP)
    if (.not. ok) return

    ! ===== Normalisation =====
    ! -- Compute energy --
    if (.not. copyStructOnly(U,BufField) ) then
      write(6,'(a,i0,a)')'[ERROR] init Plane Jet, copy struct : failed!'
    endif
    BufField%values = U%values**2.0 + V%values**2.0 + W%values**2.0
    u_prime =computeFieldAvg( BufField , spec_rank , nbcpus )
    u_prime = 0.5_WP*u_prime
!     call deleteDataLayout(BufField)
    u_prime = (2._WP/3._WP)*sqrt(u_prime)
    ! -- Normalize --
    U%values = (epsil/u_prime)*U%values
    V%values = (epsil/u_prime)*V%values
    W%values = (epsil/u_prime)*W%values

    ! ===== Jet-plane part =====
    dz = U%Lz / U%nz
    Lz  = U%Lz
    ! -- Velocity --
    do k = U%zmin,U%zmax
        zz = real(k)*dz - Lz/2.0_WP
        call computeJet(U_jet,U_coflow,nozzle,stiffness,zz, mask, vel)
        U%values(:,:,k) =       mask*U%values(:,:,k)
        V%values(:,:,k) = vel + mask*V%values(:,:,k)
        W%values(:,:,k) =       mask*W%values(:,:,k)
    end do
    ! -- Scalar solved with pseudo-spectral solver --
    if (ok .and. ASSOCIATED(ScalArray)) then
        do i=1,size(ScalArray)
            dz = ScalArray(i)%Lz / ScalArray(i)%nz
            Lz  = ScalArray(i)%Lz
            do k = ScalArray(i)%zmin,ScalArray(i)%zmax
                zz = real(k)*dz - Lz/2.0_WP
                call computeJet(1._WP,0._WP,nozzle,stiffness,zz, mask, vel)
                ScalArray(i)%values(:,:,k) = vel
            end do
        end do
    endif
    ! -- Scalar solved with mixed particle/spectral solver --
    if (ok .and. ASSOCIATED(Scal_partArray)) then
        do i=1,size(Scal_partArray)
            dz = Scal_partArray(i)%Lz / Scal_partArray(i)%nz
            Lz  = Scal_partArray(i)%Lz
            do k = Scal_partArray(i)%zmin,Scal_partArray(i)%zmax
                zz = real(k)*dz - Lz/2.0_WP
                call computeJet(1._WP,0._WP,nozzle,stiffness,zz, mask, vel)
                Scal_partArray(i)%values(:,:,k) = vel
            end do
        end do
    endif

    if (ok .and. ASSOCIATED(B)) then
       CALl parser_read('Coef Ek/Em', coef)
       if (spec_rank.EQ.0) print*, 'coef',coef
        do i=1,size(B)
            dz = B(i)%Lz / B(i)%nz
            Lz  = B(i)%Lz
            do k = B(i)%zmin,B(i)%zmax
                zz = real(k)*dz - Lz/2.0_WP
                call computeJet(1._WP,0._WP,nozzle,stiffness,zz, mask, vel)
                B(i)%values(:,:,k) = vel
            end do
        end do
    if (spec_rank.EQ.0) print*, 'tata'

    BufField%values = U%values**2.0 + V%values**2.0 + W%values**2.0
    if (spec_rank.EQ.0) print*, 'toto'

    u_prime = computeFieldAvg( BufField , spec_rank , nbcpus )
    if (spec_rank.EQ.0) print*, 'u_prime',u_prime
    BufField%values = B(1)%values**2.0 + B(2)%values**2.0 + B(3)%values**2.0
    b_prime = computeFieldAvg( BufField , spec_rank , nbcpus )
    if (spec_rank.EQ.0) print*, 'b_prime',b_prime
    coef=(u_prime/b_prime)/coef
        do i=1,size(B)
                B(i)%values(:,:,:) = coef*B(i)%values(:,:,:)
        end do
    endif
    ok = .true.

call deleteDataLayout(BufField)
end function initPlaneJet

!---------------------------------------------------------------------------------
!> For jet-plane init: the shape without turbulence is a double hyperbolic
!! tangent.
!! @autor Guillaume Balarac
!
!>    @param[in]    jetVel     = velocity in the jet center
!!    @param[in]    coflowVel  = velocity on the domain boundary
!!    @param[in]    nozzle     = how large is the plane jet ?
!!    @param[in]    stiffness  = stiffnes of the transition between jetVel et
!!                                  coflowVel
!!    @param[in]    coord      = coordinate along the jet direction
!!    @param[out]   vel        = velocity
!!    @param[out]   mask       = mask for the noise part
!---------------------------------------------------------------------------------
subroutine computeJet(jetVel,coflowVel,nozzle,stiffness,coord, mask, vel)

    !I/O data
    real(WP), intent(in):: jetVel,coflowVel,nozzle,stiffness,coord
    real(WP)            :: vel, mask

    real(WP)            :: tanh_part

    tanh_part = 0.5_WP*tanh( stiffness * (1.0_WP - 2.0_WP * abs ( coord ) / nozzle ))
    !vel = ((jetVel + coflowVel)/2.0_WP) - ((jetVel - coflowVel)) *tanh_part
    vel = ((jetVel + coflowVel)/2.0_WP) - ((-jetVel+coflowVel)) *tanh_part

    if ((0.5+tanh_part>0.05).and.(0.5+tanh_part<0.95)) then
        mask = 1
    else
        mask = 0
    end if

end subroutine computeJet

end module planejet
!> @}
