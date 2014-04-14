!> @addtogroup init
!! @{ 
!---------------------------------------------------------------------------------
!> Cosinus initialisation for velocity and scalars.
!! @autor = Guillaume Balarac and Jean-Baptiste Lagaert, LEGI
!
!
!
!---------------------------------------------------------------------------------
module KolmoFlow

  use precision_tools

  implicit none

  public  :: initKolmoFlow

contains

!---------------------------------------------------------------------------------
!> Initialisation with cosinus.
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
function initKolmoFlow(U,V,W,B,VelWN,ScalArray, ScalWN, Scal_partArray, &
        &  Scal_partWN, BfieldWN,nbcpus,spec_rank) result(ok)

    use isotropicturbulence

    implicit none

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
    real(WP)        :: zz,  dz                       
    real(WP)        :: yy,  dy                       
    real(WP)        :: xx,  dx                       
    real(WP)        :: nb
    integer         :: i,j,k
 

    ok = .false.

    ! =====
    dz = U%Lz/real(U%nz,WP)
    ! -- Velocity --
    V%values = 0.0_WP
    W%values = 0.0_WP
    do k = U%zmin,U%zmax
       zz = dz * real(k-1)
       U%values(:,:,k) = 2.0_WP * (sin(zz)+sin(3.0_WP*zz)/9.0_WP)
    end do

    ! === noise
    dx = U%Lx/real(U%nx,WP)
    dy = U%Ly/real(U%ny,WP)
    nb = 2.0
    DO k = U%zmin,U%zmax
       zz = real(k-1,WP) * dz
       if (zz.le.(acos(-1.0)+1.0).or.zz.ge.(acos(-1.0)-1.0)) then
       DO j = U%ymin,U%ymax
          yy = real(j-1,WP) * dy
          DO i = U%xmin,U%xmax
             xx = real(i-1,WP) * dx
             U%values(i,j,k) =  0.02_WP * sin(xx*nb)*sin(yy*nb)*sin(zz*nb) + U%values(i,j,k)
             V%values(i,j,k) =  0.01_WP * cos(xx*nb)*cos(yy*nb)*sin(zz*nb)
             W%values(i,j,k) =  0.01_WP * cos(xx*nb)*sin(yy*nb)*cos(zz*nb)
          ENDDO
       ENDDO
       end if
    ENDDO 
   
    ok = .true.

end function initKolmoFlow

end module KolmoFlow 
!> @}
