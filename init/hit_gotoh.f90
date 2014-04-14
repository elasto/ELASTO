!------------------------------------------------------------------------------
!
! MODULE: isotropicturbulence - initialisation from Gotoh 2012 in JCP
!
!> @author
!> Jean-Baptiste Lagaert and Patrick BEGOU, LEGI
!
!------------------------------------------------------------------------------
MODULE hit_gotoh

  use precision_tools

    USE parser_tools
    USE datalayout
    USE parallel_tools
    USE fileio
    USE transforms_tools
    USE wavenumber_tools
    USE geometric_shape_tools
    USE random_tools

    IMPLICIT NONE

    private

    public  hit_gotoh_init
    private hit_scalar
    private hit_vector

    real(WP), parameter, private :: pi=acos(-1.0_WP)

CONTAINS

!------------------------------------------------------------------------------
!> @author
!> Guillaume Balarac, Patrick BEGOU, LEGI
!
!>
!> @details
!> This subroutine initialise the velocities (U,V,W) and the scalars
!> (stored in the array of scalars called ScalArray) for starting
!> an homogeneous isotropic turbulence simulation.
!> In this first version we set the following conditions to simplify
!> the code:
!>    The computational domain is a cube: lx=ly=lz
!>    The number of points is the same in each direction
!
!> @param[in,out]   U               = longitudinal velocity storage to initialize
!> @param[in,out]   V               = vertical velocity  storage to initialize
!> @param[in,out]   W               = spanwise velocity  storage to initialize
!> @param[in]       VelWN           = the wavenumbers for U,V and W
!> @param[in,out]   B               = magnetic field (3D vector)
!> @param[in,out]   BfieldWN        = wave numbers associated to the magnetic field
!> @param[in,out]   ScalArray       = scalars solved with full pseudo-psectral method (to initialized)
!> @param[in]       ScalWN          = wavenumbers for the scalars
!> @param[in,out]   Scal_P_Array  = scalars solved with full pseudo-psectral method (to initialized)
!> @param[in]       Scal_P_WN     = wavenumbers for the scalars
!> @param[in]       mu              = viscosity
!> @param[in]       prm             = Prandt number
!> @param[in]       nbcpus          = the total amount of processes (or cpus)
!> @param[in]       spec_rank       = my MPI processe rank in the range [0:nbcpus-1]
!> @param[out]      ok              = logical set to .TRUE. if all initialization success.
!> @param[in]       visu_purpose    = logical set to true if scalars have to be initialized to a circle
!!                                    shape in order to produce "nice pictures" for communications.
!------------------------------------------------------------------------------
!=========================================================================================!
SUBROUTINE hit_gotoh_init(U,V,W,Uk,Vk,Wk,VelWN, &
        & ScalArray, ScalArrayk, ScalWN,        &
        & Scal_P_Array, Scal_P_Arrayk,Scal_P_WN,&
        & B, Bk, BfieldWN,mu,Prm,nbcpus,spec_rank,ok, visu_purpose)
!=========================================================================================!


    IMPLICIT NONE
    ! Input/ouput
    type(real_data_layout), intent(inout)           :: U, V, W
    type(complex_data_layout), intent(inout)        :: Uk, Vk, Wk
    TYPE(WaveNumbers), INTENT(IN)                   :: VelWN
    type(real_data_layout), pointer, intent(inout)  :: ScalArray(:), Scal_P_Array(:)
    type(complex_data_layout),dimension(:),intent(inout) :: ScalArrayk(:), Scal_P_Arrayk(:)
    TYPE(WaveNumbers), DIMENSION(:), INTENT(IN)     :: ScalWN, Scal_P_WN
    type(real_data_layout), pointer, intent(inout)  :: B(:)
    type(complex_data_layout),dimension(:),intent(inout) :: Bk(:)
    TYPE(WaveNumbers), INTENT(IN)                   :: BfieldWN
    integer, intent(in), optional                   :: visu_purpose
    REAL(WP), INTENT(IN)                            :: mu
    REAL(WP), INTENT(IN)                            :: Prm
    INTEGER, INTENT(IN)                             :: nbcpus,spec_rank
    LOGICAL, INTENT(OUT)                            :: ok

    ! Local variables
    REAL(WP)                                        :: eta, L_large
    REAL(WP)                                        :: k0,kv,ktheta, uprime,theta_prime
    REAL(WP)                                        :: dk, pi
    REAL(WP)                                        :: e_spec, d_spec, ndense, energy_spec
    REAL(WP)                                        :: kk

    REAL(WP):: Lx,Ly,Lz
    INTEGER :: i,j,k
    REAL(WP):: coef
    !LOGICAL                                         :: circle_init = .false. ! For product "nice" picture, scalar can be init at a circle shape
    integer                                         :: circle_init = 0 ! For product "nice" picture, scalar can be init at a circle shape

    IF (present(visu_purpose)) circle_init = visu_purpose

    pi = acos(-1.0_WP)

    ! ===== Preliminary tests and applying some restrictions on the computational domain =====
    ! Working on a cubic domain
    IF (.NOT. GetLxLyLz("Length of domain",Lx,Ly,Lz)) THEN
        IF(spec_rank .EQ.0) &
            & WRITE(6,'(a)')"[ERROR] hit_init cannot get the dimensions of the domain "
        ok=.FALSE.
        RETURN
    ENDIF
    IF(lx.NE.ly .OR. lx.NE.lz) THEN
        IF(spec_rank .EQ.0) &
        & WRITE(6,'(a)')"[ERROR] hit_init is only implemented for a cubic "// &
        & "domain in THI. "//&
        & "Additional implementation is needed is needed for your problem."
    ok=.FALSE.
    RETURN
    ENDIF    

    ! Working with same number of points in each direction for U
    IF( U%nx .NE. U%ny .OR. U%nx .NE. U%nz) THEN
    IF(spec_rank .EQ.0) &
        & WRITE(6,'(a)')"[ERROR] hit_init is only implemented for a cubic "// &
        & "domain with same number of points in each direction in THI. "//&
        & "Additional implementation is needed for your problem."
    ok=.FALSE.
    RETURN
    ENDIF   

    !  U, V and W should have the same layout
    IF(.NOT. samelayout(U,V,W)) THEN
    IF(spec_rank .EQ.0) &
        & WRITE(6,'(a)')"[INTERNAL ERROR] hit_init: U,V and W are not compatible in size "
    ok=.FALSE.
    RETURN
    ENDIF


    ! ALL seems OK, initialization can start
    dk=2.0_WP*pi/Lx
    call parser_read('k0', k0)
    call parser_read('ktheta', ktheta)
    call parser_read('theta_prime', theta_prime)
    call parser_read('k0', k0)
    call parser_read('kv', kv)
    call parser_read('uprime', uprime)

    IF (spec_rank.eq.0) THEN
        WRITE(6,'(a)')  '======================================================'
        WRITE(6,'(a)')  '|    Isotropic Turbulence Initialisation - Gotoh     |'
        WRITE(6,'(a)')  '======================================================'
        WRITE(6,'("|",a30,1x,":",6x,es14.5,"|")') ' kv ', kv
        WRITE(6,'("|",a30,1x,":",6x,es14.5,"|")') " u' ", uprime
        WRITE(6,'("|",a30,1x,":",6x,es14.5,"|")') ' k0 ', k0
        WRITE(6,'("|",a30,1x,":",6x,es14.5,"|")') ' ktheta ', ktheta
        WRITE(6,'("|",a30,1x,":",6x,es14.5,"|")') " theta_prime' ", theta_prime
        WRITE(6,'(a)')  '======================================================'
    ENDIF
    ok = hit_vector(U, V, W, Uk, Vk, Wk,VelWN, Lx, nbcpus, spec_rank)
!    ok =  hitvelocity_init(U, V, W, VelWN, mu, Lx, nbcpus, spec_rank)


    ! ===== Scalar solved with pseudo-spectral solver =====
    IF (ok .AND. ASSOCIATED(ScalArray)) THEN
        DO i=1,size(ScalArray)
            IF(ScalArray(i)%nx.NE.ScalArray(i)%ny .OR. ScalArray(i)%nx.NE.ScalArray(i)%nz) THEN
                IF(spec_rank .EQ.0) THEN
                    WRITE(6,'(a)')"[ERROR] hit_init is only implemented for a cubic "// &
                        & "domain with same number of points in each direction in THI. "//&
                        & "Additional implementation is needed for your problem."
                    WRITE(6,'(4(a,i0),a)')'Your scalar n. ',i,' is [',ScalArray(i)%nx,',',&
                        & ScalArray(i)%ny,',',ScalArray(i)%nz,'].'
                ENDIF
                ok=.FALSE.
                RETURN
            ENDIF
            IF (circle_init .eq. 1 ) then
                ok=disc_smooth(ScalArray(i), Lx/3.0)
            elseif (circle_init .eq. 2 ) then
                ok=most_smooth(ScalArray(i),spec_rank)
            else
                ok=hit_scalar(ScalArray(i),ScalArrayk(i),ScalWN(i), Lx, nbcpus, spec_rank)
            end if
            IF (.NOT. ok) RETURN
        ENDDO
    ENDIF

    ! ===== Scalar solved with mixed particles/spectral solver =====
    IF (ok .AND. ASSOCIATED(Scal_P_Array)) THEN
        DO i=1,size(Scal_P_Array)
            IF(Scal_P_Array(i)%nx.NE.Scal_P_Array(i)%ny .OR. Scal_P_Array(i)%nx.NE.Scal_P_Array(i)%nz) THEN
                IF(spec_rank .EQ.0) THEN
                    WRITE(6,'(a)')"[ERROR] hit_init is only implemented for a cubic "// &
                        & "domain with same number of points in each direction in THI. "//&
                        & "Additional implementation is needed for your problem."
                    WRITE(6,'(4(a,i0),a)')'Your scalar (solved with particles methos) n. ', &
                        & i,' is [',Scal_P_Array(i)%nx,',',&
                        & Scal_P_Array(i)%ny,',',Scal_P_Array(i)%nz,'].'
                ENDIF
                ok=.FALSE.
                RETURN
            end if
            IF (circle_init .eq. 1) then
                ok=disc_smooth(Scal_P_Array(i), Lx/3.0)
            elseif (circle_init .eq. 2 ) then
                ok=most_smooth(ScalArray(i),spec_rank)
            else
                ok=hit_scalar(Scal_P_Array(i),Scal_P_ArrayK(i),Scal_P_WN(i), Lx, nbcpus, spec_rank)
            end if
            IF (.NOT. ok) RETURN
        ENDDO
    ENDIF


    ! ===== Magnetic field (solved with mixed pseudo-spectral solver) =====
    IF (ok.AND. ASSOCIATED(B)) THEN

       dk=2.0_WP*pi/Lx

       ! Magnetic spectrum as a function of the kinetic energy spectrum
       ok = hit_vector(B(1), B(2), B(3), Bk(1), Bk(2), Bk(3), BfieldWN, &
          &  Lx, nbcpus, spec_rank)
    ENDIF


  end subroutine hit_gotoh_init

!------------------------------------------------------------------------------
!> @author
!> Guillaume Balarac, Patrick BEGOU, LEGI
!
!>
!> @details
!> This function initialize a vector for starting
!> an homogeneous isotropic turbulence simulation.
!
!> @param [in,out] Vec1 first componant of the vector storage to initialize
!> @param [in,out] Vec2 first componant of the vector storage to initialize
!> @param [in,out] Vec3 first componant of the vector storage to initialize
!> @param [in] wavenumb the waves numbers for the velocities (same discretisation for U,V and W)
!> @param [in] Lx the lenght of the domain (supposed to be the same in the 3 dimensions).
!> @param [in] nbcpus the total amount of cpus used by the application
!> @return .TRUE. if initialization is successfull.
!------------------------------------------------------------------------------
!=========================================================================================!
  function hit_vector(Vec1, Vec2, Vec3, Vec1k, Vec2k, Vec3k, wavenumb, Lx, nbcpus, spec_rank)

    implicit none
    ! Input/Ouptut
    type(real_data_layout), intent(inout) :: Vec1
    type(real_data_layout), intent(inout) :: Vec2
    type(real_data_layout), intent(inout) :: Vec3
    type(complex_data_layout), intent(inout) :: Vec1k
    type(complex_data_layout), intent(inout) :: Vec2k
    type(complex_data_layout), intent(inout) :: Vec3k
    type(WaveNumbers), intent(in)::wavenumb
    real(WP), intent(in) :: Lx
    integer, intent(in)  :: nbcpus, spec_rank
    logical :: hit_vector
    ! Local variables
    type(complex_data_layout) :: ak
    type(complex_data_layout) :: bk
    real(WP) :: dk, kc, eps
    real(WP) :: energy_spec
    real(WP) :: kk, kk2
    !real(WP) :: amp_disc
    real(WP), DIMENSION(:,:), ALLOCATABLE :: spect
    real(WP) :: ps1, ps2, psr, racsg
    integer  :: nk, ik
    integer  :: i,j,k,ierr
    complex(WP), PARAMETER :: ii=(0.0_WP,1.0_WP)
    complex(WP), DIMENSION(:,:),ALLOCATABLE :: firstslice
    logical :: res1,res2,res3
    real(WP) :: uprime, kv, k0
    integer                                 :: file_id          ! id for file io
    character(len=50)                       :: file_name        ! name of output file


    hit_vector=.FALSE.

    nk=Vec1%nx/2 + 1     ! a changer si nx != ny != nz
    dk=2.0_WP*pi/Lx
    call parser_read('k0', k0)
    call parser_read('kv', kv)
    call parser_read('uprime', uprime)
    !amp_disc = sqrt(dk)**3
    kc=pi*REAL(Vec1%nx,WP)/Lx
    eps=kc/10000000_WP
    !eps=kc/100_WP

    !--------------------------------------------------
    ! Instantiation of the arrays in the spectral domain
    !--------------------------------------------------
    IF (.NOT. initDataLayout("Ak_Spectral",ak,(Vec1%nx/2)+1,Vec1%ny,Vec1%nz,Lx,Lx,Lx,nbcpus,spec_rank,alongZ)) THEN
       WRITE(6,'(a,i0)')'[ERROR] hit_vector Failed to init ak on ', spec_rank
       RETURN
    ENDIF
    IF (.NOT. initDataLayout("Bk_Spectral",bk,(Vec1%nx/2)+1,Vec1%ny,Vec1%nz,Lx,Lx,Lx,nbcpus,spec_rank,alongZ)) THEN
       WRITE(6,'(a,i0)')'[ERROR] hit_vector Failed to init bk on ', spec_rank
       RETURN
    ENDIF

    ! Output spectrum for comparison
    allocate(spect(Vec1%nx+1,2),STAT=ierr) !WARNING spect array shape is different
    if(ierr.ne.0) then
        write(6,'()')'[ERROR] hit_vector: Not enought memory for spect on processor ',spec_rank
        return
    end if
    spect(:,1) = (/(REAL(i-1,WP)*dk,i=1,size(spect,1))/)
    spect(:,2) = 0.0_WP

    ! Compute spectrum
    do k=ak%zmin,ak%zmax
      do j=ak%ymin,ak%ymax
        do i=ak%xmin,ak%xmax
          kk=sqrt(wavenumb%kx(i)*wavenumb%kx(i)+wavenumb%ky(j)*wavenumb%ky(j)+wavenumb%kz(k)*wavenumb%kz(k))

          ! Random numbers
          call random_number(psr)
          psr= 2.0_WP*pi*psr                !2.0_WP*pi*(rand-0.5_WP)
          call random_number(ps1)
          ps1=2.0_WP*pi*ps1                 !*(rand-0.5_WP)
          call random_number(ps2)
          ps2=2.0_WP*pi*ps2                 !(rand-0.5_WP)

          ! Spectrums
          !energy_spec=spec_amp*(kk/ke)**4*exp(-2.0_WP*(kk/ke)**2)
          energy_spec=16._WP*sqrt(2._WP/pi)*(uprime**2)*((kk**4)/(k0**5))*exp(-2.0_WP*(kk/kv)**2)

          ! Coeff
          ik = 1+idint(kk/dk + 0.5_WP)
          if ((ik<=nk).and.(kk<kc).and.(kk>1e-10)) then
            ak%values(i,j,k)= sqrt(energy_spec/(4.0_WP*pi*kk**2)) &
              &                      * exp(ii*ps1) * cos(psr)
            bk%values(i,j,k)= sqrt(energy_spec/(4.0_WP*pi*kk**2)) &
              &                      * exp(ii*ps2) * sin(psr)
            spect(ik,2) = spect(ik,2) + real(ak%values(i,j,k)*conjg(ak%values(i,j,k))) &
              &                     + real(bk%values(i,j,k)*conjg(bk%values(i,j,k)))
          else
            ak%values(i,j,k) = 0._WP
            bk%values(i,j,k) = 0._WP
          end if
        end do
      end do
    end do

    ! Get imposed spectrum
    spect(:,2)=DoGlobalVectorSum(spec_rank,spect(:,2))
    if (spec_rank==0) then
      file_id = iopen()
      write(file_name,'(a)') 'spec_vel_th.table'
      file_name = trim(adjustl(file_name))
      open(unit=file_id, file=file_name, form="FORMATTED")
      write(file_id,'(a)') '# Kinetic spectrum'
      write(file_id,'(a)') '# wave_number and spectrum '
      do ik = 1, Vec1%nx+1
        write(file_id, '(e14.6,x,e14.6)') spect(ik,:)
      end do
      close(file_id)
      file_id = iclose(file_id)
    end if

    do k = Vec1k%zmin,Vec1k%zmax
      do j = Vec1k%ymin,Vec1k%ymax
        do i = Vec1k%xmin,Vec1k%xmax
          kk=sqrt(wavenumb%kx(i)**2+wavenumb%ky(j)**2+wavenumb%kz(k)**2)
          kk2=sqrt(wavenumb%kx(i)**2+wavenumb%ky(j)**2)
          if (kk2*kk>eps) then
            Vec1k%values(i,j,k)=(ak%values(i,j,k)*kk*wavenumb%ky(j)+bk%values(i,j,k)*wavenumb%kx(i) &
              &            *wavenumb%kz(k))/(kk*kk2)
            Vec2k%values(i,j,k)=(bk%values(i,j,k)*wavenumb%ky(j)*wavenumb%kz(k)-ak%values(i,j,k) &
              &            *kk*wavenumb%kx(i))/(kk*kk2)
            Vec3k%values(i,j,k)=-bk%values(i,j,k)*kk2/kk
          else
            Vec1k%values(i,j,k)=(ak%values(i,j,k)+bk%values(i,j,k))/sqrt(2.0_WP)
            Vec2k%values(i,j,k)=(bk%values(i,j,k)-ak%values(i,j,k))/sqrt(2.0_WP)
            Vec3k%values(i,j,k)= 0.0_WP
          end if
          ik = 1+idint(kk/dk + 0.5_WP)
        end do
      end do
    end do

    call btran(Vec1K,Vec1,res1)
    call btran(Vec2K,Vec2,res2)
    call btran(Vec3K,Vec3,res3)

    !--------------------------------------------------
    ! release allocated memory for local variables
    !--------------------------------------------------
    deallocate(spect)
    call deleteDataLayout(ak)
    call deleteDataLayout(bk)

    hit_vector=(res1.AND.res2.AND.res3)

  end function hit_vector


!------------------------------------------------------------------------------
!> @author
!> Guillaume Balarac, Patrick BEGOU, LEGI
!
!>
!> @details
!> This function initialise a scalar for starting
!> an homogeneous isotropic turbulence simulation.
!
!> @param[in,out] physcal   = scalar storage to initialize (in physical space)
!> @param[in,out] fouscal   = scalar storage to initialize (in Fourrier space)
!> @param[in]     wavenumb  = waves numbers for this scalar
!> @param[in]     Lx        = length of the domain. It recquieres that Lx==Ly==Lz.
!> @param [in]    nbcpus    = total amount of cpus used by the application
!> @param [in]    spec_rank = my MPI processe rank in the spectral communicator.
!> @return .TRUE. if initialization is successfull.
!------------------------------------------------------------------------------
!=========================================================================================!
  function hit_scalar(physcal, fouscal, wavenumb, lx, nbcpus, spec_rank)

    implicit none
    ! Input/Output
    type(real_data_layout), intent(inout) :: physcal
    type(complex_data_layout), intent(inout) :: fouscal
    type(wavenumbers), intent(in)::wavenumb
    integer, intent(in)  :: nbcpus, spec_rank
    real(wp), intent(in) :: lx
    logical :: hit_scalar
    ! Local variables
    real(wp) :: dk
    real(wp) :: k0, ktheta, theta_prime
    real(wp) :: kk, rand, energy_spec
    complex(wp), parameter :: ii=(0.0_wp,1.0_wp)
    integer  :: nk
    integer  :: i,j,k, ik
    logical :: res

    hit_scalar=.FALSE.

    nk=physcal%nx/2 + 1     ! a changer si nx != ny != nz
    dk=2.0_WP*pi/Lx
    call parser_read('k0', k0)
    call parser_read('ktheta', ktheta)
    call parser_read('theta_prime', theta_prime)

    fouscal%values = 0.0_WP
    do k = fouscal%zmin,fouscal%zmax
      do j = fouscal%ymin,fouscal%ymax
        do i = fouscal%xmin,fouscal%xmax
          kk=sqrt(wavenumb%kx(i)**2+wavenumb%ky(j)**2+wavenumb%kz(k)**2)
          ik = 1+idint(kk/dk + 0.5_WP)
          !if (ik<=nk) then
          if ((kk<=70).and.(kk>1e-10)) then
            !f_phi = 1.0_wp
            call random_number(rand)
            !energy_spec = (32._WP/3._WP)*sqrt(2._WP/pi)*(theta_prime**2)&    !
            !Gotoh says he uses a coefficicient equal to 32/3, but statistical values show that
            !it is a 32 coefficient.
            energy_spec = (32._WP/3._WP)*sqrt(2._WP/pi)*(theta_prime**2)&
                & *((kk**4)/(k0**5))*exp(-2.0_WP*(kk/ktheta)**2)
            fouscal%values(i,j,k) = sqrt(energy_spec/(4.0_WP*pi*kk**2))*exp(ii*2.0_wp*pi*rand)
          end if
        end do
      end do
    end do


    call btran(fouscal,physcal,res)
    if(.not. res) then
        write(6,'(a,i0)')'[ERROR] hit_scalar Failed in fouscal backsward FFT on ', spec_rank
        return
    endif

    hit_scalar=.TRUE.

  end function hit_scalar

end module hit_gotoh
