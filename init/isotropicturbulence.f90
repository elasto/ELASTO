!------------------------------------------------------------------------------
!
! MODULE: isotropicturbulence
!
!> @author
!> Patrick BEGOU, LEGI
!
!------------------------------------------------------------------------------
MODULE isotropicturbulence
    USE parser_tools
    USE datalayout
    USE parallel_tools
    USE fileio
    USE transforms_tools
    USE wavenumber_tools 
    USE geometric_shape_tools
    USE random_tools 

    IMPLICIT NONE

    PUBLIC hit_init
    PUBLIC hitvelocity_init
    PUBLIC hitscalar_init
    PUBLIC hit_vector_init

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
!> @param[in,out]   Scal_partArray  = scalars solved with full pseudo-psectral method (to initialized)
!> @param[in]       Scal_partWN     = wavenumbers for the scalars
!> @param[in]       mu              = viscosity
!> @param[in]       prm             = Prandt number
!> @param[in]       nbcpus          = the total amount of processes (or cpus)
!> @param[in]       spec_rank       = my MPI processe rank in the range [0:nbcpus-1]
!> @param[out]      ok              = logical set to .TRUE. if all initialization success.
!> @param[in]       visu_purpose    = logical set to true if scalars have to be initialized to a circle
!!                                    shape in order to produce "nice pictures" for communications.
!------------------------------------------------------------------------------
!=========================================================================================!
SUBROUTINE hit_init(U,V,W,B,VelWN,ScalArray, ScalWN, Scal_partArray, &
        &  Scal_partWN, BfieldWN,mu,Prm,nbcpus,spec_rank,ok, visu_purpose)
!=========================================================================================!


    IMPLICIT NONE
    ! Input/ouput
    TYPE(REAL_DATA_LAYOUT), INTENT(INOUT)           :: U
    TYPE(REAL_DATA_LAYOUT), INTENT(INOUT)           :: V
    TYPE(REAL_DATA_LAYOUT), INTENT(INOUT)           :: W
    TYPE(REAL_DATA_LAYOUT), POINTER, INTENT(INOUT)  :: B(:)
    TYPE(REAL_DATA_LAYOUT), POINTER, INTENT(INOUT)  :: ScalArray(:)
    TYPE(REAL_DATA_LAYOUT), POINTER, INTENT(INOUT)  :: Scal_partArray(:)
    TYPE(WaveNumbers), INTENT(IN)                   :: VelWN
    TYPE(WaveNumbers), DIMENSION(:), INTENT(IN)     :: ScalWN
    TYPE(WaveNumbers), DIMENSION(:), INTENT(IN)     :: Scal_partWN
    TYPE(WaveNumbers), INTENT(IN)                   :: BfieldWN
    !LOGICAL, intent(in), optional                   :: visu_purpose
    integer, intent(in), optional                   :: visu_purpose
    REAL(WP), INTENT(IN)                            :: mu
    REAL(WP), INTENT(IN)                            :: Prm
    INTEGER, INTENT(IN)                             :: nbcpus,spec_rank
    LOGICAL, INTENT(OUT)                            :: ok

    ! Local variables
    REAL(WP)                                        :: keta, Rlambda
    REAL(WP)                                        :: eta, L_large
    REAL(WP)                                        :: L11, ke, pi
    REAL(WP)                                        :: dk, kc, eps
    REAL(WP)                                        :: e_spec, d_spec, ndense, energy_spec
    REAL(WP)                                        :: kk
    REAL(WP)                                        :: spec_amp, amp_disc, tke, Ut, dissipation

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
    ! ===== Velocity field =====
    CALL parser_read('Kmax eta', keta)
    CALL parser_read('R Lambda', Rlambda)
    eta=(2.0_WP*keta)/REAL(U%nx,WP)
    L_large=eta*(3.0_WP/20.0_WP*Rlambda*Rlambda)**(3.0/4.0)
    L11=L_large*0.43
    ke=1.3_WP/L11
    IF (ke.lt.1) THEN
        IF (spec_rank.EQ.0) &
            & WRITE(6,'(a)')'[ERROR] Peak wavenumber not in the box, reduce R Lambda'
        ok=.FALSE.
        RETURN
    ENDIF
    dk=2.0_WP*pi/Lx
    kc=pi*REAL(U%nx,WP)/Lx
    eps=kc/10000000_WP
    e_spec=0.0_WP
    d_spec=0.0_WP
    ndense=0.0_WP

    DO k = 1,SIZE(velWN%kz)
       DO j = 1,SIZE(velWN%ky)
          DO i = 1,SIZE(velWN%kx)
             kk=sqrt(velWN%kx(i)*velWN%kx(i)+velWN%ky(j)*velWN%ky(j)+velWN%kz(k)*velWN%kz(k))
             IF ((kk.gt.eps).and.(kk.le.kc)) THEN
                ndense = ndense + 1.0_WP
                energy_spec = (kk/ke)**4*exp(-2.0_WP*(kk/ke)**2)/(4.0_WP*pi*kk**2)
                e_spec = e_spec + dk*energy_spec
                d_spec = d_spec + dk*kk*kk*energy_spec
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    !  mu = sqrt(e_spec**2/d_spec*(10.0_WP/(3.0_WP*Rlambda**2))) 
    !  tke = 1.5_WP*Ut*Ut
    spec_amp = 3.0_WP/10.0_WP*Rlambda**2*mu**2*d_spec/e_spec**2
    tke = spec_amp*e_spec
    Ut = sqrt(tke/1.5_WP)   !Rlambda*mu/(sqrt(10.0_WP)*eta**(2.0_WP/3.0_WP)*L_large**(1.0_WP/3.0_WP))
    dissipation = 2.0_WP*mu*spec_amp*d_spec
    amp_disc = sqrt(dk)**3

    IF (spec_rank.eq.0) THEN
        WRITE(6,'(a)')  '======================================================'
        WRITE(6,'(a)')  '|        Isotropic Turbulence Initialisation         |'
        WRITE(6,'(a)')  '======================================================'
        WRITE(6,'("|",a30,1x,":",6x,es14.5,"|")') ' L11 ', L11
        WRITE(6,'("|",a30,1x,":",6x,es14.5,"|")') ' L_large ', L_large
        WRITE(6,'("|",a30,1x,":",6x,es14.5,"|")') ' kc ', kc
        WRITE(6,'("|",a30,1x,":",6x,es14.5,"|")') ' ke ', ke
        WRITE(6,'("|",a30,1x,":",6x,es14.5,"|")') ' eta ', eta
        WRITE(6,'("|",a30,1x,":",6x,es14.5,"|")') ' Ut  ', Ut
        WRITE(6,'("|",a30,1x,":",6x,es14.5,"|")') ' kinetic energy ', tke
        WRITE(6,'("|",a30,1x,":",6x,es14.5,"|")') ' dissipation ', dissipation
        WRITE(6,'("|",a30,1x,":",6x,es14.5,"|")') ' spec_amp ', spec_amp
        WRITE(6,'("|",a30,1x,":",6x,es14.5,"|")') ' amp_disc ', amp_disc
        WRITE(6,'("|",a30,1x,":",6x,es14.5,"|")') ' viscosity', mu
        WRITE(6,'("|",a30,1x,":",6x,es14.5,"|")') ' e_spec ', e_spec
        WRITE(6,'("|",a30,1x,":",6x,es14.5,"|")') ' d_spec ', d_spec
        WRITE(6,'("|",a30,1x,":",6x,es14.5,"|")') ' Number of modes', ndense
        WRITE(6,'("|",a30,1x,":",6x,es14.5,"|")') ' R lambda calculated', &
                                & sqrt(20.0_WP/3.0_WP*tke**2.0/(dissipation*mu))
        WRITE(6,'(a)')  '======================================================'
    ENDIF
    ok = hit_vector_init(U, V, W, VelWN, Lx, nbcpus, spec_rank,ke,spec_amp)
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
                ok=hitscalar_init(ScalArray(i),ScalWN(i), Lx, nbcpus, spec_rank)
            end if
            IF (.NOT. ok) RETURN
        ENDDO
    ENDIF

    ! ===== Scalar solved with mixed particles/spectral solver =====
    IF (ok .AND. ASSOCIATED(Scal_partArray)) THEN
        DO i=1,size(Scal_partArray)
            IF(Scal_partArray(i)%nx.NE.Scal_partArray(i)%ny .OR. Scal_partArray(i)%nx.NE.Scal_partArray(i)%nz) THEN
                IF(spec_rank .EQ.0) THEN
                    WRITE(6,'(a)')"[ERROR] hit_init is only implemented for a cubic "// &
                        & "domain with same number of points in each direction in THI. "//&
                        & "Additional implementation is needed for your problem."
                    WRITE(6,'(4(a,i0),a)')'Your scalar (solved with particles methos) n. ', &
                        & i,' is [',Scal_partArray(i)%nx,',',&
                        & Scal_partArray(i)%ny,',',Scal_partArray(i)%nz,'].'
                ENDIF
                ok=.FALSE.
                RETURN
            ENDIF   
            IF (circle_init .eq. 1) then
                ok=disc_smooth(Scal_partArray(i), Lx/3.0)
            elseif (circle_init .eq. 2 ) then
                ok=most_smooth(ScalArray(i),spec_rank)
            else
                ok=hitscalar_init(Scal_partArray(i),Scal_partWN(i), Lx, nbcpus, spec_rank)
            end if
            IF (.NOT. ok) RETURN
        ENDDO
    ENDIF


    ! ===== Magnetic field (solved with mixed pseudo-spectral solver) =====
    IF (ok.AND. ASSOCIATED(B)) THEN
       CALL parser_read('Magnetic spectrum peak', ke)
       CALl parser_read('Coef Ek/Em', coef)

       dk=2.0_WP*pi/Lx
       kc=pi*REAL(B(1)%nx,WP)/Lx
       eps=kc/10000000_WP
       e_spec=0.0_WP

       DO k = 1,SIZE(BfieldWN%kz)
          DO j = 1,SIZE(BfieldWN%ky)
             DO i = 1,SIZE(BfieldWN%kx)
                kk=sqrt(BfieldWN%kx(i)*BfieldWN%kx(i)+BfieldWN%ky(j)*BfieldWN%ky(j)+BfieldWN%kz(k)*BfieldWN%kz(k))
                IF ((kk.gt.eps).and.(kk.le.kc)) THEN
                   energy_spec = (kk/ke)**4*exp(-2.0_WP*(kk/ke)**2)/(4.0_WP*pi*kk**2)
                   e_spec = e_spec + dk*energy_spec
                ENDIF
             ENDDO
          ENDDO
       ENDDO

       ! Magnetic spectrum as a function of the kinetic energy spectrum 
       spec_amp = tke / (coef*e_spec)
       ok = hit_vector_init(B(1), B(2), B(3), BfieldWN, Lx, nbcpus, spec_rank,ke,spec_amp)
    ENDIF

    RETURN  

END SUBROUTINE hit_init

!------------------------------------------------------------------------------
!> @author 
!> Guillaume Balarac, Patrick BEGOU, LEGI
!
!>
!> @details
!> This function initialise the velocities (U,V,W) for starting
!> an homogeneous isotropic turbulence simulation.
!
!> @param [in,out] U longitudinal velocity storage to initialize
!> @param [in,out] V vertical velocity  storage to initialize
!> @param [in,out] W spanwise velocity  storage to initialize
!> @param [in] wavenumb the waves numbers for the velocities (same discretisation for U,V and W)
!> @param [in] Lx the lenght of the domain (supposed to be the same in the 3 dimensions).
!> @param [in] mu viscosity
!> @param [in] nbcpus the total amount of cpus used by the application
!> @param [in] spec_rank my MPI processe rank in the range [0:ncpus-1]
!> @return .TRUE. if initialization is successfull. 
!------------------------------------------------------------------------------
!=========================================================================================!
FUNCTION hitvelocity_init(U, V, W, wavenumb, mu, Lx, nbcpus, spec_rank)
!=========================================================================================!
    
    IMPLICIT NONE
    TYPE(REAL_DATA_LAYOUT), INTENT(INOUT) :: U
    TYPE(REAL_DATA_LAYOUT), INTENT(INOUT) :: V
    TYPE(REAL_DATA_LAYOUT), INTENT(INOUT) :: W
    TYPE(WaveNumbers), INTENT(IN)::wavenumb
    REAL(WP), INTENT(IN) :: mu, Lx
    INTEGER, INTENT(IN)  :: nbcpus, spec_rank
    LOGICAL :: hitvelocity_init
    
    TYPE(COMPLEX_DATA_LAYOUT) :: Uk
    TYPE(COMPLEX_DATA_LAYOUT) :: Vk
    TYPE(COMPLEX_DATA_LAYOUT) :: Wk
    TYPE(COMPLEX_DATA_LAYOUT) :: ak
    TYPE(COMPLEX_DATA_LAYOUT) :: bk
    REAL(WP) :: pi
    REAL(WP) :: keta, Rlambda           !provided in input file
    REAL(WP) :: eta, L_large
    REAL(WP) :: L11, ke
    REAL(WP) :: dk, kc, eps 
    REAL(WP) :: e_spec, d_spec, ndense, energy_spec
!   REAL(WP) :: kx,ky,kz,kk, kk2
    REAL(WP) :: kk, kk2
    REAL(WP) :: spec_amp, amp_disc, tke, Ut, dissipation
    REAL(WP), DIMENSION(:,:), ALLOCATABLE :: spect
    REAL(WP), DIMENSION(:), ALLOCATABLE :: sg
    REAL(WP) :: e_total, diss_total
    REAL(WP) :: ps1, ps2, psr, racsg
    INTEGER  :: nk, ik
    INTEGER  :: i,j,k,ierr, iunit
    INTEGER, DIMENSION(:,:,:),ALLOCATABLE :: node_index 
    COMPLEX(WP), PARAMETER :: ii=(0.0_WP,1.0_WP)
    COMPLEX(WP), DIMENSION(:,:),ALLOCATABLE :: firstslice
    LOGICAL :: res1,res2,res3

    hitvelocity_init=.FALSE.
    
    !--------------------------------------------------
    ! Create pi number
    !--------------------------------------------------
    pi=acos(-1.0_WP)
    nk=U%nx/2 + 1     ! a changer si nx != ny != nz
    
    !--------------------------------------------------
    ! Instantiation of the arrays in the spectral domain
    !--------------------------------------------------
    IF (.NOT. initDataLayout("U_Spectral",Uk,(U%nx/2)+1,U%ny,U%nz,Lx,Lx,Lx,nbcpus,spec_rank,alongZ)) THEN
       WRITE(6,'(a,i0)')'[ERROR] hitvelocity_init Failed to init Uk on ', spec_rank
       RETURN
    ENDIF
    IF (.NOT. initDataLayout("V_Spectral",Vk,(V%nx/2)+1,V%ny,V%nz,Lx,Lx,Lx,nbcpus,spec_rank,alongZ)) THEN
       WRITE(6,'(a,i0)')'[ERROR] hitvelocity_init Failed to init Vk on ', spec_rank
       RETURN
    ENDIF
    IF (.NOT. initDataLayout("W_Spectral",Wk,(W%nx/2)+1,W%ny,W%nz,Lx,Lx,Lx,nbcpus,spec_rank,alongZ)) THEN
       WRITE(6,'(a,i0)')'[ERROR] hitvelocity_init Failed to init Wk on ', spec_rank
       RETURN
    ENDIF
    IF (.NOT. initDataLayout("Ak_Spectral",ak,(U%nx/2)+1,U%ny,U%nz,Lx,Lx,Lx,nbcpus,spec_rank,alongZ)) THEN
       WRITE(6,'(a,i0)')'[ERROR] hitvelocity_init Failed to init ak on ', spec_rank
       RETURN
    ENDIF
    IF (.NOT. initDataLayout("Bk_Spectral",bk,(U%nx/2)+1,U%ny,U%nz,Lx,Lx,Lx,nbcpus,spec_rank,alongZ)) THEN
       WRITE(6,'(a,i0)')'[ERROR] hitvelocity_init Failed to init bk on ', spec_rank
       RETURN
    ENDIF

    CALL parser_read('Kmax eta', keta)
    CALL parser_read('R Lambda', Rlambda)
    
    eta=(2.0_WP*keta)/REAL(U%nx,WP)
    L_large=eta*(3.0_WP/20.0_WP*Rlambda*Rlambda)**(3.0/4.0)
    L11=L_large*0.43
    ke=1.3_WP/L11
    IF (ke.lt.1) THEN
       IF (spec_rank.EQ.0) &
       & WRITE(6,'(a)')'[ERROR] Peak wavenumber not in the box, reduce R Lambda'
       RETURN
    ENDIF
  
    dk=2.0_WP*pi/Lx
    amp_disc = sqrt(dk)**3
    kc=pi*REAL(U%nx,WP)/Lx
    eps=kc/10000000_WP
    
    e_spec=0.0_WP
    d_spec=0.0_WP
    ndense=0.0_WP

    DO k = Uk%zmin,Uk%zmax
       DO j = Uk%ymin,Uk%ymax
          DO i = Uk%xmin,Uk%xmax
             kk=sqrt(wavenumb%kx(i)*wavenumb%kx(i)+wavenumb%ky(j)*wavenumb%ky(j)+wavenumb%kz(k)*wavenumb%kz(k))
             IF ((kk.gt.eps).and.(kk.le.kc)) THEN
                ndense = ndense + 1.0_WP
                energy_spec = (kk/ke)**4*exp(-2.0_WP*(kk/ke)**2)/(4.0_WP*pi*kk**2)
                e_spec = e_spec + dk*energy_spec
                d_spec = d_spec + dk*kk*kk*energy_spec
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    !sums on every nodes
    e_spec=DoGlobalSum(spec_rank,e_spec)
    d_spec=DoGlobalSum(spec_rank,d_spec)
    ndense=DoGlobalSum(spec_rank,ndense)
    
    !  mu = sqrt(e_spec**2/d_spec*(10.0_WP/(3.0_WP*Rlambda**2))) 
    !  tke = 1.5_WP*Ut*Ut
    spec_amp = 3.0_WP/10.0_WP*Rlambda**2*mu**2*d_spec/e_spec**2
    tke = spec_amp*e_spec
    Ut = sqrt(tke/1.5_WP)   !Rlambda*mu/(sqrt(10.0_WP)*eta**(2.0_WP/3.0_WP)*L_large**(1.0_WP/3.0_WP))
    dissipation = 2.0_WP*mu*spec_amp*d_spec

    IF (spec_rank.eq.0) THEN
        WRITE(6,'(a)')  '======================================================'
        WRITE(6,'(a)')  '|        Isotropic Turbulence Initialisation         |'
        WRITE(6,'(a)')  '======================================================'
        WRITE(6,'("|",a30,1x,":",6x,es14.5,"|")') ' L11 ', L11
        WRITE(6,'("|",a30,1x,":",6x,es14.5,"|")') ' L_large ', L_large
        WRITE(6,'("|",a30,1x,":",6x,es14.5,"|")') ' kc ', kc
        WRITE(6,'("|",a30,1x,":",6x,es14.5,"|")') ' ke ', ke
        WRITE(6,'("|",a30,1x,":",6x,es14.5,"|")') ' eta ', eta
        WRITE(6,'("|",a30,1x,":",6x,es14.5,"|")') ' Ut  ', Ut
        WRITE(6,'("|",a30,1x,":",6x,es14.5,"|")') ' kinetic energy ', tke
        WRITE(6,'("|",a30,1x,":",6x,es14.5,"|")') ' dissipation ', dissipation
        WRITE(6,'("|",a30,1x,":",6x,es14.5,"|")') ' spec_amp ', spec_amp
        WRITE(6,'("|",a30,1x,":",6x,es14.5,"|")') ' amp_disc ', amp_disc
        WRITE(6,'("|",a30,1x,":",6x,es14.5,"|")') ' viscosity', mu
        WRITE(6,'("|",a30,1x,":",6x,es14.5,"|")') ' e_spec ', e_spec
        WRITE(6,'("|",a30,1x,":",6x,es14.5,"|")') ' d_spec ', d_spec
        WRITE(6,'("|",a30,1x,":",6x,es14.5,"|")') ' Number of modes', ndense
        WRITE(6,'("|",a30,1x,":",6x,es14.5,"|")') ' R lambda calculated', &
                                & sqrt(20.0_WP/3.0_WP*tke**2.0/(dissipation*mu))
        WRITE(6,'(a)')  '======================================================'
    ENDIF


    ! Output spectrum for comparison
    ALLOCATE(spect(U%nx+1,2),STAT=ierr) !WARNING spect array shape is different
    IF(ierr.NE.0) THEN
        WRITE(6,'()')'[ERROR] hitvelocity_init: Not enought memory for spect on processor ',spec_rank
        RETURN
    ENDIF
    spect(:,1) = (/(REAL(i-1,WP)*dk,i=1,size(spect,1))/)
    spect(:,2) = 0.0_WP
    e_total=0.0_WP
    diss_total = 0.0_WP
    ndense = 0.0_WP

    ! Compute spectrum
    DO k=ak%zmin,ak%zmax
        DO j=ak%ymin,ak%ymax
            DO i=ak%xmin,ak%xmax
                kk=sqrt(wavenumb%kx(i)*wavenumb%kx(i)+wavenumb%ky(j)*wavenumb%ky(j)+wavenumb%kz(k)*wavenumb%kz(k))

                ! Random numbers
                CALL random_number(psr)
                psr= 2.0_WP*pi*psr                !2.0_WP*pi*(rand-0.5_WP)
                CALL random_number(ps1)
                ps1=2.0_WP*pi*ps1                 !*(rand-0.5_WP)
                CALL random_number(ps2)
                ps2=2.0_WP*pi*ps2                 !(rand-0.5_WP)

                ! Spectrums
                energy_spec=spec_amp*(kk/ke)**4*exp(-2.0_WP*(kk/ke)**2)

                ! Coeff
                ik = 1+idint(kk/dk + 0.5_WP)
                IF ((kk.gt.eps).and.(kk.le.kc)) THEN
                    ndense = ndense + 1.0_WP
                    ak%values(i,j,k)= amp_disc * sqrt(energy_spec/(4.0_WP*pi*kk**2)) &
                        &                      * exp(ii*ps1) * cos(psr)
                    bk%values(i,j,k)= amp_disc * sqrt(energy_spec/(4.0_WP*pi*kk**2)) &
                        &                      * exp(ii*ps2) * sin(psr)
                    spect(ik,2) = spect(ik,2) + REAL(ak%values(i,j,k)*conjg(ak%values(i,j,k))) &
                        &                     + REAL(bk%values(i,j,k)*conjg(bk%values(i,j,k)))
                    e_total = e_total + dk*(  REAL(ak%values(i,j,k)*conjg(ak%values(i,j,k)))   &
                        &                   + REAL(bk%values(i,j,k)*conjg(bk%values(i,j,k))))
                    diss_total = diss_total + dk*2.0_WP*mu*kk**2               &
                        &        *( REAL(ak%values(i,j,k)*conjg(ak%values(i,j,k)))  &
                        &        + REAL(bk%values(i,j,k)*conjg(bk%values(i,j,k))))
                ENDIF
            ENDDO
        ENDDO
    ENDDO
    e_total=DoGlobalSum(spec_rank,e_total)
    diss_total=DoGlobalSum(spec_rank,diss_total)
    spect(:,2)=DoGlobalVectorSum(spec_rank,spect(:,2))


    IF (spec_rank .EQ. 0) THEN
        iunit=iopen()
        OPEN(iunit,file='spectrum.analytic',form='formatted',err=100,iostat=ierr)
        REWIND(iunit)
        DO i=1,size(spect,1)
            IF ((spect(i,1).ne.0.0_WP).AND.(spect(i,2).ne.0.0_WP)) THEN
                WRITE(iunit,*,err=110,iostat=ierr) spect(i,:)  
            ENDIF
        ENDDO
        WRITE(6,'(a)')  '======================================================'
        WRITE(6,'(a)')  '|       "spectrum.analytic" file created             |'
        WRITE(6,'(a)')  '======================================================'
        WRITE(6,'("|",a30,1x,":",6x,es14.5,"|")') ' Total energy      ', e_total
        WRITE(6,'("|",a30,1x,":",6x,es14.5,"|")') ' Total dissipation ',diss_total
        WRITE(6,'(a)')  '======================================================'
110     IF (ierr .NE. 0) THEN
            WRITE(6,'(a)')'[WARNING] hitvelocity_init:' // &
                &' unable to write in spectrum.analytic file' 
            ierr=0     
        ENDIF
        CLOSE(iclose(iunit))
100     IF (ierr .NE. 0) WRITE(6,'(a)')'[WARNING] hitvelocity_init:' // &
           &' unable to create spectrum.analytic file'      
        ENDIF

    DO k = Uk%zmin,Uk%zmax
        DO j = Uk%ymin,Uk%ymax
            DO i = Uk%xmin,Uk%xmax
                kk=sqrt(wavenumb%kx(i)*wavenumb%kx(i)+wavenumb%ky(j)*wavenumb%ky(j)+wavenumb%kz(k)*wavenumb%kz(k))
                IF ((kk.GT.eps).and.(kk.LE.kc)) THEN
                    kk2=sqrt(wavenumb%kx(i)*wavenumb%kx(i)+wavenumb%ky(j)*wavenumb%ky(j))
                    if (kk2.lt.eps) then
                        Uk%values(i,j,k)=(ak%values(i,j,k)+bk%values(i,j,k))/sqrt(2.0_WP)
                        Vk%values(i,j,k)=(bk%values(i,j,k)-ak%values(i,j,k))/sqrt(2.0_WP)
                    else
                        Uk%values(i,j,k)=(ak%values(i,j,k)*kk*wavenumb%ky(j)+bk%values(i,j,k)*wavenumb%kx(i) &
                            &            *wavenumb%kz(k))/(kk*kk2)
                        Vk%values(i,j,k)=(bk%values(i,j,k)*wavenumb%ky(j)*wavenumb%kz(k)-ak%values(i,j,k) &
                            &            *kk*wavenumb%kx(i))/(kk*kk2)
                    end if
                    Wk%values(i,j,k)=-bk%values(i,j,k)*kk2/kk
                ENDIF
            ENDDO
        ENDDO
    ENDDO


    ! Oddball
    ALLOCATE(firstslice(Uk%ny,Uk%nz),stat=ierr)
    IF(ierr.NE.0) THEN
       WRITE(6,'(a,i0)')'[ERROR] hitvelocity_init Failed to allocate firstslice on ', spec_rank
       RETURN
    ENDIF
    firstslice=extractYZSlice(Uk,spec_rank,1)
    IF (Uk%xmin.LE.1 .AND. Uk%xmax.GE.1) THEN
      DO k=max(2,Uk%zmin),Uk%zmax
         DO j=max(nk+1,Uk%ymin),Uk%ymax
            Uk%values(1,j,k)=conjg(firstslice(Uk%ny+2-j,Uk%nz+2-k))
         ENDDO
      ENDDO
    ENDIF
    firstslice=extractYZSlice(Vk,spec_rank,1)
    IF (Vk%xmin.LE.1 .AND. Vk%xmax.GE.1) THEN
      DO k=max(2,Vk%zmin),Vk%zmax
         DO j=max(nk+1,Uk%ymin),Vk%ymax
            Vk%values(1,j,k)=conjg(firstslice(Vk%ny+2-j,Vk%nz+2-k))
         ENDDO
      ENDDO
    ENDIF
    firstslice=extractYZSlice(Wk,spec_rank,1)
    IF (Wk%xmin.LE.1 .AND. Wk%xmax.GE.1) THEN
      DO k=max(2,Wk%zmin),Wk%zmax
         DO j=max(nk+1,Wk%ymin),Wk%ymax
            Wk%values(1,j,k)=conjg(firstslice(Wk%ny+2-j,Wk%nz+2-k))
         ENDDO
      ENDDO
    ENDIF
    DEALLOCATE(firstslice)

!!$!!!! OLD ODDBALL NOT PROPERLY HANDLED - NEED TO BE FIXED
!!$!!!! FOR NOW - ONLY THE SIMPLE BOUNDARY (1,1,K) HANDLED
!    DO k=Uk%zmin,Uk%zmax
!         Uk%values(1,1,k)=conjg(Uk%values(1,1,Uk%nz+2-k))
!         Vk%values(1,1,k)=conjg(Vk%values(1,1,Uk%nz+2-k))
!         Wk%values(1,1,k)=conjg(Wk%values(1,1,Uk%nz+2-k))
!    ENDDO

    
    ALLOCATE(sg(U%nx+1),STAT=ierr) 
    IF(ierr.NE.0) THEN
      WRITE(6,'()')'[ERROR] hitvelocity_init: Not enought memory for '//&
      &            'sg array on processor ',spec_rank
      RETURN
    ENDIF
    sg = 0.0_WP

    ALLOCATE(node_index(Uk%xmin:Uk%xmax,Uk%ymin:Uk%ymax,Uk%zmin:Uk%zmax),STAT=ierr) 
    IF(ierr.NE.0) THEN
      WRITE(6,'()')'[ERROR] hitvelocity_init: Not enought memory for '//&
      &            'node_index array on processor ',spec_rank
      RETURN
    ENDIF
    
    DO k = Uk%zmin,Uk%zmax
        DO j = Uk%ymin,Uk%ymax
            DO i = Uk%xmin,Uk%xmax
                kk=sqrt(wavenumb%kx(i)*wavenumb%kx(i)+wavenumb%ky(j)*wavenumb%ky(j)+wavenumb%kz(k)*wavenumb%kz(k))
                ik = 1+int(kk/dk + 0.5_WP)
                node_index(i,j,k) = ik
                IF (kk.GT.eps.AND.kk.LE.kc) THEN
                    sg(ik) = sg(ik) &
                        & + REAL(Uk%values(i,j,k)*conjg(Uk%values(i,j,k))) &
                        & + REAL(Vk%values(i,j,k)*conjg(Vk%values(i,j,k))) &
                        & + REAL(Wk%values(i,j,k)*conjg(Wk%values(i,j,k)))
                ENDIF
            ENDDO
        ENDDO
    ENDDO

    sg=DoGlobalVectorSum(spec_rank,sg)

    DO i = 1,SIZE(sg)
        IF (sg(i).gt.1e-20_WP) THEN
            sg(i) = spect(i,2)/sg(i)
        ELSE
            sg(i) = 1.0_WP
        ENDIF
    ENDDO

    DO k = Uk%zmin,Uk%zmax
        DO j = Uk%ymin,Uk%ymax
            DO i = Uk%xmin,Uk%xmax
                racsg=sqrt(sg(node_index(i,j,k)))
                Uk%values(i,j,k) = Uk%values(i,j,k)*racsg
                Vk%values(i,j,k) = Vk%values(i,j,k)*racsg
                Wk%values(i,j,k) = Wk%values(i,j,k)*racsg
            ENDDO
        ENDDO
    ENDDO

    CALL btran(UK,U,res1)
    CALL btran(VK,V,res2)
    CALL btran(WK,W,res3)
    hitvelocity_init=(res1.AND.res2.AND.res3)
  
    !--------------------------------------------------
    ! release allocated memory for local variables
    !--------------------------------------------------
    DEALLOCATE(spect)
    DEALLOCATE(sg)
    DEALLOCATE(node_index)
    CALL deleteDataLayout(Uk)  
    CALL deleteDataLayout(Vk)  
    CALL deleteDataLayout(Wk)  
    CALL deleteDataLayout(ak)  
    CALL deleteDataLayout(bk)

    hitvelocity_init=(res1.AND.res2.AND.res3)
      
    RETURN
    
END FUNCTION hitvelocity_init




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
!> @param [in] ke location of the peak opf the spectrum
!> @param [in] spec_amp amplitude of the spectrum
!> @return .TRUE. if initialization is successfull. 
!------------------------------------------------------------------------------
!=========================================================================================!
FUNCTION hit_vector_init(Vec1, Vec2, Vec3, wavenumb, Lx, nbcpus, spec_rank,ke,spec_amp)
!=========================================================================================!

    IMPLICIT NONE
    TYPE(REAL_DATA_LAYOUT), INTENT(INOUT) :: Vec1
    TYPE(REAL_DATA_LAYOUT), INTENT(INOUT) :: Vec2
    TYPE(REAL_DATA_LAYOUT), INTENT(INOUT) :: Vec3
    TYPE(WaveNumbers), INTENT(IN)::wavenumb
    REAL(WP), INTENT(IN) :: Lx, ke, spec_amp
    INTEGER, INTENT(IN)  :: nbcpus, spec_rank
    LOGICAL :: hit_vector_init

    TYPE(COMPLEX_DATA_LAYOUT) :: Vec1k
    TYPE(COMPLEX_DATA_LAYOUT) :: Vec2k
    TYPE(COMPLEX_DATA_LAYOUT) :: Vec3k
    TYPE(COMPLEX_DATA_LAYOUT) :: ak
    TYPE(COMPLEX_DATA_LAYOUT) :: bk

    REAL(WP) :: pi
    REAL(WP) :: dk, kc, eps
    REAL(WP) :: energy_spec
    REAL(WP) :: kk, kk2
    REAL(WP) :: amp_disc
    REAL(WP), DIMENSION(:,:), ALLOCATABLE :: spect
    REAL(WP), DIMENSION(:), ALLOCATABLE :: sg
    REAL(WP) :: ps1, ps2, psr, racsg
    INTEGER  :: nk, ik
    INTEGER  :: i,j,k,ierr
    INTEGER, DIMENSION(:,:,:),ALLOCATABLE :: node_index 
    COMPLEX(WP), PARAMETER :: ii=(0.0_WP,1.0_WP)
    COMPLEX(WP), DIMENSION(:,:),ALLOCATABLE :: firstslice
    LOGICAL :: res1,res2,res3

    hit_vector_init=.FALSE.

    !--------------------------------------------------
    ! Create pi number
    !--------------------------------------------------
    pi=acos(-1.0_WP)
    nk=Vec1%nx/2 + 1     ! a changer si nx != ny != nz
    dk=2.0_WP*pi/Lx
    amp_disc = sqrt(dk)**3
    kc=pi*REAL(Vec1%nx,WP)/Lx
    eps=kc/10000000_WP

    !--------------------------------------------------
    ! Instantiation of the arrays in the spectral domain
    !--------------------------------------------------
    IF (.NOT. initDataLayout("Vec1_Spectral",Vec1k,(Vec1%nx/2)+1,Vec1%ny,Vec1%nz,Lx,Lx,Lx,nbcpus,spec_rank,alongZ)) THEN
       WRITE(6,'(a,i0)')'[ERROR] hit_vector_init Failed to init Vec1k on ', spec_rank
       RETURN
    ENDIF
    IF (.NOT. initDataLayout("Vec2_Spectral",Vec2k,(Vec2%nx/2)+1,Vec2%ny,Vec2%nz,Lx,Lx,Lx,nbcpus,spec_rank,alongZ)) THEN
       WRITE(6,'(a,i0)')'[ERROR] hit_vector_init Failed to init Vec2k on ', spec_rank
       RETURN
    ENDIF
    IF (.NOT. initDataLayout("Vec3_Spectral",Vec3k,(Vec3%nx/2)+1,Vec3%ny,Vec3%nz,Lx,Lx,Lx,nbcpus,spec_rank,alongZ)) THEN
       WRITE(6,'(a,i0)')'[ERROR] hit_vector_init Failed to init Vec3k on ', spec_rank
       RETURN
    ENDIF
    IF (.NOT. initDataLayout("Ak_Spectral",ak,(Vec1%nx/2)+1,Vec1%ny,Vec1%nz,Lx,Lx,Lx,nbcpus,spec_rank,alongZ)) THEN
       WRITE(6,'(a,i0)')'[ERROR] hit_vector_init Failed to init ak on ', spec_rank
       RETURN
    ENDIF
    IF (.NOT. initDataLayout("Bk_Spectral",bk,(Vec1%nx/2)+1,Vec1%ny,Vec1%nz,Lx,Lx,Lx,nbcpus,spec_rank,alongZ)) THEN
       WRITE(6,'(a,i0)')'[ERROR] hit_vector_init Failed to init bk on ', spec_rank
       RETURN
    ENDIF

    ! Output spectrum for comparison
    ALLOCATE(spect(Vec1%nx+1,2),STAT=ierr) !WARNING spect array shape is different
    IF(ierr.NE.0) THEN
        WRITE(6,'()')'[ERROR] hit_vector_init: Not enought memory for spect on processor ',spec_rank
        RETURN
    ENDIF
    spect(:,1) = (/(REAL(i-1,WP)*dk,i=1,size(spect,1))/)
    spect(:,2) = 0.0_WP

    ! Compute spectrum
    DO k=ak%zmin,ak%zmax
        DO j=ak%ymin,ak%ymax
            DO i=ak%xmin,ak%xmax
                kk=sqrt(wavenumb%kx(i)*wavenumb%kx(i)+wavenumb%ky(j)*wavenumb%ky(j)+wavenumb%kz(k)*wavenumb%kz(k))

                ! Random numbers
                CALL random_number(psr)
                psr= 2.0_WP*pi*psr                !2.0_WP*pi*(rand-0.5_WP)
                CALL random_number(ps1)
                ps1=2.0_WP*pi*ps1                 !*(rand-0.5_WP)
                CALL random_number(ps2)
                ps2=2.0_WP*pi*ps2                 !(rand-0.5_WP)

                ! Spectrums
                energy_spec=spec_amp*(kk/ke)**4*exp(-2.0_WP*(kk/ke)**2)

                ! Coeff
                ik = 1+idint(kk/dk + 0.5_WP)
                IF ((kk.gt.eps).and.(kk.le.kc)) THEN
                    ak%values(i,j,k)= amp_disc * sqrt(energy_spec/(4.0_WP*pi*kk**2)) &
                        &                      * exp(ii*ps1) * cos(psr)
                    bk%values(i,j,k)= amp_disc * sqrt(energy_spec/(4.0_WP*pi*kk**2)) &
                        &                      * exp(ii*ps2) * sin(psr)
                    spect(ik,2) = spect(ik,2) + REAL(ak%values(i,j,k)*conjg(ak%values(i,j,k))) &
                        &                     + REAL(bk%values(i,j,k)*conjg(bk%values(i,j,k)))
                ENDIF
            ENDDO
        ENDDO
    ENDDO
    spect(:,2)=DoGlobalVectorSum(spec_rank,spect(:,2))

    DO k = Vec1k%zmin,Vec1k%zmax
        DO j = Vec1k%ymin,Vec1k%ymax
            DO i = Vec1k%xmin,Vec1k%xmax
                kk=sqrt(wavenumb%kx(i)*wavenumb%kx(i)+wavenumb%ky(j)*wavenumb%ky(j)+wavenumb%kz(k)*wavenumb%kz(k))
                IF ((kk.GT.eps).and.(kk.LE.kc)) THEN
                    kk2=sqrt(wavenumb%kx(i)*wavenumb%kx(i)+wavenumb%ky(j)*wavenumb%ky(j))
                    if (kk2.lt.eps) then
                        Vec1k%values(i,j,k)=(ak%values(i,j,k)+bk%values(i,j,k))/sqrt(2.0_WP)
                        Vec2k%values(i,j,k)=(bk%values(i,j,k)-ak%values(i,j,k))/sqrt(2.0_WP)
                    else
                        Vec1k%values(i,j,k)=(ak%values(i,j,k)*kk*wavenumb%ky(j)+bk%values(i,j,k)*wavenumb%kx(i) &
                            &            *wavenumb%kz(k))/(kk*kk2)
                        Vec2k%values(i,j,k)=(bk%values(i,j,k)*wavenumb%ky(j)*wavenumb%kz(k)-ak%values(i,j,k) &
                            &            *kk*wavenumb%kx(i))/(kk*kk2)
                    end if
                    Vec3k%values(i,j,k)=-bk%values(i,j,k)*kk2/kk
                ENDIF
            ENDDO
        ENDDO
    ENDDO


    ! Oddball
    ALLOCATE(firstslice(Vec1k%ny,Vec1k%nz),stat=ierr)
    IF(ierr.NE.0) THEN
       WRITE(6,'(a,i0)')'[ERROR] hit_vector_init Failed to allocate firstslice on ', spec_rank
       RETURN
    ENDIF
    firstslice=extractYZSlice(Vec1k,spec_rank,1)
    IF (Vec1k%xmin.LE.1 .AND. Vec1k%xmax.GE.1) THEN
      DO k=max(2,Vec1k%zmin),Vec1k%zmax
         DO j=max(nk+1,Vec1k%ymin),Vec1k%ymax
            Vec1k%values(1,j,k)=conjg(firstslice(Vec1k%ny+2-j,Vec1k%nz+2-k))
         ENDDO
      ENDDO
    ENDIF
    firstslice=extractYZSlice(Vec2k,spec_rank,1)
    IF (Vec2k%xmin.LE.1 .AND. Vec2k%xmax.GE.1) THEN
      DO k=max(2,Vec2k%zmin),Vec2k%zmax
         DO j=max(nk+1,Vec1k%ymin),Vec2k%ymax
            Vec2k%values(1,j,k)=conjg(firstslice(Vec2k%ny+2-j,Vec2k%nz+2-k))
         ENDDO
      ENDDO
    ENDIF
    firstslice=extractYZSlice(Vec3k,spec_rank,1)
    IF (Vec3k%xmin.LE.1 .AND. Vec3k%xmax.GE.1) THEN
      DO k=max(2,Vec3k%zmin),Vec3k%zmax
         DO j=max(nk+1,Vec3k%ymin),Vec3k%ymax
            Vec3k%values(1,j,k)=conjg(firstslice(Vec3k%ny+2-j,Vec3k%nz+2-k))
         ENDDO
      ENDDO
    ENDIF
    DEALLOCATE(firstslice)

!!$!!!! OLD ODDBALL NOT PROPERLY HANDLED - NEED TO BE FIXED
!!$!!!! FOR NOW - ONLY THE SIMPLE BOUNDARY (1,1,K) HANDLED
!    DO k=Uk%zmin,Uk%zmax
!         Uk%values(1,1,k)=conjg(Uk%values(1,1,Uk%nz+2-k))
!         Vk%values(1,1,k)=conjg(Vk%values(1,1,Uk%nz+2-k))
!         Wk%values(1,1,k)=conjg(Wk%values(1,1,Uk%nz+2-k))
!    ENDDO


    ALLOCATE(sg(Vec1%nx+1),STAT=ierr) 
    IF(ierr.NE.0) THEN
      WRITE(6,'()')'[ERROR] hit_vector_init: Not enought memory for '//&
      &            'sg array on processor ',spec_rank
      RETURN
    ENDIF
    sg = 0.0_WP

    ALLOCATE(node_index(Vec1k%xmin:Vec1k%xmax,Vec1k%ymin:Vec1k%ymax,Vec1k%zmin:Vec1k%zmax),STAT=ierr) 
    IF(ierr.NE.0) THEN
      WRITE(6,'()')'[ERROR] hit_vector_init: Not enought memory for '//&
      &            'node_index array on processor ',spec_rank
      RETURN
    ENDIF

    DO k = Vec1k%zmin,Vec1k%zmax
        DO j = Vec1k%ymin,Vec1k%ymax
            DO i = Vec1k%xmin,Vec1k%xmax
                kk=sqrt(wavenumb%kx(i)*wavenumb%kx(i)+wavenumb%ky(j)*wavenumb%ky(j)+wavenumb%kz(k)*wavenumb%kz(k))
                ik = 1+int(kk/dk + 0.5_WP)
                node_index(i,j,k) = ik
                IF (kk.GT.eps.AND.kk.LE.kc) THEN
                    sg(ik) = sg(ik) &
                        & + REAL(Vec1k%values(i,j,k)*conjg(Vec1k%values(i,j,k))) &
                        & + REAL(Vec2k%values(i,j,k)*conjg(Vec2k%values(i,j,k))) &
                        & + REAL(Vec3k%values(i,j,k)*conjg(Vec3k%values(i,j,k)))
                ENDIF
            ENDDO
        ENDDO
    ENDDO

    sg=DoGlobalVectorSum(spec_rank,sg)

    DO i = 1,SIZE(sg)
        IF (sg(i).gt.1e-20_WP) THEN
            sg(i) = spect(i,2)/sg(i)
        ELSE
            sg(i) = 1.0_WP
        ENDIF
    ENDDO

    DO k = Vec1k%zmin,Vec1k%zmax
        DO j = Vec1k%ymin,Vec1k%ymax
            DO i = Vec1k%xmin,Vec1k%xmax
                racsg=sqrt(sg(node_index(i,j,k)))
                Vec1k%values(i,j,k) = Vec1k%values(i,j,k)*racsg
                Vec2k%values(i,j,k) = Vec2k%values(i,j,k)*racsg
                Vec3k%values(i,j,k) = Vec3k%values(i,j,k)*racsg
            ENDDO
        ENDDO
    ENDDO

    CALL btran(Vec1K,Vec1,res1)
    CALL btran(Vec2K,Vec2,res2)
    CALL btran(Vec3K,Vec3,res3)
    hit_vector_init=(res1.AND.res2.AND.res3)

    !--------------------------------------------------
    ! release allocated memory for local variables
    !--------------------------------------------------
    DEALLOCATE(spect)
    DEALLOCATE(sg)
    DEALLOCATE(node_index)
    CALL deleteDataLayout(Vec1k)  
    CALL deleteDataLayout(Vec2k)  
    CALL deleteDataLayout(Vec3k)  
    CALL deleteDataLayout(ak)  
    CALL deleteDataLayout(bk)

    hit_vector_init=(res1.AND.res2.AND.res3)

    RETURN

END FUNCTION hit_vector_init


!------------------------------------------------------------------------------
!> @author 
!> Guillaume Balarac, Patrick BEGOU, LEGI
!
!>
!> @details
!> This function initialise a scalar for starting
!> an homogeneous isotropic turbulence simulation.
!
!> @param [in,out] physcal the scalar storage to initialize (in physical space)
!> @param [in] wavenumb the waves numbers for this scalar
!> @param [in] Lx the length of the domain. It is supposed to be the same in the
!> three dimension in this version (Lx==Ly==Lz).
!> @param [in] nbcpus the total amount of cpus used by the application
!> @param [in] spec_rank my MPI processe rank in the range [0:ncpus-1]
!> @return .TRUE. if initialization is successfull. 
!------------------------------------------------------------------------------
!=========================================================================================!
  FUNCTION hitscalar_init(physcal, wavenumb, lx, nbcpus, spec_rank)
!=========================================================================================!
    IMPLICIT NONE
    TYPE(REAL_DATA_LAYOUT), INTENT(INOUT) :: physcal
    TYPE(WaveNumbers), INTENT(IN)::wavenumb
    INTEGER, INTENT(IN)  :: nbcpus, spec_rank
    REAL(WP), INTENT(IN) :: Lx
    LOGICAL :: hitscalar_init

    TYPE(COMPLEX_DATA_LAYOUT) :: fouscal
    REAL(WP) :: pi
    REAL(WP) :: dk
    REAL(WP) :: ksk0ratio, ks, kc
!   REAL(WP) :: kx,ky,kz,kk, kk2, rand,kcut
    REAL(WP) :: kk, rand,kcut
    COMPLEX(WP), PARAMETER :: ii=(0.0_WP,1.0_WP)
    INTEGER  :: nk
    INTEGER  :: i,j,k
    LOGICAL :: res
    
    hitscalar_init=.FALSE.

    IF (.NOT. initDataLayout("Scalar_Spectral",fouscal,(physcal%nx/2)+1,physcal%ny,physcal%nz, &
                    &         lx,lx,lx,nbcpus,spec_rank,alongZ)) THEN
       WRITE(6,'(a,i0)')'[ERROR] hitscalar_init Failed to init fouscal on ', spec_rank
       RETURN
    ENDIF
    
    !--------------------------------------------------
    ! Create pi number
    !--------------------------------------------------
    pi=acos(-1.0_WP)
    nk=physcal%nx/2 + 1     ! a changer si nx != ny != nz
    dk=2.0_WP*pi/Lx
    CALL parser_read('ks/ko', ksk0ratio)
    ks = ksk0ratio*dk
    kc = 2.0_WP*ks

    DO k = fouscal%zmin,fouscal%zmax
        DO j = fouscal%ymin,fouscal%ymax
            DO i = fouscal%xmin,fouscal%xmax

                kk=sqrt(wavenumb%kx(i)*wavenumb%kx(i)+wavenumb%ky(j)*wavenumb%ky(j)+wavenumb%kz(k)*wavenumb%kz(k))

                IF (kk.lt.1e-10) THEN
                    fouscal%values(i,j,k) = 0.0_WP
                ELSE
                    IF ((ks-dk/2.0_WP.le.kk).and.(kk.le.ks+dk/2.0_WP)) THEN
                        !f_phi = 1.0_WP
                        CALL random_number(rand)
                        fouscal%values(i,j,k) = &
                            & sqrt(1.0_WP/(4.0_WP*pi*kk**2))*exp(ii*2.0_WP*pi*rand)
                    ELSE
                        !f_phi = 0.0_WP
                        fouscal%values(i,j,k) = 0.0_WP
                    ENDIF
                ENDIF
            ENDDO
        ENDDO
    ENDDO
    
    CALL btran(fouscal,physcal,res)
    IF(.NOT. res) THEN
        WRITE(6,'(a,i0)')'[ERROR] hitscalar_init Failed in fouscal backsward FFT on ', spec_rank
        RETURN
    ENDIF
    
    WHERE(physcal%values <= 0.0_WP)
        physcal%values = 0.0_WP
    ELSEWHERE
        physcal%values = 1.0_WP
    ENDWHERE

    CALL ftran(physcal,fouscal,res)
    IF(.NOT. res) THEN
        WRITE(6,'(a,i0)')'[ERROR] hitscalar_init Failed in physcal forward FFT on ', spec_rank
        RETURN
    ENDIF

    kcut = REAL(physcal%nx,WP)/2.0_WP  ! Attention si nx != ny != nz
    DO k = fouscal%zmin,fouscal%zmax
        DO j = fouscal%ymin,fouscal%ymax
            DO i = fouscal%xmin,fouscal%xmax
                kk=sqrt(wavenumb%kx(i)*wavenumb%kx(i)+wavenumb%ky(j)*wavenumb%ky(j)+wavenumb%kz(k)*wavenumb%kz(k))
                IF (kk.gt.kcut) THEN
                    fouscal%values(i,j,k) = 0.0_WP
                ELSE IF (kk.GT.kc) then
                    fouscal%values(i,j,k) = fouscal%values(i,j,k)*(kc/kk)**2
                ENDIF
            ENDDO
        ENDDO
    ENDDO
    CALL btran(fouscal,physcal,res)
    IF(.NOT. res) THEN
        WRITE(6,'(a,i0)')'[ERROR] hitscalar_init Failed in fouscal second backsward FFT on ', spec_rank
        RETURN
    ENDIF
     
    IF (spec_rank .EQ.0) WRITE(6,'(a)')'[WARNING] hitscalar_init call ' //&
                         &           'to get_dns_numbers not implemented'
    !CALL get_dns_numbers(1)

    
    CALL deleteDataLayout(fouscal)  
    
    hitscalar_init=.TRUE.

    
  END FUNCTION hitscalar_init

END MODULE isotropicturbulence
