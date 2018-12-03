!> @addtogroup post_process
!! @{

!------------------------------------------------------------------------------

!
! MODULE: post_in_simu
!
!
! DESCRIPTION:
!>  This module provide a library of post-process routines. This library
!! provide all that you need to build your own post-process.
!
!> @details
!!  This module provide all post-process implement in scales. It has to
!! be used to build your own post-process traitement. Note that for post-
!! performed during computation, you have to used the module
!! post_in_simu routine (another exe will be add later to perform post-
!! process separately from simulation).
!
!> @author
!! Jean-Baptiste Lagaert, LEGI
!
!------------------------------------------------------------------------------

module post_in_simu

    use precision_tools
    use datalayout
    use post_lib
    use post_scalar
    use post_lib2D
    use avs
    use vtkxml_plus
    use post_mhd
    use post_kolmo
    use system_tools
    use stat_tools
    use parser_tools
    use post_hd
    implicit none

    private

    integer                          :: show_basic   = 1
    integer                          :: compute_spec = 1
    integer                          :: mhd_post = 1
    integer                          :: hd_post = 1
    integer                          :: iPostKolmo
    REAL(WP)                         :: mu
    REAL(WP)                         :: Prm
    !VAR LUCA
    REAL(WP)                         :: avgXY,corrXY
    integer                          :: ifreq
    !> Iteration frequency for plane jet post-processing
    integer                          :: jet_freq
    integer                          :: save_freq       ! frequency for saving field
    integer                          :: ishowdiv        ! frequency of display of velocity divergence
    integer                          :: show_scalar     ! frequency of computation of scalar PDF
    integer                          :: scal_basic_stat ! frequency of computation of scalar momentum (skew, etc)
    integer                          :: scal_moment_stat ! frequency of computation of scalar momentum (skew, etc)
    !VAR LUCA
    INTEGER                          :: freq_sigma_avgXY ! frequency of sigma average in XY plane
    INTEGER                          :: freq_sigma_corr  ! frequency of sigma correlation in XY plane
    INTEGER                          :: Nw               ! Waiting time

    integer                          :: show_mem
    character(len=16)                :: output_type     ! type of output (avs or paraview)
    integer,dimension(:),allocatable :: tagScal
    integer,dimension(:),allocatable :: tagScalPart

    logical                          :: Em_Growth_rate  ! for display or not the magnetic growth rate
! XXX Todo
!    integer,dimension(:),allocatable :: tagScalPartVector
! XXX End Todo
    integer,dimension(3)             :: tagVel
    integer,dimension(3)             :: tagB
    real(WP), dimension(:), allocatable, save :: z_pos

    ! -- Time frequency --
    !> Initial time
    real(WP)                        :: init_time
    !> Time frequency for spectrum computation
    real(WP)                        :: spec_deltaT
    !> LUCA time frequency for sigma XY average
    real(WP)                        :: sigma_avg_XY_deltaT
    !> LUCA time frequency for sigma average XY 
    integer                         :: ite_sigma_avg_XY
    !> Number of spectrum output related to real time and not iteration number
    integer                         :: ite_spec
    !> Time frequency for plane jet post-processing
    real(WP)                        :: jet_deltaT
    !> Number of jet post-processing  output related to real time and not iteration number
    integer                         :: ite_jet
    !> Time frequency for visualisation output
    real(WP)                        :: visu_deltaT
    !> Number of visualisation output related to real time and not iteration number
    integer                         :: ite_visu
    !> Time frequency for post_process
    real(WP)                        :: post_deltaT
    !> Type of paraview output
    logical                         :: paraview_3D, paraview_slice

    ! ==== Public procedures ====
    public                           :: post_simu,post_simu_luca
    public                           :: post_simu_deltaT,post_simu_deltaT_luca
    public                           :: post_in_simu_init,post_in_simu_init_luca
    public                           :: post_in_simu_del

    !> tools to compute differential diffusion
    integer                         :: tag_diff_diff
    type(REAL_DATA_LAYOUT),save     :: diff_diff
    integer                         :: ite_diff_diff, diff_diff_freq
    integer                         :: scal_1, scal_2
    real(WP)                        :: mix_th, diff_diff_deltaT
    logical                         :: differential_diff

    !variables for a priori test in simu
    integer                         :: apriori_sca
    integer                         :: apriori_velo_freq 
    integer                         :: nbAPTest,iAPTest
    integer                         :: APFiltMin,APFiltMax
    integer                         :: APFilter,ires
    character(len=str_medium),dimension(:),allocatable  :: APTest
    character(len=str_medium)       :: request 


    !> To dump Q criteria and Vorticity for paraview format
    integer                         :: tag_Qc
    integer,dimension(3)            :: tag_O
    type(REAL_DATA_LAYOUT),save     :: Qc, O1, O2, O3
    logical                         :: write_Qc_omega

contains


!------------------------------------------------------------------------------
!> Perform advanced post-process setup.
!! @author Jean-Baptiste Lagaert and Antoine Vollant
!!    @param[in]    U               = velocity along X-axis
!!    @param[in]    V               = velocity along Y-axis
!!    @param[in]    W               = velocity along Z-axis
!!    @param[in]    B               = magnetic field (3D vector)
!!    @param[in]    scalarray       = arrays of scalars in physical space
!!    @param[in]    scal_partArray  = arrays of scalars in physical space
!!    @param[in]    imhd            = 1 if MHD is ON, 0 if it is OFF
!!    @param[in]    nbscal          = the number of scalars to process
!!    @param[in]    nbscal_part     = the number of scalars to process
!!    @param[in]    me              = mpi rank
!!    @param[out]   post_freq       = frequence off post-processing
!!    @return       success         = logical equals to true is everything is right.
!------------------------------------------------------------------------------
function post_in_simu_init(U,V,W,B,scalArray,scal_partArray,imhd,nbscal,nbscal_part,me, post_freq) result(success)

    use param
    use data
    use cart_topology
    implicit none

    !I/O data
    logical,intent(in)                                          :: imhd
    integer,intent(in)                                          :: nbscal
    integer,intent(in)                                          :: nbscal_part
    integer,intent(in)                                          :: me
    type(real_data_layout),intent(in)                           :: U, V ,W
    type(real_data_layout), intent(in),pointer, dimension(:)    :: B
    type(real_data_layout), intent(in),pointer, dimension(:)    :: scalArray
    type(real_data_layout), intent(in),pointer, dimension(:)    :: scal_partArray
! XXX Todo
!    type(real_vector_data_layout),intent(in)                    :: scal_partVector
! XXX End Todo
    integer, intent(out)                                        :: post_freq
    logical :: success

    !Local data
    integer                          :: nbField,i, k, nb_z_pos
    integer                          :: spec_rank
    real(WP),dimension(3)            :: length_3D
    character (len=10)               :: type_para       ! For choosing if 3D, slice or both
    character (len=30)               :: parser_tag
! XXX Todo
!    type(real_data_layout)           :: ScalPart
! XXX End Todo

    success = .false.

    length_3D(1)=U%Lx
    length_3D(2)=U%Ly
    length_3D(3)=U%Lz

    ! Init output frequency related to iteration number
    CALL parser_read('Simulation iterations',post_freq)
    CALL parser_read('Show basic info',show_basic)
    if (show_basic>0) post_freq = min(show_basic, post_freq)
    CALL parser_read('Show Magnetic growth rate',Em_Growth_rate)
    CALL parser_read('Show scalar',show_scalar)
    if (show_scalar>0) post_freq = min(show_scalar, post_freq)
    CALL parser_read('Scalar - basic stats',scal_basic_stat)
    if (scal_basic_stat>0) post_freq = min(scal_basic_stat, post_freq)
    CALL parser_read('Scalar - moment stats',scal_moment_stat)
    if (scal_moment_stat>0) post_freq = min(scal_moment_stat, post_freq)
    CALL parser_read('Compute spectrum - ite',compute_spec)
    if (compute_spec>0) post_freq = min(compute_spec, post_freq)
    CALL parser_read('MHD post process',mhd_post)
    if (mhd_post>0) post_freq = min(mhd_post, post_freq)
    CALL parser_read('HD post process',hd_post)
    if (hd_post>0) post_freq = min(hd_post, post_freq)
    CALL parser_read('Viscosity',mu)
    CALL parser_read('Magnetic Prandl number',Prm)
    CALL parser_read('Write frequency',ifreq)
    if (ifreq>0) post_freq = min(ifreq, post_freq)
    CALL parser_read('Save frequency',save_freq)
    if (save_freq>0) post_freq = min(save_freq, post_freq)
    CALL parser_read('Show divergence',ishowdiv)
    if (ishowdiv>0) post_freq = min(ishowdiv, post_freq)
    CALL parser_read('Post Kolmo',iPostKolmo)
    if (iPostKolmo>0) post_freq = min(iPostKolmo, post_freq)

    CALL parser_read('Type of output',output_type)

    if (mhd_post>0) then
      CALL post_mhd_init(mu,Prm,Uk)
    endif
    if (hd_post>0) then
      CALL post_hd_init(mu,Uk)
    endif

    if (iPostkolmo>0) then
      CALL post_kolmo_init(U,imhd)
    endif

    call parser_read('Show memory',show_mem)
    if (show_mem>0) post_freq = min(show_mem, post_freq)
    ! Some Jet post-process
    CALL parser_read('Plane jet AVG - ite',jet_freq)
    if (jet_freq>0) post_freq = min(jet_freq, post_freq)
    call parser_read('diff diff - ite', diff_diff_freq)
    if (diff_diff_freq>0) then
        post_freq = min(diff_diff_freq, post_freq)
        call parser_read('mixing threshold', mix_th)
        call parser_read('diff diff scal1', scal_1)
        call parser_read('diff diff scal2', scal_2)
        differential_diff = .true.
    end if

    ! Init output frequency related to iteration number
    post_deltaT = 1000000._WP
    init_time = 0.0_WP
    if (parser_is_defined('Initial time')) call parser_read('Initial time', init_time)
    !CALL parser_read('Final time',post_deltaT)
!    CALL parser_read('Show basic info',show_basic)
!    if (show_basic>0) post_freq = min(show_basic, post_freq)
!    CALL parser_read('Compute scalar PDF',show_scalar)
!    if (show_scalar>0) post_freq = min(show_scalar, post_freq)
    CALL parser_read('Compute spectrum - time',spec_deltaT)
    if (spec_deltaT>0) then
        post_deltaT = min(spec_deltaT, post_deltaT)
        ite_spec = floor(init_time/spec_deltaT)
    end if
    CALL parser_read('Write output - time',visu_deltaT)
    if (visu_deltaT>0) then
        if (ifreq>0) then
            write(*,'(a,a)') '[WARNING] Ouptut are activated with both iterations', &
                & 'and time frequency - Only iterations will be considered'
            visu_deltaT = -1.0_WP
        else
            post_deltaT = min(visu_deltaT, post_deltaT)
            ite_visu = floor(init_time/visu_deltaT) - 1
        end if
    end if
    CALL parser_read('Plane jet AVG - time',jet_deltaT)
    if (jet_deltaT>0) then
        post_deltaT = min(jet_deltaT, post_deltaT)
        ite_jet = floor(init_time/jet_deltaT)
    end if

    ! differential diffusion for jet
    call parser_read('diff diff - time', diff_diff_deltaT)
    if (diff_diff_deltaT>0) then
        post_deltaT = min(diff_diff_deltaT, post_deltaT)
        call parser_read('mixing threshold', mix_th)
        call parser_read('diff diff scal1', scal_1)
        call parser_read('diff diff scal2', scal_2)
        differential_diff = .true.
    end if

    ! apriori test in simu
    call parser_read('a priori test in simu for velocity', apriori_velo_freq)
    if (apriori_sca>0) then
        post_freq = min(apriori_velo_freq, post_freq)
    endif
    call parser_read('a priori test in simu', apriori_sca)
    if (apriori_sca>0) then
        post_freq = min(apriori_sca, post_freq)
        spec_rank=getnbmycpu()
        call parser_read('number of a priori test', nbAPTest)
        if ( nbAPTest .lt. 1 ) then
          write (6,'(a)') '[ERROR] number of test is not right'
          return
        endif
        call parser_read('filter size min for a priori test',APfiltMin)
        if ( APfiltMin .lt. 2 ) then
          write (6,'(a)') '[ERROR] filter size min not correct'
          return
        endif
        call parser_read('filter size max for a priori test',APfiltMax)
        if ( APfiltMax .lt. APfiltMin ) then
          write (6,'(a)') '[ERROR] filter size max not correct'
          return
        endif
        write(request,'(a)') 'type of filter for a priori test'
        if (.not. parseFilter(APfilter,spec_rank,other=request) ) then
          write (6,'(a)') '[ERROR] reading filter in post_in_simu...'
          return
        endif
        allocate(APTest(nbAPTest),stat=ires)
        if ( ires .ne. 0 ) write (6,'(a)') '[ERROR] not enought memory to allocate APTest'
        do iAPTest= 1,nbAPTest
          write(request,'(a,1x,i0)') 'name of a priori test',iAPTest
          if (parser_is_defined(request)) then
            call parser_read(trim(adjustl(request)), APTest(iAPTest))
          else
            write (6,'(a)') '[ERROR] Name of one of a priori test is not set' 
            return
          endif
        enddo
    end if

!    CALL parser_read('MHD post process',mhd_post)
!    if (mhd_post>0) post_freq = min(mhd_post, post_freq)
!    CALL parser_read('Viscosity',mu)
!    CALL parser_read('Magnetic Prandl number',Prm)
!    CALL parser_read('Show divergence',ishowdiv)
!    if (ishowdiv>0) post_freq = min(ishowdiv, post_freq)
!    CALL parser_read('Type of output',output_type)
!    if (mhd_post>0) then
!      CALL post_mhd_init(mu,Prm,Uk)
!    endif
!    call parser_read('Show memory',show_mem)
!    if (show_mem>0) post_freq = min(show_mem, post_freq)

    ! Compute vorticity and Q criterium (the old flag is still working until no
    ! merge between paraview write and generic write)
    if(parser_is_defined('Compute Q criterium')) call parser_read('Compute Q criterium', write_Qc_Omega)

    !Init paraview context
    if (output_type .eq. 'paraview' )then
      call parser_read('slice or 3D', type_para)
      select case(trim(type_para))
      case('slice')
          paraview_slice = .true.
          paraview_3D   = .false.
      case('3D')
          paraview_3D = .true.
          paraview_slice = .false.
      case('3Dplus')
          paraview_3D = .true.
          write_Qc_omega = .true.
          paraview_slice = .false.
      case default
          paraview_slice = .true.
          paraview_3D = .true.
      end select



      if (imhd) then
        nbField = nbscal+nbscal_part+6
        !nbField = nbscal+nbscal_part+scal_partVector%nb_components+6
      else
        nbField = nbscal+nbscal_part+3
        !nbField = nbscal+nbscal_part+scal_partVector%nb_components+3
      endif
      if (differential_diff) nbField = nbField + 1
      if (write_Qc_omega) nbField = nbField + 4
      call vtkxml_init_all(nbField, nb_proc_dim , length_3D, me , coord)
      if (nbscal .gt. 0) then
        allocate(tagScal(nbscal))
        do i=1,nbscal
          call vtkxml_plus_init_field(scalArray(i), tagScal(i) )
        enddo
      endif

      if (nbscal_part .gt. 0) then
        allocate(tagScalPart(nbscal_part))
        do i=1,nbscal_part
          call vtkxml_plus_init_field(scal_partArray(i), tagScalPart(i) )
        enddo
      endif

    if (differential_diff) then
        if(scal_2<= nbscal) then
            if(.not. copyStructOnly(ScalArray(scal_2), diff_diff)) return
        else
            if(.not. copyStructOnly(Scal_partArray(scal_2-nbscal), diff_diff)) return
        end if
        diff_diff%name = 'differential_diff'
        call vtkxml_plus_init_field(diff_diff, tag_diff_diff)
    end if
    call deleteDataLayout(diff_diff)
    if (write_Qc_omega) then
       if(.not. copyStructOnly(U, Qc)) return
       if(.not. copyStructOnly(U, O1)) return
       if(.not. copyStructOnly(U, O2)) return
       if(.not. copyStructOnly(U, O3)) return
       Qc%name = 'Qcrit'
       O1%name = 'Ox'
       O2%name = 'Oy'
       O3%name = 'Oz'
       call vtkxml_plus_init_field(Qc, tag_Qc)
       call vtkxml_plus_init_field(O1, tag_O(1))
       call vtkxml_plus_init_field(O2, tag_O(2))
       call vtkxml_plus_init_field(O3, tag_O(3))
    end if
    call deleteDataLayout(Qc)
    call deleteDataLayout(O1)
    call deleteDataLayout(O2)
    call deleteDataLayout(O3)

! XXX Todo
!      if (nbscal_part .gt. 0) then
!        allocate(tagScalPartVector(scal_partVector%nb_components))
!        do i=1,nbscal_part
!          call real_vector_manipulate_comp(Scal_partVector, i, ScalPart, guru=.true.)
!          call vtkxml_plus_init_field(scalPart, tagScalPartVector(i) )
!        enddo
!      endif
! XXX End Todo

      if (imhd) then
        do i=1,3
          call vtkxml_plus_init_field(B(i), tagB(i) )
        enddo
      endif
      call vtkxml_plus_init_field(U, tagVel(1) )
      call vtkxml_plus_init_field(V, tagVel(2) )
      call vtkxml_plus_init_field(W, tagVel(3) )
    endif

    ! Init context for jet post-processing
    if ((jet_freq>0).or.(jet_deltaT>0)) then
        call parser_read('nb z position', nb_z_pos)
        allocate(z_pos(nb_z_pos))
        do k = 1, nb_z_pos
            write(parser_tag,'(a,i0)') 'z position ', k
            call parser_read(trim(parser_tag), z_pos(k))
        end do
    end if

    success = .true.

end function post_in_simu_init

!------------------------------------------------------------------------------
!! @author 
!! Luca MARRADI LIPhy 2014
!!
!! @param[in]    U               = velocity along X-axis
!! @param[in]    V               = velocity along Y-axis
!! @param[in]    W               = velocity along Z-axis
!! @param[in]    B               = magnetic field (3D vector)
!! @param[in]    scalarray       = arrays of scalars in physical space
!! @param[in]    me              = mpi rank
!! @param[out]   post_freq       = frequence off post-processing
!! @return       success         = logical equals to true is everything is right.
!! 
!! @details
!! Perform advanced post-process setup.
!------------------------------------------------------------------------------
FUNCTION post_in_simu_init_luca(sigma,me,post_freq) result(success)

    use param
    use data
    use cart_topology
    implicit none

    !== INPUT/OUTPUT DATA ==
    integer,intent(in)                                          :: me
    type(real_data_layout),intent(in)                           :: sigma
    integer, intent(out)                                        :: post_freq
    logical :: success

    !== LOCAL DATA ==
    integer                          :: nbField,i, k, nb_z_pos
    integer                          :: spec_rank
    real(WP),dimension(3)            :: length_3D
    character (len=10)               :: type_para       ! For choosing if 3D, slice or both
    character (len=30)               :: parser_tag

    success = .FALSE.

    length_3D(1)=sigma%Lx
    length_3D(2)=sigma%Ly
    length_3D(3)=sigma%Lz

    !== INIT OUTPUT FREQUENCY RELATED TO ITERATION NUMBER == 
    CALL parser_read('Simulation iterations',post_freq)
    CALL parser_read('Show basic info',show_basic)
    CALL parser_read('Sigma average XY iter',freq_sigma_avgXY)
    CALL parser_read('Sigma correlation XY iter',freq_sigma_corr)
    CALL parser_read('Waiting time iter',Nw)
    CALL parser_read('Sigma average XY time',sigma_avg_XY_deltaT)
    CALL parser_read('Save frequency',save_freq)
    CALL parser_read('Type of output',output_type)
    CALL parser_read('Write frequency',ifreq)
    if (show_basic>0) post_freq = min(show_basic, post_freq)
    if (ifreq>0) post_freq = min(ifreq, post_freq)

!!    if (hd_post>0) then
!!      CALL post_hd_init(mu,Uk)
!!    endif

    call parser_read('Show memory',show_mem)
    IF(show_mem>0) post_freq = min(show_mem, post_freq)

    ! Init output frequency related to iteration number
    post_deltaT = 1000000._WP
    init_time = 0.0_WP
    IF (parser_is_defined('Initial time')) call parser_read('Initial time', init_time)

    CALL parser_read('Final time',post_deltaT)
    CALL parser_read('Write output - time',visu_deltaT)

    IF(visu_deltaT>0) THEN
        if (ifreq>0) then
            write(*,'(a,a)') '[WARNING] Ouptut are activated with both iterations', &
                & 'and time frequency - Only iterations will be considered'
            visu_deltaT = -1.0_WP
        else
            post_deltaT = min(visu_deltaT, post_deltaT)
            ite_visu = floor(init_time/visu_deltaT) - 1
        end if
    ENDIF

    !== INIT PARAVIEW CONTEXT ==
    if (output_type .eq. 'paraview' )then
      call parser_read('slice or 3D', type_para)
      select case(trim(type_para))
      case('slice')
          paraview_slice = .true.
          paraview_3D   = .false.
      case('3D')
          paraview_3D = .true.
          paraview_slice = .false.
      case('3Dplus')
          paraview_3D = .true.
          write_Qc_omega = .true.
          paraview_slice = .false.
      case default
          paraview_slice = .true.
          paraview_3D = .true.
      end select

      call vtkxml_plus_init_field(sigma, tagVel(1) )
!!      call vtkxml_plus_init_field(V, tagVel(2) )
!!      call vtkxml_plus_init_field(W, tagVel(3) )
    endif

    success = .TRUE.

END FUNCTION post_in_simu_init_luca
!------------------------------------------------------------------------------
!> Perform custom post-process and show custom information during the run
!! Post-process are done each N iteration.
!! @author Jean-Baptiste Lagaert and Antoine Vollant
!!    @param[in]    ite             = current time iteration
!!    @param[in]    time            = current time
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in]    scalarray       = arrays of scalars in physical space
!!    @param[in]    nbscal          = the number of scalars to process
!!    @param[in]    scal_partArray  = arrays of scalars in physical space
!!    @param[in]    nbscal_part     = the number of scalars to process
!!    @param[in]    n_iter          = the number of max iteration
!!    @return       success         = logical equals to true is everything is right.
!------------------------------------------------------------------------------
function post_simu(ite,time, step_time, nbcpus,spec_rank, nbscal, scalArray, nbscal_part, scal_partArray,&
                  & U, V, W, B, imhd, n_iter) result(success)

    use data
    use post_velocity
    use differential_tools
    use rediscretization_tools
    use toolbox
    use filtering_tools , only : computeDelta,computeFilter

    ! Input/Output
    implicit none
    integer, intent(in)                                     :: ite
    real(WP), intent(in)                                    :: time, step_time
    integer, intent(in)                                     :: nbcpus
    integer, intent(in)                                     :: spec_rank
    integer, intent(in)                                     :: n_iter
    type(REAL_DATA_LAYOUT), intent(in),pointer,dimension(:) :: scalArray, scal_partArray
    type(REAL_DATA_LAYOUT), intent(in)                      :: U,V,W
    type(real_data_layout)                                  :: scalTested
    type(real_data_layout)                                  :: U_trans,V_trans,W_trans 
    type(real_data_layout),dimension(:),allocatable         :: velReDisc
    type(complex_data_layout)                               :: tempK
! XXX Todo
!    type(REAL_VECTOR_DATA_LAYOUT), intent(in)               :: scal_partVector
! XXX End Todo
    type(REAL_DATA_LAYOUT),pointer, intent(in), dimension(:):: B
    integer, intent(in)                                     :: nbscal, nbscal_part
    logical                                                 :: success
    logical                                                 :: imhd
    logical                                                 :: res
    integer                                                 :: dime,nd,filter
    real(WP)                                                :: diffu ! scalar diffusivity
    character(len=16)                                       :: ndChar

    success = .false.

    ! Show information
    if ((show_basic>0) .AND. mod(ite, show_basic).eq.0) then
        call showScalarAvg(scalarray,time,nbscal,spec_rank,'Scalar spec')
        call showScalarAvg(scal_partArray,time,nbscal_part,spec_rank,'Scalar part')
        call showVeloMinMaxAvg(U,V,W,spec_rank, nbcpus,time)
        if (.not. showEnergyAndRlambda(U,V,W,B, Uk, Vk, Wk, VelWN, mu, spec_rank,imhd,time,nbcpus,step_time,Em_Growth_rate)) then
            write (6,'(a)') 'Error on Rlambda and Energy computation' 
            return
        end if
       if (.not. compute_Integral_scale(ite, spec_rank,n_iter, time,mu)) then
           write (6,'(a)') 'Error on ReLInt and Energy computation'
           return
       end if
    end if

    ! Show memory used by the simulation
    if ((show_mem>0) .AND. mod(ite, show_mem)==0) call show_memory(spec_rank)

    ! Compute spectrum
    if ((compute_spec>0) .AND. mod(ite, compute_spec)==0) then
        if (parser_is_defined('Filter velocity for spectrum')) then
          call parser_read('Filter velocity for spectrum',nd) 
          if (.not. parseFilter(filter,spec_rank)) return
          if (.not. copyStructOnly(U,U_trans))return
          if (.not. copyStructOnly(V,V_trans))return
          if (.not. copyStructOnly(W,W_trans))return
          call computeFilter(VelWN,nd*computeDelta(U),Uk,U_trans,filter)
          call computeFilter(VelWN,nd*computeDelta(U),Vk,V_trans,filter)
          call computeFilter(VelWN,nd*computeDelta(U),Wk,W_trans,filter)
          write(ndChar,'(i3.3)') nd
          U_trans%name='velo_filt_'//trim(adjustl(ndChar))//'_'
          if (.not. compute_spectrum(U_trans,V_trans,W_trans,ite, spec_rank)) return
          call deletedatalayout(U_trans)
          call deletedatalayout(V_trans)
          call deletedatalayout(W_trans)
        endif
        if (.not. compute_all_spectrum(ite, spec_rank, n_iter,imhd)) return
    end if

    ! Compute divergence
    if (ishowdiv.gt.0 .and. mod(ite,ishowdiv).eq.0) then
       if (.not. showDivU(Uk,Vk,Wk,VelWN,spec_rank)) then
          write (6,'(a)') 'Error on divergence computation'
          return
       endif
       if (imhd) then
          if (.not. showDivB(Bk,BfieldWN,spec_rank)) then
            write (6,'(a)') 'Error on divergence computation of B field'
            return
          endif
       endif
    endif

    ! Dump fields
    if(ifreq>0) then
        if (mod(ite,ifreq).eq.0 .or. ite.eq.n_iter) then
            if ( trim(output_type) .eq. 'avs' ) then
                IF (.NOT. dataWrite(spec_rank,scalArray, scal_partArray, U, V, &
                    & W,B,imhd,ite)) THEN
                   WRITE(6,'(a)')'[ERROR] Unable to perform post-process dataWrite'
                   RETURN
                ENDIF
              !  if (.not. avsWrite(spec_rank,scalArray, scal_partArray, U, V, &
              !   & W,B,imhd,ite,n_iter)) then
              !      write (6,'(a)')'[ERROR] Failed in avs dump'
              !      return
              !  endif
            elseif ( trim(output_type) .eq. 'paraview' ) then
                if (ite == 0) then ! Save initial field in avs (not possible to restart from paraview)
                    if (.not. avsWrite(spec_rank,scalArray, scal_partArray, U, V, &
                        & W,B,imhd,ite,n_iter)) then
                        write (6,'(a)')'[ERROR] Failed in avs dump'
                        return
                    end if
                endif
                if (.not. paravieWrite(spec_rank,nbscal, scalArray, nbscal_part,scal_partArray, U, V, &
                   & W,B,imhd,tagScal,tagScalPart,tagVel,tagB, paraview_slice, paraview_3D) ) then
                    write (6,'(a)')'[ERROR] Failed in paraview dump'
                    return
                endif
            else
                write (6,'(a)')'[ERROR] Unknown type of output: abort !'
                return
            endif
            if (write_Qc_Omega) then
                if (.not. WriteQcOmega(U,Uk,Vk,Wk,VelWN,Qc,O1,O2,O3,tag_Qc,tag_O,spec_rank, ite) ) then
                    write (6,'(a)')'[ERROR] Failed in dump for Qc and Omega'
                    return
                end if
            endif
        endif
    end if

    ! Save field to be able to relaunch computations
    if((save_freq>0).and.(mod(ite, save_freq)==0)) then
        if (.not. dataWrite(spec_rank, scalArray, scal_partArray, U, V, &
            & W,B,imhd,mod(ite/save_freq,2),"sav")) then
            write (6,'(a)')'[ERROR] Failed in save field'
            return
        end if
    end if

    ! Show info of scalars
    if (show_scalar>0 .and. mod(ite, show_scalar)== 0) then
        dime = 200
        if ( size(scalarray) .gt. 1 ) then
            call showPostScalar(scalarray,nbscal,time,ite,show_scalar,dime,n_iter,spec_rank,&
                            nbcpus)  !! CTR POST
        endif
        if (.not. computeScalarsPDF(ite,scalArray,nbscal,dime,spec_rank,n_iter)) then
            write (6,'()') '[ERROR] In computing PDF'
            return
        endif
        if (.not. computeScalarsPDF(ite,scal_partArray,nbscal_part,dime,spec_rank,n_iter)) then
            write (6,'()') '[ERROR] In computing PDF'
            return
        endif
    endif

    ! Compute basic statistic on scalars
    if (scal_basic_stat>0 .and. mod(ite, scal_basic_stat)== 0) then
      do dime = 1, nbscal
        diffu = mu/schmidt(dime)
        if (diffu < 1e-8) diffu = 1._WP
        if(.not.compute_scalar_basic_stat(scalarray(dime),U, diffu, ite,n_iter,time,nbcpus,spec_rank)) then
          write (6,'(a,i0)') '[ERROR] In basic stats for scalarS ', dime
          return
        end if
      end do
      do dime = 1, nbscal_part
        diffu = mu/schmidt_part(dime)
        if (diffu < 1e-8) diffu = 1._WP
        if(.not.compute_scalar_basic_stat(scal_partArray(dime),U, diffu,ite,n_iter,time,nbcpus,spec_rank)) then
          write (6,'(a,i0)') '[ERROR] In basic stats for scalarP ', dime
          return
        end if
      end do
    endif

    ! Compute statistic moment on scalars and increment of scalar
    if (scal_moment_stat>0 .and. mod(ite, scal_moment_stat)== 0) then
      do dime = 1, nbscal
        if(.not.computeScalarMomentsStat(scalarray(dime), ite, spec_rank) ) then
          write (6,'(a,i0)') '[ERROR] In basic stats for scalarS ', dime
          return
        end if
      end do
      do dime = 1, nbscal_part
        if(.not.computeScalarMomentsStat(scal_partArray(dime), ite, spec_rank) ) then
          write (6,'(a,i0)') '[ERROR] In basic stats for scalarP ', dime
          return
        end if
      end do
    endif

    

    !Compute mean values for filtered and non filtered MHD quantities
    if ((mhd_post>0) .AND. mod(ite, mhd_post)==0) then
       call computealltheterms(U,V,W,Uk,Vk,Wk,B,Bk,VelWN,BfieldWN,nbcpus,spec_rank,time,step_time)
    end if
    if ((hd_post>0) .AND. mod(ite, hd_post)==0) then
       call computeallthetermsHD(U,V,W,Uk,Vk,Wk,VelWN,nbcpus,spec_rank,time,step_time)
    end if

    if ((iPostKolmo>0) .AND. mod(ite, iPostKolmo)==0) then
       if(imhd) then
          call compute_Mhd_stat(U,V,W,Uk,Vk,Wk,B,Bk,VelWN,spec_rank,time,step_time,ite,save_freq)
       else
          call compute_hydro_stat(U,V,W,Uk,Vk,Wk,VelWN,spec_rank,time,step_time,ite,save_freq)
       endif
    end if

    ! Jet post-process
    if ((jet_freq>0) .AND. mod(ite, jet_freq)==0) then
        dime = 200
        if(.not. compute_2D_AVG_quantities(scalArray, nbscal, scal_partArray, nbscal_part, &
            & U, V, W, ite, n_iter, spec_rank, nbcpus, time=.false.)) return
        if(spec_Rank==0) write(*,'(a)')'  [INFO] AVG2D computed'
        if (.not. compute_all_spectrum2D(ite, z_pos, spec_rank, &
            & n_iter,imhd, time=.false.)) return
        if(spec_Rank==0) write(*,'(a)')'  [INFO] Spectrum 2D computed'
        if (.not. computeScalarsPDF_2D(ite,z_pos,scalArray,nbscal,&
            & scal_partArray, nbscal_part, dime,spec_rank,n_iter, time=.false.)) return
        if(spec_Rank==0) write(*,'(a)')'  [INFO] Pdf 2D computed'
    end if

    ! Differential diffusion (special post-process for some jet cases)
    if ((diff_diff_freq>0) .AND. mod(ite, diff_diff_freq)==0) then
        if (.not.differential_diffusion(scal_1, scal_2, nbscal, nbscal_part,&
            & ScalArray, Scal_partArray, diff_diff, mix_th, tag_diff_diff,  &
            & spec_rank, ite, n_iter, time=.false.)) return
        if(spec_Rank==0) write(*,'(a)')'  [INFO] Differential diffusion computed'
    end if

    ! Apriori-test in simu for velocity
    !Compute physical quantities of filtered flow 
    if ((apriori_velo_freq>0) .AND. mod(ite, apriori_velo_freq)==0) then
       if (.not. aprioriPostVeloDNS(ite,spec_rank,U,V,W,Uk,Vk,Wk,VelWN,time,mu)) return
    endif

    ! Apriori-test in simu for a scalar
    !Compute if there is only one scalar solved in spectral or particular way
    if ((apriori_sca>0) .AND. mod(ite, apriori_sca)==0) then
       if (associated(scalArray) .and. .not. associated(scal_partArray))  then
          if (.not. copyStructOnly(scalArray(1),scalTested) ) then
             write(6,'(a)')'  [ERROR] Not able to allocate scalTested, abort!' 
             return
          endif
          scalTested%values = scalArray(1)%values
       elseif (associated(scal_partArray) .and. .not. associated(scalArray)) then
          if (.not. copyStructOnly(scal_partArray(1),scalTested) ) then
             write(6,'(a)')'  [ERROR] Not able to allocate scalTested, abort!' 
             return
          endif
          scalTested%values = scal_partArray(1)%values
       else
         if(spec_rank==0) write(6,'(a)')'  [ERROR] More than one scalar is initialized: abort!'
         return
       endif
       IF (.NOT. samelayout(scalTested,U,V,W)) THEN
         IF (.NOT. initDataLayout("tempK",tempK,(scalTested%nx/2)+1,scalTested%ny,scalTested%nz, &
            & scalTested%Lx,scalTested%Ly,scalTested%Lz,nbcpus,spec_rank,alongZ)) RETURN 
         call ftran(scalTested,tempK,res)
         IF (.NOT. res) return
         IF (.NOT. initWorkArray(scalTested,3,velReDisc)) RETURN 
         IF (.NOT. transdiscretization (tempK,scalTested,VelReDisc(1)) ) RETURN
         IF (.NOT. transdiscretization (tempK,scalTested,VelReDisc(2)) ) RETURN
         IF (.NOT. transdiscretization (tempK,scalTested,VelReDisc(3)) ) RETURN
         call deletedatalayout(tempK)
       ELSE
         IF (.NOT. initWorkArray(scalTested,3,velReDisc)) RETURN 
         velReDisc(1) = U
         velReDisc(2) = V
         velReDisc(3) = W
       END IF
       if(spec_rank==0) write(6,'(a)')'  [INFO] Apriori test are performing...'
       do iAPTest=1,nbAPTest
         if ( trim(adjustl(APTest(iAPTest))) .eq. 'scalar' ) then
           if (.not. apriori_ScalarStat(ite,velReDisc(1),velReDisc(2),velReDisc(3),scalTested,spec_rank,APfiltMin,APfiltMax,APfilter) ) then
             write (6,'(a)')'[ERROR] in computeScalarStat'
           endif
         endif
         if ( trim(adjustl(APTest(iAPTest))) .eq. 'gradient' ) then
           if (.not. apriori_Gradient(ite,velReDisc(1),velReDisc(2),velReDisc(3),scalTested,spec_rank,APfiltMin,APfiltMax,APfilter) ) then
             write (6,'(a)')'[ERROR] in apriori_Gradient'
           endif
         endif
         if ( trim(adjustl(APTest(iAPTest))) .eq. 'smagorinskydyn' ) then
           if (.not. apriori_SmagorinskyDyn(ite,velReDisc(1),velReDisc(2),velReDisc(3),scalTested,spec_rank,APfiltMin,APfiltMax,APfilter) ) then
             write (6,'(a)')'[ERROR] in apriori_SmagDyn'
           endif     
         endif
         if ( trim(adjustl(APTest(iAPTest))) .eq. 'clark' ) then
           if (.not. apriori_Clark(ite,velReDisc(1),velReDisc(2),velReDisc(3),scalTested,spec_rank,APfiltMin,APfiltMax,APfilter) ) then
             write (6,'(a)')'[ERROR] in apriori_Clark'
           endif     
         endif
         if ( trim(adjustl(APTest(iAPTest))) .eq. 'smagorinsky' ) then
           if (.not. apriori_Smagorinsky(ite,velReDisc(1),velReDisc(2),velReDisc(3),scalTested,spec_rank,APfiltMin,APfiltMax,APfilter) ) then
             write (6,'(a)')'[ERROR] in apriori_Smag'
           endif     
         endif
         if ( trim(adjustl(APTest(iAPTest))) .eq. 'clarkfabre' ) then
           if (.not. apriori_ClarkFabre(ite,velReDisc(1),velReDisc(2),velReDisc(3),scalTested,spec_rank,APfiltMin,APfiltMax,APfilter) ) then
             write (6,'(a)')'[ERROR] in apriori_ClarkFabre'
           endif     
         endif
       enddo
       call deletedatalayout(scalTested)
       res = deleteWorkArray(velReDisc)
    end if


    success = .true.

end function post_simu

!------------------------------------------------------------------------------
!> Perform custom post-process and show custom information during the run
!! Post-process are done each N iteration.
!! @author Luca MARRADI LIPhy 2014
!!    @param[in]    ite             = current time iteration
!!    @param[in]    time            = current time
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in]    scalarray       = arrays of scalars in physical space
!!    @param[in]    nbscal          = the number of scalars to process
!!    @param[in]    scal_partArray  = arrays of scalars in physical space
!!    @param[in]    nbscal_part     = the number of scalars to process
!!    @param[in]    n_iter          = the number of max iteration
!!    @return       success         = logical equals to true is everything is
!right.
!------------------------------------------------------------------------------
FUNCTION post_simu_luca(ite,time,step_time,nbcpus,spec_rank,mult,sigma,sigma_old,n_activ,n_iter,corrXYw) result(success)

    USE data
    USE post_velocity
    USE differential_tools
    USE rediscretization_tools
    USE toolbox
    USE fileio
    USE filtering_tools , only : computeDelta,computeFilter

    IMPLICIT NONE

    !== INPUT/OUTPUT DATA ==
    integer, intent(in)                       :: ite
    real(WP), intent(in)                      :: time, step_time,corrXYw
    integer, intent(in)                       :: nbcpus
    integer, intent(in)                       :: spec_rank
    integer, intent(in)                       :: n_iter
    type(REAL_DATA_LAYOUT), intent(INOUT)     :: sigma_old,sigma,n_activ,mult

    !== LOCAL DATA ==
    type(complex_data_layout)                 :: tempK
    REAL(WP)                                  :: avgXY,avg_act,corrXY
    logical                                   :: success
    logical                                   :: res
    integer                                   :: dime,nd,filter,file_id
    integer                                   :: i,j,k
    character(len=16)                         :: ndChar

    success = .FALSE.

    ! Show information
!    if ((show_basic>0) .AND. mod(ite, show_basic).eq.0) then
!        call showScalarAvg(scalarray,time,nbscal,spec_rank,'Scalar spec')
!        call showScalarAvg(scal_partArray,time,nbscal_part,spec_rank,'Scalar part')
!        call showVeloMinMaxAvg(U,V,W,spec_rank, nbcpus,time)
!        if (.not. showEnergyAndRlambda(U,V,W,B, Uk, Vk, Wk, VelWN, mu, spec_rank,imhd,time,nbcpus,step_time,Em_Growth_rate)) then
!            write (6,'(a)') 'Error on Rlambda and Energy computation' 
!            return
!        end if
!       if (.not. compute_Integral_scale(ite, spec_rank,n_iter, time,mu)) then
!           write (6,'(a)') 'Error on ReLInt and Energy computation'
!           return
!       end if
!    end if

    ! == SHOWS MEMORY USED BY THE SIMULATION ==
    IF ((show_mem>0) .AND. mod(ite, show_mem)==0) call show_memory(spec_rank)

    ! Compute spectrum
!    if ((compute_spec>0) .AND. mod(ite, compute_spec)==0) then
!        if (parser_is_defined('Filter velocity for spectrum')) then
!          call parser_read('Filter velocity for spectrum',nd) 
!          if (.not. parseFilter(filter,spec_rank)) return
!          if (.not. copyStructOnly(U,U_trans))return
!          if (.not. copyStructOnly(V,V_trans))return
!          if (.not. copyStructOnly(W,W_trans))return
!          call computeFilter(VelWN,nd*computeDelta(U),Uk,U_trans,filter)
!          call computeFilter(VelWN,nd*computeDelta(U),Vk,V_trans,filter)
!          call computeFilter(VelWN,nd*computeDelta(U),Wk,W_trans,filter)
!          write(ndChar,'(i3.3)') nd
!          U_trans%name='velo_filt_'//trim(adjustl(ndChar))//'_'
!          if (.not. compute_spectrum(U_trans,V_trans,W_trans,ite, spec_rank)) return
!          call deletedatalayout(U_trans)
!          call deletedatalayout(V_trans)
!          call deletedatalayout(W_trans)
!        endif
!        if (.not. compute_all_spectrum(ite, spec_rank, n_iter,imhd)) return
!    end if

    !== DUMP FIELDS == 
    IF(ifreq>0) THEN
        if (mod(ite,ifreq).eq.0 .or. ite.eq.n_iter) then
            if ( trim(output_type) .eq. 'avs' ) then
                if (.NOT.  dataWrite_luca(spec_rank,sigma,ite) ) then
                   write(6,'(a)')'[ERROR] Unable to perform post-process dataWrite'
                   return
                endif
            elseif ( trim(output_type) .eq. 'paraview' ) then
                if (ite == 0) then 
                    !== SAVE INITIAL FIELDS IN AVS (not possible to restart from paraview)
                    if (.not.avsWrite_luca(spec_rank,sigma,ite,n_iter)) then
                        write (6,'(a)')'[ERROR] Failed in avs dump'
                        return
                    end if
                endif
                if (.not. paravieWrite_luca(spec_rank,sigma,tagVel,paraview_slice,paraview_3D) ) then
                    write (6,'(a)')'[ERROR] Failed in paraview dump'
                    return
                endif
            else
                write (6,'(a)')'[ERROR] Unknown type of output: abort !'
                return
            endif
        endif
    ENDIF

    ! Save field to be able to relaunch computations
    IF((save_freq>0).and.(mod(ite, save_freq)==0)) THEN
        if (.not.dataWrite_luca(spec_rank,sigma,mod(ite/save_freq,2),"sav")) then
            write (6,'(a)')'[ERROR] Failed in save field'
            return
        end if
    ENDIF

    ! Show info of scalars
!    if (show_scalar>0 .and. mod(ite, show_scalar)== 0) then
!        dime = 200
!        if ( size(scalarray) .gt. 1 ) then
!            call showPostScalar(scalarray,nbscal,time,ite,show_scalar,dime,n_iter,spec_rank,&
!                            nbcpus)  !! CTR POST
!        endif
!        if (.not. computeScalarsPDF(ite,scalArray,nbscal,dime,spec_rank,n_iter)) then
!            write (6,'()') '[ERROR] In computing PDF'
!            return
!        endif
!        if (.not. computeScalarsPDF(ite,scal_partArray,nbscal_part,dime,spec_rank,n_iter)) then
!            write (6,'()') '[ERROR] In computing PDF'
!            return
!        endif
!    endif

!
   !== COMPUTES SIGMA XY CORRELATION
   IF (freq_sigma_corr>0 .AND. (mod((ite-Nw), freq_sigma_corr)== 0).AND.(ite>Nw)) THEN
       avgXY   = computeFieldAvg(sigma,spec_rank)
       corrXY  = computeFieldAvg(mult,spec_rank)- avgXY*corrXYw
     if(spec_rank==0) then
      file_id = 410
      open(unit=file_id, file='sigma_corrXY.dat',form="FORMATTED",ACCESS='APPEND')
       write(file_id,'(g15.8)') corrXY
      close(file_id)
     endif
   ENDIF

    !== COMPUTES THE SIGMA XY AVERAGE 
    IF (freq_sigma_avgXY>0 .AND. mod(ite, freq_sigma_avgXY)== 0) THEN
       avgXY   = computeFieldAvg(sigma,spec_rank)
       !avg_act = computeFieldAvg(n_activ,spec_rank)
     if(spec_rank.EQ.0) then
      print*, 'SPEC_RANK:',spec_rank
      file_id = 400!iopen()
      open(unit=file_id, file='sigma_avgXY.dat',form="FORMATTED",ACCESS='APPEND')
      open(unit=501, file='activity.dat',form="FORMATTED",ACCESS='APPEND')
       write(file_id,'(g15.8)') avgXY
       !write(501,'(g15.8)') avg_act
      close(file_id)
     endif
  ENDIF

   success = .TRUE.
END FUNCTION post_simu_luca


!------------------------------------------------------------------------------
!> Perform custom post-process and show custom information during the run
!! Post-process are done each delta(t) time (and thus does not depend of the
!! iteration number)
!! @author Jean-Baptiste Lagaert and Antoine Vollant
!!    @param[in]    ite             = current time iteration
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in]    scalarray       = arrays of scalars in physical space
!!    @param[in]    nbscal          = the number of scalars to process
!!    @param[in]    scal_partArray  = arrays of scalars in physical space
!!    @param[in]    nbscal_part     = the number of scalars to process
!!    @param[in]    n_iter          = the number of max iteration
!!    @return       success         = logical equals to true is everything is right.
!------------------------------------------------------------------------------
function post_simu_deltaT(t,nbcpus,spec_rank, nbscal, scalArray, nbscal_part, scal_partArray,&
                  & U, V, W, B, imhd, n_iter) result(success)

    use data
    use rediscretization_tools

    ! Input/Output
    implicit none
    real(WP), intent(in)                                    :: t
    integer, intent(in)                                     :: nbcpus
    integer, intent(in)                                     :: spec_rank
    integer, intent(in)                                     :: n_iter
    type(REAL_DATA_LAYOUT), intent(in),pointer,dimension(:) :: scalArray, scal_partArray
! XXX Todo
!    type(REAL_VECTOR_DATA_LAYOUT), intent(in)               :: scal_partVector
! XXX End Todo
    type(REAL_DATA_LAYOUT), intent(in)                      :: U,V,W
    type(REAL_DATA_LAYOUT),pointer, intent(in), dimension(:):: B
    integer, intent(in)                                     :: nbscal, nbscal_part
    ! Local variables
    logical                                                 :: success
    logical                                                 :: imhd
    integer                                                 :: dime
    ! For LES jet post-processing --- Implementation has to be improved !!!! But right now, I have not time for that !
    !type(real_data_layout),dimension(:),allocatable     :: scalFilDNS
!    type(COMPLEX_DATA_LAYOUT)                           :: tempK

    success = .false.

    ! Compute spectrum
    if ((spec_deltaT>0) .AND. (t>ite_spec*spec_deltaT)) then
        ite_spec = ite_spec + 1
        if (.not. compute_all_spectrum(ite_spec, spec_rank, n_iter,imhd, time=.true.)) return
    end if

    ! Differential diffusion (special post-process fo some jet cases)
    if ((diff_diff_deltaT>0) .AND. (t>ite_diff_diff*diff_diff_deltaT)) then
        ite_diff_diff = ite_diff_diff + 1
        if (.not.differential_diffusion(scal_1, scal_2, nbscal, nbscal_part,&
            & ScalArray, Scal_partArray, diff_diff, mix_th, tag_diff_diff,  &
            & spec_rank, ite_diff_diff, n_iter, time=.true.)) return
        if(spec_Rank==0) write(*,'(a)')'  [INFO] Differential diffusion computed'
    end if

    ! Jet post-process
    if ((jet_deltaT>0) .AND. (t>ite_jet*jet_deltaT)) then
        ite_jet = ite_jet + 1
        dime = 200
        if(.not. compute_2D_AVG_quantities(scalArray, nbscal, scal_partArray, nbscal_part, &
            & U, V, W, ite_jet, n_iter, spec_rank, nbcpus, time=.true.)) return
        if(spec_Rank==0) write(*,'(a)')'  [INFO] Avg 2D computed'
        if (.not. compute_all_spectrum2D(ite_jet, z_pos, spec_rank, &
            & n_iter,imhd, time=.true.)) return
        if(spec_Rank==0) write(*,'(a)')'  [INFO] Spectrum 2D computed'
        if (.not. computeScalarsPDF_2D(ite_jet,z_pos,scalArray,nbscal,&
            & scal_partArray, nbscal_part, dime,spec_rank,n_iter, time=.true.)) return
        if(spec_Rank==0) write(*,'(a)')'  [INFO] Pdf 2D computed'

        !!! LES jet post-processing --- Implementation has to be improved !!!! But right now, I have not time for that !
!        if (.not. initDataLayout("TempK", tempK,(ScalArray(1)%nx/2)+1,ScalArray(1)%ny,ScalArray(1)%nz, &
!                     & ScalArray(1)%Lx,ScalArray(1)%Ly,ScalArray(1)%Lz,nbcpus,spec_rank,alongZ)) then
!            write(6,'(a,i0,a)')'[ERROR] initDataLayout for TYPE (COMPLEX_DATA_LAYOUT) : not enought memory!'
!            return
!        endif
!        call ftran(ScalArray(1),tempK,res)
!        IF (.NOT. samelayout(ScalArray(1),ScalArray(2))) THEN
!            IF (.NOT. initWorkArray(ScalArray(2),1,scalfilDNS) ) THEN
!              IF (.NOT. transdiscretization (tempK,ScalArray(2),scalfilDNS(1)) ) RETURN
!        END IF
!        call deletedatalayout(tempK)
!        scalfilDNS(1)%name = "Filt_DNS"
!        if(.not. compute_2D_AVG_quantities(scalfilDNS, 1, scal_partArray, nbscal_part, &    !!!! ASK TO JB : JUST PUT 0 INSTEAD OF nbscal_part
!            & U, V, W, ite_jet, n_iter, spec_rank, nbcpus, time=.true.)) return
!        if (.not. computeScalarsPDF_2D(ite_jet,z_pos,scalfilDNS,1,&
!            & scal_partArray, nbscal_part, dime,spec_rank,n_iter, time=.true.)) return   !!!! ASK TO JB : JUST PUT 0 INSTEAD OF nbscal_part
!        if ( .not. deleteWorkArray(scalfilDNS) ) return
        !!! END -- LES jet post-processing --- Implementation has to be improved !!!! But right now, I have not time for that !

    end if

    ! Dump fields
    if ((visu_deltaT>0).and.(t>ite_visu*visu_deltaT)) then
        ite_visu = ite_visu + 1
        if ( trim(output_type) .eq. 'avs' ) then
                if (.not. avsWrite(spec_rank,scalArray, scal_partArray, U, V, &
                 & W,B,imhd,ite_visu,n_iter)) then
                    write (6,'(a)')'[ERROR] Failed in avs dump'
                    return
                endif
        elseif ( trim(output_type) .eq. 'paraview' ) then
            if (.not. paravieWrite(spec_rank,nbscal, scalArray, nbscal_part,scal_partArray, U, V, &
               & W,B,imhd,tagScal,tagScalPart,tagVel,tagB, paraview_slice, paraview_3D) ) then
                write (6,'(a)')'[ERROR] Failed in paraview dump'
                return
            endif
        else
            write (6,'(a)')'[ERROR] Unknown type of output: abort !'
            return
        endif
    end if

!    ! Compute divergence
!    if (ishowdiv.gt.0 .and. mod(ite,ishowdiv).eq.0) then
!       if (.not. showDivU(spec_rank)) then
!          write (6,'(a)') 'Error on divergence computation'
!          return
!       endif
!       if (imhd) then
!          if (.not. showDivB(spec_rank)) then
!            write (6,'(a)') 'Error on divergence computation of B field'
!            return
!          endif
!       endif
!    endif
!
!    ! Show info of scalars
!    if (show_scalar>0 .and. mod(ite, show_scalar)== 0) then
!        dime = 50
!        if (.not. computeScalarsPDF(ite,scalArray,nbscal,dime,spec_rank,n_iter)) then
!            write (6,'()') '[ERROR] '
!            return
!        endif
!    endif
!
!    !Compute mean values for filtered and non filtered MHD quantities
!    if ((mhd_post>0) .AND. mod(ite, mhd_post)==0) then
!       call computealltheterms(U,V,W,B,VelWN,BfieldWN,nbcpus,spec_rank,sim_time,step_time)
!    end if


    success = .true.

end function post_simu_deltaT

!------------------------------------------------------------------------------
!!  @author 
!!  Luca MARRADI LIPhy 2014
!!
!!  @param[in]    ite             = current time iteration
!!  @param[in]    spec_rank       = mpi rank inside spectral communicator
!!  @param[in]    n_iter          = the number of max iteration
!!  @return       success         = logical equals to true is everything is right.
!!  @details
!! Perform custom post-process and show custom information during the run
!! Post-process are done each delta(t) time (and thus does not depend of the
!! iteration number)
!------------------------------------------------------------------------------
FUNCTION post_simu_deltaT_luca(t,nbcpus,spec_rank,sigma, n_iter) result(success)

    use data
    use rediscretization_tools

    !== INPUT/OUTPUT ==
    implicit none
    real(WP), intent(in)                                    :: t
    integer, intent(in)                                     :: nbcpus
    integer, intent(in)                                     :: spec_rank
    integer, intent(in)                                     :: n_iter
    type(REAL_DATA_LAYOUT), intent(in)                      :: sigma

    !== LOCAL DATA ==
    logical                                                 :: success
    integer                                                 :: dime

    success = .FALSE.

    ! Compute spectrum
!    if ((spec_deltaT>0) .AND. (t>ite_spec*spec_deltaT)) then
!        ite_spec = ite_spec + 1
!        if (.not. compute_all_spectrum(ite_spec, spec_rank, n_iter,imhd, time=.true.)) return
!    end if

    ! Differential diffusion (special post-process fo some jet cases)
!    if ((diff_diff_deltaT>0) .AND. (t>ite_diff_diff*diff_diff_deltaT)) then
!        ite_diff_diff = ite_diff_diff + 1
!        if (.not.differential_diffusion(scal_1, scal_2, nbscal, nbscal_part,&
!            & ScalArray, Scal_partArray, diff_diff, mix_th, tag_diff_diff,  &
!            & spec_rank, ite_diff_diff, n_iter, time=.true.)) return
!        if(spec_Rank==0) write(*,'(a)')'  [INFO] Differential diffusion computed'
!    end if

    ! Dump fields
    if ((visu_deltaT>0).and.(t>ite_visu*visu_deltaT)) then
        ite_visu = ite_visu + 1
        if ( trim(output_type) .eq. 'avs' ) then
                if (.not. avsWrite_luca(spec_rank,sigma,ite_visu,n_iter)) then
                    write (6,'(a)')'[ERROR] Failed in avs dump'
                    return
                endif
        elseif ( trim(output_type) .eq. 'paraview' ) then
            if (.not. paravieWrite_luca(spec_rank,sigma,tagVel,paraview_slice, paraview_3D) ) then
                write (6,'(a)')'[ERROR] Failed in paraview dump'
                return
            endif
        else
            write (6,'(a)')'[ERROR] Unknown type of output: abort !'
            return
        endif
    end if
   
    !== COMPUTES THE SIGMA XY AVERAGE 
   IF((sigma_avg_XY_deltaT>0) .AND. (t>ite_sigma_avg_XY*sigma_avg_XY_deltaT)) THEN
    avgXY = computeFieldAvg(sigma,spec_rank)
    ite_sigma_avg_XY = ite_sigma_avg_XY + 1
    IF(spec_rank==0) THEN
      OPEN(unit=400,file='sigma_avgXY.dat',form="FORMATTED",ACCESS='APPEND')
      WRITE(400,'(g15.8)') avgXY
      CLOSE(400)
    ENDIF
   ENDIF

    success = .TRUE.

END FUNCTION post_simu_deltaT_luca
!------------------------------------------------------------------------------
!> Delete post-process setup.
!! @author Antoine Vollant
!!    @return       success         = logical equals to true is everything is right.
!------------------------------------------------------------------------------
function post_in_simu_del(nbscal,nbscal_part) result(success)

    use vtkxml_plus
    use io_interface
    implicit none

    !I/O data
    logical             :: success
    integer,intent(in)  :: nbscal
    integer,intent(in)  :: nbscal_part
! XXX Todo
!    integer,intent(in)  :: nbscal_partVector
! XXX End Todo

    success = .false.

    if (output_type .eq. 'paraview' )then
      call vtkxml_finish
      if (nbscal .gt. 0) deallocate(tagScal)
      if (nbscal_part .gt. 0) deallocate(tagScalPart)
! XXX Todo
!      if (nbscal_partVector .gt. 0) deallocate(tagScalPartVector)
! XXX End Todo
    endif
    if (allocated(z_pos)) deallocate(z_pos)
    if (allocated(APTest)) deallocate(APTest)
    if (.not. interface_io_close() ) write (6,'(a)') '[WARNING] HDF5 not correctly closed'

    success = .true. 

end function post_in_simu_del


end module post_in_simu
!> @}
