!> @addtogroup post_process
!! @{

!------------------------------------------------------------------------------

!
! MODULE: post_velocity
!
!
!
!> @details module which contains post-processing routines for the velocity field 
!> @author
!
!------------------------------------------------------------------------------

module post_velocity 

    use precision_tools
    use datalayout
    use toolbox
    use stat_tools
    use transforms_tools
    use post_lib
    use physical_values_tools
    use wavenumber_tools
    use velmodels

    implicit none
    private

    ! ==== Interface procedures ====

    public  :: compute_flow_stat
    public  :: showEnergyAndRlambda
    public  :: velo_spectrum_rescaled
    public  :: checkGradientModels
    public  :: aprioriRGMVel
    public  :: aprioriSikmoinsSkj
    public  :: testDissipationduidxkdujdxk
    public  :: gradientvelmodOE
    public  :: velo_test_interpol
    public  :: aprioriPostVeloDNS
    public  :: testTransfertNrjQDM 

    real(WP) :: InitialEm,Growth_B
contains


!-------------------------------------------------------------------------------------------
!> 
!! Compute quantities from filtered flow  
!!
!! @author Jean-Baptiste Lagaert, LEGI
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in]    nbcpus          = number of process in the cpu pool
!!    @param[in]    U               = Component of velocity
!!    @param[in]    V               = Component of velocity
!!    @param[in]    W               = Component of velocity
!!    @return       res             = logical value
!
!> @detail
!! Perform a lot of computation in order to characterize the flow turbulence.
!! Some of the computation are optional (and are related to optional argument)
!-------------------------------------------------------------------------------------------
function aprioriPostVeloLES(ite,spec_rank,U,Uk,Vk,Wk,WN,sim_time,visco) result(success)

    use mpilayout_tools
    use filtering_tools , only : computeDelta

    ! Input/output
    integer, intent(in)                     :: spec_rank
    type(REAL_DATA_LAYOUT), intent(in)      :: U
    type(COMPLEX_DATA_LAYOUT), intent(in)   :: Uk,Vk,Wk
    type(WAVENUMBERS), intent(in)           :: WN
    real(WP),intent(in)                     :: sim_time,visco
    integer,intent(in)                      :: ite
    logical                                 :: success

    ! Local field
    integer                                 :: nd, nbcpus
    real(WP)                                :: E_kine, lambda, R_lambda
    real(WP)                                :: epsil, eta,helicity,enstrophy
    type(COMPLEX_DATA_LAYOUT)               :: Uk_trans,Vk_trans,Wk_trans
    type(COMPLEX_DATA_LAYOUT)               :: Vort1k,Vort2k,Vort3k
    type(REAL_DATA_LAYOUT)                  :: U_trans,V_trans,W_trans
    real(WP), dimension(3)                  :: u_prime_carre, di_u_i
    real(WP)                                :: eps = 1e-10
    character(len=8)                        :: ndChar
    character(len=256)                      :: nameFic

    success = .false.
    nbcpus = getnbcpus()
    if (.not. copyStructOnly(Uk,Uk_trans))return
    if (.not. copyStructOnly(Vk,Vk_trans))return
    if (.not. copyStructOnly(Wk,Wk_trans))return
    if (.not. copyStructOnly(Uk,Vort1k))return
    if (.not. copyStructOnly(Vk,Vort2k))return
    if (.not. copyStructOnly(Wk,Vort3k))return
    if (.not. copyStructOnly(U,U_trans))return
    if (.not. copyStructOnly(U,V_trans))return
    if (.not. copyStructOnly(U,W_trans))return
    !Compute dissipation
    call computeCurlK(WN,Uk_trans,Vk_trans,Wk_trans,Vort1k,Vort2k,Vort3k) !compute vorticity
    !Compute kinetic  spectrum
    Uk_trans%name='kinetic_nd_'//trim(adjustl(ndChar))//'_'
    if (.not. compute_spectrum(Uk_trans,Vk_trans,Wk_trans,ite, spec_rank,WN)) return 
    !Compute helicity spectrum
    Uk_trans%name='helicity_nd_'//trim(adjustl(ndChar))//'_'
    if (.not. compute_spectrum(Uk_trans,Vk_trans,Wk_trans,Vort1k,Vort2k,Vort3k,ite,spec_rank,WN)) return
    !Compute enstrophy spectrum
    Vort1k%name='enstrophy_nd_'//trim(adjustl(ndChar))//'_'
    if (.not. compute_spectrum(Vort1k,Vort2k,Vort3k,ite,spec_rank,WN)) return
    !Compute global values
    helicity     = computeHelicity(WN,U_trans,V_trans,W_trans,nbcpus,spec_rank)
    enstrophy    = computeEnstrophy(WN,U_trans,V_trans,W_trans,nbcpus,spec_rank)
    ! Sum up all statistic
    write (nameFic,'(a,i3.3,a)') 'Flow_stat_nd_',nd,'.out'
    if (spec_rank .eq. 0) then
      if ( sim_time .lt. eps ) then
          open(10,file=nameFic,form='formatted')
          write(10,'(a)') &
              & "#t_simu E_kine helicity enstrophy lambda Re_lambda epsil &
              &eta  u'^2 v'^2 w'^2 <(du/dx)^2> &
              &<(dv/dy)^2> <(dw/dz)^2> dissipationSM"
      else
         open(10,file=nameFic,form='formatted',position='append')
      end if
      write(10,'(15(g15.8,x))') sim_time, E_kine, helicity ,enstrophy, lambda, R_lambda, epsil, eta, u_prime_carre, di_u_i
      close(10)
    endif

    call deletedatalayout(Uk_trans)
    call deletedatalayout(Vk_trans)
    call deletedatalayout(Wk_trans)
    call deletedatalayout(U_trans)
    call deletedatalayout(V_trans)
    call deletedatalayout(W_trans)
    call deletedatalayout(Vort1k)
    call deletedatalayout(Vort2k)
    call deletedatalayout(Vort3k)

    success = .true.

end function aprioriPostVeloLES


!-------------------------------------------------------------------------------------------
!> 
!! Compute quantities from filtered flow  
!!
!! @author Jean-Baptiste Lagaert, LEGI
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in]    nbcpus          = number of process in the cpu pool
!!    @param[in]    U               = Component of velocity
!!    @param[in]    V               = Component of velocity
!!    @param[in]    W               = Component of velocity
!!    @return       res             = logical value
!
!> @detail
!! Perform a lot of computation in order to characterize the flow turbulence.
!! Some of the computation are optional (and are related to optional argument)
!-------------------------------------------------------------------------------------------
function aprioriPostVeloDNS(ite,spec_rank,U,V,W,Uk,Vk,Wk,WN,sim_time,visco) result(success)

    use mpilayout_tools
    use filtering_tools , only : computeDelta

    ! Input/output
    integer, intent(in)                     :: spec_rank
    type(REAL_DATA_LAYOUT), intent(in)      :: U,V,W
    type(COMPLEX_DATA_LAYOUT), intent(in)   :: Uk,Vk,Wk
    type(WAVENUMBERS), intent(in)           :: WN
    real(WP),intent(in)                     :: sim_time,visco
    integer,intent(in)                      :: ite
    logical                                 :: success

    ! Local field
    integer                                 :: filter, nd, ndMin, ndMax, ndDelta, nbcpus
    real(WP)                                :: E_kine, lambda, R_lambda
    real(WP)                                :: epsil, eta,helicity,enstrophy
    type(COMPLEX_DATA_LAYOUT)               :: Uk_trans,Vk_trans,Wk_trans
    type(COMPLEX_DATA_LAYOUT)               :: Vort1k,Vort2k,Vort3k
    type(REAL_DATA_LAYOUT)                  :: U_trans,V_trans,W_trans
    type(REAL_DATA_LAYOUT)                  :: S11,S22,S33,T12
    type(REAL_DATA_LAYOUT)                  :: S12,S13,S23,T13
    type(REAL_DATA_LAYOUT)                  :: T11,T22,T33,T23
    real(WP), dimension(3)                  :: u_prime_carre, di_u_i
    real(WP)                                :: eps = 1e-10,dissipationSM
    character(len=8)                        :: ndChar,iteChar
    character(len=256)                      :: nameFic
    logical                                 :: res

    success = .false.
    nbcpus = getnbcpus()
    if (.not. parser_is_defined('Range for filtering velo')) then
      print *,'[ERROR] Range for filtering velo'
      return
    endif
    if (.not. GetNxNyNz('Range for filtering velo',ndMin,ndMax,ndDelta)) return
    nameFic="Type of filter for velocity a priori test"
    if (.not. parseFilter(filter,spec_rank,other=trim(adjustl(nameFic)))) return
!    filter = 1
    if (.not. copyStructOnly(Uk,Uk_trans))return
    if (.not. copyStructOnly(Vk,Vk_trans))return
    if (.not. copyStructOnly(Wk,Wk_trans))return
    if (.not. copyStructOnly(Uk,Vort1k))return
    if (.not. copyStructOnly(Vk,Vort2k))return
    if (.not. copyStructOnly(Wk,Vort3k))return
    if (.not. copyStructOnly(U,U_trans))return
    if (.not. copyStructOnly(U,V_trans))return
    if (.not. copyStructOnly(U,W_trans))return
    if (.not. copyStructOnly(U,S11))return
    if (.not. copyStructOnly(U,S22))return
    if (.not. copyStructOnly(U,S33))return
    if (.not. copyStructOnly(U,S12))return
    if (.not. copyStructOnly(U,S13))return
    if (.not. copyStructOnly(U,S23))return
    if (.not. copyStructOnly(U,T11))return
    if (.not. copyStructOnly(U,T22))return
    if (.not. copyStructOnly(U,T33))return
    if (.not. copyStructOnly(U,T12))return
    if (.not. copyStructOnly(U,T13))return
    if (.not. copyStructOnly(U,T23))return
    do nd=ndMin,ndMax,ndDelta
      write(ndChar,'(i3.3)') nd
      write(iteChar,'(i6.6)') ite
      call computeFilter(WN,nd*computeDelta(U),Uk,Uk_trans,filter)
      call computeFilter(WN,nd*computeDelta(U),Vk,Vk_trans,filter)
      call computeFilter(WN,nd*computeDelta(U),Wk,Wk_trans,filter)
      call btran(Uk_trans,U_trans,res)
      if (.not. res) return
      call btran(Vk_trans,V_trans,res)
      if (.not. res) return
      call btran(Wk_trans,W_trans,res)
      if (.not. res) return
      !Compute dissipation
      call computeT_ijVecA(U,V,W,U_trans,V_trans,W_trans,WN,T11,T12,T13,T22,T23,T33,filter,res,nd*computeDelta(U))
      call computeStressTensor(Uk_trans,Vk_trans,Wk_trans,WN,S11,S12,S13,S22,S23,S33,res)
      S11%values = S11%values*T11%values + S22%values*T22%values + S33%values*T33%values +&
          2.0_WP *(S12%values*T12%values + S13%values*T13%values + S23%values*T23%values )
      dissipationSM = -computeFieldAvg(S11,spec_rank)
      call computeCurlK(WN,Uk_trans,Vk_trans,Wk_trans,Vort1k,Vort2k,Vort3k) !compute vorticity
      !Compute kinetic  spectrum
      Uk_trans%name='kinetic_nd_'//trim(adjustl(ndChar))//'_'
      if (.not. compute_spectrum(Uk_trans,Vk_trans,Wk_trans,ite, spec_rank,WN)) return 
      !Compute helicity spectrum
      Uk_trans%name='helicity_nd_'//trim(adjustl(ndChar))//'_'
      if (.not. compute_spectrum(Uk_trans,Vk_trans,Wk_trans,Vort1k,Vort2k,Vort3k,ite,spec_rank,WN)) return
      !Compute enstrophy spectrum
      Vort1k%name='enstrophy_nd_'//trim(adjustl(ndChar))//'_'
      if (.not. compute_spectrum(Vort1k,Vort2k,Vort3k,ite,spec_rank,WN)) return
      !Compute global values
      helicity     = computeHelicity(WN,U_trans,V_trans,W_trans,nbcpus,spec_rank)
      enstrophy    = computeEnstrophy(WN,U_trans,V_trans,W_trans,nbcpus,spec_rank)
      call compute_flow_stat(spec_rank, nbcpus, res,U_trans,V_trans,W_trans,&
                            & Uk_trans, Vk_trans,Wk_trans, WN, visco,&
                            & E_kine, lambda, R_lambda,u_prime_carre, di_u_i,&
                            & epsil,eta)
      ! Sum up all statistic
      write (nameFic,'(a,i3.3,a)') 'Flow_stat_nd_',nd,'.out'
      if (spec_rank .eq. 0) then
        if ( sim_time .lt. eps ) then
            open(10,file=nameFic,form='formatted')
            write(10,'(a)') &
                & "#t_simu E_kine helicity enstrophy lambda Re_lambda epsil &
                &eta  u'^2 v'^2 w'^2 <(du/dx)^2> &
                &<(dv/dy)^2> <(dw/dz)^2> dissipationSM"
        else
           open(10,file=nameFic,form='formatted',position='append')
        end if
        write(10,'(15(g15.8,x))') sim_time, E_kine, helicity ,enstrophy, lambda, R_lambda, epsil, eta, u_prime_carre, di_u_i,dissipationSM 
        close(10)
      endif
    enddo

    call deletedatalayout(Uk_trans)
    call deletedatalayout(Vk_trans)
    call deletedatalayout(Wk_trans)
    call deletedatalayout(U_trans)
    call deletedatalayout(V_trans)
    call deletedatalayout(W_trans)
    call deletedatalayout(S11)
    call deletedatalayout(S22)
    call deletedatalayout(S33)
    call deletedatalayout(S12)
    call deletedatalayout(S13)
    call deletedatalayout(S23)
    call deletedatalayout(T11)
    call deletedatalayout(T22)
    call deletedatalayout(T33)
    call deletedatalayout(T12)
    call deletedatalayout(T13)
    call deletedatalayout(T23)
    call deletedatalayout(Vort1k)
    call deletedatalayout(Vort2k)
    call deletedatalayout(Vort3k)

    success = .true.

end function aprioriPostVeloDNS


!-------------------------------------------------------------------------------------------
!> Compute flow statistic : dissipation rate, Taylor micro-scale (and associated
!! Reynolds number), Kolmogorov and Batchelor WN number
!!
!! @author Jean-Baptiste Lagaert, LEGI
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in]    nbcpus          = number of process in the cpu pool
!!    @param[in]    U               = Component of velocity
!!    @param[in]    V               = Component of velocity
!!    @param[in]    W               = Component of velocity
!!    @return       res             = logical value
!
!> @detail
!! Perform a lot of computation in order to characterize the flow turbulence.
!! Some of the computation are optional (and are related to optional argument)
!-------------------------------------------------------------------------------------------

subroutine compute_flow_stat(spec_rank, nbcpus, res,    & ! Working information
                        & U,V,W, Uk, Vk,Wk, WN, visco,  & ! Input
                        & E_kine, lambda, R_lambda,     & ! basic post-process
                        & u_prime_carre, di_u_i,        & ! about isotropy
                        & epsil, eta)                     ! advanced post-process)

    use datalayout
    use differential_tools , only : computeDerivationField

    ! Input/output
    logical, intent(out)                    :: res
    integer, intent(in)                     :: spec_rank, nbcpus
    type(REAL_DATA_LAYOUT), intent(in)      :: U,V,W
    type(COMPLEX_DATA_LAYOUT), intent(in)   :: Uk,Vk,Wk
    type(WAVENUMBERS), intent(in)           :: WN
    real(WP), intent(in)                    :: visco
    real(WP), intent(out)                   :: E_kine, lambda, R_lambda
    real(WP), dimension(:), intent(out)     :: u_prime_carre, di_u_i
    real(WP), intent(out), optional         :: epsil, eta

    ! Local field
    integer                                 :: i, j, k      ! loop indices
    complex                                 :: ii=(0.0_WP,1.0_WP)
    logical                                 :: success
    type(COMPLEX_DATA_LAYOUT)               :: tmp_cmplx
    type(REAL_DATA_LAYOUT)                  :: tmp_real

    ! ===== Initialisation =====
    res = .false.
    if (.not. copystructonly(U,tmp_real) ) return
    if (.not. copystructonly(Uk,tmp_cmplx) ) return

    ! ===== Compute u' and d(u_i)/dx_i) direction by direction (to check isotropy) =====
    ! -- X component --
    tmp_real%values = U%values**2
    u_prime_carre(1) = computeFieldAvg(tmp_real,spec_rank,nbcpus)

    do k = tmp_cmplx%zmin,tmp_cmplx%zmax
        do j = tmp_cmplx%ymin,tmp_cmplx%ymax
            do i = tmp_cmplx%xmin,tmp_cmplx%xmax
                tmp_cmplx%values(i,j,k) = ii*WN%kx(i)*Uk%values(i,j,k)
            end do
        end do
    end do
    call btran(tmp_cmplx,tmp_real,success)
    if (.not.success) then
        write(6,*)'[ERROR] Unable to perform fft to compute du/dx on rank ', spec_rank
        return
    end if
!    if (.not. computeDerivationField(U,(/1/),tmp_real,spec_rank,WN)) then 
!        write(6,*)'[ERROR] Unable to perform fft to compute dw/dz on rank ', spec_rank
!        return
!    end if
    tmp_real%values = tmp_real%values**2
    di_u_i(1) = computeFieldAvg(tmp_real,spec_rank,nbcpus)

    ! -- Y component --
    tmp_real%values = V%values**2
    u_prime_carre(2) = computeFieldAvg(tmp_real,spec_rank,nbcpus)

    do k = tmp_cmplx%zmin,tmp_cmplx%zmax
        do j = tmp_cmplx%ymin,tmp_cmplx%ymax
            do i = tmp_cmplx%xmin,tmp_cmplx%xmax
                tmp_cmplx%values(i,j,k) = ii*WN%ky(j)*Vk%values(i,j,k)
            end do
        end do
    end do
    call btran(tmp_cmplx,tmp_real,success)
    if (.not.success) then
        write(6,*)'[ERROR] Unable to perform fft to compute dv/dy on rank ', spec_rank
        return
    end if
!    if (.not. computeDerivationField(V,(/2/),tmp_real,spec_rank,WN)) then 
!        write(6,*)'[ERROR] Unable to perform fft to compute dw/dz on rank ', spec_rank
!        return
!    end if
    tmp_real%values = tmp_real%values**2
    di_u_i(2) = computeFieldAvg(tmp_real,spec_rank,nbcpus)

    ! -- Z component --
    tmp_real%values = W%values**2
    u_prime_carre(3) = computeFieldAvg(tmp_real,spec_rank,nbcpus)

    do k = tmp_cmplx%zmin,tmp_cmplx%zmax
        do j = tmp_cmplx%ymin,tmp_cmplx%ymax
            do i = tmp_cmplx%xmin,tmp_cmplx%xmax
                tmp_cmplx%values(i,j,k) = ii*WN%kz(k)*Wk%values(i,j,k)
            end do
        end do
    end do
    call btran(tmp_cmplx,tmp_real,success)
    if (.not.success) then
        write(6,*)'[ERROR] Unable to perform fft to compute dw/dz on rank ', spec_rank
        return
    end if
!    if (.not. computeDerivationField(W,(/3/),tmp_real,spec_rank,WN)) then
!        write(6,*)'[ERROR] Unable to perform fft to compute dw/dz on rank ', spec_rank
!        return
!    end if
    tmp_real%values = tmp_real%values**2
    di_u_i(3) = computeFieldAvg(tmp_real,spec_rank,nbcpus)

    ! ===== Compute kinetic energy, lambda and Rlambda =====
    E_kine = sum(u_prime_carre)/3.0_WP
    !lambda = sqrt(2.0_WP*E_kine/(sum(di_u_i)/3.0_WP))
    lambda = sqrt(E_kine/(sum(di_u_i)/3.0_WP))
lambda = sqrt(u_prime_carre(1)/di_u_i(1))
    R_lambda= sqrt(E_kine)*lambda/visco
    E_kine = E_kine*3.0_WP/2.0_WP

    ! ===== Compute epsilon =====
    ! epsil = 2*viscosity*1/2*(sum on i,j)<(d(u_i)/dj+d(u_j/dj))^2>
    ! we have already computed: di_u_i = < (d(u_i)/di)^2 >
    ! it remains to compute the multiplicative constant and
    ! (sum on i,j with i/=j)<(d(u_i)/dj+d(u_j/dj))^2> = 2 (sum on i,j with i<j)<(d(u_i)/dj+d(u_j/dj))^2>
    if(present(epsil)) then
        ! -- Diagonal contribution : d(u_i)/dx_i) --
        epsil = sum(di_u_i)
        ! -- Contribution with i=x, j=y --
        ! Derivation
        do k = tmp_cmplx%zmin,tmp_cmplx%zmax
            do j = tmp_cmplx%ymin,tmp_cmplx%ymax
                do i = tmp_cmplx%xmin,tmp_cmplx%xmax
                    tmp_cmplx%values(i,j,k) = ii*(WN%ky(j)*Uk%values(i,j,k) + WN%kx(i)*Vk%values(i,j,k))
                end do
            end do
        end do
        call btran(tmp_cmplx,tmp_real,success)
        if (.not.success) then
            write(6,*)'[ERROR] Unable to perform fft to compute dw/dz on rank ', spec_rank
            return
        end if
        ! RMS
        tmp_real%values = tmp_real%values**2
        epsil = epsil + 0.5_WP*computeFieldAvg(tmp_real,spec_rank,nbcpus)

        ! -- Contribution with i=x, j=z --
        ! Derivation
        do k = tmp_cmplx%zmin,tmp_cmplx%zmax
            do j = tmp_cmplx%ymin,tmp_cmplx%ymax
                do i = tmp_cmplx%xmin,tmp_cmplx%xmax
                    tmp_cmplx%values(i,j,k) = ii*(WN%kz(k)*Uk%values(i,j,k) + WN%kx(i)*Wk%values(i,j,k))
                end do
            end do
        end do
        call btran(tmp_cmplx,tmp_real,success)
        if (.not.success) then
            write(6,*)'[ERROR] Unable to perform fft to compute dw/dz on rank ', spec_rank
            return
        end if
        ! RMS
        tmp_real%values = tmp_real%values**2
        epsil = epsil + 0.5_WP*computeFieldAvg(tmp_real,spec_rank,nbcpus)

        ! -- Contribution with i=y, j=z --
        ! Derivation
        do k = tmp_cmplx%zmin,tmp_cmplx%zmax
            do j = tmp_cmplx%ymin,tmp_cmplx%ymax
                do i = tmp_cmplx%xmin,tmp_cmplx%xmax
                    tmp_cmplx%values(i,j,k) = ii*(WN%kz(k)*Vk%values(i,j,k) + WN%ky(j)*Wk%values(i,j,k))
                end do
            end do
        end do
        call btran(tmp_cmplx,tmp_real,success)
        if (.not.success) then
            write(6,*)'[ERROR] Unable to perform fft to compute dw/dz on rank ', spec_rank
            return
        end if
        ! RMS
        tmp_real%values = tmp_real%values**2
        epsil = epsil + 0.5_WP*computeFieldAvg(tmp_real,spec_rank,nbcpus)

        ! -- Just multiply by the right constant --
        epsil = 2.0_WP*visco*epsil


    ! ===== Compute eta =====
        if(present(eta)) eta = ((visco**3)/epsil)**(1.0_WP/4.0_WP)


    end if ! if present epsil

    ! ===== Free memory and delete temporary fields =====
    call deleteDataLayout(tmp_real)
    call deleteDataLayout(tmp_cmplx)

end subroutine compute_flow_stat

!------------------------------------------------------------------------------
!> Show mean kinetic and magnetic energy.
!! @author = Guillaume BALARAC, LEGI
!!    @param[in]    U           = longitudinal velocity
!!    @param[in]    V           = vertical velocity
!!    @param[in]    W           = spanwise velocity
!!    @param[in]    B           = B vector field
!!    @param[in]    spec_rank   = mpi rank inside the spectral communicator
!!    @param[in]    imhd        = logical true if mhd is present
!!    @param[in]    nbcpus      = number of process in the pool
!!    @param[in]    sim_time    = global simulation time
!------------------------------------------------------------------------------
function showEnergyAndRlambda(U,V,W,B, Uk, Vk, Wk, WN, visco, spec_rank,imhd,sim_time,nbcpus,tstep,Emgrowth) result(success)

    type(REAL_DATA_LAYOUT), intent(in)              :: U,V,W
    type(COMPLEX_DATA_LAYOUT), intent(in)           :: Uk,Vk,Wk
    type(WAVENUMBERS), intent(in)                   :: WN
    type(REAL_DATA_LAYOUT), intent(in), dimension(:):: B
    real(WP), intent(in)                            :: visco
    integer, intent(in)                             :: spec_rank
    integer, intent(in)                             :: nbcpus
    logical, intent(in)                             :: imhd
    logical                                         :: success
    logical                                         :: Emgrowth

    real(WP)                :: E_kine, E_B, sim_time,tstep
    real(WP)                :: lambda, R_lambda
    real(WP), dimension(3)  :: u_prime_carre, di_u_i
    real(WP)    :: eps = 1e-10

    real(WP):: epsil, eta
    
    success = .false.

    ! Kinetic energy and Rlambda
    call compute_flow_stat(spec_rank, nbcpus, success,  & ! Working information
                        & U,V,W, Uk, Vk,Wk, WN, visco,  & ! Input
                        & E_kine, lambda, R_lambda,     & ! basic post-process
                        & u_prime_carre, di_u_i,        & ! about isotropy
                        & epsil, eta)                     ! advanced post-process)

    E_B = 0.0_WP
    if(imhd) then
        E_B = computeEnergy(B(1),B(2),B(3),spec_rank,nbcpus)
        if (Emgrowth) then
           if (sim_time .lt.eps) then
            InitialEm = E_B
            Growth_B = 0.0_WP
           else 
            Growth_B = ((sim_time-tstep)*Growth_B +tstep*log(E_B/InitialEm))/sim_time

           endif
        endif
    end if

    if(spec_rank.EQ.0) then
        write(6,'(a,g15.8)')'    [INFO-STAT] Mean kinetic energy= ', E_kine
        if ( sim_time .lt. eps ) then
            open(10,file='kinetic_energy.out',form='formatted')
        else
            open(10,file='kinetic_energy.out',form='formatted',position='append')
        endif
        write(10,*) sim_time, E_kine
        close(10)

        IF (imhd) THEN
            write(6,'(a,g15.8)')'    [INFO-STAT] Mean magnetic energy= ', E_B
            write(6,'(a,g15.8)')'    [INFO-STAT] Mean total energy= ', E_kine + E_B
            open(10,file='magnetic_energy.out',form='formatted',position='append')
                write(10,*) sim_time, E_B
            close(10)
            open(10,file='total_energy.out',form='formatted',position='append')
                write(10,*) sim_time, E_kine + E_B
            close(10)
            if (Emgrowth) then
               open(10,file='magnetic_growth_rate.out',form='formatted',position='append')
                    write(10,*) sim_time, Growth_B
               close(10)
            endif
        END IF

        ! Rlambda
        if ( sim_time .lt. eps ) then
            open(10,file='Re_lambda.out',form='formatted')
        else
            open(10,file='Re_lambda.out',form='formatted',position='append')
        end if
             write(10,*) sim_time, R_lambda
        close(10)

        ! Sum up all statistic
        if ( sim_time .lt. eps ) then
            open(10,file='Flow_stat.out',form='formatted')
            write(10,'(12(a10,5x))')        &
                & "#t_simu", "E_kine", "lambda", "Re_lambda", "epsil", &
                & "eta", "u'^2", "v'^2", "w'^2", "<(du/dx)^2>",        &
                & "<(dv/dy)^2>", "<(dw/dz)^2>"
        else
            open(10,file='Flow_stat.out',form='formatted',position='append')
        end if
        write(10,'(12(g14.7,x))') sim_time, E_kine, lambda, R_lambda, epsil, eta, u_prime_carre, di_u_i
        close(10)

    endif

    success = .true.

end function showEnergyAndRlambda


!------------------------------------------------------------------------------
!> Compute and save spectrum of a given field (scalars, not kinetic energy
!! spectrum).
!! @author Jean-Baptiste Lagaert
!!    @param[in]    Uk, Vk, Wk      = velocity field from wich spectrum has to been computed.
!!    @param[in]    ite             = time iteration number
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in]    nbcpus          = number of processes
!!    @return       success         = logical equals to true is everything is right.
!! @details
!!        This function compute the spectrum of a given field. It perform 
!!    a FFT in order to obtain it in the Fourrier space (wich is required
!!    to computing the spectrum). It could be perform "out a simulation"
!!    (ie from a field write in an output).
!------------------------------------------------------------------------------
!function velo_spectrum_rescaled(Uk, Vk, Wk, WN, epsil, eta_k, nbcpus, visco, k_lambda, spec_rank) result(success)
function velo_spectrum_rescaled(Uk, Vk, Wk, WN, epsil, eta_k, k_lambda, spec_rank) result(success)

    use datalayout
    use fileio
    use parallel_tools
    use transforms_tools

    ! Input/Output
    type(COMPLEX_DATA_LAYOUT), intent(in)   :: Uk, Vk, WK
    type(WaveNumbers), intent(in)           :: WN               ! associated wave numbers
    real(WP), intent(in)                    :: epsil, eta_k, k_lambda
!    real(WP), intent(in)                    :: visco
    integer, intent(in)                     :: spec_rank
!    integer, intent(in)                     :: nbcpus
    logical                                 :: success
    ! Others variables
    real(WP), dimension(:), allocatable     :: WN_norm          ! to know the different wave number modules
    real(WP), dimension(:), allocatable     :: spectrum         ! to save spectrum
    integer                                 :: size_spec        ! spectrum size
    integer                                 :: i, j ,k, ik      ! indices
    !integer                                 :: nk               ! indice for maximal wave number
    real(WP)                                :: kk               ! norm of wave number
    real(WP)                                :: dk               ! space step in Fourrier space
    real(WP)                                :: kc, eps          ! to filter spectrum
    integer                                 :: file_id          ! id for file io
    character(len=50)                       :: file_name        ! name of output file
    !integer                                 :: ires             ! error code for allocation
    real(WP), allocatable, dimension(:)     :: rescaling        ! to compute rescaling 
    real(WP), allocatable, dimension(:)     :: weight           ! to apply mask for the different resaling
!    real(WP)                                :: cst_B,           ! rescaling constants
    real(WP)                                :: cst_OC    ! rescaling constants
    real(WP)                                :: max_sc           ! spectrum range associated to the different rescaling
!    real(WP)                                :: min_sc   ! spectrum range associated to the different rescaling

    ! Init
    success = .false.

    ! ===== Check if spectral space step is the same above all direction =====
    if((Uk%Lx/=Uk%Ly).or.(Uk%Lx/=Uk%Lz)) then
        write(*,'(a)') '[warning] Spectrum not yet implemented if not(Lx=Ly=Lz)'
        return
    end if


    ! ===== Initialisation =====
    success = .false.
    dk=2.0_WP*acos(-1.0_WP)/Uk%Lx
    kc = 2.0_WP*acos(-1.0_WP)*(Uk%nx-1)/Uk%Lx
    eps = dk/2.0


    ! ===== Compute spectrum =====
    size_spec = size(WN%kx)
    allocate(spectrum(size_spec))
    allocate(WN_norm(size_spec))
    WN_norm = -1
    spectrum = 0
    ! Compute
    do k = Uk%zmin,Uk%zmax
        do j =  Uk%ymin,Uk%ymax
            do i =  Uk%xmin,Uk%xmax
                ! Compute norm of wave number
                kk = WN%kx(i)*WN%kx(i) + &
                    & WN%ky(j)*WN%ky(j) + &
                    & WN%kz(k)*WN%kz(k)
                kk = sqrt(real(kk))
                ! Compute indice
                ik = nint(kk/dk) 
                if ((kk>eps).and.(kk<=kc)) then
                    ! And then update spectrum
                    spectrum(ik) = spectrum(ik) + Uk%values(i,j,k)*conjg(Uk%values(i,j,k)) &
                        & + Vk%values(i,j,k)*conjg(Vk%values(i,j,k)) &
                        & + Wk%values(i,j,k)*conjg(Wk%values(i,j,k))
                end if
            end do
        end do
    end do

    ! ===== Sum values on all processes =====
    spectrum = DoGlobalVectorSum(spec_rank, spectrum)
    do ik = 1, size_spec
        WN_norm(ik) = dk*real(ik)
    end do

    ! ===== Rescaling =====
    allocate(weight(size_spec))
    allocate(rescaling(size_spec))
    ! -- Intertial range --
    rescaling = (epsil**(2.0_WP/3.0_WP))*(WN_norm**(-5._WP/3._WP))
    ! Compute constant coefficient
    max_sc = k_lambda
    where((WN_norm<max_sc).and.(spectrum>0.0_WP))
        weight = 1.0_WP
    elsewhere 
        weight = 0.0_WP
    end where
    cst_OC = sum(weight)
    if((cst_OC<1.0_WP).and.(spec_rank==0)) print*, 'poids = 0 et k_lambda = ,', k_lambda
    where(weight>0.0_WP)
        weight = weight*spectrum/rescaling
    end where
    cst_OC = sum(weight)/cst_OC
    !if (cst_OC>0.0_WP) then
    !    rescaling = cst_OC*rescaling
    !end if

    ! ===== Write output =====
    if (spec_rank==0) then
        file_id = iopen()
        write(file_name,'(a)') 'rescaled_velo_spectrum.table'
        open(unit=file_id, file=trim(adjustl(file_name)), form="FORMATTED", ERR=300)
        write(file_id, '(3(a,e14.6),a)') '#### Rescaled velo spectrum - espilon = ',  &
            & epsil, ' eta_k = ', eta_k, 'cst_OC',cst_OC, ' ####'
        write(file_id, '(a,a,a)') &
            & '# Fields are : k // k.eta_B // E // '           ,&
            & 'chi*(epsil**(-1./3.))*(k**(-5./**3.) // ',&
            & 'cst_OC*chi*(epsil**(-1./3.))*(k**(-5./**3.) // '
        do ik = 1, size_spec
            if (WN_norm(ik)>=0) write(file_id, '(5(e14.6,x))') &
                & WN_norm(ik), WN_norm(ik)*eta_k, spectrum(ik), rescaling(ik), cst_OC*rescaling(ik)
        end do
        close(file_id)
        file_id = iclose(file_id)
    end if

    ! ===== Free memory =====
    deallocate(spectrum)
    deallocate(rescaling)
    deallocate(weight)
    deallocate(WN_norm)

    ! ===== Return sucess =====
    success = .true.
    return

    ! ===== IO error =====
300 write(6,'(a)') 'Unable to open file "'//trim(file_name)//'" for spectrum'
    return

end function velo_spectrum_rescaled

!-------------------------------------------------------------------------------------------
!> @detail
!> check the equality between gradient model and gradient model with incompressibility 
!!
!! @author Antoine Vollant, LEGI
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in]    nbcpus          = number of process in the cpu pool
!!    @param[in]    U               = Component of velocity
!!    @param[in]    V               = Component of velocity
!!    @param[in]    W               = Component of velocity
!!    @return       res             = logical value
!
!-------------------------------------------------------------------------------------------
subroutine testTransfertNrjQDM(spec_rank,res,U,V,W)

    use filtering_tools
    use datalayout
    use differential_tools
    use mpilayout_tools

    ! Input/output
    logical, intent(out)                    :: res
    integer, intent(in)                     :: spec_rank
    type(REAL_DATA_LAYOUT), intent(in)      :: U,V,W

    ! Local field
    type(WAVENUMBERS)                       :: WN
    type(REAL_DATA_LAYOUT)   :: resS,UBar,VBar,WBar
    type(REAL_DATA_LAYOUT)   :: SBarChap11,SBarChap12,SBarChap13,SBarChap22,SBarChap23,SBarChap33
    type(REAL_DATA_LAYOUT)   :: SmBarChap11,SmBarChap12,SmBarChap13,SmBarChap22,SmBarChap23,SmBarChap33
    type(REAL_DATA_LAYOUT)   :: S11,S12,S13,S22,S23,S33
    type(REAL_DATA_LAYOUT)   :: SBar11,SBar12,SBar13,SBar22,SBar23,SBar33
    type(REAL_DATA_LAYOUT)   :: SmBar11,SmBar12,SmBar13,SmBar22,SmBar23,SmBar33
    type(REAL_DATA_LAYOUT)   :: T11,T12,T13,T22,T23,T33,temp
    type(REAL_DATA_LAYOUT)   :: L11,L12,L13,L22,L23,L33
    type(REAL_DATA_LAYOUT)   :: tau11,tau12,tau13,tau22,tau23,tau33
    type(REAL_DATA_LAYOUT)   :: tauChap11,tauChap12,tauChap13,tauChap22,tauChap23,tauChap33,tau_SBar_toutChap
    type(REAL_DATA_LAYOUT)   :: tau_SBar,coef1,L_SBarChap,coef2,T_SBarChap,coef3,tauChap_SBarChap,coef4
    type(REAL_DATA_LAYOUT)   :: Unuldiv,Vnuldiv,Wnuldiv
    type(REAL_DATA_LAYOUT)   :: UBarChap,VBarChap,WBarChap,epsNu
    type(COMPLEX_DATA_LAYOUT):: UBark,VBark,WBark,UBarChapk,VBarChapk,WBarChapk
    type(COMPLEX_DATA_LAYOUT):: Uk,Vk,Wk
    integer                  :: nbcpus,filter,ii,jj,jjmin,jjmax
    real(WP)                 :: ndTestFilt,nd,AvgtauChap_SBarChap,AvgT_SBarChap,AvgL_SBarChap,avgcoef1,avgcoef2,avgcoef3,avgcoef4
    real(WP)                 :: avgtau_SBar,nu,avgcoef5,avgcoef6,avg_epsNu,Avgtau_Sbar_toutChap,avgcoef7,avgcoef8
    real(wp)                 :: SikmoinsSjkchapSijchap,SikmoinsSjkSijtoutchap  
    real(wp)                 :: minUbar,maxUbar,minVbar,maxVbar,minWbar,maxWbar
    character(len=16)        :: charNd
    logical                  :: existe,varTestFilt

    ! ===== Initialisation =====
    res = .false.
    varTestFilt = .false.
    nbcpus = getnbcpus()
    if (spec_rank .eq. 0 ) print *,'testTransfertNrjQDM'
    if (.not. samelayout(U,V,W)) return
    if (.not. initDataLayout("UFilk", UBark,(U%nx/2)+1,U%ny,U%nz, &
       & U%Lx,U%Ly,U%Lz,nbcpus,spec_rank,alongZ)) then 
      write(6,'(a,i0,a)')'[ERROR] Not able to create Exk'
      return
    endif
    CALL parser_read('Viscosity',nu)
    if ( parser_is_defined('Variation of test filter size')) then
      call parser_read('Variation of test filter size',varTestFilt)
    endif
    if (.not. copystructonly(U,Unuldiv) ) return
    if (.not. copystructonly(U,Vnuldiv) ) return
    if (.not. copystructonly(U,Wnuldiv) ) return
    if (.not. copystructonly(U,tau12) ) return
    if (.not. copystructonly(U,tau13) ) return
    if (.not. copystructonly(U,tau11) ) return
    if (.not. copystructonly(U,tau12) ) return
    if (.not. copystructonly(U,tau13) ) return
    if (.not. copystructonly(U,tau22) ) return
    if (.not. copystructonly(U,tau23) ) return
    if (.not. copystructonly(U,tau33) ) return
    if (.not. copystructonly(U,tauChap11) ) return
    if (.not. copystructonly(U,tauChap12) ) return
    if (.not. copystructonly(U,tauChap13) ) return
    if (.not. copystructonly(U,tauChap22) ) return
    if (.not. copystructonly(U,tauChap23) ) return
    if (.not. copystructonly(U,tauChap33) ) return
    if (.not. copystructonly(U,S11) ) return
    if (.not. copystructonly(U,S12) ) return
    if (.not. copystructonly(U,S13) ) return
    if (.not. copystructonly(U,S22) ) return
    if (.not. copystructonly(U,S23) ) return
    if (.not. copystructonly(U,S33) ) return
    if (.not. copystructonly(U,SBar11) ) return
    if (.not. copystructonly(U,SBar12) ) return
    if (.not. copystructonly(U,SBar13) ) return
    if (.not. copystructonly(U,SBar22) ) return
    if (.not. copystructonly(U,SBar23) ) return
    if (.not. copystructonly(U,SBar33) ) return
    if (.not. copystructonly(U,SmBar11) ) return
    if (.not. copystructonly(U,SmBar12) ) return
    if (.not. copystructonly(U,SmBar13) ) return
    if (.not. copystructonly(U,SmBar22) ) return
    if (.not. copystructonly(U,SmBar23) ) return
    if (.not. copystructonly(U,SmBar33) ) return
    if (.not. copystructonly(U,SBarChap11) ) return
    if (.not. copystructonly(U,SBarChap12) ) return
    if (.not. copystructonly(U,SBarChap13) ) return
    if (.not. copystructonly(U,SBarChap22) ) return
    if (.not. copystructonly(U,SBarChap23) ) return
    if (.not. copystructonly(U,SBarChap33) ) return
    if (.not. copystructonly(U,SmBarChap11) ) return
    if (.not. copystructonly(U,SmBarChap12) ) return
    if (.not. copystructonly(U,SmBarChap13) ) return
    if (.not. copystructonly(U,SmBarChap22) ) return
    if (.not. copystructonly(U,SmBarChap23) ) return
    if (.not. copystructonly(U,SmBarChap33) ) return
    if (.not. copystructonly(U,L11) ) return
    if (.not. copystructonly(U,L12) ) return
    if (.not. copystructonly(U,L13) ) return
    if (.not. copystructonly(U,L22) ) return
    if (.not. copystructonly(U,L23) ) return
    if (.not. copystructonly(U,L33) ) return
    if (.not. copystructonly(U,T11) ) return
    if (.not. copystructonly(U,T12) ) return
    if (.not. copystructonly(U,T13) ) return
    if (.not. copystructonly(U,T22) ) return
    if (.not. copystructonly(U,T23) ) return
    if (.not. copystructonly(U,T33) ) return
    if (.not. copystructonly(U,resS) ) return
    if (.not. copystructonly(U,UBar) ) return
    if (.not. copystructonly(U,VBar) ) return
    if (.not. copystructonly(U,WBar) ) return
    if (.not. copystructonly(U,UBarChap) ) return
    if (.not. copystructonly(U,VBarChap) ) return
    if (.not. copystructonly(U,WBarChap) ) return
    if (.not. copystructonly(UBark,VBark) ) return
    if (.not. copystructonly(UBark,WBark) ) return
    if (.not. copystructonly(UBark,VBarChapk) ) return
    if (.not. copystructonly(UBark,WBarChapk) ) return
    if (.not. copystructonly(UBark,UBarChapk) ) return
    if (.not. copystructonly(UBark,Uk) ) return
    if (.not. copystructonly(UBark,Vk) ) return
    if (.not. copystructonly(UBark,Wk) ) return
    if (.not. copystructonly(U,temp) ) return
    if (.not. copystructonly(U,tau_SBar) ) return
    if (.not. copystructonly(U,coef1) ) return
    if (.not. copystructonly(U,L_SBarChap) ) return
    if (.not. copystructonly(U,coef2) ) return
    if (.not. copystructonly(U,T_SBarChap) ) return
    if (.not. copystructonly(U,tau_SBar_toutChap) ) return
    if (.not. copystructonly(U,coef3) ) return
    if (.not. copystructonly(U,epsNu) ) return
    if (.not. copystructonly(U,tauChap_SBarChap) ) return
    if (.not. copystructonly(U,coef4) ) return
    !Calcul des nombres d'onde
    if (.not. initWN(WN,spec_rank,U%nx,U%ny,U%nz) .or. &
       &.not.  computeWN(WN,spec_rank,U%Lx,U%Ly,U%Lz,U%nx,U%ny,U%nz) ) then 
      write(6,'(a,i0,a)')'[ERROR] Cannot compute WN'
      return
    endif
    
    !cutoff -> 1
    filter=1
    ndTestFilt = 2.0_WP
    if (.not. diff_tool_initBuff(UBark)) return
    call ftran (U,Uk,res)
    call ftran (V,Vk,res)
    call ftran (W,Wk,res)
    IF (.NOT. showDivU(Uk,Vk,Wk,WN,spec_rank)) return
    IF (.NOT. NullifyDiv(Uk,Vk,WK,WN,spec_rank)) return 
    IF (.NOT. showDivU(Uk,Vk,Wk,WN,spec_rank)) return
    call btran (Uk,Unuldiv,res)
    call btran (Vk,Vnuldiv,res)
    call btran (Wk,Wnuldiv,res)
    if (varTestFilt) then
      jjmin = 15
      jjmax = 45
    else
      jjmin = 1
      jjmax = 1
    endif
    !do ii=2,U%nx/16,2
    ii=1
      nd=real(ii,WP)
      do jj=jjmin,jjmax
        if (varTestFilt) then
          ndTestFilt = real(jj,WP)/10.0_WP
          write (charNd,'(f4.1,a,f4.1)') nd,'_',ndTestFilt
          if (spec_rank .eq. 0) print *,'[INFO testTransfertNrjQDM] nd =',nd,'ndTestFilt =',ndTestFilt
        else
          if (spec_rank .eq. 0) print *,'[INFO testTransfertNrjQDM] nd =',nd
          write (charNd,'(f4.1)') nd
        endif
        call computeFilter(WN,computeDelta(U)*nd,Uk,UBark,filter)
        call computeFilter(WN,computeDelta(U)*nd,Vk,VBark,filter)
        call computeFilter(WN,computeDelta(U)*nd,Wk,WBark,filter)
!        UBark%values = Uk%values
!        VBark%values = Vk%values
!        WBark%values = Wk%values
    call btran (UBark,UBarChap,res)
    if (.not. equalToVal(UBarChap,Unuldiv)) print *,'UBarChap != U'
        call computeFilter(WN,ndTestFilt*computeDelta(U)*nd,Uk,UBarChapk,filter)
        call computeFilter(WN,ndTestFilt*computeDelta(U)*nd,Vk,VBarChapk,filter)
        call computeFilter(WN,ndTestFilt*computeDelta(U)*nd,Wk,WBarChapk,filter)
!        call btran(UBark,UBar,res)
!        call btran(VBark,VBar,res)
!        call btran(WBark,WBar,res)
        UBar%values = Unuldiv%values
        VBar%values = Vnuldiv%values
        WBar%values = Wnuldiv%values
        if (.not. equalToVal(Unuldiv,UBar)) print *,' U != Ubar'
        UBar%name="UBar"
        VBar%name="VBar"
        WBar%name="WBar"
        if (.not. computeFieldPDF(int(nd),UBar,300,spec_rank,10000)) return
        if (.not. computeFieldPDF(int(nd),VBar,300,spec_rank,10000)) return
        if (.not. computeFieldPDF(int(nd),WBar,300,spec_rank,10000)) return
        call btran(UBarChapk,UBarChap,res)
        call btran(VBarChapk,VBarChap,res)
        call btran(WBarChapk,WBarChap,res)
        UBarChap%name="UBarChap"
        VBarChap%name="VBarChap"
        WBarChap%name="WBarChap"
        if (.not. computeFieldPDF(int(nd),UBarChap,300,spec_rank,10000)) return
        if (.not. computeFieldPDF(int(nd),VBarChap,300,spec_rank,10000)) return
        if (.not. computeFieldPDF(int(nd),WBarChap,300,spec_rank,10000)) return
  
        !Sij Bar
        call computeStressTensor(UBark,VBark,WBark,WN,SBar11,SBar12,SBar13,SBar22,SBar23,SBar33,res)
        !Sij Bar Chapeau
        call computeStressTensor(UBarchapk,VBarchapk,WBarchapk,WN,SBarChap11,SBarChap12,SBarChap13,SBarChap22,SBarChap23,SBarChap33,res)
        !Tau et Tau Chapeau
        call computeT_ijVecA(Unuldiv,Vnuldiv,Wnuldiv,UBar,VBar,WBar,WN,tau11,tau12,tau13,tau22,tau23,tau33,filter,res,computeDelta(U)*nd)
        call computeFilter(WN,ndTestFilt*computeDelta(U)*nd,tau11,tauChap11,filter)
        call computeFilter(WN,ndTestFilt*computeDelta(U)*nd,tau12,tauChap12,filter)
        call computeFilter(WN,ndTestFilt*computeDelta(U)*nd,tau13,tauChap13,filter)
        call computeFilter(WN,ndTestFilt*computeDelta(U)*nd,tau22,tauChap22,filter)
        call computeFilter(WN,ndTestFilt*computeDelta(U)*nd,tau23,tauChap23,filter)
        call computeFilter(WN,ndTestFilt*computeDelta(U)*nd,tau33,tauChap33,filter)
        
        !T
        T11%values=UBar%values ** 2
        call computeFilter(WN,ndTestFilt*computeDelta(U)*nd,T11,temp,filter)
        T11%values=temp%values - UBarChap%values ** 2
        T22%values=VBar%values ** 2
        call computeFilter(WN,ndTestFilt*computeDelta(U)*nd,T22,temp,filter)
        T22%values=temp%values - VBarChap%values ** 2 
        T33%values=WBar%values ** 2
        call computeFilter(WN,ndTestFilt*computeDelta(U)*nd,T33,temp,filter)
        T33%values=temp%values - WBarChap%values ** 2
  
        T12%values=UBar%values * VBar%values
        call computeFilter(WN,ndTestFilt*computeDelta(U)*nd,T12,temp,filter)
        T12%values=temp%values - UBarChap%values * VBarChap%values
        T13%values=UBar%values * WBar%values
        call computeFilter(WN,ndTestFilt*computeDelta(U)*nd,T13,temp,filter)
        T13%values=temp%values - UBarChap%values * WBarChap%values
        T23%values=VBar%values * WBar%values
        call computeFilter(WN,ndTestFilt*computeDelta(U)*nd,T23,temp,filter)
        T23%values=temp%values - VBarChap%values * WBarChap%values
        !les procedures sont equivalentes
        !call computeT_ijVecA(UBar,VBar,WBar,UBarChap,VBarChap,WBarChap,WN,T11,T12,T13,T22,T23,T33,filter,res,ndTestFilt*computeDelta(U)*nd)
        
        !L
        L11%values=Unuldiv%values ** 2
        call computeFilter(WN,computeDelta(U)*nd,L11,temp,filter)
        call computeFilter(WN,ndTestFilt*computeDelta(U)*nd,temp,L11,filter)
        L11%values=L11%values - UBarChap%values ** 2
        L22%values=Vnuldiv%values ** 2
        call computeFilter(WN,computeDelta(U)*nd,L22,temp,filter)
        call computeFilter(WN,ndTestFilt*computeDelta(U)*nd,temp,L22,filter)
        L22%values=L22%values - VBarChap%values ** 2 
        L33%values=Wnuldiv%values ** 2
        call computeFilter(WN,computeDelta(U)*nd,L33,temp,filter)
        call computeFilter(WN,ndTestFilt*computeDelta(U)*nd,temp,L33,filter)
        L33%values=L33%values - WBarChap%values ** 2
  
        L12%values=Unuldiv%values * Vnuldiv%values
        call computeFilter(WN,computeDelta(U)*nd,L12,temp,filter)
        call computeFilter(WN,ndTestFilt*computeDelta(U)*nd,temp,L12,filter)
        L12%values=L12%values - UBarChap%values * VBarChap%values
        L13%values=Unuldiv%values * Wnuldiv%values
        call computeFilter(WN,computeDelta(U)*nd,L13,temp,filter)
        call computeFilter(WN,ndTestFilt*computeDelta(U)*nd,temp,L13,filter)
        L13%values=L13%values - UBarChap%values * WBarChap%values
        L23%values=Vnuldiv%values * Wnuldiv%values
        call computeFilter(WN,computeDelta(U)*nd,L23,temp,filter)
        call computeFilter(WN,ndTestFilt*computeDelta(U)*nd,temp,L23,filter)
        L23%values=L23%values - VBarChap%values * WBarChap%values
  
        !Compute Transfert terms 
        ! tau_SBar         = tauij SijBar         = [(Ui Uj)Bar - UiBar UjBar] SijBar
        ! tauChap_SBarChap = tauijChap SijBarChap = [(Ui Uj)BarChap - (UiBar UjBar)Chap] SijBarChap
        ! T_SBarChap       = TijChap SijBarChap   = [(UiBar UjBar)Chap - UiBarChap UjBarChap] SijBarChap
        ! L_SBarChap       = LijChap SijBarChap   = [(Ui Uj)BarChap - UiBarChap UjBarChap] SijBarChap
        tau_SBar%values = tau11%values * SBar11%values + tau22%values * SBar22%values + tau33%values * SBar33%values +&
                               & 2.0_WP * (tau12%values * SBar12%values + tau13%values * SBar13%values + tau23%values * SBar23%values )
        call computeFilter(WN,ndTestFilt*computeDelta(U)*nd,tau_SBar,tau_SBar_toutChap,filter)
        tauChap_SBarChap%values = tauChap11%values * SBarChap11%values + tauChap22%values * SBarChap22%values + tauChap33%values * SBarChap33%values +&
                               & 2.0_WP * (tauChap12%values * SBarChap12%values + tauChap13%values * SBarChap13%values + tauChap23%values * SBarChap23%values )
        T_SBarChap%values = T11%values * SBarChap11%values + T22%values * SBarChap22%values + T33%values * SBarChap33%values +&
                               & 2.0_WP * (T12%values * SBarChap12%values + T13%values * SBarChap13%values + T23%values * SBarChap23%values )
        L_SBarChap%values = L11%values * SBarChap11%values + L22%values * SBarChap22%values + L33%values * SBarChap33%values +&
                               & 2.0_WP * (L12%values * SBarChap12%values + L13%values * SBarChap13%values + L23%values * SBarChap23%values )
        tauChap_SBarChap%name='tauChap_SBarChap_'//trim(adjustl(charNd))
        T_SBarChap%name='T_SBarChap_'//trim(adjustl(charNd))
        L_SBarChap%name='L_SBarChap_'//trim(adjustl(charNd))
        AvgtauChap_SBarChap = computeFieldAvg (tauChap_SBarChap,spec_rank)
        AvgT_SBarChap= computeFieldAvg (T_SBarChap,spec_rank)
        AvgL_SBarChap= computeFieldAvg (L_SBarChap,spec_rank)
        AvgTau_SBar = computeFieldAvg (tau_SBar,spec_rank)
        Avgtau_Sbar_toutChap = computeFieldAvg (tau_Sbar_toutChap,spec_rank)
  
        !U,V,W -> Ubar,Vbar,Wbar -> S-ij Bar
        if (.not. computeSplitingSij(SBar11,SBar12,SBar13,SBar22,SBar23,SBar33,SmBar11,SmBar12,SmBar13,SmBar22,SmBar23,SmBar33,'moins')) return
        if (.not. computeSplitingSij(SBarChap11,SBarChap12,SBarChap13,SBarChap22,SBarChap23,SBarChap33,&
                                   & SmBarChap11,SmBarChap12,SmBarChap13,SmBarChap22,SmBarChap23,SmBarChap33,'moins')) return
  
        !(Sik-Bar SjkBar)Chap -> Lij
!        temp%values = SmBar11%values*SBar11%values+SmBar12%values*SBar12%values+SmBar13%values*SBar13%values
!        call computeFilter(WN,ndTestFilt*computeDelta(U)*nd,temp,L11,filter)
!        temp%values = SmBar12%values*SBar12%values+SmBar22%values*SBar22%values+SmBar23%values*SBar23%values
!        call computeFilter(WN,ndTestFilt*computeDelta(U)*nd,temp,L22,filter)
!        temp%values = SmBar13%values*SBar13%values+SmBar23%values*SBar23%values+SmBar33%values*SBar33%values
!        call computeFilter(WN,ndTestFilt*computeDelta(U)*nd,temp,L33,filter)
!        temp%values = SmBar11%values*SBar12%values+SmBar12%values*SBar22%values+SmBar13%values*SBar23%values
!        call computeFilter(WN,ndTestFilt*computeDelta(U)*nd,temp,L12,filter)
!        temp%values = SmBar11%values*SBar13%values+SmBar12%values*SBar23%values+SmBar13%values*SBar33%values
!        call computeFilter(WN,ndTestFilt*computeDelta(U)*nd,temp,L13,filter)
!        temp%values = SmBar12%values*SBar13%values+SmBar22%values*SBar23%values+SmBar23%values*SBar33%values
!        call computeFilter(WN,ndTestFilt*computeDelta(U)*nd,temp,L23,filter)
!MODIF
        call RGMSikmoinsSkj_tij(tauChap11,tauChap12,tauChap13,tauChap22,tauChap23,tauChap33,Ubar,Ubark,Vbark,Wbark,WN,spec_rank,1.0_WP,1.0_WP)
!        call computeFilter(WN,2.0_WP*computeDelta(U),tauChap11,L11,filter)
!        call computeFilter(WN,2.0_WP*computeDelta(U),tauChap22,L22,filter)
!        call computeFilter(WN,2.0_WP*computeDelta(U),tauChap33,L33,filter)
!        call computeFilter(WN,2.0_WP*computeDelta(U),tauChap12,L12,filter)
!        call computeFilter(WN,2.0_WP*computeDelta(U),tauChap13,L13,filter)
!        call computeFilter(WN,2.0_WP*computeDelta(U),tauChap23,L23,filter)
        call computeFilter(WN,ndTestFilt*computeDelta(U)*nd,tauChap11,L11,filter)
        call computeFilter(WN,ndTestFilt*computeDelta(U)*nd,tauChap22,L22,filter)
        call computeFilter(WN,ndTestFilt*computeDelta(U)*nd,tauChap33,L33,filter)
        call computeFilter(WN,ndTestFilt*computeDelta(U)*nd,tauChap12,L12,filter)
        call computeFilter(WN,ndTestFilt*computeDelta(U)*nd,tauChap13,L13,filter)
        call computeFilter(WN,ndTestFilt*computeDelta(U)*nd,tauChap23,L23,filter)
        !Sik-BarChap SjkBarChap -> Tij
        T11%values = SmBarChap11%values*SBarChap11%values+SmBarChap12%values*SBarChap12%values+SmBarChap13%values*SBarChap13%values
        T22%values = SmBarChap12%values*SBarChap12%values+SmBarChap22%values*SBarChap22%values+SmBarChap23%values*SBarChap23%values
        T33%values = SmBarChap13%values*SBarChap13%values+SmBarChap23%values*SBarChap23%values+SmBarChap33%values*SBarChap33%values
        T12%values = SmBarChap11%values*SBarChap12%values+SmBarChap12%values*SBarChap22%values+SmBarChap13%values*SBarChap23%values
        T13%values = SmBarChap11%values*SBarChap13%values+SmBarChap12%values*SBarChap23%values+SmBarChap13%values*SBarChap33%values
        T23%values = SmBarChap12%values*SBarChap13%values+SmBarChap22%values*SBarChap23%values+SmBarChap23%values*SBarChap33%values
        !Sik-Bar SjkBar -> tauChapij
!MODIF
!        tauChap11%values = SmBar11%values*SBar11%values+SmBar12%values*SBar12%values+SmBar13%values*SBar13%values
!        tauChap22%values = SmBar12%values*SBar12%values+SmBar22%values*SBar22%values+SmBar23%values*SBar23%values
!        tauChap33%values = SmBar13%values*SBar13%values+SmBar23%values*SBar23%values+SmBar33%values*SBar33%values
!        tauChap12%values = SmBar11%values*SBar12%values+SmBar12%values*SBar22%values+SmBar13%values*SBar23%values
!        tauChap13%values = SmBar11%values*SBar13%values+SmBar12%values*SBar23%values+SmBar13%values*SBar33%values
!        tauChap23%values = SmBar12%values*SBar13%values+SmBar22%values*SBar23%values+SmBar23%values*SBar33%values
  
        !C1 = <tauij SijBar> / <DeltaBar**2 S-ikBar SjkBar SijBar>
        !C2 = <Lij SijBariChap> / <DeltaBarChap**2 S-ikBarChap SjkBarChap SijBarChap>
        !C3 = <Tij SijBarChap> / < (DeltaBarChap**2 S-ikBarChap SjkBarChap - DeltaBar**2 S-ikBar SjkBar ) SijBarChap >
        !C4 = <tauijChap SijBarChap> / <DeltaBar**2 (S-ikBar SjkBar)Chap SijBarChap>
        !C5 = <eps_nu + Tij SijBarChap> / < DeltaBar**2 [ ( S-ikBar SjkBar SijBar )Chap - (S-ikBar SjkBar)Chap SijBarChap) >
        !C6 = <(tauij SijBar)Chap> / <DeltaBar**2 (S-ikBar SjkBar SijBar )Chap>
        !C7 = < (Lij Hij>/<Hij Hij> avec Lij=(u|_i u|_j)^ - u|^_i u|^_j et Hij = delta|**2/12 S|^moinsik S|^jk 
        !C8 = < (Lij Hij>/<Hij Hij> avec Lij=(u|_i u|_j)^ - u|^_i u|^_j et Hij = delta|^**2/12 S|^moinsik S|^jk - delta|**2/12 (S|moinsik S|kj)^

        !C7
        !    computeT_ijVecA(in1 ,in2 ,in3 ,in1f,in2f,in3f            ,wave,out11,out12,out13,out22,out23,out33,filter,success,delta)
        call computeT_ijVecA(Ubar,Vbar,Wbar,UbarChap,VbarChap,WbarChap,WN  ,S11  ,S12  ,S13  ,S22  ,S23  ,S33  ,filter,res    ,ndTestFilt*computeDelta(U)*nd)
        coef1%values = (computeDelta(U)*nd)**2 *( S11%values * T11%values + S22%values * T22%values + S33%values * T33%values + & 
                                         2.0_WP*( S13%values * T13%values + S12%values * T12%values + S23%values * T23%values ))
        coef2%values = (computeDelta(U)*nd)**4 *( T11%values * T11%values + T22%values * T22%values + T33%values * T33%values + & 
                                         2.0_WP*( T13%values * T13%values + T12%values * T12%values + T23%values * T23%values ))
        Avgcoef7 = computeFieldAvg(coef1,spec_rank) / computeFieldAvg(coef2,spec_rank)
        !C8
        coef1%values = S11%values * ( (ndTestFilt*computeDelta(U)*nd)**2*T11%values - (computeDelta(U)*nd)**2 * L11%values ) + &
                     & S22%values * ( (ndTestFilt*computeDelta(U)*nd)**2*T22%values - (computeDelta(U)*nd)**2 * L22%values ) + &
                     & S33%values * ( (ndTestFilt*computeDelta(U)*nd)**2*T33%values - (computeDelta(U)*nd)**2 * L33%values ) + &
             2.0_WP*(  S12%values * ( (ndTestFilt*computeDelta(U)*nd)**2*T12%values - (computeDelta(U)*nd)**2 * L12%values ) + &
                     & S13%values * ( (ndTestFilt*computeDelta(U)*nd)**2*T13%values - (computeDelta(U)*nd)**2 * L13%values ) + &
                     & S23%values * ( (ndTestFilt*computeDelta(U)*nd)**2*T23%values - (computeDelta(U)*nd)**2 * L23%values ) )
        coef2%values = ( (ndTestFilt*computeDelta(U)*nd)**2*T11%values - (computeDelta(U)*nd)**2 * L11%values )**2 + &
                     & ( (ndTestFilt*computeDelta(U)*nd)**2*T22%values - (computeDelta(U)*nd)**2 * L22%values )**2 + &
                     & ( (ndTestFilt*computeDelta(U)*nd)**2*T33%values - (computeDelta(U)*nd)**2 * L33%values )**2 + &
             2.0_WP*(  ( (ndTestFilt*computeDelta(U)*nd)**2*T12%values - (computeDelta(U)*nd)**2 * L12%values )**2 + &
                     & ( (ndTestFilt*computeDelta(U)*nd)**2*T13%values - (computeDelta(U)*nd)**2 * L13%values )**2 + &
                     & ( (ndTestFilt*computeDelta(U)*nd)**2*T23%values - (computeDelta(U)*nd)**2 * L23%values )**2 )
        Avgcoef8 = computeFieldAvg(coef1,spec_rank) / computeFieldAvg(coef2,spec_rank) 

        coef1%values = (computeDelta(U)*nd)**2*( tauChap11%values * SBar11%values + tauChap22%values * SBar22%values + tauChap33%values * SBar33%values +&
                                        2.0_WP*( tauChap12%values * SBar12%values + tauChap13%values * SBar13%values + tauChap23%values * SBar23%values ))
        coef2%values = (ndTestFilt*computeDelta(U)*nd)**2*( T11%values * SBarChap11%values + T22%values * SBarChap22%values + T33%values * SBarChap33%values +&
                                                   2.0_WP*( T12%values * SBarChap12%values + T13%values * SBarChap13%values + T23%values * SBarChap23%values) )
        coef3%values =  ((ndTestFilt*computeDelta(U)*nd)**2*T11%values - (computeDelta(U)*nd)**2*L11%values) * SBarChap11%values +&
                        ((ndTestFilt*computeDelta(U)*nd)**2*T22%values - (computeDelta(U)*nd)**2*L22%values) * SBarChap22%values +&
                        ((ndTestFilt*computeDelta(U)*nd)**2*T33%values - (computeDelta(U)*nd)**2*L33%values) * SBarChap33%values +&
                 2.0_WP*((ndTestFilt*computeDelta(U)*nd)**2*T12%values - (computeDelta(U)*nd)**2*L12%values) * SBarChap12%values +&
                 2.0_WP*((ndTestFilt*computeDelta(U)*nd)**2*T13%values - (computeDelta(U)*nd)**2*L13%values) * SBarChap13%values +&
                 2.0_WP*((ndTestFilt*computeDelta(U)*nd)**2*T23%values - (computeDelta(U)*nd)**2*L23%values) * SBarChap23%values
        coef4%values =(computeDelta(U)*nd)**2 * (L11%values*SBarChap11%values + L22%values*SBarChap22%values + L33%values*SBarChap33%values + &
                                     & 2.0_WP * (L12%values*SBarChap12%values + L13%values*SBarChap13%values + L23%values*SBarChap23%values))
  
        avgCoef1 = computeFieldAvg (coef1,spec_rank)
        avgCoef2 = computeFieldAvg (coef2,spec_rank)
        avgCoef3 = computeFieldAvg (coef3,spec_rank)
        avgCoef4 = computeFieldAvg (coef4,spec_rank)
        avgCoef1 = avgtau_SBar/avgCoef1
        avgCoef2 = avgL_SBarChap / avgCoef2 
        avgCoef3 = avgT_SBarChap / avgCoef3
        avgCoef4 = avgtauChap_SBarChap / avgCoef4

       ! avg_epsNu = nu * computeFieldAvg(epsNu,spec_rank)
        call computeViscousDissResolFiltered(epsNu,UBar,UBark,VBark,WBark,WN,nu,filter,ndTestFilt*computeDelta(U)*nd,res)
        avg_epsNu = computeFieldAvg(epsNu,spec_rank)

        !C6
        coef3%values =   tauChap11%values * SBar11%values + tauChap22%values * SBar22%values + tauChap33%values * SBar33%values +& 
                2.0_WP*( tauChap12%values * SBar12%values + tauChap13%values * SBar13%values + tauChap23%values * SBar23%values )
!MODIF
        call computeFilter(WN,ndTestFilt*computeDelta(U)*nd,coef3,coef1,filter)
!        call computeFilter(WN,2.0_WP*computeDelta(U),coef3,coef1,filter)
        avgCoef6 = (computeDelta(U)*nd)**2*computeFieldAvg (coef1,spec_rank)
        avgCoef6 = avgtau_Sbar_toutChap / avgCoef6
        
        !C5
        coef1%values = (computeDelta(U)*nd)**2*coef1%values
        SikmoinsSjkchapSijchap = computeFieldAvg(coef4,spec_rank)
        SikmoinsSjkSijtoutchap = computeFieldAvg(coef1,spec_rank)
        coef3%values = coef1%values - coef4%values 
        avgCoef5 = computeFieldAvg (coef3,spec_rank)
        if (spec_rank .eq. 0) print *,'[avg_epsNu]                         ',avg_epsNu  
        if (spec_rank .eq. 0) print *,'[avgT_SBarChap]                     ',avgT_SBarChap
        if (spec_rank .eq. 0) print *,'[(S-S)^S^]                          ',SikmoinsSjkchapSijchap
        if (spec_rank .eq. 0) print *,'[(S-SS)^]                           ',SikmoinsSjkSijtoutchap
        if (spec_rank .eq. 0) print *,'[denom]                             ',avgCoef5
        avgCoef5 = ( avg_epsNu + avgT_SBarChap ) / avgCoef5
        if (spec_rank .eq. 0) print *,'[coef]                              ',avgCoef5

  
        if (spec_rank .eq. 0) then
          inquire(file='TransfertNrjQDM_coef.table', exist=existe)
          open(unit=50,file='TransfertNrjQDM_coef.table', form="FORMATTED",position="append")
          if ( .not. existe) then
            write(50,'(a)') '#nd ndTestFilt AvgtauChap_SBarChap AvgT_SBarChap AvgL_SBarChap Avg_epsNu Avgtau_Sbar_toutChap avgcoef1 avgcoef2 avgcoef3 avgcoef4 avgCoef5 avgCoef6&
                            & avgCoef7 avgCoef8'
          endif
          write(50,'(f4.1,1x,f4.1,1x,13(g15.8,1x))') nd,ndTestFilt,AvgtauChap_SBarChap,AvgT_SBarChap,AvgL_SBarChap,avg_epsNu,avgtau_Sbar_toutChap,&
                                                    &avgcoef1,avgcoef2,avgcoef3,avgcoef4,avgCoef5,avgCoef6,avgcoef7,avgcoef8
          close(50)
        endif
      enddo
      if (varTestFilt) then
        if (spec_rank .eq. 0) then 
          open(unit=50,file='TransfertNrjQDM_coef.table', form="FORMATTED",position="append")
          write (50,'(a)') ''
          close(50)
        endif
      endif
!    enddo

    ! ===== Free memory and delete temporary fields =====
    if( .not. deleteWN(WN)) return
    call deleteDataLayout(UBark)
    call deleteDataLayout(tau11)
    call deleteDataLayout(tau12)
    call deleteDataLayout(tau13)
    call deleteDataLayout(tau22)
    call deleteDataLayout(tau23)
    call deleteDataLayout(tau33)
    call deleteDataLayout(tauChap11)
    call deleteDataLayout(tauChap12)
    call deleteDataLayout(tauChap13)
    call deleteDataLayout(tauChap22)
    call deleteDataLayout(tauChap23)
    call deleteDataLayout(tauChap33)
    call deleteDataLayout(SBar11)
    call deleteDataLayout(SBar12)
    call deleteDataLayout(SBar13)
    call deleteDataLayout(SBar22)
    call deleteDataLayout(SBar23)
    call deleteDataLayout(SBar33)
    call deleteDataLayout(SmBar11)
    call deleteDataLayout(SmBar12)
    call deleteDataLayout(SmBar13)
    call deleteDataLayout(SmBar22)
    call deleteDataLayout(SmBar23)
    call deleteDataLayout(SmBar33)
    call deleteDataLayout(SBarChap11)
    call deleteDataLayout(SBarChap12)
    call deleteDataLayout(SBarChap13)
    call deleteDataLayout(SBarChap22)
    call deleteDataLayout(SBarChap23)
    call deleteDataLayout(SBarChap33)
    call deleteDataLayout(SmBarChap11)
    call deleteDataLayout(SmBarChap12)
    call deleteDataLayout(SmBarChap13)
    call deleteDataLayout(SmBarChap22)
    call deleteDataLayout(SmBarChap23)
    call deleteDataLayout(SmBarChap33)
    call deleteDataLayout(L11)
    call deleteDataLayout(L12)
    call deleteDataLayout(L13)
    call deleteDataLayout(L22)
    call deleteDataLayout(L23)
    call deleteDataLayout(L33)
    call deleteDataLayout(S11)
    call deleteDataLayout(S12)
    call deleteDataLayout(S13)
    call deleteDataLayout(S22)
    call deleteDataLayout(S23)
    call deleteDataLayout(S33)
    call deleteDataLayout(T11)
    call deleteDataLayout(T12)
    call deleteDataLayout(T13)
    call deleteDataLayout(T22)
    call deleteDataLayout(T23)
    call deleteDataLayout(T33)
    call deleteDataLayout(resS)
    call deleteDataLayout(Uk)
    call deleteDataLayout(Vk)
    call deleteDataLayout(Wk)
    call deleteDataLayout(UBar)
    call deleteDataLayout(VBar)
    call deleteDataLayout(WBar)
    call deleteDataLayout(Unuldiv)
    call deleteDataLayout(Vnuldiv)
    call deleteDataLayout(Wnuldiv)
    call deleteDataLayout(UBarChap)
    call deleteDataLayout(VBarChap)
    call deleteDataLayout(WBarChap)
    call deleteDataLayout(VBarChapk)
    call deleteDataLayout(WBarChapk)
    call deleteDataLayout(UBarChapk)
    call deleteDataLayout(temp)
    call deleteDataLayout(tau_SBar)
    call deleteDataLayout(coef1)
    call deleteDataLayout(L_SBarChap)
    call deleteDataLayout(coef2)
    call deleteDataLayout(T_SBarChap)
    call deleteDataLayout(coef3)
    call deleteDataLayout(tauChap_SBarChap)
    call deleteDataLayout(coef4)
    call deleteDataLayout(epsNu)

    res = .true.

end subroutine testTransfertNrjQDM

!-------------------------------------------------------------------------------------------
!> @detail
!> check the equality between gradient model and gradient model with incompressibility 
!!
!! @author Antoine Vollant, LEGI
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in]    nbcpus          = number of process in the cpu pool
!!    @param[in]    U               = Component of velocity
!!    @param[in]    V               = Component of velocity
!!    @param[in]    W               = Component of velocity
!!    @return       res             = logical value
!
!-------------------------------------------------------------------------------------------
subroutine testDissipationduidxkdujdxk(spec_rank,res,U,V,W)

    use filtering_tools
    use datalayout
    use differential_tools
    use mpilayout_tools

    ! Input/output
    logical, intent(out)                    :: res
    integer, intent(in)                     :: spec_rank
    type(REAL_DATA_LAYOUT), intent(inout)   :: U,V,W

    ! Local field
    type(WAVENUMBERS)                       :: WN
    type(REAL_DATA_LAYOUT)   :: resS,Ufil,Vfil,Wfil,nrjDiss
    type(REAL_DATA_LAYOUT),allocatable,dimension(:,:)   :: Sm,Spp
    type(REAL_DATA_LAYOUT),allocatable,dimension(:,:)   :: O,S
    type(COMPLEX_DATA_LAYOUT)               :: Ufilk,Vfilk,Wfilk,T1k,T2k,T3k
    integer                                 :: i,j,k,l,m,n,nbcpus,nd,filter
    character(len=16)                       ::charNd

    ! ===== Initialisation =====
    res = .false.
    nbcpus = getnbcpus()
    if (spec_rank .eq. 0 ) print *,'testsikmoinskjsij'
    if (.not. initDataLayout("UFilk", UFilk,(U%nx/2)+1,U%ny,U%nz, &
       & U%Lx,U%Ly,U%Lz,nbcpus,spec_rank,alongZ)) then 
      write(6,'(a,i0,a)')'[ERROR] Not able to create Exk'
      return
    endif
    allocate(S(3,3))
    allocate(O(3,3))
    allocate(Sm(3,3))
    allocate(Spp(3,3))
    if (.not. copystructonly(U,S(1,1)) ) return
    if (.not. copystructonly(U,S(1,2)) ) return
    if (.not. copystructonly(U,S(1,3)) ) return
    if (.not. copystructonly(U,S(2,1)) ) return
    if (.not. copystructonly(U,S(2,2)) ) return
    if (.not. copystructonly(U,S(2,3)) ) return
    if (.not. copystructonly(U,S(3,1)) ) return
    if (.not. copystructonly(U,S(3,2)) ) return
    if (.not. copystructonly(U,S(3,3)) ) return
    if (.not. copystructonly(U,Spp(1,1)) ) return
    if (.not. copystructonly(U,Spp(1,2)) ) return
    if (.not. copystructonly(U,Spp(1,3)) ) return
    if (.not. copystructonly(U,Spp(2,1)) ) return
    if (.not. copystructonly(U,Spp(2,2)) ) return
    if (.not. copystructonly(U,Spp(2,3)) ) return
    if (.not. copystructonly(U,Spp(3,1)) ) return
    if (.not. copystructonly(U,Spp(3,2)) ) return
    if (.not. copystructonly(U,Spp(3,3)) ) return
    if (.not. copystructonly(U,Sm(1,1)) ) return
    if (.not. copystructonly(U,Sm(1,2)) ) return
    if (.not. copystructonly(U,Sm(1,3)) ) return
    if (.not. copystructonly(U,Sm(2,1)) ) return
    if (.not. copystructonly(U,Sm(2,2)) ) return
    if (.not. copystructonly(U,Sm(2,3)) ) return
    if (.not. copystructonly(U,Sm(3,1)) ) return
    if (.not. copystructonly(U,Sm(3,2)) ) return
    if (.not. copystructonly(U,Sm(3,3)) ) return
    if (.not. copystructonly(U,O(1,1)) ) return
    if (.not. copystructonly(U,O(2,2)) ) return
    if (.not. copystructonly(U,O(3,3)) ) return
    if (.not. copystructonly(U,O(1,2)) ) return
    if (.not. copystructonly(U,O(1,3)) ) return
    if (.not. copystructonly(U,O(2,1)) ) return
    if (.not. copystructonly(U,O(2,3)) ) return
    if (.not. copystructonly(U,O(3,1)) ) return
    if (.not. copystructonly(U,O(3,2)) ) return
    if (.not. copystructonly(U,resS) ) return
    if (.not. copystructonly(U,nrjDiss) ) return
    if (.not. copystructonly(U,Ufil) ) return
    if (.not. copystructonly(U,Vfil) ) return
    if (.not. copystructonly(U,Wfil) ) return
    if (.not. copystructonly(Ufilk,Vfilk) ) return
    if (.not. copystructonly(Ufilk,Wfilk) ) return
    if (.not. copystructonly(Ufilk,T1k) ) return
    if (.not. copystructonly(Ufilk,T2k) ) return
    if (.not. copystructonly(Ufilk,T3k) ) return
    !Calcul des nombres d'onde
    if (.not. initWN(WN,spec_rank,U%nx,U%ny,U%nz) .or. &
       &.not.  computeWN(WN,spec_rank,U%Lx,U%Ly,U%Lz,U%nx,U%ny,U%nz) ) then 
      write(6,'(a,i0,a)')'[ERROR] Cannot compute WN'
      return
    endif
    
    !cutoff -> 1
    filter=1
    if (.not. diff_tool_initBuff(Ufilk)) return
    do nd=2,U%nx/8,2
      write (charNd,'(i3.3)') nd
      call computeFilter(WN,computeDelta(U)*nd,U,Ufilk,filter)
      call computeFilter(WN,computeDelta(V)*nd,V,Vfilk,filter)
      call computeFilter(WN,computeDelta(W)*nd,W,Wfilk,filter)
      IF (.NOT. showDivU(Ufilk,Vfilk,Wfilk,WN,spec_rank)) return
      IF (.NOT. NullifyDiv(Ufilk, Vfilk, WfilK, WN, spec_rank)) return 
      IF (.NOT. showDivU(Ufilk,Vfilk,Wfilk,WN,spec_rank)) return
      call btran(Ufilk,Ufil,res)
      call btran(Vfilk,Vfil,res)
      call btran(Wfilk,Wfil,res)

      !Dissipation Exact
      call computeT_ijVecA(U,V,W,Ufil,Vfil,Wfil,WN,S(1,1),S(1,2),S(1,3),S(2,2),S(2,3),S(3,3),filter,res,computeDelta(U)*nd)
      call ComputeDivSymTijVec1Vec2(T1k,T2k,T3k,S(1,2),S(1,3),S(2,3),S(1,1),S(2,2),S(3,3),WN,res)
      call btran(T1k,S(1,1),res)
      call btran(T2k,S(1,2),res)
      call btran(T3K,S(1,3),res)
      nrjDiss%values = - ( S(1,1)%values * Ufil%values + S(1,2)%values * Vfil%values + S(1,3)%values * Wfil%values )
      nrjDiss%name='Exact_'//trim(adjustl(charNd))
      if (.not. computeFieldPDF(1,nrjDiss,300,spec_rank,1)) return

      if (.not. computeFieldGradient(WN,Ufil,S(1,1),S(1,2),S(1,3),nbcpus,spec_rank)) return
      if (.not. computeFieldGradient(WN,Vfil,S(2,1),S(2,2),S(2,3),nbcpus,spec_rank)) return
      if (.not. computeFieldGradient(WN,Wfil,S(3,1),S(3,2),S(3,3),nbcpus,spec_rank)) return
      O(1,2)%values = 0.5_WP * ( S(1,2)%values - S(2,1)%values )
      O(1,3)%values = 0.5_WP * ( S(1,3)%values - S(3,1)%values )
      O(2,3)%values = 0.5_WP * ( S(2,3)%values - S(3,2)%values )
      O(2,1)%values = - O(1,2)%values 
      O(3,1)%values = - O(1,3)%values  
      O(3,2)%values = - O(2,3)%values  
      O(1,1)%values = 0.0_WP 
      O(2,2)%values = 0.0_WP
      O(3,3)%values = 0.0_WP
      call computeStressTensor(Ufilk,Vfilk,Wfilk,WN,S(1,1),S(1,2),S(1,3),S(2,2),S(2,3),S(3,3),res)
      S(2,1)%values = S(1,2)%values
      S(3,1)%values = S(1,3)%values
      S(3,2)%values = S(2,3)%values
      if (.not. computeSplitingSij(S(1,1),S(1,2),S(1,3),S(2,2),S(2,3),S(3,3),Sm(1,1),Sm(1,2),Sm(1,3),Sm(2,2),Sm(2,3),Sm(3,3),'moins')) return
      Sm(2,1)%values = Sm(1,2)%values
      Sm(3,1)%values = Sm(1,3)%values
      Sm(3,2)%values = Sm(2,3)%values
      if (.not. computeSplitingSij(S(1,1),S(1,2),S(1,3),S(2,2),S(2,3),S(3,3),Spp(1,1),Spp(1,2),Spp(1,3),Spp(2,2),Spp(2,3),Spp(3,3),'plus')) return
      Spp(2,1)%values = Spp(1,2)%values
      Spp(3,1)%values = Spp(1,3)%values
      Spp(3,2)%values = Spp(2,3)%values
  

      !SikmoinsSkjSij
      resS%values=0.0_WP
      do l=U%zmin,U%zmax
        do m=U%ymin,U%ymax
          do n=U%xmin,U%xmax
            do i=1,3
              do j=1,3
                do k=1,3
                  resS%values(n,m,l) = resS%values(n,m,l) + Sm(i,k)%values(n,m,l) * S(k,j)%values(n,m,l) * S(i,j)%values(n,m,l)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      resS%name='sikmoinsskjsij_'//trim(adjustl(charNd))
      if (.not. computeFieldPDF(1,resS,300,spec_rank,1)) return
      if (.not. computeTwoFieldsPDF(1,nrjDiss,resS,(/300,300/),spec_rank,1)) return
      !SikplusSkjSij
      resS%values=0.0_WP
      do l=U%zmin,U%zmax
        do m=U%ymin,U%ymax
          do n=U%xmin,U%xmax
            do i=1,3
              do j=1,3
                do k=1,3
                  resS%values(n,m,l) = resS%values(n,m,l) + Spp(i,k)%values(n,m,l) * S(k,j)%values(n,m,l) * S(i,j)%values(n,m,l)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      resS%name='sikplusskjsij_'//trim(adjustl(charNd))
      if (.not. computeFieldPDF(1,resS,300,spec_rank,1)) return
      if (.not. computeTwoFieldsPDF(1,nrjDiss,resS,(/300,300/),spec_rank,1)) return
      !OikOjkSij
      resS%values=0.0_WP
      do l=U%zmin,U%zmax
        do m=U%ymin,U%ymax
          do n=U%xmin,U%xmax
            do i=1,3
              do j=1,3
                do k=1,3
                    resS%values(n,m,l) = resS%values(n,m,l) + O(i,k)%values(n,m,l) * O(j,k)%values(n,m,l) * S(i,j)%values(n,m,l)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      resS%name='OikOjkSij_'//trim(adjustl(charNd))
      if (.not. computeFieldPDF(1,resS,300,spec_rank,1)) return
      if (.not. computeTwoFieldsPDF(1,nrjDiss,resS,(/300,300/),spec_rank,1)) return
      !(Sikmoins-Sikplus)SkjSij
      resS%values=0.0_WP
      do l=U%zmin,U%zmax
        do m=U%ymin,U%ymax
          do n=U%xmin,U%xmax
            do i=1,3
              do j=1,3
                do k=1,3
                    resS%values(n,m,l) = resS%values(n,m,l) + ( Sm(i,k)%values(n,m,l) - Spp(i,k)%values(n,m,l) ) * S(k,j)%values(n,m,l) * S(i,j)%values(n,m,l)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      resS%name='(Sikmoins-Sikplus)SkjSij_'//trim(adjustl(charNd))
      if (.not. computeFieldPDF(1,resS,300,spec_rank,1)) return
      if (.not. computeTwoFieldsPDF(1,nrjDiss,resS,(/300,300/),spec_rank,1)) return
      !(Sikmoins+Sikplus)SkjSij
      resS%values=0.0_WP
      do l=U%zmin,U%zmax
        do m=U%ymin,U%ymax
          do n=U%xmin,U%xmax
            do i=1,3
              do j=1,3
                do k=1,3
                    resS%values(n,m,l) = resS%values(n,m,l) + ( Sm(i,k)%values(n,m,l) + Spp(i,k)%values(n,m,l) ) * S(k,j)%values(n,m,l) * S(i,j)%values(n,m,l)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      resS%name='(Sikmoins+Sikplus)SkjSij_'//trim(adjustl(charNd))
      if (.not. computeFieldPDF(1,resS,300,spec_rank,1)) return
      if (.not. computeTwoFieldsPDF(1,nrjDiss,resS,(/300,300/),spec_rank,1)) return
      !(-SikplusSkj+OikOjk)Sij
      resS%values=0.0_WP
      do l=U%zmin,U%zmax
        do m=U%ymin,U%ymax
          do n=U%xmin,U%xmax
            do i=1,3
              do j=1,3
                do k=1,3
                    resS%values(n,m,l) = resS%values(n,m,l) + ( -Spp(i,k)%values(n,m,l) * S(k,j)%values(n,m,l) + O(i,k)%values(n,m,l) * O(j,k)%values(n,m,l) ) * S(i,j)%values(n,m,l)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      resS%name='(-SikplusSkj+OikOjk)Sij_'//trim(adjustl(charNd))
      if (.not. computeFieldPDF(1,resS,300,spec_rank,1)) return
      if (.not. computeTwoFieldsPDF(1,nrjDiss,resS,(/300,300/),spec_rank,1)) return
      !(SikplusSkj+OikOjk)Sij
      resS%values=0.0_WP
      do l=U%zmin,U%zmax
        do m=U%ymin,U%ymax
          do n=U%xmin,U%xmax
            do i=1,3
              do j=1,3
                do k=1,3
                    resS%values(n,m,l) = resS%values(n,m,l) + ( Spp(i,k)%values(n,m,l) * S(k,j)%values(n,m,l) + O(i,k)%values(n,m,l) * O(j,k)%values(n,m,l) ) * S(i,j)%values(n,m,l)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      resS%name='(SikplusSkj+OikOjk)Sij_'//trim(adjustl(charNd))
      if (.not. computeFieldPDF(1,resS,300,spec_rank,1)) return
      if (.not. computeTwoFieldsPDF(1,nrjDiss,resS,(/300,300/),spec_rank,1)) return
      !(SikSkj+OikOjk)Sij
      resS%values=0.0_WP
      do l=U%zmin,U%zmax
        do m=U%ymin,U%ymax
          do n=U%xmin,U%xmax
            do i=1,3
              do j=1,3
                do k=1,3
                    resS%values(n,m,l) = resS%values(n,m,l) + ( S(i,k)%values(n,m,l) * S(k,j)%values(n,m,l) + O(i,k)%values(n,m,l) * O(j,k)%values(n,m,l) ) * S(i,j)%values(n,m,l)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      resS%name='(SikSkj+OikOjk)Sij_'//trim(adjustl(charNd))
      if (.not. computeFieldPDF(1,resS,300,spec_rank,1)) return
      if (.not. computeTwoFieldsPDF(1,nrjDiss,resS,(/300,300/),spec_rank,1)) return
      !(SikmoinsSkj+OikOjk)Sij
      resS%values=0.0_WP
      do l=U%zmin,U%zmax
        do m=U%ymin,U%ymax
          do n=U%xmin,U%xmax
            do i=1,3
              do j=1,3
                do k=1,3
                    resS%values(n,m,l) = resS%values(n,m,l) + ( Sm(i,k)%values(n,m,l) * S(k,j)%values(n,m,l) + O(i,k)%values(n,m,l) * O(j,k)%values(n,m,l) ) * S(i,j)%values(n,m,l)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      resS%name='(SikmoinsSkj+OikOjk)Sij_'//trim(adjustl(charNd))
      if (.not. computeFieldPDF(1,resS,300,spec_rank,1)) return
      if (.not. computeTwoFieldsPDF(1,nrjDiss,resS,(/300,300/),spec_rank,1)) return
      !test
      O(1,1)%values = Sm(1,1)%values*S(1,1)%values+ Sm(1,2)%values*S(1,2)%values+ Sm(1,3)%values*S(1,3)%values
      O(1,2)%values = Sm(1,1)%values*S(1,2)%values+ Sm(1,2)%values*S(2,2)%values+ Sm(1,3)%values*S(2,3)%values
      O(1,3)%values = Sm(1,1)%values*S(1,3)%values+ Sm(1,2)%values*S(2,3)%values+ Sm(1,3)%values*S(3,3)%values
      O(2,2)%values = Sm(1,2)%values*S(1,2)%values+ Sm(2,2)%values*S(2,2)%values+ Sm(2,3)%values*S(2,3)%values
      O(2,3)%values = Sm(1,2)%values*S(1,3)%values+ Sm(2,2)%values*S(2,3)%values+ Sm(2,3)%values*S(3,3)%values
      O(3,3)%values = Sm(1,3)%values*S(1,3)%values+ Sm(2,3)%values*S(2,3)%values+ Sm(3,3)%values*S(3,3)%values
      resS%values = O(1,1)%values*S(1,1)%values + O(2,2)%values*S(2,2)%values + O(3,3)%values*S(3,3)%values +&
               & 2.0_WP * (O(1,2)%values*S(1,2)%values + O(1,3)%values*S(1,3)%values + O(2,3)%values*S(2,3)%values)
      resS%name='(testsik-sjk)Sij_'//trim(adjustl(charNd))
      if (.not. computeFieldPDF(1,resS,300,spec_rank,1)) return
      if (.not. computeTwoFieldsPDF(1,nrjDiss,resS,(/300,300/),spec_rank,1)) return
    enddo

    ! ===== Free memory and delete temporary fields =====
    if( .not. deleteWN(WN)) return
    call deleteDataLayout(S(1,1))
    call deleteDataLayout(S(1,2))
    call deleteDataLayout(S(1,3))
    call deleteDataLayout(S(2,1))
    call deleteDataLayout(S(3,1))
    call deleteDataLayout(S(2,3))
    call deleteDataLayout(S(3,2))
    call deleteDataLayout(S(2,2))
    call deleteDataLayout(S(3,3))
    call deleteDataLayout(O(1,1))
    call deleteDataLayout(O(2,2))
    call deleteDataLayout(O(3,3))
    call deleteDataLayout(O(1,2))
    call deleteDataLayout(O(1,3))
    call deleteDataLayout(O(2,3))
    call deleteDataLayout(O(2,1))
    call deleteDataLayout(O(3,1))
    call deleteDataLayout(O(3,2))
    call deleteDataLayout(Spp(1,1))
    call deleteDataLayout(Spp(2,2))
    call deleteDataLayout(Spp(3,3))
    call deleteDataLayout(Spp(1,3))
    call deleteDataLayout(Spp(2,3))
    call deleteDataLayout(Spp(1,2))
    call deleteDataLayout(Spp(3,1))
    call deleteDataLayout(Spp(3,2))
    call deleteDataLayout(Spp(2,1))
    call deleteDataLayout(Sm(1,1))
    call deleteDataLayout(Sm(2,2))
    call deleteDataLayout(Sm(3,3))
    call deleteDataLayout(Sm(1,3))
    call deleteDataLayout(Sm(2,3))
    call deleteDataLayout(Sm(1,2))
    call deleteDataLayout(Sm(3,1))
    call deleteDataLayout(Sm(3,2))
    call deleteDataLayout(Sm(2,1))
    call deleteDataLayout(resS)
    call deleteDataLayout(Ufilk)
    call deleteDataLayout(Vfilk)
    call deleteDataLayout(Wfilk)
    call deleteDataLayout(Ufil)
    call deleteDataLayout(Vfil)
    call deleteDataLayout(Wfil)
    call deleteDataLayout(T1k)
    call deleteDataLayout(T2k)
    call deleteDataLayout(T3k)
    call deleteDataLayout(nrjDiss)

    res = .true.

end subroutine testDissipationduidxkdujdxk



!-------------------------------------------------------------------------------------------
!> @detail
!> check the equality between gradient model and gradient model with incompressibility 
!!
!! @author Antoine Vollant, LEGI
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in]    nbcpus          = number of process in the cpu pool
!!    @param[in]    U               = Component of velocity
!!    @param[in]    V               = Component of velocity
!!    @param[in]    W               = Component of velocity
!!    @return       res             = logical value
!
!-------------------------------------------------------------------------------------------
subroutine checkGradientModels(spec_rank,nbcpus,res,U,V,W)

    use filtering_tools
    use datalayout
    use differential_tools , only : computeDerivationField

    ! Input/output
    logical, intent(out)                    :: res
    integer, intent(in)                     :: spec_rank, nbcpus
    type(REAL_DATA_LAYOUT), intent(in)      :: U,V,W

    ! Local field
    type(WAVENUMBERS)                       :: WN
    integer                                 :: filter,nd 
    real(WP)                                :: delta,t1,t2,t3
    logical                                 :: success
    type(COMPLEX_DATA_LAYOUT)               :: T1GM1K,T2GM1K,T3GM1K,Uk 
    type(COMPLEX_DATA_LAYOUT)               :: T1GM2K,T2GM2K,T3GM2K
    type(REAL_DATA_LAYOUT)                  :: T1GM1,T2GM1,T3GM1 
    type(REAL_DATA_LAYOUT)                  :: T1GM2,T2GM2,T3GM2
!    type(COMPLEX_DATA_LAYOUT)               :: Ukfil,Vkfil,Wkfil 
    type(REAL_DATA_LAYOUT)                  :: Ufil,Vfil,Wfil

    ! ===== Initialisation =====
    res = .false.
    if (.not. initDataLayout("Uk", Uk,(U%nx/2)+1,U%ny,U%nz, &
       & U%Lx,U%Ly,U%Lz,nbcpus,spec_rank,alongZ)) then 
      write(6,'(a,i0,a)')'[ERROR] Not able to create Exk'
      return
    endif
    if (.not. copystructonly(U,T1GM1) ) return
    if (.not. copystructonly(U,T2GM1) ) return
    if (.not. copystructonly(U,T3GM1) ) return
    if (.not. copystructonly(U,T1GM2) ) return
    if (.not. copystructonly(U,T2GM2) ) return
    if (.not. copystructonly(U,T3GM2) ) return
    if (.not. copystructonly(Uk,T1GM1K) ) return
    if (.not. copystructonly(Uk,T2GM1K) ) return
    if (.not. copystructonly(Uk,T3GM1K) ) return
    if (.not. copystructonly(Uk,T1GM2K) ) return
    if (.not. copystructonly(Uk,T2GM2K) ) return
    if (.not. copystructonly(Uk,T3GM2K) ) return
    if (.not. copystructonly(U,Ufil) ) return
    if (.not. copystructonly(U,Vfil) ) return
    if (.not. copystructonly(U,Wfil) ) return
    !Calcul des nombres d'onde
    if (.not. initWN(WN,spec_rank,U%nx,U%ny,U%nz) .or. &
       &.not.  computeWN(WN,spec_rank,U%Lx,U%Ly,U%Lz,U%nx,U%ny,U%nz) ) then 
      write(6,'(a,i0,a)')'[ERROR] Cannot compute WN'
      return
    endif
    nd = 8
    filter = 1
    delta=nd*computeDelta(U)
    call computeFilter(WN,delta,U,Ufil,filter)
    call computeFilter(WN,delta,V,Vfil,filter)
    call computeFilter(WN,delta,W,Wfil,filter)
    call cpu_time(t1)
    call gradientVelIncomp(T1GM1K,T2GM1K,T3GM1K,Ufil,Vfil,Wfil,WN,spec_rank,delta=delta)
    call cpu_time(t2)
    call gradientVelDIV(T1GM2K,T2GM2K,T3GM2K,Ufil,Vfil,Wfil,WN,delta,spec_rank,nbcpus)
    call cpu_time(t3)
    call btran(T1GM1K,T1GM1,success)
    if (.not. success) return
    call btran(T2GM1K,T2GM1,success)
    if (.not. success) return
    call btran(T3GM1K,T3GM1,success)
    if (.not. success) return
    call btran(T1GM2K,T1GM2,success)
    if (.not. success) return
    call btran(T2GM2K,T2GM2,success)
    if (.not. success) return
    call btran(T3GM2K,T3GM2,success)
    if (.not. success) return

    if (spec_rank .eq. 0) then
      write (6,'(a,1x,g15.8,1x,a)') '[INFO] gradientVelIncomp',t2-t1,'s'
      write (6,'(a,1x,g15.8,1x,a)') '[INFO] gradientVel      ',t3-t2,'s'
    endif
    if (.not. equalToVal(T1GM1,T1GM2) .or. &
       &.not. equalToVal(T2GM1,T2GM2) .or. &
       &.not. equalToVal(T3GM1,T3GM2)) then
     write (6,'(a,1x,3(g15.8,1x))') '[INFO] MODELS ARE NOT EQUIVALENT'
     write (6,'(a,1x,3(g15.8,1x))') '[INFO]',T1GM1%values(U%nx/2,U%ymin,U%zmin),&
                                            &T2GM1%values(U%nx/2,U%ymin,U%zmin),&
                                            &T3GM1%values(U%nx/2,U%ymin,U%zmin)
     write (6,'(a,1x,3(g15.8,1x))') '[INFO]',T1GM2%values(U%nx/2,U%ymin,U%zmin),&
                                            &T2GM2%values(U%nx/2,U%ymin,U%zmin),&
                                            &T3GM2%values(U%nx/2,U%ymin,U%zmin)
    else
     write (6,'(a,1x,3(g15.8,1x))') '[INFO] MODELS ARE EQUIVALENT !!!!'
    endif

    ! ===== Free memory and delete temporary fields =====
    if( .not. deleteWN(WN)) return
    call deleteDataLayout(Uk)
    call deleteDataLayout(T1GM1K)
    call deleteDataLayout(T2GM1K)
    call deleteDataLayout(T3GM1K)
    call deleteDataLayout(T1GM2K)
    call deleteDataLayout(T2GM2K)
    call deleteDataLayout(T3GM2K)
    call deleteDataLayout(T1GM1)
    call deleteDataLayout(T2GM1)
    call deleteDataLayout(T3GM1)
    call deleteDataLayout(T1GM2)
    call deleteDataLayout(T2GM2)
    call deleteDataLayout(T3GM2)
    call deleteDataLayout(Ufil)
    call deleteDataLayout(Vfil)
    call deleteDataLayout(Wfil)

    res = .true.

end subroutine checkGradientModels 



!-------------------------------------------------------------------------------------------
!> @detail
!> check the equality between gradient model and gradient model with incompressibility 
!!
!! @author Antoine Vollant, LEGI
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in]    nbcpus          = number of process in the cpu pool
!!    @param[in]    U               = Component of velocity
!!    @param[in]    V               = Component of velocity
!!    @param[in]    W               = Component of velocity
!!    @return       res             = logical value
!
!-------------------------------------------------------------------------------------------
subroutine aprioriRGMVel(iFic,spec_rank,nbcpus,res,U,V,W)

    use filtering_tools
    use datalayout
    use differential_tools , only : computeDerivationField
    use conditional_mean_tools
    use fileio 

    ! Input/output
    logical, intent(out)                    :: res
    integer, intent(in)                     :: spec_rank,nbcpus,iFic
    type(REAL_DATA_LAYOUT), intent(in)      :: U,V,W

    ! Local field
    type(WAVENUMBERS)                       :: WN
    integer                                 :: filter,nd,file_id
    real(WP)                                :: delta
    logical                                 :: success,existe
    type(COMPLEX_DATA_LAYOUT)               :: T1RGMK,T2RGMK,T3RGMK
    type(COMPLEX_DATA_LAYOUT)               :: Ukfil,Vkfil,Wkfil 
    type(REAL_DATA_LAYOUT)                  :: T1GM,T2GM,T3GM 
    type(REAL_DATA_LAYOUT)                  :: T1DNS,T2DNS,T3DNS
    type(REAL_DATA_LAYOUT)                  :: T1RGM,T2RGM,T3RGM
    type(REAL_DATA_LAYOUT)                  :: Ufil,Vfil,Wfil
    type(REAL_DATA_LAYOUT)                  :: dudx,dudy,dudz 
    type(REAL_DATA_LAYOUT),dimension(:),allocatable :: uvw 
    INTEGER                                 :: bin, typeDisc
    INTEGER                                 :: i,j,k 
    REAL(WP)                                :: T1GMminus,T1RGMminus,T2RGMminus,T2GMminus,T3RGMminus,T3GMminus
    REAL(WP)                                :: mseQuadGM1,mseQuadGM2,mseQuadGM3
    REAL(WP)                                :: mseQuadRGM1,mseQuadRGM2,mseQuadRGM3
    REAL(WP)                                :: mseQuadRGM1_plus,mseQuadRGM2_plus,mseQuadRGM3_plus
    REAL(WP)                                :: mseIrrGM1,mseIrrGM2,mseIrrGM3
    REAL(WP)                                :: mseQuadSmag1,mseQuadSmag2,mseQuadSmag3
    REAL(WP)                                :: mseQuadSmagStat1,mseQuadSmagStat2,mseQuadSmagStat3
    REAL(WP)                                :: mseIrrSmag1,mseIrrSmag2,mseIrrSmag3
    REAL(WP)                                :: mseIrrRGM1,mseIrrRGM2,mseIrrRGM3
    REAL(WP)                                :: mseIrrGMVort1,mseIrrGMVort2,mseIrrGMVort3
    REAL(WP)                                :: mseQuadGMVort1,mseQuadGMVort2,mseQuadGMVort3
    REAL(WP)                                :: mseIrrRGM1_plus,mseIrrRGM2_plus,mseIrrRGM3_plus
    REAL(WP)                                :: mseIrrRGM1_2var,mseIrrRGM2_2var,mseIrrRGM3_2var
    REAL(WP)                                :: mseQuadRGM1_2var,mseQuadRGM2_2var,mseQuadRGM3_2var
    REAL(WP)                                :: mseQuadRGM1dyn,mseQuadRGM2dyn,mseQuadRGM3dyn
    REAL(WP)                                :: mseIrrSikmoinsSkj1,mseIrrSikmoinsSkj2,mseIrrSikmoinsSkj3
    REAL(WP)                                :: mseIrrSikmoinsSkjOikOjk1,mseIrrSikmoinsSkjOikOjk2,mseIrrSikmoinsSkjOikOjk3
    REAL(WP)                                :: mseQuadSikmoinsSkj1,mseQuadSikmoinsSkj2,mseQuadSikmoinsSkj3
    REAL(WP)                                :: mseIrrSikmoinsSkjSikOjkOikSjk1,mseIrrSikmoinsSkjSikOjkOikSjk2,mseIrrSikmoinsSkjSikOjkOikSjk3 
    REAL(WP)                                :: mseIrrSikOjkOikSjk1,mseIrrSikOjkOikSjk2,mseIrrSikOjkOikSjk3
    REAL(WP)                                :: mseIrrSikmoinsSkj_OikOjk1,mseIrrSikmoinsSkj_OikOjk2,mseIrrSikmoinsSkj_OikOjk3
    REAL(WP)                                :: mseIrrOikOjk1,mseIrrOikOjk2,mseIrrOikOjk3
    REAL(WP)                                :: mseIrrSikmoinsSkj_SikOjkOikSjk1,mseIrrSikmoinsSkj_SikOjkOikSjk2,mseIrrSikmoinsSkj_SikOjkOikSjk3
    REAL(WP)                                :: mseIrrgraddOkdxldUjdxl1,mseIrrgraddOkdxldUjdxl2,mseIrrgraddOkdxldUjdxl3
    REAL(WP)                                :: mseIrrgraddOkdxlSjlmoins1,mseIrrgraddOkdxlSjlmoins2,mseIrrgraddOkdxlSjlmoins3
    REAL(WP)                                :: mseIrrgraddOkdxlSjlplus1,mseIrrgraddOkdxlSjlplus2,mseIrrgraddOkdxlSjlplus3
    REAL(WP)                                :: mseIrrgraddOkdxlOmegajl1,mseIrrgraddOkdxlOmegajl2,mseIrrgraddOkdxlOmegajl3
    REAL(WP)                                :: mseRGMSikmoinsSkjSikOjkOikSjkDyn1,mseRGMSikmoinsSkjSikOjkOikSjkDyn2,mseRGMSikmoinsSkjSikOjkOikSjkDyn3 
    REAL(WP)                                :: mseRGMSikmoinsSkjSikOjkOikSjkDyn21,mseRGMSikmoinsSkjSikOjkOikSjkDyn22,mseRGMSikmoinsSkjSikOjkOikSjkDyn23
    REAL(WP)                                :: mseRGMSikmoinsSkjDyn1,mseRGMSikmoinsSkjDyn2,mseRGMSikmoinsSkjDyn3 
    REAL(WP)                                :: dissExact,dissRGM,dissGM,diffdissRGM
    REAL(WP)                                :: diffExact,diffRGM,diffGM,diffdissGM,diffdissExacte 
    REAL(WP),dimension(3)                   :: dissExactComp,dissRGMComp,dissGMComp 
    REAL(WP),dimension(3)                   :: diffExactComp,diffRGMComp,diffGMComp
    REAL(WP)                                :: dissIrrgraddOkdxldUjdxl,dissIrrgraddOkdxlOmegajl,dissIrrgraddOkdxlSjlmoins
    REAL(WP)                                :: dissIrrgraddOkdxlSjlplus,dissIrrOikOjk,dissIrrQDMvort,dissIrrSikmoinsSkj
    REAL(WP)                                :: dissIrrSikmoinsSkj_OikOjk,dissIrrSikmoinsSkj_SikOjkOikSjk,dissIrrSikmoinsSkjDyn
    REAL(WP)                                :: dissIrrSikmoinsSkjOikOjk,dissIrrSikmoinsSkjSikOjkOikSjk,dissIrrSikOjkOikSjk,dissQDMvort
    REAL(WP)                                :: dissSikmoinsSkjStat,dissSmagDyna,dissSmagIrr,dissSmagStat,dissIrrSikmoinsSkjSikOjkOikSjkDyn,dissIrrSikmoinsSkjSikOjkOikSjkDyn2
    REAL(WP)                                :: mseRGMSikmoinsSkjDynFab1,mseRGMSikmoinsSkjDynFab2,mseRGMSikmoinsSkjDynFab3
    REAL(WP)                                :: mseSikmoinsSkjSikOjkOikSjkDynFab1,mseSikmoinsSkjSikOjkOikSjkDynFab2,mseSikmoinsSkjSikOjkOikSjkDynFab3
    REAL(WP)                                :: mseSikmoinsSkjSikOjkOikSjkDyn2Fab1,mseSikmoinsSkjSikOjkOikSjkDyn2Fab2,mseSikmoinsSkjSikOjkOikSjkDyn2Fab3
    REAL(WP)                                :: dissSikmoinsSkjSikOjkOikSjkDyn2Fab,dissSikmoinsSkjSikOjkOikSjkDynFab,dissSikmoinsSkjDynFab,dissSikmoinsSkjDyn

    ! ===== Initialisation =====
    res = .false.
    if (.not. initDataLayout("T1RGMK",T1RGMK,(U%nx/2)+1,U%ny,U%nz, &
       & U%Lx,U%Ly,U%Lz,nbcpus,spec_rank,alongZ)) then 
      write(6,'(a,i0,a)')'[ERROR] Not able to create Exk'
      return
    endif
    if (.not. copystructonly(T1RGMK,T2RGMK) ) return
    if (.not. copystructonly(T1RGMK,T3RGMK) ) return
    if (.not. copystructonly(T1RGMK,Ukfil) ) return
    if (.not. copystructonly(T1RGMK,Vkfil) ) return
    if (.not. copystructonly(T1RGMK,Wkfil) ) return
    if (.not. copystructonly(U,T1GM) ) return
    if (.not. copystructonly(U,T2GM) ) return
    if (.not. copystructonly(U,T3GM) ) return
    if (.not. copystructonly(U,T1RGM) ) return
    if (.not. copystructonly(U,T2RGM) ) return
    if (.not. copystructonly(U,T3RGM) ) return
    if (.not. copystructonly(U,T1DNS) ) return
    if (.not. copystructonly(U,T2DNS) ) return
    if (.not. copystructonly(U,T3DNS) ) return
    if (.not. copystructonly(U,Ufil) ) return
    if (.not. copystructonly(U,Vfil) ) return
    if (.not. copystructonly(U,Wfil) ) return
    if (.not. copystructonly(U,dudx) ) return
    if (.not. copystructonly(U,dudy) ) return
    if (.not. copystructonly(U,dudz) ) return
    if (.not. initWorkArray(U,3,uvw)) return
    !Calcul des nombres d'onde
    if (.not. initWN(WN,spec_rank,U%nx,U%ny,U%nz) .or. &
       &.not.  computeWN(WN,spec_rank,U%Lx,U%Ly,U%Lz,U%nx,U%ny,U%nz) ) then 
      write(6,'(a,i0,a)')'[ERROR] Cannot compute WN'
      return
    endif
    ! cutoff -> 1
    ! box -> 2
    filter = 1
    if ( parser_is_defined('Filter type for post process')) then
      if (.not. parseFilter(filter,spec_rank,other='Filter type for post process')) return
    endif
    bin=250
    typeDisc=1
    do nd=2,U%nx/16,2
      delta=nd*computeDelta(U)
      call computeFilter(WN,delta,U,Ukfil,filter)
      call computeFilter(WN,delta,V,Vkfil,filter)
      call computeFilter(WN,delta,W,Wkfil,filter)
      IF (.NOT. showDivU(Ukfil,Vkfil,Wkfil,WN,spec_rank)) return
      IF (.NOT. NullifyDiv(Ukfil,Vkfil,Wkfil, WN, spec_rank)) return 
      IF (.NOT. showDivU(Ukfil,Vkfil,Wkfil,WN,spec_rank)) return
      call btran(Ukfil,Ufil,res)
      call btran(Vkfil,Vfil,res)
      call btran(Wkfil,Wfil,res)
      !DNS
      call computeT_ijVecA(U,V,W,Ufil,Vfil,Wfil,WN,T1RGM,T2RGM,T3RGM,T1GM,T2GM,T3GM,filter,success,delta)
      call ComputeDivSymTijVec1Vec2(T1RGMK,T2RGMK,T3RGMK,T2RGM,T3RGM,T2GM,T1RGM,T1GM,T3GM,WN,res)
      call btran(T1RGMK,T1DNS,success)
      if (.not. success) return
      call btran(T2RGMK,T2DNS,success)
      if (.not. success) return
      call btran(T3RGMK,T3DNS,success)
      if (.not. success) return
      !RGM
      call RGMVel(T1RGMK,T2RGMk,T3RGMk,U,Ufil,Vfil,Wfil,WN,spec_rank,'moins',delta=delta)
      call btran(T1RGMK,T1RGM,success)
      if (.not. success) return
      call btran(T2RGMK,T2RGM,success)
      if (.not. success) return
      call btran(T3RGMK,T3RGM,success)
      if (.not. success) return
      !GM
      !call gradientVelIncomp(T1RGMK,T2RGMK,T3RGMK,Ufil,Vfil,Wfil,WN,spec_rank,delta=delta,coef=1.0_WP)
      call gradientVelDIV(T1RGMK,T2RGMK,T3RGMK,Ufil,Vfil,Wfil,WN,delta,spec_rank,nbcpus)
      call btran(T1RGMK,T1GM,success)
      if (.not. success) return
      call btran(T2RGMK,T2GM,success)
      if (.not. success) return
      call btran(T3RGMK,T3GM,success)
      if (.not. success) return
      !Calcul transfert totaux
      dudx%values    = T1DNS%values * Ufil%values + T2DNS%values * Vfil%values + T3DNS%values * Wfil%values
      diffdissExacte = computeFieldAvg(dudx,spec_rank)
      dudy%values    = T1RGM%values * Ufil%values + T2RGM%values * Vfil%values + T3RGM%values * Wfil%values
      diffdissRGM    = computeFieldAvg(dudy,spec_rank)
      dudz%values    = T1GM%values * Ufil%values + T2GM%values * Vfil%values + T3GM%values * Wfil%values
      diffdissGM     = computeFieldAvg(dudz,spec_rank)
      !Calcul des erreurs irr et quad GM
      if (.not. computeEONVar( spec_rank,(/T1GM/),T1DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudx) ) return 
      if (.not. computeEONVar( spec_rank,(/T2GM/),T2DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudy) ) return 
      if (.not. computeEONVar( spec_rank,(/T3GM/),T3DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudz) ) return 
      mseIrrGM1  = computeMSE(T1DNS,dudx,spec_rank,.true.)
      mseIrrGM2  = computeMSE(T2DNS,dudy,spec_rank,.true.)
      mseIrrGM3  = computeMSE(T3DNS,dudz,spec_rank,.true.)
      mseQuadGM1 = computeMSE(T1DNS,T1GM,spec_rank,.true.)
      mseQuadGM2 = computeMSE(T2DNS,T2GM,spec_rank,.true.)
      mseQuadGM3 = computeMSE(T3DNS,T3GM,spec_rank,.true.)
      !Calcul des erreurs irr et quad RGM ( C delta^2 d^2ui/dxkdxj S^-_jk )
      if (.not. computeEONVar( spec_rank,(/T1RGM/),T1DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudx) ) return 
      if (.not. computeEONVar( spec_rank,(/T2RGM/),T2DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudy) ) return 
      if (.not. computeEONVar( spec_rank,(/T3RGM/),T3DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudz) ) return 
      mseIrrRGM1  = computeMSE(T1DNS,dudx,spec_rank,.true.)
      mseIrrRGM2  = computeMSE(T2DNS,dudy,spec_rank,.true.)
      mseIrrRGM3  = computeMSE(T3DNS,dudz,spec_rank,.true.)
      mseQuadRGM1 = computeMSE(T1DNS,T1RGM,spec_rank,.true.)
      mseQuadRGM2 = computeMSE(T2DNS,T2RGM,spec_rank,.true.)
      mseQuadRGM3 = computeMSE(T3DNS,T3RGM,spec_rank,.true.)
      !Calcul de l'erreur quad du RGM Plus Moins Dyn ( C delta^2 d^2ui/dxkdxj S^-_jk + Cd delta^2 d^2ui/dxkdxj S^+_jk )
      call RGMDivDynV2(T1RGMK,T2RGMk,T3RGMk,U,Ufil,Vfil,Wfil,WN,spec_rank,filter)
      call btran(T1RGMK,T1GM,success)
      if (.not. success) return
      call btran(T2RGMK,T2GM,success)
      if (.not. success) return
      call btran(T3RGMK,T3GM,success)
      if (.not. success) return
      mseQuadRGM1_2var  = computeMSE(T1DNS,T1GM,spec_rank,.true.)
      mseQuadRGM2_2var  = computeMSE(T2DNS,T2GM,spec_rank,.true.)
      mseQuadRGM3_2var  = computeMSE(T3DNS,T3GM,spec_rank,.true.)
      !Calcul de l'erreur quad du RGM Moins Dyn ( Cd delta^2 d^2ui/dxkdxj S^-_jk )
      call RGMDivDynV1(T1RGMK,T2RGMk,T3RGMk,U,Ufil,Vfil,Wfil,WN,spec_rank,filter)
      call btran(T1RGMK,T1GM,success)
      if (.not. success) return
      call btran(T2RGMK,T2GM,success)
      if (.not. success) return
      call btran(T3RGMK,T3GM,success)
      if (.not. success) return
      mseQuadRGM1Dyn= computeMSE(T1DNS,T1GM,spec_rank,.true.)
      mseQuadRGM2Dyn= computeMSE(T2DNS,T2GM,spec_rank,.true.)
      mseQuadRGM3Dyn= computeMSE(T3DNS,T3GM,spec_rank,.true.)
      !Calcul des erreurs irr et quad RGM plus ( C delta^2 d^2ui/dxkdxj S^+_jk )
      call RGMVel(T1RGMK,T2RGMk,T3RGMk,U,Ufil,Vfil,Wfil,WN,spec_rank,'plus',delta=delta)
      call btran(T1RGMK,T1GM,success)
      if (.not. success) return
      call btran(T2RGMK,T2GM,success)
      if (.not. success) return
      call btran(T3RGMK,T3GM,success)
      if (.not. success) return
      if (.not. computeEONVar( spec_rank,(/T1GM/),T1DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudx) ) return 
      if (.not. computeEONVar( spec_rank,(/T2GM/),T2DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudy) ) return 
      if (.not. computeEONVar( spec_rank,(/T3GM/),T3DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudz) ) return 
      mseIrrRGM1_plus  = computeMSE(T1DNS,dudx,spec_rank,.true.)
      mseIrrRGM2_plus  = computeMSE(T2DNS,dudy,spec_rank,.true.)
      mseIrrRGM3_plus  = computeMSE(T3DNS,dudz,spec_rank,.true.)
      mseQuadRGM1_plus = computeMSE(T1DNS,T1GM,spec_rank,.true.)
      mseQUadRGM2_plus = computeMSE(T2DNS,T2GM,spec_rank,.true.)
      mseQuadRGM3_plus = computeMSE(T3DNS,T3GM,spec_rank,.true.)
      !Calcul des erreurs irr et quad RGM 2 variables
      if (.not. computeEONVar( spec_rank,(/T1RGM,T1GM/),T1DNS,&
                      & (/typeDisc,typeDisc/),(/bin,bin/),&
                      & fieldOpt = dudx) ) return 
      if (.not. computeEONVar( spec_rank,(/T2RGM,T2GM/),T2DNS,&
                      & (/typeDisc,typeDisc/),(/bin,bin/),&
                      & fieldOpt = dudy) ) return 
      if (.not. computeEONVar( spec_rank,(/T3RGM,T3GM/),T3DNS,&
                      & (/typeDisc,typeDisc/),(/bin,bin/),&
                      & fieldOpt = dudz) ) return 
      mseIrrRGM1_2var  = computeMSE(T1DNS,dudx,spec_rank,.true.)
      mseIrrRGM2_2var  = computeMSE(T2DNS,dudy,spec_rank,.true.)
      mseIrrRGM3_2var  = computeMSE(T3DNS,dudz,spec_rank,.true.)
      !Verif 
      call gradientVelDIV(T1RGMK,T2RGMK,T3RGMK,Ufil,Vfil,Wfil,WN,delta,spec_rank,nbcpus)
      call btran(T1RGMK,dudx,success)
      if (.not. success) return
      call btran(T2RGMK,dudy,success)
      if (.not. success) return
      call btran(T3RGMK,dudz,success)
      if (.not. success) return
      T1GM%values = T1GM%values + T1RGM%values
      T2GM%values = T2GM%values + T2RGM%values
      T3GM%values = T3GM%values + T3RGM%values
      if (equalToVal(T1GM,dudx) .and. &
         &equalToVal(T2GM,dudy) .and. &
         &equalToVal(T3GM,dudz) ) then
        if (spec_rank .eq. 0) print *,'[INFO] RGM^plus + RGM^moins == GM'
      else
        if (spec_rank .eq. 0) print *,'[INFO] RGM^plus + RGM^moins != GM'
      endif
      !Calcul erreur irr et quad Smag
      !Calcul de Smag statique pour l'erreur irr
      call smagorinskyVel(T1RGMK,T2RGMk,T3RGMk,Ufil,Vfil,Wfil,Ukfil,Vkfil,Wkfil,WN,1,spec_rank,res)
      if (.not. res) return
      call btran(T1RGMK,T1GM,success)
      if (.not. success) return
      call btran(T2RGMk,T2GM,success)
      if (.not. success) return
      call btran(T3RGMk,T3GM,success)
      if (.not. success) return
      mseQuadSmagStat1 = computeMSE(T1DNS,T1GM,spec_rank,.true.)
      mseQuadSmagStat2 = computeMSE(T2DNS,T2GM,spec_rank,.true.)
      mseQuadSmagStat3 = computeMSE(T3DNS,T3GM,spec_rank,.true.)
      uvw(1)%values = - ( T1GM%values * Ufil%values + T2GM%values * Vfil%values + T3GM%values * Wfil%values )
      dissSmagStat = computeFieldAvg(uvw(1),spec_rank)
      if (.not. computeEONVar( spec_rank,(/T1GM/),T1DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudx) ) return 
      if (.not. computeEONVar( spec_rank,(/T2GM/),T2DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudy) ) return 
      if (.not. computeEONVar( spec_rank,(/T3GM/),T3DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudz) ) return 
      mseIrrSmag1  = computeMSE(T1DNS,dudx,spec_rank,.true.)
      mseIrrSmag2  = computeMSE(T2DNS,dudy,spec_rank,.true.)
      mseIrrSmag3  = computeMSE(T3DNS,dudz,spec_rank,.true.)
      uvw(1)%values = - ( dudx%values * Ufil%values + dudy%values * Vfil%values + dudz%values * Wfil%values )
      dissSmagIrr = computeFieldAvg(uvw(1),spec_rank)
      !Calcul de Smag dynamique pour l'erreur quad 
      call smagorinskyVel(T1RGMK,T2RGMk,T3RGMk,Ufil,Vfil,Wfil,Ukfil,Vkfil,Wkfil,WN,2,spec_rank,res,numVelFilter=1,type_avg=0)
      if (.not. res ) return
      call btran(T1RGMK,T1GM,success)
      if (.not. success) return
      call btran(T2RGMk,T2GM,success)
      if (.not. success) return
      call btran(T3RGMk,T3GM,success)
      if (.not. success) return
      mseQuadSmag1 = computeMSE(T1DNS,T1GM,spec_rank,.true.)
      mseQuadSmag2 = computeMSE(T2DNS,T2GM,spec_rank,.true.)
      mseQuadSmag3 = computeMSE(T3DNS,T3GM,spec_rank,.true.)
      uvw(1)%values = - ( T1GM%values * Ufil%values + T2GM%values * Vfil%values + T3GM%values * Wfil%values )
      dissSmagDyna = computeFieldAvg(uvw(1),spec_rank)
      !Calcul de irr et quad gradient QDM vort
      call gradientQDMvort(T1RGMK,T2RGMK,T3RGMK,Ufil,Vfil,Wfil,WN,spec_rank)
      call btran(T1RGMK,T1GM,success)
      if (.not. success) return
      call btran(T2RGMK,T2GM,success)
      if (.not. success) return
      call btran(T3RGMK,T3GM,success)
      if (.not. success) return
      if (.not. computeEONVar( spec_rank,(/T1GM/),T1DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudx) ) return 
      if (.not. computeEONVar( spec_rank,(/T2GM/),T2DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudy) ) return 
      if (.not. computeEONVar( spec_rank,(/T3GM/),T3DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudz) ) return 
      mseIrrGMVort1 = computeMSE(T1DNS,dudx,spec_rank,.true.)
      mseIrrGMVort2 = computeMSE(T2DNS,dudy,spec_rank,.true.)
      mseIrrGMVort3 = computeMSE(T3DNS,dudz,spec_rank,.true.)
      uvw(1)%values = - ( dudx%values * Ufil%values + dudy%values * Vfil%values + dudz%values * Wfil%values )
      dissIrrQDMvort = computeFieldAvg(uvw(1),spec_rank)
      mseQuadGMVort1 = computeMSE(T1DNS,T1GM,spec_rank,.true.)
      mseQuadGMVort2 = computeMSE(T2DNS,T2GM,spec_rank,.true.)
      mseQuadGMVort3 = computeMSE(T3DNS,T3GM,spec_rank,.true.)
      uvw(1)%values = - ( T1GM%values * Ufil%values + T2GM%values * Vfil%values + T3GM%values * Wfil%values )
      dissQDMvort = computeFieldAvg(uvw(1),spec_rank)
      !Calcul de irr RGMSIKmoinsSKJOikOjk(T1,T2,T3,VecX,VecY,VecZ,Wave,spec_rank,delta,coef,ite)
      call RGMSikmoinsSkjOikOjk(T1RGMK,T2RGMK,T3RGMK,Ufil,Ukfil,Vkfil,Wkfil,WN,spec_rank,delta=delta)
      call btran(T1RGMK,T1GM,success)
      if (.not. success) return
      call btran(T2RGMK,T2GM,success)
      if (.not. success) return
      call btran(T3RGMK,T3GM,success)
      if (.not. success) return
      if (.not. computeEONVar( spec_rank,(/T1GM/),T1DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudx) ) return 
      if (.not. computeEONVar( spec_rank,(/T2GM/),T2DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudy) ) return 
      if (.not. computeEONVar( spec_rank,(/T3GM/),T3DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudz) ) return 
      mseIrrSikmoinsSkjOikOjk1 = computeMSE(T1DNS,dudx,spec_rank,.true.)
      mseIrrSikmoinsSkjOikOjk2 = computeMSE(T2DNS,dudy,spec_rank,.true.)
      mseIrrSikmoinsSkjOikOjk3 = computeMSE(T3DNS,dudz,spec_rank,.true.)
      uvw(1)%values = - ( dudx%values * Ufil%values + dudy%values * Vfil%values + dudz%values * Wfil%values )
      dissIrrSikmoinsSkjOikOjk = computeFieldAvg(uvw(1),spec_rank)
      !Calcul de irr et quad RGMSIKmoinsSKJ
      call RGMSikmoinsSkj(T1RGMK,T2RGMK,T3RGMK,Ufil,Ukfil,Vkfil,Wkfil,WN,spec_rank,delta=delta)
      call btran(T1RGMK,T1GM,success)
      if (.not. success) return
      call btran(T2RGMK,T2GM,success)
      if (.not. success) return
      call btran(T3RGMK,T3GM,success)
      if (.not. success) return
      if (.not. computeEONVar( spec_rank,(/T1GM/),T1DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudx) ) return 
      if (.not. computeEONVar( spec_rank,(/T2GM/),T2DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudy) ) return 
      if (.not. computeEONVar( spec_rank,(/T3GM/),T3DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudz) ) return 
      mseIrrSikmoinsSkj1 = computeMSE(T1DNS,dudx,spec_rank,.true.)
      mseIrrSikmoinsSkj2 = computeMSE(T2DNS,dudy,spec_rank,.true.)
      mseIrrSikmoinsSkj3 = computeMSE(T3DNS,dudz,spec_rank,.true.)
      uvw(1)%values = - ( dudx%values * Ufil%values + dudy%values * Vfil%values + dudz%values * Wfil%values )
      dissIrrSikmoinsSkj = computeFieldAvg(uvw(1),spec_rank)
      mseQuadSikmoinsSkj1 = computeMSE(T1DNS,T1GM,spec_rank,.true.)
      mseQuadSikmoinsSkj2 = computeMSE(T2DNS,T2GM,spec_rank,.true.)
      mseQuadSikmoinsSkj3 = computeMSE(T3DNS,T3GM,spec_rank,.true.)
      uvw(1)%values = - ( T1GM%values * Ufil%values + T2GM%values * Vfil%values + T3GM%values * Wfil%values )
      dissSikmoinsSkjStat= computeFieldAvg(uvw(1),spec_rank)
      !Calcul du coef
      call RGMSikmoinsSkj(T1RGMK,T2RGMK,T3RGMK,Ufil,Ukfil,Vkfil,Wkfil,WN,spec_rank,delta=delta)
      call btran(T1RGMK,T1GM,success)
      if (.not. success) return
      call btran(T2RGMK,T2GM,success)
      if (.not. success) return
      call btran(T3RGMK,T3GM,success)
      if (.not. success) return
      uvw(1)%values = T1GM%values / T1DNS%values 
!      uvw(2)%values = T2GM%values / T2DNS%values 
!      uvw(3)%values = T3GM%values / T3DNS%values 
      if (.not.computeFieldPDF(1,uvw(1),51,spec_rank,10,minDelta=-2.0_WP,maxDelta=2.0_WP)) return
      !Calcul de irr RGMSIKSKJ-OikOjk
      call RGMOikOjk(T1RGMK,T2RGMK,T3RGMK,Ufil,Ukfil,Vkfil,Wkfil,WN,spec_rank,delta=delta)
      call btran(T1RGMK,T1RGM,success)
      if (.not. success) return
      call btran(T2RGMK,T2RGM,success)
      if (.not. success) return
      call btran(T3RGMK,T3RGM,success)
      if (.not. success) return
      if (.not. computeEONVar( spec_rank,(/T1GM,T1RGM/),T1DNS,&
                      & (/typeDisc,typeDisc/),(/bin,bin/),&
                      & fieldOpt = dudx) ) return 
      if (.not. computeEONVar( spec_rank,(/T2GM,T2RGM/),T2DNS,&
                      & (/typeDisc,typeDisc/),(/bin,bin/),&
                      & fieldOpt = dudy) ) return 
      if (.not. computeEONVar( spec_rank,(/T3GM,T3RGM/),T3DNS,&
                      & (/typeDisc,typeDisc/),(/bin,bin/),&
                      & fieldOpt = dudz) ) return 
      mseIrrSikmoinsSkj_OikOjk1 = computeMSE(T1DNS,dudx,spec_rank,.true.)
      mseIrrSikmoinsSkj_OikOjk2 = computeMSE(T2DNS,dudy,spec_rank,.true.)
      mseIrrSikmoinsSkj_OikOjk3 = computeMSE(T3DNS,dudz,spec_rank,.true.)
      uvw(1)%values = - ( dudx%values * Ufil%values + dudy%values * Vfil%values + dudz%values * Wfil%values )
      dissIrrSikmoinsSkj_OikOjk = computeFieldAvg(uvw(1),spec_rank)
      !Calcul de irr RGMOikOjk
      if (.not. computeEONVar( spec_rank,(/T1RGM/),T1DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudx) ) return 
      if (.not. computeEONVar( spec_rank,(/T2RGM/),T2DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudy) ) return 
      if (.not. computeEONVar( spec_rank,(/T3RGM/),T3DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudz) ) return 
      mseIrrOikOjk1 = computeMSE(T1DNS,dudx,spec_rank,.true.)
      mseIrrOikOjk2 = computeMSE(T2DNS,dudy,spec_rank,.true.)
      mseIrrOikOjk3 = computeMSE(T3DNS,dudz,spec_rank,.true.)
      uvw(1)%values = - ( dudx%values * Ufil%values + dudy%values * Vfil%values + dudz%values * Wfil%values )
      dissIrrOikOjk= computeFieldAvg(uvw(1),spec_rank)
      !Calcul de irr RGMSikmoinsSkjSikOjkOikSjk
      call RGMSikmoinsSkjSikOjkOikSjk(T1RGMK,T2RGMK,T3RGMK,Ufil,Ukfil,Vkfil,Wkfil,WN,spec_rank,delta=delta)
      call btran(T1RGMK,T1RGM,success)
      if (.not. success) return
      call btran(T2RGMK,T2RGM,success)
      if (.not. success) return
      call btran(T3RGMK,T3RGM,success)
      if (.not. success) return
      if (.not. computeEONVar( spec_rank,(/T1RGM/),T1DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudx) ) return 
      if (.not. computeEONVar( spec_rank,(/T2RGM/),T2DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudy) ) return 
      if (.not. computeEONVar( spec_rank,(/T3RGM/),T3DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudz) ) return 
      mseIrrSikmoinsSkjSikOjkOikSjk1 = computeMSE(T1DNS,dudx,spec_rank,.true.)
      mseIrrSikmoinsSkjSikOjkOikSjk2 = computeMSE(T2DNS,dudy,spec_rank,.true.)
      mseIrrSikmoinsSkjSikOjkOikSjk3 = computeMSE(T3DNS,dudz,spec_rank,.true.)
      uvw(1)%values = - ( dudx%values * Ufil%values + dudy%values * Vfil%values + dudz%values * Wfil%values )
      dissIrrSikmoinsSkjSikOjkOikSjk= computeFieldAvg(uvw(1),spec_rank)
      !Calcul de irr RGMSikOjkOikSjk
      call RGMSikOjkOikSjk(T1RGMK,T2RGMK,T3RGMK,Ufil,Ukfil,Vkfil,Wkfil,WN,spec_rank,delta=delta)
      call btran(T1RGMK,T1RGM,success)
      if (.not. success) return
      call btran(T2RGMK,T2RGM,success)
      if (.not. success) return
      call btran(T3RGMK,T3RGM,success)
      if (.not. success) return
      if (.not. computeEONVar( spec_rank,(/T1RGM/),T1DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudx) ) return 
      if (.not. computeEONVar( spec_rank,(/T2RGM/),T2DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudy) ) return 
      if (.not. computeEONVar( spec_rank,(/T3RGM/),T3DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudz) ) return 
      mseIrrSikOjkOikSjk1 = computeMSE(T1DNS,dudx,spec_rank,.true.)
      mseIrrSikOjkOikSjk2 = computeMSE(T2DNS,dudy,spec_rank,.true.)
      mseIrrSikOjkOikSjk3 = computeMSE(T3DNS,dudz,spec_rank,.true.)
      uvw(1)%values = - ( dudx%values * Ufil%values + dudy%values * Vfil%values + dudz%values * Wfil%values )
      dissIrrSikOjkOikSjk= computeFieldAvg(uvw(1),spec_rank)
      !Calcul de irr et quad RGMSIKmoinsSKJ
      !call RGMSIKmoinsSKJ(T1RGMK,T2RGMK,T3RGMK,Ufil,Vfil,Wfil,WN,spec_rank,delta=delta)
      call btran(T1RGMK,T1GM,success)
      if (.not. success) return
      call btran(T2RGMK,T2GM,success)
      if (.not. success) return
      call btran(T3RGMK,T3GM,success)
      if (.not. success) return
      if (.not. computeEONVar( spec_rank,(/T1GM,T1RGM/),T1DNS,&
                      & (/typeDisc,typeDisc/),(/bin,bin/),&
                      & fieldOpt = dudx) ) return 
      if (.not. computeEONVar( spec_rank,(/T2GM,T2RGM/),T2DNS,&
                      & (/typeDisc,typeDisc/),(/bin,bin/),&
                      & fieldOpt = dudy) ) return 
      if (.not. computeEONVar( spec_rank,(/T3GM,T3RGM/),T3DNS,&
                      & (/typeDisc,typeDisc/),(/bin,bin/),&
                      & fieldOpt = dudz) ) return 
      mseIrrSikmoinsSkj_SikOjkOikSjk1 = computeMSE(T1DNS,dudx,spec_rank,.true.)
      mseIrrSikmoinsSkj_SikOjkOikSjk2 = computeMSE(T2DNS,dudy,spec_rank,.true.)
      mseIrrSikmoinsSkj_SikOjkOikSjk3 = computeMSE(T3DNS,dudz,spec_rank,.true.)
      uvw(1)%values = - ( dudx%values * Ufil%values + dudy%values * Vfil%values + dudz%values * Wfil%values )
      dissIrrSikmoinsSkj_SikOjkOikSjk= computeFieldAvg(uvw(1),spec_rank)
      !Calcul de irr graddOkdxldUjdxl
      call graddOkdxldUjdxl(T1RGMK,T2RGMK,T3RGMK,Ufil,Vfil,Wfil,WN,spec_rank,delta=delta)
      call btran(T1RGMK,T1GM,success)
      if (.not. success) return
      call btran(T2RGMK,T2GM,success)
      if (.not. success) return
      call btran(T3RGMK,T3GM,success)
      if (.not. success) return
      if (.not. computeEONVar( spec_rank,(/T1GM/),T1DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudx) ) return 
      if (.not. computeEONVar( spec_rank,(/T2GM/),T2DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudy) ) return 
      if (.not. computeEONVar( spec_rank,(/T3GM/),T3DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudz) ) return 
      mseIrrgraddOkdxldUjdxl1 = computeMSE(T1DNS,dudx,spec_rank,.true.)
      mseIrrgraddOkdxldUjdxl2 = computeMSE(T2DNS,dudy,spec_rank,.true.)
      mseIrrgraddOkdxldUjdxl3 = computeMSE(T3DNS,dudz,spec_rank,.true.)
      uvw(1)%values = - ( dudx%values * Ufil%values + dudy%values * Vfil%values + dudz%values * Wfil%values )
      dissIrrgraddOkdxldUjdxl= computeFieldAvg(uvw(1),spec_rank)
      !Calcul de irr graddOkdxlSjlmoins 
      call graddOkdxlSjl(T1RGMK,T2RGMK,T3RGMK,Ufil,Vfil,Wfil,WN,spec_rank,'moins',delta=delta)
      call btran(T1RGMK,T1GM,success)
      if (.not. success) return
      call btran(T2RGMK,T2GM,success)
      if (.not. success) return
      call btran(T3RGMK,T3GM,success)
      if (.not. success) return
      if (.not. computeEONVar( spec_rank,(/T1GM/),T1DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudx) ) return 
      if (.not. computeEONVar( spec_rank,(/T2GM/),T2DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudy) ) return 
      if (.not. computeEONVar( spec_rank,(/T3GM/),T3DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudz) ) return 
      mseIrrgraddOkdxlSjlmoins1 = computeMSE(T1DNS,dudx,spec_rank,.true.)
      mseIrrgraddOkdxlSjlmoins2 = computeMSE(T2DNS,dudy,spec_rank,.true.)
      mseIrrgraddOkdxlSjlmoins3 = computeMSE(T3DNS,dudz,spec_rank,.true.)
      uvw(1)%values = - ( dudx%values * Ufil%values + dudy%values * Vfil%values + dudz%values * Wfil%values )
      dissIrrgraddOkdxlSjlmoins= computeFieldAvg(uvw(1),spec_rank)
      !Calcul de irr graddOkdxlSjlplus
      call graddOkdxlSjl(T1RGMK,T2RGMK,T3RGMK,Ufil,Vfil,Wfil,WN,spec_rank,'plus',delta=delta)
      call btran(T1RGMK,T1GM,success)
      if (.not. success) return
      call btran(T2RGMK,T2GM,success)
      if (.not. success) return
      call btran(T3RGMK,T3GM,success)
      if (.not. success) return
      if (.not. computeEONVar( spec_rank,(/T1GM/),T1DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudx) ) return 
      if (.not. computeEONVar( spec_rank,(/T2GM/),T2DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudy) ) return 
      if (.not. computeEONVar( spec_rank,(/T3GM/),T3DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudz) ) return 
      mseIrrgraddOkdxlSjlplus1 = computeMSE(T1DNS,dudx,spec_rank,.true.)
      mseIrrgraddOkdxlSjlplus2 = computeMSE(T2DNS,dudy,spec_rank,.true.)
      mseIrrgraddOkdxlSjlplus3 = computeMSE(T3DNS,dudz,spec_rank,.true.)
      uvw(1)%values = - ( dudx%values * Ufil%values + dudy%values * Vfil%values + dudz%values * Wfil%values )
      dissIrrgraddOkdxlSjlplus = computeFieldAvg(uvw(1),spec_rank)
      !Calcul de irr graddOkdxlOmegajl 
      call graddOkdxlOmegajl(T1RGMK,T2RGMK,T3RGMK,Ufil,Vfil,Wfil,WN,spec_rank,delta=delta)
      call btran(T1RGMK,T1GM,success)
      if (.not. success) return
      call btran(T2RGMK,T2GM,success)
      if (.not. success) return
      call btran(T3RGMK,T3GM,success)
      if (.not. success) return
      if (.not. computeEONVar( spec_rank,(/T1GM/),T1DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudx) ) return 
      if (.not. computeEONVar( spec_rank,(/T2GM/),T2DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudy) ) return 
      if (.not. computeEONVar( spec_rank,(/T3GM/),T3DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudz) ) return 
      mseIrrgraddOkdxlOmegajl1 = computeMSE(T1DNS,dudx,spec_rank,.true.)
      mseIrrgraddOkdxlOmegajl2 = computeMSE(T2DNS,dudy,spec_rank,.true.)
      mseIrrgraddOkdxlOmegajl3 = computeMSE(T3DNS,dudz,spec_rank,.true.)
      uvw(1)%values = - ( dudx%values * Ufil%values + dudy%values * Vfil%values + dudz%values * Wfil%values )
      dissIrrgraddOkdxlOmegajl = computeFieldAvg(uvw(1),spec_rank)
      !Calcul de quad RGMSikmoinsSkjSikOjkOikSjkDyn 
      call RGMSikmoinsSkjSikOjkOikSjkDyn(T1RGMK,T2RGMK,T3RGMK,Ufil,Vfil,Wfil,WN,spec_rank,filter,delta=delta)
      call btran(T1RGMK,T1GM,success)
      if (.not. success) return
      call btran(T2RGMK,T2GM,success)
      if (.not. success) return
      call btran(T3RGMK,T3GM,success)
      if (.not. success) return
      mseRGMSikmoinsSkjSikOjkOikSjkDyn1 = computeMSE(T1DNS,T1GM,spec_rank,.true.)
      mseRGMSikmoinsSkjSikOjkOikSjkDyn2 = computeMSE(T2DNS,T2GM,spec_rank,.true.)
      mseRGMSikmoinsSkjSikOjkOikSjkDyn3 = computeMSE(T3DNS,T3GM,spec_rank,.true.)
      uvw(1)%values = - ( T1GM%values * Ufil%values + T2GM%values * Vfil%values + T3GM%values * Wfil%values )
      dissIrrSikmoinsSkjSikOjkOikSjkDyn = computeFieldAvg(uvw(1),spec_rank)
      !Calcul de quad RGMSikmoinsSkjSikOjkOikSjkDyn2 
      call RGMSikmoinsSkjSikOjkOikSjkDyn2(T1RGMK,T2RGMK,T3RGMK,Ufil,Vfil,Wfil,WN,spec_rank,filter,delta=delta)
      call btran(T1RGMK,T1GM,success)
      if (.not. success) return
      call btran(T2RGMK,T2GM,success)
      if (.not. success) return
      call btran(T3RGMK,T3GM,success)
      if (.not. success) return
      mseRGMSikmoinsSkjSikOjkOikSjkDyn21 = computeMSE(T1DNS,T1GM,spec_rank,.true.)
      mseRGMSikmoinsSkjSikOjkOikSjkDyn22 = computeMSE(T2DNS,T2GM,spec_rank,.true.)
      mseRGMSikmoinsSkjSikOjkOikSjkDyn23 = computeMSE(T3DNS,T3GM,spec_rank,.true.)
      uvw(1)%values = - ( T1GM%values * Ufil%values + T2GM%values * Vfil%values + T3GM%values * Wfil%values )
      dissIrrSikmoinsSkjSikOjkOikSjkDyn2 = computeFieldAvg(uvw(1),spec_rank)
      !Calcul de quad RGMSikmoinsSkjDyn2 
      call RGMSikmoinsSkjDyn(T1RGMK,T2RGMK,T3RGMK,Ufil,Vfil,Wfil,WN,spec_rank,filter,delta=delta)
      call btran(T1RGMK,T1GM,success)
      if (.not. success) return
      call btran(T2RGMK,T2GM,success)
      if (.not. success) return
      call btran(T3RGMK,T3GM,success)
      if (.not. success) return
      mseRGMSikmoinsSkjDyn1 = computeMSE(T1DNS,T1GM,spec_rank,.true.)
      mseRGMSikmoinsSkjDyn2 = computeMSE(T2DNS,T2GM,spec_rank,.true.)
      mseRGMSikmoinsSkjDyn3 = computeMSE(T3DNS,T3GM,spec_rank,.true.)
      uvw(1)%values = - ( T1GM%values * Ufil%values + T2GM%values * Vfil%values + T3GM%values * Wfil%values )
      dissSikmoinsSkjDyn = computeFieldAvg(uvw(1),spec_rank)
      !Calcul de quad RGMSikmoinsSkjDynFab 
      call RGMSikmoinsSkjDynFab(T1RGMK,T2RGMK,T3RGMK,Ufil,Ukfil,Vkfil,Wkfil,WN,spec_rank,filter,delta=delta)
      call btran(T1RGMK,T1GM,success)
      if (.not. success) return
      call btran(T2RGMK,T2GM,success)
      if (.not. success) return
      call btran(T3RGMK,T3GM,success)
      if (.not. success) return
      mseRGMSikmoinsSkjDynFab1 = computeMSE(T1DNS,T1GM,spec_rank,.true.)
      mseRGMSikmoinsSkjDynFab2 = computeMSE(T2DNS,T2GM,spec_rank,.true.)
      mseRGMSikmoinsSkjDynFab3 = computeMSE(T3DNS,T3GM,spec_rank,.true.)
      uvw(1)%values = - ( T1GM%values * Ufil%values + T2GM%values * Vfil%values + T3GM%values * Wfil%values )
      dissSikmoinsSkjDynFab = computeFieldAvg(uvw(1),spec_rank)
      !Calcul de quad RGMSikmoinsSkjSikOjkOikSjkDynFab  
      call RGMSikmoinsSkjSikOjkOikSjkDynFab(T1RGMK,T2RGMK,T3RGMK,Ufil,Vfil,Wfil,WN,spec_rank,filter,delta=delta)
      call btran(T1RGMK,T1GM,success)
      if (.not. success) return
      call btran(T2RGMK,T2GM,success)
      if (.not. success) return
      call btran(T3RGMK,T3GM,success)
      if (.not. success) return
      mseSikmoinsSkjSikOjkOikSjkDynFab1 = computeMSE(T1DNS,T1GM,spec_rank,.true.)
      mseSikmoinsSkjSikOjkOikSjkDynFab2 = computeMSE(T2DNS,T2GM,spec_rank,.true.)
      mseSikmoinsSkjSikOjkOikSjkDynFab3 = computeMSE(T3DNS,T3GM,spec_rank,.true.)
      uvw(1)%values = - ( T1GM%values * Ufil%values + T2GM%values * Vfil%values + T3GM%values * Wfil%values )
      dissSikmoinsSkjSikOjkOikSjkDynFab = computeFieldAvg(uvw(1),spec_rank)
      !Calcul de quad RGMSikmoinsSkjSikOjkOikSjkDyn2Fab
      call RGMSikmoinsSkjSikOjkOikSjkDyn2Fab(T1RGMK,T2RGMK,T3RGMK,Ufil,Vfil,Wfil,WN,spec_rank,filter,delta=delta)
      call btran(T1RGMK,T1GM,success)
      if (.not. success) return
      call btran(T2RGMK,T2GM,success)
      if (.not. success) return
      call btran(T3RGMK,T3GM,success)
      if (.not. success) return
      mseSikmoinsSkjSikOjkOikSjkDyn2Fab1 = computeMSE(T1DNS,T1GM,spec_rank,.true.)
      mseSikmoinsSkjSikOjkOikSjkDyn2Fab2 = computeMSE(T2DNS,T2GM,spec_rank,.true.)
      mseSikmoinsSkjSikOjkOikSjkDyn2Fab3 = computeMSE(T3DNS,T3GM,spec_rank,.true.)
      uvw(1)%values = - ( T1GM%values * Ufil%values + T2GM%values * Vfil%values + T3GM%values * Wfil%values )
      dissSikmoinsSkjSikOjkOikSjkDyn2Fab = computeFieldAvg(uvw(1),spec_rank)

      !diffExact= - < d/dxj(tauij ui)>
         ! computeT_ijVecA(in1,in2,in3,......,wave,out11,out12,out13,out22,out23,out33,filter,success,delta)
      call computeT_ijVecA(U,V,W,Ufil,Vfil,Wfil,WN,T1RGM,T2RGM,T3RGM,T1GM ,T2GM ,T3GM ,filter,success,delta)
      T1DNS%values=Ufil%values * T1RGM%values
      T2DNS%values=Ufil%values * T2RGM%values
      T3DNS%values=Ufil%values * T3RGM%values
      if (.not. computeFieldDivergence(WN,T1DNS,T2DNS,T3DNS,dudx,nbcpus,spec_rank)) return
      diffExactComp(1)= - computeFieldAvg(dudx,spec_rank)
      T1DNS%values=Vfil%values * T2RGM%values
      T2DNS%values=Vfil%values * T1GM%values
      T3DNS%values=Vfil%values * T2GM%values
      if (.not. computeFieldDivergence(WN,T1DNS,T2DNS,T3DNS,dudy,nbcpus,spec_rank)) return
      diffExactComp(2)= - computeFieldAvg(dudy,spec_rank)
      T1DNS%values=Wfil%values * T3RGM%values
      T2DNS%values=Wfil%values * T2GM%values
      T3DNS%values=Wfil%values * T3GM%values
      if (.not. computeFieldDivergence(WN,T1DNS,T2DNS,T3DNS,dudz,nbcpus,spec_rank)) return
      diffExactComp(3)= - computeFieldAvg(dudz,spec_rank)
      diffExact       = diffExactComp(1) + diffExactComp(2) + diffExactComp(3) 

      !dissExact=  < tauij dui/dxj >
      if (.not. computeFieldGradient(WN,Ufil,dudx,dudy,dudz,nbcpus,spec_rank) ) return
      T1RGM%values=T1RGM%values*dudx%values + T2RGM%values*dudy%values + T3RGM%values*dudz%values
      dissExactComp(1)=   computeFieldAvg(T1RGM,spec_rank)
      if (.not. computeFieldGradient(WN,Vfil,dudx,dudy,dudz,nbcpus,spec_rank) ) return
      T2RGM%values=T2RGM%values*dudx%values + T1GM%values*dudy%values  + T2GM%values*dudz%values
      dissExactComp(2)=   computeFieldAvg(T2RGM,spec_rank)
      if (.not. computeFieldGradient(WN,Wfil,dudx,dudy,dudz,nbcpus,spec_rank) ) return
      T3RGM%values=T3RGM%values*dudx%values + T2GM%values*dudy%values  + T3GM%values*dudz%values
      dissExactComp(3)=   computeFieldAvg(T3RGM,spec_rank)
      T1RGM%values= T1RGM%values + T2RGM%values + T3RGM%values
      dissExact=computeFieldAvg(T1RGM,spec_rank)

      uvw(1)%values=Ufil%values
      uvw(2)%values=Vfil%values
      uvw(3)%values=Wfil%values
      !diffGM   = - < C Delta^2 d/dxj( duj/dxk dui/dxk ui ) >
      !         = - < C Delta^2 duj/dxk ( d^2ui/dxjdxk ui + dui/dxk dui/dxj ) >
      diffGM=0.0_WP
      do i=1,3
        Vfil%values=0.0_WP
        do j=1,3
          do k=1,3
            if (.not. computeDerivationField(uvw(j),(/k/),T2DNS,spec_rank,WN)) return
            if (.not. computeDerivationField(uvw(i),(/k/),T3DNS,spec_rank,WN)) return
            T2DNS%values= T2DNS%values * T3DNS%values*uvw(i)%values
            if (.not. computeDerivationField(T2DNS,(/j/),T1DNS,spec_rank,WN)) return
            Vfil%values = Vfil%values + T1DNS%values
          enddo
        enddo
        diffGMComp(i) = - delta**2/12.0_WP*computeFieldAvg(Vfil,spec_rank)
        diffGM = diffGM + diffGMComp(i) 
      enddo

      !dissGM =  < C Delta^2  duj/dxk dui/dxk dui/dxj  >
      dissGM=0.0_WP
      do i=1,3
        Vfil%values=0.0_WP
        do j=1,3
          if (.not. computeDerivationField(uvw(i),(/j/),T1DNS,spec_rank,WN)) return
          do k=1,3
            if (.not. computeDerivationField(uvw(j),(/k/),T2DNS,spec_rank,WN)) return
            if (.not. computeDerivationField(uvw(i),(/k/),T3DNS,spec_rank,WN)) return
            Vfil%values = T2DNS%values * T3DNS%values * T1DNS%values  + Vfil%values
          enddo
        enddo
        dissGMComp(i) = delta**2/12.0_WP*computeFieldAvg(Vfil,spec_rank)
        dissGM = dissGM + dissGMComp(i) 
      enddo


      !echange < d/dxj tau_ij ui > pour RGM
      call RGMVel(T1RGMK,T2RGMk,T3RGMk,U,uvw(1),uvw(2),uvw(3),WN,spec_rank,'moins',delta=delta)
      call btran(T1RGMK,T1RGM,success)
      if (.not. success) return
      call btran(T2RGMK,T2RGM,success)
      if (.not. success) return
      call btran(T3RGMK,T3RGM,success)
      if (.not. success) return
      call ftran(uvw(1),T1RGMK,success) ! -> Ukbar
      if (.not. success) return
      call ftran(uvw(2),T2RGMK,success) ! -> Vkbar
      if (.not. success) return
      call ftran(uvw(3),T3RGMK,success) ! -> Wkbar
      if (.not. success) return
      call computeStressTensor(T1RGMK,T2RGMK,T3RGMK,WN,T1GM,T1RGM,T2RGM,T2GM,T3RGM,T3GM,res)
      do k=U%zmin,U%zmax
        do j=U%ymin,U%ymax
          do i=U%xmin,U%xmax
            if (.not. computeSplitingSij(T1GM%values(i,j,k),T1RGM%values(i,j,k),&
                                        &T2RGM%values(i,j,k),T2GM%values(i,j,k),&
                                        &T3RGM%values(i,j,k),T3GM%values(i,j,k),&
                                        &T1GMminus,T1RGMminus,T2RGMminus,T2GMminus,T3RGMminus,T3GMminus,&
                                        &"moins")) print *,'[ERROR]' 
            T1GM%values(i,j,k) = T1GMminus   !S11minus
            T2GM%values(i,j,k) = T2GMminus   !S22minus
            T3GM%values(i,j,k) = T3GMminus   !S33minus
            T1RGM%values(i,j,k) = T1RGMminus !S12minus
            T2RGM%values(i,j,k) = T2RGMminus !S13minus
            T3RGM%values(i,j,k) = T3RGMminus !S23minus
          enddo
        enddo
      enddo

      !diffRGM  = - < C Delta^2 Sjk^moins d/dxj( dui/dxk ui  )>
      diffRGM=0.0_WP
      do i=1,3
        Vfil%values=0.0_WP
        if (.not. computeFieldGradient(WN,uvw(i),dudx,dudy,dudz,nbcpus,spec_rank) ) return
        !S11^minus*d/dx(dui/dx ui)
        T1DNS%values=dudx%values*uvw(i)%values
        if (.not. computeDerivationField(T1DNS,(/1/),T3DNS,spec_rank,WN)) return
!        if (.not. computeFieldGradient(WN,T1DNS,T3DNS,Ufil,Ufil,nbcpus,spec_rank) ) return
        Vfil%values=T1GM%values*T3DNS%values + Vfil%values
        !S22^minus*d/dy(dui/dy ui)
        T1DNS%values=dudy%values*uvw(i)%values
        if (.not. computeDerivationField(T1DNS,(/2/),T3DNS,spec_rank,WN)) return
!        if (.not. computeFieldGradient(WN,T1DNS,Ufil,T3DNS,Ufil,nbcpus,spec_rank) ) return
        Vfil%values=T2GM%values*T3DNS%values + Vfil%values
        !S33^minus*d/dz(dui/dz ui)
        T1DNS%values=dudz%values*uvw(i)%values
        if (.not. computeDerivationField(T1DNS,(/3/),T3DNS,spec_rank,WN)) return
!        if (.not. computeFieldGradient(WN,T1DNS,Ufil,Ufil,T3DNS,nbcpus,spec_rank) ) return
        Vfil%values=T3GM%values*T3DNS%values + Vfil%values

        !2*S12^minus*(d^2ui/dxdy ui + dui/dx dui/dy)
        if (.not. computeDerivationField(uvw(i),(/1,2/),T3DNS,spec_rank,WN)) return
        T1DNS%values=2.0_WP*T1RGM%values* (uvw(i)%values*T3DNS%values + dudy%values * dudx%values )
        Vfil%values=Vfil%values + T1DNS%values
        !2*S13^minus*(d^2ui/dxdz ui + dui/dx dui/dz)
        if (.not. computeDerivationField(uvw(i),(/1,3/),T3DNS,spec_rank,WN)) return
        T1DNS%values=2.0_WP*T2RGM%values* (uvw(i)%values*T3DNS%values + dudx%values * dudz%values )
        Vfil%values=Vfil%values + T1DNS%values
        !2*S23^minus*(d^2ui/dzdy ui + dui/dz dui/dy)
        if (.not. computeDerivationField(uvw(i),(/2,3/),T3DNS,spec_rank,WN)) return
        T1DNS%values=2.0_WP*T3RGM%values* (uvw(i)%values*T3DNS%values + dudy%values * dudz%values )
        Vfil%values=Vfil%values + T1DNS%values
        diffRGMComp(i) = - delta**2/12.0_WP * computeFieldAvg(Vfil,spec_rank)
        diffRGM = diffRGMComp(i) + diffRGM
      enddo

      !dissRGM  =  < C Delta^2 Sjk^moins dui/dxk dui/dxj  >
      !         =    C Delta^2 < S11^moins (dui/dx)^2 + S22^moins (dui/dy)^2 + S33^moins (dui/dz)^2 
      !           +  2 (S12^moins dui/dx dui/dy + S13^moins dui/dx dui/dz + S23^moins dui/dy dui/dz)
      dissRGM=0.0_WP
      do i=1,3
        if (.not. computeFieldGradient(WN,uvw(i),dudx,dudy,dudz,nbcpus,spec_rank) ) return
        Vfil%values = T1GM%values * dudx%values**2 + T2GM%values * dudy%values**2 + T3GM%values * dudz%values**2 + &
                    & 2.0_WP*T1RGM%values * dudx%values*dudy%values +& 
                    & 2.0_WP*T2RGM%values * dudx%values*dudz%values +&
                    & 2.0_WP*T3RGM%values * dudy%values*dudz%values
        dissRGMComp(i) = delta**2/12.0_WP * computeFieldAvg(Vfil,spec_rank)
        dissRGM = dissRGMComp(i) + dissRGM
      enddo

      !dump
      if (spec_rank .EQ. 0) write (6,'(a)')                         '################    ERREUR IRREDUCTIBLE    #######################' 
      if (spec_rank .EQ. 0) write (6,'(a)')                         '###########           Modele de base       #######################' 
      if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO IRR]  GM de base                  :',nd,mseIrrGM1,mseIrrGM2,mseIrrGM3
      if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO IRR]  SMAG                        :',nd,mseIrrSmag1,mseIrrSmag2,mseIrrSmag3
      if (spec_rank .EQ. 0) write (6,'(a)')                         '###########Modele sur la divergence de t_ij#######################' 
      if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO IRR]  RGM -                       :',nd,mseIrrRGM1,mseIrrRGM2,mseIrrRGM3
      if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO IRR]  RGM +/- 2vars               :',nd,mseIrrRGM1_2var,mseIrrRGM2_2var,mseIrrRGM3_2var
      if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO IRR]  RGM +                       :',nd,mseIrrRGM1_plus,mseIrrRGM2_plus,mseIrrRGM3_plus
      if (spec_rank .EQ. 0) write (6,'(a)')                         '###########       Modele sur la vorticite  #######################' 
      if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO IRR]  VORT                        :',nd,mseIrrGMVort1,mseIrrGMVort2,mseIrrGMVort3
      if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO IRR] graddOkdxldUjdxl             :',nd,mseIrrgraddOkdxldUjdxl1,mseIrrgraddOkdxldUjdxl2,mseIrrgraddOkdxldUjdxl3
      if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO IRR] graddOkdxlSjl-               :',nd,mseIrrgraddOkdxlSjlmoins1,mseIrrgraddOkdxlSjlmoins2,&
                                                                                                             &mseIrrgraddOkdxlSjlmoins3
      if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO IRR] graddOkdxlSjl+               :',nd,mseIrrgraddOkdxlSjlplus1,mseIrrgraddOkdxlSjlplus2,mseIrrgraddOkdxlSjlplus3
      if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO IRR] graddOkdxlOmegajl            :',nd,mseIrrgraddOkdxlOmegajl1,mseIrrgraddOkdxlOmegajl2,mseIrrgraddOkdxlOmegajl3
      if (spec_rank .EQ. 0) write (6,'(a)')                         '###########       Modele sur t_ij          #######################' 
      if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO IRR] Sik-Skj                      :',nd,mseIrrSikmoinsSkj1,mseIrrSikmoinsSkj2,mseIrrSikmoinsSkj3
      if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO IRR] Sik-SkjOikOjk                :',nd,mseIrrSikmoinsSkjOikOjk1,mseIrrSikmoinsSkjOikOjk2,mseIrrSikmoinsSkjOikOjk3
      if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO IRR] Sik-Skj_OikOjk               :',nd,mseIrrSikmoinsSkj_OikOjk1,mseIrrSikmoinsSkj_OikOjk2,&
                                                                                                             &mseIrrSikmoinsSkj_OikOjk3
      if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO IRR] OikOjk                       :',nd,mseIrrOikOjk1,mseIrrOikOjk2,mseIrrOikOjk3
      if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO IRR] Sik-SkjSikOjkOikSjk          :',nd,mseIrrSikmoinsSkjSikOjkOikSjk1,mseIrrSikmoinsSkjSikOjkOikSjk2,&
                                                                                                             &mseIrrSikmoinsSkjSikOjkOikSjk3
      if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO IRR] Sik-Skj_SikOjkOikSjk         :',nd,mseIrrSikmoinsSkj_SikOjkOikSjk1,mseIrrSikmoinsSkj_SikOjkOikSjk2,&
                                                                                                             &mseIrrSikmoinsSkj_SikOjkOikSjk3
      if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO IRR] SikOjkOikSjk                 :',nd,mseIrrSikOjkOikSjk1,mseIrrSikOjkOikSjk2,mseIrrSikOjkOikSjk3
      if (spec_rank .EQ. 0) write (6,'(a)')                         '################    ERREUR QUADRATIQUE     #######################' 
      if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO QUAD] GM                          :',nd,mseQuadGM1,mseQuadGM2,mseQuadGM3
      if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO QUAD] RGM -                       :',nd,mseQuadRGM1,mseQuadRGM2,mseQuadRGM3
      if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO QUAD] RGM +                       :',nd,mseQuadRGM1_plus,mseQuadRGM2_plus,mseQuadRGM3_plus
      if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO QUAD] SMAG dyn                    :',nd,mseQuadSmag1,mseQuadSmag2,mseQuadSmag3
      if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO QUAD] SMAG stat                   :',nd,mseQuadSmagStat1,mseQuadSmagStat2,mseQuadSmagStat3
      if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO QUAD] RGMdyn -                    :',nd,mseQuadRGM1_2var,mseQuadRGM2_2var,mseQuadRGM3_2var
      if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO QUAD] RGMdyn+/-                   :',nd,mseQuadRGM1dyn,mseQuadRGM2dyn,mseQuadRGM3dyn
      if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO QUAD] Vort                        :',nd,mseQuadGMVort1,mseQuadGMVort2,mseQuadGMVort3
      if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO QUAD] Sik-Skj    Statique         :',nd,mseQuadSikmoinsSkj1,mseQuadSikmoinsSkj2,mseQuadSikmoinsSkj3
      if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO QUAD] Sik-Skj    Dyn              :',nd,mseRGMSikmoinsSkjDyn1,mseRGMSikmoinsSkjDyn2,mseRGMSikmoinsSkjDyn3
      if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO QUAD] Sik-Skj    DynFab           :',nd,mseRGMSikmoinsSkjDynFab1,mseRGMSikmoinsSkjDynFab2,mseRGMSikmoinsSkjDynFab3
      if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO QUAD] Sik-SkjSikOjkOikSjk dyn1    :',nd,mseRGMSikmoinsSkjSikOjkOikSjkDyn1,mseRGMSikmoinsSkjSikOjkOikSjkDyn2,&
                                                                                                             &mseRGMSikmoinsSkjSikOjkOikSjkDyn3
      if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO QUAD] Sik-SkjSikOjkOikSjk dyn2    :',nd,mseRGMSikmoinsSkjSikOjkOikSjkDyn21,mseRGMSikmoinsSkjSikOjkOikSjkDyn22,&
                                                                                                             &mseRGMSikmoinsSkjSikOjkOikSjkDyn23
      if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO QUAD] Sik-SkjSikOjkOikSjk dyn1Fab :',nd,mseSikmoinsSkjSikOjkOikSjkDynFab1,mseSikmoinsSkjSikOjkOikSjkDynFab2,&
                                                                                                             &mseSikmoinsSkjSikOjkOikSjkDynFab3
      if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO QUAD] Sik-SkjSikOjkOikSjk dynFab2 :',nd,mseSikmoinsSkjSikOjkOikSjkDyn2Fab1,mseSikmoinsSkjSikOjkOikSjkDyn2Fab2,&
                                                                                                             &mseSikmoinsSkjSikOjkOikSjkDyn2Fab3
      if (spec_rank .EQ. 0) write (6,'(a)')                         '################    DIFFUSION DISSIPATION  #######################' 
      if (spec_rank .EQ. 0) write (6,'(3(a,1x,i0,1x,2(g15.8,1x)))') '[INFO diff + diss] GM          :',nd,diffdissGM,diffGM+dissGM
      if (spec_rank .EQ. 0) write (6,'(3(a,1x,i0,1x,2(g15.8,1x)))') '[INFO diff + diss]RGM          :',nd,diffdissRGM,diffRGM+dissRGM
      if (spec_rank .EQ. 0) write (6,'(3(a,1x,i0,1x,2(g15.8,1x)))') '[INFO diff + diss]DNS          :',nd,diffdissExacte,diffExact+dissExact
      if (spec_rank .EQ. 0) write (6,'(a,1x,i0,1x,4(g15.8,1x))')    '[INFO diff]     DNS            :',nd,diffExactComp(1),diffExactComp(2),diffExactComp(3),diffExact
      if (spec_rank .EQ. 0) write (6,'(a,1x,i0,1x,4(g15.8,1x))')    '[INFO diss]     DNS            :',nd,dissExactComp(1),dissExactComp(2),dissExactComp(3),dissExact
      if (spec_rank .EQ. 0) write (6,'(a,1x,i0,1x,4(g15.8,1x))')    '[INFO diff]     RGM            :',nd,diffRGMComp(1),diffRGMComp(2),diffRGMComp(3),diffRGM
      if (spec_rank .EQ. 0) write (6,'(a,1x,i0,1x,4(g15.8,1x))')    '[INFO diss]     RGM            :',nd,dissRGMComp(1),dissRGMComp(2),dissRGMComp(3),dissRGM
      if (spec_rank .EQ. 0) write (6,'(a,1x,i0,1x,4(g15.8,1x))')    '[INFO diff]      GM            :',nd,diffGMComp(1),diffGMComp(2),diffGMComp(3),diffGM
      if (spec_rank .EQ. 0) write (6,'(a,1x,i0,1x,4(g15.8,1x))')    '[INFO diss]      GM            :',nd,dissGMComp(1),dissGMComp(2),dissGMComp(3),dissGM
      if (spec_rank .EQ. 0) then
        file_id=iopen()
        inquire(file='bilanMSE_GM_VS_RGM.table', exist=existe)
        open(unit=file_id, file='bilanMSE_GM_VS_RGM.table', form="FORMATTED",position="append")
        if ( .not. existe) then
          write(file_id, '(a)') '#iFic nd mseQuadGM1 mseQuadGM2 mseQuadGM3 &
                                         &mseQuadRGM1 mseQuadRGM2 mseQuadRGM3 &
                                         &mseIrrGM1 mseIrrGM2 mseIrrGM3 &
                                         &mseIrrRGM1 mseIrrRGM2 mseIrrRGM3 &
                                         &mseIrrRGM1_2var mseIrrRGM2_2var mseIrrRGM3_2var &
                                         &mseQuadRGM1_plus mseQuadRGM2_plus mseQuadRGM3_plus &
                                         &mseIrrRGM1_plus mseIrrRGM2_plus mseIrrRGM3_plus &
                                         &mseQuadSmag1 mseQuadSmag2 mseQuadSmag3 &
                                         &mseIrrSmag1 mseIrrSmag2 mseIrrSmag3 &
                                         &mseQuadSmagStat1 mseQuadSmagStat2 mseQuadSmagStat3 &
                                         &mseQuadRGM1_2var mseQuadRGM2_2var mseQuadRGM3_2var &
                                         &mseQuadRGM1dyn mseQuadRGM2dyn mseQuadRGM3dyn &
                                         &mseQuadGMVort1 mseQuadGMVort2 mseQuadGMVort3 &
                                         &mseIrrGMVort1 mseIrrGMVort2 mseIrrGMVort3 &
                                         &mseIrrSikmoinsSkj1 mseIrrSikmoinsSkj2 mseIrrSikmoinsSkj3 &
                                         &mseIrrSikmoinsSkjOikOjk1 mseIrrSikmoinsSkjOikOjk2 mseIrrSikmoinsSkjOikOjk3 &
                                         &mseIrrSikmoinsSkj_OikOjk1 mseIrrSikmoinsSkj_OikOjk2 mseIrrSikmoinsSkj_OikOjk3 &
                                         &mseIrrOikOjk1 mseIrrOikOjk2 mseIrrOikOjk3 &
                                         &mseIrrgraddOkdxldUjdxl1 mseIrrgraddOkdxldUjdxl2 mseIrrgraddOkdxldUjdxl3 &
                                         &mseIrrgraddOkdxlSjlmoins1 mseIrrgraddOkdxlSjlmoins2 mseIrrgraddOkdxlSjlmoins3 &
                                         &mseIrrgraddOkdxlSjlplus1 mseIrrgraddOkdxlSjlplus2 mseIrrgraddOkdxlSjlplus3 &
                                         &mseIrrgraddOkdxlOmegajl1 mseIrrgraddOkdxlOmegajl2 mseIrrgraddOkdxlOmegajl3 &
                                         &mseRGMSikmoinsSkjSikOjkOikSjkDyn1 mseRGMSikmoinsSkjSikOjkOikSjkDyn2 mseRGMSikmoinsSkjSikOjkOikSjkDyn3 &
                                         &mseIrrSikmoinsSkjSikOjkOikSjk1 mseIrrSikmoinsSkjSikOjkOikSjk2 mseIrrSikmoinsSkjSikOjkOikSjk3 &
                                         &mseIrrSikmoinsSkj_SikOjkOikSjk1 mseIrrSikmoinsSkj_SikOjkOikSjk2 mseIrrSikmoinsSkj_SikOjkOikSjk3 &
                                         &mseIrrSikOjkOikSjk1 mseIrrSikOjkOikSjk2 mseIrrSikOjkOikSjk3 &
                                         &mseQuadSikmoinsSkj1 mseQuadSikmoinsSkj2 mseQuadSikmoinsSkj3 &
                              &mseRGMSikmoinsSkjSikOjkOikSjkDyn21 mseRGMSikmoinsSkjSikOjkOikSjkDyn22 mseRGMSikmoinsSkjSikOjkOikSjkDyn23 &
                                              &mseRGMSikmoinsSkjDyn1 mseRGMSikmoinsSkjDyn2 mseRGMSikmoinsSkjDyn3'
        endif
        write(file_id,'(i0,1x,i0,1x,78(g15.8,1x))') iFic,nd,mseQuadGM1,mseQuadGM2,mseQuadGM3,&
                                                          &mseQuadRGM1,mseQuadRGM2,mseQuadRGM3,&
                                                          &mseIrrGM1,mseIrrGM2,mseIrrGM3,&
                                                          &mseIrrRGM1,mseIrrRGM2,mseIrrRGM3,&
                                                          &mseIrrRGM1_2var,mseIrrRGM2_2var,mseIrrRGM3_2var,&
                                                          &mseQuadRGM1_plus,mseQuadRGM2_plus,mseQuadRGM3_plus,&
                                                          &mseIrrRGM1_plus,mseIrrRGM2_plus,mseIrrRGM3_plus,&
                                                          &mseQuadSmag1,mseQuadSmag2,mseQuadSmag3,&
                                                          &mseIrrSmag1,mseIrrSmag2,mseIrrSmag3,&
                                                          &mseQuadSmagStat1,mseQuadSmagStat2,mseQuadSmagStat3,&
                                                          &mseQuadRGM1_2var,mseQuadRGM2_2var,mseQuadRGM3_2var,&
                                                          &mseQuadRGM1dyn,mseQuadRGM2dyn,mseQuadRGM3dyn,&
                                                          &mseQuadGMVort1,mseQuadGMVort2,mseQuadGMVort3,&
                                                          &mseIrrGMVort1,mseIrrGMVort2,mseIrrGMVort3,&
                                                          &mseIrrSikmoinsSkj1,mseIrrSikmoinsSkj2,mseIrrSikmoinsSkj3,&
                                                          &mseIrrSikmoinsSkjOikOjk1,mseIrrSikmoinsSkjOikOjk2,mseIrrSikmoinsSkjOikOjk3,&
                                                          &mseIrrSikmoinsSkj_OikOjk1,mseIrrSikmoinsSkj_OikOjk2,mseIrrSikmoinsSkj_OikOjk3,&
                                                          &mseIrrOikOjk1,mseIrrOikOjk2,mseIrrOikOjk3,&
                                                          &mseIrrgraddOkdxldUjdxl1,mseIrrgraddOkdxldUjdxl2,mseIrrgraddOkdxldUjdxl3,&
                                                          &mseIrrgraddOkdxlSjlmoins1,mseIrrgraddOkdxlSjlmoins2,mseIrrgraddOkdxlSjlmoins3,&
                                                          &mseIrrgraddOkdxlSjlplus1,mseIrrgraddOkdxlSjlplus2,mseIrrgraddOkdxlSjlplus3,&
                                                          &mseIrrgraddOkdxlOmegajl1,mseIrrgraddOkdxlOmegajl2,mseIrrgraddOkdxlOmegajl3,&
                                                     &mseRGMSikmoinsSkjSikOjkOikSjkDyn1,mseRGMSikmoinsSkjSikOjkOikSjkDyn2,mseRGMSikmoinsSkjSikOjkOikSjkDyn3,&
                                                 &mseIrrSikmoinsSkjSikOjkOikSjk1,mseIrrSikmoinsSkjSikOjkOikSjk2,mseIrrSikmoinsSkjSikOjkOikSjk3,&
                                                 &mseIrrSikmoinsSkj_SikOjkOikSjk1,mseIrrSikmoinsSkj_SikOjkOikSjk2,mseIrrSikmoinsSkj_SikOjkOikSjk3,&
                                                          &mseIrrSikOjkOikSjk1,mseIrrSikOjkOikSjk2,mseIrrSikOjkOikSjk3,&
                                                          &mseQuadSikmoinsSkj1,mseQuadSikmoinsSkj2,mseQuadSikmoinsSkj3,&
                                                &mseRGMSikmoinsSkjSikOjkOikSjkDyn21,mseRGMSikmoinsSkjSikOjkOikSjkDyn22,mseRGMSikmoinsSkjSikOjkOikSjkDyn23,&
                                              &mseRGMSikmoinsSkjDyn1,mseRGMSikmoinsSkjDyn2,mseRGMSikmoinsSkjDyn3
        close(file_id)
        file_id=iclose(file_id)

        file_id=iopen()
        inquire(file='bilanMSEIRR.table', exist=existe)
        open(unit=file_id, file='bilanMSEIRR.table', form="FORMATTED",position="append")
        if ( .not. existe) then
          write(file_id, '(a)') '#iFic nd mseIrrGM1 mseIrrGM2 mseIrrGM3 &
                                         &mseIrrRGM1 mseIrrRGM2 mseIrrRGM3 &
                                         &mseIrrRGM1_2var mseIrrRGM2_2var mseIrrRGM3_2var &
                                         &mseIrrRGM1_plus mseIrrRGM2_plus mseIrrRGM3_plus &
                                         &mseIrrSmag1 mseIrrSmag2 mseIrrSmag3 &
                                         &mseIrrGMVort1 mseIrrGMVort2 mseIrrGMVort3 &
                                         &mseIrrSikmoinsSkj1 mseIrrSikmoinsSkj2 mseIrrSikmoinsSkj3 &
                                         &mseIrrSikmoinsSkjOikOjk1 mseIrrSikmoinsSkjOikOjk2 mseIrrSikmoinsSkjOikOjk3 &
                                         &mseIrrSikmoinsSkj_OikOjk1 mseIrrSikmoinsSkj_OikOjk2 mseIrrSikmoinsSkj_OikOjk3 &
                                         &mseIrrOikOjk1 mseIrrOikOjk2 mseIrrOikOjk3 &
                                         &mseIrrgraddOkdxldUjdxl1 mseIrrgraddOkdxldUjdxl2 mseIrrgraddOkdxldUjdxl3 &
                                         &mseIrrgraddOkdxlSjlmoins1 mseIrrgraddOkdxlSjlmoins2 mseIrrgraddOkdxlSjlmoins3 &
                                         &mseIrrgraddOkdxlSjlplus1 mseIrrgraddOkdxlSjlplus2 mseIrrgraddOkdxlSjlplus3 &
                                         &mseIrrgraddOkdxlOmegajl1 mseIrrgraddOkdxlOmegajl2 mseIrrgraddOkdxlOmegajl3 &
                                         &mseIrrSikmoinsSkjSikOjkOikSjk1 mseIrrSikmoinsSkjSikOjkOikSjk2 mseIrrSikmoinsSkjSikOjkOikSjk3 &
                                         &mseIrrSikmoinsSkj_SikOjkOikSjk1 mseIrrSikmoinsSkj_SikOjkOikSjk2 mseIrrSikmoinsSkj_SikOjkOikSjk3 &
                                         &mseIrrSikOjkOikSjk1 mseIrrSikOjkOikSjk2 mseIrrSikOjkOikSjk3'
        endif
        write(file_id,'(i0,1x,i0,1x,51(g15.8,1x))') iFic,nd,mseIrrGM1,mseIrrGM2,mseIrrGM3,&
                                                          &mseIrrRGM1,mseIrrRGM2,mseIrrRGM3,&
                                                          &mseIrrRGM1_2var,mseIrrRGM2_2var,mseIrrRGM3_2var,&
                                                          &mseIrrRGM1_plus,mseIrrRGM2_plus,mseIrrRGM3_plus,&
                                                          &mseIrrSmag1,mseIrrSmag2,mseIrrSmag3,&
                                                          &mseIrrGMVort1,mseIrrGMVort2,mseIrrGMVort3,&
                                                          &mseIrrSikmoinsSkj1,mseIrrSikmoinsSkj2,mseIrrSikmoinsSkj3,&
                                                          &mseIrrSikmoinsSkjOikOjk1,mseIrrSikmoinsSkjOikOjk2,mseIrrSikmoinsSkjOikOjk3,&
                                                          &mseIrrSikmoinsSkj_OikOjk1,mseIrrSikmoinsSkj_OikOjk2,mseIrrSikmoinsSkj_OikOjk3,&
                                                          &mseIrrOikOjk1,mseIrrOikOjk2,mseIrrOikOjk3,&
                                                          &mseIrrgraddOkdxldUjdxl1,mseIrrgraddOkdxldUjdxl2,mseIrrgraddOkdxldUjdxl3,&
                                                          &mseIrrgraddOkdxlSjlmoins1,mseIrrgraddOkdxlSjlmoins2,mseIrrgraddOkdxlSjlmoins3,&
                                                          &mseIrrgraddOkdxlSjlplus1,mseIrrgraddOkdxlSjlplus2,mseIrrgraddOkdxlSjlplus3,&
                                                          &mseIrrgraddOkdxlOmegajl1,mseIrrgraddOkdxlOmegajl2,mseIrrgraddOkdxlOmegajl3,&
                                                 &mseIrrSikmoinsSkjSikOjkOikSjk1,mseIrrSikmoinsSkjSikOjkOikSjk2,mseIrrSikmoinsSkjSikOjkOikSjk3,&
                                                 &mseIrrSikmoinsSkj_SikOjkOikSjk1,mseIrrSikmoinsSkj_SikOjkOikSjk2,mseIrrSikmoinsSkj_SikOjkOikSjk3,&
                                                           &mseIrrSikOjkOikSjk1,mseIrrSikOjkOikSjk2,mseIrrSikOjkOikSjk3
        close(file_id)
        file_id=iclose(file_id)

        file_id=iopen()
        inquire(file='bilanMSEQUAD.table', exist=existe)
        open(unit=file_id, file='bilanMSEQUAD.table', form="FORMATTED",position="append")
        if ( .not. existe) then
          write(file_id, '(a)') '#iFic nd mseQuadGM1 mseQuadGM2 mseQuadGM3 &
                                         &mseQuadRGM1 mseQuadRGM2 mseQuadRGM3 &
                                         &mseQuadRGM1_plus mseQuadRGM2_plus mseQuadRGM3_plus mseIrrRGM1_plus &
                                         &mseQuadSmag1 mseQuadSmag2 mseQuadSmag3 &
                                         &mseQuadSmagStat1 mseQuadSmagStat2 mseQuadSmagStat3 &
                                         &mseQuadRGM1_2var mseQuadRGM2_2var mseQuadRGM3_2var &
                                         &mseQuadRGM1dyn mseQuadRGM2dyn mseQuadRGM3dyn &
                                         &mseQuadGMVort1 mseQuadGMVort2 mseQuadGMVort3 &
                                         &mseRGMSikmoinsSkjSikOjkOikSjkDyn1 mseRGMSikmoinsSkjSikOjkOikSjkDyn2 mseRGMSikmoinsSkjSikOjkOikSjkDyn3 &
                                         &mseQuadSikmoinsSkj1 mseQuadSikmoinsSkj2 mseQuadSikmoinsSkj3 &
                                         &mseRGMSikmoinsSkjSikOjkOikSjkDyn21 mseRGMSikmoinsSkjSikOjkOikSjkDyn22 mseRGMSikmoinsSkjSikOjkOikSjkDyn23 &
                                         &mseRGMSikmoinsSkjDyn1 mseRGMSikmoinsSkjDyn2 mseRGMSikmoinsSkjDyn3 &
                                         &mseRGMSikmoinsSkjDynFab1 mseRGMSikmoinsSkjDynFab2 mseRGMSikmoinsSkjDynFab3 &
                                         &mseSikmoinsSkjSikOjkOikSjkDynFab1 mseSikmoinsSkjSikOjkOikSjkDynFab2 mseSikmoinsSkjSikOjkOikSjkDynFab3 &
                                         &mseSikmoinsSkjSikOjkOikSjkDyn2Fab1 mseSikmoinsSkjSikOjkOikSjkDyn2Fab2 mseSikmoinsSkjSikOjkOikSjkDyn2Fab3'
        endif
        write(file_id,'(i0,1x,i0,1x,45(g15.8,1x))') iFic,nd,mseQuadGM1,mseQuadGM2,mseQuadGM3,&
                                                          &mseQuadRGM1,mseQuadRGM2,mseQuadRGM3,&
                                                          &mseQuadRGM1_plus,mseQuadRGM2_plus,mseQuadRGM3_plus,&
                                                          &mseQuadSmag1,mseQuadSmag2,mseQuadSmag3,&
                                                          &mseQuadSmagStat1,mseQuadSmagStat2,mseQuadSmagStat3,&
                                                          &mseQuadRGM1_2var,mseQuadRGM2_2var,mseQuadRGM3_2var,&
                                                          &mseQuadRGM1dyn,mseQuadRGM2dyn,mseQuadRGM3dyn,&
                                                          &mseQuadGMVort1,mseQuadGMVort2,mseQuadGMVort3,&
                                           &mseRGMSikmoinsSkjSikOjkOikSjkDyn1,mseRGMSikmoinsSkjSikOjkOikSjkDyn2,mseRGMSikmoinsSkjSikOjkOikSjkDyn3,&
                                                          &mseQuadSikmoinsSkj1,mseQuadSikmoinsSkj2,mseQuadSikmoinsSkj3,&
                                    &mseRGMSikmoinsSkjSikOjkOikSjkDyn21,mseRGMSikmoinsSkjSikOjkOikSjkDyn22,mseRGMSikmoinsSkjSikOjkOikSjkDyn23,&
                                              &mseRGMSikmoinsSkjDyn1,mseRGMSikmoinsSkjDyn2,mseRGMSikmoinsSkjDyn3,&
                                              &mseRGMSikmoinsSkjDynFab1,mseRGMSikmoinsSkjDynFab2,mseRGMSikmoinsSkjDynFab3,&
                                              &mseSikmoinsSkjSikOjkOikSjkDynFab1,mseSikmoinsSkjSikOjkOikSjkDynFab2,mseSikmoinsSkjSikOjkOikSjkDynFab3,&
                                              &mseSikmoinsSkjSikOjkOikSjkDyn2Fab1,mseSikmoinsSkjSikOjkOikSjkDyn2Fab2,mseSikmoinsSkjSikOjkOikSjkDyn2Fab3
        close(file_id)
        file_id=iclose(file_id)

        file_id=iopen()
        inquire(file='bilanDIFFDISS_GM_VS_RGM.table', exist=existe)
        open(unit=file_id, file='bilanDIFFDISS_GM_VS_RGM.table', form="FORMATTED",position="append")
        if ( .not. existe) then
          write(file_id,'(a)') '#iFic nd diffExactComp(1) diffExactComp(2) diffExactComp(3) diffExact &
                                          &dissExactComp(1) dissExactComp(2) dissExactComp(3) dissExact &
                                          &diffRGMComp(1) diffRGMComp(2) diffRGMComp(3) diffRGM &
                                          &dissRGMComp(1) dissRGMComp(2) dissRGMComp(3) dissRGM &
                                          &diffGMComp(1) diffGMComp(2) diffGMComp(3) diffGM &
                                          &dissGMComp(1) dissGMComp(2) dissGMComp(3) dissGM &
                                          &dissIrrgraddOkdxldUjdxl dissIrrgraddOkdxlOmegajl dissIrrgraddOkdxlSjlmoins &
                                          &dissIrrgraddOkdxlSjlplus dissIrrOikOjk dissIrrQDMvort &
                                          &dissIrrSikmoinsSkj dissIrrSikmoinsSkj_OikOjk dissIrrSikmoinsSkj_SikOjkOikSjk &
                                          &dissIrrSikmoinsSkjDyn dissIrrSikmoinsSkjOikOjk dissIrrSikmoinsSkjSikOjkOikSjk &
                                          &dissIrrSikOjkOikSjk dissQDMvort dissSikmoinsSkjStat &
                                          &dissSmagDyna dissSmagIrr dissSmagStat &
                                          &dissIrrSikmoinsSkjSikOjkOikSjkDyn dissIrrSikmoinsSkjSikOjkOikSjkDyn2 dissSikmoinsSkjSikOjkOikSjkDyn2Fab &
                                          &dissSikmoinsSkjSikOjkOikSjkDynFab dissSikmoinsSkjDynFab dissSikmoinsSkjDyn'
        endif
        write(file_id,'(i0,1x,i0,1x,48(g15.8,1x))') iFic,nd,diffExactComp(1),diffExactComp(2),diffExactComp(3),diffExact,&
                                          &dissExactComp(1),dissExactComp(2),dissExactComp(3),dissExact,&
                                          &diffRGMComp(1),diffRGMComp(2),diffRGMComp(3),diffRGM,&
                                          &dissRGMComp(1),dissRGMComp(2),dissRGMComp(3),dissRGM,&
                                          &diffGMComp(1),diffGMComp(2),diffGMComp(3),diffGM,&
                                          &dissGMComp(1),dissGMComp(2),dissGMComp(3),dissGM,&
                                          &dissIrrgraddOkdxldUjdxl,dissIrrgraddOkdxlOmegajl,dissIrrgraddOkdxlSjlmoins,&
                                          &dissIrrgraddOkdxlSjlplus,dissIrrOikOjk,dissIrrQDMvort,&
                                          &dissIrrSikmoinsSkj,dissIrrSikmoinsSkj_OikOjk,dissIrrSikmoinsSkj_SikOjkOikSjk,&
                                          &dissIrrSikmoinsSkjDyn,dissIrrSikmoinsSkjOikOjk,dissIrrSikmoinsSkjSikOjkOikSjk,&
                                          &dissIrrSikOjkOikSjk,dissQDMvort,dissSikmoinsSkjStat,&
                                          &dissSmagDyna,dissSmagIrr,dissSmagStat,&
                                          &dissIrrSikmoinsSkjSikOjkOikSjkDyn,dissIrrSikmoinsSkjSikOjkOikSjkDyn2,dissSikmoinsSkjSikOjkOikSjkDyn2Fab,&
                                          &dissSikmoinsSkjSikOjkOikSjkDynFab,dissSikmoinsSkjDynFab,dissSikmoinsSkjDyn
        close(file_id)
        file_id=iclose(file_id)
      endif
    enddo

    ! ===== Free memory and delete temporary fields =====
    if( .not. deleteWN(WN)) return
    call deleteDataLayout(Ufil)
    call deleteDataLayout(Vfil)
    call deleteDataLayout(Wfil)
    call deleteDataLayout(T1DNS)
    call deleteDataLayout(T2DNS)
    call deleteDataLayout(T3DNS)
    call deleteDataLayout(T1GM)
    call deleteDataLayout(T2GM)
    call deleteDataLayout(T3GM)
    call deleteDataLayout(T1RGM)
    call deleteDataLayout(T2RGM)
    call deleteDataLayout(T3RGM)
    call deleteDataLayout(T1RGMK)
    call deleteDataLayout(T2RGMK)
    call deleteDataLayout(T3RGMK)
    call deleteDataLayout(Ukfil)
    call deleteDataLayout(Vkfil)
    call deleteDataLayout(Wkfil)
    call deleteDataLayout(dudx)
    call deleteDataLayout(dudy)
    call deleteDataLayout(dudz)
    if (.not. deleteWorkArray(uvw)) return
    res = .true.

end subroutine aprioriRGMVel

!-------------------------------------------------------------------------------------------
!> Post-processing for the 2012 CTR summer program VELOCITY
!!
!! @author Guillaume Balarac
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in]    nbcpus          = number of process in the cpu pool
!!    @param[in]    U               = Component of velocity
!!    @param[in]    V               = Component of velocity
!!    @param[in]    W               = Component of velocity
!!    @param[in]    Scal            = Scalar to compute
!!    @param[in]    ite             = the number of file
!!    @return       res             = logical value
!-------------------------------------------------------------------------------------------
function gradientvelmodOE(U,V,W,Scal,nbcpus,spec_rank) result (success)

    use conditional_mean_tools
    use subgrid_tools , only : computeSplitingSij

#ifdef LAPACK
  external DGEEV
#endif
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !I/O data
  type(real_data_layout),intent(inout)     :: U,V,W,Scal
  integer,intent(in)                       :: spec_rank,nbcpus
  logical                                  :: success

  !Local data
  type(real_data_layout)                             :: S11, S12, S13, S22, S23, S33
  type(real_data_layout)                             :: Sm11, Sm12, Sm13, Sm22, Sm23, Sm33
  type(real_data_layout)                             :: T11, T12, T13, T22, T23, T33
  type(real_data_layout)                             :: Z11, Z12, Z13, Z22, Z23, Z33
  type(real_data_layout),dimension(:),allocatable    :: phi, phidiss, TiDNS, TiM, SijmTj
  type(real_data_layout),dimension(:),allocatable    :: temp, temp1, temp2, temp3, temp4 , temp5
  type(real_data_layout),dimension(:),allocatable    :: VelScal_ft, Vel_ft, dZdxi_ft, SdZdxi_ft, Li, Mi
  type(real_data_layout)                             :: S11_ft,S12_ft,S13_ft, S22_ft, S23_ft, S33_ft
  type(real_data_layout)                             :: Uf, Vf, Wf

  type(complex_data_layout),dimension(:),allocatable :: tempK, tempKf
  type(complex_data_layout)                          :: ExK

  type(wavenumbers)                                  :: ScalWaveN
  real(WP)                                           :: qe_gm, ie_gm, deltx, DynCoef, filtersize
  logical                                            :: res
  integer                                            :: bin1,bin2,bin3,typeDisc,n_iter
  integer                                            :: filter
  integer                                            :: nd,ii
  integer                                            :: i,j,k
  character(len=64)                                  :: numFil


  !! for eigenvalues computation
  real(WP), dimension(3,3) :: E1,E2,E3
#ifdef LAPACK
  integer                  :: a, b, info
  real(WP), dimension(3,3) :: tab,VecP,VL,TvecP
  real(WP), dimension(3,3) :: M1,M2,M3
  real(WP), dimension(3)   :: valR,valI
  real(WP)                 :: tabwork(102)
#endif
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  success = .false.

  bin1 = 100
  bin2 = bin1
  bin3 = bin1
  typeDisc  = 1
  filter = 1

  n_iter = 100

  deltx = computeDelta(U)

  !allocation

  if (.not. initDataLayout("ExK", ExK,(U%nx/2)+1,U%ny,U%nz, &
     & U%Lx,U%Ly,U%Lz,nbcpus,spec_rank,alongZ)) then
     write(6,'(a,i0,a)')'[ERROR] initDataLayout for TYPE (COMPLEX_DATA_LAYOUT) : not enought memory!'
     return
  endif
  res = initWorkArray(ExK,4,tempK)
  res = initWorkArray(ExK,4,tempKf)
#ifndef DEBUGFONCTION
  call deletedatalayout(ExK)
#endif
  if (.not. copyStructOnly(U,Sm11,"Sm11") .or. & 
     &.not. copyStructOnly(U,Sm12,"Sm12") .or. &
     &.not. copyStructOnly(U,Sm13,"Sm13") .or. &
     &.not. copyStructOnly(U,Sm22,"Sm22") .or. &
     &.not. copyStructOnly(U,Sm23,"Sm23") .or. &
     &.not. copyStructOnly(U,Sm33,"Sm33") .or. &
     &.not. copyStructOnly(U,S11,"S11") .or. &
     &.not. copyStructOnly(U,S12,"S12") .or. &
     &.not. copyStructOnly(U,S13,"S13") .or. &
     &.not. copyStructOnly(U,S22,"S22") .or. &
     &.not. copyStructOnly(U,S23,"S23") .or. &
     &.not. copyStructOnly(U,S33,"S33") .or. &
     &.not. copyStructOnly(U,T11,"T11") .or. &
     &.not. copyStructOnly(U,T12,"T12") .or. &
     &.not. copyStructOnly(U,T13,"T13") .or. &
     &.not. copyStructOnly(U,T22,"T22") .or. &
     &.not. copyStructOnly(U,T23,"T23") .or. &
     &.not. copyStructOnly(U,T33,"T33") .or. &
     &.not. copyStructOnly(U,Z11,"Z11") .or. &
     &.not. copyStructOnly(U,Z12,"Z12") .or. &
     &.not. copyStructOnly(U,Z13,"Z13") .or. &
     &.not. copyStructOnly(U,Z22,"Z22") .or. &
     &.not. copyStructOnly(U,Z23,"Z23") .or. &
     &.not. copyStructOnly(U,Z33,"Z33") .or. &
     &.not. copyStructOnly(U,S11_ft,"S11_ft") .or. &
     &.not. copyStructOnly(U,S12_ft,"S12_ft") .or. &
     &.not. copyStructOnly(U,S13_ft,"S13_ft") .or. &
     &.not. copyStructOnly(U,S22_ft,"S22_ft") .or. &
     &.not. copyStructOnly(U,S23_ft,"S23_ft") .or. &
     &.not. copyStructOnly(U,S33_ft,"S33_ft") .or. &
     &.not. copyStructOnly(U,Uf,"Uf") .or. &
     &.not. copyStructOnly(U,Vf,"Vf") .or. &
     &.not. copyStructOnly(U,Wf,"Wf") ) then
     write(6,'(a,i0,a)')'[ERROR] copyStructOnly for TYPE (REAL_DATA_LAYOUT) : not enought memory!'
     return
  endif

  if (.not. initWorkArray(U,3,temp) .or. &
     &.not. initWorkArray(U,3,temp1) .or. &
     &.not. initWorkArray(U,3,temp2) .or. &
     &.not. initWorkArray(U,3,temp3) .or. &
     &.not. initWorkArray(U,3,temp4) .or. &
     &.not. initWorkArray(U,3,temp5) .or. &
     &.not. initWorkArray(U,4,phi) .or. &
     &.not. initWorkArray(U,3,TiDNS) .or. &
     &.not. initWorkArray(U,3,TiM) .or. &
     &.not. initWorkArray(U,3,SijmTj) .or. &
     &.not. initWorkArray(U,4,phidiss) .or. &
     &.not. initWorkArray(U,3,VelScal_ft) .or. &
     &.not. initWorkArray(U,3,Vel_ft) .or. &
     &.not. initWorkArray(U,3,dZdxi_ft) .or. &
     &.not. initWorkArray(U,3,SdZdxi_ft) .or. &
     &.not. initWorkArray(U,3,Li) .or. &
     &.not. initWorkArray(U,3,Mi)) then
     write(6,'(a,i0,a)')'[ERROR] initWorkArray for TYPE (REAL_DATA_LAYOUT) : not enought memory!'
     return
  endif


  res = initWN(ScalWaveN,spec_rank,Scal%nx,Scal%ny,Scal%nz)
  res = computeWN(ScalWaveN,spec_rank,Scal%Lx,Scal%Ly,Scal%Lz,Scal%nx,Scal%ny,Scal%nz)

  if (spec_rank.EQ.0) then
     open(20,file='IE_QE_T1_DS.out',form='formatted')
     open(21,file='IE_QE_T2_DS.out',form='formatted')
     open(22,file='IE_QE_T3_DS.out',form='formatted')
     open(23,file='mean_diss_DS.out',form='formatted')

     open(30,file='IE_QE_T1_S-.out',form='formatted')
     open(31,file='IE_QE_T2_S-.out',form='formatted')
     open(32,file='IE_QE_T3_S-.out',form='formatted')
     open(33,file='mean_diss_S-.out',form='formatted')
  endif

!  if (.not. diff_tool_initBuff(tempK(1))) return

  call ftran(U,tempK(1),res)
  call ftran(V,tempK(2),res)
  call ftran(W,tempK(3),res)
  ie_gm = computeDiv( spec_rank,tempK(1),tempK(2),tempK(3),tempK(4),ScalWaveN )
  IF (spec_rank.EQ.0) WRITE(6,'(a,g15.8)')'  [INFO] Divergence max of U field is ',ie_gm
  IF (.NOT. NullifyDiv(tempK(1),tempK(2),tempK(3), ScalWaveN, spec_rank)) RETURN
  ie_gm = computeDiv( spec_rank,tempK(1),tempK(2),tempK(3),tempK(4),ScalWaveN )
  IF (spec_rank.EQ.0) WRITE(6,'(a,g15.8)')'  [INFO] Divergence max of U field is ',ie_gm

  !! -> loop on filter size
  do ii = 1,8

    nd = 2*ii
    filtersize = real(nd)*deltx
    if (spec_rank .eq. 0)  print*,'[INFO] ',ii,' TAILLE de FILTRE = ',filtersize

    ! Compute Tij DNS
    !!!!!!!!!!!!!!!!!
    ! -> Compute the filtered velocities and the filtered scalar
    !!!!!!!!!!!!!!!!!
    call ftran(U,tempK(1),res)
    call ftran(V,tempK(2),res)
    call ftran(W,tempK(3),res)
    call computeFilter(ScalWaveN,filtersize,tempK(1),tempKf(1),filter) !U bar
    call computeFilter(ScalWaveN,filtersize,tempK(2),tempKf(2),filter) !V bar
    call computeFilter(ScalWaveN,filtersize,tempK(3),tempKf(3),filter) !W bar
    call btran(tempKf(1),Uf,res) !U bar
    call btran(tempKf(2),Vf,res) !V bar
    call btran(tempKf(3),Wf,res) !W bar

    !!!!!!!!!!!!!!!!!
    ! -> Compute TiDNS
    !!!!!!!!!!!!!!!!!
    T11%values = U%values * U%values
    T12%values = V%values * U%values
    T13%values = W%values * U%values
    T22%values = V%values * V%values
    T23%values = V%values * W%values
    T33%values = W%values * W%values
    call ftran(T11,tempK(1),res)
    call ftran(T12,tempK(2),res)
    call ftran(T13,tempK(3),res)
    call computeFilter(ScalWaveN,filtersize,tempK(1),tempKf(1),filter) !U bar
    call computeFilter(ScalWaveN,filtersize,tempK(2),tempKf(2),filter) !V bar
    call computeFilter(ScalWaveN,filtersize,tempK(3),tempKf(3),filter) !W bar
    call btran(tempKf(1),T11,res) !U bar
    call btran(tempKf(2),T12,res) !V bar
    call btran(tempKf(3),T13,res) !W bar
    call ftran(T22,tempK(1),res)
    call ftran(T23,tempK(2),res)
    call ftran(T33,tempK(3),res)
    call computeFilter(ScalWaveN,filtersize,tempK(1),tempKf(1),filter) !U bar
    call computeFilter(ScalWaveN,filtersize,tempK(2),tempKf(2),filter) !V bar
    call computeFilter(ScalWaveN,filtersize,tempK(3),tempKf(3),filter) !W bar
    call btran(tempKf(1),T22,res) !U bar
    call btran(tempKf(2),T23,res) !V bar
    call btran(tempKf(3),T33,res) !W bar
    T11%values = T11%values - Uf%values * Uf%values
    T12%values = T12%values - Uf%values * Vf%values
    T13%values = T13%values - Uf%values * Wf%values
    T22%values = T22%values - Vf%values * Vf%values
    T23%values = T23%values - Vf%values * Wf%values
    T33%values = T33%values - Wf%values * Wf%values

    ! compute duidxj
    ! temp1 <- dUfdxi
    ! temp2 <- dVfdxi
    ! temp3 <- dWfdxi
    res = computeFieldGradient(ScalWaveN,Uf,temp1(1),temp1(2),temp1(3),nbcpus,spec_rank)
    res = computeFieldGradient(ScalWaveN,Vf,temp2(1),temp2(2),temp2(3),nbcpus,spec_rank)
    res = computeFieldGradient(ScalWaveN,Wf,temp3(1),temp3(2),temp3(3),nbcpus,spec_rank)

    ! compute Sij = 1/2 * ( duidxj + dujdxi) and Oij = 1/2 * ( duidxj - dujdxi )
    S11%values = temp1(1)%values
    S12%values = 0.5 * ( temp1(2)%values + temp2(1)%values )
    S13%values = 0.5 * ( temp1(3)%values + temp3(1)%values )
    S22%values = temp2(2)%values
    S23%values = 0.5 * ( temp2(3)%values + temp3(2)%values )
    S33%values = temp3(3)%values

    !! Smagorinsky model

    ! compute |Sij| <-- temp5(1)
    call computeTensorNorme1(S11,S12,S13,S22,S23,S33,temp5(1),res)

    !Dynamic Procedure
    DynCoef = - 2.0_WP * (0.18_WP)**2
    ! TO DO HERE

    ! QE and IE for model Tij = C . deltx^2 . |Sij| . Sij
    ! for dT1jdxj
    Z11%values = DynCoef * real(filtersize)**2.0*temp5(1)%values*S11%values
    Z12%values = DynCoef * real(filtersize)**2.0*temp5(1)%values*S12%values
    Z13%values = DynCoef * real(filtersize)**2.0*temp5(1)%values*S13%values
    Z22%values = DynCoef * real(filtersize)**2.0*temp5(1)%values*S22%values
    Z23%values = DynCoef * real(filtersize)**2.0*temp5(1)%values*S23%values
    Z33%values = DynCoef * real(filtersize)**2.0*temp5(1)%values*S33%values

    ! take divergence -- this is the parameter  <-- phi(1)
    res = computeFieldDivergence(ScalWaveN,Z11,Z12,Z13,phi(1),nbcpus,spec_rank)
    ! take divergence of exact term -- this is the objective <-- phi(2)
    res = computeFieldDivergence(ScalWaveN,T11,T12,T13,phi(2),nbcpus,spec_rank)
    ! Compute the optimal estimator of GM
    res = computeEONVar( spec_rank,(/phi(1)/),& !Les variables
                       & phi(2),&                            !Le champ à estimer
                       & (/typeDisc/),&     !Les discrétisations par variable
                       & (/bin1/),&                 !le nombre d'intervalles par variable
                       & fieldOpt=temp(1))                    !champ optimal en sortie
    ! IE of GM
    ie_gm = computeMSE(phi(2),temp(1),spec_rank,.true.)     ! true = adim par la variance
    ! QE of GM
    qe_gm = computeMSE(phi(2),phi(1),spec_rank,.true.)     ! true = adim par la variance
    if (spec_rank.EQ.0) then
       write(20,*) nd, ie_gm, qe_gm
    end if

    ! take divergence -- this is the parameter  <-- phi(1)
    res = computeFieldDivergence(ScalWaveN,Z12,Z22,Z23,phi(1),nbcpus,spec_rank)
    ! take divergence of exact term -- this is the objective <-- phi(2)
    res = computeFieldDivergence(ScalWaveN,T12,T22,T23,phi(2),nbcpus,spec_rank)
    ! Compute the optimal estimator of GM
    res = computeEONVar( spec_rank,(/phi(1)/),& !Les variables
                       & phi(2),&                            !Le champ à estimer
                       & (/typeDisc/),&     !Les discrétisations par variable
                       & (/bin1/),&                 !le nombre d'intervalles par variable
                       & fieldOpt=temp(1))                    !champ optimal en sortie
    ! IE of GM
    ie_gm = computeMSE(phi(2),temp(1),spec_rank,.true.)     ! true = adim par la variance
    ! QE of GM
    qe_gm = computeMSE(phi(2),phi(1),spec_rank,.true.)     ! true = adim par la variance
    if (spec_rank.EQ.0) then
       write(21,*) nd, ie_gm, qe_gm
    end if

    ! take divergence -- this is the parameter  <-- phi(1)
    res = computeFieldDivergence(ScalWaveN,Z13,Z23,Z33,phi(1),nbcpus,spec_rank)
    ! take divergence of exact term -- this is the objective <-- phi(2)
    res = computeFieldDivergence(ScalWaveN,T13,T23,T33,phi(2),nbcpus,spec_rank)
    ! Compute the optimal estimator of GM
    res = computeEONVar( spec_rank,(/phi(1)/),& !Les variables
                       & phi(2),&                            !Le champ à estimer
                       & (/typeDisc/),&     !Les discrétisations par variable
                       & (/bin1/),&                 !le nombre d'intervalles par variable
                       & fieldOpt=temp(1))                    !champ optimal en sortie
    ! IE of GM
    ie_gm = computeMSE(phi(2),temp(1),spec_rank,.true.)     ! true = adim par la variance
    ! QE of GM
    qe_gm = computeMSE(phi(2),phi(1),spec_rank,.true.)     ! true = adim par la variance
    if (spec_rank.EQ.0) then
       write(22,*) nd, ie_gm, qe_gm
    end if

    ! SGS dissipation
    temp5(1)%values = T11%values * S11%values +  T22%values * S22%values +  T33%values * S33%values + 2.0_WP*(T12%values * S12%values +  T13%values * S13%values +  T23%values * S23%values)
    temp5(2)%values = Z11%values * S11%values +  Z22%values * S22%values +  Z33%values * S33%values + 2.0_WP*(Z12%values * S12%values +  Z13%values * S13%values +  Z23%values * S23%values)
    ie_gm = computeFieldAvg( temp5(1) , spec_rank , nbcpus )
    qe_gm = computeFieldAvg( temp5(2) , spec_rank , nbcpus )
    if (spec_rank.EQ.0) then
       write(23,*) nd, ie_gm, qe_gm
    end if


    ! Proposed model : Tij = C . deltx^2 . Sik^- . Sjk
    ! Compute Sik^- <-- Smik
    E1 = 0
    E2 = 0
    E3 = 0
    E1(1,1) = 1._WP
    E2(2,2) = 1._WP
    E3(3,3) = 1._WP
#ifdef LAPACK
    do k = U%zmin, U%zmax
       do j = U%ymin, U%ymax
          do i = U%xmin, U%xmax
             tab(1,1) = S11%values(i,j,k)
             tab(1,2) = S12%values(i,j,k)
             tab(1,3) = S13%values(i,j,k)
             tab(2,1) = S12%values(i,j,k)
             tab(2,2) = S22%values(i,j,k)
             tab(2,3) = S23%values(i,j,k)
             tab(3,1) = S13%values(i,j,k)
             tab(3,2) = S23%values(i,j,k)
             tab(3,3) = S33%values(i,j,k)
             call DGEEV('N','V',3,tab,3,valR,valI,VL,3,VecP,3,tabwork,102,INFO)
             do a=1,3
               do b=1,3
                 TvecP(a,b)=VecP(b,a)
               enddo
             enddo
             M1 = matmul(E1,TvecP)
             M1 = matmul(vecP,M1)
             M2 = matmul(E2,TvecP)
             M2 = matmul(vecP,M2)
             M3 = matmul(E3,TvecP)
             M3 = matmul(vecP,M3)

             if (M1(1,2).ne.M1(2,1)) print*,'M1_12 ne M1_21', M1(1,2), M1(2,1)
             if (M1(1,3).ne.M1(3,1)) print*,'M1_13 ne M1_31', M1(1,3), M1(3,1)
             if (M1(2,3).ne.M1(3,2)) print*,'M1_23 ne M1_32', M1(2,3), M1(3,2)
             if (M2(1,2).ne.M2(2,1)) print*,'M2_12 ne M2_21', M2(1,2), M2(2,1)
             if (M2(1,3).ne.M2(3,1)) print*,'M2_13 ne M2_31', M2(1,3), M2(3,1)
             if (M2(2,3).ne.M2(3,2)) print*,'M2_23 ne M2_32', M2(2,3), M2(3,2)
             if (M3(1,2).ne.M3(2,1)) print*,'M3_12 ne M3_21', M3(1,2), M3(2,1)
             if (M3(1,3).ne.M3(3,1)) print*,'M3_13 ne M3_31', M3(1,3), M3(3,1)
             if (M3(2,3).ne.M3(3,2)) print*,'M3_23 ne M3_32', M3(2,3), M3(3,2)

             Sm11%values(i,j,k)     =  MIN(valR(1),0.0_WP)*M1(1,1) + MIN(valR(2),0.0_WP)*M2(1,1) + MIN(valR(3),0.0_WP)*M3(1,1)
             Sm12%values(i,j,k)     =  MIN(valR(1),0.0_WP)*M1(1,2) + MIN(valR(2),0.0_WP)*M2(1,2) + MIN(valR(3),0.0_WP)*M3(1,2)
             Sm13%values(i,j,k)     =  MIN(valR(1),0.0_WP)*M1(1,3) + MIN(valR(2),0.0_WP)*M2(1,3) + MIN(valR(3),0.0_WP)*M3(1,3)
             Sm22%values(i,j,k)     =  MIN(valR(1),0.0_WP)*M1(2,2) + MIN(valR(2),0.0_WP)*M2(2,2) + MIN(valR(3),0.0_WP)*M3(2,2)
             Sm23%values(i,j,k)     =  MIN(valR(1),0.0_WP)*M1(2,3) + MIN(valR(2),0.0_WP)*M2(2,3) + MIN(valR(3),0.0_WP)*M3(2,3)
             Sm33%values(i,j,k)     =  MIN(valR(1),0.0_WP)*M1(3,3) + MIN(valR(2),0.0_WP)*M2(3,3) + MIN(valR(3),0.0_WP)*M3(3,3)

          end do
       end do
    end do
#endif
#ifndef LAPACK
  !print *,'[WARNING] LAPACK not compiled [WARNING]'
  !success=.false.
  !return
  if (.not. computeSplitingSij(S11,S12,S13,S22,S23,S33,&
                          & Sm11,Sm12,Sm13,Sm22,Sm23,Sm33,&
                          & 'moins')) return
#endif

    ! Dynamic procedure
    DynCoef = 1.0_WP/12.0_WP
    ! TO DO HERE

    ! Compute Zij = C . deltx^2 . Sik^- . Sjk   (i.e. Smik * Sjk)  !! FINISH %VALUES
    Z11%values = DynCoef * (filtersize)**2.0 * (Sm11%values*S11%values + Sm12%values*S12%values + Sm13%values*S13%values)
    Z12%values = DynCoef * (filtersize)**2.0 * (Sm11%values*S12%values + Sm12%values*S22%values + Sm13%values*S23%values)
    Z13%values = DynCoef * (filtersize)**2.0 * (Sm11%values*S13%values + Sm12%values*S23%values + Sm13%values*S33%values)
    Z22%values = DynCoef * (filtersize)**2.0 * (Sm12%values*S12%values + Sm22%values*S22%values + Sm23%values*S23%values)
    Z23%values = DynCoef * (filtersize)**2.0 * (Sm12%values*S13%values + Sm22%values*S23%values + Sm23%values*S33%values)
    Z33%values = DynCoef * (filtersize)**2.0 * (Sm13%values*S13%values + Sm23%values*S23%values + Sm33%values*S33%values)

    ! take divergence -- this is the parameter  <-- phi(1)
    res = computeFieldDivergence(ScalWaveN,Z11,Z12,Z13,phi(1),nbcpus,spec_rank)
    ! take divergence of exact term -- this is the objective <-- phi(2)
    res = computeFieldDivergence(ScalWaveN,T11,T12,T13,phi(2),nbcpus,spec_rank)
    ! Compute the optimal estimator of GM
    res = computeEONVar( spec_rank,(/phi(1)/),& !Les variables
                       & phi(2),&                            !Le champ à estimer
                       & (/typeDisc/),&     !Les discrétisations par variable
                       & (/bin1/),&                 !le nombre d'intervalles par variable
                       & fieldOpt=temp(1))                    !champ optimal en sortie
    ! IE of GM
    ie_gm = computeMSE(phi(2),temp(1),spec_rank,.true.)     ! true = adim par la variance
    ! QE of GM
    qe_gm = computeMSE(phi(2),phi(1),spec_rank,.true.)     ! true = adim par la variance
    if (spec_rank.EQ.0) then
       write(30,*) nd, ie_gm, qe_gm
    end if

    ! take divergence -- this is the parameter  <-- phi(1)
    res = computeFieldDivergence(ScalWaveN,Z12,Z22,Z23,phi(1),nbcpus,spec_rank)
    ! take divergence of exact term -- this is the objective <-- phi(2)
    res = computeFieldDivergence(ScalWaveN,T12,T22,T23,phi(2),nbcpus,spec_rank)
    ! Compute the optimal estimator of GM
    res = computeEONVar( spec_rank,(/phi(1)/),& !Les variables
                       & phi(2),&                            !Le champ à estimer
                       & (/typeDisc/),&     !Les discrétisations par variable
                       & (/bin1/),&                 !le nombre d'intervalles par variable
                       & fieldOpt=temp(1))                    !champ optimal en sortie
    ! IE of GM
    ie_gm = computeMSE(phi(2),temp(1),spec_rank,.true.)     ! true = adim par la variance
    ! QE of GM
    qe_gm = computeMSE(phi(2),phi(1),spec_rank,.true.)     ! true = adim par la variance
    if (spec_rank.EQ.0) then
       write(31,*) nd, ie_gm, qe_gm
    end if

    ! take divergence -- this is the parameter  <-- phi(1)
    res = computeFieldDivergence(ScalWaveN,Z13,Z23,Z33,phi(1),nbcpus,spec_rank)
    ! take divergence of exact term -- this is the objective <-- phi(2)
    res = computeFieldDivergence(ScalWaveN,T13,T23,T33,phi(2),nbcpus,spec_rank)
    ! Compute the optimal estimator of GM
    res = computeEONVar( spec_rank,(/phi(1)/),& !Les variables
                       & phi(2),&                            !Le champ à estimer
                       & (/typeDisc/),&     !Les discrétisations par variable
                       & (/bin1/),&                 !le nombre d'intervalles par variable
                       & fieldOpt=temp(1))                    !champ optimal en sortie
    ! IE of GM
    ie_gm = computeMSE(phi(2),temp(1),spec_rank,.true.)     ! true = adim par la variance
    ! QE of GM
    qe_gm = computeMSE(phi(2),phi(1),spec_rank,.true.)     ! true = adim par la variance
    if (spec_rank.EQ.0) then
       write(32,*) nd, ie_gm, qe_gm
    end if

    ! SGS dissipation
    temp5(1)%values = T11%values * S11%values +  T22%values * S22%values +  T33%values * S33%values + 2.0_WP*(T12%values * S12%values +  T13%values * S13%values +  T23%values * S23%values)
    temp5(2)%values = Z11%values * S11%values +  Z22%values * S22%values +  Z33%values * S33%values + 2.0_WP*(Z12%values * S12%values +  Z13%values * S13%values +  Z23%values * S23%values)
    write(numFil,'(i2.2)') nd
    temp5(1)%name='dissDivTijExact_'//trim(numFil)
    if (.not. computeFieldPDF(1,temp5(1),300,spec_rank,10000)) return
    ie_gm = computeFieldAvg( temp5(1) , spec_rank , nbcpus )
    qe_gm = computeFieldAvg( temp5(2) , spec_rank , nbcpus )
    if (spec_rank.EQ.0) then
       write(33,*) nd, ie_gm, qe_gm
    end if
    do k = U%zmin, U%zmax
       do j = U%ymin, U%ymax
          do i = U%xmin, U%xmax
             if (temp(2)%values(i,j,k).gt.0.0_WP) print*,'Something wrong for velocity model :', temp(2)%values(i,j,k)
          end do
       end do
    end do

  end do

  success = .true.

end function gradientvelmodOE


function velo_test_interpol(U, V, W, scal, nbcpus, spec_rank) result(success)

  use post_lib, only : compute_velo_spectrum
  use Interpolation_velo, only : Interpol_3D, Interpol_init
  use data
  use datalayout
  USE param
  use stat_tools, only : computeFieldMin, computeFieldMax

  type(REAL_DATA_LAYOUT), intent(inout) :: U, V, W, scal
  integer, intent(in)               :: nbcpus, spec_rank
  ! Local variable
  real(WP), dimension(3)            :: dx_f, dx_c
  real(WP)                          :: mini, maxi
  logical                           :: success
  CHARACTER(len=4)        :: advec_interp! choose advection solver (p_ON = particle method at order N, N = 2 or 4)

  ! ===== Initialisation =====
  success = .false.
  if(spec_rank == 0) write(*,'(a)') '  [STATUS] Start interpolation'
  ! Init interpolation context
  call parser_read('advection interpol', advec_interp)
  call interpol_init(advec_interp, .true.)
  ! Compute space step
  dx_c = (/ U%Lx/U%Nx, U%Ly/U%Ny, U%Lz/U%Nz /)
  dx_f = (/ scal%Lx/scal%Nx, scal%Ly/scal%Ny, scal%Lz/scal%Nz /)
  if(spec_rank == 0) write(*,'(a)') '  [STATUS] Init ok'

  ! ====== Compute spectrum before interpolation =====
  IF (.NOT. initDataLayout("U_Spectral", &
          & Uk,(U%nx/2)+1,U%ny,U%nz,U%Lx,U%Ly,U%Lz,nbcpus,spec_rank,alongZ)) THEN
      WRITE(6,'(a,i0,a)')'[ERROR] velo_test_interpol on process ',spec_rank,&
          & ': not enought memory for U velocities in fourier space.'
      RETURN
  ENDIF
  IF (.NOT. initDataLayout("V_Spectral", &
          & Vk,(V%nx/2)+1,V%ny,V%nz,V%Lx,V%Ly,V%Lz,nbcpus,spec_rank,alongZ)) THEN
      WRITE(6,'(a,i0,a)')'[ERROR] velo_test_interpol on process ',spec_rank,&
          & ': not enought memory for V velocities in fourier space.'
      RETURN
  ENDIF
  IF (.NOT. initDataLayout("W_Spectral", &
          & Wk,(W%nx/2)+1,W%ny,W%nz,W%Lx,W%Ly,W%Lz,nbcpus,spec_rank,alongZ)) THEN
      WRITE(6,'(a,i0,a)')'[ERROR] velo_test_interpol on process ',spec_rank,&
          & ': not enought memory for W velocities in fourier space.'
      RETURN
  ENDIF
  IF(.not. wavesNumber_init(VelWN,spec_rank,U%Lx,U%Ly,U%Lz,U%nx,U%ny,U%nz)) THEN
  write(*,'(a)') '[ERROR] velo_test_interpol : fail to compute WN'
  return
  end if
  call ftran(U, Uk,success)
  call ftran(V, Vk,success)
  call ftran(W, Wk,success)
  if(.not.compute_velo_spectrum(0, spec_rank,1)) then
    write(*,'(a)') '[ERROR] velo_test_interpol : fail to compute velo spectrum'
    return
  end if
  if(spec_rank == 0) write(*,'(a)') '  [STATUS] Spectrum before interpol OK'
  call deleteDatalayout(Uk)
  call deleteDatalayout(Vk)
  call deleteDatalayout(Wk)
  if(.not. deleteWN(VelWN)) write(*,'(a)') '[ERROR] velo_test_interpol : fail to delete WN'

  ! ===== Perform interpolation without permutation =====
  if(spec_rank==0) write(*,'(a, 3g15.8)') '      [INFO] velo_test_interpol : dx_f/dx_c = ', dx_f/dx_c
  call showVeloMinMaxAvg(U,V,W,spec_rank,nbcpus,0._WP)
  
  call Interpol_3D(U%values, dx_c, scal%values, dx_f)
  maxi = computeFieldMax(scal, spec_rank)
  mini = computeFieldMin(scal, spec_rank)
  if(spec_rank == 0) write(*,'(2(a,g15.8))') '    [INFO] Interpol U: max = ', maxi, ' min = ', mini
  call deleteDatalayout(U)
  scal%name = "U"
  U = scal
  call Interpol_3D(V%values, dx_c, scal%values, dx_f)
  maxi = computeFieldMax(scal, spec_rank)
  mini = computeFieldMin(scal, spec_rank)
  if(spec_rank == 0) write(*,'(2(a,g15.8))') '    [INFO] Interpol V: max = ', maxi, ' min = ', mini
  call deleteDatalayout(V)
  V = scal
  V%name = "V"
  call Interpol_3D(W%values, dx_c, scal%values, dx_f)
  maxi = computeFieldMax(scal, spec_rank)
  mini = computeFieldMin(scal, spec_rank)
  if(spec_rank == 0) write(*,'(2(a,g15.8))') '    [INFO] Interpol W: max = ', maxi, ' min = ', mini
  call deleteDatalayout(W)
  W = scal
  W%name = "W"
  call deleteDatalayout(scal)
  if(spec_rank == 0) write(*,'(a)') '  [STATUS] Interpol ok'
  call showVeloMinMaxAvg(U,V,W,spec_rank,nbcpus,0._WP)


  !  ===== Allocate storage for velocities in fourier space =====
  IF (.NOT. initDataLayout("U_Spectral", &
          & Uk,(U%nx/2)+1,U%ny,U%nz,U%Lx,U%Ly,U%Lz,nbcpus,spec_rank,alongZ)) THEN
      WRITE(6,'(a,i0,a)')'[ERROR] velo_test_interpol on process ',spec_rank,&
          & ': not enought memory for U velocities in fourier space.'
      RETURN
  ENDIF
  IF (.NOT. initDataLayout("V_Spectral", &
          & Vk,(V%nx/2)+1,V%ny,V%nz,V%Lx,V%Ly,V%Lz,nbcpus,spec_rank,alongZ)) THEN
      WRITE(6,'(a,i0,a)')'[ERROR] velo_test_interpol on process ',spec_rank,&
          & ': not enought memory for V velocities in fourier space.'
      RETURN
  ENDIF
  IF (.NOT. initDataLayout("W_Spectral", &
          & Wk,(W%nx/2)+1,W%ny,W%nz,W%Lx,W%Ly,W%Lz,nbcpus,spec_rank,alongZ)) THEN
      WRITE(6,'(a,i0,a)')'[ERROR] velo_test_interpol on process ',spec_rank,&
          & ': not enought memory for W velocities in fourier space.'
      RETURN
  ENDIF
  IF(.not. wavesNumber_init(VelWN,spec_rank,U%Lx,U%Ly,U%Lz,U%nx,U%ny,U%nz)) THEN
  write(*,'(a)') '[ERROR] velo_test_interpol : fail to compute WN'
  return
  end if
  if(spec_rank == 0) write(*,'(a)') '  [STATUS] Init spec ok'





  ! ===== Go in spectral space =====
  call ftran(U, Uk,success)
  call ftran(V, Vk,success)
  call ftran(W, Wk,success)
  if (.not. success) then
    write(*,'(a)') '[ERROR] velo_test_interpol : fail to compute FFT'
    return
  end if
  success = .false.

  ! ===== Compute stats =====
  ! Kinetic spectrum
  if(.not.compute_velo_spectrum(1, spec_rank,1)) then
    write(*,'(a)') '[ERROR] velo_test_interpol : fail to compute velo spectrum'
    return
  end if
  if(spec_rank == 0) write(*,'(a)') '  [STATUS] Spectrum ok'
  call deleteDatalayout(Vk)
  call deleteDatalayout(Wk)
  ! Pdf
  if(.not.computeFieldPDF(1,U,400,spec_rank,1)) then
    write(*,'(a)') '[ERROR] velo_test_interpol : fail to compute U pdf'
    return
  end if
  if(spec_rank == 0) write(*,'(a)') '  [STATUS] PDF ok'
  ! Gradient pdf
  if(.not.computeDerivationField_Spe(Uk,(/1/),U,spec_rank,VelWN)) then
    write(*,'(a)') '[ERROR] velo_test_interpol : fail to compute dU/dx'
    return
  end if
  U%name = "U_dx"
  if(.not.computeFieldPDF(1,U,400,spec_rank,1)) then
    write(*,'(a)') '[ERROR] velo_test_interpol : fail to compute pdf of dU/dx'
    return
  end if
  if(spec_rank == 0) write(*,'(a)') '  [STATUS] Grad PDF ok'

  success = .true.

end function velo_test_interpol

!-------------------------------------------------------------------------------------------
!> @detail
!> check the equality between gradient model and gradient model with incompressibility 
!!
!! @author Antoine Vollant, LEGI
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in]    nbcpus          = number of process in the cpu pool
!!    @param[in]    U               = Component of velocity
!!    @param[in]    V               = Component of velocity
!!    @param[in]    W               = Component of velocity
!!    @return       res             = logical value
!
!-------------------------------------------------------------------------------------------
subroutine aprioriSikmoinsSkj(iFic,spec_rank,nbcpus,res,U,V,W)

    use filtering_tools
    use datalayout
    use differential_tools , only : computeDerivationField
    use conditional_mean_tools
    use fileio 

    ! Input/output
    logical, intent(out)                    :: res
    integer, intent(in)                     :: spec_rank,nbcpus,iFic
    type(REAL_DATA_LAYOUT), intent(in)      :: U,V,W

    ! Local field
    type(WAVENUMBERS)                       :: WN
    integer                                 :: filter,nd,file_id
    real(WP)                                :: delta
    logical                                 :: success,existe
    type(COMPLEX_DATA_LAYOUT)               :: T1RGMK,T2RGMK,T3RGMK
    type(COMPLEX_DATA_LAYOUT)               :: Ukfil,Vkfil,Wkfil 
    type(REAL_DATA_LAYOUT)                  :: T1GM,T2GM,T3GM 
    type(REAL_DATA_LAYOUT)                  :: T1DNS,T2DNS,T3DNS
    type(REAL_DATA_LAYOUT)                  :: T1RGM,T2RGM,T3RGM
    type(REAL_DATA_LAYOUT)                  :: Ufil,Vfil,Wfil
    type(REAL_DATA_LAYOUT)                  :: dudx,dudy,dudz 
    type(REAL_DATA_LAYOUT)                  :: S11,S12,S13
    type(REAL_DATA_LAYOUT)                  :: S22,S23,S33,norm
    type(REAL_DATA_LAYOUT)                  :: nrjDiss,nrjDissExact
    INTEGER                                 :: bin, typeDisc
    REAL(WP)                                :: mseQuadGM1,mseQuadGM2,mseQuadGM3
    REAL(WP)                                :: mseIrrGM1,mseIrrGM2,mseIrrGM3
    REAL(WP)                                :: mseQuadSmagDyn1,mseQuadSmagDyn2,mseQuadSmagDyn3
    REAL(WP)                                :: mseQuadSmagStat1,mseQuadSmagStat2,mseQuadSmagStat3
    REAL(WP)                                :: mseIrrSmag1,mseIrrSmag2,mseIrrSmag3
    REAL(WP)                                :: mseIrrSikmoinsSkj1,mseIrrSikmoinsSkj2,mseIrrSikmoinsSkj3
    REAL(WP)                                :: mseQuadSikmoinsSkj1,mseQuadSikmoinsSkj2,mseQuadSikmoinsSkj3
    REAL(WP)                                :: mseRGMSikmoinsSkjDyn1,mseRGMSikmoinsSkjDyn2,mseRGMSikmoinsSkjDyn3 
    REAL(WP)                                :: dissExact,dissIrrSikmoinsSkj
    REAL(WP)                                :: dissSmagDyna,dissSmagIrr,dissSmagStat
    REAL(WP)                                :: mseRGMSikmoinsSkjDynFab1,mseRGMSikmoinsSkjDynFab2,mseRGMSikmoinsSkjDynFab3
    REAL(WP)                                :: dissIrrGrad,dissQuadGrad,dissQuadSikmoinsSkjDyn,dissQuadSikmoinsSkjDynFab,dissQuadSikmoinsSkjStat
    REAL(WP)                                :: mseIrrSikmoinsSkjOikOjk1,mseIrrSikmoinsSkjOikOjk2,mseIrrSikmoinsSkjOikOjk3,dissIrrSikmoinsSkjOikOjk
    REAL(WP)                                :: mseQuadSikmoinsSkjOikOjk1,mseQuadSikmoinsSkjOikOjk2,mseQuadSikmoinsSkjOikOjk3,dissQuadSikmoinsSkjOikOjk
    REAL(WP)                                :: coefExactDiss,mseQuadSikmoinsSkjCoefDissExact1,mseQuadSikmoinsSkjCoefDissExact2,mseQuadSikmoinsSkjCoefDissExact3
    REAL(WP)                                :: dissQuadSikmoinsSkjOikOjkDynFab,dissQuadSikmoinsSkjOikOjkDyn2Fab
    REAL(WP)                                :: dissRGMcoefDissExact,mseQuadSikmoinsSkjOikOjkDynFab1,mseQuadSikmoinsSkjOikOjkDynFab2,mseQuadSikmoinsSkjOikOjkDynFab3  
    REAL(WP)                                :: mseQuadSikmoinsSkjOikOjkDyn2Fab1,mseQuadSikmoinsSkjOikOjkDyn2Fab2,mseQuadSikmoinsSkjOikOjkDyn2Fab3  
    REAL(WP)                                :: mseQuadSikmoinsSkjDynYouMoin1,mseQuadSikmoinsSkjDynYouMoin2,mseQuadSikmoinsSkjDynYouMoin3,dissQuadSikmoinsSkjDynYouMoin
    REAL(WP)                                :: mseQuadDissGMOpt,mseQuadDissGM,mseQuadDissSmag,mseQuadDissSmagOpt,mseQuadDissSmagDyn,mseQuadDissSikmoinsSkjOpt
    REAL(WP)                                :: mseQuadDissSikmoinsSkjStat,mseQuadDissSikmoinsSkjDyn,mseQuadDissSikmoinsSkjDynFab,mseQuadDissSikmoinsSkjDynYouMoin

    ! ===== Initialisation =====
    res = .false.
    if (.not. initDataLayout("T1RGMK",T1RGMK,(U%nx/2)+1,U%ny,U%nz, &
       & U%Lx,U%Ly,U%Lz,nbcpus,spec_rank,alongZ)) then 
      write(6,'(a,i0,a)')'[ERROR] Not able to create Exk'
      return
    endif
    if (.not. copystructonly(T1RGMK,T2RGMK) ) return
    if (.not. copystructonly(T1RGMK,T3RGMK) ) return
    if (.not. copystructonly(U,T1GM) ) return
    if (.not. copystructonly(U,T2GM) ) return
    if (.not. copystructonly(U,T3GM) ) return
    if (.not. copystructonly(U,T1RGM) ) return
    if (.not. copystructonly(U,T2RGM) ) return
    if (.not. copystructonly(U,T3RGM) ) return
    if (.not. copystructonly(U,T1DNS) ) return
    if (.not. copystructonly(U,T2DNS) ) return
    if (.not. copystructonly(U,T3DNS) ) return
    if (.not. copystructonly(U,S11) ) return
    if (.not. copystructonly(U,S12) ) return
    if (.not. copystructonly(U,S13) ) return
    if (.not. copystructonly(U,S22) ) return
    if (.not. copystructonly(U,S23) ) return
    if (.not. copystructonly(U,S33) ) return
    if (.not. copystructonly(U,norm) ) return
    if (.not. copystructonly(U,Ufil) ) return
    if (.not. copystructonly(U,Vfil) ) return
    if (.not. copystructonly(U,Wfil) ) return
    if (.not. copystructonly(U,dudx) ) return
    if (.not. copystructonly(U,dudy) ) return
    if (.not. copystructonly(U,dudz) ) return
    if (.not. copystructonly(T1RGMK,Ukfil) ) return
    if (.not. copystructonly(T1RGMK,Vkfil) ) return
    if (.not. copystructonly(T1RGMK,Wkfil) ) return
    if (.not. copystructonly(U,nrjDiss) ) return
    if (.not. copystructonly(U,nrjDissExact) ) return
    !Calcul des nombres d'onde
    if (.not. initWN(WN,spec_rank,U%nx,U%ny,U%nz) .or. &
       &.not.  computeWN(WN,spec_rank,U%Lx,U%Ly,U%Lz,U%nx,U%ny,U%nz) ) then 
      write(6,'(a,i0,a)')'[ERROR] Cannot compute WN'
      return
    endif
    ! cutoff -> 1
    ! box -> 2
    filter = 1
    if ( parser_is_defined('Filter type for post process')) then
      if (.not. parseFilter(filter,spec_rank,other='Filter type for post process')) return
    endif
    bin=250
    typeDisc=1
    if (.not. diff_tool_initBuff(T1RGMK)) return
    do nd=2,U%nx/16,2
    !do nd=1,2
      delta=nd*computeDelta(U)
      call computeFilter(WN,delta,U,Ukfil,filter)
      call computeFilter(WN,delta,V,Vkfil,filter)
      call computeFilter(WN,delta,W,Wkfil,filter)
      IF (.NOT. showDivU(Ukfil,Vkfil,Wkfil,WN,spec_rank)) return
      IF (.NOT. NullifyDiv(Ukfil,Vkfil,Wkfil, WN, spec_rank)) return 
      IF (.NOT. showDivU(Ukfil,Vkfil,Wkfil,WN,spec_rank)) return
      call btran(Ukfil,Ufil,res)
      call btran(Vkfil,Vfil,res)
      call btran(Wkfil,Wfil,res)
      !DNS
      call computeT_ijVecA(U,V,W,Ufil,Vfil,Wfil,WN,T1RGM,T2RGM,T3RGM,T1GM,T2GM,T3GM,filter,success,delta)
      call ComputeDivSymTijVec1Vec2(T1RGMK,T2RGMK,T3RGMK,T2RGM,T3RGM,T2GM,T1RGM,T1GM,T3GM,WN,res)
      call btran(T1RGMK,T1DNS,success)
      if (.not. success) return
      call btran(T2RGMK,T2DNS,success)
      if (.not. success) return
      call btran(T3RGMK,T3DNS,success)
      if (.not. success) return
      nrjDissExact%name='dissipationExact'
      nrjDissExact%values = - ( T1DNS%values * Ufil%values + T2DNS%values * Vfil%values + T3DNS%values * Wfil%values )
      dissExact= computeFieldAvg(nrjDiss,spec_rank)
      !GM
      !Calcul des erreurs irr et quad GM
      call gradientVelDIV(T1RGMK,T2RGMK,T3RGMK,Ufil,Vfil,Wfil,WN,delta,spec_rank,nbcpus)
      call btran(T1RGMK,T1GM,success)
      if (.not. success) return
      call btran(T2RGMK,T2GM,success)
      if (.not. success) return
      call btran(T3RGMK,T3GM,success)
      if (.not. success) return
      if (.not. computeEONVar( spec_rank,(/T1GM/),T1DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudx) ) return 
      if (.not. computeEONVar( spec_rank,(/T2GM/),T2DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudy) ) return 
      if (.not. computeEONVar( spec_rank,(/T3GM/),T3DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudz) ) return 
      mseIrrGM1  = computeMSE(T1DNS,dudx,spec_rank,.true.)
      mseIrrGM2  = computeMSE(T2DNS,dudy,spec_rank,.true.)
      mseIrrGM3  = computeMSE(T3DNS,dudz,spec_rank,.true.)
      nrjDiss%values = - ( dudx%values * Ufil%values + dudy%values * Vfil%values + dudz%values * Wfil%values )
      dissIrrGrad= computeFieldAvg(nrjDiss,spec_rank)
      mseQuadDissGMOpt = computeMSE(nrjDissExact,nrjDiss,spec_rank,.true.)
      nrjDiss%name="GMOpt"
      if (.not. computeTwoFieldsPDF(nd,nrjDissExact,nrjDiss,(/300,300/),spec_rank,100)) return
      mseQuadGM1 = computeMSE(T1DNS,T1GM,spec_rank,.true.)
      mseQuadGM2 = computeMSE(T2DNS,T2GM,spec_rank,.true.)
      mseQuadGM3 = computeMSE(T3DNS,T3GM,spec_rank,.true.)
      nrjDiss%values = - ( T1GM%values * Ufil%values + T2GM%values * Vfil%values + T3GM%values * Wfil%values )
      mseQuadDissGM = computeMSE(nrjDissExact,nrjDiss,spec_rank,.true.)
      nrjDiss%name="GM"
      if (.not. computeTwoFieldsPDF(nd,nrjDissExact,nrjDiss,(/300,300/),spec_rank,100)) return
      dissQuadGrad = computeFieldAvg(nrjDiss,spec_rank)
      !Calcul erreur irr et quad Smag
      !Calcul de Smag statique pour l'erreur irr
      call smagorinskyVel(T1RGMK,T2RGMK,T3RGMK,Ufil,Vfil,Wfil,Ukfil,Vkfil,Wkfil,WN,1,spec_rank,res,delta=delta)
      if (.not. res) return
      call btran(T1RGMK,T1GM,success)
      if (.not. success) return
      call btran(T2RGMK,T2GM,success)
      if (.not. success) return
      call btran(T3RGMK,T3GM,success)
      if (.not. success) return
      mseQuadSmagStat1 = computeMSE(T1DNS,T1GM,spec_rank,.true.)
      mseQuadSmagStat2 = computeMSE(T2DNS,T2GM,spec_rank,.true.)
      mseQuadSmagStat3 = computeMSE(T3DNS,T3GM,spec_rank,.true.)
      nrjDiss%values = - ( T1GM%values * Ufil%values + T2GM%values * Vfil%values + T3GM%values * Wfil%values )
      nrjDiss%name="Smag"
      if (.not. computeTwoFieldsPDF(nd,nrjDissExact,nrjDiss,(/300,300/),spec_rank,100)) return
      mseQuadDissSmag = computeMSE(nrjDissExact,nrjDiss,spec_rank,.true.)
      dissSmagStat = computeFieldAvg(nrjDiss,spec_rank)
      if (.not. computeEONVar( spec_rank,(/T1GM/),T1DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudx) ) return 
      if (.not. computeEONVar( spec_rank,(/T2GM/),T2DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudy) ) return 
      if (.not. computeEONVar( spec_rank,(/T3GM/),T3DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudz) ) return 
      mseIrrSmag1  = computeMSE(T1DNS,dudx,spec_rank,.true.)
      mseIrrSmag2  = computeMSE(T2DNS,dudy,spec_rank,.true.)
      mseIrrSmag3  = computeMSE(T3DNS,dudz,spec_rank,.true.)
      nrjDiss%values = - ( dudx%values * Ufil%values + dudy%values * Vfil%values + dudz%values * Wfil%values )
      nrjDiss%name="SmagOpt"
      if (.not. computeTwoFieldsPDF(nd,nrjDissExact,nrjDiss,(/300,300/),spec_rank,100)) return
      mseQuadDissSmagOpt = computeMSE(nrjDissExact,nrjDiss,spec_rank,.true.)
      dissSmagIrr = computeFieldAvg(nrjDiss,spec_rank)
      !Calcul de Smag dynamique pour l'erreur quad 
      call smagorinskyVel(T1RGMK,T2RGMK,T3RGMK,Ufil,Vfil,Wfil,Ukfil,Vkfil,Wkfil,WN,2,spec_rank,res,numVelFilter=1,type_avg=0,delta=delta)
      if (.not. res ) return
      call btran(T1RGMK,T1GM,success)
      if (.not. success) return
      call btran(T2RGMK,T2GM,success)
      if (.not. success) return
      call btran(T3RGMK,T3GM,success)
      if (.not. success) return
      mseQuadSmagDyn1 = computeMSE(T1DNS,T1GM,spec_rank,.true.)
      mseQuadSmagDyn2 = computeMSE(T2DNS,T2GM,spec_rank,.true.)
      mseQuadSmagDyn3 = computeMSE(T3DNS,T3GM,spec_rank,.true.)
      nrjDiss%values = - ( T1GM%values * Ufil%values + T2GM%values * Vfil%values + T3GM%values * Wfil%values )
      dissSmagDyna = computeFieldAvg(nrjDiss,spec_rank)
      mseQuadDissSmagDyn= computeMSE(nrjDissExact,nrjDiss,spec_rank,.true.)
      nrjDiss%name="SmagDyn"
      if (.not. computeTwoFieldsPDF(nd,nrjDissExact,nrjDiss,(/300,300/),spec_rank,100)) return
      !Calcul de irr et quad RGMSIKmoinsSKJ
      call RGMSikmoinsSkj(T1RGMK,T2RGMK,T3RGMK,Ufil,Ukfil,Vkfil,Wkfil,WN,spec_rank,delta=delta)
      call btran(T1RGMK,T1GM,success)
      if (.not. success) return
      call btran(T2RGMK,T2GM,success)
      if (.not. success) return
      call btran(T3RGMK,T3GM,success)
      if (.not. success) return
      if (.not. computeEONVar( spec_rank,(/T1GM/),T1DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudx) ) return 
      if (.not. computeEONVar( spec_rank,(/T2GM/),T2DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudy) ) return 
      if (.not. computeEONVar( spec_rank,(/T3GM/),T3DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudz) ) return 
      mseIrrSikmoinsSkj1 = computeMSE(T1DNS,dudx,spec_rank,.true.)
      mseIrrSikmoinsSkj2 = computeMSE(T2DNS,dudy,spec_rank,.true.)
      mseIrrSikmoinsSkj3 = computeMSE(T3DNS,dudz,spec_rank,.true.)
      nrjDiss%values = - ( dudx%values * Ufil%values + dudy%values * Vfil%values + dudz%values * Wfil%values )
      mseQuadDissSikmoinsSkjOpt= computeMSE(nrjDissExact,nrjDiss,spec_rank,.true.)
      nrjDiss%name="RGMSik-Sjk_opt"
      if (.not. computeTwoFieldsPDF(nd,nrjDissExact,nrjDiss,(/300,300/),spec_rank,100)) return
      dissIrrSikmoinsSkj = computeFieldAvg(nrjDiss,spec_rank)
      mseQuadSikmoinsSkj1 = computeMSE(T1DNS,T1GM,spec_rank,.true.)
      mseQuadSikmoinsSkj2 = computeMSE(T2DNS,T2GM,spec_rank,.true.)
      mseQuadSikmoinsSkj3 = computeMSE(T3DNS,T3GM,spec_rank,.true.)
      nrjDiss%values = - ( T1GM%values * Ufil%values + T2GM%values * Vfil%values + T3GM%values * Wfil%values )
      dissQuadSikmoinsSkjStat= computeFieldAvg(nrjDiss,spec_rank)
      mseQuadDissSikmoinsSkjStat= computeMSE(nrjDissExact,nrjDiss,spec_rank,.true.)
      nrjDiss%name="RGMSik-Sjk_stat"
      if (.not. computeTwoFieldsPDF(nd,nrjDissExact,nrjDiss,(/300,300/),spec_rank,100)) return
      !Calcul de quad RGMSikmoinsSjkDyn  Procédure dynamique usuelle
      call RGMSikmoinsSkjDyn(T1RGMK,T2RGMK,T3RGMK,Ufil,Vfil,Wfil,WN,spec_rank,filter,delta=delta)
      call btran(T1RGMK,T1GM,success)
      if (.not. success) return
      call btran(T2RGMK,T2GM,success)
      if (.not. success) return
      call btran(T3RGMK,T3GM,success)
      if (.not. success) return
      mseRGMSikmoinsSkjDyn1 = computeMSE(T1DNS,T1GM,spec_rank,.true.)
      mseRGMSikmoinsSkjDyn2 = computeMSE(T2DNS,T2GM,spec_rank,.true.)
      mseRGMSikmoinsSkjDyn3 = computeMSE(T3DNS,T3GM,spec_rank,.true.)
      nrjDiss%values = - ( T1GM%values * Ufil%values + T2GM%values * Vfil%values + T3GM%values * Wfil%values )
      dissQuadSikmoinsSkjDyn = computeFieldAvg(nrjDiss,spec_rank)
      mseQuadDissSikmoinsSkjDyn= computeMSE(nrjDissExact,nrjDiss,spec_rank,.true.)
      nrjDiss%name="RGMSik-Sjk_dyn"
      if (.not. computeTwoFieldsPDF(nd,nrjDissExact,nrjDiss,(/300,300/),spec_rank,100)) return
      !Calcul de quad RGMSikmoinsSkjDynFab Procédure dynamique sans la partie sur-filtrée
      call RGMSikmoinsSkjDynFab(T1RGMK,T2RGMK,T3RGMK,Ufil,Ukfil,Vkfil,Wkfil,WN,spec_rank,filter,delta=delta)
      call btran(T1RGMK,T1GM,success)
      if (.not. success) return
      call btran(T2RGMK,T2GM,success)
      if (.not. success) return
      call btran(T3RGMK,T3GM,success)
      if (.not. success) return
      mseRGMSikmoinsSkjDynFab1 = computeMSE(T1DNS,T1GM,spec_rank,.true.)
      mseRGMSikmoinsSkjDynFab2 = computeMSE(T2DNS,T2GM,spec_rank,.true.)
      mseRGMSikmoinsSkjDynFab3 = computeMSE(T3DNS,T3GM,spec_rank,.true.)
      nrjDiss%values = - ( T1GM%values * Ufil%values + T2GM%values * Vfil%values + T3GM%values * Wfil%values )
      dissQuadSikmoinsSkjDynFab = computeFieldAvg(nrjDiss,spec_rank)
      mseQuadDissSikmoinsSkjDynFab= computeMSE(nrjDissExact,nrjDiss,spec_rank,.true.)
      nrjDiss%name="RGMSik-Sjk_dynFab"
      if (.not. computeTwoFieldsPDF(nd,nrjDissExact,nrjDiss,(/300,300/),spec_rank,100)) return
      !Calcul de irr et quad RGMSIKmoinsSKJOikOjk
      call RGMSikmoinsSkjOikOjk(T1RGMK,T2RGMK,T3RGMK,Ufil,Ukfil,Vkfil,Wkfil,WN,spec_rank,delta=delta)
      call btran(T1RGMK,T1GM,success)
      if (.not. success) return
      call btran(T2RGMK,T2GM,success)
      if (.not. success) return
      call btran(T3RGMK,T3GM,success)
      if (.not. success) return
      if (.not. computeEONVar( spec_rank,(/T1GM/),T1DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudx) ) return 
      if (.not. computeEONVar( spec_rank,(/T2GM/),T2DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudy) ) return 
      if (.not. computeEONVar( spec_rank,(/T3GM/),T3DNS,&
                      & (/typeDisc/),(/bin/),&
                      & fieldOpt = dudz) ) return 
      mseIrrSikmoinsSkjOikOjk1 = computeMSE(T1DNS,dudx,spec_rank,.true.)
      mseIrrSikmoinsSkjOikOjk2 = computeMSE(T2DNS,dudy,spec_rank,.true.)
      mseIrrSikmoinsSkjOikOjk3 = computeMSE(T3DNS,dudz,spec_rank,.true.)
      nrjDiss%values = - ( dudx%values * Ufil%values + dudy%values * Vfil%values + dudz%values * Wfil%values )
      dissIrrSikmoinsSkjOikOjk = computeFieldAvg(nrjDiss,spec_rank)
      mseQuadSikmoinsSkjOikOjk1 = computeMSE(T1DNS,T1GM,spec_rank,.true.)
      mseQuadSikmoinsSkjOikOjk2 = computeMSE(T2DNS,T2GM,spec_rank,.true.)
      mseQuadSikmoinsSkjOikOjk3 = computeMSE(T3DNS,T3GM,spec_rank,.true.)
      nrjDiss%values = - ( T1GM%values * Ufil%values + T2GM%values * Vfil%values + T3GM%values * Wfil%values )
      dissQuadSikmoinsSkjOikOjk= computeFieldAvg(nrjDiss,spec_rank)
      !Calcul de irr et quad RGMSIKmoinsSKJOikOjkdynFab
      call RGMSikmoinsSkjOikOjkDynFab(T1RGMK,T2RGMK,T3RGMK,Ufil,Vfil,Wfil,WN,spec_rank,filter,delta=delta)
      call btran(T1RGMK,T1GM,success)
      if (.not. success) return
      call btran(T2RGMK,T2GM,success)
      if (.not. success) return
      call btran(T3RGMK,T3GM,success)
      if (.not. success) return
      mseQuadSikmoinsSkjOikOjkDynFab1 = computeMSE(T1DNS,T1GM,spec_rank,.true.)
      mseQuadSikmoinsSkjOikOjkDynFab2 = computeMSE(T2DNS,T2GM,spec_rank,.true.)
      mseQuadSikmoinsSkjOikOjkDynFab3 = computeMSE(T3DNS,T3GM,spec_rank,.true.)
      nrjDiss%values = - ( T1GM%values * Ufil%values + T2GM%values * Vfil%values + T3GM%values * Wfil%values )
      dissQuadSikmoinsSkjOikOjkDynFab= computeFieldAvg(nrjDiss,spec_rank)
      !Calcul de irr et quad RGMSIKmoinsSKJOikOjkdyn2Fab
      call RGMSikmoinsSkjOikOjkDyn2Fab(T1RGMK,T2RGMK,T3RGMK,Ufil,Vfil,Wfil,WN,spec_rank,filter,delta=delta)
      call btran(T1RGMK,T1GM,success)
      if (.not. success) return
      call btran(T2RGMK,T2GM,success)
      if (.not. success) return
      call btran(T3RGMK,T3GM,success)
      if (.not. success) return
      mseQuadSikmoinsSkjOikOjkDyn2Fab1 = computeMSE(T1DNS,T1GM,spec_rank,.true.)
      mseQuadSikmoinsSkjOikOjkDyn2Fab2 = computeMSE(T2DNS,T2GM,spec_rank,.true.)
      mseQuadSikmoinsSkjOikOjkDyn2Fab3 = computeMSE(T3DNS,T3GM,spec_rank,.true.)
      nrjDiss%values = - ( T1GM%values * Ufil%values + T2GM%values * Vfil%values + T3GM%values * Wfil%values )
      dissQuadSikmoinsSkjOikOjkDyn2Fab= computeFieldAvg(nrjDiss,spec_rank)
      !Calcul de irr et quad RGMSIKmoinsSKJOikOjkdyn2Fab
      call RGMSikmoinsSkjDynYouMoin(T1RGMK,T2RGMK,T3RGMK,Ufil,Vfil,Wfil,Ukfil,Vkfil,Wkfil,WN,spec_rank,filter,delta=delta)
      call btran(T1RGMK,T1GM,success)
      if (.not. success) return
      call btran(T2RGMK,T2GM,success)
      if (.not. success) return
      call btran(T3RGMK,T3GM,success)
      if (.not. success) return
      mseQuadSikmoinsSkjDynYouMoin1 = computeMSE(T1DNS,T1GM,spec_rank,.true.)
      mseQuadSikmoinsSkjDynYouMoin2 = computeMSE(T2DNS,T2GM,spec_rank,.true.)
      mseQuadSikmoinsSkjDynYouMoin3 = computeMSE(T3DNS,T3GM,spec_rank,.true.)
      nrjDiss%values = - ( T1GM%values * Ufil%values + T2GM%values * Vfil%values + T3GM%values * Wfil%values )
      dissQuadSikmoinsSkjDynYouMoin = computeFieldAvg(nrjDiss,spec_rank)
      mseQuadDissSikmoinsSkjDynYouMoin= computeMSE(nrjDissExact,nrjDiss,spec_rank,.true.)
      nrjDiss%name="RGMSik-Sjk_dynYouMoin"
      if (.not. computeTwoFieldsPDF(nd,nrjDissExact,nrjDiss,(/300,300/),spec_rank,100)) return

      !Test coefs
      call RGMSikmoinsSkj(T1RGMK,T2RGMK,T3RGMK,Ufil,Ukfil,Vkfil,Wkfil,WN,spec_rank,delta=delta)
      call btran(T1RGMK,T1GM,success)
      if (.not. success) return
      call btran(T2RGMK,T2GM,success)
      if (.not. success) return
      call btran(T3RGMK,T3GM,success)
      if (.not. success) return
      nrjDiss%values = - ( T1GM%values * Ufil%values + T2GM%values * Vfil%values + T3GM%values * Wfil%values )
      coefExactDiss = dissExact / computeFieldAvg(nrjDiss,spec_rank) 
      if (spec_rank .EQ. 0) write (6,'(a,1x,g15.8)') '[INFO coefExactDiss RGM]  coef dyn rgm       :',coefExactDiss
      call RGMSikmoinsSkj(T1RGMK,T2RGMK,T3RGMK,Ufil,Ukfil,Vkfil,Wkfil,WN,spec_rank,delta=delta)
      call btran(T1RGMK,T1GM,success)
      if (.not. success) return
      call btran(T2RGMK,T2GM,success)
      if (.not. success) return
      call btran(T3RGMK,T3GM,success)
      if (.not. success) return
      mseQuadSikmoinsSkjCoefDissExact1 = computeMSE(T1DNS,T1GM,spec_rank,.true.)
      mseQuadSikmoinsSkjCoefDissExact2 = computeMSE(T2DNS,T2GM,spec_rank,.true.)
      mseQuadSikmoinsSkjCoefDissExact3 = computeMSE(T3DNS,T3GM,spec_rank,.true.)
      nrjDiss%values = - ( T1GM%values * Ufil%values + T2GM%values * Vfil%values + T3GM%values * Wfil%values )
      dissRGMcoefDissExact = computeFieldAvg(nrjDiss,spec_rank) 

      !dump
      if (spec_rank .EQ. 0) write (6,'(a)')                         '################    ERREUR IRREDUCTIBLE    #######################' 
      if (spec_rank .EQ. 0) write (6,'(a)')                         '###########           Modele de base       #######################' 
      if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO IRR]  GM de base                  :',nd,mseIrrGM1,mseIrrGM2,mseIrrGM3
      if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO IRR]  SMAG                        :',nd,mseIrrSmag1,mseIrrSmag2,mseIrrSmag3
      if (spec_rank .EQ. 0) write (6,'(a)')                         '###########       Modele sur t_ij          #######################' 
      if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO IRR] Sik-Skj                      :',nd,mseIrrSikmoinsSkj1,mseIrrSikmoinsSkj2,mseIrrSikmoinsSkj3
      if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO IRR] Sik-Skj + OikOjk             :',nd,mseIrrSikmoinsSkjOikOjk1,mseIrrSikmoinsSkjOikOjk2,mseIrrSikmoinsSkjOikOjk3
      if (spec_rank .EQ. 0) write (6,'(a)')                         '################    ERREUR QUADRATIQUE     #######################' 
      if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO QUAD] GM                          :',nd,mseQuadGM1,mseQuadGM2,mseQuadGM3
      if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO QUAD] SMAG dyn                    :',nd,mseQuadSmagDyn1,mseQuadSmagDyn2,mseQuadSmagDyn3
      if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO QUAD] SMAG stat                   :',nd,mseQuadSmagStat1,mseQuadSmagStat2,mseQuadSmagStat3
      if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO QUAD] Sik-Skj    Statique         :',nd,mseQuadSikmoinsSkj1,mseQuadSikmoinsSkj2,mseQuadSikmoinsSkj3
      if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO QUAD] Sik-Skj    Dyn              :',nd,mseRGMSikmoinsSkjDyn1,mseRGMSikmoinsSkjDyn2,mseRGMSikmoinsSkjDyn3
      if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO QUAD] Sik-Skj    DynFab           :',nd,mseRGMSikmoinsSkjDynFab1,mseRGMSikmoinsSkjDynFab2,mseRGMSikmoinsSkjDynFab3
     if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO QUAD] Sik-Skj + OikOjk            :',nd,mseQuadSikmoinsSkjOikOjk1,mseQuadSikmoinsSkjOikOjk2,mseQuadSikmoinsSkjOikOjk3
     if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO QUAD] Sik-Skj  Dyncoef diss exact :',nd,mseQuadSikmoinsSkjCoefDissExact1,mseQuadSikmoinsSkjCoefDissExact2,mseQuadSikmoinsSkjCoefDissExact3
     if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO QUAD] Sik-Skj + OikOjk  dynFab    :',nd,mseQuadSikmoinsSkjOikOjkDynFab1,mseQuadSikmoinsSkjOikOjkDynFab2,mseQuadSikmoinsSkjOikOjkDynFab3
     if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO QUAD] Sik-Skj + OikOjk  dyn2Fab   :',nd,mseQuadSikmoinsSkjOikOjkDyn2Fab1,mseQuadSikmoinsSkjOikOjkDyn2Fab2,mseQuadSikmoinsSkjOikOjkDyn2Fab3
     if (spec_rank .EQ. 0) write (6,'(1(a,1x,i0,1x,3(g15.8,1x)))') '[INFO QUAD] Sik-Skj  dyn You and Moin   :',nd,mseQuadSikmoinsSkjDynYouMoin1,mseQuadSikmoinsSkjDynYouMoin2,mseQuadSikmoinsSkjDynYouMoin3
      if (spec_rank .EQ. 0) then
        file_id=iopen()
        inquire(file='bilanMSEIRR.table', exist=existe)
        open(unit=file_id, file='bilanMSEIRR.table', form="FORMATTED",position="append")
        if ( .not. existe) then
          write(file_id, '(a)') '#iFic nd mseIrrGM1 mseIrrGM2 mseIrrGM3 &
                                         &mseIrrSmag1 mseIrrSmag2 mseIrrSmag3 &
                                         &mseIrrSikmoinsSkj1 mseIrrSikmoinsSkj2 mseIrrSikmoinsSkj3 &
                                         &mseIrrSikmoinsSkjOikOjk1 mseIrrSikmoinsSkjOikOjk2 mseIrrSikmoinsSkjOikOjk3'
        endif
        write(file_id,'(i0,1x,i0,1x,12(g15.8,1x))') iFic,nd,mseIrrGM1,mseIrrGM2,mseIrrGM3,&
                                                          &mseIrrSmag1,mseIrrSmag2,mseIrrSmag3,&
                                                          &mseIrrSikmoinsSkj1,mseIrrSikmoinsSkj2,mseIrrSikmoinsSkj3,&
                                                          &mseIrrSikmoinsSkjOikOjk1,mseIrrSikmoinsSkjOikOjk2,mseIrrSikmoinsSkjOikOjk3
        close(file_id)
        file_id=iclose(file_id)

        file_id=iopen()
        inquire(file='bilanMSEQUAD.table', exist=existe)
        open(unit=file_id, file='bilanMSEQUAD.table', form="FORMATTED",position="append")
        if ( .not. existe) then
          write(file_id, '(a)') '#iFic nd mseQuadGM1 mseQuadGM2 mseQuadGM3 &
                                         &mseQuadSmagDyn1 mseQuadSmagDyn2 mseQuadSmagDyn3 &
                                         &mseQuadSmagStat1 mseQuadSmagStat2 mseQuadSmagStat3 &
                                         &mseQuadSikmoinsSkj1 mseQuadSikmoinsSkj2 mseQuadSikmoinsSkj3 &
                                         &mseRGMSikmoinsSkjDyn1 mseRGMSikmoinsSkjDyn2 mseRGMSikmoinsSkjDyn3 &
                                         &mseRGMSikmoinsSkjDynFab1 mseRGMSikmoinsSkjDynFab2 mseRGMSikmoinsSkjDynFab3 &
                                         &mseQuadSikmoinsSkjOikOjk1 mseQuadSikmoinsSkjOikOjk2 mseQuadSikmoinsSkjOikOjk3 &
                                         &mseQuadSikmoinsSkjCoefDissExact1 mseQuadSikmoinsSkjCoefDissExact2 mseQuadSikmoinsSkjCoefDissExact3 &
                                         &mseQuadSikmoinsSkjOikOjkDynFab1 mseQuadSikmoinsSkjOikOjkDynFab2 mseQuadSikmoinsSkjOikOjkDynFab3 &
                                         &mseQuadSikmoinsSkjOikOjkDyn2Fab1 mseQuadSikmoinsSkjOikOjkDyn2Fab2 mseQuadSikmoinsSkjOikOjkDyn2Fab3 &
                                         &mseQuadSikmoinsSkjDynYouMoin1 mseQuadSikmoinsSkjDynYouMoin2 mseQuadSikmoinsSkjDynYouMoin3'
        endif
        write(file_id,'(i0,1x,i0,1x,33(g15.8,1x))') iFic,nd,mseQuadGM1,mseQuadGM2,mseQuadGM3,&
                                                          &mseQuadSmagDyn1,mseQuadSmagDyn2,mseQuadSmagDyn3,&
                                                          &mseQuadSmagStat1,mseQuadSmagStat2,mseQuadSmagStat3,&
                                                          &mseQuadSikmoinsSkj1,mseQuadSikmoinsSkj2,mseQuadSikmoinsSkj3,&
                                              &mseRGMSikmoinsSkjDyn1,mseRGMSikmoinsSkjDyn2,mseRGMSikmoinsSkjDyn3,&
                                              &mseRGMSikmoinsSkjDynFab1,mseRGMSikmoinsSkjDynFab2,mseRGMSikmoinsSkjDynFab3,&
                                              &mseQuadSikmoinsSkjOikOjk1,mseQuadSikmoinsSkjOikOjk2,mseQuadSikmoinsSkjOikOjk3,&
                                         &mseQuadSikmoinsSkjCoefDissExact1,mseQuadSikmoinsSkjCoefDissExact2,mseQuadSikmoinsSkjCoefDissExact3,&
                                         &mseQuadSikmoinsSkjOikOjkDynFab1,mseQuadSikmoinsSkjOikOjkDynFab2,mseQuadSikmoinsSkjOikOjkDynFab3,&
                                         &mseQuadSikmoinsSkjOikOjkDyn2Fab1,mseQuadSikmoinsSkjOikOjkDyn2Fab2,mseQuadSikmoinsSkjOikOjkDyn2Fab3,&
                                         &mseQuadSikmoinsSkjDynYouMoin1,mseQuadSikmoinsSkjDynYouMoin2,mseQuadSikmoinsSkjDynYouMoin3
        close(file_id)
        file_id=iclose(file_id)

        file_id=iopen()
        inquire(file='bilanDIFFDISS_GM_VS_RGM.table', exist=existe)
        open(unit=file_id, file='bilanDIFFDISS_GM_VS_RGM.table', form="FORMATTED",position="append")
        if ( .not. existe) then
          write(file_id,'(a)') '#iFic nd  dissExact dissIrrGrad dissQuadGrad &
                                &dissSmagIrr dissSmagStat dissSmagDyna &
                                &dissIrrSikmoinsSkj dissQuadSikmoinsSkjStat dissQuadSikmoinsSkjDyn dissQuadSikmoinsSkjDynFab &
                                &dissIrrSikmoinsSkjOikOjk dissQuadSikmoinsSkjOikOjk dissRGMcoefDissExact dissQuadSikmoinsSkjOikOjkDynFab dissQuadSikmoinsSkjOikOjkDyn2Fab &
                                &dissQuadSikmoinsSkjDynYouMoin'
        endif
        write(file_id,'(i0,1x,i0,1x,16(g15.8,1x))') iFic,nd,dissExact,&
                                          &dissIrrGrad,dissQuadGrad,&
                                          &dissSmagIrr,dissSmagStat,dissSmagDyna,&
                                          &dissIrrSikmoinsSkj,dissQuadSikmoinsSkjStat,dissQuadSikmoinsSkjDyn,dissQuadSikmoinsSkjDynFab,&
                                          &dissIrrSikmoinsSkjOikOjk,dissQuadSikmoinsSkjOikOjk,dissRGMcoefDissExact,&
                                          &dissQuadSikmoinsSkjOikOjkDynFab,dissQuadSikmoinsSkjOikOjkDyn2Fab,dissQuadSikmoinsSkjDynYouMoin
        close(file_id)
        file_id=iclose(file_id)

        file_id=iopen()
        inquire(file='bilanDISS_MSE.table', exist=existe)
        open(unit=file_id, file='bilanDISS_MSE.table', form="FORMATTED",position="append")
        if ( .not. existe) then
          write(file_id,'(a)') '#iFic nd mseQuadDissGMOpt mseQuadDissGM mseQuadDissSmag mseQuadDissSmagOpt & 
               &mseQuadDissSmagDyn mseQuadDissSikmoinsSkjOpt mseQuadDissSikmoinsSkjStat mseQuadDissSikmoinsSkjDyn &
               &mseQuadDissSikmoinsSkjDynFab mseQuadDissSikmoinsSkjDynYouMoin'
        endif
        write(file_id,'(i0,1x,i0,1x,10(g15.8,1x))') iFic,nd,mseQuadDissGMOpt,mseQuadDissGM,mseQuadDissSmag,&
                                                    mseQuadDissSmagOpt,mseQuadDissSmagDyn,mseQuadDissSikmoinsSkjOpt,&
                                                    mseQuadDissSikmoinsSkjStat,mseQuadDissSikmoinsSkjDyn,mseQuadDissSikmoinsSkjDynFab,&
                                                    mseQuadDissSikmoinsSkjDynYouMoin
        close(file_id)
        file_id=iclose(file_id)
      endif
    enddo

    ! ===== Free memory and delete temporary fields =====
    if( .not. deleteWN(WN)) return
    call deleteDataLayout(Ufil)
    call deleteDataLayout(Vfil)
    call deleteDataLayout(Wfil)
    call deleteDataLayout(T1DNS)
    call deleteDataLayout(T2DNS)
    call deleteDataLayout(T3DNS)
    call deleteDataLayout(T1GM)
    call deleteDataLayout(T2GM)
    call deleteDataLayout(T3GM)
    call deleteDataLayout(T1RGM)
    call deleteDataLayout(T2RGM)
    call deleteDataLayout(T3RGM)
    call deleteDataLayout(T1RGMK)
    call deleteDataLayout(T2RGMK)
    call deleteDataLayout(T3RGMK)
    call deleteDataLayout(dudx)
    call deleteDataLayout(dudy)
    call deleteDataLayout(dudz)
    call deleteDataLayout(nrjDiss)
    call deleteDataLayout(nrjDissExact)
    call deleteDataLayout(Ukfil)
    call deleteDataLayout(Vkfil)
    call deleteDataLayout(Wkfil)
    call deleteDataLayout(S11)
    call deleteDataLayout(S12)
    call deleteDataLayout(S13)
    call deleteDataLayout(S22)
    call deleteDataLayout(S23)
    call deleteDataLayout(S33)
    call deleteDataLayout(norm)
    res = .true.

end subroutine aprioriSikmoinsSkj



end module post_velocity
!> @}
