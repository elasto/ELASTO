!USEFORTEST toolbox
!USEFORTEST postprocess
!USEFORTEST avgcond 
!> @addtogroup toolbox
!! @{
!------------------------------------------------------------------------------
!
! MODULE: toolbox
!
!> @author
!> Antoine Vollant
!
! DESCRIPTION:
!> The aim of this module is to provide elementaries piece of code which can be use elsewhere.
!> BEWARE there are no systematic tests in these routines in order to not degrade performances
!>
!------------------------------------------------------------------------------


module stat_tools

 use precision_tools
 use parallel_tools
 use datalayout
 use mpilayout_tools
 use wavenumber_tools
 use transforms_tools 
! use toolbox
 
 implicit none

 public

 interface equalToVal 
   module procedure equalToVal_real_sca 
   module procedure equalToVal_real_wp_sca 
   module procedure equalToVal_wp_sca 
   module procedure equalToVal_realDataLayout 
 end interface equalToVal 

 interface computeFieldAvg2D 
   module procedure computeFieldAvg2DSelectDir
   module procedure computeFieldAvg2D
 end interface computeFieldAvg2D

 interface compute_spectrum 
   module procedure compute_spectrum_sca
   module procedure compute_spectrum_vec
   module procedure compute_spectrum_vec_k
   module procedure compute_spectrum_vecDotvec_k
 end interface compute_spectrum

 contains


!---------------------------------------------------------------------------
!> @details 
!> Compute maximum of a scalar fields.
!! @author = Antoine Vollant, LEGI
!!    @param[in]    field   = scalar field in physical space
!!    @param[in]    spec_rank   = mpi rank inside the spectral communicator
!!    @return       maxi = the maximum of scalar field
!---------------------------------------------------------------------------
function computeFieldMax(field,spec_rank) result(maxi)

     ! Input/Output
     type(REAL_DATA_LAYOUT), intent(in) :: field
     integer, intent(in) :: spec_rank
     real(WP)::maxi

     !Local data
     real(WP)::maxiL

     ! Compute
     maxiL=maxval(field%values)
     maxi =DoGlobalMax(spec_rank, maxiL)

end function computeFieldMax

!---------------------------------------------------------------------------
!> @details 
!> Compute minimum of a scalar fields.
!! @author = Antoine Vollant, LEGI
!!    @param[in]    field   = scalar field in physical space
!!    @param[in]    spec_rank   = mpi rank inside the spectral communicator
!!    @return       mini = the minimum of scalar field
!---------------------------------------------------------------------------
function computeFieldMin(field,spec_rank) result(mini)

     ! Input/Output
     type(REAL_DATA_LAYOUT), intent(in) :: field
     integer, intent(in) :: spec_rank
     real(WP)::mini

     !Local data
     real(WP)::miniL

     ! Compute
     miniL=minval( field%values )
     mini =DoGlobalMin( spec_rank, miniL )

end function computeFieldMin

!---------------------------------------------------------------------------
!> @details 
!> Compute maximum of a complex scalar fields.
!! @author = Antoine Vollant, LEGI
!!    @param[in]    field   = scalar field in physical space
!!    @param[in]    spec_rank   = mpi rank inside the spectral communicator
!!    @return       maxi = the maximum of scalar field
!---------------------------------------------------------------------------
function computeFieldCMax(field,spec_rank) result(maxi)

     ! Input/Output
     type(COMPLEX_DATA_LAYOUT), intent(in) :: field
     integer, intent(in) :: spec_rank
     real(WP)::maxi

     !Local data
     real(WP)::maxiL

     ! Compute
     maxiL=maxval(abs(field%values))
     maxi =DoGlobalMax(spec_rank, maxiL)

end function computeFieldCMax

!---------------------------------------------------------------------------
!> @details 
!> Compute minimum of a scalar fields.
!! @author = Antoine Vollant, LEGI
!!    @param[in]    field   = scalar field in physical space
!!    @param[in]    spec_rank   = mpi rank inside the spectral communicator
!!    @return       mini = the minimum of scalar field
!---------------------------------------------------------------------------
function computeFieldCMin(field,spec_rank) result(mini)

     ! Input/Output
     type(COMPLEX_DATA_LAYOUT), intent(in) :: field
     integer, intent(in) :: spec_rank
     real(WP)::mini

     !Local data
     real(WP)::miniL

     ! Compute
     miniL=minval(abs(field%values ))
     mini =DoGlobalMin( spec_rank, miniL )

end function computeFieldCMin

!---------------------------------------------------------------------------
!> @details 
!> Compute average of a scalar fields.
!! @author = Antoine Vollant, LEGI
!!    @param[in]    field   = scalar field in physical space
!!    @param[in]    spec_rank   = mpi rank inside the spectral communicator
!!    @param[in]    nbcpus   = number of processing unit 
!!    @return       avg = the average of scalar field
!---------------------------------------------------------------------------
function computeFieldAvg(field,spec_rank,nbcpus) result(avg)

     ! Input/Output
     type(REAL_DATA_LAYOUT), intent(in) :: field
     integer, intent(in)                :: spec_rank
     integer, intent(in),optional       :: nbcpus 
     real(WP)                           :: avg

     ! Local data
     real(WP)                           :: somme
!     real(WP)                           :: temp
!     integer                            :: nxLoc,nyLoc,nzLoc
!     integer                            :: i,j,k,l
!     real(WP)                           :: cpu
!     real(WP)                           :: cpu,t1,t2,t3,avg2

!Second way -> different values according to nbcpus
!-test     if (present(nbcpus)) then
!-test       cpu = nbcpus
!-test     else
!-test       cpu = getnbcpus()
!-test     endif
!-test
!-test     ! Compute
!-test     nxLoc = field%xmax - field%xmin + 1
!-test     nyLoc = field%ymax - field%ymin + 1
!-test     nzLoc = field%zmax - field%zmin + 1
!-test     somme = sum(field%values)
!-test     somme = somme / ( nxLoc * nyLoc * nzLoc )
!-test     temp=DoGlobalSum(spec_rank, somme)
!-test     avg=temp / cpu

!-test!Third Way
!-test     if (present(nbcpus)) then
!-test       cpu = nbcpus
!-test    else
!-test       cpu = getnbcpus()
!-test     endif
!-test!call cpu_time(t1)
!-test     l = 1
!-test     temp = 0.0_WP
!-test     do k=field%zmin,field%zmax
!-test       do j=field%ymin,field%ymax
!-test         do i=field%xmin,field%xmax
!-test           temp = real(l,WP)/real(l+1,WP) * temp + field%values(i,j,k)/real(l+1,WP)
!-test           l = l + 1
!-test         enddo
!-test       enddo
!-test     enddo
!-test     avg=DoGlobalSum(spec_rank, temp )/real(cpu,WP)
!call cpu_time(t2)

!First way -> overflow
       somme = sum(field%values)
       avg=DoGlobalSum(spec_rank, somme)/field%nx
       avg=avg/field%ny
       avg=avg/field%nz

!call cpu_time(t3)
!     if (spec_rank .eq. 0) write (6,'(3(a,1x,e15.8,1x,a,1x,e15.8,1x))') '    [INFO AVG] Moyenne recurrence',avg2,'temps calcul'&
!                                                                &,t2-t1,'Moyenne brut',avg,'temps calcul',t3-t2,'delta de',avg2-avg


end function computeFieldAvg

!---------------------------------------------------------------------------
!> @details 
!> Compute average of a scalar fields depending of X,Y or Z position according to dir choosen
!> Please note that size of array avg2D = number of points in domain in z direction
!> @author = Antoine Vollant, LEGI
!>    @param[in]    field   = scalar field in physical space
!>    @param[in]    spec_rank   = mpi rank inside the spectral communicator
!>    @param[in]    nbcpus   = number of processing unit
!>    @param[in]    dir = direction to compute average 
!>    @param[out]   avg2D = array of averages of scalar field
!>    @return       success = true if all work fine
!---------------------------------------------------------------------------
function computeFieldAvg2DSelectDir(field,avg2D,spec_rank,dir) result (success)

     ! Input/Output
     type(REAL_DATA_LAYOUT), intent(in)     :: field
     real(WP), dimension(:), intent(inout)  :: avg2D
     integer, intent(in)                    :: spec_rank
     integer, intent(in)                    :: dir
     logical                                :: success

     ! Local data
     real(WP), dimension(:), allocatable    :: somme
     integer                                :: cpuY, cpuZ , ires
     integer                                :: nxLoc,nyLoc,nzLoc, k!, i, j

     success = .false.

     if ( (size(avg2D,1) .ne. field%nx .and. dir .eq. 1) .or. &
        & (size(avg2D,1) .ne. field%ny .and. dir .eq. 2) .or. &
        & (size(avg2D,1) .ne. field%nz .and. dir .eq. 3)) then
       write(*,'(a)') '[ERROR] In computeField2DAvg : wrong size of avg2D'
       return
     endif
     allocate(somme(size(avg2D)),stat=ires)
     if(ires .ne. 0) then
       write(*,'(a)') '[ERROR] In computeField2DAvg : not enought memory'
       return
     end if
     somme = 0.0_WP
     call getdecompcpus(cpuY,cpuZ)
     nxLoc = field%xmax - field%xmin + 1
     nyLoc = field%ymax - field%ymin + 1
     nzLoc = field%zmax - field%zmin + 1

     !Compute average in (Y,Z) plans
     if (dir .eq. 1 ) then  
       do k = field%xmin, field%xmax
          somme(k) = sum(field%values(k,:,:))
       end do
       somme = somme / ( nzLoc * nyLoc )
       avg2D = DoGlobalVectorSum(spec_rank, somme)
       avg2D = avg2D / (cpuY*cpuZ)
     endif
     !Compute average in (X,Z) plans
     if (dir .eq. 2 ) then
       do k = field%ymin, field%ymax
          somme(k) = sum(field%values(:,k,:))
       end do
       somme = somme / ( nzLoc * nxLoc )
       avg2D = DoGlobalVectorSum(spec_rank, somme)
       avg2D = avg2D / (cpuZ)
     endif
     !Compute average in (X,Y) plans
     if (dir .eq. 3 ) then
       do k = field%zmin, field%zmax
          somme(k) = sum(field%values(:,:,k))
       end do
       somme = somme / ( nyLoc * nxLoc )
       avg2D = DoGlobalVectorSum(spec_rank, somme)
       avg2D = avg2D / (cpuY)
     endif
     ! Free memory
     deallocate(somme)

     success = .true.

end function computeFieldAvg2DSelectDir

!---------------------------------------------------------------------------
!> @details 
!> Compute average of a scalar fields depending of Z position - useful for plane
!> jet
!> Please note that size of array avg2D = number of points in domain in z direction
!> @author = Antoine Vollant, LEGI
!>    @param[in]    field   = scalar field in physical space
!>    @param[in]    spec_rank   = mpi rank inside the spectral communicator
!>    @param[in]    nbcpus   = number of processing unit
!>    @param[out]   avg = the average of scalar field
!>    @return       success = true if all work fine
!---------------------------------------------------------------------------
function computeFieldAvg2D(field,avg2D,spec_rank) result (success)

     ! Input/Output
     type(REAL_DATA_LAYOUT), intent(in)     :: field
     real(WP), dimension(:), intent(inout)  :: avg2D
     integer, intent(in)                    :: spec_rank
     logical                                :: success

     ! Local data
     real(WP), dimension(:), allocatable    :: somme
     integer                                :: cpuY, cpuZ
     integer                                :: nxLoc,nyLoc, k!, i, j

     success = .false.
     call getdecompcpus(cpuY,cpuZ)

     ! Check size
     !if ((field%zmax-field%zmin+1)> size(avg2D)) then
     if ((field%nz) .ne. size(avg2D)) then
         if(spec_rank==0) write(*,'(a)') '[ERROR] In computeField2DAvg : wrong size for avg2D'
         return
     end if

    ! Compute
     allocate(somme(size(avg2D)))
     nxLoc = field%xmax - field%xmin + 1
     nyLoc = field%ymax - field%ymin + 1
     somme = 0.0
     do k = field%zmin, field%zmax
        somme(k) = sum(field%values(:,:,k))
        !do j = field%ymin, field%ymax
        !    do i = field%xmin, field%xmax
        !        somme(k) = somme(k) + field%values(i,j,k)
        !    end do
        !end do
    end do
    somme = somme / ( nxLoc * nyLoc )
    avg2D = DoGlobalVectorSum(spec_rank, somme)
    avg2D = avg2D / cpuY

    ! Free memory
    deallocate(somme)

    success = .true.

end function computeFieldAvg2D


!---------------------------------------------------------------------------
!> @details
!> Show average of a scalar fields from its spectral values
!! @author = Jean-Baptiste Lagaert, LEGI
!!    @param[in]    field   = scalar field in physical space
!!    @param[in]    spec_rank   = mpi rank inside the spectral communicator
!!    @param[in]    nbcpus   = number of processing unit
!---------------------------------------------------------------------------
subroutine showFieldAvg_spec(field,WN,tag,spec_rank)

     ! Input/Output
     type(COMPLEX_DATA_LAYOUT), intent(in)  :: field
     type(waveNumbers), intent(in)          :: WN
     integer, intent(in)                    :: spec_rank
     character(len=*)                       :: tag

     integer                            :: i, j, k

     if (spec_rank==0) then
         i = field%xmin
         j = field%ymin
         k = field%zmin
         write(*,'(a,f8.4,a,e12.4)') '[info] moyenne (ie k = ', sqrt(WN%kx(i)**2 + WN%ky(j)**2 + WN%kz(k)**2), &
            & ') '//trim(tag), abs(field%values(i,j,k))
     end if



end subroutine showFieldAvg_spec


!---------------------------------------------------------------------------
!> @details
!> Compute variance of a scalar fields.
!! @author = Antoine Vollant, LEGI
!!    @param[in]    field   = scalar field in physical space
!!    @param[in]    spec_rank   = mpi rank inside the spectral communicator
!!    @param[in]    nbcpus   = number of processing unit
!!    @return       var = variance of scalar field
!---------------------------------------------------------------------------
function computeFieldVar(field,spec_rank) result(var)

     ! Input/Output
     type(REAL_DATA_LAYOUT), intent(in) :: field
     integer, intent(in)                :: spec_rank
     real(WP)                           :: var

     var = computeFieldCoVar(field,field,spec_rank)

end function computeFieldVar

!---------------------------------------------------------------------------
!> @details
!> Compute covariance of two scalars fields.
!! @author = Antoine Vollant, LEGI
!!    @param[in]    field1   = scalar field in physical space
!!    @param[in]    field2   = scalar field in physical space
!!    @param[in]    nbcpus   = number of processing unit
!!    @param[in]    spec_rank   = mpi rank inside the spectral communicator
!!    @return       covar = variance of scalar field
!---------------------------------------------------------------------------
function computeFieldCoVar(field1,field2,spec_rank) result(covar)

     ! Input/Output
     type(REAL_DATA_LAYOUT), intent(in) :: field1,field2
     integer, intent(in)                :: spec_rank
     real(WP)                           :: covar

     ! Local data
     real(WP):: avgSquarreSca , squarreAvgSca, tmp1, tmp2
     type(REAL_DATA_LAYOUT)             :: temp

     covar = 0.0_WP
     ! Compute
     temp = field1
     temp%values = field1%values*field2%values
     avgSquarreSca = computeFieldAvg( temp , spec_rank )
     tmp1 = computeFieldAvg( field1 , spec_rank )
     tmp2 = computeFieldAvg( field2 , spec_rank )
     squarreAvgSca = tmp1 * tmp2
     covar = avgSquarreSca - squarreAvgSca

     call deleteDataLayout(temp)

end function computeFieldCoVar

!------------------------------------------------------------------------------
!> Compute correlation of two scalars field.
!! @author Antoine Vollant, LEGI
!!    @param[in]    inField1        = scalar field 1
!!    @param[in]    inField2        = scalar field 2
!!    @param[in]    nbcpus          = number of process in the pool
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[inout] Correlation     = correlation
!!    @return       success         = logical equals to true is everything is right.
!! @details
!!        This function compute the divergence of an array field. 
!------------------------------------------------------------------------------
function computeCorr(Correlation,inField1,inField2,spec_rank) result(success)

    ! Input/Output
    integer, intent(in)                       :: spec_rank
    real(WP),intent(inout)                    :: Correlation
    type(REAL_DATA_LAYOUT), intent(in)        :: inField1,inField2
    logical                                   :: success

    !Local datas
    real(WP)                                  :: covariance
    real(WP)                                  :: variance1,variance2

    success = .false.

    if (.not. samelayout(inField1,inField2) ) then
      write(6,'(a,i0,a)')'[ERROR] Fields in computeCorr are not the same!'
      return
    endif

    covariance = computeFieldCoVar(inField1,inField2,spec_rank) 
    variance1  = computeFieldVar(inField1,spec_rank) 
    variance2  = computeFieldVar(inField2,spec_rank)
    Correlation = covariance / sqrt( variance1 * variance2 )
 
    success = .true.

end function computeCorr

!---------------------------------------------------------------------------
!> @details
!> Compute kurtosis of a scalar fields. if kurtosis equal 0, the distribution has a gaussian shape
!! @author = Antoine Vollant, LEGI
!!    @param[in]    field   = scalar field in physical space
!!    @param[inout] kurt is the kurtosis of data field
!!    @param[in]    spec_rank   = mpi rank inside the spectral communicator
!!    @return       success equal to .true. if everything is OK and .false. otherwize
!---------------------------------------------------------------------------
function computeFieldKurt(field,spec_rank) result(kurt)

     ! I/O data
     type(REAL_DATA_LAYOUT), intent(in) :: field
     real(WP)                           :: kurt
     integer, intent(in)                :: spec_rank

     !Local data
     real(WP)                :: average,variance
     type(REAL_DATA_LAYOUT)  :: temp

     kurt = 0.0_WP
     if (.not. copyStructOnly(field,temp) ) then
       if (spec_rank .eq. 0) write(6,'(a)') '[ERROR] in computeFieldSkew, copy structure failed'
       return
     endif
     average = computeFieldAvg(field,spec_rank)
     variance = computeFieldVar(field,spec_rank)
     temp%values = (field%values - average)**4
     kurt = computeFieldAvg(temp,spec_rank)
     kurt = kurt / variance**(2.0_WP) - 3.0_WP
     call deletedatalayout(temp)

end function computeFieldKurt


!---------------------------------------------------------------------------
!> @details
!> Compute skewness of a scalar fields.
!! @author = Antoine Vollant, LEGI
!!    @param[in]    field   = scalar field in physical space
!!    @param[inout] skew is the skewness of data field
!!    @param[in]    spec_rank   = mpi rank inside the spectral communicator
!!    @return       success equal to .true. if everything is OK and .false. otherwize
!---------------------------------------------------------------------------
function computeFieldSkew(field,spec_rank) result(skew)

     ! I/O data
     type(REAL_DATA_LAYOUT), intent(in) :: field
     real(WP)                           :: skew
     integer, intent(in)                :: spec_rank

     !Local data
     real(WP)                :: average,variance
     type(REAL_DATA_LAYOUT)  :: temp

     skew = 0.0_WP
     if (.not. copyStructOnly(field,temp) ) then
       if (spec_rank .eq. 0) write(6,'(a)') '[ERROR] in computeFieldSkew, copy structure failed'
       return
     endif
     average = computeFieldAvg(field,spec_rank)
     variance = computeFieldVar(field,spec_rank)
     temp%values = (field%values - average)**3
     skew = computeFieldAvg(temp,spec_rank)
     skew = skew / variance**(1.5_WP)
     call deletedatalayout(temp)

end function computeFieldSkew


!---------------------------------------------------------------------------
!> @details
!> Compute mixed skewness of velocity and a scalar fields.
!! @author = Jean-Baptiste Lagaert, LEGI
!!    @param[in]    scal      = scalar field in physical space
!!    @param[in]    U         = velocity first component in physical space
!!    @param[inout] mix_skew  = mixed skewness of velocity and data field
!!    @param[in]    spec_rank = mpi rank inside the spectral communicator
!!    @param[in]    sc_deriv  = derivate of scalar field along X, in physical space - OPTIONAL
!!    @return       success equal to .true. if everything is OK and .false. otherwize
!---------------------------------------------------------------------------
function computeMixedSkew(mix_skew,U,scal,spec_rank,sc_deriv) result(success)

  use differential_tools

  ! I/O data
  type(REAL_DATA_LAYOUT), intent(in):: scal, U
  real(WP), intent(out)             :: mix_skew
  integer, intent(in)               :: spec_rank
  logical                           :: success
  type(REAL_DATA_LAYOUT), intent(in), optional:: sc_deriv

  !Local data
  real(WP)                :: dsc_variance,dU_variance
  integer,dimension(1)    :: numComp
  type(REAL_DATA_LAYOUT)  :: temp, temp_sc

  success = .false.

  ! Init datalayout
  if (.not. copyStructOnly(scal,temp) ) then
    if (spec_rank .eq. 0) write(6,'(a)') '[ERROR] in computeMixedSkew, copy structure 1 failed'
    return
  endif
  ! Compute derivative of velocity
  numComp = (/1/)
  if (.not. computeDerivationField(U,(/1/),temp,spec_rank)) then
    if (spec_rank .eq. 0) write(6,'(a)') '[ERROR] in computeMixedSkew, derivation of U failed'
    return
  end if
  ! Compute variance of du/dx
  dU_variance = computeFieldVar(temp,spec_rank)
  ! Compute scalar derivative if needed and after compute avg[(dsc/dx)*{(dU/dx)**2}]
  if (present(sc_deriv)) then
    temp%values = (sc_deriv%values**2)*temp%values
    mix_skew = computeFieldAvg(temp,spec_rank)
    temp%values = (sc_deriv%values**2)
    dsc_variance = computeFieldAvg(temp,spec_rank)
  else
    if (.not. copyStructOnly(scal,temp_sc) ) then
      if (spec_rank .eq. 0) write(6,'(a)') '[ERROR] in computeMixedSkew, copy structure 2 failed'
      return
    endif
    if (.not. computeDerivationField(scal,(/1/),temp_sc,spec_rank)) then
      if (spec_rank .eq. 0) write(6,'(a)') '[ERROR] in computeMixedSkew, derivation of scal failed'
      return
    end if
    temp_sc%values = (temp_sc%values**2)
    dsc_variance = computeFieldAvg(temp_sc,spec_rank)
    temp%values = temp_sc%values*temp%values
    mix_skew = computeFieldAvg(temp,spec_rank)
    call deletedatalayout(temp_sc)
  end if
  ! Compute mixed skewness
  mix_skew = mix_skew / (sqrt(dU_variance)*dsc_variance)
  ! Free memory
  call deletedatalayout(temp)
  success = .true.

end function computeMixedSkew

!---------------------------------------------------------------------------
!> @details
!> Compute ordinary moment of a scalar fields.
!> m_r = Esperance (X^r)
!! @author = Antoine Vollant, LEGI
!!    @param[in]    field   = scalar field in physical space
!!    @param[in]    order   = scalar field in physical space
!!    @param[inout] moment  = the statistical value (output) 
!!    @param[in]    spec_rank   = mpi rank inside the spectral communicator
!!    @return       success equal to .true. if everything is OK and .false. otherwize
!---------------------------------------------------------------------------
function computeOrdinaryMoment(field,order,spec_rank) result (moment)

  ! I/O data
  type(REAL_DATA_LAYOUT), intent(in) :: field
  real(WP)                           :: moment 
  integer, intent(in)                :: order
  integer, intent(in)                :: spec_rank 

  !Local data
  type(REAL_DATA_LAYOUT)  :: temp

  moment = 0.0_WP
  if (.not. copyStructOnly(field,temp) ) then
    if (spec_rank .eq. 0) write(6,'(a)') '[ERROR] in computeOrdinaryMoment, copy structure 1 failed'
    return
  endif

  if (order .lt. 1) return
  temp%values = field%values ** order
  moment = computeFieldAvg(temp,spec_rank)

  call deletedatalayout(temp)

end function computeOrdinaryMoment

!---------------------------------------------------------------------------
!> @details
!> Compute centered moment of a scalar fields.
!> mu_r = Esperance[(X-Esperance(X))^r]
!! @author = Antoine Vollant, LEGI
!!    @param[in]    field   = scalar field in physical space
!!    @param[in]    order   = scalar field in physical space
!!    @param[inout] moment  = the statistical value (output) 
!!    @param[in]    spec_rank   = mpi rank inside the spectral communicator
!!    @return       success equal to .true. if everything is OK and .false. otherwize
!---------------------------------------------------------------------------
function computeCenteredMoment(field,order,spec_rank) result (moment)

  ! I/O data
  type(REAL_DATA_LAYOUT), intent(in) :: field
  real(WP)                           :: moment 
  integer, intent(in)                :: order
  integer, intent(in)                :: spec_rank 

  !Local data
  real(WP)                :: avg
  integer                 :: i,j,k
  type(REAL_DATA_LAYOUT)  :: temp

  moment = 0.0_WP
  if (.not. copyStructOnly(field,temp) ) then
    if (spec_rank .eq. 0) write(6,'(a)') '[ERROR] in computeOrdinaryMoment, copy structure 1 failed'
    return
  endif
  if (order .lt. 1) return
  avg = computeFieldAvg(field,spec_rank)
  do k=field%zmin,field%zmax
    do j=field%ymin,field%ymax
      do i=field%xmin,field%xmax
        temp%values(i,j,k) = (field%values(i,j,k) - avg) ** order
      enddo
    enddo
  enddo
  moment = computeFieldAvg(temp,spec_rank)
  call deletedatalayout(temp)

end function computeCenteredMoment

!---------------------------------------------------------------------------
!> @details
!> Compute centered moment of a scalar fields.
!> Beta_(r-2) = Esperance[(X-mu)^r sigma^-r] avec mu = m_1 et sigma = sqrt(mu_2)
!! @author = Antoine Vollant, LEGI
!!    @param[in]    field   = scalar field in physical space
!!    @param[in]    order   = scalar field in physical space
!!    @param[inout] moment  = the statistical value (output) 
!!    @param[in]    spec_rank   = mpi rank inside the spectral communicator
!!    @return       success equal to .true. if everything is OK and .false. otherwize
!---------------------------------------------------------------------------
function computeReduceCenteredMoment(field,order,spec_rank) result (moment)

  ! I/O data
  type(REAL_DATA_LAYOUT), intent(in) :: field
  real(WP)                           :: moment 
  integer, intent(in)                :: order
  integer, intent(in)                :: spec_rank 

  !Local data
  real(WP)                :: avg,sig
  integer                 :: i,j,k
  type(REAL_DATA_LAYOUT)  :: temp

  moment = 0.0_WP
  if (.not. copyStructOnly(field,temp) ) then
    if (spec_rank .eq. 0) write(6,'(a)') '[ERROR] in computeOrdinaryMoment, copy structure 1 failed'
    return
  endif

  if (order .lt. 0) return
  avg = computeFieldAvg(field,spec_rank)
  sig = computeCenteredMoment(field,2,spec_rank)
  sig = sqrt(sig)
  do k=field%zmin,field%zmax
    do j=field%ymin,field%ymax
      do i=field%xmin,field%xmax
        temp%values(i,j,k) = ((field%values(i,j,k) - avg)/sig) ** (order+2)
      enddo
    enddo
  enddo
  moment = computeFieldAvg(temp,spec_rank)
  call deletedatalayout(temp)

end function computeReduceCenteredMoment

!------------------------------------------------------------------------------
!> Compute MSE
!! @author = Antoine Vollant, LEGI
!!    @param[in]    U           = longitudinal velocity
!!    @param[in]    V           = vertical velocity
!!    @param[in]    W           = spanwise velocity
!!    @param[in]    spec_rank   = mpi rank inside the spectral communicator
!!    @return       mse
!------------------------------------------------------------------------------
function computeMSE(FieldObj,FieldComputed,spec_rank,adim) result(mse)

    use datalayout
    
    !I/O data
    type(REAL_DATA_LAYOUT), intent(in)  :: FieldObj,FieldComputed 
    integer, intent(in)                 :: spec_rank
    logical, intent(in)                 :: adim 
    real(WP)                            :: mse 
 
    !Local data
    type(REAL_DATA_LAYOUT)  :: temp 
    real(WP)                :: var 
    integer                 :: i,j,k
    

    mse = 0.0_WP
    if (.not. samelayout(FieldObj,FieldComputed) ) then
      write(6,'(a)')'[ERROR] Fields in computeMSE are not the same!'
      return
    endif
    if (.not. copyStructOnly(FieldObj,temp)) then 
      write(6,'(a)')'[ERROR] computeMSE , copyStruct : failed!'
    endif
    do k=FieldObj%zmin,FieldObj%zmax
      do j=FieldObj%ymin,FieldObj%ymax
        do i=FieldObj%xmin,FieldObj%xmax
          temp%values(i,j,k) = (FieldObj%values(i,j,k) - FieldComputed%values(i,j,k))**2 
        enddo
      enddo
    enddo
    mse = computeFieldAvg(temp,spec_rank)
!    if (spec_rank .eq. 0 ) then
!      write(6,'(a,a,a,f10.8)')'    [INFO-STAT] mse for ',trim(adjustl(FieldComputed%name)),' is ',mse 
!    endif
    if ( adim ) then
      var = computeFieldVar(FieldObj,spec_rank)
!      if (spec_rank .eq. 0 ) then
!        write(6,'(a,a,a,f10.8)')'   [INFO-STAT] variance for ',trim(adjustl(FieldObj%name)),' is ',var
!      endif
      mse = mse/var
    endif

    call deleteDataLayout(temp)

end function computeMSE

!---------------------------------------------------------------------------
!> @details 
!> Compare two values 
!! @author  Antoine Vollant, LEGI
!! @param[in]    valIn 
!! @param[in]    valToCompare
!! @return       posix equal -1 if valIn<valToCompare,
!!                               0 if valIn=valToCompare and
!!                               1 if valIn>valToCompare
!---------------------------------------------------------------------------
 function equalToVal_realdatalayout(valIn,valToCompare) result(success)

  type(real_data_layout), intent(in)   :: valIn,valToCompare
  logical                              :: success

  !Local data
  integer   :: i,j,k

  success = .false.
  do k = valIn%zmin,valIn%zmax
    do j = valIn%ymin,valIn%ymax
      do i = valIn%xmin,valIn%xmax
        if (.not. equalToVal_wp_sca(valIn%values(i,j,k),valToCompare%values(i,j,k))) return
      enddo
    enddo
  enddo
  success = .true.

 end function

!---------------------------------------------------------------------------
!> @details
!> Compare two values
!! @author  Antoine Vollant, LEGI
!! @param[in]    valIn
!! @param[in]    valToCompare
!! @return       posix equal -1 if valIn<valToCompare,
!!                               0 if valIn=valToCompare and
!!                               1 if valIn>valToCompare
!---------------------------------------------------------------------------
 function equalToVal_wp_Sca(valIn,valToCompare) result(success)

  real(WP), intent(in)   :: valIn,valToCompare
  logical                :: success

  !Local data
  real(WP)               :: eps = 1e-7_WP

  success = .false.
  if ( valIn .lt. valToCompare + eps ) then
    if ( valIn .gt. valToCompare - eps ) then
      success = .true.
    endif
  endif

 end function equalToVal_wp_Sca

!---------------------------------------------------------------------------
!> @details
!> Compare two values
!! @author  Antoine Vollant, LEGI
!! @param[in]    valIn
!! @param[in]    valToCompare
!! @return       posix equal -1 if valIn<valToCompare,
!!                               0 if valIn=valToCompare and
!!                               1 if valIn>valToCompare
!---------------------------------------------------------------------------
 function equalToVal_real_wp_Sca(valIn,valToCompare) result(success)

  real, intent(in)       :: valIn
  real(WP), intent(in)   :: valToCompare
  logical                :: success

  !Local data
  real                   :: eps = 1e-7_WP

  success = .false.
  if ( valIn .lt. valToCompare + eps ) then
    if ( valIn .gt. valToCompare - eps ) then
      success = .true.
    endif
  endif

 end function equalToVal_real_wp_Sca


!---------------------------------------------------------------------------
!> @details
!> Compare two values
!! @author  Antoine Vollant, LEGI
!! @param[in]    valIn
!! @param[in]    valToCompare
!! @return       posix equal -1 if valIn<valToCompare,
!!                               0 if valIn=valToCompare and
!!                               1 if valIn>valToCompare
!---------------------------------------------------------------------------
 function equalToVal_real_Sca(valIn,valToCompare) result(success)

  real, intent(in)   :: valIn,valToCompare
  logical                :: success

  !Local data
  real                   :: eps = 1e-4_WP

  success = .false.
  if ( valIn .lt. valToCompare + eps ) then
    if ( valIn .gt. valToCompare - eps ) then
      success = .true.
    endif
  endif

 end function equalToVal_real_Sca

!------------------------------------------------------------------------------
!> Compute and print PDF of two scalars field.
!> To plot the pdf, use 'set pm3d map' in gnuplot. then use splot 'pdf_field1_field2.table'
!! @author Antoine Vollant, LEGI
!!    @param[in]    ite             = current time iteration
!!    @param[in]    n_iter          = number max of iteration
!!    @param[in]    field1          = field1 to compute PDF
!!    @param[in]    field2          = field2 to compute PDF
!!    @param[in]    dime            = number of samples for field 1 and 2 (array)
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in]    n_iter          = number of max iteration
!!    @return       success         = logical equals to true is everything is right.
!! @details
!!        This function compute the probability density function of a scalar field.
!------------------------------------------------------------------------------
function computeTwoFieldsPDF(ite,field1,field2,dime,spec_rank,n_iter,minDelta,maxDelta) result(success)

    use datalayout
    use fileio
    use toolbox

    ! Input/Output
    TYPE(REAL_DATA_LAYOUT),intent(in)         :: field1
    TYPE(REAL_DATA_LAYOUT),intent(in)         :: field2
    integer,intent(in)                        :: ite
    integer,intent(in)                        :: n_iter
    integer,intent(in),dimension(2)           :: dime
    integer,intent(in)                        :: spec_rank
    real(WP),intent(in),dimension(2),optional :: minDelta,maxDelta
    logical                                   :: success
    ! Others variables
    real(WP),dimension(:,:), allocatable    :: yPdfLoc          ! to save y coordinates of PDF
    real(WP),dimension(:,:), allocatable    :: yPdf             ! to save y coordinates of PDF
    integer                                 :: i, j ,k ,ll1,ll2 ! indices
    real(WP)                                :: cptLoc,cpt       ! count
    integer                                 :: ires
    real(WP)                                :: mini1,maxi1,taux
    real(WP)                                :: mini2,maxi2
    real(WP)                                :: delta1,delta2    ! size of one PDF's sample
    integer                                 :: file_id          ! id for file io
    character(len=str_long)                       :: file_name        ! name of output file

    ! ===== Init =====
    success = .false.
    if (.not.(size(dime,1) .eq. 2) ) return
    if (.not. samelayout(field1,field2)) return
    allocate ( yPdf(dime(1),dime(2)),stat=ires )
    if (ires .ne. 0) return
    allocate ( yPdfLoc(dime(1),dime(2)),stat=ires )
    if (ires .ne. 0) return
    yPdfLoc = 0.0_WP
    ypdf = 0.0_WP
    cptLoc = 0.0_WP

    if (present(minDelta) .and. present(maxDelta) ) then
      if (.not. (size(minDelta,1) .eq. 2) .or. .not. (size(maxDelta,1) .eq. 2)) return
      mini1 = minDelta(1)
      maxi1 = maxDelta(1)
      mini2 = minDelta(2)
      maxi2 = maxDelta(2)
      do k = field1%zmin,field1%zmax
        do j = field1%ymin,field1%ymax
          do i = field1%xmin,field1%xmax
            ll1 = computeIndice(minDelta(1) , maxDelta(1) , field1%values(i,j,k) , dime(1))
            ll2 = computeIndice(minDelta(2) , maxDelta(2) , field2%values(i,j,k) , dime(2))
            if ( (ll1 .ge. 1) .and. (ll1 .le. dime(1)) .and. &
               & (ll2 .ge. 1) .and. (ll2 .le. dime(2)) ) then 
              yPdfLoc(ll1,ll2)=yPdfLoc(ll1,ll2)+1.0_WP
              cptLoc = cptLoc + 1.0_WP
            endif
          enddo
        enddo
      enddo
      do i=1,dime(1)
        yPdf(i,:) = DoGlobalVectorSum(spec_rank,yPdfLoc(i,:))
      enddo
      cpt  = DoGlobalSum(spec_rank,cptLoc)
      if ( spec_rank .eq. 0 )  then
        taux = cpt  / real(field1%nx,WP)
        taux = taux / real(field1%ny,WP)
        taux = taux / real(field1%nz,WP)
        write (6,'(a,1x,a,1x,a,2(1x,f7.3,1x,f7.3,1x,a),1x,f7.3,1x,a)') '    [INFO PDF] Field',trim(adjustl(field1%name)),&
        'between',mini1,mini2,'and',maxi1,maxi2,'represent',taux*100_WP,'% of total distribution'
      endif
      yPdf = yPdf / cpt
    else
      mini1 = computeFieldMin( field1 , spec_rank )
      maxi1 = computeFieldMax( field1 , spec_rank )
      mini2 = computeFieldMin( field2 , spec_rank )
      maxi2 = computeFieldMax( field2 , spec_rank )
      do k = field1%zmin,field1%zmax
        do j = field1%ymin,field1%ymax
          do i = field1%xmin,field1%xmax
            ll1 = computeIndice(mini1, maxi1, field1%values(i,j,k) , dime(1))
            ll2 = computeIndice(mini2, maxi2, field2%values(i,j,k) , dime(2))
            yPdfLoc(ll1,ll2)=yPdfLoc(ll1,ll2)+1.0_WP
          enddo
        enddo
      enddo
      do i=1,dime(1)
        yPdf(i,:) = DoGlobalVectorSum(spec_rank,yPdfLoc(i,:))
      enddo
      yPdf = yPdf/real(field1%nx,WP)
      yPdf = yPdf/real(field1%ny,WP)
      yPdf = yPdf/real(field1%nz,WP)
    endif
    delta1=(maxi1-mini1)/real((dime(1)),WP)
    delta2=(maxi2-mini2)/real((dime(2)),WP)
    ! Write output
    if (spec_rank==0) then
      file_id = iopen()
      if (n_iter .lt. 1000) then
        write(file_name,'(a,a,a,a,a,i3.3,a)') 'pdf_',trim(adjustl(field1%name)),'_',trim(adjustl(field2%name)),'_at_',ite,'.table'
      elseif (n_iter .lt. 1000000) then
        write(file_name,'(a,a,a,a,a,i6.6,a)') 'pdf_',trim(adjustl(field1%name)),'_',trim(adjustl(field2%name)),'_at_',ite,'.table'
      else
        write(file_name,'(a,a,a,a,a,i0,a)') 'pdf_',trim(adjustl(field1%name)),'_',trim(adjustl(field2%name)),'_at_',ite,'.table'
      endif
      file_name = trim(adjustl(file_name))
      open(unit=file_id, file=file_name, form="FORMATTED")
      do i = 1,dime(1) 
        do j = 1,dime(2) 
          write(file_id, '(g15.8,1x,g15.8,1x,g15.8)') mini1+(i-0.5)*delta1, mini2+(j-0.5)*delta2, yPdf(i,j)
        end do
        write(file_id, '(a)') '' 
      end do
      close(file_id)
      file_id = iclose(file_id)
    end if
    ! Free memory for next scalar field
    deallocate(yPdf)
    deallocate(yPdfLoc)

    ! ===== Return sucess =====
    success = .true.
    return

end function computeTwoFieldsPDF


!------------------------------------------------------------------------------
!> Compute and print PDF of a scalars field.
!! @author Antoine Vollant, LEGI
!!    @param[in]    ite             = current time iteration
!!    @param[in]    n_iter          = number max of iteration
!!    @param[in]    field           = field to compute PDF
!!    @param[in]    dime            = number of samples in PDF
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in]    n_iter          = number of max iteration
!!    @return       success         = logical equals to true is everything is right.
!! @details
!!        This function compute the probability density function of a scalar field.
!------------------------------------------------------------------------------
function computeFieldPDF(ite,field,dime,spec_rank,n_iter,minDelta,maxDelta) result(success)

    use datalayout
    use fileio
    use toolbox

    ! Input/Output
    TYPE(REAL_DATA_LAYOUT), intent(in) :: field
    integer, intent(in)          :: ite
    integer, intent(in)          :: n_iter
    integer, intent(in)          :: dime
    integer, intent(in)          :: spec_rank
    real(WP),intent(in),optional :: minDelta,maxDelta
    logical             :: success
    ! Others variables
    real(WP), dimension(:), allocatable     :: yPdfLoc  ! to save y coordinates of PDF
    real(WP), dimension(:), allocatable     :: yPdf     ! to save y coordinates of PDF
    real(WP), dimension(:), allocatable     :: xPdf     ! to save x coordinates of PDF
    integer                                 :: i, j ,k ,ll ! indices
    real(WP)                                :: cptLoc,cpt ! count
    integer                                 :: ires
    real(WP)                                :: mini,maxi,taux
    real(WP)                                :: delta     ! size of one PDF's sample
    integer                                 :: file_id   ! id for file io
    character(len=50)                       :: file_name ! name of output file

    ! ===== Init =====
    success = .false.
    allocate ( yPdf(dime),stat=ires )
    if (ires .ne. 0) return
    allocate ( yPdfLoc(dime),stat=ires )
    if (ires .ne. 0) return
    allocate ( xPdf(dime),stat=ires )
    if (ires .ne. 0) return
    yPdfLoc = 0.0_WP
    ypdf = 0.0_WP
    cptLoc = 0.0_WP

    if (present(minDelta) .and. present(maxDelta) ) then
      mini = minDelta
      maxi = maxDelta
      do k = field%zmin,field%zmax
        do j = field%ymin,field%ymax
          do i = field%xmin,field%xmax
            ll = computeIndice(mini , maxi , field%values(i,j,k) , dime)
            if ( (ll .ge. 1) .and. (ll .le. dime) ) then 
              yPdfLoc(ll)=yPdfLoc(ll)+1.0_WP
              cptLoc = cptLoc + 1.0_WP
            endif
          enddo
        enddo
      enddo
      yPdf = DoGlobalVectorSum(spec_rank,yPdfLoc)
      cpt  = DoGlobalSum(spec_rank,cptLoc)
      if ( spec_rank .eq. 0 )  then
        taux = cpt  / real(field%nx,WP)
        taux = taux / real(field%ny,WP)
        taux = taux / real(field%nz,WP)
        write (6,'(a,1x,a,1x,a,2(1x,f7.3,1x,a),1x,f7.3,1x,a)') '    [INFO PDF] Field',trim(adjustl(field%name)),&
        'between',mini,'and',maxi,'represent',taux*100_WP,'% of total distribution'
      endif
      yPdf = yPdf / cpt
    else
      mini = computeFieldMin( field , spec_rank )
      maxi = computeFieldMax( field , spec_rank )
      do k = field%zmin,field%zmax
        do j = field%ymin,field%ymax
          do i = field%xmin,field%xmax
            ll = computeIndice(mini , maxi , field%values(i,j,k) , dime)
            yPdfLoc(ll)=yPdfLoc(ll)+1.0_WP
          enddo
        enddo
      enddo
      yPdf = DoGlobalVectorSum(spec_rank,yPdfLoc)
      yPdf = yPdf/real(field%nx,WP)
      yPdf = yPdf/real(field%ny,WP)
      yPdf = yPdf/real(field%nz,WP)
    endif
    delta=(maxi-mini)/real((dime),WP)
!    yPdf = yPdf/delta
    do ll = 1,dime
      xPdf(ll)=mini+(ll-0.5)*delta
    enddo
    ! Write output
    if (spec_rank==0) then
      file_id = iopen()
      if (n_iter .lt. 1000) then
        write(file_name,'(a,a,a,i3.3,a)') 'pdf_',trim(adjustl(field%name)),'_at_',ite,'.table'
      elseif (n_iter .lt. 1000000) then
        write(file_name,'(a,a,a,i6.6,a)') 'pdf_',trim(adjustl(field%name)),'_at_',ite,'.table'
      else
        write(file_name,'(a,a,a,i0,a)') 'pdf_',trim(adjustl(field%name)),'_at_',ite,'.table'
      endif
      file_name = trim(adjustl(file_name))
      open(unit=file_id, file=file_name, form="FORMATTED")
      do i = 1,dime 
        write(file_id, '(g15.8,1x,g15.8)') xPdf(i), yPdf(i)
      end do
      close(file_id)
      file_id = iclose(file_id)
    end if
    ! Free memory for next scalar field
    deallocate(yPdf)
    deallocate(xPdf)
    deallocate(yPdfLoc)

    ! ===== Return sucess =====
    success = .true.
    return

end function computeFieldPDF

!------------------------------------------------------------------------------
!> Compute and save spectrum of a given field (scalars, not kinetic energy
!! spectrum).
!! @author Jean-Baptiste Lagaert
!!    @param[in]    field           = field from wich spectrum has to been computed.
!!    @param[in]    ite             = time iteration number
!!    @param[in]    spec_rank   = mpi rank inside spectral communicator
!!    @param[in]    nbcpus          = number of processes
!!    @return       success         = logical equals to true is everything is right.
!! @details
!!        This function compute the spectrum of a given field. It performs
!!    a FFT in order to obtain it in the Fourrier space (wich is required
!!    to computing the spectrum). It could be perform "out a simulation"
!!    (ie from a field write in an output).
!------------------------------------------------------------------------------
function compute_spectrum_sca(field, ite, spec_rank) result(success)

    use datalayout
    use fileio
    use wavenumber_tools
    use transforms_tools

    ! Input/Output
    type(REAL_DATA_LAYOUT), intent(in)      :: field
    integer, intent(in)                    :: ite
    integer, intent(in)                     :: spec_rank
    logical                                 :: success
    ! Others variables
    type(COMPLEX_DATA_LAYOUT)               :: spec_field       ! to store spectral values of "field"
    type(WaveNumbers)                       :: WN               ! associated wave numbers
    real(WP), dimension(:), allocatable     :: WN_norm          ! to know the different wave number modules
    real(WP), dimension(:), allocatable     :: spectrum         ! to save spectrum
    integer                                 :: size_spec        ! spectrum size
    integer                                 :: i, j ,k, ik      ! indices
    real(WP)                                :: kk               ! norm of wave number
    real(WP)                                :: dk               ! space step in Fourrier space
    real(WP)                                :: kc, eps          ! to filter spectrum
    integer                                 :: file_id          ! id for file io
    character(len=256)                       :: file_name        ! name of output file
    integer                                 :: ires             ! error code for allocation
    integer                                 :: nbcpus

    ! Init
    success = .false.

    ! ===== Check if spectral space step is the same above all direction =====
    if((field%Lx/=field%Ly).or.(field%Lx/=field%Lz)) then
        write(*,'(a)') '[warning] Spectrum not yet implemented if not(Lx=Ly=Lz)'
        return
    else
        dk=2.0_WP*acos(-1.0_WP)/field%Lx
    end if
    kc = 2.0_WP*acos(-1.0_WP)*(field%nx-1)/field%Lx
    !kc = acos(-1.0_WP)*(field%nx-1)/field%Lx
    eps = dk/2.0
    nbcpus=getnbcpus()

    ! ===== Go in spectral space =====
    ! Allocate storage in fourrier space
    if (.not. initDataLayout("scalar_spectral", spec_field,(field%nx/2)+1, &
            & field%ny,field%nz,field%Lx,field%Ly,field%Lz,nbcpus,spec_rank,alongZ)) then
        write(6,'(a,i0,a,a)')'[ERROR] initSolver on process ',spec_rank,&
            & ': not enought memory to transform field in fourier space in order ', &
            & 'to compute spectrum.'
        return
    endif
    ! Perform FFT
    call ftran(field,spec_field,success)
    if (.not. success) then
        write(6,'(a)')'[ERROR] spectrum computation: cannot perform FFT'
        return
    end if

    ! ===== Create WN =====
    if (.not. initWN(WN,spec_rank,field%nx,field%ny,field%nz)) then
      write (6,'(a)') '[ERROR] failed in init WN'
      return
    endif
    if (.not. computeWN(WN,spec_rank,field%Lx,field%Ly,field%Lz,field%nx,field%ny,field%nz)) then
      write (6,'(a)') '[ERROR] failed in init WN'
      return
    endif

    ! ===== Compute spectrum =====
    size_spec = size(WN%kx)
    allocate(spectrum(size_spec),stat=ires)
    if (ires/=0) then
      write(6,'(a,i0,a)')'[ERROR] on process ',spec_rank,': not enought memory .'
      return
    endif
    allocate(WN_norm(size_spec),stat=ires)
    if (ires/=0) then
      write(6,'(a,i0,a)')'[ERROR] on process ',spec_rank,': not enought memory .'
      return
    endif
    WN_norm = -1
    spectrum = 0.0_WP
    ! Compute
    do k = spec_field%zmin,spec_field%zmax
        do j =  spec_field%ymin,spec_field%ymax
            do i =  spec_field%xmin,spec_field%xmax
                ! Compute norm of wave number
                kk = WN%kx(i)*WN%kx(i) + &
                   & WN%ky(j)*WN%ky(j) + &
                   & WN%kz(k)*WN%kz(k)
                kk = sqrt(real(kk,WP))
                if ((kk>eps).and.(kk<=kc)) then
                  ! Compute indice
                  ik = nint(kk/dk) + 1
                  ! And then update spectrum
      !            if ( spec_rank .eq. 0 .and. (ik .gt. size_spec)) print *,'Nombre onde',ik 
                  !if ( ik .gt. size_spec ) ik = size_spec
                  if ( ik .le. size_spec ) then 
                  spectrum(ik) = spectrum(ik) + spec_field%values(i,j,k)*conjg(spec_field%values(i,j,k))
                  endif
                endif
            end do
        end do
    end do

    ! ===== Sum values on all processes =====
    spectrum = DoGlobalVectorSum(spec_rank, spectrum)

    ! ===== Write output =====
    if (spec_rank==0) then
        do ik = 1, size_spec
            WN_norm(ik) = dk*real(ik-1)
        end do
        file_id = iopen()
        write(file_name,'(a,a,i6.6,a)') 'spec_',trim(field%name),ite,'.table'
        open(unit=file_id, file=trim(adjustl(file_name)), form="FORMATTED", ERR=300)
        do ik = 1, size_spec
            if (WN_norm(ik)>=0) write(file_id, '(e14.6,x,e14.6)') WN_norm(ik), spectrum(ik)
        end do
        close(file_id)
        file_id = iclose(file_id)
    end if

    ! ===== Free memory =====
    deallocate(spectrum)
    deallocate(WN_norm)
    if (associated(WN%kx)) deallocate(WN%kx)
    if (associated(WN%ky)) deallocate(WN%ky)
    if (associated(WN%kz)) deallocate(WN%kz)
    WN%kx=>  null()
    WN%ky=>  null()
    WN%kz=>  null()
    call deleteDataLayout(spec_field)

    ! ===== Return sucess =====
    success = .true.
    return

    ! ===== IO error =====
300 write(6,'(a)') 'Unable to open file "'//trim(file_name)//'" for spectrum'
    return

end function compute_spectrum_sca

!------------------------------------------------------------------------------
!> Compute and save spectrum of a given field (scalars, not kinetic energy
!! spectrum).
!! @author Jean-Baptiste Lagaert
!!    @param[in]    field           = field from wich spectrum has to been computed.
!!    @param[in]    ite             = time iteration number
!!    @param[in]    spec_rank   = mpi rank inside spectral communicator
!!    @param[in]    nbcpus          = number of processes
!!    @return       success         = logical equals to true is everything is right.
!! @details
!!        This function compute the spectrum of a given field. It performs
!!    a FFT in order to obtain it in the Fourrier space (wich is required
!!    to computing the spectrum). It could be perform "out a simulation"
!!    (ie from a field write in an output).
!------------------------------------------------------------------------------
function compute_spectrum_vecDotvec_k(Fx,Fy,Fz,Gx,Gy,Gz,ite,spec_rank,WN) result(success)

    use datalayout
    use fileio
    use wavenumber_tools
    use transforms_tools

    ! Input/Output
    type(COMPLEX_DATA_LAYOUT), intent(in)   :: Fx,Fy,Fz 
    type(COMPLEX_DATA_LAYOUT), intent(in)   :: Gx,Gy,Gz 
    integer, intent(in)                     :: ite
    type(WaveNumbers),intent(in)            :: WN               ! associated wave numbers
    integer, intent(in)                     :: spec_rank
    logical                                 :: success
    ! Others variables
    real(WP), dimension(:), allocatable     :: WN_norm          ! to know the different wave number modules
    real(WP), dimension(:), allocatable     :: spectrum         ! to save spectrum
    integer                                 :: size_spec        ! spectrum size
    integer                                 :: i, j ,k, ik      ! indices
    real(WP)                                :: kk               ! norm of wave number
    real(WP)                                :: dk               ! space step in Fourrier space
    real(WP)                                :: kc, eps          ! to filter spectrum
    integer                                 :: file_id          ! id for file io
    character(len=256)                      :: file_name        ! name of output file
    integer                                 :: ires             ! error code for allocation
    integer                                 :: nbcpus

    ! Init
    success = .false.

    ! ===== Check if spectral space step is the same above all direction =====
    if (.not. samelayout(Fx,Fy,Fz,Gx)) return
    if (.not. samelayout(Gx,Gy,Gz)) return
    if((fx%Lx/=fx%Ly).or.(fx%Lx/=fx%Lz)) then
        write(*,'(a)') '[warning] Spectrum not yet implemented if not(Lx=Ly=Lz)'
        return
    else
        dk=2.0_WP*acos(-1.0_WP)/fx%Lx
    end if
    kc = 2.0_WP*acos(-1.0_WP)*(fx%nx-1)/fx%Lx
    eps = dk/2.0
    nbcpus = getnbcpus()

    ! ===== Compute spectrum =====
    size_spec = size(WN%kx)
    allocate(spectrum(size_spec),stat=ires)
    if (ires/=0) then
      write(6,'(a,i0,a)')'[ERROR] on process ',spec_rank,': not enought memory .'
      return
    endif
    allocate(WN_norm(size_spec),stat=ires)
    if (ires/=0) then
      write(6,'(a,i0,a)')'[ERROR] on process ',spec_rank,': not enought memory .'
      return
    endif
    WN_norm = -1
    spectrum = 0.0_WP
    ! Compute
    do k = fx%zmin,fx%zmax
        do j =  fx%ymin,fx%ymax
            do i =  fx%xmin,fx%xmax
                ! Compute norm of wave number
                kk = WN%kx(i)*WN%kx(i) + &
                   & WN%ky(j)*WN%ky(j) + &
                   & WN%kz(k)*WN%kz(k)
                kk = sqrt(real(kk,WP))
                if ((kk>eps).and.(kk<=kc)) then
                  ! Compute indice
                  ik = nint(kk/dk) + 1
                  ! And then update spectrum
                  if ( ik .le. size_spec ) then 
                  spectrum(ik) = spectrum(ik) + 0.5_WP*(fx%values(i,j,k)*conjg(gx%values(i,j,k)) &
                                            & +fy%values(i,j,k)*conjg(gy%values(i,j,k)) &
                                            & +fz%values(i,j,k)*conjg(gz%values(i,j,k)) &
                                            & +gx%values(i,j,k)*conjg(gx%values(i,j,k)) &
                                            & +gy%values(i,j,k)*conjg(gy%values(i,j,k)) &
                                            & +gz%values(i,j,k)*conjg(gz%values(i,j,k)))
                  endif
                endif
            end do
        end do
    end do

    ! ===== Sum values on all processes =====
    spectrum = DoGlobalVectorSum(spec_rank, spectrum)

    ! ===== Write output =====
    if (spec_rank==0) then
        do ik = 1, size_spec
            WN_norm(ik) = dk*real(ik-1)
        end do
        file_id = iopen()
        write(file_name,'(a,a,i6.6,a)') 'spec_',trim(fx%name),ite,'.table'
        open(unit=file_id, file=trim(adjustl(file_name)), form="FORMATTED", ERR=300)
        do ik = 1, size_spec
            if (WN_norm(ik)>=0) write(file_id, '(e14.6,x,e14.6)') WN_norm(ik), spectrum(ik)
        end do
        close(file_id)
        file_id = iclose(file_id)
    end if

    ! ===== Return sucess =====
    success = .true.
    return

    ! ===== IO error =====
300 write(6,'(a)') 'Unable to open file "'//trim(file_name)//'" for spectrum'
    return

end function compute_spectrum_vecDotvec_k




!------------------------------------------------------------------------------
!> Compute and save spectrum of a given field (scalars, not kinetic energy
!! spectrum).
!! @author Jean-Baptiste Lagaert
!!    @param[in]    field           = field from wich spectrum has to been computed.
!!    @param[in]    ite             = time iteration number
!!    @param[in]    spec_rank   = mpi rank inside spectral communicator
!!    @param[in]    nbcpus          = number of processes
!!    @return       success         = logical equals to true is everything is right.
!! @details
!!        This function compute the spectrum of a given field. It performs
!!    a FFT in order to obtain it in the Fourrier space (wich is required
!!    to computing the spectrum). It could be perform "out a simulation"
!!    (ie from a field write in an output).
!------------------------------------------------------------------------------
function compute_spectrum_vec_k(Fx,Fy,Fz,ite,spec_rank,WN) result(success)

    use datalayout
    use fileio
    use wavenumber_tools
    use transforms_tools

    ! Input/Output
    type(COMPLEX_DATA_LAYOUT), intent(in)   :: Fx,Fy,Fz 
    integer, intent(in)                     :: ite
    type(WaveNumbers),intent(in)            :: WN               ! associated wave numbers
    integer, intent(in)                     :: spec_rank
    logical                                 :: success
    ! Others variables
    real(WP), dimension(:), allocatable     :: WN_norm          ! to know the different wave number modules
    real(WP), dimension(:), allocatable     :: spectrum         ! to save spectrum
    integer                                 :: size_spec        ! spectrum size
    integer                                 :: i, j ,k, ik      ! indices
    real(WP)                                :: kk               ! norm of wave number
    real(WP)                                :: dk               ! space step in Fourrier space
    real(WP)                                :: kc, eps          ! to filter spectrum
    integer                                 :: file_id          ! id for file io
    character(len=256)                      :: file_name        ! name of output file
    integer                                 :: ires             ! error code for allocation
    integer                                 :: nbcpus

    ! Init
    success = .false.

    ! ===== Check if spectral space step is the same above all direction =====
    if (.not. samelayout(Fx,Fy,Fz)) return
    if((fx%Lx/=fx%Ly).or.(fx%Lx/=fx%Lz)) then
        write(*,'(a)') '[warning] Spectrum not yet implemented if not(Lx=Ly=Lz)'
        return
    else
        dk=2.0_WP*acos(-1.0_WP)/fx%Lx
    end if
    kc = 2.0_WP*acos(-1.0_WP)*(fx%nx-1)/fx%Lx
    eps = dk/2.0
    nbcpus = getnbcpus()

    ! ===== Compute spectrum =====
    size_spec = size(WN%kx)
    allocate(spectrum(size_spec),stat=ires)
    if (ires/=0) then
      write(6,'(a,i0,a)')'[ERROR] on process ',spec_rank,': not enought memory .'
      return
    endif
    allocate(WN_norm(size_spec),stat=ires)
    if (ires/=0) then
      write(6,'(a,i0,a)')'[ERROR] on process ',spec_rank,': not enought memory .'
      return
    endif
    WN_norm = -1
    spectrum = 0.0_WP
    ! Compute
    do k = fx%zmin,fx%zmax
        do j =  fx%ymin,fx%ymax
            do i =  fx%xmin,fx%xmax
                ! Compute norm of wave number
                kk = WN%kx(i)*WN%kx(i) + &
                   & WN%ky(j)*WN%ky(j) + &
                   & WN%kz(k)*WN%kz(k)
                kk = sqrt(real(kk,WP))
                if ((kk>eps).and.(kk<=kc)) then
                  ! Compute indice
                  ik = nint(kk/dk) + 1
                  ! And then update spectrum
                  if ( ik .le. size_spec ) then 
                  spectrum(ik) = spectrum(ik) +fx%values(i,j,k)*conjg(fx%values(i,j,k)) &
                                            & +fy%values(i,j,k)*conjg(fy%values(i,j,k)) &
                                            & +fz%values(i,j,k)*conjg(fz%values(i,j,k))
                  endif
                endif
            end do
        end do
    end do

    ! ===== Sum values on all processes =====
    spectrum = DoGlobalVectorSum(spec_rank, spectrum)

    ! ===== Write output =====
    if (spec_rank==0) then
        do ik = 1, size_spec
            WN_norm(ik) = dk*real(ik-1)
        end do
        file_id = iopen()
        write(file_name,'(a,a,i6.6,a)') 'spec_',trim(fx%name),ite,'.table'
        open(unit=file_id, file=trim(adjustl(file_name)), form="FORMATTED", ERR=300)
        do ik = 1, size_spec
            if (WN_norm(ik)>=0) write(file_id, '(e14.6,x,e14.6)') WN_norm(ik), spectrum(ik)
        end do
        close(file_id)
        file_id = iclose(file_id)
    end if

    ! ===== Return sucess =====
    success = .true.
    return

    ! ===== IO error =====
300 write(6,'(a)') 'Unable to open file "'//trim(file_name)//'" for spectrum'
    return

end function compute_spectrum_vec_k


!------------------------------------------------------------------------------
!> Compute and save spectrum of a given field (scalars, not kinetic energy
!! spectrum).
!! @author Jean-Baptiste Lagaert
!!    @param[in]    field           = field from wich spectrum has to been computed.
!!    @param[in]    ite             = time iteration number
!!    @param[in]    spec_rank   = mpi rank inside spectral communicator
!!    @param[in]    nbcpus          = number of processes
!!    @return       success         = logical equals to true is everything is right.
!! @details
!!        This function compute the spectrum of a given field. It performs
!!    a FFT in order to obtain it in the Fourrier space (wich is required
!!    to computing the spectrum). It could be perform "out a simulation"
!!    (ie from a field write in an output).
!------------------------------------------------------------------------------
function compute_spectrum_vec(Fx,Fy,Fz,ite,spec_rank) result(success)

    use datalayout
    use fileio
    use wavenumber_tools
    use transforms_tools

    ! Input/Output
    type(REAL_DATA_LAYOUT), intent(in)      :: Fx,Fy,Fz 
    integer, intent(in)                     :: ite
    integer, intent(in)                     :: spec_rank
    logical                                 :: success
    ! Others variables
    type(COMPLEX_DATA_LAYOUT)               :: fxk,fyk,fzk      ! to store spectral values of "field"
    type(WaveNumbers)                       :: WN               ! associated wave numbers
    real(WP), dimension(:), allocatable     :: WN_norm          ! to know the different wave number modules
    real(WP), dimension(:), allocatable     :: spectrum         ! to save spectrum
    integer                                 :: size_spec        ! spectrum size
!    integer                                 :: i, j ,k, ik      ! indices
!    real(WP)                                :: kk               ! norm of wave number
    real(WP)                                :: dk               ! space step in Fourrier space
    real(WP)                                :: kc, eps          ! to filter spectrum
    integer                                 :: file_id          ! id for file io
    character(len=256)                      :: file_name        ! name of output file
    integer                                 :: ires             ! error code for allocation
    integer                                 :: nbcpus

    ! Init
    success = .false.

    ! ===== Check if spectral space step is the same above all direction =====
    if (.not. samelayout(Fx,Fy,Fz)) return
    if((fx%Lx/=fx%Ly).or.(fx%Lx/=fx%Lz)) then
        write(*,'(a)') '[warning] Spectrum not yet implemented if not(Lx=Ly=Lz)'
        return
    else
        dk=2.0_WP*acos(-1.0_WP)/fx%Lx
    end if
    kc = 2.0_WP*acos(-1.0_WP)*(fx%nx-1)/fx%Lx
    eps = dk/2.0
    nbcpus = getnbcpus()

    ! ===== Go in spectral space =====
    ! Allocate storage in fourrier space
    if (.not. initDataLayout("fxk",fxk,(fx%nx/2)+1, &
            & fx%ny,fx%nz,fx%Lx,fx%Ly,fx%Lz,nbcpus,spec_rank,alongZ)) then
        write(6,'(a,i0,a,a)')'[ERROR] initSolver on process ',spec_rank,&
            & ': not enought memory to transform field in fourier space in order ', &
            & 'to compute spectrum.'
        return
    endif
    if (.not. copyStructOnly(fxk,fyk,'fyk')) return
    if (.not. copyStructOnly(fxk,fzk,'fzk')) return
    ! Perform FFT
    call ftran(fx,fxk,success)
    if (.not. success) then
        write(6,'(a)')'[ERROR] spectrum computation: cannot perform FFT'
        return
    end if
    call ftran(fy,fyk,success)
    if (.not. success) then
        write(6,'(a)')'[ERROR] spectrum computation: cannot perform FFT'
        return
    end if
    call ftran(fz,fzk,success)
    if (.not. success) then
        write(6,'(a)')'[ERROR] spectrum computation: cannot perform FFT'
        return
    end if

    ! ===== Create WN =====
    if (.not. initWN(WN,spec_rank,fx%nx,fx%ny,fx%nz)) then
      write (6,'(a)') '[ERROR] failed in init WN'
      return
    endif
    if (.not. computeWN(WN,spec_rank,fx%Lx,fx%Ly,fx%Lz,fx%nx,fx%ny,fx%nz)) then
      write (6,'(a)') '[ERROR] failed in init WN'
      return
    endif


    if(.not. compute_spectrum_vec_k(Fxk,Fyk,Fzk,ite,spec_rank,WN)) return

!    ! ===== Compute spectrum =====
!    size_spec = size(WN%kx)
!    allocate(spectrum(size_spec),stat=ires)
!    if (ires/=0) then
!      write(6,'(a,i0,a)')'[ERROR] on process ',spec_rank,': not enought memory .'
!      return
!    endif
!    allocate(WN_norm(size_spec),stat=ires)
!    if (ires/=0) then
!      write(6,'(a,i0,a)')'[ERROR] on process ',spec_rank,': not enought memory .'
!      return
!    endif
!    WN_norm = -1
!    spectrum = 0.0_WP
!    ! Compute
!    do k = fxk%zmin,fxk%zmax
!        do j =  fxk%ymin,fxk%ymax
!            do i =  fxk%xmin,fxk%xmax
!                ! Compute norm of wave number
!                kk = WN%kx(i)*WN%kx(i) + &
!                   & WN%ky(j)*WN%ky(j) + &
!                   & WN%kz(k)*WN%kz(k)
!                kk = sqrt(real(kk,WP))
!                if ((kk>eps).and.(kk<=kc)) then
!                  ! Compute indice
!                  ik = nint(kk/dk) + 1
!                  ! And then update spectrum
!                  if ( ik .le. size_spec ) then 
!                  spectrum(ik) = spectrum(ik) +fxk%values(i,j,k)*conjg(fxk%values(i,j,k)) &
!                                            & +fyk%values(i,j,k)*conjg(fyk%values(i,j,k)) &
!                                            & +fzk%values(i,j,k)*conjg(fzk%values(i,j,k))
!                  endif
!                endif
!            end do
!        end do
!    end do
!
!    ! ===== Sum values on all processes =====
!    spectrum = DoGlobalVectorSum(spec_rank, spectrum)
!
!    ! ===== Write output =====
!    if (spec_rank==0) then
!        do ik = 1, size_spec
!            WN_norm(ik) = dk*real(ik-1)
!        end do
!        file_id = iopen()
!        write(file_name,'(a,a,i4.4,a)') 'spec_',trim(fx%name),ite,'.table'
!        open(unit=file_id, file=trim(adjustl(file_name)), form="FORMATTED", ERR=300)
!        do ik = 1, size_spec
!            if (WN_norm(ik)>=0) write(file_id, '(e14.6,x,e14.6)') WN_norm(ik), spectrum(ik)
!        end do
!        close(file_id)
!        file_id = iclose(file_id)
!    end if
!
    ! ===== Free memory =====
    deallocate(spectrum)
    deallocate(WN_norm)
    if (associated(WN%kx)) deallocate(WN%kx)
    if (associated(WN%ky)) deallocate(WN%ky)
    if (associated(WN%kz)) deallocate(WN%kz)
    WN%kx=>  null()
    WN%ky=>  null()
    WN%kz=>  null()
    call deleteDataLayout(fxk)
    call deleteDataLayout(fyk)
    call deleteDataLayout(fzk)
!
    ! ===== Return sucess =====
    success = .true.
    return

!    ! ===== IO error =====
!300 write(6,'(a)') 'Unable to open file "'//trim(file_name)//'" for spectrum'
!    return

end function compute_spectrum_vec


!------------------------------------------------------------------------------
!> Compute and print PDF of a Scalars.
!! @author Antoine Vollant, LEGI
!!    @param[in]    ite             = current time iteration
!!    @param[in]    n_iter          = number max of iteration
!!    @param[in]    scalarArray     = array of scalars
!!    @param[in]    nbscalar        = number of scalars in the array
!!    @param[in]    dime            = number of samples in PDF
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in]    n_iter          = number of max iteration
!!    @return       success         = logical equals to true is everything is right.
!! @details
!!        This function compute the probability density function of all scalars fields.
!------------------------------------------------------------------------------
function computeScalarsPDF(ite,scalArray,nbscal,dime,spec_rank,n_iter) result (success)

    use datalayout

    ! Input/Output
    type(REAL_DATA_LAYOUT), intent(in), dimension(:)    :: scalArray
    integer, intent(in) :: ite
    integer, intent(in) :: nbscal
    integer, intent(in) :: dime
    integer, intent(in) :: n_iter
    integer, intent(in) :: spec_rank
    logical             :: success
    ! Others variables
    integer :: sca

    ! Init
    success = .false.

    ! Compute pdf
    do sca=1,nbscal
      if(.not.computeFieldPDF(ite,scalArray(sca),dime,spec_rank,n_iter)) then
        write(*,'(a,i3.3)') '[WARNING] Failed to compute PDF of scalars fiels ',sca
        return
      end if
    enddo
    ! Return sucess
    success = .true.

end function computeScalarsPDF


end module stat_tools
!> @}
