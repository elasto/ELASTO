!> @addtogroup init
!! @{
!------------------------------------------------------------------------------
!
! MODULE: geometric_shape
!
!
! DESCRIPTION: 
!> The module provides tool to init a field to a given geometric shape.
!
!> @details
!!     This module contains the generic procedure to initialize and parametrise the
!! advection solver based on particles method. It also contain the subroutine
!! "advec_step" wich solve the equation for a given time step. It is the only one
!! module which is supposed to be included by a code using this library of
!! particle methods.
!
!> @author
!! Jean-Baptiste Lagaert, LEGI
!
!------------------------------------------------------------------------------
module geometric_shape_tools

    use precision_tools
    use realdatalayout


    implicit none
    private

    type point 
     real(WP) :: X 
     real(WP) :: Y 
     real(WP) :: Z 
    end type point 

    interface computeCoordsSphere
      module procedure computeCoordsSphere_pareto
      module procedure computeCoordsSphere_coef
    end interface computeCoordsSphere

    ! ===== Public procedures =====
    public          :: disc_smooth
    public          :: most_smooth
    public          :: computeCoordsSphere
      

contains

!------------------------------------------------------------------------------
!> Initialise the field to a "smoothed" disc
!!    @param[in,out]    field       = field to initialized
!!    @param[in]        rayon       = size of the disc
!!    @param[in]        center      = center coordinate (optional, if not precise, then the domain center)
!!    @param[in]        smoothness  = stifness coefficient (optional)
!!    @return           res         = success (logical)
!! @details
!! Initialize the field to a regularised heaviside matching to a circle: 
!! if distance between X and the center << rayon then field(X) = 1
!! if distance between X and the center >> rayon then field(X) = 0
!------------------------------------------------------------------------------
function disc_smooth(field, rayon, center, smoothness) result(res)

    ! Input/Output
    type(REAL_DATA_LAYOUT), intent(inout)   :: field
    real(WP)                                :: rayon
    real(WP), dimension(3), optional        :: center
    real(WP), optional                      :: smoothness
    logical                                 :: res
    ! Other local variables
    integer                                 :: i, j, k      ! mesh indices
    real(WP)                                :: rx, ry, rz,rr! rayon
    real(WP)                                :: dx, dy, dz   ! space step
    real(WP), dimension(3)                  :: disc_center
    real(WP)                                :: stifness_coeff

    res = .false.

    if (present(center)) then
        disc_center = center
    else
        disc_center(1) = field%Lx/2.0
        disc_center(2) = field%Ly/2.0
        disc_center(3) = field%Lz/2.0
    end if

    stifness_coeff = 10.0
    if (present(smoothness)) stifness_coeff = smoothness

    dx = field%Lx/field%nx
    dy = field%Lx/field%ny
    dz = field%Lx/field%nz

    do k = field%zmin,field%zmax 
        rz = (k*dz-disc_center(3))**2
        do j = field%ymin,field%ymax
            ry = (j*dy-disc_center(2))**2
            do i = field%xmin,field%xmax
                rx = (i*dx-disc_center(1))**2
                rr = sqrt(rz+ry+rx)
                field%values(i,j,k) = tanh(stifness_coeff*(rayon - rr))
                
            end do
        end do
    end do

    res = .true.


end function



!------------------------------------------------------------------------------
!> Initialise the field to "MoST" letters
!!    @param[in,out]    field       = field to initialized
!!    @param[in]        spec_rank   = rank of processus in cpu pool
!!    @return           res         = success (logical)
!! @details
!! Initialize the field to a regularised heaviside matching to MOST lettres: 
!! If the distance between a letter and a point is greater than epLettre/3 field = 0
!! and 0 otherwize
!------------------------------------------------------------------------------
function most_smooth(field,spec_rank) result(success)

    use stat_tools

    !I/O data
    type(real_data_layout),intent(inout) :: field
    integer,intent(in)                   :: spec_rank
    logical                              :: success

    !Local data
    type(point),dimension(:),allocatable :: M1,M2,O,S,T1,T2
    type(point),dimension(:),allocatable :: coord 
    real(WP)   :: epLettre,finLettre 
    real(WP)   :: lBoite
    real(WP)   :: margeV
    real(WP)   :: dx,dy,dz 
    real(WP)   :: rx,t 
    real(WP)   :: stifness,pi 
    real(WP),dimension(6)   :: ry,rz 
    real(WP)   :: pas,dist,distCurrent 
    integer    :: i,j,k,l,n,nbInc
    integer    :: numLettre 
    logical    :: res,premPassage 


    success = .false.

    if (spec_rank .eq. 0) write (6,'(a)') '[INFO] Initializing scalar in progress...'

    pi = real(acos(-1.0_WP),WP)
    stifness = 1.0_WP 
    lBoite=field%Lx 
    epLettre = 0.05_WP*lBoite
    margeV = 0.2_WP*lBoite
    !nbInc = 100field
    nbInc = field%nx
    allocate (M1(3))
    allocate (M2(3))
    allocate (O(6))
    allocate (S(6))
    allocate (T1(2))
    allocate (T2(2))
    allocate (coord(6))

    numLettre = 0 
    M1(1)%Y=numLettre*lBoite/4.0_WP+epLettre
    M1(1)%Z=margeV
    M1(2)%Y=numLettre*lBoite/4.0_WP+epLettre
    M1(2)%Z=1.2_WP*lBoite
    M1(3)%Y=numLettre*lBoite/4.0_WP+lBoite/8.0_WP
    M1(3)%Z=0.5_WP*lBoite
    M2(1)%Y=numLettre*lBoite/4.0_WP+lBoite/8.0_WP
    M2(1)%Z=0.5_WP*lBoite
    M2(2)%Y=numLettre*lBoite/4.0_WP+lBoite/4.0_WP-epLettre
    M2(2)%Z=1.2_WP*lBoite
    M2(3)%Y=numLettre*lBoite/4.0_WP+lBoite/4.0_WP-epLettre
    M2(3)%Z=margeV
    do i=1,3
      M1(i)%X = 0.5_WP*lBoite
      M2(i)%X = 0.5_WP*lBoite
    enddo
    numLettre = 1; 
    O(1)%Y=numLettre*lBoite/4.0_WP+lBoite/8.0_WP
    O(1)%Z=0.5_WP*lBoite
    O(2)%Y=numLettre*lBoite/4.0_WP+epLettre
    O(2)%Z=0.5_WP*lBoite
    O(3)%Y=numLettre*lBoite/4.0_WP+epLettre
    !O(3)%Z=margeV
    O(3)%Z=0.0_WP
    O(4)%Y=numLettre*lBoite/4.0_WP+lBoite/4.0_WP-epLettre
    !O(4)%Z=margeV
    O(4)%Z=0.0_WP
    O(5)%Y=numLettre*lBoite/4.0_WP+lBoite/4.0_WP-epLettre
    O(5)%Z=0.5_WP*lBoite
    O(6)%Y=numLettre*lBoite/4.0_WP+lBoite/8.0_WP
    O(6)%Z=0.5_WP*lBoite
    do i=1,6
      O(i)%X = 0.5_WP*lBoite
    enddo
    numLettre = 2 
    S(1)%Y=numLettre*lBoite/4.0_WP+lBoite/4.0_WP-2*epLettre
    S(1)%Z=lBoite-margeV
    S(2)%Y=numLettre*lBoite/4.0_WP+epLettre
    S(2)%Z=lBoite-margeV
    S(3)%Y=numLettre*lBoite/4.0_WP+epLettre
    S(3)%Z=0.5_WP*lBoite
    S(4)%Y=numLettre*lBoite/4.0_WP+lBoite/4.0_WP-epLettre
    S(4)%Z=0.5_WP*lBoite
    S(5)%Y=numLettre*lBoite/4.0_WP+lBoite/4.0_WP-epLettre
    S(5)%Z=margeV
    S(6)%Y=numLettre*lBoite/4.0_WP+epLettre
    S(6)%Z=margeV
    do i=1,6
      S(i)%X = 0.5_WP*lBoite
    enddo
    numLettre = 3 
    T1(1)%Y=numLettre*lBoite/4.0_WP+epLettre
    T1(1)%Z=lBoite-margeV
    T1(2)%Y=numLettre*lBoite/4.0_WP+lBoite/4.0_WP-epLettre
    T1(2)%Z=lBoite-margeV
    T2(1)%Y=numLettre*lBoite/4.0_WP+lBoite/8.0_WP
    T2(1)%Z=lBoite-margeV
    T2(2)%Y=numLettre*lBoite/4.0_WP+lBoite/8.0_WP
    T2(2)%Z=margeV
    do i=1,2
      T1(i)%X = 0.5_WP*lBoite
      T2(i)%X = 0.5_WP*lBoite
    enddo

    if (spec_rank .eq. 0) then
    endif

    dx = field%Lx/field%nx
    dy = field%Ly/field%ny
    dz = field%Lz/field%nz
 
    pas = 1.0_WP/real(nbInc,WP)
    premPassage = .true.

    !R)echerche des points les plus proches
    do k = field%zmin,field%zmax
      do j = field%ymin,field%ymax
        do i = field%xmin,field%xmax
          rx = (i*dx-0.5_WP*field%Lx)**2
          dist = sqrt( field%Lx**2 + field%Ly**2 + field%Lz**2 )
          do l=0,nbInc
            t = real(l,WP)*pas
            res = computeBezierCurve( M1 , t ,3, coord(1) )
            res = computeBezierCurve( M2 , t ,3, coord(6) )
            res = computeBezierCurve( O , t ,6, coord(2) )
            res = computeBezierCurve( S , t ,6, coord(3) )
            res = computeBezierCurve( T1 , t ,2, coord(4) )
            res = computeBezierCurve( T2 , t ,2, coord(5) )
            if (.not. res ) return
            do n=1,6
              rz(n) = (real(k,WP)*dz - coord(n)%Z)**2
              ry(n) = (real(j,WP)*dy - coord(n)%Y)**2
              distCurrent = sqrt(rx + ry(n) + rz(n) )
              if ( distCurrent .lt. dist ) dist = distCurrent 
            enddo 
          enddo
          premPassage = .false.
          finLettre = epLettre/0.90_WP
          if ( dist .le. finLettre ) then
            field%values(i,j,k) = 0.5_WP*(1.0_WP+tanh(stifness*tan(pi*(0.5_WP-dist/finLettre))))
          else
            field%values(i,j,k) = 0.0_WP
          endif
        enddo
      enddo
    enddo

    !To set average of scalar init to 0.5
    field%values = 0.5_WP/computeFieldAvg(field,spec_rank) * field%values

    deallocate (M1)
    deallocate (M2)
    deallocate (O)
    deallocate (S)
    deallocate (T1)
    deallocate (T2)
    deallocate (coord)
    success = .true.

end function most_smooth


!------------------------------------------------------------------------------
!> computes coordinates of o bezier curbe in z fixed plan 
!!    @param[in,out]    coords = field to initialized
!!    @param[in]        points = size of the disc
!!    @param[in]        param = center coordinate (optional, if not precise, then the domain center)
!!    @return           success (logical)
!------------------------------------------------------------------------------
function computeBezierCurve( points , param , nbPoints, coords  ) result(success)
  
    use physical_values_tools , only : computeBernsteinPol

    !I/O data
    type(point),dimension(:),intent(in)      :: points 
    type(point),intent(inout)                :: coords 
    integer,intent(in)                       :: nbPoints
    real(WP),intent(in)                      :: param
    logical                                  :: success

    !Local data
    integer :: i,m

    success = .false.

    m = nbPoints-1 
    coords%X = 0.0_WP
    coords%Y = 0.0_WP
    coords%Z = 0.0_WP
    do i=0,m
      coords%Y = computeBernsteinPol( m , i , param ) * points(i+1)%Y + coords%Y       
      coords%Z = computeBernsteinPol( m , i , param ) * points(i+1)%Z + coords%Z       
    enddo

    success = .true.

end function

!---------------------------------------------------------------------------
!> @details 
!> compute coordinates of all points on the grid which are closer from sphere
!> of radius Rsphere
!> @author = Antoine Vollant, LEGI
!>    @param[in]    field for grids information 
!>    @param[in]    me is the number of processus in spectral communicator
!>    @param[in]    Rsphere is the radius of sphere
!>    @param[in]    coordsPoints is an array which contains all cartesians coordinates
!>                  of points to look at
!>    @param[in]    countPts is the number of points found  
!>    @return       success logical .true. or .false 
!---------------------------------------------------------------------------
function computeCoordsSphere_pareto(field,Rsphere,coordsPoints,countPts) result(success)

    use parallel_tools

    !I/O data
    type(real_data_layout),intent(in)               :: field 
    integer,dimension(:,:),allocatable,intent(inout)    :: coordsPoints 
    real(WP),intent(in)                             :: Rsphere
    integer,intent(inout)                           :: countPts
    logical                                         :: success

    !Local data
    integer                                         :: i,j,k,m,n
    integer                                         :: ires,sizeCoord
    integer                                         :: dimax,djmax,dkmax
    real(WP)                                        :: dx,dy,dz
    real(WP)                                        :: Rcurrent,delta
    integer,dimension(:,:),allocatable              :: temp
    !character(len=str_long)                         :: nameite,filename
    !integer                                         :: file_id
    integer,dimension(:),allocatable                :: pareto
    real(WP)                                        :: xProj,yProj,zProj,coefK
    real(WP)                                        :: distCur 

    success = .false.
    countPts = 0
    allocate(coordsPoints(3,100),stat=ires)
    if ( ires .ne. 0) then 
      print *,'ERROR allocation'
      return
    endif
    dx = field%Lx/real(field%nx,WP)
    dy = field%Ly/real(field%ny,WP)
    dz = field%Lz/real(field%nz,WP)
    delta = sqrt(3.0_WP)*(dx*dy*dz)**(1.0_WP/3.0_WP)
    dimax=nint(1.15*Rsphere/dx)
    djmax=nint(1.15*Rsphere/dy)
    dkmax=nint(1.15*Rsphere/dz)
    do i=-dimax , dimax
      do j =-djmax , djmax
        do k =-dkmax , dkmax
          Rcurrent =sqrt( (real(i,WP)*dx)**2 + &
                         &(real(j,WP)*dy)**2 + &
                         &(real(k,WP)*dz)**2  )
          if ( (Rcurrent .lt. Rsphere + delta) .and. (Rcurrent .gt. Rsphere - delta) ) then
            countPts = countPts+1
            if ( countPts .gt. size(coordsPoints,2) ) then
              sizeCoord = size(coordsPoints,2)
              allocate(temp(3,sizeCoord),stat=ires)
              if ( ires .ne. 0) then 
                print *,'ERROR allocation'
                return
              endif
              temp = coordsPoints
              deallocate(coordsPoints)
              allocate(coordsPoints(3,sizeCoord+100),stat=ires)
              if ( ires .ne. 0) then
                print *,'ERROR allocation'
                return
              endif
              do m = 1,sizeCoord
                coordsPoints(:,m) = temp(:,m)
              enddo              
              if ( ires .ne. 0) then
                print *,'ERROR allocation'
                return
              endif
              deallocate(temp)
            endif
            coordsPoints(1,countPts)=i
            coordsPoints(2,countPts)=j
            coordsPoints(3,countPts)=k
          endif
        enddo
      enddo
    enddo
    allocate(pareto(countPts),stat=ires)
    if ( ires .ne. 0) then
      print *,'ERROR allocation'
      return
    endif
    pareto=1
    do i = 1 , countPts
      coefk = Rsphere/sqrt((coordsPoints(1,i)*dx)**2+(coordsPoints(2,i)*dy)**2+(coordsPoints(3,i)*dz)**2)
      xProj = coordsPoints(1,i)*dx*coefk
      yProj = coordsPoints(2,i)*dy*coefk
      zProj = coordsPoints(3,i)*dz*coefk
      distCur = (coordsPoints(1,i)*dx*(1.0_WP-coefk)) ** 2.0_WP + &
                (coordsPoints(2,i)*dy*(1.0_WP-coefk)) ** 2.0_WP + &
                (coordsPoints(3,i)*dz*(1.0_WP-coefk)) ** 2.0_WP
      j = 1
      do while ( j.le.countPts .and. pareto(i) .eq. 1 ) 
        if ( j .ne. i ) then
          if ( (coordsPoints(1,j)*dx - xProj)**2 + &
               (coordsPoints(2,j)*dy - yProj)**2 + &
               (coordsPoints(3,j)*dz - zProj)**2 .lt. &
                distCur ) then
                pareto(i) = 0
          endif
        endif
        j=j+1
      enddo
    enddo
    sizeCoord = size(coordsPoints,2)
    allocate(temp(3,sizeCoord),stat=ires)
    if ( ires .ne. 0) then 
      print *,'ERROR allocation'
      return
    endif
    temp = coordsPoints
    deallocate(coordsPoints)
    allocate(coordsPoints(3,sum(pareto)),stat=ires)
    if ( ires .ne. 0) then
      print *,'ERROR allocation'
      return
    endif
    n=1
    do m = 1,countPts
      if ( pareto(m) == 1 ) then
        coordsPoints(:,n) = temp(:,m)
        n = n+1
      endif
    enddo              
    countPts = sum(pareto)
!    if (me .eq. 0 ) then
!      file_id = iopen()
!      write(filename,'(a)') 'deltaCoordSpherePareto.table'
!      open(unit=file_id, file=trim(adjustl(filename)), form="FORMATTED")
!      write(file_id, '(a)')  '#numpoint dX dY dZ X Y Z'
!      do i=1,countPts
!        write (file_id,'(i0,6(1x,g20.8))') &
!            & i,abs(coordsPoints(1,i)*dx*(1.0_WP-Rsphere/sqrt((coordsPoints(1,i)*dx)**2+(coordsPoints(2,i)*dy)**2+(coordsPoints(3,i)*dz)**2)))&
!               ,abs(coordsPoints(2,i)*dy*(1.0_WP-Rsphere/sqrt((coordsPoints(1,i)*dx)**2+(coordsPoints(2,i)*dy)**2+(coordsPoints(3,i)*dz)**2)))&
!               ,abs(coordsPoints(3,i)*dz*(1.0_WP-Rsphere/sqrt((coordsPoints(1,i)*dx)**2+(coordsPoints(2,i)*dy)**2+(coordsPoints(3,i)*dz)**2)))&
!               ,coordsPoints(1,i)*dx,coordsPoints(2,i)*dy,coordsPoints(3,i)*dz
!      enddo
!      close(file_id)
!      file_id = iclose(file_id)
!    endif

    deallocate(pareto)

    success = .true.

end function computeCoordsSphere_pareto 

!---------------------------------------------------------------------------
!> @details 
!> compute the distances between the point P (coordPoint) and all the others
!> and the normp (normL2**ordre) between the vector at coordPoint (VecRef) and  
!> the other vectors (vecField). Note you can compute normp between vectors of dimension 1 to N
!> @author = Antoine Vollant, LEGI
!>    @param[in]    coordPoint = scalar field in physical space
!>    @param[in]    VecRef 
!>    @param[in]    vecField
!>    @param[in]    ordre 
!>    @param[in]    Rsphere
!>    @param[in]    Rboule
!>    @param[in]    me   = mpi rank inside the spectral communicator
!>    @param[inout] normp 
!>    @param[inout] Dist 
!>    @return       success 
!---------------------------------------------------------------------------
function computeCoordsSphere_coef(field,Rsphere,coordsPoints,countPts,coef,RsphereMax) result(success)

    use parallel_tools

    !I/O data
    type(real_data_layout),intent(in)               :: field 
    integer,dimension(:,:),allocatable,intent(inout)    :: coordsPoints 
    real(WP),intent(in)                             :: Rsphere
    real(WP),intent(in)                             :: coef
    integer,intent(inout)                           :: countPts
    real(WP),intent(in),optional                    :: RsphereMax
    logical                                         :: success

    !Local data
    integer                                         :: i,j,k,m
    integer                                         :: ires,sizeCoord
    integer                                         :: dimax,djmax,dkmax
    real(WP)                                        :: dx,dy,dz
    real(WP)                                        :: Rmin,Rmax,Rcurrent,delta
    integer,dimension(:,:),allocatable              :: temp
!    character(len=str_long)                         :: nameite,filename
!    integer                                         :: file_id

    success = .false.
    countPts = 0
    allocate(coordsPoints(3,100),stat=ires)
    if ( ires .ne. 0) then 
      print *,'ERROR allocation'
      return
    endif
    if ( present (RsphereMax) ) then
      Rmin = Rsphere
      Rmax = RsphereMax
    else
      Rmin = Rsphere
      Rmax = Rsphere
    endif
    dx = field%Lx/real(field%nx,WP)
    dy = field%Ly/real(field%ny,WP)
    dz = field%Lz/real(field%nz,WP)
    delta = sqrt(3.0_WP)/2.0_WP*(dx*dy*dz)**(1.0_WP/3.0_WP)
    delta = coef * delta
    dimax=nint(Rsphere/dx)
    djmax=nint(Rsphere/dy)
    dkmax=nint(Rsphere/dz)
    do i=-dimax , dimax
      do j =-djmax , djmax
        do k =-dkmax , dkmax
          Rcurrent =sqrt( (real(i,WP)*dx)**2 + &
                         &(real(j,WP)*dy)**2 + &
                         &(real(k,WP)*dz)**2  )
          if ( (Rcurrent .lt. Rmax + delta) .and. (Rcurrent .gt. Rmin - delta) ) then
            countPts = countPts+1
            if ( countPts .gt. size(coordsPoints,2) ) then
              sizeCoord = size(coordsPoints,2)
              allocate(temp(3,sizeCoord),stat=ires)
              if ( ires .ne. 0) then 
                print *,'ERROR allocation'
                return
              endif
              temp = coordsPoints
              deallocate(coordsPoints)
              allocate(coordsPoints(3,sizeCoord+100),stat=ires)
              if ( ires .ne. 0) then
                print *,'ERROR allocation'
                return
              endif
              do m = 1,sizeCoord
                coordsPoints(:,m) = temp(:,m)
              enddo              
              deallocate(temp)
            endif
            coordsPoints(1,countPts)=i
            coordsPoints(2,countPts)=j
            coordsPoints(3,countPts)=k
          endif
        enddo
      enddo
    enddo
!    if (me .eq. 0 ) then
!      file_id = iopen()
!      write(nameite,'(f6.3)') delta
!      write(filename,'(a)') 'deltaCoordSphere'//trim(adjustl(nameite))//'.table'
!      open(unit=file_id, file=trim(adjustl(filename)), form="FORMATTED")
!      write(file_id, '(a)')  '#numpoint dX dY dZ X Y Z'
!      do i=1,countPts
!        write (file_id,'(i0,6(1x,g20.8))') &
!            & i,abs(coordsPoints(1,i)*dx*(1.0_WP-Rsphere/sqrt((coordsPoints(1,i)*dx)**2+(coordsPoints(2,i)*dy)**2+(coordsPoints(3,i)*dz)**2)))&
!               ,abs(coordsPoints(2,i)*dy*(1.0_WP-Rsphere/sqrt((coordsPoints(1,i)*dx)**2+(coordsPoints(2,i)*dy)**2+(coordsPoints(3,i)*dz)**2)))&
!               ,abs(coordsPoints(3,i)*dz*(1.0_WP-Rsphere/sqrt((coordsPoints(1,i)*dx)**2+(coordsPoints(2,i)*dy)**2+(coordsPoints(3,i)*dz)**2)))&
!               ,coordsPoints(1,i)*dx,coordsPoints(2,i)*dy,coordsPoints(3,i)*dz
!      enddo
!      close(file_id)
!      file_id = iclose(file_id)
!    endif

    success = .true.

end function computeCoordsSphere_coef



end module geometric_shape_tools

!> @}
