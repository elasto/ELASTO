!> @addtogroup post_process
!! @{

!------------------------------------------------------------------------------
!
! MODULE: post_scalar
!
!
!
!> @details
!> @author
!
!------------------------------------------------------------------------------

module post_scalar

    use precision_tools
    use datalayout
    use toolbox
    use stat_tools
    use transforms_tools
    use post_lib
    use scamodels
    use wavenumber_tools
    use conditional_mean_tools
    use parallel_tools
    use subgrid_tools
    use dealias_tools
    use fileio
    use parser_tools
    use differential_tools
    use filtering_tools
    use physical_values_tools

    implicit none
    private

    ! ==== Interface procedures ====
#ifdef POSTCTR
    public :: gradientmodOE
#endif
    public :: computedTidxi
    public :: compute_scalar_basic_stat
    public :: computeScalarMomentsStat 
    public :: apriori_ScalarStat
    public :: compute_scalar_spectrum_rescaled
    public :: apriori_Clark
    public :: apriori_ClarkFabre
    public :: apriori_Gradient
    public :: apriori_Smagorinsky
    public :: apriori_SmagorinskyDyn
    public :: apriori_RGM
    public :: apriori_RGM_analyt
    public :: computeFunctionStructureV1


contains


!---------------------------------------------------------------------------
!> @details 
!> compute the distances between the point P (coordPoint) and all the others
!> and the field (normL2**ordre) between the vector at coordPoint (VecRef) and  
!> the other vectors (vecField). Note you can compute field between vectors of dimension 1 to N
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
function computeIncrementSca(coordPoint,VecRef,vecField,normp,Dist,ordre,me,Rsphere,Rboule) result(success)

    use parallel_tools

    !I/O data
    type(real_data_layout),intent(in)               :: vecField 
    integer,dimension(3),intent(in)                 :: coordPoint 
    real(WP),intent(in)                             :: VecRef
    type(real_data_layout),intent(inout)            :: normp 
    type(real_data_layout),intent(inout)            :: Dist 
    real(WP),intent(in)                             :: ordre
    integer,intent(in)                              :: me
    real(WP),intent(in),optional                    :: Rsphere
    real(WP),intent(in),optional                    :: Rboule
    logical                                         :: success

    !Local data
    integer                                         :: i,j,k,m
    integer                                         :: dimax,djmax,dkmax
!    integer                                         :: countPts
    real(WP)                                        :: Rdisque
    integer                                         :: bufI,bufJ,bufK,dxsave,dysave,dzsave
    real(WP)                                        :: dx,dy,dz,t1,t2,distsave
    integer,dimension(27,3)                         :: signe 

    success = .false.
!call CPU_time(t1)
    if (.not. samelayout(vecField,normp) .or. &
       &.not. samelayout(vecField,Dist)  ) then
      write(6,'(a,i0,a)')'[ERROR] Fields are not the same!'
      return
    endif
    dx = vecField%Lx/real(vecField%nx,WP)
    dy = vecField%Ly/real(vecField%ny,WP)
    dz = vecField%Lz/real(vecField%nz,WP)
!    countPts = 0

    if ( .not. present (Rsphere) .and. .not. present (Rboule) ) then
!      if (me .eq. 0) write (6,'(a)') '[INFO] Calcul de la distance pour tous les points '
      signe = reshape ((/-1,-1,-1, 0,-1,-1, 1,-1,-1, -1,0,-1, 0,0,-1, 1,0,-1, -1,1,-1, 0,1,-1, 1,1,-1,&
                        &-1,-1,0,  0,-1,0,  1,-1,0,  -1,0,0,  0,0,0,  1,0,0,  -1,1,0,  0,1,0,  1,1,0,&
                        &-1,-1,1,  0,-1,1,  1,-1,1,  -1,0,1,  0,0,1,  1,0,1,  -1,1,1,  0,1,1,  1,1,1&
                        &/),shape(signe),order=(/2,1/))
      do k=normp%zmin,normp%zmax
        do j=normp%ymin,normp%ymax
          do i=normp%xmin,normp%xmax
            distSave = 2.0_WP*vecField%Lx 
            do m=1,size(signe,1)
              bufI = signe(m,1) * vecField%nx + coordPoint(1)
              bufJ = signe(m,2) * vecField%ny + coordPoint(2)
              bufK = signe(m,3) * vecField%nz + coordPoint(3)
              Dist%values(i,j,k) = sqrt( (real(i - bufI,WP)*dx)**2 + &
                                       & (real(j - bufJ,WP)*dy)**2 + &
                                       & (real(k - bufK,WP)*dz)**2   )
              if ( Dist%values(i,j,k) .lt. distSave ) distSave = Dist%values(i,j,k)
            enddo
            Dist%values(i,j,k) = distSave
            normp%values(i,j,k) = ( VecRef - vecField%values(i,j,k) ) ** ordre
          enddo
        enddo
      enddo
    elseif ( present (Rsphere) .and. .not. present (Rboule) ) then
!      if (me .eq. 0) write (6,'(a)') '[INFO] Calcul a R fixe'
      normp%values =  0.0_WP 
      Dist%values  = -1.0_WP 
!      dkmax=floor(Rsphere/dz)
      dkmax=nint(Rsphere/dz)
      do k=coordPoint(3)-dkmax,coordPoint(3)+dkmax
        dzsave = abs( coordPoint(3) - k )
        if ( k .gt. normp%nz ) then
          bufK = k - normp%nz
        elseif ( k .lt. 1 ) then
          bufK = normp%nz +k
        else
          bufK = k
        endif
        if ( (bufK .ge. normp%zmin) .and. (bufK .le. normp%zmax)) then
          Rdisque = sqrt(Rsphere**2-(real(dzsave,WP)*dz)**2)
!          djmax   = floor(Rdisque/dy)
          djmax   = nint(Rdisque/dy)
          do j=coordPoint(2)-djmax,coordPoint(2)+djmax
            dysave = abs( coordPoint(2) - j )
            if ( j .gt. normp%nz ) then
              bufJ = j - normp%ny
            elseif ( j .lt. 1 ) then
              bufJ = normp%ny + j
            else
              bufJ = j
            endif
            if ( (bufJ .ge. normp%ymin) .and. (bufJ .le. normp%ymax)) then
              dimax=floor(sqrt(Rdisque**2-(real(dysave,WP)*dy)**2)/dx)
              !dimax=nint(sqrt(Rdisque**2-(real(dysave,WP)*dy)**2)/dx)
              do  m =1,2
                i = coordPoint(1) + (-1.0_WP)**m*dimax
                dxsave = abs (i - coordPoint(1))
                if ( i .gt. normp%nx) then
                  bufI = i - normp%nx
                elseif ( i .lt. 1 ) then
                  bufI = normp%nx + i
                else
                  bufI = i
                endif
                if ((bufI .ge. normp%xmin) .and. (bufI .le. normp%xmax)) then
!                  countPts = countPts+1
                  normp%values(bufI,bufJ,bufK) = ( VecRef - vecField%values(bufI,bufJ,bufK) ) ** ordre
                  Dist%values(bufI,bufJ,bufK) = sqrt( (real(dxsave,WP)*dx)**2 + &
                                                    & (real(dysave,WP)*dy)**2 + &
                                                    & (real(dzsave,WP)*dz)**2   )
                endif
              enddo
            endif
          enddo
        endif
      enddo
!      countPts = DoGlobalSum(me,countPts)
!      if ( me .eq. 0) write(6,'(a,1x,i0)'),'[INFO computeIncrement] countPts = ',countPts
    elseif ( .not. present (Rsphere) .and. present (Rboule) ) then
      !de 0 à R
      if (me .eq. 0) write (6,'(a,1x,f5.2)') '[INFO] Calcul de 0 à',Rboule
      normp%values =  0.0_WP
      Dist%values  = -1.0_WP 
      dkmax=floor(Rboule/dz)
      do k=coordPoint(3)-dkmax,coordPoint(3)+dkmax
        dzsave = abs( coordPoint(3) - k )
        bufK = k
        if ( k .gt. normp%nz ) bufK = k - normp%nz
        if ( k .lt. 1 )        bufK = normp%nz +k
        if ( (bufK .ge. normp%zmin) .and. (bufK .le. normp%zmax)) then
          Rdisque = sqrt(Rboule**2-(real(dzsave,WP)*dz)**2)
          djmax   = floor(Rdisque/dy)
          do j=coordPoint(2)-djmax,coordPoint(2)+djmax
            dysave = abs( coordPoint(2) - j )
            bufJ = j
            if ( j .gt. normp%nz ) bufJ = j - normp%ny
            if ( j .lt. 1 )        bufJ = normp%ny + j
            if ( (bufJ .ge. normp%ymin) .and. (bufJ .le. normp%ymax)) then
              dimax=floor(sqrt(Rdisque**2-(real(dysave,WP)*dy)**2)/dx)
              do i=coordPoint(1)-dimax,coordPoint(1)+dimax 
                dxsave = abs (i - coordPoint(1))
                bufI = i
                if ( i .gt. normp%nx) bufI = i - normp%nx
                if ( i .lt. 1 )       bufI = normp%nx + i
                if ((bufI .ge. normp%xmin) .and. (bufI .le. normp%xmax)) then
!                  countPts = countPts+1
                  normp%values(bufI,bufJ,bufK) = ( VecRef - vecField%values(bufI,bufJ,bufK) ) ** ordre
                  Dist%values(bufI,bufJ,bufK) = sqrt( (real(dxsave,WP)*dx)**2 + &
                                                    & (real(dysave,WP)*dy)**2 + &
                                                    & (real(dzsave,WP)*dz)**2   )
                endif
              enddo
            endif
          enddo
        endif
      enddo
!     signe = reshape ((/-1,-1,-1, 0,-1,-1, 1,-1,-1, -1,0,-1, 0,0,-1, 1,0,-1, -1,1,-1, 0,1,-1, 1,1,-1,&
!                       &-1,-1,0,  0,-1,0,  1,-1,0,  -1,0,0,  0,0,0,  1,0,0,  -1,1,0,  0,1,0,  1,1,0,&
!                       &-1,-1,1,  0,-1,1,  1,-1,1,  -1,0,1,  0,0,1,  1,0,1,  -1,1,1,  0,1,1,  1,1,1&
!                       &/),shape(signe),order=(/2,1/))
!     bufI=coordPoint(1)
!     bufJ=coordPoint(2)
!     bufK=coordPoint(3)
!     m = 1
!     !Pour r allant de 0 à taille de boite sur deux
!     do k=normp%zmin,normp%zmax
!       do j=normp%ymin,normp%ymax
!         do i=normp%xmin,normp%xmax
!           n = 1
!           do while ( ((2.0_WP*real(i - bufI,WP)*dx/vecField%Lx)**2 + &
!                    &  (2.0_WP*real(j - bufJ,WP)*dy/vecField%Ly)**2 + &
!                    &  (2.0_WP*real(k - bufK,WP)*dz/vecField%Lz)**2 .gt. 1.0_WP) .and. (n .le. size(signe,1)) )
!             bufI = signe(m,1) * vecField%nx + coordPoint(1)
!             bufJ = signe(m,2) * vecField%ny + coordPoint(2)
!             bufK = signe(m,3) * vecField%nz + coordPoint(3)
!             m = m + 1
!             n = n + 1
!             if (m .eq. size(signe,1)+1) m = 1
!           enddo
!           if  ( (2.0_WP*real(i - bufI,WP)*dx/vecField%Lx)**2 + &
!               & (2.0_WP*real(j - bufJ,WP)*dy/vecField%Ly)**2 + &
!               & (2.0_WP*real(k - bufK,WP)*dz/vecField%Lz)**2 .lt. 1.0_WP + 1e-15 ) then
!             Dist%values(i,j,k) = sqrt( (real(i - bufI,WP)*dx)**2 + &
!                                      & (real(j - bufJ,WP)*dy)**2 + &
!                                      & (real(k - bufK,WP)*dz)**2   )
!             normp%values(bufI,bufJ,bufK) = ( VecRef - vecField%values(bufI,bufJ,bufK) ) ** ordre
!             countPts = countPts + 1
!           else
!             normp%values(i,j,k) =  0.0_WP 
!             Dist%values(i,j,k)  = -1.0_WP 
!           endif
!         enddo
!       enddo
!     enddo
!      countPts = DoGlobalSum(me,countPts)
!      if ( me .eq. 0) write(6,'(a,1x,i0)'),'[INFO computeIncrement] countPts = ',countPts
    else
      if ( me .eq. 0) print *,'[ERROR]'
      return
    endif
!call CPU_time(t2)
!    if (me == 0) print *,'[INFO] tps calcul NormpDist',t2-t1    

    success = .true.

end function computeIncrementSca



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
function computeIncrementVec(coordPoint,VecRef,vecField,normp,Dist,ordre,me,Rsphere,Rboule) result(success)

    use parallel_tools

    !I/O data
    type(real_data_layout),dimension(:),intent(in)  :: vecField 
    integer,dimension(3),intent(in)                 :: coordPoint 
    real(WP),dimension(:),intent(in)                :: VecRef
    type(real_data_layout),intent(inout)            :: normp
    type(real_data_layout),intent(inout)            :: Dist 
    real(WP),intent(in)                             :: ordre
    integer,intent(in)                              :: me
    real(WP),intent(in),optional                    :: Rsphere
    real(WP),intent(in),optional                    :: Rboule
    logical                                         :: success

    !Local data
    integer                                         :: i,j,k,l,m,sizeOfVecRef 
    integer                                         :: dimax,djmax,dkmax
!    integer                                         :: countPts
    real(WP)                                        :: Rdisque
    integer                                         :: bufI,bufJ,bufK,dxsave,dysave,dzsave
    real(WP)                                        :: dx,dy,dz,distsave
    integer,dimension(27,3)                         :: signe 

    success = .false.
    sizeOfVecRef = size(VecRef,1)
    if ( sizeOfVecRef .le. 0) then
      write (6,'(a)') '[ERROR] sizeOfVecRef <= 0'
      return 
    endif
    if (.not. samelayout(vecField(1),normp) .or. &
       &.not. samelayout(vecField(1),Dist)  ) then
      write(6,'(a,i0,a)')'[ERROR] Fields are not the same!'
      return
    endif
    dx = vecField(1)%Lx/real(vecField(1)%nx,WP)
    dy = vecField(1)%Ly/real(vecField(1)%ny,WP)
    dz = vecField(1)%Lz/real(vecField(1)%nz,WP)
!    countPts = 0

    if ( .not. present (Rsphere) .and. .not. present (Rboule) ) then
!      if (me .eq. 0) write (6,'(a)') '[INFO] Calcul de la distance pour tous les points '
      signe = reshape ((/-1,-1,-1, 0,-1,-1, 1,-1,-1, -1,0,-1, 0,0,-1, 1,0,-1, -1,1,-1, 0,1,-1, 1,1,-1,&
                        &-1,-1,0,  0,-1,0,  1,-1,0,  -1,0,0,  0,0,0,  1,0,0,  -1,1,0,  0,1,0,  1,1,0,&
                        &-1,-1,1,  0,-1,1,  1,-1,1,  -1,0,1,  0,0,1,  1,0,1,  -1,1,1,  0,1,1,  1,1,1&
                        &/),shape(signe),order=(/2,1/))
      do k=normp%zmin,normp%zmax
        do j=normp%ymin,normp%ymax
          do i=normp%xmin,normp%xmax
            distSave = 2.0_WP*vecField(1)%Lx 
            do m=1,size(signe,1)
              bufI = signe(m,1) * vecField(1)%nx + coordPoint(1)
              bufJ = signe(m,2) * vecField(1)%ny + coordPoint(2)
              bufK = signe(m,3) * vecField(1)%nz + coordPoint(3)
              Dist%values(i,j,k) = sqrt( (real(i - bufI,WP)*dx)**2 + &
                                       & (real(j - bufJ,WP)*dy)**2 + &
                                       & (real(k - bufK,WP)*dz)**2   )
              if ( Dist%values(i,j,k) .lt. distSave ) distSave = Dist%values(i,j,k)
            enddo
            Dist%values(i,j,k) = distSave
            normp%values(i,j,k) = 0.0_WP
            do l=1,sizeOfVecRef
              normp%values(i,j,k) = normp%values(i,j,k) + ( VecRef(l) - vecField(l)%values(i,j,k) ) **2
            enddo
            normp%values(i,j,k) = sqrt(normp%values(i,j,k)) ** ordre
          enddo
        enddo
      enddo
    elseif ( present (Rsphere) .and. .not. present (Rboule) ) then
!      if (me .eq. 0) write (6,'(a)') '[INFO] Calcul a R fixe'
      normp%values =  0.0_WP 
      Dist%values  = -1.0_WP 
!      dkmax=floor(Rsphere/dz)
      dkmax=nint(Rsphere/dz)
      do k=coordPoint(3)-dkmax,coordPoint(3)+dkmax
        dzsave = abs( coordPoint(3) - k )
        if ( k .gt. normp%nz ) then
          bufK = k - normp%nz
        elseif ( k .lt. 1 ) then
          bufK = normp%nz +k
        else
          bufK = k
        endif
        if ( (bufK .ge. normp%zmin) .and. (bufK .le. normp%zmax)) then
          Rdisque = sqrt(Rsphere**2-(real(dzsave,WP)*dz)**2)
!          djmax   = floor(Rdisque/dy)
          djmax   = nint(Rdisque/dy)
          do j=coordPoint(2)-djmax,coordPoint(2)+djmax
            dysave = abs( coordPoint(2) - j )
            if ( j .gt. normp%nz ) then
              bufJ = j - normp%ny
            elseif ( j .lt. 1 ) then
              bufJ = normp%ny + j
            else
              bufJ = j
            endif
            if ( (bufJ .ge. normp%ymin) .and. (bufJ .le. normp%ymax)) then
              dimax=floor(sqrt(Rdisque**2-(real(dysave,WP)*dy)**2)/dx)
              !dimax=nint(sqrt(Rdisque**2-(real(dysave,WP)*dy)**2)/dx)
              do  m =1,2
                i = coordPoint(1) + (-1.0_WP)**m*dimax
                dxsave = abs (i - coordPoint(1))
                if ( i .gt. normp%nx) then
                  bufI = i - normp%nx
                elseif ( i .lt. 1 ) then
                  bufI = normp%nx + i
                else
                  bufI = i
                endif
                if ((bufI .ge. normp%xmin) .and. (bufI .le. normp%xmax)) then
!                  countPts = countPts+1
                  do l=1,sizeOfVecRef
                    normp%values(bufI,bufJ,bufK) = normp%values(bufI,bufJ,bufK) + ( VecRef(l) - vecField(l)%values(bufI,bufJ,bufK) ) **2
                  enddo
                  normp%values(bufI,bufJ,bufK) = sqrt(normp%values(bufI,bufJ,bufK)) ** ordre
                  Dist%values(bufI,bufJ,bufK) = sqrt( (real(dxsave,WP)*dx)**2 + &
                                                    & (real(dysave,WP)*dy)**2 + &
                                                    & (real(dzsave,WP)*dz)**2   )
                endif
              enddo
            endif
          enddo
        endif
      enddo
!      countPts = DoGlobalSum(me,countPts)
!      if ( me .eq. 0) write(6,'(a,1x,i0)'),'[INFO computeIncrement] countPts = ',countPts
    elseif ( .not. present (Rsphere) .and. present (Rboule) ) then
      !de 0 à R
      if (me .eq. 0) write (6,'(a,1x,f5.2)') '[INFO] Calcul de 0 à',Rboule
      normp%values =  0.0_WP
      Dist%values  = -1.0_WP 
      dkmax=floor(Rboule/dz)
      do k=coordPoint(3)-dkmax,coordPoint(3)+dkmax
        dzsave = abs( coordPoint(3) - k )
        bufK = k
        if ( k .gt. normp%nz ) bufK = k - normp%nz
        if ( k .lt. 1 )        bufK = normp%nz +k
        if ( (bufK .ge. normp%zmin) .and. (bufK .le. normp%zmax)) then
          Rdisque = sqrt(Rboule**2-(real(dzsave,WP)*dz)**2)
          djmax   = floor(Rdisque/dy)
          do j=coordPoint(2)-djmax,coordPoint(2)+djmax
            dysave = abs( coordPoint(2) - j )
            bufJ = j
            if ( j .gt. normp%nz ) bufJ = j - normp%ny
            if ( j .lt. 1 )        bufJ = normp%ny + j
            if ( (bufJ .ge. normp%ymin) .and. (bufJ .le. normp%ymax)) then
              dimax=floor(sqrt(Rdisque**2-(real(dysave,WP)*dy)**2)/dx)
              do i=coordPoint(1)-dimax,coordPoint(1)+dimax 
                dxsave = abs (i - coordPoint(1))
                bufI = i
                if ( i .gt. normp%nx) bufI = i - normp%nx
                if ( i .lt. 1 )       bufI = normp%nx + i
                if ((bufI .ge. normp%xmin) .and. (bufI .le. normp%xmax)) then
!                  countPts = countPts+1
                  do l=1,sizeOfVecRef
                    normp%values(bufI,bufJ,bufK) = normp%values(bufI,bufJ,bufK) + ( VecRef(l) - vecField(l)%values(bufI,bufJ,bufK) ) **2
                  enddo
                  normp%values(bufI,bufJ,bufK) = sqrt(normp%values(bufI,bufJ,bufK)) ** ordre
                  Dist%values(bufI,bufJ,bufK) = sqrt( (real(dxsave,WP)*dx)**2 + &
                                                    & (real(dysave,WP)*dy)**2 + &
                                                    & (real(dzsave,WP)*dz)**2   )
                endif
              enddo
            endif
          enddo
        endif
      enddo
!     signe = reshape ((/-1,-1,-1, 0,-1,-1, 1,-1,-1, -1,0,-1, 0,0,-1, 1,0,-1, -1,1,-1, 0,1,-1, 1,1,-1,&
!                       &-1,-1,0,  0,-1,0,  1,-1,0,  -1,0,0,  0,0,0,  1,0,0,  -1,1,0,  0,1,0,  1,1,0,&
!                       &-1,-1,1,  0,-1,1,  1,-1,1,  -1,0,1,  0,0,1,  1,0,1,  -1,1,1,  0,1,1,  1,1,1&
!                       &/),shape(signe),order=(/2,1/))
!     bufI=coordPoint(1)
!     bufJ=coordPoint(2)
!     bufK=coordPoint(3)
!     m = 1
!     !Pour r allant de 0 à taille de boite sur deux
!     do k=normp%zmin,normp%zmax
!       do j=normp%ymin,normp%ymax
!         do i=normp%xmin,normp%xmax
!           n = 1
!           do while ( ((2.0_WP*real(i - bufI,WP)*dx/vecField(1)%Lx)**2 + &
!                    &  (2.0_WP*real(j - bufJ,WP)*dy/vecField(1)%Ly)**2 + &
!                    &  (2.0_WP*real(k - bufK,WP)*dz/vecField(1)%Lz)**2 .gt. 1.0_WP) .and. (n .le. size(signe,1)) )
!             bufI = signe(m,1) * vecField(1)%nx + coordPoint(1)
!             bufJ = signe(m,2) * vecField(1)%ny + coordPoint(2)
!             bufK = signe(m,3) * vecField(1)%nz + coordPoint(3)
!             m = m + 1
!             n = n + 1
!             if (m .eq. size(signe,1)+1) m = 1
!           enddo
!           if  ( (2.0_WP*real(i - bufI,WP)*dx/vecField(1)%Lx)**2 + &
!               & (2.0_WP*real(j - bufJ,WP)*dy/vecField(1)%Ly)**2 + &
!               & (2.0_WP*real(k - bufK,WP)*dz/vecField(1)%Lz)**2 .lt. 1.0_WP + 1e-15 ) then
!             Dist%values(i,j,k) = sqrt( (real(i - bufI,WP)*dx)**2 + &
!                                      & (real(j - bufJ,WP)*dy)**2 + &
!                                      & (real(k - bufK,WP)*dz)**2   )
!             normp%values(bufI,bufJ,bufK) = 0.0_WP
!             do l=1,sizeOfVecRef
!               normp%values(bufI,bufJ,bufK) = normp%values(bufI,bufJ,bufK) + ( VecRef(l) - vecField(l)%values(bufI,bufJ,bufK) ) **2
!             enddo
!             normp%values(bufI,bufJ,bufK) = sqrt(normp%values(bufI,bufJ,bufK)) ** ordre
!             countPts = countPts + 1
!           else
!             normp%values(i,j,k) =  0.0_WP 
!             Dist%values(i,j,k)  = -1.0_WP 
!           endif
!         enddo
!       enddo
!     enddo
!      countPts = DoGlobalSum(me,countPts)
!      if ( me .eq. 0) write(6,'(a,1x,i0)'),'[INFO computeIncrement] countPts = ',countPts
    else
      if ( me .eq. 0) print *,'[ERROR]'
      return
    endif

    success = .true.

end function computeIncrementVec

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
!>    @param[in]    me   = mpi rank inside the spectral communicator
!>    @param[inout] normp 
!>    @param[inout] Dist 
!>    @return       success 
!---------------------------------------------------------------------------
function computeFunctionStructureVec(field1,field2,field3,spec_rank,iteration) result(success)

    use geometric_shape_tools , only : computeCoordsSphere
    use mpilayout_tools
    use avs

    !I/O data
    type(real_data_layout),intent(in)               :: field1,field2,field3
    integer,intent(in)                              :: spec_rank
    integer,intent(in)                              :: iteration
    logical                                         :: success

    !Local data
    type(real_data_layout)                          :: increment
    real(WP)                                        :: ordre
    integer                                         :: bufI,bufJ,bufK
    integer                                         :: i,j,k,l,m
    integer                                         :: countPts
    real(WP)                                        :: lCorr,avgSca,varSca,skewSca,kurtSca
    real(WP)                                        :: t1,t2
    integer,dimension(:,:),allocatable              :: coordsPoints 
    real(wp),dimension(:),allocatable               :: currentInc1,currentInc2,currentInc3
    real(WP),dimension(:),allocatable               :: distTot
    real(WP),dimension(:),allocatable               :: avg
    !character(len=str_long)                         :: namelayout

    success = .false.

    !Lecture input
    if (parser_is_defined('Fonction structure longueur correlation')) then
      call parser_read('Fonction structure longueur correlation',lCorr)
    else
      lCorr = field1%Lx/10.0_WP
    endif
    if (parser_is_defined('Fonction structure ordre')) then
      call parser_read('Fonction structure ordre',ordre)
    else
      ordre=2
    endif
    if ( iteration == 0 .or. iteration == 1 ) then
      if (spec_rank .eq. 0 ) then
        write(6,'(a,1x,f5.2)')     '[INFO FS] ordre                    :',ordre
        write(6,'(a,1x,f5.2)')     '[INFO FS] longueur correlation     :',lCorr
      endif
    endif

    !allocataion mémoire
    if (.not. copyStructOnly(field1,increment,'increment')) return
    allocate (currentInc1(getnbcpus()),stat=i)
    if (i.ne. 0) return
    allocate (currentInc2(getnbcpus()),stat=i)
    if (i.ne. 0) return
    allocate (currentInc3(getnbcpus()),stat=i)
    if (i.ne. 0) return
    allocate (distTot(getnbcpus()),stat=i)
    if (i.ne. 0) return
    allocate (avg(getnbcpus()),stat=i)
    if (i.ne. 0) return
    if (.not. computeCoordsSphere(field1,lCorr,coordsPoints,countPts)) return
!    if (.not. computeCoordsSphere(field,lCorr,coordsPoints,countPts,0.05_WP)) return
    increment%name='increment'
    if ( spec_rank .eq. 0) write (6,'(a)') '[INFO FS] computation in progress...'
    call CPU_time(t1)
    do k=field1%zmin,field1%zmax
      do j=field1%ymin,field1%ymax
        do i=field1%xmin,field1%xmax
          if (.not. DoAllGatherScal(field1%values(i,j,k),currentInc1,spec_rank)) return
          if (.not. DoAllGatherScal(field2%values(i,j,k),currentInc2,spec_rank)) return
          if (.not. DoAllGatherScal(field3%values(i,j,k),currentInc3,spec_rank)) return
          avg = 0.0_WP
          distTot = 0.0_WP
          do l=1,countPts
            bufI = i + coordsPoints(1,l)
            bufJ = j + coordsPoints(2,l)
            bufK = k + coordsPoints(3,l)
            if ( bufI .gt. field1%nx) then
              bufI = bufI - field1%nx
            elseif ( bufI .lt. 1 ) then
              bufI = field1%nx + bufI
            endif
            if ( bufJ .gt. field1%ny) then
              bufJ = bufJ - field1%ny
            elseif ( bufJ .lt. 1 ) then
              bufJ = field1%ny + bufJ
            endif
            if ( bufK .gt. field1%nz) then
              bufK = bufK - field1%nz
            elseif ( bufK .lt. 1 ) then
              bufK = field1%nz+ bufK
            endif
            distTot(foundprocess(bufI,bufJ,bufK,field1)+1) = 1.0_WP + distTot(foundprocess(bufI,bufJ,bufK,field1)+1) 
            if (spec_rank .eq. foundprocess(bufI,bufJ,bufK,field1) ) then
              avg(foundprocess(i,j,k,field1)+1) =( (currentInc1(foundprocess(i,j,k,field1)+1) - field1%values(bufI,bufJ,bufK) ) ** ordre & 
                                               +  (currentInc2(foundprocess(i,j,k,field1)+1) - field2%values(bufI,bufJ,bufK) ) ** ordre & 
                                               +  (currentInc3(foundprocess(i,j,k,field1)+1) - field3%values(bufI,bufJ,bufK) ) ** ordre) ** (1.0_WP/ordre) & 
                                               + avg(foundprocess(i,j,k,field1)+1)
            endif
          enddo
          avg = doglobalVectorSum(spec_rank,avg)
          distTot = doglobalVectorSum(spec_rank,distTot)
          increment%values(i,j,k) = avg(spec_rank+1)/distTot(spec_rank+1)
        enddo
      enddo
    enddo
    call CPU_time(t2)
    if (spec_rank == 0) write(6,'(a,1x,f15.5)') '[INFO FS] tps calcul',t2-t1    
    avgSca  = computeFieldAvg(increment,spec_rank)
    varSca  = computeFieldVar(increment,spec_rank)
    skewSca = computeFieldSkew(increment,spec_rank)
    kurtSca = computeFieldKurt(increment,spec_rank)
    if (spec_rank .eq. 0) write(6,'(5(a,1x,f15.10,1x),a,1x,i0,1x,a,1x,f8.4)') '[INFO FS] ordre de increment',ordre,'moyenne=',avgSca,'variance=',varSca,&
                                        &'etalement',skewSca,'aplatissement',kurtSca,'nbpts',countPts,'coef',real(m,WP)*0.05_WP
    if (.not. computeFieldPDF(iteration,increment,100,spec_rank,1000,-0.8_WP,0.8_WP)) return
    deallocate(coordsPoints)
    deallocate(avg)
    deallocate(distTot)
    deallocate(currentInc1)
    deallocate(currentInc2)
    deallocate(currentInc3)
    call deletedatalayout(increment) 

    success = .true.

end function computeFunctionStructureVec


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
!>    @param[in]    me   = mpi rank inside the spectral communicator
!>    @param[inout] normp 
!>    @param[inout] Dist 
!>    @return       success 
!---------------------------------------------------------------------------
function computeFunctionStructureSca(field,spec_rank,iteration) result(success)

    use geometric_shape_tools , only : computeCoordsSphere
    use mpilayout_tools
    use avs

    !I/O data
    type(real_data_layout),intent(in)               :: field 
    integer,intent(in)                              :: spec_rank
    integer,intent(in)                              :: iteration
    logical                                         :: success

    !Local data
    type(real_data_layout)                          :: increment
    real(WP)                                        :: ordre
    integer                                         :: bufI,bufJ,bufK
    integer                                         :: i,j,k,l,m
    integer                                         :: countPts
    real(WP)                                        :: lCorr,avgSca,varSca,skewSca,kurtSca
    real(WP)                                        :: t1,t2
    integer,dimension(:,:),allocatable              :: coordsPoints 
    real(wp),dimension(:),allocatable               :: currentInc
    real(WP),dimension(:),allocatable               :: distTot
    real(WP),dimension(:),allocatable               :: avg
    !character(len=str_long)                         :: namelayout

    success = .false.

    !Lecture input
    if (parser_is_defined('Fonction structure longueur correlation')) then
      call parser_read('Fonction structure longueur correlation',lCorr)
    else
      lCorr = field%Lx/10.0_WP
    endif
    if (parser_is_defined('Fonction structure ordre')) then
      call parser_read('Fonction structure ordre',ordre)
    else
      ordre=2
    endif
    if ( iteration == 0 .or. iteration == 1 ) then
      if (spec_rank .eq. 0 ) then
        write(6,'(a,1x,f5.2)')     '[INFO FS] ordre                    :',ordre
        write(6,'(a,1x,f5.2)')     '[INFO FS] longueur correlation     :',lCorr
      endif
    endif

    !allocataion mémoire
    if (.not. copyStructOnly(field,increment,'increment')) return
    allocate (currentInc(getnbcpus()),stat=i)
    if (i.ne. 0) return
    allocate (distTot(getnbcpus()),stat=i)
    if (i.ne. 0) return
    allocate (avg(getnbcpus()),stat=i)
    if (i.ne. 0) return
    if (.not. computeCoordsSphere(field,lCorr,coordsPoints,countPts)) return
!    if (.not. computeCoordsSphere(field,lCorr,coordsPoints,countPts,0.05_WP)) return
    increment%name='increment'
    if ( spec_rank .eq. 0) write (6,'(a)') '[INFO FS] computation in progress...'
    call CPU_time(t1)
    do k=field%zmin,field%zmax
      do j=field%ymin,field%ymax
        do i=field%xmin,field%xmax
          if (.not. DoAllGatherScal(field%values(i,j,k),currentInc,spec_rank)) return
          avg = 0.0_WP
          distTot = 0.0_WP
          do l=1,countPts
            bufI = i + coordsPoints(1,l)
            bufJ = j + coordsPoints(2,l)
            bufK = k + coordsPoints(3,l)
            if ( bufI .gt. field%nx) then
              bufI = bufI - field%nx
            elseif ( bufI .lt. 1 ) then
              bufI = field%nx + bufI
            endif
            if ( bufJ .gt. field%ny) then
              bufJ = bufJ - field%ny
            elseif ( bufJ .lt. 1 ) then
              bufJ = field%ny + bufJ
            endif
            if ( bufK .gt. field%nz) then
              bufK = bufK - field%nz
            elseif ( bufK .lt. 1 ) then
              bufK = field%nz+ bufK
            endif
            distTot(foundprocess(bufI,bufJ,bufK,field)+1) = 1.0_WP + distTot(foundprocess(bufI,bufJ,bufK,field)+1) 
            if (spec_rank .eq. foundprocess(bufI,bufJ,bufK,field) ) then
              avg(foundprocess(i,j,k,field)+1) = (currentInc(foundprocess(i,j,k,field)+1) - field%values(bufI,bufJ,bufK) ) ** ordre + avg(foundprocess(i,j,k,field)+1)
            endif
          enddo
          avg = doglobalVectorSum(spec_rank,avg)
          distTot = doglobalVectorSum(spec_rank,distTot)
          increment%values(i,j,k) = avg(spec_rank+1)/distTot(spec_rank+1)
        enddo
      enddo
    enddo
    call CPU_time(t2)
    if (spec_rank == 0) write(6,'(a,1x,f15.5)') '[INFO FS] tps calcul',t2-t1    
    avgSca  = computeFieldAvg(increment,spec_rank)
    varSca  = computeFieldVar(increment,spec_rank)
    skewSca = computeFieldSkew(increment,spec_rank)
    kurtSca = computeFieldKurt(increment,spec_rank)
    if (spec_rank .eq. 0) write(6,'(5(a,1x,f15.10,1x),a,1x,i0,1x,a,1x,f8.4)') '[INFO FS] ordre de increment',ordre,'moyenne=',avgSca,'variance=',varSca,&
                                        &'etalement',skewSca,'aplatissement',kurtSca,'nbpts',countPts,'coef',real(m,WP)*0.05_WP
    if (.not. computeFieldPDF(iteration,increment,100,spec_rank,1000,-0.8_WP,0.8_WP)) return
    deallocate(coordsPoints)
    deallocate(avg)
    deallocate(distTot)
    deallocate(currentInc)
    call deletedatalayout(increment) 

    success = .true.

end function computeFunctionStructureSca


!!---------------------------------------------------------------------------
!!> @details 
!!> compute the distances between the point P (coordPoint) and all the others
!!> and the normp (normL2**ordre) between the vector at coordPoint (VecRef) and  
!!> the other vectors (vecField). Note you can compute normp between vectors of dimension 1 to N
!!> @author = Antoine Vollant, LEGI
!!>    @param[in]    coordPoint = scalar field in physical space
!!>    @param[in]    VecRef 
!!>    @param[in]    vecField
!!>    @param[in]    ordre 
!!>    @param[in]    me   = mpi rank inside the spectral communicator
!!>    @param[inout] normp 
!!>    @param[inout] Dist 
!!>    @return       success 
!!---------------------------------------------------------------------------
!function computeFunctionStructureSca(field,spec_rank,iteration) result(success)
!
!    use mpilayout_tools
!    use avs
!    use conditional_mean_tools , only : computeEONVar
!
!    !I/O data
!    type(real_data_layout),intent(in)               :: field 
!    integer,intent(in)                              :: spec_rank
!    integer,intent(in)                              :: iteration
!    logical                                         :: success
!
!    !Local data
!    type(real_data_layout)                          :: dist 
!    type(real_data_layout)                          :: normp 
!!    type(real_data_layout)                          :: increment1
!    type(real_data_layout)                          :: increment2
!    real(WP)                                        :: ordre
!    integer                                         :: rootSum,root
!    integer                                         :: i,j,k,ires 
!    integer                                         :: l,m,n
!    real(WP),dimension(1)                           :: VecRefBuf
!    integer                                         :: countPts
!    real(WP)                                        :: lCorr,avg,distTot
!    real(WP)                                        :: var,skew,kurt,t1,t2
!
!    success = .false.
!
!    !allocataion mémoire
!    if (.not. copyStructOnly(field,dist,'dist') .or. &
!       &.not. copyStructOnly(field,normp,'normp') .or. &
!!       &.not. copyStructOnly(field,increment1,'increment1') .or. &
!       &.not. copyStructOnly(field,increment2,'increment2')) then 
!      write(6,'(a)') '[ERROR] in computeFunctionStructure : copy structure 1 failed'
!      return
!    endif
!
!    !Lecture input
!    if (parser_is_defined('Fonction structure longueur correlation')) then
!      call parser_read('Fonction structure longueur correlation',lCorr)
!    else
!      lCorr = field%Lx/2.0_WP
!    endif
!    if (parser_is_defined('Fonction structure ordre')) then
!      call parser_read('Fonction structure ordre',ordre)
!    else
!      ordre=2
!    endif
!    if ( iteration == 0 ) then
!      if (spec_rank .eq. 0 ) then
!        write(6,'(a,1x,f5.2)')     '[INFO FS] ordre                    :',ordre
!        write(6,'(a,1x,f5.2)')     '[INFO FS] longueur correlation     :',lCorr
!      endif
!    endif
!
!    countPts = 0
!    if ( spec_rank .eq. 0) write (6,'(a)') '[INFO FS] computation in progress...'
!    call CPU_time(t1)
!    do k=1,field%nz
!      do j=1,field%ny
!        do i=1,field%nx
!          root = 0
!          if ( (i .ge. field%xmin) .and. (i .le. field%xmax) .and. &
!             & (j .ge. field%ymin) .and. (j .le. field%ymax) .and. &
!             & (k .ge. field%zmin) .and. (k .le. field%zmax)) then
!            root = spec_rank
!            VecRefBuf = field%values(i,j,k)
!          endif
!          rootSum = doglobalsum(spec_rank,root)
!          if (.not. DoBcastVector(VecRefBuf,1,rootSum,spec_rank)) then 
!            write (6,'(a)') '[ERROR] in computeFunctionStructure : DoBcastVector'
!            return
!          endif
!
!!          if ( .not. computeIncrementSca((/i,j,k/),VecRefBuf(1),field,normp,dist,ordre,spec_rank)) then
!!            write (6,'(a)') '[ERROR] in computeFunctionStructure : computeIncrement'
!!            return
!!          endif
!!          avg = 0.0_WP
!!          distTot = 0.0_WP
!!          do l=normp%zmin,normp%zmax
!!            do m=normp%ymin,normp%ymax
!!              do n=normp%xmin,normp%xmax
!!                if ( dist%values(n,m,l) .gt. lCorr - sqrt(3.0_WP)*normp%Lx/normp%nx .and. dist%values(n,m,l) .lt. lCorr + sqrt(3.0_WP)*normp%Lx/normp%nx ) then
!!                  avg = normp%values(n,m,l)/(abs(lCorr-dist%values(n,m,l))) + avg
!!                  distTot = 1.0_WP/(abs(lCorr-dist%values(n,m,l))) + distTot 
!!                endif
!!              enddo
!!            enddo
!!          enddo
!!          avg = doglobalsum(spec_rank,avg)
!!          distTot = doglobalsum(spec_rank,distTot)
!!          if ( (i .ge. field%xmin) .and. (i .le. field%xmax) .and. &
!!             & (j .ge. field%ymin) .and. (j .le. field%ymax) .and. &
!!             & (k .ge. field%zmin) .and. (k .le. field%zmax)) then
!!            increment1%values(i,j,k) = avg/distTot
!!          endif
!
!          if ( .not. computeIncrementSca((/i,j,k/),VecRefBuf(1),field,normp,dist,ordre,spec_rank,Rsphere=lCorr)) then
!            write (6,'(a)') '[ERROR] in computeFunctionStructure : computeIncrement'
!            return
!          endif
!          avg = 0.0_WP
!          distTot = 0.0_WP
!          do l=normp%zmin,normp%zmax
!            do m=normp%ymin,normp%ymax
!              do n=normp%xmin,normp%xmax
!                if ( dist%values(n,m,l) .gt. -1e-15 ) then
!                  avg = normp%values(n,m,l)/(abs(lCorr-dist%values(n,m,l))) + avg
!                  distTot = 1.0_WP/(abs(lCorr-dist%values(n,m,l))) + distTot 
!                endif
!              enddo
!            enddo
!          enddo
!          avg = doglobalsum(spec_rank,avg)
!          distTot = doglobalsum(spec_rank,distTot)
!          if ( (i .ge. field%xmin) .and. (i .le. field%xmax) .and. &
!             & (j .ge. field%ymin) .and. (j .le. field%ymax) .and. &
!             & (k .ge. field%zmin) .and. (k .le. field%zmax)) then
!            increment2%values(i,j,k) = avg/distTot
!          endif
!
!          countPts=countPts+1
!        enddo
!      enddo
!    enddo
!    call CPU_time(t2)
!    if (spec_rank == 0) write(6,'(a,1x,f15.5)') '[INFO FS] tps calcul',t2-t1
!
!!    avg  = computeFieldAvg(increment1,spec_rank)
!!    var  = computeFieldVar(increment1,spec_rank)
!!    skew = computeFieldSkew(increment1,spec_rank)
!!    kurt = computeFieldKurt(increment1,spec_rank)
!!    if (spec_rank .eq. 0) write(6,'(5(a,1x,f15.10,1x))') '[INFO FS] ordre de increment',ordre,'moyenne=',avg,'variance=',var,&
!!                                        &'etalement',skew,'aplatissement',kurt
!!    if (.not. computeFieldPDF(iteration,increment1,100,spec_rank,1000)) return
!
!    avg  = computeFieldAvg(increment2,spec_rank)
!    var  = computeFieldVar(increment2,spec_rank)
!    skew = computeFieldSkew(increment2,spec_rank)
!    kurt = computeFieldKurt(increment2,spec_rank)
!    if (.not. computeFieldPDF(iteration,increment2,100,spec_rank,1000)) return
!    if (spec_rank .eq. 0) write(6,'(5(a,1x,f15.10,1x))') '[INFO FS] ordre de increment',ordre,'moyenne=',avg,'variance=',var,&
!                                        &'etalement',skew,'aplatissement',kurt
!
!    call deletedatalayout(normp)
!    call deletedatalayout(dist)
!!    call deletedatalayout(increment1)
!    call deletedatalayout(increment2)
!
!    success = .true.
!
!end function computeFunctionStructureSca

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
!>    @param[in]    me   = mpi rank inside the spectral communicator
!>    @param[inout] normp
!>    @param[inout] Dist
!>    @return       success
!---------------------------------------------------------------------------
function computeFunctionStructureV1(field,spec_rank,iteration) result(success)

    use mpilayout_tools
    use avs

    !I/O data
    type(real_data_layout),dimension(:),intent(in)  :: field 
    integer,intent(in)                              :: spec_rank
    integer,intent(in)                              :: iteration
    logical                                         :: success

    !Local data
    type(real_data_layout)                          :: dist 
    type(real_data_layout)                          :: normp 
    type(real_data_layout)                          :: increment 
    real(WP)                                        :: ordre
    integer                                         :: rootSum,root,sampling
    integer                                         :: i,j,k,l,m,n,ires 
    integer,dimension(6)                            :: coordRang
    integer                                         :: sizeOfVecRef
    real(WP),dimension(:),allocatable               :: VecRefBuf
    character(len=6)                                :: nameordre,nameite
    character(len=str_long)                         :: filename
    integer                                         :: file_id,indice,countPts
    logical                                         :: logLongRang,existe
    real(WP)                                        :: lCorr,avg,distTot
    integer,dimension(:),allocatable                :: nombreLoc,nombreglob
    real(WP),dimension(:),allocatable               :: valueLoc,valueglob
    logical                                         :: dumpFile,dumpField

    success = .false.

    coordRang = 0
    logLongRang = .false.
    sizeOfVecRef = size(field,1)
    dumpFile = .true.
    dumpField = .true.
    !allocataion mémoire
    allocate(VecRefBuf(sizeOfVecRef),stat=ires)
    if (.not. copyStructOnly(field(1),dist,'dist') .or. &
       &.not. copyStructOnly(field(1),normp,'normp') .or. &
       &.not. copyStructOnly(field(1),increment,'increment')) then 
      write(6,'(a)') '[ERROR] in computeFunctionStructure : copy structure 1 failed'
      return
    endif

    !Lecture input
    if (parser_is_defined('Fonction structure limite points inf') .and. &
       &parser_is_defined('Fonction structure limite points sup')) then
      if (.not. GetNxNyNz('Fonction structure limite points inf',coordRang(1),coordRang(2),coordRang(3))) return
      if (.not. GetNxNyNz('Fonction structure limite points sup',coordRang(4),coordRang(5),coordRang(6))) return
    else
      coordRang(1) = (field(1)%xmax+field(1)%xmin)/2 
      coordRang(2) = (field(1)%ymax+field(1)%ymin)/2
      coordRang(3) = (field(1)%zmax+field(1)%zmin)/2
      coordRang(4) = coordRang(1)
      coordRang(5) = coordRang(2)
      coordRang(6) = coordRang(3)
    endif
    if (parser_is_defined('Fonction structure longueur correlation')) then
      call parser_read('Fonction structure longueur correlation',lCorr)
    elseif (parser_is_defined('Fonction structure longueur correlation max')) then
      call parser_read('Fonction structure longueur correlation max',lCorr)
      logLongRang = .true.
    else
      lCorr = field(1)%Lx/2.0_WP
    endif
    if (parser_is_defined('Fonction structure ordre')) then
      call parser_read('Fonction structure ordre',ordre)
    else
      ordre=2
    endif
    if ( iteration == 0 ) then
      if (spec_rank .eq. 0 ) then
        write(6,'(a,1x,f5.2)')     '[INFO FS] ordre                    :',ordre
        write(6,'(a,1x,2(i0,1x))') '[INFO FS] coords points X          :',coordRang(1),coordRang(4)
        write(6,'(a,1x,2(i0,1x))') '[INFO FS] coords points Y          :',coordRang(2),coordRang(5)
        write(6,'(a,1x,2(i0,1x))') '[INFO FS] coords points Z          :',coordRang(3),coordRang(6)
        write(6,'(a,1x,f5.2)')     '[INFO FS] longueur correlation     :',lCorr
      endif
    endif

    sampling = floor(lCorr*field(1)%nx/field(1)%Lx) 
    if ( logLongRang ) then 
      allocate (nombreLoc(sampling),stat=ires)
      if (ires .ne. 0 ) return
      allocate (valueLoc(sampling),stat=ires)
      if (ires .ne. 0 ) return
      allocate (nombreglob(sampling),stat=ires)
      if (ires .ne. 0 ) return
      allocate (valueglob(sampling),stat=ires)
      if (ires .ne. 0 ) return
    endif
    
    do i=0,1
      if ( (coordRang(1+i*3) .lt. 1) .or. (coordRang(1+i*3) .gt. field(1)%nx) .or. &
         & (coordRang(2+i*3) .lt. 1) .or. (coordRang(2+i*3) .gt. field(1)%ny) .or. &
         & (coordRang(3+i*3) .lt. 1) .or. (coordRang(3+i*3) .gt. field(1)%nz)) then
        write (6,'(a)') '[ERROR] in computeFunctionStructure : points asked is out of domain'
        return
      endif
    enddo

    if (logLongRang ) then
      nombreglob = 0
      valueglob = 0.0_WP
    endif
    countPts = 0
    do k=coordRang(3),coordRang(6)
      do j=coordRang(2),coordRang(5)
        do i=coordRang(1),coordRang(4)
          root = 0
          if ( (i .ge. field(1)%xmin) .and. (i .le. field(1)%xmax) .and. &
             & (j .ge. field(1)%ymin) .and. (j .le. field(1)%ymax) .and. &
             & (k .ge. field(1)%zmin) .and. (k .le. field(1)%zmax)) then
            root = spec_rank
            do l=1,sizeOfVecRef
              VecRefBuf(l) = field(l)%values(i,j,k)
            enddo
          endif
          rootSum = doglobalsum(spec_rank,root)
          if (.not. DoBcastVector(VecRefBuf,sizeOfVecRef,rootSum,spec_rank)) then 
            write (6,'(a)') '[ERROR] in computeFunctionStructure : DoBcastVector'
            return
          endif
          if ( logLongRang ) then
            if ( .not. computeIncrementVec((/i,j,k/),VecRefBuf,field,normp,dist,ordre,spec_rank,Rboule=lCorr)) then
              write (6,'(a)') '[ERROR] in computeFunctionStructure : computeIncrement'
              return
            endif
            nombreloc = 0
            valueloc = 0.0_WP
            do l=normp%zmin,normp%zmax
              do m=normp%ymin,normp%ymax
                do n=normp%xmin,normp%xmax
                  if ( dist%values(n,m,l) .gt. 0.0_WP -1e-15 .and. dist%values(n,m,l) .lt. lCorr + 1e-15 ) then
                    indice = computeIndice(0.0_WP,lCorr,dist%values(n,m,l),sampling)
                    nombreloc(indice) = nombreloc(indice) + 1
                    valueloc(indice) = valueloc(indice) + normp%values(n,m,l)
                  endif
                enddo
              enddo
            enddo
            nombreloc = doglobalVectorsum(spec_rank,nombreloc)
            valueloc = doglobalVectorsum(spec_rank,valueloc)
            nombreglob = nombreglob + nombreloc
            valueglob = valueglob + valueloc
          else 

            !if ( .not. computeIncrementVec((/i,j,k/),VecRefBuf,(/field/),normp,dist,ordre,spec_rank,Rsphere=lCorr)) then
            if ( .not. computeIncrementVec((/i,j,k/),VecRefBuf,(/field/),normp,dist,ordre,spec_rank)) then
              write (6,'(a)') '[ERROR] in computeFunctionStructure : computeIncrement'
              return
            endif
            if (.not. computeFieldPDF(iteration,dist,100,spec_rank,1000,0.3_WP,0.8_WP)) return
            avg = 0.0_WP
            distTot = 0.0_WP
            do l=normp%zmin,normp%zmax
              do m=normp%ymin,normp%ymax
                do n=normp%xmin,normp%xmax
                  if ( dist%values(n,m,l) .gt. lCorr - sqrt(3.0_WP)*normp%Lx/normp%nx .and. dist%values(n,m,l) .lt. lCorr + sqrt(3.0_WP)*normp%Lx/normp%nx ) then
                  !if ( dist%values(n,m,l) .gt. -1e-15 ) then
                    avg = normp%values(n,m,l)/(abs(lCorr-dist%values(n,m,l))) + avg
                    distTot = 1.0_WP/(abs(lCorr-dist%values(n,m,l))) + distTot 
                  endif
                enddo
              enddo
            enddo
            avg = doglobalsum(spec_rank,avg)
            distTot = doglobalsum(spec_rank,distTot)
            avg = avg/distTot
            print *,'Fonction structure 1',avg

            if ( .not. computeIncrementVec((/i,j,k/),VecRefBuf,(/field/),normp,dist,ordre,spec_rank,Rsphere=lCorr)) then
              write (6,'(a)') '[ERROR] in computeFunctionStructure : computeIncrement'
              return
            endif
            if (.not. computeFieldPDF(iteration,dist,100,spec_rank,1000,0.3_WP,0.8_WP)) return
            avg = 0.0_WP
            distTot = 0.0_WP
            do l=normp%zmin,normp%zmax
              do m=normp%ymin,normp%ymax
                do n=normp%xmin,normp%xmax
                  if ( dist%values(n,m,l) .gt. -1e-15 ) then
                    avg = normp%values(n,m,l)/(abs(lCorr-dist%values(n,m,l))) + avg
                    distTot = 1.0_WP/(abs(lCorr-dist%values(n,m,l))) + distTot 
                  endif
                enddo
              enddo
            enddo
            avg = doglobalsum(spec_rank,avg)
            distTot = doglobalsum(spec_rank,distTot)
            avg = avg/distTot
            print *,'Fonction structure 2',avg


          endif
          if ( dumpFile ) then
            if (spec_rank .eq. 0 ) then
              file_id = iopen()
              write(nameordre,'(f5.2)') ordre
              write(nameite,'(i6.6)') iteration
              write(filename,'(a)') 'FSloc_'//trim(field(1)%name)//'_ordre-'//trim(adjustl(nameordre))//'_ite-'//trim(adjustl(nameite))//'.table'
              inquire( file=trim(adjustl(filename)), exist=existe)
              open(unit=file_id, file=trim(adjustl(filename)), form="FORMATTED",position="append")
              if ( logLongRang ) then
                if (.not. existe) then
                  write(file_id, '(a)')  '#i j k distance Sp(X,distance) pdf(Sp)'
                endif
                do l=1,sampling
                  write (file_id,'(4(i0,1x),2(f20.8,1x),1x,i0)') countPts,i,j,k,lCorr/(sampling)*real((l-0.5),WP),valueloc(l)/real(nombreloc(l),WP),nombreloc(l)
                enddo
                write (file_id,'(a)') '' 
              else
                if (.not. existe) then
                  write(file_id, '(a)')  '#i j k distance Sp(X,distance)'
                endif
                write (file_id,'(4(i0,1x),2(f20.8,1x),1x,i0)') countPts,i,j,k,lCorr,avg
              endif
              close(file_id)
              file_id = iclose(file_id)
            endif
          endif
          countPts=countPts+1
        enddo
      enddo
    enddo
    if ( dumpFile .and. logLongRang) then
      if (spec_rank .eq. 0 ) then
        file_id = iopen()
        write(nameordre,'(f5.2)') ordre
        write(nameite,'(i6.6)') iteration
        write(filename,'(a)') 'FSglob_'//trim(field(1)%name)//'_ordre-'//trim(adjustl(nameordre))//'_ite-'//trim(adjustl(nameite))//'.table'
        inquire( file=trim(adjustl(filename)), exist=existe)
        open(unit=file_id, file=trim(adjustl(filename)), form="FORMATTED",position="append")
        if (.not. existe) then
          write(file_id, '(a)')  '#distance Sp(X,distance) pdf(Sp)'
        endif
        do l=1,sampling
          write (file_id,'(2(f20.8,1x),1x,i0)') lCorr/(sampling)*real(l-0.5,WP),valueglob(l)/real(nombreglob(l),WP),nombreglob(l)
        enddo
        close(file_id)
        file_id = iclose(file_id)
      endif
    endif

    if ( dumpField ) then
      if (.not. dump_UVWP_ForAVS(trim(normp%name),spec_rank,normp%lx,normp%ly,normp%lz,normp,usename=.true.)) then
        write (6,'(a)') '[ERROR] in computeFunctionStructure : dump_UVWP_ForAVS'
        return
      endif
      if (.not. dump_UVWP_ForAVS(trim(dist%name),spec_rank,dist%lx,dist%ly,dist%lz,dist,usename=.true.)) then
        write (6,'(a)') '[ERROR] in computeFunctionStructure : dump_UVWP_ForAVS'
        return
      endif
    endif

    deallocate(VecRefBuf)
    if (logLongRang ) then
      deallocate(nombreloc)
      deallocate(valueloc)
      deallocate(nombreglob)
      deallocate(valueglob)
    endif
    call deletedatalayout(normp)    
    call deletedatalayout(dist)    
    call deletedatalayout(increment) 

    success = .true.

end function computeFunctionStructureV1

!------------------------------------------------------------------------------
!> Compute several moments of fields and its derivation field
!! @author Antoine Vollant
!!    @param[in]    field           = field from wich spectrum has to been computed.
!!    @param[in]    iteration       = time iteration number
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @return       success         = logical equals to true is everything is right.
!------------------------------------------------------------------------------
function computeScalarMomentsStat(field, iteration, spec_rank) result(success)


  ! Input/Output
  type(REAL_DATA_LAYOUT), intent(in):: field
  integer, intent(in)               :: iteration
  integer, intent(in)               :: spec_rank
  logical                           :: success

  type(REAL_DATA_LAYOUT)            :: deriv
  integer                           :: file_id
  character(len=str_long)           :: filename

  real(WP),dimension(:),allocatable :: ordiField,centField,cRedField
  real(WP),dimension(:),allocatable :: ordiDeriv,centDeriv,cRedDeriv
  integer                           :: i,order,maxOrder,ires
  character(len=7)                  :: nameite

  success = .false.

  ! -- More advanced stats --
  ! Init datalayout
  if (.not. copyStructOnly(field,deriv) ) then
    if (spec_rank .eq. 0) write(6,'(a)') '[ERROR] Post-process: in compute_scalar_Basic_Stat copy structure failed'
    return
  endif
  deriv%name=trim(field%name)//'_dX'
  ! Compute scalar first derivate along X
  if (.not. computeDerivationField(field,(/1/),deriv,spec_rank)) then
    if (spec_rank .eq. 0) write(6,'(a)') '[ERROR] Post-process: in compute_scalar_Basic_Stat scalar derivation failed'
    return
  end if

  maxorder = 40
  allocate(ordiField(maxorder),stat=ires)
  if (ires .ne. 0) return
  allocate(centField(maxorder),stat=ires)
  if (ires .ne. 0) return
  allocate(cRedField(maxorder),stat=ires)
  if (ires .ne. 0) return
  allocate(ordiDeriv(maxorder),stat=ires)
  if (ires .ne. 0) return
  allocate(centDeriv(maxorder),stat=ires)
  if (ires .ne. 0) return
  allocate(cRedDeriv(maxorder),stat=ires)
  if (ires .ne. 0) return
  cRedField = 0.0_WP
  cRedDeriv = 0.0_WP
  do order = 1,maxorder
    ordiField(order) = computeOrdinaryMoment(field,order,spec_rank)
    centField(order) = computeCenteredMoment(field,order,spec_rank)
    cRedField(order) = computeReduceCenteredMoment(field,order,spec_rank)
    ordiDeriv(order) = computeOrdinaryMoment(deriv,order,spec_rank)
    centDeriv(order) = computeCenteredMoment(deriv,order,spec_rank)
    cRedDeriv(order) = computeReduceCenteredMoment(deriv,order,spec_rank)
  enddo

  if (spec_rank .eq. 0 ) then
    file_id = iopen()
    write(nameite,'(i6.6)') iteration
    write(filename,'(a)') 'statMoment_'//trim(adjustl(field%name))//'_'//trim(adjustl(nameite))//'.table'
    open(unit=file_id, file=trim(adjustl(filename)), form="FORMATTED")
    write(file_id, '(a84)')  '#order ordiSca centerSca reduced_centerSca ordiDeriv centerDeriv reduced_centerSca'
    do i=1,maxorder
      write (file_id,'(i0,6(1x,g20.8))') &
          & i,ordiField(i),centField(i),cRedField(i),ordiDeriv(i),centDeriv(i),cRedDeriv(i)
    enddo
    close(file_id)
    file_id = iclose(file_id)
  endif

  if (.not. computeFunctionStructureSca(field,spec_rank,iteration)) then
    write (6,'(a)') '[ERROR] in computeFunctionStructure'
    return
  endif

  ! Free memory
  call deletedatalayout(deriv)
  deallocate(ordiField)
  deallocate(centField)
  deallocate(cRedField)
  deallocate(ordiDeriv)
  deallocate(centDeriv)
  deallocate(cRedDeriv)

  success = .true.

end function computeScalarMomentsStat


!------------------------------------------------------------------------------
!> Compute basic statististic on scalar distribution (momentum of different orders) and save
!! their evolution into a field.
!! spectrum).
!! @author Jean-Baptiste Lagaert
!!    @param[in]    field           = field from wich spectrum has to been computed.
!!    @param[in]    ite             = time iteration number
!!    @param[in]    spec_rank   = mpi rank inside spectral communicator
!!    @param[in]    nbcpus          = number of processes
!!    @return       success         = logical equals to true is everything is right.
!! @details
!!        This function compute the spectrum of a given field. It perform
!!    a FFT in order to obtain it in the Fourrier space (wich is required
!!    to computing the spectrum). It could be perform "out a simulation"
!!    (ie from a field write in an output).
!------------------------------------------------------------------------------
function compute_scalar_basic_stat(field, U, diffu, iteration, ite_max, time, nbcpus, spec_rank) result(success)


  ! Input/Output
  type(REAL_DATA_LAYOUT), intent(in):: field, U
  real(WP), intent(in)              :: time
  real(WP), intent(in)              :: diffu
  integer, intent(in)               :: iteration, ite_max
  integer, intent(in)               :: nbcpus
  integer, intent(in)               :: spec_rank
  logical                           :: success

  type(REAL_DATA_LAYOUT)            :: deriv
  real(WP)                          :: energ, mini, maxi, avg, var, skew, kurt
  real(WP)                          :: sk_grad, flat_grad,mix_skew, dissip_rate
  real(WP)                          :: lambda, grad_var
  integer                           :: file_id
  logical                           :: existe
  character(len=str_long)           :: filename

  success = .false.

  ! -- First stats --
  energ   = computeEnergy( field ,spec_rank,nbcpus)
  mini    = computeFieldMin( field , spec_rank )
  maxi    = computeFieldMax( field , spec_rank )

  ! -- 1-4 order moment --
  avg     = computeFieldAvg (field,spec_rank)
  var     = computeFieldVar (field,spec_rank)
  skew    = computeFieldSkew(field,spec_rank)
  kurt    = computeFieldKurt(field,spec_rank)

  ! -- More advanced stats --
  ! Init datalayout
  if (.not. copyStructOnly(field,deriv) ) then
    if (spec_rank .eq. 0) write(6,'(a)') '[ERROR] Post-process: in compute_scalar_Basic_Stat copy structure failed'
    return
  endif
  deriv%name=trim(field%name)//'_dX'
  ! Compute scalar first derivate along X
  if (.not. computeDerivationField(field,(/1/),deriv,spec_rank)) then
    if (spec_rank .eq. 0) write(6,'(a)') '[ERROR] Post-process: in compute_scalar_Basic_Stat scalar derivation failed'
    return
  end if
  ! Compute skewness of scalar gradient (flatness = kurtosis)
  sk_grad  = computeFieldSkew(deriv,spec_rank)
  ! Compute flatness of scalar gradient (flatness = kurtosis)
  !flat_grad = 3._WP + computeFieldKurt(deriv,spec_rank)
  ! Compute mixed skewness of scalar gradient
  if (.not. computeMixedSkew(mix_skew,U,field,spec_rank,deriv) ) then
    write (6,'(a)') '[ERROR] failed to compute mixed skewness of scalar gradient'
    return
  endif

  ! -- Compute some pdf -- (before changing deriv values)
  ! Scalar pdf
  if (.not.computeFieldPDF(iteration,field,300,spec_rank,ite_max)) then
    write(*,'(a)') '[ERROR] Failed to compute PDF of scalars fiels'
    return
  end if
  ! Pdf of d(scalar)/dx
  if (.not.computeFieldPDF(iteration,deriv,600,spec_rank,ite_max)) then
    write(*,'(a)') '[ERROR] Failed to compute PDF of scalars first derivate along X'
    return
  end if

  ! -- Some other stat -- (now we can change deriv%values)
  ! Compute gradient variance
  grad_var     = computeFieldVar (deriv,spec_rank)
  deriv%values = deriv%values**2
  !grad_var = computeFieldAvg(deriv , spec_rank , nbcpus )
  ! Compute flatness of scalar gradient (flatness = kurtosis)
  deriv%values = deriv%values**2  ! thanks to line 1461, deriv =(dU/dx)**4
  flat_grad = computeFieldAvg(deriv, spec_rank, nbcpus)/(grad_var**2)
  ! Compute lambda, the taylor microscale for scalar field
  lambda = sqrt(var/grad_var)

  ! -- Compute dissipation rate --
  ! For memory optimisation, "deriv" is re-use to save sum[(d(field)/dx_i), i=1,3]
  if(.not. computeSumSquareFirstDerivateField(field,deriv,spec_rank)) then
    write (6,'(a)') '[ERROR] failed to compute mixed skewness of scalar gradient'
    return
  endif
  deriv%name=trim(field%name)//'_dissip'
  ! Then compute the average of the square power
  dissip_rate = computeFieldAvg(deriv, spec_rank )
  ! Compute the pdf
  deriv%values = deriv%values/dissip_rate
  dissip_rate = 2._WP*diffu*dissip_rate
  if (.not.computeFieldPDF(iteration,deriv,300,spec_rank,ite_max,0._WP,150._WP)) then
    write(*,'(a)') '[ERROR] Failed to compute PDF of scalars fiels'
    return
  end if


  ! Free memory
  call deletedatalayout(deriv)

  !if (.not. compute_spectrum(field, iFic spec_rank) ) then
  !  write(*,'(a)') '[ERROR] Failed to compute spectrum of scalars fiels'
  !  return
  !end if

  if (spec_rank .eq. 0 ) then
    file_id = iopen()
    write(filename,'(a)') 'stat_'//trim(field%name)//'.table'
    inquire( file=trim(adjustl(filename)), exist=existe)
    open(unit=file_id, file=trim(adjustl(filename)), form="FORMATTED",position="append")
    if (.not. existe) then
      write(file_id, '(a6,x,14(a14,1x))')  '#1-ite', '2-time', '3-min', &
                                & '4-avg', '5-max', '6-zz/2',           &
                                & '6-var', '7-skew', '8-kurt',          &
                                & '9-SkGrad','10-GradFlat','11-MixSkew',&
                                & '12-DissipRate', '13-lambda', '14-VarGrad'
    endif
    if (iteration<1e6) then
      write (file_id,'(i6,1x,14(f14.7,1x))') &
          & iteration,time,mini,        &
          & avg,maxi,energ,             &
          & var,skew,kurt,              &
          & sk_grad, flat_grad,mix_skew,&
          & dissip_rate, lambda, grad_var
    else
      write (file_id,'(i0,1x,14(f14.7,1x))') &
          & iteration,time,mini,        &
          & avg,maxi,energ,             &
          & var,skew,kurt,              &
          & sk_grad, flat_grad,mix_skew,&
          & dissip_rate, lambda, grad_var
    end if
    close(file_id)
    file_id = iclose(file_id)
  endif

  success = .true.

end function compute_scalar_basic_stat


!-------------------------------------------------------------------------------------------
!> Compute the SGS term and the viscous term of scalar equation and it dumps results in ASCII
!! files
!! @author Antoine Volllant
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in]    nbcpus          = number of process in the cpu pool
!!    @param[in]    U               = Component of velocity in physical space not filtered
!!    @param[in]    V               = Component of velocity in physical space not filtered
!!    @param[in]    W               = Component of velocity in physical space not filtered
!!    @param[in]    Scal            = Scalar to compute  in physical space not filtered
!!    @param[inout] dTidxi          = Exact SGS terme for convection/diffusion equation of scalar
!!    @return       res             = logical value
!-------------------------------------------------------------------------------------------
 function computedTidxi(nbcpus,nd,filter,spec_rank, Scal, U, V ,W , dTidxi, Ti ) result (res)


  !I/O data
  type(real_data_layout),intent(in)      :: U, V ,W, Scal
  type(real_data_layout),intent(inout)   :: dTidxi
  type(real_data_layout),intent(inout),dimension(:),optional :: Ti
  integer,intent(in)                     :: nbcpus
  integer,intent(in)                     :: spec_rank
  integer,intent(in)                     :: filter,nd
  logical                                :: res

  !Local data
  type(real_data_layout),dimension(:),allocatable    :: temp
  type(wavenumbers)                                  :: ScalWaveN
  real(WP)                                           :: deltx

  res = .false.

  !Initialisation des tableaux de travail
  if (.not. initWorkArray(U,7,temp)) then
     write(6,'(a,i0,a)')'[ERROR] in computedTidxi : not enought memory!'
     return
  endif

  if (.not. initWN(ScalWaveN,spec_rank,Scal%nx,Scal%ny,Scal%nz)) return
  if (.not. computeWN(ScalWaveN,spec_rank,Scal%Lx,Scal%Ly,Scal%Lz,Scal%nx,Scal%ny,Scal%nz)) return
  deltx = computeDelta(U)

  temp(1)%values=U%values * Scal%values
  temp(2)%values=V%values * Scal%values
  temp(3)%values=W%values * Scal%values
  call computeFilter(ScalWaveN,nd*deltx,temp(1),temp(5),filter) !UZ bar
  call computeFilter(ScalWaveN,nd*deltx,temp(2),temp(6),filter) !VZ bar
  call computeFilter(ScalWaveN,nd*deltx,temp(3),temp(7),filter) !WZ bar
  call computeFilter(ScalWaveN,nd*deltx,U,temp(1),filter) !U bar
  call computeFilter(ScalWaveN,nd*deltx,V,temp(2),filter) !V bar
  call computeFilter(ScalWaveN,nd*deltx,W,temp(3),filter) !W bar
  call computeFilter(ScalWaveN,nd*deltx,Scal,temp(4),filter) !Z bar

       
  temp(1)%values  = temp(5)%values - temp(1)%values * temp(4)%values
  temp(2)%values  = temp(6)%values - temp(2)%values * temp(4)%values
  temp(3)%values  = temp(7)%values - temp(3)%values * temp(4)%values

  if ( present(Ti)) then
    if ( size(Ti,1) .eq. 3) then
      if ( samelayout(Ti(1),U) ) then
        Ti(1)%values = temp(1)%values
        Ti(2)%values = temp(2)%values
        Ti(3)%values = temp(3)%values
      else
        write(6,'(a,i0,a)')'[ERROR] in computedTidxi : Ti is different from U!'
        return
      endif
    else
      write(6,'(a,i0,a)')'[ERROR] in computedTidxi : Ti is not a vector!'
      return
    endif
  endif

  if (.not. computeFieldDivergence(ScalWaveN,temp(1),temp(2),temp(3),dTidxi,nbcpus,spec_rank)) return

  if ( .not. deleteWorkArray(temp) .or. &
     & .not. deleteWN(ScalWaveN)) then
   write (6,'(a)') '[ERROR] Memory leackage'
   return
  endif

  res = .true.

 end function computedTidxi


!-------------------------------------------------------------------------------------------
!> Compute the SGS term and the viscous term of scalar equation and it dumps results in ASCII
!! files
!! @author Antoine Volllant
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in]    nbcpus          = number of process in the cpu pool
!!    @param[in]    U               = Component of velocity
!!    @param[in]    V               = Component of velocity
!!    @param[in]    W               = Component of velocity
!!    @param[in]    Scal            = Scalar to compute
!!    @param[in]    ite             = the number of file
!!    @return       res             = logical value
!-------------------------------------------------------------------------------------------
 function computeClarkTerms(U,V,W,ScalIn,filter,nd,nbcpus,spec_rank,Te1,Te2) result(success)

  !I/O data
  type(real_data_layout),intent(in)        :: U,V,W,ScalIn
  type(real_data_layout),intent(inout)     :: Te1,Te2
  integer,intent(in)                       :: nbcpus,spec_rank
  integer,intent(in)                       :: filter,nd
  logical                                  :: success

  !Local data
  type(real_data_layout),dimension(:),allocatable    :: Ubar
  type(real_data_layout),dimension(:),allocatable    :: Uhatbar
  type(real_data_layout),dimension(:),allocatable    :: L
  type(real_data_layout),dimension(:),allocatable    :: UbarZbar
  type(real_data_layout),dimension(:),allocatable    :: UbarZbarTilde
  type(real_data_layout),dimension(:),allocatable    :: K
  type(real_data_layout),dimension(:),allocatable    :: K1
  type(real_data_layout),dimension(:),allocatable    :: N
  type(real_data_layout),dimension(:),allocatable    :: M
  type(real_data_layout),dimension(:),allocatable    :: Scal
  type(real_data_layout),dimension(:),allocatable    :: SdZ
  type(complex_data_layout),dimension(:),allocatable :: K1K
  type(complex_data_layout),dimension(:),allocatable :: K2K
  type(complex_data_layout),dimension(:),allocatable :: UK
  type(complex_data_layout),dimension(:),allocatable :: UKbar
  type(complex_data_layout),dimension(:),allocatable :: UbarZbarK
  type(complex_data_layout),dimension(:),allocatable :: UbarZbarKTilde
  type(complex_data_layout),dimension(:),allocatable :: UKhatbar
  type(complex_data_layout),dimension(:),allocatable :: ScalK
  type(real_data_layout),dimension(:),allocatable    :: temp
  type(complex_data_layout)                          :: ExK
  type(wavenumbers)                                  :: WaveNum
  logical                                            :: res
  real(WP)                                           :: delta,coefClark
  real(WP)                                           :: coef

  success = .false.
  delta = computeDelta(U)

  if (.not. initDataLayout("ExK", ExK,(U%nx/2)+1,U%ny,U%nz, &
     & U%Lx,U%Ly,U%Lz,nbcpus,spec_rank,alongZ)) then
     write(6,'(a,i0,a)')'[ERROR] initDataLayout for TYPE (COMPLEX_DATA_LAYOUT) : not enought memory!'
     return
  endif
  !Init des tableaux de travail
  if (.not. initWorkArray(U,12,temp)) then
    write(6,'(a,i0,a)')'[ERROR] initWorkArray : not enought memory!'
    return
  endif
  res = initWorkArray(U,3,Ubar)
  res = initWorkArray(U,3,SdZ)
  res = initWorkArray(U,3,Uhatbar)
  res = initWorkArray(U,3,L)
  res = initWorkArray(U,3,K)
  res = initWorkArray(U,3,N)
  res = initWorkArray(U,3,UbarZbar)
  res = initWorkArray(U,3,UbarZbarTilde)
  res = initWorkArray(U,3,K1)
  res = initWorkArray(U,3,M)
  res = initWorkArray(U,4,Scal)
  res = initWorkArray(ExK,3,K1K)
  res = initWorkArray(ExK,3,UbarZbarK)
  res = initWorkArray(ExK,3,UbarZbarKTilde)
  res = initWorkArray(ExK,3,K2K)
  res = initWorkArray(ExK,3,UK)
  res = initWorkArray(ExK,3,UKbar)
  res = initWorkArray(ExK,3,UKhatbar)
  res = initWorkArray(ExK,3,ScalK)

  !Calcul des nombres d'onde
  if (.not. initWN(WaveNum,spec_rank,U%nx,U%ny,U%nz) ) then
    write(6,'(a,i0,a)')'[ERROR] Cannot init WN'
    return
  endif
  if (.not.  computeWN(WaveNum,spec_rank,U%Lx,U%Ly,U%Lz,U%nx,U%ny,U%nz) ) then
    write(6,'(a,i0,a)')'[ERROR] Cannot compute WN'
    return
  endif

 !Passage dans l'espace spectrale
 call ftran(U,UK(1),res)
 call ftran(V,UK(2),res)
 call ftran(W,UK(3),res)
 call ftran(ScalIn,ScalK(1),res)

 !Filtrage de tous les termes par le filtre 'filtre' à la longueur de coupure nd*delta
 call computeFilter(WaveNum,nd*delta,UK(1),UKbar(1),filter)
 call computeFilter(WaveNum,nd*delta,UK(2),UKbar(2),filter)
 call computeFilter(WaveNum,nd*delta,UK(3),UKbar(3),filter)
 call computeFilter(WaveNum,nd*delta,ScalK(1),ScalK(2),filter)
 !Retour espace physique des quantitées filtrées
 call btran(UKbar(1),Ubar(1),res)   !Ubar
 call btran(UKbar(2),Ubar(2),res)   !Vbar
 call btran(UKbar(3),Ubar(3),res)  !Wbar
 call btran(ScalK(2),Scal(1),res)  !Zbar

  !Calcul du terme du gradient d/dxi( dUbar_i/dxj dZbar/dxj )
  !coef = 1.0_WP/computeDelta(U)**2.0
  coef = 1.0_WP/12.0_WP
  !call gradientSca(ExK,Ubar(1),UK(1),UK(2),UK(3),ScalK(2),ScalIn,WaveNum,res,coef,coef*(nd*delta)**2)
  call gradientSca(ExK,Ubar(1),UK(1),UK(2),UK(3),ScalK(2),ScalIn,WaveNum,res,1,coef*(nd*delta)**2)
  call btran(ExK,Te2,res)

  !Calcul du coef dynamique de Clark
  !calcul des Li
! DO NOT WORK  ????
!  call computeGermanoSca(WaveNum,Ubar(1),Ubar(2),Ubar(3),Scal(1),&
!                       &UKbar(1),UKbar(2),UKbar(3),ScalK(2),filter,delta,&
!                       &L(1),L(2),L(3))
!  if (spec_rank .eq. 0 ) write(6,'(a,1x,f10.5)')'[INFO] L1',L(1)%values(L(1)%xmin,L(1)%ymin,L(1)%zmin)
! DO NOT WORK  ????

  !filtre test pour le calcul des coefs
  !Filtrage de tous les termes par le filtre 'filtre' à la longueur de coupure nd*delta
  call computeFilter(WaveNum,2.0_WP*nd*delta,UKbar(1),UKhatbar(1),filter)
  call computeFilter(WaveNum,2.0_WP*nd*delta,UKbar(2),UKhatbar(2),filter)
  call computeFilter(WaveNum,2.0_WP*nd*delta,UKbar(3),UKhatbar(3),filter)
  call computeFilter(WaveNum,2.0_WP*nd*delta,ScalK(2),ScalK(3),filter)
  !Retour espace physique des quantitées filtrées
  call btran(UKhatbar(1),Uhatbar(1),res)
  call btran(UKhatbar(2),Uhatbar(2),res)
  call btran(UKhatbar(3),Uhatbar(3),res)
  call btran(ScalK(3),Scal(2),res)  !Zhatbar

  !calcul des Li
  UbarZbar(1)%values = Ubar(1)%values * Scal(1)%values
  UbarZbar(2)%values = Ubar(2)%values * Scal(1)%values
  UbarZbar(3)%values = Ubar(3)%values * Scal(1)%values
  call ftran(UbarZbar(1),UbarZbarK(1),res)
  call ftran(UbarZbar(2),UbarZbarK(2),res)
  call ftran(UbarZbar(3),UbarZbarK(3),res)
  call computeFilter(WaveNum,2.0_WP*nd*delta,UbarZbarK(1),UbarZbarKTilde(1),filter)
  call computeFilter(WaveNum,2.0_WP*nd*delta,UbarZbarK(2),UbarZbarKTilde(2),filter)
  call computeFilter(WaveNum,2.0_WP*nd*delta,UbarZbarK(3),UbarZbarKTilde(3),filter)
  call btran(UbarZbarKTilde(1),UbarZbarTilde(1),res)
  call btran(UbarZbarKTilde(2),UbarZbarTilde(2),res)
  call btran(UbarZbarKTilde(3),UbarZbarTilde(3),res)
  L(1)%values=UbarZbarTilde(1)%values - Scal(2)%values * Uhatbar(1)%values
  L(2)%values=UbarZbarTilde(2)%values - Scal(2)%values * Uhatbar(2)%values
  L(3)%values=UbarZbarTilde(3)%values - Scal(2)%values * Uhatbar(3)%values

  !calcul des Ki
  res = computeFieldGradient(WaveNum,Uhatbar(1),temp(1),temp(2),temp(3),nbcpus,spec_rank)
  res = computeFieldGradient(WaveNum,Uhatbar(2),temp(4),temp(5),temp(6),nbcpus,spec_rank)
  res = computeFieldGradient(WaveNum,Uhatbar(3),temp(7),temp(8),temp(9),nbcpus,spec_rank)
  res = computeFieldGradient(WaveNum,Scal(2),temp(10),temp(11),temp(12),nbcpus,spec_rank)
  K(1)%values = (2*nd*delta)**2/12.0_WP * (temp(1)%values * temp(10)%values + &
                                   &temp(2)%values * temp(11)%values + &
                                   &temp(3)%values * temp(12)%values)
  K(2)%values = (2*nd*delta)**2/12.0_WP * (temp(4)%values * temp(10)%values + &
                                   &temp(5)%values * temp(11)%values + &
                                   &temp(6)%values * temp(12)%values)
  K(3)%values = (2*nd*delta)**2/12.0_WP * (temp(7)%values * temp(10)%values + &
                                   &temp(8)%values * temp(11)%values + &
                                   &temp(9)%values * temp(12)%values)

  res = computeFieldGradient(WaveNum,Ubar(1),temp(1),temp(2),temp(3),nbcpus,spec_rank)
  res = computeFieldGradient(WaveNum,Ubar(2),temp(4),temp(5),temp(6),nbcpus,spec_rank)
  res = computeFieldGradient(WaveNum,Ubar(3),temp(7),temp(8),temp(9),nbcpus,spec_rank)
  res = computeFieldGradient(WaveNum,Scal(1),temp(10),temp(11),temp(12),nbcpus,spec_rank)
  K1(1)%values = (nd*delta)**2/12.0_WP * (temp(1)%values * temp(10)%values + &
                                   &temp(2)%values * temp(11)%values + &
                                   &temp(3)%values * temp(12)%values)
  K1(2)%values = (nd*delta)**2/12.0_WP * (temp(4)%values * temp(10)%values + &
                                   &temp(5)%values * temp(11)%values + &
                                   &temp(6)%values * temp(12)%values)
  K1(3)%values = (nd*delta)**2/12.0_WP * (temp(7)%values * temp(10)%values + &
                                   &temp(8)%values * temp(11)%values + &
                                   &temp(9)%values * temp(12)%values)
  call ftran(K1(1),K1K(1),res)
  call ftran(K1(2),K1K(2),res)
  call ftran(K1(3),K1K(3),res)
  call computeFilter(WaveNum,2.0_WP*nd*delta,K1K(1),K2K(1),filter)
  call computeFilter(WaveNum,2.0_WP*nd*delta,K1K(2),K2K(2),filter)
  call computeFilter(WaveNum,2.0_WP*nd*delta,K1K(3),K2K(3),filter)
  call btran(K2K(1),K1(1),res)
  call btran(K2K(2),K1(2),res)
  call btran(K2K(3),K1(3),res)

  K(1)%values = K(1)%values - K1(1)%values
  K(2)%values = K(2)%values - K1(2)%values
  K(3)%values = K(3)%values - K1(3)%values

  !calcul des Ni
  call computeSdZ(Ubar(1),UKbar(1),UKbar(2),UKbar(3),ScalK(2),&
                  &WaveNum,N(1),N(2),N(3),res)
  call ftran(N(1),K1K(1),res)
  call ftran(N(2),K1K(2),res)
  call ftran(N(3),K1K(3),res)
  call computeFilter(WaveNum,2.0_WP*nd*delta,K1K(1),K2K(1),filter)
  call computeFilter(WaveNum,2.0_WP*nd*delta,K1K(2),K2K(2),filter)
  call computeFilter(WaveNum,2.0_WP*nd*delta,K1K(3),K2K(3),filter)
  call btran(K2K(1),M(1),res)
  call btran(K2K(2),M(2),res)
  call btran(K2K(3),M(3),res)
  call computeSdZ(Uhatbar(1),UKhatbar(1),UKhatbar(2),UKhatbar(3),ScalK(3),&
                  &WaveNum,N(1),N(2),N(3),res)

  N(1)%values = (2.0_WP*nd*delta)**2*N(1)%values - (nd*delta)**2*M(1)%values
  N(2)%values = (2.0_WP*nd*delta)**2*N(2)%values - (nd*delta)**2*M(2)%values
  N(3)%values = (2.0_WP*nd*delta)**2*N(3)%values - (nd*delta)**2*M(3)%values

  temp(1)%values = (L(1)%values-K(1)%values)*N(1)%values + (L(2)%values-K(2)%values)*N(2)%values +&
                 & (L(3)%values-K(3)%values)*N(3)%values
  temp(2)%values = N(1)%values*N(1)%values + N(2)%values*N(2)%values + N(3)%values*N(3)%values

  CoefClark = computeFieldAvg(temp(1),spec_rank) / computeFieldAvg(temp(2),spec_rank)
  if (spec_rank .eq. 0) write (6,'(a,1x,f9.6)') '[INFO] CoefDynClark=',CoefClark

  !Smagorinsky
!  call smagorinskySca(ExK,Ubar(1),UKbar(1),UKbar(2),UKbar(3),ScalK(2),WaveNum,res,spec_rank,&
!                     &statcoeff=CoefClark*(nd*delta)**2)
!  call btran(ExK,Te1,res)
  call computeSdZ(Ubar(1),UKbar(1),UKbar(2),UKbar(3),ScalK(2),WaveNum,SdZ(1),SdZ(2),SdZ(3),res)
  SdZ(1)%values = (nd*delta)**2 * SdZ(1)%values
  SdZ(2)%values = (nd*delta)**2 * SdZ(2)%values
  SdZ(3)%values = (nd*delta)**2 * SdZ(3)%values
  if (.not. computeFieldDivergence(WaveNum,SdZ(1),SdZ(2),SdZ(3),Te1,&
     & nbcpus,spec_rank)) then
    write (6,'(a)') '[ERROR] debugcomputeMseV0, divergence failed'
    return
  endif

  !Désallocation des tableaux de travail
  res = deleteWorkArray(Ubar)
  res = deleteWorkArray(SdZ)
  res = deleteWorkArray(Uhatbar)
  res = deleteWorkArray(L)
  res = deleteWorkArray(K)
  res = deleteWorkArray(N)
  res = deleteWorkArray(UbarZbar)
  res = deleteWorkArray(UbarZbarTilde)
  res = deleteWorkArray(K1)
  res = deleteWorkArray(M)
  res = deleteWorkArray(Scal)
  res = deleteWorkArray(K1K)
  res = deleteWorkArray(UbarZbarK)
  res = deleteWorkArray(UbarZbarKTilde)
  res = deleteWorkArray(K2K)
  res = deleteWorkArray(UK)
  res = deleteWorkArray(UKbar)
  res = deleteWorkArray(UKhatbar)
  res = deleteWorkArray(ScalK)
  res = deleteWorkArray(temp)
  res = deleteWN(WaveNum)
  call deletedataLayout(ExK)

    success = .true.

 end function computeClarkTerms


!-------------------------------------------------------------------------------------------
!> Compute the SGS term and the viscous term of scalar equation and it dumps results in ASCII
!! files
!! @author Antoine Volllant
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in]    nbcpus          = number of process in the cpu pool
!!    @param[in]    U               = Component of velocity
!!    @param[in]    V               = Component of velocity
!!    @param[in]    W               = Component of velocity
!!    @param[in]    Scal            = Scalar to compute
!!    @param[in]    ite             = the number of file
!!    @return       res             = logical value
!-------------------------------------------------------------------------------------------
 function computeClarkFabreTerms(U,V,W,ScalIn,filter,nd,nbcpus,spec_rank,Te1,Te2) result(success)

  !I/O data
  type(real_data_layout),intent(in)        :: U,V,W,ScalIn
  type(real_data_layout),intent(inout)     :: Te1,Te2
  integer,intent(in)                       :: nbcpus,spec_rank
  integer,intent(in)                       :: filter,nd
  logical                                  :: success

  !Local data
  type(real_data_layout),dimension(:),allocatable    :: Ubar
  type(real_data_layout),dimension(:),allocatable    :: SdZ
  type(real_data_layout),dimension(:),allocatable    :: Uhatbar
  type(real_data_layout),dimension(:),allocatable    :: L
  type(real_data_layout),dimension(:),allocatable    :: UbarZbar
  type(real_data_layout),dimension(:),allocatable    :: UbarZbarTilde
  type(real_data_layout),dimension(:),allocatable    :: K
  type(real_data_layout),dimension(:),allocatable    :: N
  type(real_data_layout),dimension(:),allocatable    :: Scal
  type(complex_data_layout),dimension(:),allocatable :: UK
  type(complex_data_layout),dimension(:),allocatable :: UKbar
  type(complex_data_layout),dimension(:),allocatable :: UbarZbarK
  type(complex_data_layout),dimension(:),allocatable :: UbarZbarKTilde
  type(complex_data_layout),dimension(:),allocatable :: UKhatbar
  type(complex_data_layout),dimension(:),allocatable :: ScalK
  type(real_data_layout),dimension(:),allocatable    :: temp
  type(complex_data_layout)                          :: ExK
  type(wavenumbers)                                  :: WaveNum
  logical                                            :: res
  real(WP)                                           :: delta,coefFabre
  real(WP)                                           :: coef

  success = .false.
  delta = computeDelta(U)

  if (.not. initDataLayout("ExK", ExK,(U%nx/2)+1,U%ny,U%nz, &
     & U%Lx,U%Ly,U%Lz,nbcpus,spec_rank,alongZ)) then
     write(6,'(a,i0,a)')'[ERROR] initDataLayout for TYPE (COMPLEX_DATA_LAYOUT) : not enought memory!'
     return
  endif
  !Init des tableaux de travail
  if (.not. initWorkArray(U,12,temp)) then
    write(6,'(a,i0,a)')'[ERROR] initWorkArray : not enought memory!'
    return
  endif
  res = initWorkArray(U,3,Ubar)
  res = initWorkArray(U,3,Uhatbar)
  res = initWorkArray(U,3,L)
  res = initWorkArray(U,3,K)
  res = initWorkArray(U,3,N)
  res = initWorkArray(U,3,UbarZbar)
  res = initWorkArray(U,3,UbarZbarTilde)
  res = initWorkArray(U,4,Scal)
  res = initWorkArray(U,3,SdZ)
  res = initWorkArray(ExK,3,UbarZbarK)
  res = initWorkArray(ExK,3,UbarZbarKTilde)
  res = initWorkArray(ExK,3,UK)
  res = initWorkArray(ExK,3,UKbar)
  res = initWorkArray(ExK,3,UKhatbar)
  res = initWorkArray(ExK,3,ScalK)

  !Calcul des nombres d'onde
  if (.not. initWN(WaveNum,spec_rank,U%nx,U%ny,U%nz) ) then
    write(6,'(a,i0,a)')'[ERROR] Cannot init WN'
    return
  endif
  if (.not.  computeWN(WaveNum,spec_rank,U%Lx,U%Ly,U%Lz,U%nx,U%ny,U%nz) ) then
    write(6,'(a,i0,a)')'[ERROR] Cannot compute WN'
    return
  endif

  !Passage dans l'espace spectrale
  call ftran(U,UK(1),res)
  call ftran(V,UK(2),res)
  call ftran(W,UK(3),res)
  call ftran(ScalIn,ScalK(1),res)

  !Filtrage de tous les termes par le filtre 'filtre' à la longueur de coupure nd*delta
  call computeFilter(WaveNum,nd*delta,UK(1),UKbar(1),filter)
  call computeFilter(WaveNum,nd*delta,UK(2),UKbar(2),filter)
  call computeFilter(WaveNum,nd*delta,UK(3),UKbar(3),filter)
  call computeFilter(WaveNum,nd*delta,ScalK(1),ScalK(2),filter)
  !Retour espace physique des quantitées filtrées
  call btran(UKbar(1),Ubar(1),res)   !Ubar
  if (.not. res ) print *,'btran 1'
  call btran(UKbar(2),Ubar(2),res)   !Vbar
  if (.not. res ) print *,'btran 2'
  call btran(UKbar(3),Ubar(3),res)  !Wbar
  if (.not. res ) print *,'btran 3'
  call btran(ScalK(2),Scal(1),res)  !Zbar
  if (.not. res ) print *,'btran 4'

  !Calcul du terme du gradient d/dxi( dUbar_i/dxj dZbar/dxj )
  !coef = 1.0_WP/computeDelta(U)**2.0
  coef = 1.0_WP/12.0_WP
  !call gradientSca(ExK,Ubar(1),UK(1),UK(2),UK(3),ScalK(2),ScalIn,WaveNum,res,coef,coef*(nd*delta)**2)
  call gradientSca(ExK,Ubar(1),UK(1),UK(2),UK(3),ScalK(2),ScalIn,WaveNum,res,1,coef*(nd*delta)**2)
  call btran(ExK,Te2,res)

  !Calcul du coef dynamique de Y.Fabre
  !calcul des Li
! DO NOT WORK  ????
!  call computeGermanoSca(WaveNum,Ubar(1),Ubar(2),Ubar(3),Scal(1),&
!                       &UKbar(1),UKbar(2),UKbar(3),ScalK(2),filter,delta,&
!                       &L(1),L(2),L(3))
!  if (spec_rank .eq. 0 ) write(6,'(a,1x,f10.5)')'[INFO] L1',L(1)%values(L(1)%xmin,L(1)%ymin,L(1)%zmin)
! DO NOT WORK  ????

  !filtre test pour le calcul des coefs
  !Filtrage de tous les termes par le filtre 'filtre' à la longueur de coupure nd*delta
  call computeFilter(WaveNum,2.0_WP*nd*delta,UKbar(1),UKhatbar(1),filter)
  call computeFilter(WaveNum,2.0_WP*nd*delta,UKbar(2),UKhatbar(2),filter)
  call computeFilter(WaveNum,2.0_WP*nd*delta,UKbar(3),UKhatbar(3),filter)
  call computeFilter(WaveNum,2.0_WP*nd*delta,ScalK(2),ScalK(3),filter)
  !Retour espace physique des quantitées filtrées
  call btran(UKhatbar(1),Uhatbar(1),res)
  call btran(UKhatbar(2),Uhatbar(2),res)
  call btran(UKhatbar(3),Uhatbar(3),res)
  call btran(ScalK(3),Scal(2),res)  !Zhatbar

  !calcul des Li
  UbarZbar(1)%values = Ubar(1)%values * Scal(1)%values
  UbarZbar(2)%values = Ubar(2)%values * Scal(1)%values
  UbarZbar(3)%values = Ubar(3)%values * Scal(1)%values
  call ftran(UbarZbar(1),UbarZbarK(1),res)
  call ftran(UbarZbar(2),UbarZbarK(2),res)
  call ftran(UbarZbar(3),UbarZbarK(3),res)
  call computeFilter(WaveNum,2.0_WP*nd*delta,UbarZbarK(1),UbarZbarKTilde(1),filter)
  call computeFilter(WaveNum,2.0_WP*nd*delta,UbarZbarK(2),UbarZbarKTilde(2),filter)
  call computeFilter(WaveNum,2.0_WP*nd*delta,UbarZbarK(3),UbarZbarKTilde(3),filter)
  call btran(UbarZbarKTilde(1),UbarZbarTilde(1),res)
  call btran(UbarZbarKTilde(2),UbarZbarTilde(2),res)
  call btran(UbarZbarKTilde(3),UbarZbarTilde(3),res)
  L(1)%values=UbarZbarTilde(1)%values - Scal(2)%values * Uhatbar(1)%values
  L(2)%values=UbarZbarTilde(2)%values - Scal(2)%values * Uhatbar(2)%values
  L(3)%values=UbarZbarTilde(3)%values - Scal(2)%values * Uhatbar(3)%values

  !calcul des Ki = DeltaBarTilde^2 / 12 x dUiBarTilde/dxj x dZBarTilde/dxj
  res = computeFieldGradient(WaveNum,Uhatbar(1),temp(1),temp(2),temp(3),nbcpus,spec_rank)
  res = computeFieldGradient(WaveNum,Uhatbar(2),temp(4),temp(5),temp(6),nbcpus,spec_rank)
  res = computeFieldGradient(WaveNum,Uhatbar(3),temp(7),temp(8),temp(9),nbcpus,spec_rank)
  res = computeFieldGradient(WaveNum,Scal(2),temp(10),temp(11),temp(12),nbcpus,spec_rank)
  K(1)%values = (2*nd*delta)**2/12.0_WP * (temp(1)%values * temp(10)%values + &
                                   &temp(2)%values * temp(11)%values + &
                                   &temp(3)%values * temp(12)%values)
  K(2)%values = (2*nd*delta)**2/12.0_WP * (temp(4)%values * temp(10)%values + &
                                   &temp(5)%values * temp(11)%values + &
                                   &temp(6)%values * temp(12)%values)
  K(3)%values = (2*nd*delta)**2/12.0_WP * (temp(7)%values * temp(10)%values + &
                                   &temp(8)%values * temp(11)%values + &
                                   &temp(9)%values * temp(12)%values)

  !calcul des Ni = DeltaBarTilde^2 x |SbarTilde| x dZbarTilde / dxi
  !    computeSdZ(U,Uk,Vk,Wk,ScalK,ScalWN,outSdZ1,outSdZ2,outSdZ3,res)
  call computeSdZ(Uhatbar(1),UKhatbar(1),UKhatbar(2),UKhatbar(3),ScalK(3),&
                  &WaveNum,N(1),N(2),N(3),res)
  N(1)%values = (2*nd*delta)**2*N(1)%values
  N(2)%values = (2*nd*delta)**2*N(2)%values
  N(3)%values = (2*nd*delta)**2*N(3)%values

  temp(1)%values = (L(1)%values-K(1)%values)*N(1)%values + (L(2)%values-K(2)%values)*N(2)%values +&
                 & (L(3)%values-K(3)%values)*N(3)%values
  temp(2)%values = N(1)%values*N(1)%values + N(2)%values*N(2)%values + N(3)%values*N(3)%values

  CoefFabre = computeFieldAvg(temp(1),spec_rank) / computeFieldAvg(temp(2),spec_rank)
  if (spec_rank .eq. 0) write (6,'(a,1x,f9.6)') '[INFO] CoefFabre=',CoefFabre

!  call smagorinskySca(ExK,Ubar(1),UKbar(1),UKbar(2),UKbar(3),ScalK(2),WaveNum,res,spec_rank,&
!                     &statcoeff=CoefFabre*(nd*delta)**2)
!  call btran(ExK,Te1,res)
  call computeSdZ(Ubar(1),UKbar(1),UKbar(2),UKbar(3),ScalK(2),WaveNum,SdZ(1),SdZ(2),SdZ(3),res)
  SdZ(1)%values = (nd*delta)**2 * SdZ(1)%values
  SdZ(2)%values = (nd*delta)**2 * SdZ(2)%values
  SdZ(3)%values = (nd*delta)**2 * SdZ(3)%values
  if (.not. computeFieldDivergence(WaveNum,SdZ(1),SdZ(2),SdZ(3),Te1,&
     & nbcpus,spec_rank)) then
    write (6,'(a)') '[ERROR] debugcomputeMseV0, divergence failed'
    return
  endif

  !Désallocation des tableaux de travail
  res = deleteWorkArray(Ubar)
  res = deleteWorkArray(SdZ)
  res = deleteWorkArray(Uhatbar)
  res = deleteWorkArray(L)
  res = deleteWorkArray(K)
  res = deleteWorkArray(N)
  res = deleteWorkArray(UbarZbar)
  res = deleteWorkArray(UbarZbarTilde)
  res = deleteWorkArray(Scal)
  res = deleteWorkArray(UbarZbarK)
  res = deleteWorkArray(UbarZbarKTilde)
  res = deleteWorkArray(UK)
  res = deleteWorkArray(UKbar)
  res = deleteWorkArray(UKhatbar)
  res = deleteWorkArray(ScalK)
  res = deleteWorkArray(temp)
  res = deleteWN(WaveNum)
  call deletedataLayout(ExK)

    success = .true.

 end function computeClarkFabreTerms



 !-------------------------------------------------------------------------------------------
 !> Compute the SGS term and the viscous term of scalar equation and it dumps results in ASCII
 !! files 
 !! @author Antoine Volllant 
 !!    @param[in]    spec_rank       = mpi rank inside spectral communicator
 !!    @param[in]    nbcpus          = number of process in the cpu pool 
 !!    @param[in]    U               = Component of velocity in physical space
 !!    @param[in]    V               = Component of velocity in physical space
 !!    @param[in]    W               = Component of velocity in physical space
 !!    @param[in]    Scal            = Scalar in physical space 
 !!    @param[in]    Exk             = field in spectral space 
 !!    @param[in]    ScalBar        = Scalar filtered in physical space 
 !!    @param[in]    ite             = the number of file
 !!    @return       res             = logical value 
 !-------------------------------------------------------------------------------------------
  function computeAprioriTest(iFic,U,V,W,Scal,ScalBar,Exk,spec_rank,nd,filter,phi,&
                             & dTidxiMod,nameMod,computeModBiObj,dumpModOpt,Timod) result(success)

 
   !I/O data
   type(real_data_layout),intent(in)              :: U, V ,W, Scal
   type(complex_data_layout),intent(in)           :: Exk 
   type(real_data_layout),intent(inout)           :: ScalBar
   type(real_data_layout),dimension(:),intent(in),optional  :: Timod
   integer,intent(in)                             :: spec_rank
   integer,intent(in)                             :: iFic,filter
   integer,intent(in)                             :: nd
   type(real_data_layout),dimension(:),intent(in) :: phi
   type(real_data_layout),intent(in)              :: dTidxiMod
   character(len=*),intent(in)                    :: nameMod
   logical,intent(in)                             :: computeModBiObj,dumpModOpt
   logical                                        :: success
 
   !Local data
   type(real_data_layout),dimension(:),allocatable  :: TiDNS
   real(WP)                           :: delta
   integer                            :: inc,sizeOfPhi
   integer                            :: nbcpus 
   integer                            :: ires
   logical                            :: res
   character(len=64)                  :: filtre,file_name,write_alf
   character(len=64)                  :: fmtwrite1,fmtwrite2
   character(len=64)                  :: write_iFic 
   character(len=64)                  :: nameModShort 
   integer                            :: file_id,i,binTot
   integer,dimension(:),allocatable   :: bin,typeDisc 
   logical                            :: existe
   real(WP)                           :: eps = 1e-10_WP
   real(WP)                           :: time1,time2
   real(WP) ::  moyDissModOpt,moyDissMod,moyDissOpt,moyDissExact
   real(WP) ::  errTotSgs,errIrrSgs,errIrrDiss,errQuadDiss,errTotDiss
   real(WP) ::  alf,moyDissModOpt1,moyDissModOpt2,moyZbar
   real(WP) ::  errQuadDiss1,errQuadDiss2,errIrrSgs1,errIrrSgs2
   real(WP),dimension(:,:),allocatable :: modBiObj1
   real(WP),dimension(:,:),allocatable :: modBiObj2
   real(WP),dimension(:,:),allocatable :: modOptDiss
   real(WP),dimension(:,:),allocatable :: modOptSgs
   real(WP),dimension(:,:),allocatable :: modOptZbar
   type(real_data_layout)              :: dTidxiDNS
   type(real_data_layout)              :: angle 
   type(real_data_layout)              :: norm 
   type(complex_data_layout)           :: dealiasingK 
   type(real_data_layout)              :: dTidxiModOpt
   type(real_data_layout)              :: fieldDissOpt
   type(real_data_layout)              :: fieldDissModOpt
   type(real_data_layout)              :: fieldDissMod
   type(real_data_layout)              :: fieldDissExact
   type(real_data_layout)              :: dTidxiBiObj1
   type(real_data_layout)              :: dTidxiBiObj2
   type(real_data_layout)              :: dissBiObj1
   type(real_data_layout)              :: dissBiObj2
   type(wavenumbers)                   :: WaveNum
   real(WP)                            :: minNorm,moyNorm,maxNorm,varNorm,kurtNorm
   real(WP)                            :: minAngle,moyAngle,maxAngle,varAngle,kurtAngle 
 
   success = .false.
   nbcpus = getnbcpus()
   delta = computeDelta(U)

   sizeOfPhi=size(phi,1)
   allocate(typeDisc(sizeofPhi),stat=ires)
   if (ires .ne. 0) return
   allocate(bin(sizeofPhi),stat=ires)
   if (ires .ne. 0) return
   bin = 300
   typeDisc = 1
   binTot = bin(1)**sizeofPhi
   write(write_iFic,'(i5.5)') iFic
   write_iFic = trim(adjustl(write_iFic))
   write(filtre,'(i2.2)') nd
   filtre = trim(adjustl(filtre))
   write(fmtwrite1,'(a,i0,a)') '(',sizeOfPhi+3,'(f15.8,1x))'
   fmtwrite1=trim(adjustl(fmtwrite1))
   write(fmtwrite2,'(a,i0,a)') '(',sizeOfPhi+2,'(f15.8,1x))'
   fmtwrite2=trim(adjustl(fmtwrite2))
   nameModShort=trim(adjustl(nameMod))
   !if (spec_rank .eq. 0) then
   !  write(6,'(a,i0)')'[INFO test apriori] sizeofPhi   = ',sizeOfPhi
   !  write(6,'(a,i0)')'[INFO test apriori] binTot      = ',binTot
   !  write(6,'(a,a)')'[INFO test apriori] write_iFic   = ',write_iFic
   !  write(6,'(a,a)')'[INFO test apriori] filtre       = ',filtre
   !  write(6,'(a,a)')'[INFO test apriori] fmtwrite1    = ',fmtwrite1
   !  write(6,'(a,a)')'[INFO test apriori] fmtwrite2    = ',fmtwrite2
   !  write(6,'(a,a)')'[INFO test apriori] nameModShort = ',nameModShort
   !endif
 
   !Init des tableau de travail
   if (.not. copyStructOnly(U,dTidxiDNS) .or. &
      &.not. copyStructOnly(Exk,dealiasingK) .or. &
      &.not. copyStructOnly(U,dTidxiModOpt) .or. &
      &.not. copyStructOnly(U,fieldDissOpt) .or. &
      &.not. copyStructOnly(U,fieldDissModOpt) .or. &
      &.not. copyStructOnly(U,fieldDissMod) .or. &
      &.not. copyStructOnly(U,fieldDissExact)) then
     write(6,'(a,i0,a)')'[ERROR] Not enought memory 1'
     return
   endif
   allocate(modOptDiss(sizeofPhi+1,binTot),stat=ires)
   if (ires .ne. 0) then
     write(6,'(a,i0,a)')'[ERROR] not enought memory 2' 
     return
   endif
   allocate(modOptSgs(sizeofPhi+1,binTot),stat=ires)
   if (ires .ne. 0) then
     write(6,'(a,i0,a)')'[ERROR] not enought memory 3' 
     return
   endif
   allocate(modOptZbar(sizeofPhi+1,binTot),stat=ires)
   if (ires .ne. 0) then
     write(6,'(a,i0,a)')'[ERROR] not enought memory 4' 
     return
   endif

   fieldDissOpt%name   = trim(nameModShort)//'_fieldDissOpt_nd_'//trim(filtre)//'_'
   fieldDissMod%name   = trim(nameModShort)//'_fieldDissMod_nd_'//trim(filtre)//'_'
   fieldDissModOpt%name= trim(nameModShort)//'_fieldDissModOpt_nd_'//trim(filtre)//'_'

   !Calcul des nombres d'onde
   if (.not. initWN(WaveNum,spec_rank,U%nx,U%ny,U%nz) .or. &
      &.not.  computeWN(WaveNum,spec_rank,U%Lx,U%Ly,U%Lz,U%nx,U%ny,U%nz) ) then
     write(6,'(a,i0,a)')'[ERROR] Cannot compute WN'
     return
   endif

   !Filtrage des quantitées
   !Calcul du terme sous maille exact
   if ( present(Timod)) then
     if (.not. initWorkArray(U,3,TiDNS)) return
     if (.not. copyStructOnly(U,angle)) return
     if (.not. copyStructOnly(U,norm)) return
     if (.not. computedTidxi(nbcpus,nd,filter,spec_rank, Scal, U, V ,W , dTidxiDNS,TiDNS) ) then
       write(6,'(a,i0,a)')'[ERROR] in computedTidxi' 
       return
     endif
     if (.not. computeVectorSimilarity(TiDNS,Timod,norm,angle,.true.,ific) ) then
       write(6,'(a,i0,a)')'[ERROR] in computeVectorSimilarity' 
       return
     endif
     minAngle  = computeFieldMin(angle,spec_rank)
     maxAngle  = computeFieldMax(angle,spec_rank)
     moyAngle  = computeFieldAvg(angle,spec_rank)
     varAngle  = computeFieldVar(angle,spec_rank)
     kurtAngle = computeFieldKurt(angle,spec_rank)
     minNorm   = computeFieldMin(norm,spec_rank)
     maxNorm   = computeFieldMax(norm,spec_rank)
     moyNorm   = computeFieldAvg(norm,spec_rank)
     varNorm   = computeFieldVar(norm,spec_rank)
     kurtNorm  = computeFieldKurt(norm,spec_rank)
     res = deleteWorkArray(TiDNS)
     call deletedatalayout(angle)
     call deletedatalayout(norm)
     if (spec_rank .eq. 0) then
       write(6,'(5(a,1x,f10.4,1x))')'[INFO ANGLE TiDNS/TiMod] min=',minAngle,'moy=',moyAngle,'max=',maxAngle,'var=',varAngle,'kurt=',kurtAngle
       write(6,'(5(a,1x,f10.4,1x))')'[INFO NORM  TiDNS/TiMod] min=',minNorm,'moy=',moyNorm,'max=',maxNorm,'var=',varNorm,'kurt=',kurtNorm 
       !dump angle's stat
       file_id = iopen()
       !write(file_name,'(a)') trim(nameModShort)//'_angle_'//trim(filtre)//'_'//trim(write_iFic)//'.table'
       write(file_name,'(a)') trim(nameModShort)//'_angle_'//trim(filtre)//'.table'
       inquire( file=trim(adjustl(file_name)), exist=existe)
       open(unit=file_id, file=trim(adjustl(file_name)), form="FORMATTED",position="append")
       if ( .not. existe) then
         write(file_id, '(a)') '#fichier min moy max var kurt' 
       endif
       write(file_id,'(i0,5(1x,f10.4))') iFic,minAngle,moyAngle,maxAngle,varAngle,kurtAngle 
       close(file_id)
       !dump norms's stat
       !write(file_name,'(a)') trim(nameModShort)//'_norm_'//trim(filtre)//'_'//trim(write_iFic)//'.table'
       write(file_name,'(a)') trim(nameModShort)//'_norm_'//trim(filtre)//'.table'
       inquire( file=trim(adjustl(file_name)), exist=existe)
       open(unit=file_id, file=trim(adjustl(file_name)), form="FORMATTED",position="append")
       if ( .not. existe) then
         write(file_id, '(a)') '#fichier min moy max var kurt' 
       endif
       write(file_id,'(i0,5(1x,f10.4))') iFic,minNorm,moyNorm,maxNorm,varNorm,kurtNorm
       close(file_id)
       file_id=iclose(file_id)
     endif
   else
     if (.not. computedTidxi(nbcpus,nd,filter,spec_rank, Scal, U, V ,W , dTidxiDNS) ) then
       write(6,'(a,i0,a)')'[ERROR] in computedTidxi' 
       return
     endif
   endif
   !Calcul de la moyenne de Zbar 
   moyZbar = computeFieldAvg( ScalBar , spec_rank )
   !Champ de dissipation
   fieldDissMod%values = ScalBar%values * dTidxiMod%values
   fieldDissExact%values = ScalBar%values * dTidxiDNS%values
   !Calcul de l'estimateur optimale du flux scalaire sous maille 
   CALL cpu_time(time1)
   if (.not. computeEONVar( spec_rank,phi,dTidxiDNS,&
                      & typeDisc,bin,&
                      & estimOpt = modOptSgs,&
                      & fieldOpt = dTidxiModOpt) ) return 
   CALL cpu_time(time2)
   if (spec_rank .eq. 0) then
     write(6,'(a,1x,f6.3,1x,a)')'[INFO test apriori] estimateur optimale dTidxiModOpt',time2-time1,'secondes'
   endif
   !Calcul de l'estimateur optimal du scalaire 
   CALL cpu_time(time1)
   if (.not. computeEONVar( spec_rank,phi,ScalBar,&
                      & typeDisc,bin,&
                      & estimOpt = modOptZbar)) return
   CALL cpu_time(time2)
   if (spec_rank .eq. 0) then
     write(6,'(a,1x,f6.3,1x,a)')'[INFO test apriori] estimateur optimale modOptZbar',time2-time1,'secondes'
   endif
   !Calcul de l'estimateur optimal de la dissipation 
   CALL cpu_time(time1)
   if (.not. computeEONVar( spec_rank,phi,fieldDissExact,&
                      & typeDisc,bin,&
                      & estimOpt = modOptDiss,&
                      & fieldOpt = fieldDissOpt) ) return 
   CALL cpu_time(time2)
   if (spec_rank .eq. 0) then
     write(6,'(a,1x,f6.3,1x,a)')'[INFO test apriori] estimateur optimale fieldDissOpt',time2-time1,'secondes'
   endif
   !Champ de dissipation
   fieldDissModOpt%values = ScalBar%values * dTidxiModOpt%values
   !Dealiasing
   CALL cpu_time(time1)
   call ftran(fieldDissModOpt,dealiasingK,res)
   call dealiasIsotrope(dealiasingK,WaveNum)
   call btran(dealiasingK,fieldDissModOpt,res)
   call ftran(fieldDissMod,dealiasingK,res)
   call dealiasIsotrope(dealiasingK,WaveNum)
   call btran(dealiasingK,fieldDissMod,res)
   CALL cpu_time(time2)
   if (spec_rank .eq. 0) then
     write(6,'(a,1x,f6.3,1x,a)')'[INFO test apriori] Dealiasing',time2-time1,'secondes'
   endif
   !Ecriture des spectres de dissipation / du scalaire
   CALL cpu_time(time1)
   if (.not. compute_spectrum(fieldDissOpt, iFic ,  spec_rank) ) return
   if (.not. compute_spectrum(fieldDissModOpt, iFic , spec_rank) ) return
   if (.not. compute_spectrum(fieldDissMod, iFic , spec_rank) ) return
   CALL cpu_time(time2)
   if (spec_rank .eq. 0) then
     write(6,'(a,1x,f6.3,1x,a)')'[INFO test apriori] Spectrum',time2-time1,'secondes'
   endif
   !Moyenne dissipation
   CALL cpu_time(time1)
   moyDissModOpt = computeFieldAvg( fieldDissModOpt , spec_rank ) 
   moyDissMod    = computeFieldAvg( fieldDissMod , spec_rank ) 
   moyDissOpt    = computeFieldAvg( fieldDissOpt , spec_rank ) 
   moyDissExact  = computeFieldAvg( fieldDissExact , spec_rank ) 
   CALL cpu_time(time2)
   if (spec_rank .eq. 0) then
     write(6,'(a,1x,f6.3,1x,a)')'[INFO test apriori] Averages',time2-time1,'secondes'
   endif
   !Erreurs
   CALL cpu_time(time1)
   errTotSgs   = computeMSE(dTidxiDNS,dTidxiMod,spec_rank,.true.)
   errIrrSgs   = computeMSE(dTidxiDNS,dTidxiModOpt,spec_rank,.true.)
   errIrrDiss  = computeMSE(fieldDissExact,fieldDissOpt,spec_rank,.true.)
   errQuadDiss = computeMSE(fieldDissExact,fieldDissModOpt,spec_rank,.true.)
   errTotDiss  = computeMSE(fieldDissExact,fieldDissMod,spec_rank,.true.)
   CALL cpu_time(time2)
   if (spec_rank .eq. 0) then
     write(6,'(a)')'[INFO test apriori] Start mse...'
   endif
   !Dump datas
   CALL cpu_time(time1)
   if (spec_rank .eq. 0) then
     write(6,'(a)')'[INFO test apriori] Start dumping datas...'
   endif
   if (spec_rank .eq. 0) then
     if (dumpModOpt) then
       !ecrit les modelles optimaux
       file_id=iopen()
       write(file_name,'(a)') trim(nameModShort)//'_modOpt_'//trim(filtre)//'_'//trim(write_iFic)//'.table'
       inquire( file=trim(adjustl(file_name)), exist=existe)
       open(unit=file_id, file=trim(adjustl(file_name)), form="FORMATTED",position="append")
       if ( .not. existe) then
         write(file_id, '(a)') '#phi(s) dTidxiOpt ZbardTidxiOpt ZbarOpt' 
       endif
       do i =1,binTot
         write(file_id,fmtwrite1) modOptSgs(:,i),modOptDiss(sizeOfPhi+1,i),modOptZbar(sizeOfPhi+1,i)
         enddo
       close(file_id)
       file_id=iclose(file_id)
     endif

     !ecrit les resultats d erreurs 
     file_id=iopen()
     write(file_name,'(a)') trim(nameModShort)//'_aprioriTestMod_nd_'//trim(filtre)//'.table'
     inquire( file=trim(adjustl(file_name)), exist=existe)
     open(unit=file_id, file=trim(adjustl(file_name)), form="FORMATTED",position="append")
     if ( .not. existe) then
       write(file_id, '(a)') '#numFic -moyDissModOpt -moyDissMod -moyDissOpt -moyDissExact &
                            &errTotSgs errIrrSgs errIrrDiss errQuadDiss errTotDiss errRelDissModOpt errRelDissMod'
     endif
     write(file_id,'(i0,1x,11(f10.6,1x))') iFic,-moyDissModOpt,-moyDissMod,-moyDissOpt,-moyDissExact,errTotSgs,&
                        &errIrrSgs,errIrrDiss,errQuadDiss,errTotDiss,(moyDissModOpt-moyDissExact)/moyDissExact,&
                        &(moyDissMod-moyDissExact)/moyDissExact
     close(file_id)
     file_id=iclose(file_id)
   end if
   CALL cpu_time(time2)
   if (spec_rank .eq. 0) then
     write(6,'(a,1x,f6.3,1x,a)')'[INFO test apriori] Dumping datas',time2-time1,'secondes'
   endif
   !Release memory
   call deletedatalayout(dTidxiModOpt)
   call deletedatalayout(fieldDissOpt)
   call deletedatalayout(fieldDissModOpt)
   call deletedatalayout(fieldDissMod)

   if ( computeModBiObj ) then
     if (.not. copyStructOnly(U,dTidxiBiObj1) .or. &
        &.not. copyStructOnly(U,dTidxiBiObj2) .or. &
        &.not. copyStructOnly(U,dissBiObj1) .or. &
        &.not. copyStructOnly(U,dissBiObj2)) then
       write(6,'(a,i0,a)')'[ERROR] Not enought memory 1'
       return
     endif
     allocate(modBiObj1(sizeofPhi+1,binTot),stat=ires)
     if (ires .ne. 0) then
       write(6,'(a,i0,a)')'[ERROR] not enought memory 2' 
       return
     endif
     allocate(modBiObj2(sizeofPhi+1,binTot),stat=ires)
     if (ires .ne. 0) then
       write(6,'(a,i0,a)')'[ERROR] not enought memory 3' 
       return
     endif
     !Calcul des modèles bi optim
     do i=1,sizeOfPhi
       modBiObj1(i,:)=modOptSgs(i,:)
       modBiObj2(i,:)=modOptSgs(i,:)
     enddo
     do inc = 1,6
       !Creation des modèles Bi-optimaux
       alf = 0.20_WP*(inc-1)
       if (spec_rank .eq. 0) then
         write(6,'(a,1x,f6.3,1x,a)')'[INFO test apriori] Start alpha',alf,'...'
       endif
       write(write_alf,'(f6.3)') alf
       write_alf = trim(adjustl(write_alf))
       if (spec_rank .eq. 0) then
         write(6,'(a)')'[INFO test apriori]          Building bi-objective model...'
       endif
       do i = 1,binTot
         if ( modOptZbar(sizeOfPhi+1,i) .gt. eps ) then 
           modBiObj1(sizeOfPhi+1,i) = (1.0_WP - alf ) * modOptSgs(sizeOfPhi+1,i) + &
                       & alf * modOptDiss(sizeOfPhi+1,i) / modOptZbar(sizeOfPhi+1,i)
         else
           modBiObj1(sizeOfPhi+1,i) = (1.0_WP - alf ) * modOptSgs(sizeOfPhi+1,i)
         endif
       enddo
       modBiObj2(sizeOfPhi+1,:) = (1.0_WP - alf ) * modOptSgs(sizeOfPhi+1,:) + &
                                       & alf * modOptDiss(sizeOfPhi+1,:) / moyZbar
       !reconstruction des champs
       CALL cpu_time(time1)
       if (.not. buildEOFieldFromEO(phi,modBiObj1,dTidxiBiObj1)) return
       if (.not. buildEOFieldFromEO(phi,modBiObj2,dTidxiBiObj2)) return
       CALL cpu_time(time2)
       if (spec_rank .eq. 0) then
         write(6,'(a,1x,f8.3,1x,a)')'[INFO test apriori]          Building fields...',time2-time1,'secondes'
       endif
       !Calcul des dissipations
       dissBiObj1%values = ScalBar%values * dTidxiBiObj1%values
       dissBiObj2%values = ScalBar%values * dTidxiBiObj2%values
       !Dealiasing
       if (spec_rank .eq. 0) then
         write(6,'(a)')'[INFO test apriori]          Starting dealiasing...'
       endif
       call ftran(dissBiObj1,dealiasingK,res)
       call dealiasIsotrope(dealiasingK,WaveNum)
       call btran(dealiasingK,dissBiObj1,res)
       call ftran(dissBiObj2,dealiasingK,res)
       call dealiasIsotrope(dealiasingK,WaveNum)
       call btran(dealiasingK,dissBiObj2,res)
       !Moyenne des dissipations
       moyDissModOpt1 = computeFieldAvg( dissBiObj1 , spec_rank ) 
       moyDissModOpt2 = computeFieldAvg( dissBiObj2 , spec_rank ) 
       !Erreurs
       errQuadDiss1 = computeMSE(fieldDissExact,dissBiObj1,spec_rank,.true.)
       errQuadDiss2 = computeMSE(fieldDissExact,dissBiObj2,spec_rank,.true.)
       errIrrSgs1   = computeMSE(dTidxiDNS,dTidxiBiObj1,spec_rank,.true.)
       errIrrSgs2   = computeMSE(dTidxiDNS,dTidxiBiObj2,spec_rank,.true.)
       !Etude spectrale
       dissBiObj1%name=trim(nameModShort)//'_dissBiObj1_nd_'//trim(filtre)//'_alf_'//trim(write_alf)//'_'
       dissBiObj2%name=trim(nameModShort)//'_dissBiObj2_nd_'//trim(filtre)//'_alf_'//trim(write_alf)//'_'
       if (.not. compute_spectrum(dissBiObj1, iFic, spec_rank) ) return
       if (.not. compute_spectrum(dissBiObj2, iFic, spec_rank) ) return
       !Dump datas
       if (spec_rank .eq. 0) then
         file_id=iopen()
         write(file_name,'(a)') trim(nameModShort)//'_aprioriTestBiObj_nd_'//trim(filtre)//&
                                &'_fic_'//trim(write_iFic)//'.table'
         inquire( file=trim(adjustl(file_name)), exist=existe)
         open(unit=file_id, file=trim(adjustl(file_name)), form="FORMATTED",position="append")
         if ( .not. existe) then
           write(file_id, '(a)') '#alf -moyDissMod1 -moyDissMod2 errQuadDiss1 errQuadDiss2 & 
                                 & errIrrSgsBiObj1 errIrrSgsBiObj2 errRelDissMod1 errRelDissMod2' 
         endif
         write(file_id,'(f5.2,1x,8(f10.6,1x))') alf,-moyDissModOpt1,-moyDissModOpt2,errQuadDiss1,errQuadDiss2, &
                                           & errIrrSgs1,errIrrSgs2,(moyDissModOpt1-moyDissExact)/moyDissExact, & 
                                           & (moyDissModOpt2-moyDissExact)/moyDissExact
         close(file_id)
         !Ecriture des modèles bi_obj
         
         write(file_name,'(a)') trim(nameModShort)//'_modBiObj_nd_'//trim(filtre)//'_alf_'//trim(write_alf)&
                                     &//'_fic_'//trim(write_iFic)//'.table'
         inquire( file=trim(adjustl(file_name)), exist=existe)
         open(unit=file_id, file=trim(adjustl(file_name)), form="FORMATTED",position="append")
         if ( .not. existe) then
           write(file_id,'(a)') '#phiGrad modBiObj1 modBiObj2'
         endif
         do i=1,binTot
           write(file_id,fmtwrite2) modBiObj1(:,i),modBiObj2(sizeOfPhi+1,i)
         enddo
         close(file_id)
         file_id=iclose(file_id)
       end if
     enddo
     !Release memory
     deallocate(modBiObj1)
     deallocate(modBiObj2)
     call deletedatalayout(dTidxiBiObj1)
     call deletedatalayout(dTidxiBiObj2)
     call deletedatalayout(dissBiObj1)
     call deletedatalayout(dissBiObj2)
   endif

   !Release memory
   call deletedatalayout(dTidxiDNS)
   call deletedatalayout(fieldDissExact)
   deallocate(modOptDiss)
   deallocate(modOptZbar)
   call deletedatalayout(dealiasingK)
   deallocate(modOptSgs)
   deallocate(typeDisc)
   deallocate(bin)

 
   success = .true.
 
  end function computeAprioriTest



!-------------------------------------------------------------------------------------------
!> Compute the SGS term and the viscous term of scalar equation and it dumps results in ASCII
!! files
!! @author Antoine Volllant
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in]    nbcpus          = number of process in the cpu pool
!!    @param[in]    U               = Component of velocity
!!    @param[in]    V               = Component of velocity
!!    @param[in]    W               = Component of velocity
!!    @param[in]    Scal            = Scalar to compute
!!    @param[in]    ite             = the number of file
!!    @return       res             = logical value
!-------------------------------------------------------------------------------------------
 function apriori_ClarkFabre(iFic,U,V,W,Scal,spec_rank,filtMin,filtMax,filter,dTidxi) result(success)

  !I/O data
  type(real_data_layout),intent(in)      :: U, V ,W, Scal
  integer,intent(in)                     :: spec_rank
  integer,intent(in)                     :: iFic,filter
  integer,intent(in)                     :: filtMin,filtMax
  logical                                :: success
  type(real_data_layout),optional,intent(inout) :: dTidxi

  !Local data
  real(WP)                          :: delta,coeff
  integer                           :: nbcpus
  integer                           :: nd
  logical                           :: res,boucle
  type(wavenumbers)                 :: WaveNum
  type(complex_data_layout)         :: ExK
  type(real_data_layout),dimension(:),allocatable :: SdZ
  type(complex_data_layout),dimension(:),allocatable :: VelKBar
  type(complex_data_layout),dimension(:),allocatable :: VelK
  type(real_data_layout)            :: ScalBar
  type(complex_data_layout)         :: ScalK
  type(complex_data_layout)         :: ScalKBar
  type(complex_data_layout)         :: fieldPhiModK
  type(real_data_layout)         :: fieldPhiMod1
  type(real_data_layout)         :: fieldPhiMod2
  type(real_data_layout)         :: fieldPhiMod3
  type(real_data_layout)            :: dTidxiMod

  success = .false.
  delta = computeDelta(U)
  nbcpus = getnbcpus()

  !Init des tableau de travail
  if (.not. initDataLayout("ExK", ExK,(U%nx/2)+1,U%ny,U%nz, &
     & U%Lx,U%Ly,U%Lz,nbcpus,spec_rank,alongZ)) then
    write(6,'(a,i0,a)')'[ERROR] Not able to create Exk'
    return
  endif
  if (.not. initWorkArray(ExK,3,VelKBar) .or. &
     &.not. initWorkArray(ExK,3,VelK) .or. &
     &.not. initWorkArray(U,3,SdZ) .or. &
     &.not. copyStructOnly(Exk,ScalK) .or. &
     &.not. copyStructOnly(Exk,ScalKBar) .or. &
     &.not. copyStructOnly(Exk,fieldPhiModK) .or. &
     &.not. copyStructOnly(U,fieldPhiMod1) .or. &
     &.not. copyStructOnly(U,fieldPhiMod2) .or. &
     &.not. copyStructOnly(U,fieldPhiMod3) .or. &
     &.not. copyStructOnly(U,ScalBar) .or. &
     &.not. copyStructOnly(U,dTidxiMod) ) then
    write(6,'(a,i0,a)')'[ERROR] Not enought memory 1'
    return
  endif

  !Calcul des nombres d'onde
  if (.not. initWN(WaveNum,spec_rank,U%nx,U%ny,U%nz) .or. &
     &.not. computeWN(WaveNum,spec_rank,U%Lx,U%Ly,U%Lz,U%nx,U%ny,U%nz) ) then
    write(6,'(a,i0,a)')'[ERROR] Cannot compute WN'
    return
  endif

  nd = filtMin
  boucle = .true.
  do while (boucle)
    !passage dans l'espace spectrale
    call ftran(U,VelK(1),res)
    if (.not. res) return
    call ftran(V,VelK(2),res)
    if (.not. res) return
    call ftran(W,VelK(3),res)
    if (.not. res) return
    call ftran(Scal,ScalK,res)
    if (.not. res) return
    !Filtrage de tous les termes par le filtre 'filtre' à la longueur de coupure nd*delta
    call computeFilter(WaveNum,nd*delta,VelK(1),VelKBar(1),filter)
    call computeFilter(WaveNum,nd*delta,VelK(2),VelKBar(2),filter)
    call computeFilter(WaveNum,nd*delta,VelK(3),VelKBar(3),filter)
    call computeFilter(WaveNum,nd*delta,ScalK,ScalKBar,filter)
    !Retour dans l'espace physique de z
    call btran(ScalKBar,ScalBar,res)
    if (.not. res) return
    !Calcul des terms des modèles sous maille
    !Gradient
    call gradientSca(fieldPhiModK,U,VelKBar(1),VelKBar(2),VelKBar(3),ScalKBar,Scal,WaveNum,res,nd)
    if (.not. res) return
    call btran(fieldPhiModK,fieldPhiMod1,res)
    if (.not. res) return
    !Smagorinsky
    call computeSdZ(U,VelKBar(1),VelKBar(2),VelKBar(3),ScalKBar,WaveNum,SdZ(1),SdZ(2),SdZ(3),res)
    if (.not. res) return
    if (.not. coefDynClarkFabreSca(U,V,W,VelKBar(1),VelKBar(2),VelKBar(3),Scal,ScalK,WaveNum,filter,nd,spec_rank,&
       & CoeffDynScalar=coeff) ) then
      write (6,'(a)') '[ERROR] in DynClarkFabreSca'
      return
    endif
    SdZ(1)%values = (nd*delta)**2 * SdZ(1)%values
    SdZ(2)%values = (nd*delta)**2 * SdZ(2)%values
    SdZ(3)%values = (nd*delta)**2 * SdZ(3)%values
    if (.not. computeFieldDivergence(WaveNum,SdZ(1),SdZ(2),SdZ(3),fieldPhiMod2,nbcpus,spec_rank)) return
    SdZ(1)%values = coeff * SdZ(1)%values
    SdZ(2)%values = coeff * SdZ(2)%values
    SdZ(3)%values = coeff * SdZ(3)%values
    if (.not. computeFieldDivergence(WaveNum,SdZ(1),SdZ(2),SdZ(3),fieldPhiMod3,nbcpus,spec_rank)) return
    !--SGS Clark--
    dTidxiMod%values = fieldPhiMod3%values + fieldPhiMod1%values
    if (present (dTidxi) ) then
      if ( samelayout(dTidxiMod,dTidxi) ) dTidxi%values = dTidxiMod%values
    else
      if (.not. computeAprioriTest(iFic,U,V,W,Scal,ScalBar,Exk,spec_rank,nd,filter,&
                               & (/fieldPhiMod1,fieldPhiMod2/),dTidxiMod,'ClarkFabre',.true.,.true.) ) then
        write (6,'(a)') '[ERROR] computeAprioriTest failed'
        return
      endif
    endif
    !--Incrementation boucle--
    nd = nd + 2
    if ( nd .gt. filtMax) then
      boucle = .false.
    endif
  enddo

  !Release memory
  call deletedataLayout(ExK)
  res = deleteWorkArray(SdZ)
  res = deleteWorkArray(VelKBar)
  res = deleteWorkArray(VelK)
  res = deleteWN(WaveNum)
  call deletedatalayout(ScalBar)
  call deletedatalayout(ScalK)
  call deletedatalayout(ScalKBar)
  call deletedatalayout(fieldPhiModK)
  call deletedatalayout(fieldPhiMod1)
  call deletedatalayout(fieldPhiMod2)
  call deletedatalayout(fieldPhiMod3)
  call deletedatalayout(dTidxiMod )

  success = .true.

 end function apriori_ClarkFabre


!-------------------------------------------------------------------------------------------
!> Compute the SGS term and the viscous term of scalar equation and it dumps results in ASCII
!! files
!! @author Antoine Volllant
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in]    nbcpus          = number of process in the cpu pool
!!    @param[in]    U               = Component of velocity
!!    @param[in]    V               = Component of velocity
!!    @param[in]    W               = Component of velocity
!!    @param[in]    Scal            = Scalar to compute
!!    @param[in]    ite             = the number of file
!!    @return       res             = logical value
!-------------------------------------------------------------------------------------------
 function apriori_Clark(iFic,U,V,W,Scal,spec_rank,filtMin,filtMax,filter,dTidxi) result(success)

  !I/O data
  type(real_data_layout),intent(in)      :: U, V ,W, Scal
  integer,intent(in)                     :: spec_rank
  integer,intent(in)                     :: iFic,filter
  integer,intent(in)                     :: filtMin,filtMax
  logical                                :: success
  type(real_data_layout),optional,intent(inout) :: dTidxi

  !Local data
  real(WP)                          :: delta,coeff
  integer                           :: nbcpus
  integer                           :: nd
  logical                           :: res,boucle
  type(wavenumbers)                 :: WaveNum
  type(complex_data_layout)         :: ExK
  type(real_data_layout),dimension(:),allocatable :: SdZ
  type(complex_data_layout),dimension(:),allocatable :: VelKBar
  type(complex_data_layout),dimension(:),allocatable :: VelK
  type(real_data_layout)            :: ScalBar
  type(complex_data_layout)         :: ScalK
  type(complex_data_layout)         :: ScalKBar
  type(complex_data_layout)         :: fieldPhiModK
  type(real_data_layout)         :: fieldPhiMod1
  type(real_data_layout)         :: fieldPhiMod2
  type(real_data_layout)         :: fieldPhiMod3
  type(real_data_layout)            :: dTidxiMod

  success = .false.
  delta = computeDelta(U)
  nbcpus = getnbcpus()

  !Init des tableau de travail
  if (.not. initDataLayout("ExK", ExK,(U%nx/2)+1,U%ny,U%nz, &
     & U%Lx,U%Ly,U%Lz,nbcpus,spec_rank,alongZ)) then
    write(6,'(a,i0,a)')'[ERROR] Not able to create Exk'
    return
  endif
  if (.not. initWorkArray(ExK,3,VelKBar) .or. &
     &.not. initWorkArray(ExK,3,VelK) .or. &
     &.not. initWorkArray(U,3,SdZ) .or. &
     &.not. copyStructOnly(Exk,ScalK) .or. &
     &.not. copyStructOnly(Exk,ScalKBar) .or. &
     &.not. copyStructOnly(Exk,fieldPhiModK) .or. &
     &.not. copyStructOnly(U,fieldPhiMod1) .or. &
     &.not. copyStructOnly(U,fieldPhiMod2) .or. &
     &.not. copyStructOnly(U,fieldPhiMod3) .or. &
     &.not. copyStructOnly(U,ScalBar) .or. &
     &.not. copyStructOnly(U,dTidxiMod) ) then
    write(6,'(a,i0,a)')'[ERROR] Not enought memory 1'
    return
  endif

  !Calcul des nombres d'onde
  if (.not. initWN(WaveNum,spec_rank,U%nx,U%ny,U%nz) .or. &
     &.not.  computeWN(WaveNum,spec_rank,U%Lx,U%Ly,U%Lz,U%nx,U%ny,U%nz) ) then
    write(6,'(a,i0,a)')'[ERROR] Cannot compute WN'
    return
  endif

  nd = filtMin
  boucle = .true.
  do while (boucle)
    !passage dans l'espace spectrale
    call ftran(U,VelK(1),res)
    if (.not. res) return
    call ftran(V,VelK(2),res)
    if (.not. res) return
    call ftran(W,VelK(3),res)
    if (.not. res) return
    call ftran(Scal,ScalK,res)
    if (.not. res) return
    !Filtrage de tous les termes par le filtre 'filtre' à la longueur de coupure nd*delta
    call computeFilter(WaveNum,nd*delta,VelK(1),VelKBar(1),filter)
    call computeFilter(WaveNum,nd*delta,VelK(2),VelKBar(2),filter)
    call computeFilter(WaveNum,nd*delta,VelK(3),VelKBar(3),filter)
    call computeFilter(WaveNum,nd*delta,ScalK,ScalKBar,filter)
    !Retour dans l'espace physique de z
    call btran(ScalKBar,ScalBar,res)
    if (.not. res) return
    !Calcul des terms des modèles sous maille
    !Gradient
    call gradientSca(fieldPhiModK,U,VelKBar(1),VelKBar(2),VelKBar(3),ScalKBar,Scal,WaveNum,res,nd)
    if (.not. res) return
    call btran(fieldPhiModK,fieldPhiMod1,res)
    if (.not. res) return
    !Smagorinsky
    call computeSdZ(U,VelKBar(1),VelKBar(2),VelKBar(3),ScalKBar,WaveNum,SdZ(1),SdZ(2),SdZ(3),res)
    if (.not. res) return
    if (.not. coefDynClarkSca(U,V,W,VelKBar(1),VelKBar(2),VelKBar(3),Scal,ScalK,WaveNum,filter,nd,spec_rank,&
       & CoeffDynScalar=coeff) ) then
      write (6,'(a)') '[ERROR] in DynClarkFabreSca'
      return
    endif
    SdZ(1)%values = (nd*delta)**2 * SdZ(1)%values
    SdZ(2)%values = (nd*delta)**2 * SdZ(2)%values
    SdZ(3)%values = (nd*delta)**2 * SdZ(3)%values
    if (.not. computeFieldDivergence(WaveNum,SdZ(1),SdZ(2),SdZ(3),fieldPhiMod2,nbcpus,spec_rank)) return
    SdZ(1)%values = coeff * SdZ(1)%values
    SdZ(2)%values = coeff * SdZ(2)%values
    SdZ(3)%values = coeff * SdZ(3)%values
    if (.not. computeFieldDivergence(WaveNum,SdZ(1),SdZ(2),SdZ(3),fieldPhiMod3,nbcpus,spec_rank)) return
    !--SGS Clark--
    dTidxiMod%values = fieldPhiMod3%values + fieldPhiMod1%values
    if (present (dTidxi) ) then
      if ( samelayout(dTidxiMod,dTidxi) ) dTidxi%values = dTidxiMod%values
    else
      if (.not. computeAprioriTest(iFic,U,V,W,Scal,ScalBar,Exk,spec_rank,nd,filter,&
                               & (/fieldPhiMod1,fieldPhiMod2/),dTidxiMod,'Clark',.true.,.true.) ) then
        write (6,'(a)') '[ERROR] computeAprioriTest failed'
        return
      endif
    endif
    !--Incrementation boucle--
    nd = nd + 2
    if ( nd .gt. filtMax) then
      boucle = .false.
    endif
  enddo

  !Release memory
  call deletedataLayout(ExK)
  res = deleteWorkArray(SdZ)
  res = deleteWorkArray(VelKBar)
  res = deleteWorkArray(VelK)
  res = deleteWN(WaveNum)
  call deletedatalayout(ScalBar)
  call deletedatalayout(ScalK)
  call deletedatalayout(ScalKBar)
  call deletedatalayout(fieldPhiModK)
  call deletedatalayout(fieldPhiMod1)
  call deletedatalayout(fieldPhiMod2)
  call deletedatalayout(fieldPhiMod3)
  call deletedatalayout(dTidxiMod )

  success = .true.

 end function apriori_Clark

!-------------------------------------------------------------------------------------------
!> Compute the SGS term and the viscous term of scalar equation and it dumps results in ASCII
!! files
!! @author Antoine Volllant
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in]    nbcpus          = number of process in the cpu pool
!!    @param[in]    U               = Component of velocity
!!    @param[in]    V               = Component of velocity
!!    @param[in]    W               = Component of velocity
!!    @param[in]    Scal            = Scalar to compute
!!    @param[in]    ite             = the number of file
!!    @return       res             = logical value
!-------------------------------------------------------------------------------------------
 function apriori_Smagorinsky(iFic,U,V,W,Scal,spec_rank,filtMin,filtMax,filter,dTidxi) result(success)

  !I/O data
  type(real_data_layout),intent(in)      :: U, V ,W, Scal
  integer,intent(in)                     :: spec_rank
  integer,intent(in)                     :: iFic,filter
  integer,intent(in)                     :: filtMin,filtMax
  logical                                :: success
  type(real_data_layout),optional,intent(inout) :: dTidxi

  !Local data
  real(WP)                          :: delta
  integer                           :: nbcpus
  integer                           :: nd
  logical                           :: res,boucle
  type(wavenumbers)                 :: WaveNum
  type(complex_data_layout)         :: ExK
  type(real_data_layout),dimension(:),allocatable :: SdZ
  type(complex_data_layout),dimension(:),allocatable :: VelKBar
  type(complex_data_layout),dimension(:),allocatable :: VelK
  type(real_data_layout)            :: ScalBar
  type(complex_data_layout)         :: ScalK
  type(complex_data_layout)         :: ScalKBar
  type(complex_data_layout)         :: fieldPhiModK
  type(real_data_layout)            :: fieldPhiMod
  type(real_data_layout)            :: dTidxiMod

  success = .false.
  delta = computeDelta(U)
  nbcpus = getnbcpus()

  !Init des tableau de travail
  if (.not. initDataLayout("ExK", ExK,(U%nx/2)+1,U%ny,U%nz, &
     & U%Lx,U%Ly,U%Lz,nbcpus,spec_rank,alongZ)) then
    write(6,'(a,i0,a)')'[ERROR] Not able to create Exk'
    return
  endif
  if (.not. initWorkArray(ExK,3,VelKBar) .or. &
     &.not. initWorkArray(ExK,3,VelK) .or. &
     &.not. initWorkArray(U,3,SdZ) .or. &
     &.not. copyStructOnly(Exk,ScalK) .or. &
     &.not. copyStructOnly(Exk,ScalKBar) .or. &
     &.not. copyStructOnly(Exk,fieldPhiModK) .or. &
     &.not. copyStructOnly(U,ScalBar) .or. &
     &.not. copyStructOnly(U,fieldPhiMod) .or. &
     &.not. copyStructOnly(U,dTidxiMod) ) then
    write(6,'(a,i0,a)')'[ERROR] Not enought memory 1'
    return
  endif

  !Calcul des nombres d'onde
  if (.not. initWN(WaveNum,spec_rank,U%nx,U%ny,U%nz) .or. &
     &.not.  computeWN(WaveNum,spec_rank,U%Lx,U%Ly,U%Lz,U%nx,U%ny,U%nz) ) then
    write(6,'(a,i0,a)')'[ERROR] Cannot compute WN'
    return
  endif

  nd = filtMin
  boucle = .true.
  do while (boucle)
    !passage dans l'espace spectrale
    call ftran(U,VelK(1),res)
    if (.not. res) return
    call ftran(V,VelK(2),res)
    if (.not. res) return
    call ftran(W,VelK(3),res)
    if (.not. res) return
    call ftran(Scal,ScalK,res)
    if (.not. res) return
    !Filtrage de tous les termes par le filtre 'filtre' à la longueur de coupure nd*delta
    call computeFilter(WaveNum,nd*delta,VelK(1),VelKBar(1),filter)
    call computeFilter(WaveNum,nd*delta,VelK(2),VelKBar(2),filter)
    call computeFilter(WaveNum,nd*delta,VelK(3),VelKBar(3),filter)
    call computeFilter(WaveNum,nd*delta,ScalK,ScalKBar,filter)
    !Retour dans l'espace physique de z
    call btran(ScalKBar,ScalBar,res)       !Zbar
    if (.not. res) return
    call computeSdZ(U,VelKBar(1),VelKBar(2),VelKBar(3),ScalKBar,WaveNum,SdZ(1),SdZ(2),SdZ(3),res)
    SdZ(1)%values = (nd*delta)**2 * SdZ(1)%values
    SdZ(2)%values = (nd*delta)**2 * SdZ(2)%values
    SdZ(3)%values = (nd*delta)**2 * SdZ(3)%values
    if (.not. computeFieldDivergence(WaveNum,SdZ(1),SdZ(2),SdZ(3),fieldPhiMod,nbcpus,spec_rank)) return
    SdZ(1)%values = 0.18_WP**2/0.6_WP * SdZ(1)%values
    SdZ(2)%values = 0.18_WP**2/0.6_WP * SdZ(2)%values
    SdZ(3)%values = 0.18_WP**2/0.6_WP * SdZ(3)%values
    if (.not. computeFieldDivergence(WaveNum,SdZ(1),SdZ(2),SdZ(3),dTidxiMod,nbcpus,spec_rank)) return
    if (present (dTidxi) ) then
      if ( samelayout(dTidxiMod,dTidxi) ) dTidxi%values = dTidxiMod%values
    else
      if (.not. computeAprioriTest(iFic,U,V,W,Scal,ScalBar,Exk,spec_rank,nd,filter,(/fieldPhiMod/),&
                               & dTidxiMod,'Smag',.true.,.true.) ) then
        write (6,'(a)') '[ERROR] computeAprioriTest failed'
        return
      endif
    endif
    nd = nd + 2
    if ( nd .gt. filtMax) then
      boucle = .false.
    endif
  enddo

  !Release memory
  call deletedataLayout(ExK)
  res = deleteWorkArray(VelKBar)
  res = deleteWorkArray(SdZ)
  res = deleteWorkArray(VelK)
  res = deleteWN(WaveNum)
  call deletedatalayout(ScalBar)
  call deletedatalayout(ScalK)
  call deletedatalayout(ScalKBar)
  call deletedatalayout(fieldPhiModK)
  call deletedatalayout(fieldPhiMod)
  call deletedatalayout(dTidxiMod )

  success = .true.

 end function apriori_Smagorinsky

!-------------------------------------------------------------------------------------------
!> Compute the SGS term and the viscous term of scalar equation and it dumps results in ASCII
!! files
!! @author Antoine Volllant
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in]    nbcpus          = number of process in the cpu pool
!!    @param[in]    U               = Component of velocity
!!    @param[in]    V               = Component of velocity
!!    @param[in]    W               = Component of velocity
!!    @param[in]    Scal            = Scalar to compute
!!    @param[in]    ite             = the number of file
!!    @return       res             = logical value
!-------------------------------------------------------------------------------------------
 function apriori_SmagorinskyDyn(iFic,U,V,W,Scal,spec_rank,filtMin,filtMax,filter,dTidxi) result(success)

  !I/O data
  type(real_data_layout),intent(in)      :: U, V ,W, Scal
  integer,intent(in)                     :: spec_rank
  integer,intent(in)                     :: iFic,filter
  integer,intent(in)                     :: filtMin,filtMax
  logical                                :: success
  type(real_data_layout),optional,intent(inout) :: dTidxi

  !Local data
  real(WP)                          :: delta,coeff
  integer                           :: nbcpus
  integer                           :: nd
  logical                           :: res,boucle
  type(wavenumbers)                 :: WaveNum
  type(complex_data_layout)         :: ExK
  type(real_data_layout),dimension(:),allocatable :: SdZ
  type(complex_data_layout),dimension(:),allocatable :: VelKBar
  type(complex_data_layout),dimension(:),allocatable :: VelK
  type(real_data_layout)            :: ScalBar
  type(complex_data_layout)         :: ScalK
  type(complex_data_layout)         :: ScalKBar
  type(complex_data_layout)         :: fieldPhiModK
  type(real_data_layout)            :: fieldPhiMod
  type(real_data_layout)            :: dTidxiMod

  success = .false.
  delta = computeDelta(U)
  nbcpus = getnbcpus()

  !Init des tableau de travail
  if (.not. initDataLayout("ExK", ExK,(U%nx/2)+1,U%ny,U%nz, &
     & U%Lx,U%Ly,U%Lz,nbcpus,spec_rank,alongZ)) then
    write(6,'(a,i0,a)')'[ERROR] Not able to create Exk'
    return
  endif
  if (.not. initWorkArray(ExK,3,VelKBar) .or. &
     &.not. initWorkArray(ExK,3,VelK) .or. &
     &.not. initWorkArray(U,3,SdZ) .or. &
     &.not. copyStructOnly(Exk,ScalK) .or. &
     &.not. copyStructOnly(Exk,ScalKBar) .or. &
     &.not. copyStructOnly(Exk,fieldPhiModK) .or. &
     &.not. copyStructOnly(U,ScalBar) .or. &
     &.not. copyStructOnly(U,fieldPhiMod) .or. &
     &.not. copyStructOnly(U,dTidxiMod) ) then
    write(6,'(a,i0,a)')'[ERROR] Not enought memory 1'
    return
  endif

  !Calcul des nombres d'onde
  if (.not. initWN(WaveNum,spec_rank,U%nx,U%ny,U%nz) .or. &
     &.not.  computeWN(WaveNum,spec_rank,U%Lx,U%Ly,U%Lz,U%nx,U%ny,U%nz) ) then
    write(6,'(a,i0,a)')'[ERROR] Cannot compute WN'
    return
  endif

  nd = filtMin
  boucle = .true.
  do while (boucle)
    !passage dans l'espace spectrale
    call ftran(U,VelK(1),res)
    if (.not. res) return
    call ftran(V,VelK(2),res)
    if (.not. res) return
    call ftran(W,VelK(3),res)
    if (.not. res) return
    call ftran(Scal,ScalK,res)
    if (.not. res) return
    !Filtrage de tous les termes par le filtre 'filtre' à la longueur de coupure nd*delta
    call computeFilter(WaveNum,nd*delta,VelK(1),VelKBar(1),filter)
    call computeFilter(WaveNum,nd*delta,VelK(2),VelKBar(2),filter)
    call computeFilter(WaveNum,nd*delta,VelK(3),VelKBar(3),filter)
    call computeFilter(WaveNum,nd*delta,ScalK,ScalKBar,filter)
    !Retour dans l'espace physique de z
    call btran(ScalKBar,ScalBar,res)       !Zbar
    if (.not. res) return
    call computeSdZ(U,VelKBar(1),VelKBar(2),VelKBar(3),ScalKBar,WaveNum,SdZ(1),SdZ(2),SdZ(3),res)
    if (.not. coefDynSmagSca(U,V,W,VelKBar(1),VelKBar(2),VelKBar(3),Scal,ScalK,WaveNum,filter,nd,spec_rank,&
       & CoeffDynScalar=coeff) ) then
      write (6,'(a)') '[ERROR] in DynClarkFabreSca'
      res = .false.
    endif
    SdZ(1)%values = (nd*delta)**2 * SdZ(1)%values
    SdZ(2)%values = (nd*delta)**2 * SdZ(2)%values
    SdZ(3)%values = (nd*delta)**2 * SdZ(3)%values
    if (.not. computeFieldDivergence(WaveNum,SdZ(1),SdZ(2),SdZ(3),fieldPhiMod,nbcpus,spec_rank)) return
    SdZ(1)%values = coeff * SdZ(1)%values
    SdZ(2)%values = coeff * SdZ(2)%values
    SdZ(3)%values = coeff * SdZ(3)%values
    if (.not. computeFieldDivergence(WaveNum,SdZ(1),SdZ(2),SdZ(3),dTidxiMod,nbcpus,spec_rank)) return
    if (present (dTidxi) ) then
      if ( samelayout(dTidxiMod,dTidxi) ) dTidxi%values = dTidxiMod%values
    else
      if (.not. computeAprioriTest(iFic,U,V,W,Scal,ScalBar,Exk,spec_rank,nd,filter,(/fieldPhiMod/),&
                               & dTidxiMod,'SmagDyn',.true.,.true.) ) then
        write (6,'(a)') '[ERROR] computeAprioriTest failed'
        return
      endif
    endif
    nd = nd + 2
    if ( nd .gt. filtMax) then
      boucle = .false.
    endif
  enddo

  !Release memory
  call deletedataLayout(ExK)
  res = deleteWorkArray(VelKBar)
  res = deleteWorkArray(SdZ)
  res = deleteWorkArray(VelK)
  res = deleteWN(WaveNum)
  call deletedatalayout(ScalBar)
  call deletedatalayout(ScalK)
  call deletedatalayout(ScalKBar)
  call deletedatalayout(fieldPhiModK)
  call deletedatalayout(fieldPhiMod)
  call deletedatalayout(dTidxiMod )

  success = .true.

 end function apriori_SmagorinskyDyn

!-------------------------------------------------------------------------------------------
!> Compute the SGS term and the viscous term of scalar equation and it dumps results in ASCII
!! files
!! @author Antoine Volllant
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in]    nbcpus          = number of process in the cpu pool
!!    @param[in]    U               = Component of velocity
!!    @param[in]    V               = Component of velocity
!!    @param[in]    W               = Component of velocity
!!    @param[in]    Scal            = Scalar to compute
!!    @param[in]    ite             = the number of file
!!    @return       res             = logical value
!-------------------------------------------------------------------------------------------
 function apriori_Gradient(iFic,U,V,W,Scal,spec_rank,filtMin,filtMax,filter,dTidxi) result(success)

  !I/O data
  type(real_data_layout),intent(in)      :: U, V ,W, Scal
  integer,intent(in)                     :: spec_rank
  integer,intent(in)                     :: iFic,filter
  integer,intent(in)                     :: filtMin,filtMax
  logical                                :: success
  type(real_data_layout),optional,intent(inout) :: dTidxi

  !Local data
  real(WP)                          :: delta
  integer                           :: nbcpus
  integer                           :: nd
  logical                           :: res,boucle
  type(wavenumbers)                 :: WaveNum
  type(complex_data_layout)         :: ExK
  type(complex_data_layout),dimension(:),allocatable :: VelKBar
  type(complex_data_layout),dimension(:),allocatable :: VelK
  type(real_data_layout)            :: ScalBar
  type(complex_data_layout)         :: ScalK
  type(complex_data_layout)         :: ScalKBar
  type(complex_data_layout)         :: fieldPhiModK
  type(real_data_layout)            :: fieldPhiMod
  type(real_data_layout)            :: dTidxiMod

  success = .false.
  delta = computeDelta(U)
  nbcpus = getnbcpus()

  !Init des tableau de travail
  if (.not. initDataLayout("ExK", ExK,(U%nx/2)+1,U%ny,U%nz, &
     & U%Lx,U%Ly,U%Lz,nbcpus,spec_rank,alongZ)) then
    write(6,'(a,i0,a)')'[ERROR] Not able to create Exk'
    return
  endif
  if (.not. initWorkArray(ExK,3,VelKBar) .or. &
     &.not. initWorkArray(ExK,3,VelK) .or. &
     &.not. copyStructOnly(Exk,ScalK) .or. &
     &.not. copyStructOnly(Exk,ScalKBar) .or. &
     &.not. copyStructOnly(Exk,fieldPhiModK) .or. &
     &.not. copyStructOnly(U,ScalBar) .or. &
     &.not. copyStructOnly(U,fieldPhiMod) .or. &
     &.not. copyStructOnly(U,dTidxiMod) ) then
    write(6,'(a,i0,a)')'[ERROR] Not enought memory 1'
    return
  endif

  !Calcul des nombres d'onde
  if (.not. initWN(WaveNum,spec_rank,U%nx,U%ny,U%nz) .or. &
     &.not.  computeWN(WaveNum,spec_rank,U%Lx,U%Ly,U%Lz,U%nx,U%ny,U%nz) ) then
    write(6,'(a,i0,a)')'[ERROR] Cannot compute WN'
    return
  endif

  nd = filtMin
  boucle = .true.
  do while (boucle)
    !passage dans l'espace spectrale
    call ftran(U,VelK(1),res)
    if (.not. res) return
    call ftran(V,VelK(2),res)
    if (.not. res) return
    call ftran(W,VelK(3),res)
    if (.not. res) return
    call ftran(Scal,ScalK,res)
    if (.not. res) return
    !Filtrage de tous les termes par le filtre 'filtre' à la longueur de coupure nd*delta
    call computeFilter(WaveNum,nd*delta,VelK(1),VelKBar(1),filter)
    call computeFilter(WaveNum,nd*delta,VelK(2),VelKBar(2),filter)
    call computeFilter(WaveNum,nd*delta,VelK(3),VelKBar(3),filter)
    call computeFilter(WaveNum,nd*delta,ScalK,ScalKBar,filter)
    !Retour dans l'espace physique de z
    call btran(ScalKBar,ScalBar,res)       !Zbar
    if (.not. res) return
    !Calcul des terms des modèles sous maille
    call GradientSca(fieldPhiModK,U,VelKBar(1),VelKBar(2),VelKBar(3),ScalKBar,Scal,WaveNum,res,nd,1.0_WP)
    if (.not. res) return
    call btran(fieldPhiModK,fieldPhiMod,res)
    if (.not. res) return
    !Calcul du modèle du gradient
    dTidxiMod%values = 1.0_WP/12.0_WP * fieldPhiMod%values
    if (present (dTidxi) ) then
      if ( samelayout(dTidxiMod,dTidxi) ) dTidxi%values = dTidxiMod%values
    else
      if (.not. computeAprioriTest(iFic,U,V,W,Scal,ScalBar,Exk,spec_rank,nd,filter,(/fieldPhiMod/),&
                               & dTidxiMod,'Grad',.true.,.true.) ) then
        write (6,'(a)') '[ERROR] computeAprioriTest failed'
        return
      endif
    endif
    nd = nd + 2
    if ( nd .gt. filtMax) then
      boucle = .false.
    endif
  enddo

  !Release memory
  call deletedataLayout(ExK)
  res = deleteWorkArray(VelKBar)
  res = deleteWorkArray(VelK)
  res = deleteWN(WaveNum)
  call deletedatalayout(ScalBar)
  call deletedatalayout(ScalK)
  call deletedatalayout(ScalKBar)
  call deletedatalayout(fieldPhiModK)
  call deletedatalayout(fieldPhiMod)
  call deletedatalayout(dTidxiMod )

  success = .true.

 end function apriori_Gradient


!-------------------------------------------------------------------------------------------
!> Compute the SGS term and the viscous term of scalar equation and it dumps results in ASCII
!! files
!! @author Antoine Volllant
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in]    nbcpus          = number of process in the cpu pool
!!    @param[in]    U               = Component of velocity
!!    @param[in]    V               = Component of velocity
!!    @param[in]    W               = Component of velocity
!!    @param[in]    Scal            = Scalar to compute
!!    @param[in]    ite             = the number of file
!!    @return       res             = logical value
!-------------------------------------------------------------------------------------------
 function apriori_RGM_analyt(iFic,U,V,W,Scal,spec_rank,filtMin,filtMax,filter,dTidxi) result(success)

  !I/O data
  type(real_data_layout),intent(in)      :: U, V ,W, Scal
  integer,intent(in)                     :: spec_rank
  integer,intent(in)                     :: iFic,filter
  integer,intent(in)                     :: filtMin,filtMax
  logical                                :: success
  type(real_data_layout),optional,intent(inout) :: dTidxi

  !Local data
  real(WP)                          :: delta
  integer                           :: nbcpus
  integer                           :: nd
  logical                           :: res,boucle
  type(wavenumbers)                 :: WaveNum
  type(complex_data_layout)         :: ExK
  type(complex_data_layout),dimension(:),allocatable :: VelKBar
  type(complex_data_layout),dimension(:),allocatable :: VelK
  type(real_data_layout)            :: ScalBar
  type(complex_data_layout)         :: ScalK
  type(complex_data_layout)         :: ScalKBar
  type(complex_data_layout)         :: fieldPhiModK
  type(real_data_layout)            :: fieldPhiMod
  type(real_data_layout)            :: dTidxiMod
  type(real_data_layout),dimension(:),allocatable :: TiMod

  success = .false.
  delta = computeDelta(U)
  nbcpus = getnbcpus()

  !Init des tableau de travail
  if (.not. initDataLayout("ExK", ExK,(U%nx/2)+1,U%ny,U%nz, &
     & U%Lx,U%Ly,U%Lz,nbcpus,spec_rank,alongZ)) then
    write(6,'(a,i0,a)')'[ERROR] Not able to create Exk'
    return
  endif
  if (.not. initWorkArray(ExK,3,VelKBar) .or. &
     &.not. initWorkArray(ExK,3,VelK) .or. &
     &.not. initWorkArray(U,3,TiMod) .or. &
     &.not. copyStructOnly(Exk,ScalK) .or. &
     &.not. copyStructOnly(Exk,ScalKBar) .or. &
     &.not. copyStructOnly(Exk,fieldPhiModK) .or. &
     &.not. copyStructOnly(U,ScalBar) .or. &
     &.not. copyStructOnly(U,fieldPhiMod) .or. &
     &.not. copyStructOnly(U,dTidxiMod) ) then
    write(6,'(a,i0,a)')'[ERROR] Not enought memory 1'
    return
  endif

  !Calcul des nombres d'onde
  if (.not. initWN(WaveNum,spec_rank,U%nx,U%ny,U%nz) .or. &
     &.not.  computeWN(WaveNum,spec_rank,U%Lx,U%Ly,U%Lz,U%nx,U%ny,U%nz) ) then
    write(6,'(a,i0,a)')'[ERROR] Cannot compute WN'
    return
  endif

  nd = filtMin
  boucle = .true.
  do while (boucle)
    !passage dans l'espace spectrale
    call ftran(U,VelK(1),res)
    if (.not. res) return
    call ftran(V,VelK(2),res)
    if (.not. res) return
    call ftran(W,VelK(3),res)
    if (.not. res) return
    call ftran(Scal,ScalK,res)
    if (.not. res) return
    !Filtrage de tous les termes par le filtre 'filtre' à la longueur de coupure nd*delta
    call computeFilter(WaveNum,nd*delta,VelK(1),VelKBar(1),filter)
    call computeFilter(WaveNum,nd*delta,VelK(2),VelKBar(2),filter)
    call computeFilter(WaveNum,nd*delta,VelK(3),VelKBar(3),filter)
    call computeFilter(WaveNum,nd*delta,ScalK,ScalKBar,filter)
    !Retour dans l'espace physique de z
    call btran(ScalKBar,ScalBar,res)       !Zbar
    if (.not. res) return
    !Calcul des terms des modèles sous maille
    if (present (dTidxi) ) then
      call RGMSca_analyt(fieldPhiModK,U,VelKBar(1),VelKBar(2),VelKBar(3),ScalKBar,WaveNum,res,spec_rank,1.0_WP)
    else
      call RGMSca_analyt(fieldPhiModK,U,VelKBar(1),VelKBar(2),VelKBar(3),ScalKBar,WaveNum,res,spec_rank,1.0_WP,Timod)
    endif
    if (.not. res) return
    call btran(fieldPhiModK,fieldPhiMod,res)
    if (.not. res) return
    !Calcul du modèle du gradient
    dTidxiMod%values = 1.0_WP/12.0_WP * fieldPhiMod%values
    if (present (dTidxi) ) then
      if ( samelayout(dTidxiMod,dTidxi) ) dTidxi%values = dTidxiMod%values
    else
      if (.not. computeAprioriTest(iFic,U,V,W,Scal,ScalBar,Exk,spec_rank,nd,filter,(/fieldPhiMod/),&
                               & dTidxiMod,'RGM_mod',.false.,.false.,Timod) ) then
        write (6,'(a)') '[ERROR] computeAprioriTest failed'
        return
      endif
    endif
    nd = nd + 2
    if ( nd .gt. filtMax) then
      boucle = .false.
    endif
  enddo

  !Release memory
  call deletedataLayout(ExK)
  res = deleteWorkArray(VelKBar)
  res = deleteWorkArray(VelK)
  res = deleteWorkArray(TiMod)
  res = deleteWN(WaveNum)
  call deletedatalayout(ScalBar)
  call deletedatalayout(ScalK)
  call deletedatalayout(ScalKBar)
  call deletedatalayout(fieldPhiModK)
  call deletedatalayout(fieldPhiMod)
  call deletedatalayout(dTidxiMod )

  success = .true.

 end function apriori_RGM_analyt

!-------------------------------------------------------------------------------------------
!> Compute the SGS term and the viscous term of scalar equation and it dumps results in ASCII
!! files
!! @author Antoine Volllant
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in]    nbcpus          = number of process in the cpu pool
!!    @param[in]    U               = Component of velocity
!!    @param[in]    V               = Component of velocity
!!    @param[in]    W               = Component of velocity
!!    @param[in]    Scal            = Scalar to compute
!!    @param[in]    ite             = the number of file
!!    @return       res             = logical value
!-------------------------------------------------------------------------------------------
 function apriori_RGM(iFic,U,V,W,Scal,spec_rank,filtMin,filtMax,filter) result(success)

  !I/O data
  type(real_data_layout),intent(in)      :: U, V ,W, Scal
  integer,intent(in)                     :: spec_rank
  integer,intent(in)                     :: iFic,filter
  integer,intent(in)                     :: filtMin,filtMax
  logical                                :: success

  !Local data
  real(WP)                          :: delta
  integer                           :: nbcpus
  integer                           :: nd
  logical                           :: res,boucle
  type(wavenumbers)                 :: WaveNum
  type(complex_data_layout)         :: ExK
  type(complex_data_layout),dimension(:),allocatable :: VelKBar
  type(complex_data_layout),dimension(:),allocatable :: VelK
  type(real_data_layout)            :: ScalBar
  type(complex_data_layout)         :: ScalK
  type(complex_data_layout)         :: ScalKBar
  type(complex_data_layout)         :: fieldPhiModK
  type(real_data_layout)            :: fieldPhiMod
  type(real_data_layout)            :: dTidxiMod
  type(real_data_layout),dimension(:),allocatable :: TiMod

  success = .false.
  delta = computeDelta(U)
  nbcpus = getnbcpus()

  !Init des tableau de travail
  if (.not. initDataLayout("ExK", ExK,(U%nx/2)+1,U%ny,U%nz, &
     & U%Lx,U%Ly,U%Lz,nbcpus,spec_rank,alongZ)) then
    write(6,'(a,i0,a)')'[ERROR] Not able to create Exk'
    return
  endif
  if (.not. initWorkArray(ExK,3,VelKBar) .or. &
     &.not. initWorkArray(ExK,3,VelK) .or. &
     &.not. initWorkArray(U,3,TiMod) .or. &
     &.not. copyStructOnly(Exk,ScalK) .or. &
     &.not. copyStructOnly(Exk,ScalKBar) .or. &
     &.not. copyStructOnly(Exk,fieldPhiModK) .or. &
     &.not. copyStructOnly(U,ScalBar) .or. &
     &.not. copyStructOnly(U,fieldPhiMod) .or. &
     &.not. copyStructOnly(U,dTidxiMod) ) then
    write(6,'(a,i0,a)')'[ERROR] Not enought memory 1'
    return
  endif

  !Calcul des nombres d'onde
  if (.not. initWN(WaveNum,spec_rank,U%nx,U%ny,U%nz) .or. &
     &.not.  computeWN(WaveNum,spec_rank,U%Lx,U%Ly,U%Lz,U%nx,U%ny,U%nz) ) then
    write(6,'(a,i0,a)')'[ERROR] Cannot compute WN'
    return
  endif

  nd = filtMin
  boucle = .true.
  do while (boucle)
    !passage dans l'espace spectrale
    call ftran(U,VelK(1),res)
    if (.not. res) return
    call ftran(V,VelK(2),res)
    if (.not. res) return
    call ftran(W,VelK(3),res)
    if (.not. res) return
    call ftran(Scal,ScalK,res)
    if (.not. res) return
    !Filtrage de tous les termes par le filtre 'filtre' à la longueur de coupure nd*delta
    call computeFilter(WaveNum,nd*delta,VelK(1),VelKBar(1),filter)
    call computeFilter(WaveNum,nd*delta,VelK(2),VelKBar(2),filter)
    call computeFilter(WaveNum,nd*delta,VelK(3),VelKBar(3),filter)
    call computeFilter(WaveNum,nd*delta,ScalK,ScalKBar,filter)
    !Retour dans l'espace physique de z
    call btran(ScalKBar,ScalBar,res)       !Zbar
    if (.not. res) return
    !Calcul des terms des modèles sous maille
    call RGMScaForTest(fieldPhiModK,U,VelKBar(1),VelKBar(2),VelKBar(3),ScalKBar,WaveNum,res,spec_rank,1.0_WP,Timod)
    if (.not. res) return
    call btran(fieldPhiModK,fieldPhiMod,res)
    if (.not. res) return
    !Calcul du modèle du gradient
    dTidxiMod%values = 1.0_WP/12.0_WP * fieldPhiMod%values
    if (.not. computeAprioriTest(iFic,U,V,W,Scal,ScalBar,Exk,spec_rank,nd,filter,(/fieldPhiMod/),&
                             & dTidxiMod,'RGM_mod',.false.,.false.,Timod) ) then
      write (6,'(a)') '[ERROR] computeAprioriTest failed'
      return
    endif
    nd = nd + 2
    if ( nd .gt. filtMax) then
      boucle = .false.
    endif
  enddo

  !Release memory
  call deletedataLayout(ExK)
  res = deleteWorkArray(VelKBar)
  res = deleteWorkArray(VelK)
  res = deleteWorkArray(TiMod)
  res = deleteWN(WaveNum)
  call deletedatalayout(ScalBar)
  call deletedatalayout(ScalK)
  call deletedatalayout(ScalKBar)
  call deletedatalayout(fieldPhiModK)
  call deletedatalayout(fieldPhiMod)
  call deletedatalayout(dTidxiMod )

  success = .true.

 end function apriori_RGM


  !-------------------------------------------------------------------------------------------
  !> Compute the SGS term and the viscous term of scalar equation and it dumps results in ASCII
  !! files
  !! @author Antoine Volllant
  !!    @param[in]    spec_rank       = mpi rank inside spectral communicator
  !!    @param[in]    nbcpus          = number of process in the cpu pool
  !!    @param[in]    U               = Component of velocity
  !!    @param[in]    V               = Component of velocity
  !!    @param[in]    W               = Component of velocity
  !!    @param[in]    Scal            = Scalar to compute
  !!    @param[in]    ite             = the number of file
  !!    @return       res             = logical value
  !-------------------------------------------------------------------------------------------
 function apriori_ScalarStat(iFic,U,V,W,Scal,spec_rank,filtMin,filtMax,filter) result(success)

  !I/O data
  type(real_data_layout),intent(in)      :: U, V ,W, Scal
  integer,intent(in)                     :: spec_rank
  integer,intent(in)                     :: iFic
  integer,intent(in)                     :: filter
  integer,intent(in)                     :: filtMin,filtMax
  logical                                :: success

  !Local data
  type(real_data_layout)                 :: UBar, VBar ,WBar, ScalBar
  type(wavenumbers)                      :: WaveNum
  integer                                :: nd
  real(WP)                               :: delta
  logical                                :: boucle

  success = .false.

  if (.not. copyStructOnly(U,UBar) .or. &
     &.not. copyStructOnly(U,VBar) .or. &
     &.not. copyStructOnly(U,WBar) .or. &
     &.not. copyStructOnly(U,ScalBar)) then
    write(6,'(a)') '[ERROR] in apriori_ScalarStat: not enought memory'
    return
  endif

  !Calcul des nombres d'onde
  if (.not. initWN(WaveNum,spec_rank,U%nx,U%ny,U%nz) .or. &
     &.not.  computeWN(WaveNum,spec_rank,U%Lx,U%Ly,U%Lz,U%nx,U%ny,U%nz) ) then
    write(6,'(a,i0,a)')'[ERROR] Cannot compute WN'
    return
  endif

  delta = computeDelta(U)
  nd = filtMin
  boucle = .true.
  do while (boucle)
    call computeFilter(WaveNum,nd*delta,U,UBar,filter)
    call computeFilter(WaveNum,nd*delta,V,VBar,filter)
    call computeFilter(WaveNum,nd*delta,W,WBar,filter)
    call computeFilter(WaveNum,nd*delta,Scal,ScalBar,filter)
    if (.not. computeScalarStat_Apriori(iFic,U,V,W,Scal,UBar,VBar,WBar,ScalBar,spec_rank, &
                             & nd,filter)) return
    nd = nd + 2
    if ( nd .gt. filtMax) then
      boucle = .false.
    endif
  enddo


  if(.not. deleteWN(WaveNum)) return
  call deletedatalayout(Ubar)
  call deletedatalayout(Vbar)
  call deletedatalayout(Wbar)
  call deletedatalayout(ScalBar)

  success = .true.

 end function apriori_ScalarStat



  !-------------------------------------------------------------------------------------------
  !> Compute the SGS term and the viscous term of scalar equation and it dumps results in ASCII
  !! files 
  !! @author Antoine Volllant 
  !!    @param[in]    spec_rank       = mpi rank inside spectral communicator
  !!    @param[in]    nbcpus          = number of process in the cpu pool 
  !!    @param[in]    U               = Component of velocity 
  !!    @param[in]    V               = Component of velocity 
  !!    @param[in]    W               = Component of velocity
  !!    @param[in]    Scal            = Scalar to compute 
  !!    @param[in]    ite             = the number of file
  !!    @return       res             = logical value 
  !-------------------------------------------------------------------------------------------
 function computeScalarStat_Apriori(iFic,U,V,W,Scal,UBar,VBar,WBar,ScalBar,spec_rank,nd,&
                           &filter) result(success)

  !I/O data
  type(real_data_layout),intent(in)      :: U, V ,W, Scal
  type(real_data_layout),intent(in)      :: UBar, VBar ,WBar, ScalBar
  integer,intent(in)                     :: spec_rank
  integer,intent(in)                     :: iFic
  integer,intent(in)                     :: nd
  integer,intent(in)                     :: filter
  logical                                :: success
  
  !Local data
  real(WP)                               :: mini,maxi,avg,var,nu
  real(WP)                               :: kurt,energ,delta,skew
  real(WP)                               :: Sc,moyNrjSGS,moyNrjDissVisq
  real(WP)                               :: moyConvNrjScal,moyNrjDiffvisq,moyConvNrjScalSGS 
  real(WP)                               :: moyNrjDiffVisqSGS,moyNrjDissVisqSGS
  type(real_data_layout)                 :: dTidxi
  type(real_data_layout)                 :: lapNrjScal
  type(real_data_layout)                 :: fieldDissExact
  type(real_data_layout)                 :: nrjScalBar
  type(real_data_layout)                 :: convNrjScal
  type(real_data_layout)                 :: lapNrjScalSGS
  type(real_data_layout)                 :: divNrjDissVisqSGS
  type(real_data_layout)                 :: nrjScalSGS2
  type(real_data_layout)                 :: nrjScalSGS1
  type(real_data_layout)                 :: nrjScal
  type(real_data_layout)                 :: nrjDissvisq
  type(real_data_layout)                 :: lapNrjScalVisqSGS 
  type(real_data_layout),dimension(:),allocatable :: uBarZZBar
  type(real_data_layout),dimension(:),allocatable :: gradScal
  type(real_data_layout),dimension(:),allocatable :: gradScalBar
  type(real_data_layout),dimension(:),allocatable :: gradBarScal
  type(real_data_layout),dimension(:),allocatable :: gradNrjScal
  type(real_data_layout),dimension(:),allocatable :: nrjDissVisqSGS
  type(complex_data_layout),dimension(:),allocatable :: VelKBar 
  type(complex_data_layout)              :: dealiasingK
  type(wavenumbers)                       :: WaveNum
  integer                                :: nbcpus
  logical                                :: res,existe
  character(len=64)                      :: filtre,file_name
  
  success = .false.
  nbcpus  = getnbcpus()
  delta   = computeDelta(U)
  CALL parser_read('Schmidt number for scalar 1',Sc)
  CALL parser_read('Viscosity',nu)

  if (.not. initDataLayout("dealiasingK",dealiasingK,(U%nx/2)+1,U%ny,U%nz, &
     & U%Lx,U%Ly,U%Lz,nbcpus,spec_rank,alongZ) .or. &
     &.not. initWorkArray(dealiasingK,3,VelKBar) .or. &
     &.not. initWorkArray(U,3,uBarZZBar) .or. &
     &.not. initWorkArray(U,3,gradScal) .or. &
     &.not. initWorkArray(U,3,gradScalBar) .or. &
     &.not. initWorkArray(U,3,gradBarScal) .or. &
     &.not. initWorkArray(U,3,gradNrjScal) .or. &
     &.not. initWorkArray(U,3,nrjDissVisqSGS) .or. &
     &.not. copyStructOnly(U,lapNrjScal) .or. &
     &.not. copyStructOnly(U,fieldDissExact) .or. &
     &.not. copyStructOnly(U,nrjScalBar) .or. &
     &.not. copyStructOnly(U,convNrjScal) .or. &
     &.not. copyStructOnly(U,nrjScalSGS1) .or. &
     &.not. copyStructOnly(U,lapNrjScalSGS) .or. &
     &.not. copyStructOnly(U,divNrjDissVisqSGS) .or. &
     &.not. copyStructOnly(U,nrjScalSGS2) .or. &
     &.not. copyStructOnly(U,nrjScalSGS1) .or. &
     &.not. copyStructOnly(U,nrjScal) .or. &
     &.not. copyStructOnly(U,nrjDissvisq) .or. &
     &.not. copyStructOnly(U,lapNrjScalVisqSGS) .or. &
     &.not. copyStructOnly(U,dTidxi)) then
    write(6,'(a)') '[ERROR] in computeScalarStat_Apriori : not enought memory'
    return
  endif 

  write(filtre,'(i2.2)') nd
  filtre = trim(adjustl(filtre))

  !Calcul des nombres d'onde   
  if (.not. initWN(WaveNum,spec_rank,U%nx,U%ny,U%nz) .or. &
     &.not.  computeWN(WaveNum,spec_rank,U%Lx,U%Ly,U%Lz,U%nx,U%ny,U%nz) ) then
    write(6,'(a,i0,a)')'[ERROR] Cannot compute WN'
    return
  endif

  energ   = computeEnergy   ( ScalBar ,spec_rank,nbcpus)
  mini    = computeFieldMin ( ScalBar , spec_rank )
  maxi    = computeFieldMax ( ScalBar , spec_rank )
  avg     = computeFieldAvg ( ScalBar , spec_rank )
  var     = computeFieldVar ( ScalBar , spec_rank ) 
  skew    = computeFieldSkew(ScalBar,spec_rank)
  kurt    = computeFieldKurt(ScalBar,spec_rank)

  if (.not.computeFieldPDF(iFic,ScalBar,300,spec_rank,1000)) then
    write(*,'(a)') '[ERROR] Failed to compute PDF of scalars fiels'
    return
  end if
  if (.not. compute_spectrum(ScalBar, iFic , spec_rank) ) then
    write(*,'(a)') '[ERROR] Failed to compute PDF of scalars fiels'
    return
  end if
  
  if (.not. computedTidxi(nbcpus,nd,filter,spec_rank, Scal, U, V ,W , dTidxi) ) return
  fieldDissExact%name = 'fieldDissExact_nd_'//trim(filtre)//'_' 
  fieldDissExact%values = dTidxi%values * ScalBar%values
  call ftran(fieldDissExact,dealiasingK,res)
  call dealiasIsotrope(dealiasingK,WaveNum)
  call btran(dealiasingK,fieldDissExact,res)
  if (.not. compute_spectrum(fieldDissExact, iFic , spec_rank) ) return
  
  !Calcul des termes moyens de l'equation de transport de l'nrj scalaire resolue
    !terme convectif
  nrjScalBar%values = 0.5_WP*(ScalBar%values)**2
  if (.not. computeFieldGradient(WaveNum,nrjScalBar,gradNrjScal(1),gradNrjScal(2),gradNrjScal(3),nbcpus,spec_rank)) return 
  convNrjScal%values = UBar%values * gradNrjScal(1)%values +&
                     & VBar%values * gradNrjScal(2)%values +&
                     & WBar%values * gradNrjScal(3)%values
  moyConvNrjScal = computeFieldAvg( convNrjScal , spec_rank )
    !terme de diffusion visqueuse
  if (.not. computeFieldLaplacian(WaveNum,nrjScalBar,lapNrjScal,nbcpus,spec_rank)) return
  moyNrjDiffvisq = computeFieldAvg( lapNrjScal , spec_rank )
  moyNrjDiffvisq = nu/Sc * moyNrjDiffvisq 
    !terme de dissipation visqueuse
  if (.not. computeFieldGradient(WaveNum,ScalBar,gradScal(1),gradScal(2),gradScal(3),nbcpus,spec_rank)) return 
  nrjDissvisq%values = gradScal(1)%values * gradScal(1)%values + &
                     & gradScal(2)%values * gradScal(2)%values + &
                     & gradScal(3)%values * gradScal(3)%values 
  moyNrjDissVisq = computeFieldAvg( nrjDissVisq , spec_rank )
  moyNrjDissVisq = nu/Sc * moyNrjDissVisq
  !Calcul des termes moyens de l'equation de transport de l'nrj scalaire sous maille
   !calcul de l'nrj scalaire sgs facon 1 ((zz)Bar-zBarzBar)/2
  nrjScal%values = (Scal%values)**2 
  call computeFilter(WaveNum,nd*delta,nrjScal,nrjScalBar,filter)
  nrjScalSGS1%values = 0.5_WP * nrjScalBar%values - nrjScalBar%values
   !calcul de l'nrj scalaire sgs facon 2 (z-zbar)**2/2
  nrjScalSGS2%values = 0.5_WP * (Scal%values - ScalBar%values) **2
    !terme convectif +uiBar d(Esgs)/dxi
  if (.not. computeFieldGradient(WaveNum,nrjScalBar,gradNrjScal(1),gradNrjScal(2),gradNrjScal(3),nbcpus,spec_rank)) return 
  convNrjScal%values = UBar%values * gradNrjScal(1)%values +&
                     & VBar%values * gradNrjScal(2)%values +&
                     & WBar%values * gradNrjScal(3)%values
  moyConvNrjScalSGS = computeFieldAvg( convNrjScal , spec_rank )
    !terme de dissipation visqueuse SGS -0.5*d/dxi((uiZZ)bar - uiBarZbarZbar)
  nrjDissVisqSGS(1)%values = U%values * Scal%values * Scal%values
  nrjDissVisqSGS(2)%values = V%values * Scal%values * Scal%values
  nrjDissVisqSGS(3)%values = W%values * Scal%values * Scal%values
  uBarZZBar(1)%values = U%values * 2.0_WP * nrjScalBar%values
  uBarZZBar(2)%values = V%values * 2.0_WP * nrjScalBar%values
  uBarZZBar(3)%values = W%values * 2.0_WP * nrjScalBar%values
  nrjDissVisqSGS(1)%values = nrjDissVisqSGS(1)%values - uBarZZBar(1)%values
  nrjDissVisqSGS(2)%values = nrjDissVisqSGS(2)%values - uBarZZBar(2)%values
  nrjDissVisqSGS(3)%values = nrjDissVisqSGS(3)%values - uBarZZBar(3)%values
  if (.not. computeFieldDivergence(WaveNum,nrjDissVisqSGS(1),nrjDissVisqSGS(2),nrjDissVisqSGS(3),&
                                           &divNrjDissVisqSGS,nbcpus,spec_rank)) return
  divNrjDissVisqSGS%values = - 0.5_WP*divNrjDissVisqSGS%values
  moyNrjDissVisqSGS = computeFieldAvg( divNrjDissVisqSGS , spec_rank ) 
    !terme de diffusion visqueuse SGS +D d^2Esgs/dxi^2
  if (.not. computeFieldLaplacian(WaveNum,nrjScalSGS1,lapNrjScalSGS,nbcpus,spec_rank)) return
  moyNrjDiffVisqSGS = computeFieldAvg( lapNrjScalVisqSGS , spec_rank ) 
  moyNrjDiffVisqSGS = nu/Sc * moyNrjDiffVisqSGS 
    !terme +D( dz/dxi**2Bar - dZbar/dxi**2)
  if (.not. computeFieldGradient(WaveNum,ScalBar,gradScalBar(1),gradScalBar(2),gradScalBar(3),nbcpus,spec_rank)) return 
  if (.not. computeFieldGradient(WaveNum,Scal,gradScal(1),gradScal(2),gradScal(3),nbcpus,spec_rank)) return 
  gradScalBar(1)%values = gradScalBar(1)%values * gradScalBar(1)%values
  gradScalBar(2)%values = gradScalBar(2)%values * gradScalBar(2)%values
  gradScalBar(3)%values = gradScalBar(3)%values * gradScalBar(3)%values
  gradScal(1)%values = gradScal(1)%values * gradScal(1)%values
  gradScal(2)%values = gradScal(2)%values * gradScal(2)%values
  gradScal(3)%values = gradScal(3)%values * gradScal(3)%values
  call computeFilter(WaveNum,nd*delta,gradScal(1),gradBarScal(1),filter)
  call computeFilter(WaveNum,nd*delta,gradScal(2),gradBarScal(2),filter)
  call computeFilter(WaveNum,nd*delta,gradScal(3),gradBarScal(3),filter)
  gradScalBar(1)%values = gradBarScal(1)%values - gradScalBar(1)%values 
  gradScalBar(2)%values = gradBarScal(2)%values - gradScalBar(2)%values
  gradScalBar(3)%values = gradBarScal(3)%values - gradScalBar(3)%values
  gradBarScal(1)%values = gradScalBar(1)%values + gradScalBar(2)%values + gradScalBar(3)%values
  moyNrjSGS = computeFieldAvg( gradBarScal(1) , spec_rank ) 
  moyNrjSGS = nu/Sc * moyNrjSGS
  
  !Termes calcule
  if (spec_rank .eq. 0 ) then
    write(file_name,'(a)') 'scalaire_aprioriTest_nd_'//trim(filtre)//'.table'
    inquire( file=trim(adjustl(file_name)), exist=existe)
    open(unit=66, file=trim(adjustl(file_name)), form="FORMATTED",position="append")
    if ( .not. existe) then
      write(66, '(a)') '#numFic min(z) <z> max(z) var(z) skew(z) kurt(z) zz/2  &
                            & u_id(E_z)/dxi Dd^2(E_z)/dxi^2  D(d(E_z)/dxi)^2 u_id(Esgs)/dxi &
                            & Dd^2(Esgs)/dxi^2 -0.5*d/dxi((u_izz)bar-u_iBarzbarzbar)&
                            & +D(d(z)/dxi^2Bar-di(zbar)/dxi^2)'
    endif
    write (66,'(i0,1x,14(f15.8,1x))') iFic,mini,avg,maxi,var,skew,kurt,energ,&
                                          &moyConvNrjScal,moyNrjDiffvisq,moyNrjDissVisq,&
                                  &moyConvNrjScalSGS, moyNrjDiffVisqSGS,moyNrjDissVisqSGS,moyNrjSGS
    close(66)
  endif

  call deletedatalayout(convNrjScal)
  call deletedatalayout(dealiasingK)
  call deletedatalayout(lapNrjScal)
  call deletedatalayout(fieldDissExact)
  call deletedatalayout(nrjScalBar)
  call deletedatalayout(nrjScalSGS1)
  call deletedatalayout(lapNrjScalSGS)
  call deletedatalayout(divNrjDissVisqSGS)
  call deletedatalayout(nrjScalSGS2)
  call deletedatalayout(nrjScalSGS1)
  call deletedatalayout(nrjScal)
  call deletedatalayout(nrjDissvisq)
  call deletedatalayout(lapNrjScalVisqSGS)
  call deletedatalayout(dTidxi)
  if (.not. deleteWorkArray(VelKBar) .or. &
     &.not. deleteWorkArray(gradNrjScal) .or. &
     &.not. deleteWorkArray(gradScal) .or. &
     &.not. deleteWorkArray(uBarZZBar) .or. &
     &.not. deleteWorkArray(nrjDissVisqSGS) .or. &
     &.not. deleteWorkArray(gradScalBar) .or. &
     &.not. deleteWN(WaveNum) .or. &
     &.not. deleteWorkArray(gradBarScal)) then
    write(6,'(a)') '[ERROR] memory leackage'
  endif
  success = .true.
  
 end function computeScalarStat_Apriori


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
!!        This function compute the spectrum of a given field. It perform
!!    a FFT in order to obtain it in the Fourrier space (wich is required
!!    to computing the spectrum). It could be perform "out a simulation"
!!    (ie from a field write in an output).
!------------------------------------------------------------------------------
function compute_scalar_spectrum_rescaled(field, WN, chi, epsil, eta_k, eta_b, nbcpus, visco, k_lambda, spec_rank) result(success)


    ! Input/Output
    type(REAL_DATA_LAYOUT), intent(in)      :: field
    type(WaveNumbers), intent(in)           :: WN               ! associated wave numbers
    real(WP), intent(in)                    :: chi, epsil, eta_k, eta_b, visco, k_lambda
    integer, intent(in)                     :: spec_rank
    integer, intent(in)                     :: nbcpus
    logical                                 :: success
    ! Others variables
    type(COMPLEX_DATA_LAYOUT)               :: spec_field       ! to store spectral values of "field"
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
    real(WP), allocatable, dimension(:,:)   :: rescaling        ! to compute rescaling
    real(WP), allocatable, dimension(:)     :: weight           ! to apply mask for the different resaling
    real(WP)                                :: cst_B, cst_OC    ! rescaling constants
    real(WP)                                :: cst_Kr, cst_Ba   ! rescaling constants for dissipation
    real(WP)                                :: max_sc, min_sc   ! spectrum range associated to the different rescaling

    ! Init
    success = .false.

    ! ===== Check if spectral space step is the same above all direction =====
    if((field%Lx/=field%Ly).or.(field%Lx/=field%Lz)) then
        write(*,'(a)') '[warning] Spectrum not yet implemented if not(Lx=Ly=Lz)'
        return
    end if



    ! ===== Go in spectral space =====
    ! Allocate storage in fourrier space
    if (.not. initDataLayout("scalar_rescaled_spectral", spec_field,(field%nx/2)+1, &
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
    success = .false.
    dk=2.0_WP*acos(-1.0_WP)/spec_field%Lx
    kc = 2.0_WP*acos(-1.0_WP)*(spec_field%nx-1)/spec_field%Lx
    eps = dk/2.0


    ! ===== Compute spectrum =====
    size_spec = size(WN%kx)
    allocate(spectrum(size_spec))
    allocate(WN_norm(size_spec))
    WN_norm = -1
    spectrum = 0
    ! Compute
    do k = spec_field%zmin,spec_field%zmax
        do j =  spec_field%ymin,spec_field%ymax
            do i =  spec_field%xmin,spec_field%xmax
                ! Compute norm of wave number
                kk = WN%kx(i)*WN%kx(i) + &
                    & WN%ky(j)*WN%ky(j) + &
                    & WN%kz(k)*WN%kz(k)
                kk = sqrt(real(kk))
                ! Compute indice
                ik = nint(kk/dk)
                if ((kk>eps).and.(kk<=kc)) then
                ! And then update spectrum
                spectrum(ik) = spectrum(ik) &
                    & + spec_field%values(i,j,k)*conjg(spec_field%values(i,j,k))
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
    allocate(rescaling(size_spec,4))
    ! -- Intertial range --
    rescaling(:,1) = chi*(epsil**(-1._WP/3._WP))*(WN_norm**(-5._WP/3._WP))
    ! Compute constant coefficient
    max_sc = k_lambda
    where((WN_norm<max_sc).and.(spectrum>0.0_WP))
        weight = 1.0_WP
    elsewhere
        weight = 0.0_WP
    end where
    cst_OC = max(1.0_WP,sum(weight))
    where(weight>0.0_WP)
        weight = weight*spectrum/rescaling(:,1)
    end where
    cst_OC = sum(weight)/cst_OC
    !if (cst_OC>0.0_WP) then
    rescaling(:,1) = cst_OC*rescaling(:,1)
    !end if

    ! -- If Schmidt bigger than one, k^-1 rescaling --
    !if(eta_B>eta_k) then
        rescaling(:,2) =  (chi*sqrt(visco/epsil)/WN_norm)
! XXX Todo : parser si présent
        min_sc = 2*k_lambda
        max_sc = min(8.0_WP*k_lambda, k_lambda + ((1.0_WP/eta_B)+k_lambda)/2.0_WP)
        where((WN_norm>min_sc).and.(WN_norm<max_sc))
            weight = 1.0_WP
        elsewhere
            weight = 0.0_WP
        end where
        cst_B = sum(weight)
        !if (spec_rank==0) then
        !    if(cst_B<=1) then
        !    else
        !    end if
        !
        if(spec_rank==0) print*, 'interpol range for scalar = ', min_sc, ' to', max_sc
        cst_B = max(1.0_WP, cst_B)
        where(weight>0.0_WP)
            weight = weight*spectrum/rescaling(:,2)
        end where
        cst_B = sum(weight)/cst_B
        !if (cst_B>0.0_WP) then
        rescaling(:,2) = cst_B*rescaling(:,2)
        !end if

    ! -- If Schmidt bigger than one, decay rescaling --
        cst_Kr = cst_B
        cst_Ba = cst_B
        if(parser_is_defined('cst_Kr')) call parser_read('cst_Kr', cst_Kr)
        if(parser_is_defined('cst_Ba')) call parser_read('cst_Ba', cst_Ba)

        rescaling(:,3) = rescaling(:,2)*exp(-cst_Ba*(eta_B*WN_norm)**2)
        rescaling(:,4) = sqrt(6.0_WP*cst_Kr)*eta_B*WN_norm
        rescaling(:,4) = rescaling(:,2)*(1.0_WP + rescaling(:,4))*exp(-rescaling(:,4))
        !end if


    ! ===== Write output =====
    if (spec_rank==0) then
        file_id = iopen()
        write(file_name,'(a)') 'rescaled_scalar_spectrum.table'
        open(unit=file_id, file=trim(adjustl(file_name)), form="FORMATTED", ERR=300)
        write(file_id, '(8(a,e14.6),a)') '#### Rescaled scalar spectrum - espilon = ',  &
            & epsil, ' chi = ', chi, ' eta_k = ', eta_k, 'eta_b = ', eta_b, 'cst_OC',   &
            & cst_OC, 'cst_B', cst_B, 'cst_Ba ', cst_Ba, 'cst_Kr ', cst_Kr, ' ####'
        write(file_id, '(a,a,a,a)') &
            & '# Fields are : k // k.eta_B // E // '           ,&
            & 'cst_OC*chi*(epsil**(-1./3.))*(k**(-5./**3.) // ',&
            & 'cst_B*(chi*sqrt(visco/epsil)/k)u // '           ,&
            & ' dissipation gaussienne et exponentielle'
        do ik = 1, size_spec
            if (WN_norm(ik)>=0) write(file_id, '(7(e14.6,x))') &
                & WN_norm(ik), WN_norm(ik)*eta_b, spectrum(ik), rescaling(ik,:)
        end do
        close(file_id)
        file_id = iclose(file_id)
    end if

    ! ===== Free memory =====
    deallocate(spectrum)
    deallocate(WN_norm)
    deallocate(rescaling)
    deallocate(weight)
    call deleteDataLayout(spec_field)

    ! ===== Return sucess =====
    success = .true.
    return

    ! ===== IO error =====
300 write(6,'(a)') 'Unable to open file "'//trim(file_name)//'" for spectrum'
    return

end function compute_scalar_spectrum_rescaled



#ifdef POSTCTR
!-------------------------------------------------------------------------------------------
!> Post-processing for the 2012 CTR summer program
!!
!! @author Antoine Vollant and Guillaume Balarac
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in]    nbcpus          = number of process in the cpu pool
!!    @param[in]    U               = Component of velocity
!!    @param[in]    V               = Component of velocity
!!    @param[in]    W               = Component of velocity
!!    @param[in]    Scal            = Scalar to compute
!!    @param[in]    ite             = the number of file
!!    @return       res             = logical value
!-------------------------------------------------------------------------------------------
 function gradientmodOE(U,V,W,Scal,nbcpus,spec_rank) result (success)



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
  type(real_data_layout)                             :: T11, T12, T13, T22, T23, T33
  type(real_data_layout)                             :: Z11, Z12, Z13, Z22, Z23, Z33
  type(real_data_layout)                             :: O12, O13, O23
  type(real_data_layout),dimension(:),allocatable    :: phi, phidiss, TiDNS, TiM, SijmTj
  type(real_data_layout),dimension(:),allocatable    :: temp, temp1, temp2, temp3, temp4 , temp5, temp6
  type(real_data_layout),dimension(:),allocatable    :: VelScal_ft, Vel_ft, dZdxi_ft, SdZdxi_ft, Li, Mi, Pi, Qi
  type(real_data_layout)                             :: Scal_ft
  type(real_data_layout)                             :: S11_ft,S12_ft,S13_ft, S22_ft, S23_ft, S33_ft
  type(real_data_layout)                             :: dTidxi_DNS
  type(real_data_layout)                             :: dTidxi_M
  type(real_data_layout)                             :: Uf, Vf, Wf,Scalf

  type(complex_data_layout),dimension(:),allocatable :: tempK, tempKf
  type(complex_data_layout)                          :: ExK

  type(wavenumbers)                                  :: ScalWaveN
  real(WP)                                           :: qe_gm, ie_gm, deltx, DynCoef, filtersize
  !real(WP)                                           :: Coef1, Coef2
  logical                                            :: res
  integer                                            :: bin1,bin2,bin3,typeDisc,n_iter
  integer                                            :: filter
  integer                                            :: nd,ii
  integer                                            :: i,j,k

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
  if (.not. copyStructOnly(U,S11,"S11") .or. & 
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
     &.not. copyStructOnly(U,O12,"O12") .or. &
     &.not. copyStructOnly(U,O13,"O13") .or. &
     &.not. copyStructOnly(U,O23,"O23") .or. &
     &.not. copyStructOnly(U,dTidxi_M,"dTidxi_M") .or. &
     &.not. copyStructOnly(U,dTidxi_DNS,"dTidxi_DNS") .or. &
     &.not. copyStructOnly(U,Uf,"Uf") .or. &
     &.not. copyStructOnly(U,Vf,"Vf") .or. &
     &.not. copyStructOnly(U,Wf,"Wf") .or. &
     &.not. copyStructOnly(U,Scalf,"Scalf") .or. &
     &.not. copyStructOnly(U,Scal_ft,"Scal_ft")) then
     write(6,'(a,i0,a)')'[ERROR] copyStructOnly for TYPE (REAL_DATA_LAYOUT) : not enought memory!'
     return
  endif
  if (.not. initWorkArray(U,3,temp) .or. &
     &.not. initWorkArray(U,3,temp1) .or. &
     &.not. initWorkArray(U,3,temp2) .or. &
     &.not. initWorkArray(U,3,temp3) .or. &
     &.not. initWorkArray(U,3,temp4) .or. &
     &.not. initWorkArray(U,3,temp5) .or. &
     &.not. initWorkArray(U,3,temp6) .or. &
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
     &.not. initWorkArray(U,3,Mi) .or. &
     &.not. initWorkArray(U,3,Pi) .or. &
     &.not. initWorkArray(U,3,Qi)) then
     write(6,'(a,i0,a)')'[ERROR] initWorkArray for TYPE (REAL_DATA_LAYOUT) : not enought memory!'
     return
  endif

  res = initWN(ScalWaveN,spec_rank,Scal%nx,Scal%ny,Scal%nz)
  res = computeWN(ScalWaveN,spec_rank,Scal%Lx,Scal%Ly,Scal%Lz,Scal%nx,Scal%ny,Scal%nz)

  if (spec_rank.EQ.0) then
     open(14,file='IE_QE_GM.out',form='formatted')
     open(11,file='IE_QE_diff_GM.out',form='formatted')
     open(12,file='mean_diff_GM.out',form='formatted')
     open(13,file='mean_diss_GM.out',form='formatted')

     open(20,file='IE_QE_GM2.out',form='formatted')
     open(21,file='IE_QE_diff_GM2.out',form='formatted')
     open(22,file='mean_diff_GM2.out',form='formatted')

     open(30,file='IE_QE_S-O2.out',form='formatted')
     open(31,file='IE_QE_diff_S-O2.out',form='formatted')
     open(32,file='mean_diff_S-O2.out',form='formatted')

     open(40,file='IE_QE_S+O2.out',form='formatted')
     open(41,file='IE_QE_diff_S+O2.out',form='formatted')
     open(42,file='mean_diff_S+O2.out',form='formatted')

     open(50,file='IE_QE_S-2.out',form='formatted')
     open(51,file='IE_QE_diff_S-2.out',form='formatted')
     open(52,file='mean_diff_S-2.out',form='formatted')
     open(53,file='mean_diss_S-2.out',form='formatted')

     open(60,file='IE_QE_S-O.out',form='formatted')
     open(61,file='IE_QE_diff_S-O.out',form='formatted')
     open(62,file='mean_diff_S-O.out',form='formatted')

     open(70,file='IE_QE_S+O.out',form='formatted')
     open(71,file='IE_QE_diff_S+O.out',form='formatted')
     open(72,file='mean_diff_S+O.out',form='formatted')

     open(80,file='IE_QE_DynSmag.out',form='formatted')
     open(81,file='IE_QE_diff_DynSmag.out',form='formatted')
     open(82,file='mean_diff_DynSmag.out',form='formatted')
     open(83,file='mean_diss_DynSmag.out',form='formatted')

     open(90,file='IE_QE_S+2.out',form='formatted')
     open(91,file='IE_QE_diss_S+2.out',form='formatted')
     open(92,file='mean_diss_S+2.out',form='formatted')

     open(100,file='IE_QE_O2.out',form='formatted')
     open(101,file='IE_QE_diss_O2.out',form='formatted')
     open(102,file='mean_diss_O2.out',form='formatted')

     open(110,file='IE_QE_wang.out',form='formatted')
     open(111,file='IE_QE_diss_wang.out',form='formatted')
     open(112,file='mean_diss_wang.out',form='formatted')
  endif

  call ftran(U,tempK(1),res)
  call ftran(V,tempK(2),res)
  call ftran(W,tempK(3),res)
  ie_gm = computeDiv( spec_rank,tempK(1),tempK(2),tempK(3),tempK(4),ScalWaveN )
  IF (spec_rank.EQ.0) WRITE(6,'(a,g14.8)')'  [INFO] Divergence max of U field is ',ie_gm
  IF (.NOT. NullifyDiv(tempK(1),tempK(2),tempK(3), ScalWaveN, spec_rank)) RETURN
  ie_gm = computeDiv( spec_rank,tempK(1),tempK(2),tempK(3),tempK(4),ScalWaveN )
  IF (spec_rank.EQ.0) WRITE(6,'(a,g14.8)')'  [INFO] Divergence max of U field is ',ie_gm

  !! -> loop on filter size
  do ii = 1,8

    nd = 2*ii
    filtersize = real(nd)*deltx
    if (spec_rank .eq. 0)  print*,'[INFO] ',ii,' TAILLE de FILTRE = ',filtersize

    !!!!!!!!!!!!!!!!!
    ! -> Exact divergence of the SGS scalar flux
    !!!!!!!!!!!!!!!!!
    !Calcul du terme de flux scalaire sous maille filtré exact
    res = computedTidxi(nbcpus,nd,filter,spec_rank, Scal, U, V ,W , dTidxi_DNS )

    !!!!!!!!!!!!!!!!!
    ! -> Compute the filtered velocities and the filtered scalar
    !!!!!!!!!!!!!!!!!
    call ftran(U,tempK(1),res)
    call ftran(V,tempK(2),res)
    call ftran(W,tempK(3),res)
    call ftran(Scal,tempK(4),res)
    call computeFilter(ScalWaveN,filtersize,tempK(1),tempKf(1),filter) !U bar
    call computeFilter(ScalWaveN,filtersize,tempK(2),tempKf(2),filter) !V bar
    call computeFilter(ScalWaveN,filtersize,tempK(3),tempKf(3),filter) !W bar
    call computeFilter(ScalWaveN,filtersize,tempK(4),tempKf(4),filter) !Z bar
    call btran(tempKf(1),Uf,res) !U bar
    call btran(tempKf(2),Vf,res) !V bar
    call btran(tempKf(3),Wf,res) !W bar
    call btran(tempKf(4),Scalf,res) !Z bar




    !!!!!!!!!!!!!!!!!
    ! -> Compute TiDNS
    !!!!!!!!!!!!!!!!!
    temp5(1)%values = U%values * Scal%values
    temp5(2)%values = V%values * Scal%values
    temp5(3)%values = W%values * Scal%values
    call ftran(temp5(1),tempK(1),res)
    call ftran(temp5(2),tempK(2),res)
    call ftran(temp5(3),tempK(3),res)
    call computeFilter(ScalWaveN,filtersize,tempK(1),tempKf(1),filter) !U bar
    call computeFilter(ScalWaveN,filtersize,tempK(2),tempKf(2),filter) !V bar
    call computeFilter(ScalWaveN,filtersize,tempK(3),tempKf(3),filter) !W bar
    call btran(tempKf(1),TiDNS(1),res) !U bar
    call btran(tempKf(2),TiDNS(2),res) !V bar
    call btran(tempKf(3),TiDNS(3),res) !W bar

    !!!!!!!!!!!!!!!!!
    !-> IE and QE for the full gradient model: dTidxi = delta^2 / 12 * d/dxi ( dui/dxj*dZdxj )  => Phi1 = d/dxi ( dui/dxj*dZdxj )
    !!!!!!!!!!!!!!!!!
    ! compute duidxj and dzdxj
    ! temp1 <- dUfdxi
    ! temp2 <- dVfdxi
    ! temp3 <- dWfdxi
    ! temp4 <- dScalfdxi
    res = computeFieldGradient(ScalWaveN,Uf,temp1(1),temp1(2),temp1(3),nbcpus,spec_rank)
    res = computeFieldGradient(ScalWaveN,Vf,temp2(1),temp2(2),temp2(3),nbcpus,spec_rank)
    res = computeFieldGradient(ScalWaveN,Wf,temp3(1),temp3(2),temp3(3),nbcpus,spec_rank)
    res = computeFieldGradient(ScalWaveN,Scalf,temp4(1),temp4(2),temp4(3),nbcpus,spec_rank)

    ! product : temp <- duidxj*dzdxj
    TiM(1)%values = real(filtersize)**2.0 / 12.0 * (temp1(1)%values*temp4(1)%values + temp1(2)%values*temp4(2)%values + temp1(3)%values*temp4(3)%values)
    TiM(2)%values = real(filtersize)**2.0 / 12.0 * (temp2(1)%values*temp4(1)%values + temp2(2)%values*temp4(2)%values + temp2(3)%values*temp4(3)%values)
    TiM(3)%values = real(filtersize)**2.0 / 12.0 * (temp3(1)%values*temp4(1)%values + temp3(2)%values*temp4(2)%values + temp3(3)%values*temp4(3)%values)

    ! take divergence -- this is the Phi1 parameter
    res = computeFieldDivergence(ScalWaveN,TiM(1),TiM(2),TiM(3),phi(1),nbcpus,spec_rank)
    dTidxi_M%values = phi(1)%values

    ! Compute the optimal estimator of GM
    res = computeEONVar( spec_rank,(/phi(1)/),& !Les variables
                       & dTidxi_DNS,&                            !Le champ à estimer
                       & (/typeDisc/),&     !Les discrétisations par variable
                       & (/bin1/),&                 !le nombre d'intervalles par variable
                       & fieldOpt=temp(1))                    !champ optimal en sortie

    ! IE of GM
    ie_gm = computeMSE(dTidxi_DNS,temp(1),spec_rank,.true.)     ! true = adim par la variance

    ! QE of GM
    qe_gm = computeMSE(dTidxi_DNS,dTidxi_M,spec_rank,.true.)     ! true = adim par la variance

    if (spec_rank.EQ.0) then
       write(14,*) nd, ie_gm, qe_gm
    end if

    ! Dissipation term: !!! warning !!! Dissipation is defined as zdTidxi here

    ! DNS dissipation <- temp(2)
    temp(2)%values = Scalf%values*dTidxi_DNS%values
    ! model parameter
    phidiss(1)%values = Scalf%values*phi(1)%values
    ! optimal estimator dissipation <- temp(1)
    res = computeEONVar( spec_rank,(/phidiss(1)/),& !Les variables
                       & temp(2),&                            !Le champ à estimer
                       & (/typeDisc/),&     !Les discrétisations par variable
                       & (/bin1/),&                 !le nombre d'intervalles par variable
                       & fieldOpt=temp(1))                    !champ optimal en sortie
    ! Model dissipation <- temp(3)
    temp(3)%values = Scalf%values*dTidxi_M%values
    ! IE et QE of GM for dissipation
    ie_gm = computeMSE(temp(2),temp(1),spec_rank,.true.)     ! true = adim par la variance
    qe_gm = computeMSE(temp(2),temp(3),spec_rank,.true.)     ! true = adim par la variance
    if (spec_rank.EQ.0) then
       write(11,*) nd, ie_gm, qe_gm
    end if
    ! global averaging of dissipation: ie_gm <- DNS dissipation and qe_gm <- model dissipation
    ie_gm = computeFieldAvg( temp(2) , spec_rank , nbcpus )
    qe_gm = computeFieldAvg( temp(3) , spec_rank , nbcpus )
    if (spec_rank.EQ.0) then
       write(12,*) nd, ie_gm, qe_gm
    end if

    temp(2)%values = TiDNS(1)%values*temp4(1)%values + TiDNS(2)%values*temp4(2)%values + TiDNS(2)%values*temp4(2)%values
    temp(2)%name = "DNS_Diss"
    if (.not. computeFieldPDF(nd,temp(2),bin1,spec_rank,n_iter)) then
         write (6,'()') '[ERROR] '
         return
    endif
    ie_gm = computeFieldAvg( temp(2) , spec_rank , nbcpus )

    temp(2)%values = TiM(1)%values*temp4(1)%values + TiM(2)%values*temp4(2)%values + TiM(3)%values*temp4(3)%values
    temp(2)%name = "GM_Diss"
    if (.not. computeFieldPDF(nd,temp(2),bin1,spec_rank,n_iter)) then
         write (6,'()') '[ERROR] '
         return
    endif
    qe_gm = computeFieldAvg( temp(2) , spec_rank , nbcpus )
    if (spec_rank.EQ.0) then
       write(13,*) nd, ie_gm, qe_gm
    end if

    ! compute Sij = 1/2 * ( duidxj + dujdxi) and Oij = 1/2 * ( duidxj - dujdxi )
    S11%values = temp1(1)%values
    S12%values = 0.5 * ( temp1(2)%values + temp2(1)%values )
    S13%values = 0.5 * ( temp1(3)%values + temp3(1)%values )
    S22%values = temp2(2)%values
    S23%values = 0.5 * ( temp2(3)%values + temp3(2)%values )
    S33%values = temp3(3)%values
    O12%values = 0.5 * ( temp1(2)%values - temp2(1)%values )
    O13%values = 0.5 * ( temp1(3)%values - temp3(1)%values )
    O23%values = 0.5 * ( temp2(3)%values - temp3(2)%values )

!    do k = U%zmin, U%zmax
!       do j = U%ymin, U%ymax
!          do i = U%xmin, U%xmax
!             qe_gm = S12%values(i,j,k) + O12%values(i,j,k)
!             if (qe_gm.ne.temp1(2)%values(i,j,k)) print*,'12 something wrong',qe_gm,temp1(2)%values(i,j,k)
!             qe_gm = S12%values(i,j,k) - O12%values(i,j,k)
!             if (qe_gm.ne.temp2(1)%values(i,j,k)) print*,'21 something wrong',qe_gm,temp2(1)%values(i,j,k)
!
!             qe_gm = S13%values(i,j,k) + O13%values(i,j,k)
!             if (qe_gm.ne.temp1(3)%values(i,j,k)) print*,'13 something wrong',qe_gm,temp1(3)%values(i,j,k)
!             qe_gm = S13%values(i,j,k) - O13%values(i,j,k)
!             if (qe_gm.ne.temp3(1)%values(i,j,k)) print*,'31 something wrong',qe_gm,temp3(1)%values(i,j,k)
!
!             qe_gm = S23%values(i,j,k) + O23%values(i,j,k)
!             if (qe_gm.ne.temp2(3)%values(i,j,k)) print*,'23 something wrong',qe_gm,temp2(3)%values(i,j,k)
!             qe_gm = S23%values(i,j,k) - O23%values(i,j,k)
!             if (qe_gm.ne.temp3(2)%values(i,j,k)) print*,'32 something wrong',qe_gm,temp3(2)%values(i,j,k)
!          enddo
!       enddo
!    enddo

    !!!!!!!!!!!!!!!!!
    !-> IE and QE for the Smagorinsky model: dTidxi = C * delta^2 * d/dxi ( |Sij| * dZdxi )
    !   => Phi1 = d/dxi ( |Sij| * dZdxi )
    !!!!!!!!!!!!!!!!!

    ! compute |Sij| = sqrt(2.0*Sij*Sij)   <-- temp5(1)
    call computeTensorNorme1(S11,S12,S13,S22,S23,S33,temp5(1),res)

    ! compute dynamic coef and model
    !DynCoef = - 0.141
    call ftran(Uf,tempK(1),res)
    call ftran(Vf,tempK(2),res)
    call ftran(Wf,tempK(3),res)
    call ftran(Scalf,tempK(4),res)
    ! Compute Li = TestFilter (uiz) - TestFilter(ui)*TestFilter(z)
    ! Compute TestFilter(ui) and TestFilter(z)
    call ftran(Uf,tempK(1),res)
    call ftran(Vf,tempK(2),res)
    call ftran(Wf,tempK(3),res)
    call ftran(Scalf,tempK(4),res)
    call computeFilter(ScalWaveN,2.0*filtersize,tempK(1),tempKf(1),filter) !U bar
    call computeFilter(ScalWaveN,2.0*filtersize,tempK(2),tempKf(2),filter) !V bar
    call computeFilter(ScalWaveN,2.0*filtersize,tempK(3),tempKf(3),filter) !W bar
    call computeFilter(ScalWaveN,2.0*filtersize,tempK(4),tempKf(4),filter) !Z bar
    call btran(tempKf(1),Vel_ft(1),res) !U bar
    call btran(tempKf(2),Vel_ft(2),res) !V bar
    call btran(tempKf(3),Vel_ft(3),res) !W bar
    call btran(tempKf(4),Scal_ft,res) !Z bar
    VelScal_ft(1)%values = Uf%values * Scalf%values
    VelScal_ft(2)%values = Vf%values * Scalf%values
    VelScal_ft(3)%values = Wf%values * Scalf%values
    call ftran(VelScal_ft(1),tempK(1),res)
    call ftran(VelScal_ft(2),tempK(2),res)
    call ftran(VelScal_ft(3),tempK(3),res)
    call computeFilter(ScalWaveN,2.0*filtersize,tempK(1),tempKf(1),filter) !U bar
    call computeFilter(ScalWaveN,2.0*filtersize,tempK(2),tempKf(2),filter) !V bar
    call computeFilter(ScalWaveN,2.0*filtersize,tempK(3),tempKf(3),filter) !W bar
    call btran(tempKf(1),VelScal_ft(1),res) !U bar
    call btran(tempKf(2),VelScal_ft(2),res) !V bar
    call btran(tempKf(3),VelScal_ft(3),res) !W bar
    ! Compute Li
    Li(1)%values = VelScal_ft(1)%values - Vel_ft(1)%values * Scal_ft%values
    Li(2)%values = VelScal_ft(2)%values - Vel_ft(2)%values * Scal_ft%values
    Li(3)%values = VelScal_ft(3)%values - Vel_ft(3)%values * Scal_ft%values

    ! Compute Mi = (2 * filtersize)^2* |Sij_ft|*dZdxi_ft - filtersize^2 * SdZdxi_ft
    ! Compute Sij_ft by filtering Sij
    call ftran(S11,tempK(1),res)
    call ftran(S12,tempK(2),res)
    call ftran(S13,tempK(3),res)
    call computeFilter(ScalWaveN,2.0*filtersize,tempK(1),tempKf(1),filter) !U bar
    call computeFilter(ScalWaveN,2.0*filtersize,tempK(2),tempKf(2),filter) !V bar
    call computeFilter(ScalWaveN,2.0*filtersize,tempK(3),tempKf(3),filter) !W bar
    call btran(tempKf(1),S11_ft,res) !U bar
    call btran(tempKf(2),S12_ft,res) !V bar
    call btran(tempKf(3),S13_ft,res) !W bar
    call ftran(S22,tempK(1),res)
    call ftran(S23,tempK(2),res)
    call ftran(S33,tempK(3),res)
    call computeFilter(ScalWaveN,2.0*filtersize,tempK(1),tempKf(1),filter) !U bar
    call computeFilter(ScalWaveN,2.0*filtersize,tempK(2),tempKf(2),filter) !V bar
    call computeFilter(ScalWaveN,2.0*filtersize,tempK(3),tempKf(3),filter) !W bar
    call btran(tempKf(1),S22_ft,res) !U bar
    call btran(tempKf(2),S23_ft,res) !V bar
    call btran(tempKf(3),S33_ft,res) !W bar
    ! Compute |Sij_ft|  <-- temp5(2)
    call computeTensorNorme1(S11_ft,S12_ft,S13_ft,S22_ft,S23_ft,S33_ft,temp5(2),res)
    ! Compute dZdxi_ft by filtering of dZdxi (which is in temp4)
    call ftran(temp4(1),tempK(1),res)
    call ftran(temp4(2),tempK(2),res)
    call ftran(temp4(3),tempK(3),res)
    call computeFilter(ScalWaveN,2.0*filtersize,tempK(1),tempKf(1),filter) !U bar
    call computeFilter(ScalWaveN,2.0*filtersize,tempK(2),tempKf(2),filter) !V bar
    call computeFilter(ScalWaveN,2.0*filtersize,tempK(3),tempKf(3),filter) !W bar
    call btran(tempKf(1),dZdxi_ft(1),res) !U bar
    call btran(tempKf(2),dZdxi_ft(2),res) !V bar
    call btran(tempKf(3),dZdxi_ft(3),res) !W bar

    ! Compute SdZdxi_ft
    SdZdxi_ft(1)%values = temp5(1)%values * temp4(1)%values
    SdZdxi_ft(2)%values = temp5(1)%values * temp4(2)%values
    SdZdxi_ft(3)%values = temp5(1)%values * temp4(3)%values
    call ftran(SdZdxi_ft(1),tempK(1),res)
    call ftran(SdZdxi_ft(2),tempK(2),res)
    call ftran(SdZdxi_ft(3),tempK(3),res)
    call computeFilter(ScalWaveN,2.0*filtersize,tempK(1),tempKf(1),filter) !U bar
    call computeFilter(ScalWaveN,2.0*filtersize,tempK(2),tempKf(2),filter) !V bar
    call computeFilter(ScalWaveN,2.0*filtersize,tempK(3),tempKf(3),filter) !W bar
    call btran(tempKf(1),SdZdxi_ft(1),res) !U bar
    call btran(tempKf(2),SdZdxi_ft(2),res) !V bar
    call btran(tempKf(3),SdZdxi_ft(3),res) !W bar
    ! Compute Mi = (2 * filtersize)^2* |Sij_ft|*dZdxi_ft - filtersize^2 * SdZdxi_ft
    Mi(1)%values = (2.0 * filtersize)**2.0 * temp5(2)%values*dZdxi_ft(1)%values - filtersize**2.0 * SdZdxi_ft(1)%values
    Mi(2)%values = (2.0 * filtersize)**2.0 * temp5(2)%values*dZdxi_ft(2)%values - filtersize**2.0 * SdZdxi_ft(2)%values
    Mi(3)%values = (2.0 * filtersize)**2.0 * temp5(2)%values*dZdxi_ft(3)%values - filtersize**2.0 * SdZdxi_ft(3)%values

    ! Compute Li*Mi <-- temp5(1) and Mi*Mi <-- temp5(2)
    temp5(3)%values = Li(1)%values * Mi(1)%values + Li(2)%values * Mi(2)%values + Li(3)%values * Mi(3)%values
    temp5(2)%values = Mi(1)%values * Mi(1)%values + Mi(2)%values * Mi(2)%values + Mi(3)%values * Mi(3)%values
    ! Compute DynCoef = <LiMi> / <LiMi>
    ie_gm = computeFieldAvg( temp5(3) , spec_rank , nbcpus )
    qe_gm = computeFieldAvg( temp5(2) , spec_rank , nbcpus )
    DynCoef = ie_gm / qe_gm

    !!TO TEST!!call ftran(U,tempK(1),res)
    !!TO TEST!!call ftran(V,tempK(2),res)
    !!TO TEST!!call ftran(W,tempK(3),res)
    !!TO TEST!!call ftran(Scal,tempK(4),res)
    !!TO TEST!!call computeFilter(ScalWaveN,filtersize,tempK(1),tempKf(1),filter) !U bar
    !!TO TEST!!call computeFilter(ScalWaveN,filtersize,tempK(2),tempKf(2),filter) !V bar
    !!TO TEST!!call computeFilter(ScalWaveN,filtersize,tempK(3),tempKf(3),filter) !W bar
    !!TO TEST!!call computeFilter(ScalWaveN,filtersize,tempK(4),tempKf(4),filter) !Z bar
    !!TO TEST!!call coefDynSca(Uf,Vf,Wf,tempKf(1),tempKf(2),tempKf(3),Scalf,tempKf(4),ScalWaveN,filter,DynCoef,spec_rank,res,filterSize=filtersize)
    if (spec_rank.eq.0) print*,'DynCoef Smag: ',nd,DynCoef

    ! compute vector |Sij| * dZdxi
    TiM(1)%values = DynCoef *real(filtersize)**2.0 * ( temp5(1)%values * temp4(1)%values )
    TiM(2)%values = DynCoef *real(filtersize)**2.0 * ( temp5(1)%values * temp4(2)%values )
    TiM(3)%values = DynCoef *real(filtersize)**2.0 * ( temp5(1)%values * temp4(3)%values )

    ! compute d/dxi ( |Sij| * dZdxi )
    res = computeFieldDivergence(ScalWaveN,TiM(1),TiM(2),TiM(3),phi(1),nbcpus,spec_rank)

    dTidxi_M%values = phi(1)%values

    ! Compute the optimal estimator of Model
    res = computeEONVar( spec_rank,(/phi(1)/),& !Les variables
                       & dTidxi_DNS,&                            !Le champ à estimer
                       & (/typeDisc/),&     !Les discrétisations par variable
                       & (/bin2/),&                 !le nombre d'intervalles par variable
                       & fieldOpt=temp(1))                    !champ optimal en sortie

    ! IE of GM
    ie_gm = computeMSE(dTidxi_DNS,temp(1),spec_rank,.true.)     ! true = adim par la variance

    ! QE of GM
    qe_gm = computeMSE(dTidxi_DNS,dTidxi_M,spec_rank,.true.)     ! true = adim par la variance

    if (spec_rank.EQ.0) then
       write(80,*) nd, ie_gm, qe_gm
    end if

    ! Dissipation term: !!! warning !!! Dissipation is defined as zdTidxi here

    ! DNS dissipation <- temp(2)
    temp(2)%values = Scalf%values*dTidxi_DNS%values
    ! model parameter
    phidiss(1)%values = Scalf%values*phi(1)%values
    ! optimal estimator dissipation <- temp(1)
    res = computeEONVar( spec_rank,(/phidiss(1)/),& !Les variables
                       & temp(2),&                            !Le champ à estimer
                       & (/typeDisc/),&     !Les discrétisations par variable
                       & (/bin3/),&                 !le nombre d'intervalles par variable
                       & fieldOpt=temp(1))                    !champ optimal en sortie
    ! Model dissipation <- temp(3)
    temp(3)%values = Scalf%values*dTidxi_M%values
    temp(3)%name = "GM_Diss"
    if (.not. computeFieldPDF(nd,temp(3),bin1,spec_rank,n_iter)) then
         write (6,'()') '[ERROR] '
         return
    endif
    ! IE et QE of GM for dissipation
    ie_gm = computeMSE(temp(2),temp(1),spec_rank,.true.)     ! true = adim par la variance
    qe_gm = computeMSE(temp(2),temp(3),spec_rank,.true.)     ! true = adim par la variance
    if (spec_rank.EQ.0) then
       write(81,*) nd, ie_gm, qe_gm
    end if
    ! global averaging of dissipation: ie_gm <- DNS dissipation and qe_gm <- model dissipation
    ie_gm = computeFieldAvg( temp(2) , spec_rank , nbcpus )
    qe_gm = computeFieldAvg( temp(3) , spec_rank , nbcpus )
    if (spec_rank.EQ.0) then
       write(82,*) nd, ie_gm, qe_gm
    end if

    temp(2)%values = TiM(1)%values*temp4(1)%values + TiM(2)%values*temp4(2)%values + TiM(3)%values*temp4(3)%values
    temp(2)%name = "DynSmag_Diss"
    if (.not. computeFieldPDF(nd,temp(2),bin1,spec_rank,n_iter)) then
         write (6,'()') '[ERROR] '
         return
    endif
    qe_gm = computeFieldAvg( temp(2) , spec_rank , nbcpus )
    if (spec_rank.EQ.0) then
       write(83,*) nd, ie_gm, qe_gm
    end if

    ! Model 3 in Wang et al., Heat and Fluid Flow, 2007 (model 3)

    ! Dynamic coefficient :
    Qi(1)%values = (2.0*filtersize)**2.0 * ( S11_ft%values*dZdxi_ft(1)%values + S12_ft%values*dZdxi_ft(2)%values + S13_ft%values*dZdxi_ft(3)%values )
    Qi(2)%values = (2.0*filtersize)**2.0 * ( S12_ft%values*dZdxi_ft(1)%values + S22_ft%values*dZdxi_ft(2)%values + S23_ft%values*dZdxi_ft(3)%values )
    Qi(3)%values = (2.0*filtersize)**2.0 * ( S13_ft%values*dZdxi_ft(1)%values + S23_ft%values*dZdxi_ft(2)%values + S33_ft%values*dZdxi_ft(3)%values )
    Pi(1)%values = S11%values*temp4(1)%values + S12%values*temp4(2)%values + S13%values*temp4(3)%values
    Pi(2)%values = S12%values*temp4(1)%values + S22%values*temp4(2)%values + S23%values*temp4(3)%values
    Pi(3)%values = S13%values*temp4(1)%values + S23%values*temp4(2)%values + S33%values*temp4(3)%values

    call ftran(Pi(1),tempK(1),res)
    call ftran(Pi(2),tempK(2),res)
    call ftran(Pi(3),tempK(3),res)
    call computeFilter(ScalWaveN,2.0*filtersize,tempK(1),tempKf(1),filter)
    call computeFilter(ScalWaveN,2.0*filtersize,tempK(2),tempKf(2),filter)
    call computeFilter(ScalWaveN,2.0*filtersize,tempK(3),tempKf(3),filter)
    call btran(tempKf(1),Pi(1),res)
    call btran(tempKf(2),Pi(2),res)
    call btran(tempKf(3),Pi(3),res)
    Qi(1)%values = Qi(1)%values - (filtersize)**2.0 * Pi(1)%values
    Qi(2)%values = Qi(2)%values - (filtersize)**2.0 * Pi(2)%values
    Qi(3)%values = Qi(3)%values - (filtersize)**2.0 * Pi(3)%values


    ! Compute Coef1 (for Smag part)
    temp5(3)%values =  (Li(1)%values*Mi(1)%values + Li(2)%values*Mi(2)%values + Li(3)%values*Mi(3)%values)*(Qi(1)%values*Qi(1)%values + Qi(2)%values*Qi(2)%values + Qi(3)%values*Qi(3)%values) &
                 &  -  (Li(1)%values*Qi(1)%values + Li(2)%values*Qi(2)%values + Li(3)%values*Qi(3)%values)*(Mi(1)%values*Qi(1)%values + Mi(2)%values*Qi(2)%values + Mi(3)%values*Qi(3)%values)
    temp5(2)%values =  (Qi(1)%values*Qi(1)%values + Qi(2)%values*Qi(2)%values + Qi(3)%values*Qi(3)%values)*(Mi(1)%values*Mi(1)%values + Mi(2)%values*Mi(2)%values + Mi(3)%values*Mi(3)%values) &
                 &  -  (Mi(1)%values*Qi(1)%values + Mi(2)%values*Qi(2)%values + Mi(3)%values*Qi(3)%values)*(Mi(1)%values*Qi(1)%values + Mi(2)%values*Qi(2)%values + Mi(3)%values*Qi(3)%values)
    !ie_gm = computeFieldAvg( temp5(3) , spec_rank , nbcpus )
    !qe_gm = computeFieldAvg( temp5(2) , spec_rank , nbcpus )
    !Coef1 = ie_gm / qe_gm
    !if (spec_rank.eq.0) print*,' Wang model, coef1 =', Coef1
    temp6(1)%values = temp5(3)%values / temp5(2)%values
    ie_gm = computeFieldAvg( temp6(1) , spec_rank , nbcpus )
    if (spec_rank.eq.0) print*,' Wang model, mean of coef1 - before clipping =', ie_gm
    do k = U%zmin, U%zmax
       do j = U%ymin, U%ymax
          do i = U%xmin, U%xmax
             if (temp6(1)%values(i,j,k).lt.-0.1) temp6(1)%values(i,j,k) = -0.1
             if (temp6(1)%values(i,j,k).gt.0.1) temp6(1)%values(i,j,k) = 0.1
          end do
       end do
    end do
    ie_gm = computeFieldAvg( temp6(1) , spec_rank , nbcpus )
    if (spec_rank.eq.0) print*,' Wang model, mean of coef1 - after clipping =', ie_gm

    ! Compute Coef2 (for other part)
    temp5(3)%values =  (Li(1)%values*Mi(1)%values + Li(2)%values*Mi(2)%values + Li(3)%values*Mi(3)%values)*(Mi(1)%values*Qi(1)%values + Mi(2)%values*Qi(2)%values + Mi(3)%values*Qi(3)%values) &
                 &  -  (Li(1)%values*Qi(1)%values + Li(2)%values*Qi(2)%values + Li(3)%values*Qi(3)%values)*(Mi(1)%values*Mi(1)%values + Mi(2)%values*Mi(2)%values + Mi(3)%values*Mi(3)%values)
    temp5(2)%values =  (Mi(1)%values*Qi(1)%values + Mi(2)%values*Qi(2)%values + Mi(3)%values*Qi(3)%values)*(Qi(1)%values*Mi(1)%values + Qi(2)%values*Mi(2)%values + Qi(3)%values*Mi(3)%values) &
                 &  -  (Mi(1)%values*Mi(1)%values + Mi(2)%values*Mi(2)%values + Mi(3)%values*Mi(3)%values)*(Qi(1)%values*Qi(1)%values + Qi(2)%values*Qi(2)%values + Qi(3)%values*Qi(3)%values)
    !ie_gm = computeFieldAvg( temp5(3) , spec_rank , nbcpus )
    !qe_gm = computeFieldAvg( temp5(2) , spec_rank , nbcpus )
    !Coef2 = ie_gm / qe_gm
    !if (spec_rank.eq.0) print*,' Wang model, coef2 =', Coef2
    temp6(2)%values = temp5(3)%values / temp5(2)%values
    ie_gm = computeFieldAvg( temp6(2) , spec_rank , nbcpus )
    if (spec_rank.eq.0) print*,' Wang model, mean of coef2 - before clipping =', ie_gm
    do k = U%zmin, U%zmax
       do j = U%ymin, U%ymax
          do i = U%xmin, U%xmax
             if (temp6(2)%values(i,j,k).lt.-0.1) temp6(2)%values(i,j,k) = -0.1
             if (temp6(2)%values(i,j,k).gt.0.1) temp6(2)%values(i,j,k) = 0.1
          end do
       end do
    end do
    ie_gm = computeFieldAvg( temp6(2) , spec_rank , nbcpus )
    if (spec_rank.eq.0) print*,' Wang model, mean of coef1 - after clipping =', ie_gm


    Pi(1)%values = S11%values*temp4(1)%values + S12%values*temp4(2)%values + S13%values*temp4(3)%values
    Pi(2)%values = S12%values*temp4(1)%values + S22%values*temp4(2)%values + S23%values*temp4(3)%values
    Pi(3)%values = S13%values*temp4(1)%values + S23%values*temp4(2)%values + S33%values*temp4(3)%values
    Pi(1)%values = temp6(2)%values * (filtersize)**2.0 * Pi(1)%values
    Pi(2)%values = temp6(2)%values * (filtersize)**2.0 * Pi(2)%values
    Pi(3)%values = temp6(2)%values * (filtersize)**2.0 * Pi(3)%values
    res = computeFieldDivergence(ScalWaveN,Pi(1),Pi(2),Pi(3),phi(1),nbcpus,spec_rank)

    TiM(1)%values = temp6(1)%values *real(filtersize)**2.0 * ( temp5(1)%values * temp4(1)%values )
    TiM(2)%values = temp6(1)%values *real(filtersize)**2.0 * ( temp5(1)%values * temp4(2)%values )
    TiM(3)%values = temp6(1)%values *real(filtersize)**2.0 * ( temp5(1)%values * temp4(3)%values )
    res = computeFieldDivergence(ScalWaveN,TiM(1),TiM(2),TiM(3),phi(2),nbcpus,spec_rank)

    TiM(1)%values = TiM(1)%values + Pi(1)%values
    TiM(2)%values = TiM(2)%values + Pi(2)%values
    TiM(3)%values = TiM(3)%values + Pi(3)%values
    res = computeFieldDivergence(ScalWaveN,TiM(1),TiM(2),TiM(3),dTidxi_M,nbcpus,spec_rank)

   ! Compute the optimal estimator of Model
    res = computeEONVar( spec_rank,(/phi(1),phi(2)/),& !Les variables
                       & dTidxi_DNS,&                            !Le champ à estimer
                       & (/typeDisc,typeDisc/),&     !Les discrétisations par variable
                       & (/bin1,bin2/),&                 !le nombre d'intervalles par variable
                       & fieldOpt=temp(1))                    !champ optimal en sortie

    ! IE of GM
    ie_gm = computeMSE(dTidxi_DNS,temp(1),spec_rank,.true.)     ! true = adim par la variance

    ! QE of GM
    qe_gm = computeMSE(dTidxi_DNS,dTidxi_M,spec_rank,.true.)     ! true = adim par la variance

    if (spec_rank.EQ.0) then
       write(110,*) nd, ie_gm, qe_gm
    end if

    ! DNS dissipation <- temp(2)
    temp(2)%values = Scalf%values*dTidxi_DNS%values
    ! model parameter
    phidiss(1)%values = Scalf%values*phi(1)%values
    phidiss(2)%values = Scalf%values*phi(2)%values
    ! optimal estimator dissipation <- temp(1)
    res = computeEONVar( spec_rank,(/phidiss(1),phidiss(2)/),& !Les variables
                       & temp(2),&                            !Le champ à estimer
                       & (/typeDisc,typeDisc/),&     !Les discrétisations par variable
                       & (/bin1,bin2/),&                 !le nombre d'intervalles par variable
                       & fieldOpt=temp(1))                    !champ optimal en sortie
    ! Model dissipation <- temp(3)
    temp(3)%values = Scalf%values*dTidxi_M%values
    ! IE et QE of GM for dissipation
    ie_gm = computeMSE(temp(2),temp(1),spec_rank,.true.)     ! true = adim par la variance
    qe_gm = computeMSE(temp(2),temp(3),spec_rank,.true.)     ! true = adim par la variance
    if (spec_rank.EQ.0) then
       write(111,*) nd, ie_gm, qe_gm
    end if
    ! global averaging of dissipation: ie_gm <- DNS dissipation and qe_gm <- model dissipation
    ie_gm = computeFieldAvg( temp(2) , spec_rank , nbcpus )
    qe_gm = computeFieldAvg( temp(3) , spec_rank , nbcpus )
    if (spec_rank.EQ.0) then
       write(112,*) nd, ie_gm, qe_gm
    end if










    !!!!!!!!!!!!!!!!!
    !-> IE and QE for the decompose gradient model: dTidxi = delta^2 / 12 * ( d/dxi ( Sij^+ * dZdxj ) + d/dxi ( Sij^- * dZdxj ) + d/dxi ( Oij * dZdxj) )
    !   => Phi1 = d/dxi ( Sij^+ * dZdxj )
    !   => Phi2 = d/dxi ( Sij^- * dZdxj )
    !   => Phi3 = d/dxi ( Oij * dZdxj )
    !!!!!!!!!!!!!!!!!

    ! compute Phi3
    temp(1)%values =                                O12%values*temp4(2)%values + O13%values*temp4(3)%values
    temp(2)%values = - O12%values*temp4(1)%values                              + O23%values*temp4(3)%values
    temp(3)%values = - O13%values*temp4(1)%values - O23%values*temp4(2)%values
    res = computeFieldDivergence(ScalWaveN,temp(1),temp(2),temp(3),phi(3),nbcpus,spec_rank)

    ! compute Sij+ and Sij- and Phi1 and Phi2
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

             temp(1)%values(i,j,k) = ( MIN(valR(1),0.0_WP)* ( M1(1,1)*temp4(1)%values(i,j,k) + M1(1,2)*temp4(2)%values(i,j,k) + M1(1,3)*temp4(3)%values(i,j,k)) &
                                   & + MIN(valR(2),0.0_WP)* ( M2(1,1)*temp4(1)%values(i,j,k) + M2(1,2)*temp4(2)%values(i,j,k) + M2(1,3)*temp4(3)%values(i,j,k)) &
                                   & + MIN(valR(3),0.0_WP)* ( M3(1,1)*temp4(1)%values(i,j,k) + M3(1,2)*temp4(2)%values(i,j,k) + M3(1,3)*temp4(3)%values(i,j,k)))

             temp(2)%values(i,j,k) = ( MIN(valR(1),0.0_WP)* ( M1(2,1)*temp4(1)%values(i,j,k) + M1(2,2)*temp4(2)%values(i,j,k) + M1(2,3)*temp4(3)%values(i,j,k)) &
                                   & + MIN(valR(2),0.0_WP)* ( M2(2,1)*temp4(1)%values(i,j,k) + M2(2,2)*temp4(2)%values(i,j,k) + M2(2,3)*temp4(3)%values(i,j,k)) &
                                   & + MIN(valR(3),0.0_WP)* ( M3(2,1)*temp4(1)%values(i,j,k) + M3(2,2)*temp4(2)%values(i,j,k) + M3(2,3)*temp4(3)%values(i,j,k)))

             temp(3)%values(i,j,k) = ( MIN(valR(1),0.0_WP)* ( M1(3,1)*temp4(1)%values(i,j,k) + M1(3,2)*temp4(2)%values(i,j,k) + M1(3,3)*temp4(3)%values(i,j,k)) &
                                   & + MIN(valR(2),0.0_WP)* ( M2(3,1)*temp4(1)%values(i,j,k) + M2(3,2)*temp4(2)%values(i,j,k) + M2(3,3)*temp4(3)%values(i,j,k)) &
                                   & + MIN(valR(3),0.0_WP)* ( M3(3,1)*temp4(1)%values(i,j,k) + M3(3,2)*temp4(2)%values(i,j,k) + M3(3,3)*temp4(3)%values(i,j,k)))

             temp5(1)%values(i,j,k) = ( MAX(valR(1),0.0_WP)* ( M1(1,1)*temp4(1)%values(i,j,k) + M1(1,2)*temp4(2)%values(i,j,k) + M1(1,3)*temp4(3)%values(i,j,k)) &
                                    & + MAX(valR(2),0.0_WP)* ( M2(1,1)*temp4(1)%values(i,j,k) + M2(1,2)*temp4(2)%values(i,j,k) + M2(1,3)*temp4(3)%values(i,j,k)) &
                                    & + MAX(valR(3),0.0_WP)* ( M3(1,1)*temp4(1)%values(i,j,k) + M3(1,2)*temp4(2)%values(i,j,k) + M3(1,3)*temp4(3)%values(i,j,k)))

             temp5(2)%values(i,j,k) = ( MAX(valR(1),0.0_WP)* ( M1(2,1)*temp4(1)%values(i,j,k) + M1(2,2)*temp4(2)%values(i,j,k) + M1(2,3)*temp4(3)%values(i,j,k)) &
                                    & + MAX(valR(2),0.0_WP)* ( M2(2,1)*temp4(1)%values(i,j,k) + M2(2,2)*temp4(2)%values(i,j,k) + M2(2,3)*temp4(3)%values(i,j,k)) &
                                    & + MAX(valR(3),0.0_WP)* ( M3(2,1)*temp4(1)%values(i,j,k) + M3(2,2)*temp4(2)%values(i,j,k) + M3(2,3)*temp4(3)%values(i,j,k)))

             temp5(3)%values(i,j,k) = ( MAX(valR(1),0.0_WP)* ( M1(3,1)*temp4(1)%values(i,j,k) + M1(3,2)*temp4(2)%values(i,j,k) + M1(3,3)*temp4(3)%values(i,j,k)) &
                                    & + MAX(valR(2),0.0_WP)* ( M2(3,1)*temp4(1)%values(i,j,k) + M2(3,2)*temp4(2)%values(i,j,k) + M2(3,3)*temp4(3)%values(i,j,k)) &
                                    & + MAX(valR(3),0.0_WP)* ( M3(3,1)*temp4(1)%values(i,j,k) + M3(3,2)*temp4(2)%values(i,j,k) + M3(3,3)*temp4(3)%values(i,j,k)))

          end do
       end do
    end do
#endif

    SijmTj(1)%values = temp(1)%values
    SijmTj(2)%values = temp(2)%values
    SijmTj(3)%values = temp(3)%values
    res = computeFieldDivergence(ScalWaveN,temp(1),temp(2),temp(3),phi(2),nbcpus,spec_rank)
    res = computeFieldDivergence(ScalWaveN,temp5(1),temp5(2),temp5(3),phi(1),nbcpus,spec_rank)

    dTidxi_M%values = real(filtersize)**2.0 / 12.0 * (phi(1)%values + phi(2)%values + phi(3)%values)

    ! Compute the optimal estimator of Model
    res = computeEONVar( spec_rank,(/phi(1),phi(2),phi(3)/),& !Les variables
                       & dTidxi_DNS,&                            !Le champ à estimer
                       & (/typeDisc,typeDisc,typeDisc/),&     !Les discrétisations par variable
                       & (/bin1,bin2,bin3/),&                 !le nombre d'intervalles par variable
                       & fieldOpt=temp(1))                    !champ optimal en sortie

    ! IE of GM
    ie_gm = computeMSE(dTidxi_DNS,temp(1),spec_rank,.true.)     ! true = adim par la variance

    ! QE of GM
    qe_gm = computeMSE(dTidxi_DNS,dTidxi_M,spec_rank,.true.)     ! true = adim par la variance

    if (spec_rank.EQ.0) then
       write(20,*) nd, ie_gm, qe_gm
    end if

    ! Dissipation term: !!! warning !!! Dissipation is defined as zdTidxi here

    ! DNS dissipation <- temp(2)
    temp(2)%values = Scalf%values*dTidxi_DNS%values
    ! model parameter
    phidiss(1)%values = Scalf%values*phi(1)%values
    phidiss(2)%values = Scalf%values*phi(2)%values
    phidiss(3)%values = Scalf%values*phi(3)%values
    ! optimal estimator dissipation <- temp(1)
    res = computeEONVar( spec_rank,(/phidiss(1),phidiss(2),phidiss(3)/),& !Les variables
                       & temp(2),&                            !Le champ à estimer
                       & (/typeDisc,typeDisc,typeDisc/),&     !Les discrétisations par variable
                       & (/bin1,bin2,bin3/),&                 !le nombre d'intervalles par variable
                       & fieldOpt=temp(1))                    !champ optimal en sortie
    ! Model dissipation <- temp(3)
    temp(3)%values = Scalf%values*dTidxi_M%values
    ! IE et QE of GM for dissipation
    ie_gm = computeMSE(temp(2),temp(1),spec_rank,.true.)     ! true = adim par la variance
    qe_gm = computeMSE(temp(2),temp(3),spec_rank,.true.)     ! true = adim par la variance
    if (spec_rank.EQ.0) then
       write(21,*) nd, ie_gm, qe_gm
    end if
    ! global averaging of dissipation: ie_gm <- DNS dissipation and qe_gm <- model dissipation
    ie_gm = computeFieldAvg( temp(2) , spec_rank , nbcpus )
    qe_gm = computeFieldAvg( temp(3) , spec_rank , nbcpus )
    if (spec_rank.EQ.0) then
       write(22,*) nd, ie_gm, qe_gm
    end if

    !!!!!!!!!!!!!!!!!
    !-> IE and QE for the regularized gradient model: dTidxi = delta^2 / 12 * ( d/dxi ( Sij^- * dZdxj ) + d/dxi ( Oij * dZdxj) )
    !   => Phi1 = d/dxi ( Sij^- * dZdxj )
    !   => Phi2 = d/dxi ( Oij * dZdxj )
    !!!!!!!!!!!!!!!!!
    dTidxi_M%values = real(filtersize)**2.0 / 12.0 * (phi(2)%values + phi(3)%values)

    ! Compute the optimal estimator of Model
    res = computeEONVar( spec_rank,(/phi(2),phi(3)/),& !Les variables
                       & dTidxi_DNS,&                            !Le champ à estimer
                       & (/typeDisc,typeDisc/),&     !Les discrétisations par variable
                       & (/bin2,bin3/),&                 !le nombre d'intervalles par variable
                       & fieldOpt=temp(1))                    !champ optimal en sortie

    ! IE of GM
    ie_gm = computeMSE(dTidxi_DNS,temp(1),spec_rank,.true.)     ! true = adim par la variance

    ! QE of GM
    qe_gm = computeMSE(dTidxi_DNS,dTidxi_M,spec_rank,.true.)     ! true = adim par la variance

    if (spec_rank.EQ.0) then
       write(30,*) nd, ie_gm, qe_gm
    end if

    ! Dissipation term: !!! warning !!! Dissipation is defined as zdTidxi here

    ! DNS dissipation <- temp(2)
    temp(2)%values = Scalf%values*dTidxi_DNS%values
    ! model parameter -- already computed
    ! optimal estimator dissipation <- temp(1)
    res = computeEONVar( spec_rank,(/phidiss(2),phidiss(3)/),& !Les variables
                       & temp(2),&                            !Le champ à estimer
                       & (/typeDisc,typeDisc/),&     !Les discrétisations par variable
                       & (/bin2,bin3/),&                 !le nombre d'intervalles par variable
                       & fieldOpt=temp(1))                    !champ optimal en sortie
    ! Model dissipation <- temp(3)
    temp(3)%values = Scalf%values*dTidxi_M%values
    ! IE et QE of GM for dissipation
    ie_gm = computeMSE(temp(2),temp(1),spec_rank,.true.)     ! true = adim par la variance
    qe_gm = computeMSE(temp(2),temp(3),spec_rank,.true.)     ! true = adim par la variance
    if (spec_rank.EQ.0) then
       write(31,*) nd, ie_gm, qe_gm
    end if
    ! global averaging of dissipation: ie_gm <- DNS dissipation and qe_gm <- model dissipation
    ie_gm = computeFieldAvg( temp(2) , spec_rank , nbcpus )
    qe_gm = computeFieldAvg( temp(3) , spec_rank , nbcpus )
    if (spec_rank.EQ.0) then
       write(32,*) nd, ie_gm, qe_gm
    end if

    !!!!!!!!!!!!!!!!!
    !-> IE and QE for the regularized gradient model: dTidxi = delta^2 / 12 * ( d/dxi ( Sij^+ * dZdxj ) + d/dxi ( Oij * dZdxj) )
    !   => Phi1 = d/dxi ( Sij^+ * dZdxj )
    !   => Phi2 = d/dxi ( Oij * dZdxj )
    !!!!!!!!!!!!!!!!!
    dTidxi_M%values = real(filtersize)**2.0 / 12.0 * (phi(1)%values + phi(3)%values)

    ! Compute the optimal estimator of Model
    res = computeEONVar( spec_rank,(/phi(1),phi(3)/),& !Les variables
                       & dTidxi_DNS,&                            !Le champ à estimer
                       & (/typeDisc,typeDisc/),&     !Les discrétisations par variable
                       & (/bin2,bin3/),&                 !le nombre d'intervalles par variable
                       & fieldOpt=temp(1))                    !champ optimal en sortie

    ! IE of GM
    ie_gm = computeMSE(dTidxi_DNS,temp(1),spec_rank,.true.)     ! true = adim par la variance

    ! QE of GM
    qe_gm = computeMSE(dTidxi_DNS,dTidxi_M,spec_rank,.true.)     ! true = adim par la variance

    if (spec_rank.EQ.0) then
       write(40,*) nd, ie_gm, qe_gm
    end if

    ! Dissipation term: !!! warning !!! Dissipation is defined as zdTidxi here

    ! DNS dissipation <- temp(2)
    temp(2)%values = Scalf%values*dTidxi_DNS%values
    ! model parameter -- already computed
    ! optimal estimator dissipation <- temp(1)
    res = computeEONVar( spec_rank,(/phidiss(1),phidiss(3)/),& !Les variables
                       & temp(2),&                            !Le champ à estimer
                       & (/typeDisc,typeDisc/),&     !Les discrétisations par variable
                       & (/bin2,bin3/),&                 !le nombre d'intervalles par variable
                       & fieldOpt=temp(1))                    !champ optimal en sortie
    ! Model dissipation <- temp(3)
    temp(3)%values = Scalf%values*dTidxi_M%values
    ! IE et QE of GM for dissipation
    ie_gm = computeMSE(temp(2),temp(1),spec_rank,.true.)     ! true = adim par la variance
    qe_gm = computeMSE(temp(2),temp(3),spec_rank,.true.)     ! true = adim par la variance
    if (spec_rank.EQ.0) then
       write(41,*) nd, ie_gm, qe_gm
    end if
    ! global averaging of dissipation: ie_gm <- DNS dissipation and qe_gm <- model dissipation
    ie_gm = computeFieldAvg( temp(2) , spec_rank , nbcpus )
    qe_gm = computeFieldAvg( temp(3) , spec_rank , nbcpus )
    if (spec_rank.EQ.0) then
       write(42,*) nd, ie_gm, qe_gm
    end if

    !!!!!!!!!!!!!!!!!
    !-> IE and QE for the regularized gradient model: dTidxi = delta^2 / 12 * ( d/dxi ( Sij^- * dZdxj ) )
    !   => Phi2 = d/dxi ( Sij^- * dZdxj )
    !!!!!!!!!!!!!!!!!

    ! Dynamic procedure 4 (cf. RGM-(VII) ) ::
    ! C * (2.0*filtersize)^2 * Sij^-_ft * dZdxj_ft * dZdxi_ft = Li * dZdxi_ft
    !
    !  C =  <  Li * dZdxi_ft  >  /  <  (2.0*filtersize)^2 * Sij^-_ft * dZdxj_ft * dZdxi_ft  >
    !
    ! Li * dZdxi_ft   <-- Li(1)
    ! (2.0*filtersize)^2 * Sij^-_ft * dZdxj_ft * dZdxi_ft    <-- Mi(1)
    !
    !     C = < Li(1) > / < Mi(1) >
    !
    ! Li(1) computation
    Li(1)%values = Li(1)%values*dZdxi_ft(1)%values + Li(2)%values*dZdxi_ft(2)%values + Li(3)%values*dZdxi_ft(3)%values

    ! Mi(1) computation
    ! Sij^-_ft * dZdxj_ft <-- temp
    call ftran(S11,tempK(1),res)
    call ftran(S12,tempK(2),res)
    call ftran(S13,tempK(3),res)
    call computeFilter(ScalWaveN,2.0*filtersize,tempK(1),tempKf(1),filter) !U bar
    call computeFilter(ScalWaveN,2.0*filtersize,tempK(2),tempKf(2),filter) !V bar
    call computeFilter(ScalWaveN,2.0*filtersize,tempK(3),tempKf(3),filter) !W bar
    call btran(tempKf(1),S11_ft,res) !U bar
    call btran(tempKf(2),S12_ft,res) !V bar
    call btran(tempKf(3),S13_ft,res) !W bar
    call ftran(S22,tempK(1),res)
    call ftran(S23,tempK(2),res)
    call ftran(S33,tempK(3),res)
    call computeFilter(ScalWaveN,2.0*filtersize,tempK(1),tempKf(1),filter) !U bar
    call computeFilter(ScalWaveN,2.0*filtersize,tempK(2),tempKf(2),filter) !V bar
    call computeFilter(ScalWaveN,2.0*filtersize,tempK(3),tempKf(3),filter) !W bar
    call btran(tempKf(1),S22_ft,res) !U bar
    call btran(tempKf(2),S23_ft,res) !V bar
    call btran(tempKf(3),S33_ft,res) !W bar
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
             tab(1,1) = S11_ft%values(i,j,k)
             tab(1,2) = S12_ft%values(i,j,k)
             tab(1,3) = S13_ft%values(i,j,k)
             tab(2,1) = S12_ft%values(i,j,k)
             tab(2,2) = S22_ft%values(i,j,k)
             tab(2,3) = S23_ft%values(i,j,k)
             tab(3,1) = S13_ft%values(i,j,k)
             tab(3,2) = S23_ft%values(i,j,k)
             tab(3,3) = S33_ft%values(i,j,k)
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
             temp(1)%values(i,j,k) = ( MIN(valR(1),0.0_WP)* ( M1(1,1)*dZdxi_ft(1)%values(i,j,k) + M1(1,2)*dZdxi_ft(2)%values(i,j,k) + M1(1,3)*dZdxi_ft(3)%values(i,j,k)) &
                                   & + MIN(valR(2),0.0_WP)* ( M2(1,1)*dZdxi_ft(1)%values(i,j,k) + M2(1,2)*dZdxi_ft(2)%values(i,j,k) + M2(1,3)*dZdxi_ft(3)%values(i,j,k)) &
                                   & + MIN(valR(3),0.0_WP)* ( M3(1,1)*dZdxi_ft(1)%values(i,j,k) + M3(1,2)*dZdxi_ft(2)%values(i,j,k) + M3(1,3)*dZdxi_ft(3)%values(i,j,k)))

             temp(2)%values(i,j,k) = ( MIN(valR(1),0.0_WP)* ( M1(2,1)*dZdxi_ft(1)%values(i,j,k) + M1(2,2)*dZdxi_ft(2)%values(i,j,k) + M1(2,3)*dZdxi_ft(3)%values(i,j,k)) &
                                   & + MIN(valR(2),0.0_WP)* ( M2(2,1)*dZdxi_ft(1)%values(i,j,k) + M2(2,2)*dZdxi_ft(2)%values(i,j,k) + M2(2,3)*dZdxi_ft(3)%values(i,j,k)) &
                                   & + MIN(valR(3),0.0_WP)* ( M3(2,1)*dZdxi_ft(1)%values(i,j,k) + M3(2,2)*dZdxi_ft(2)%values(i,j,k) + M3(2,3)*dZdxi_ft(3)%values(i,j,k)))

             temp(3)%values(i,j,k) = ( MIN(valR(1),0.0_WP)* ( M1(3,1)*dZdxi_ft(1)%values(i,j,k) + M1(3,2)*dZdxi_ft(2)%values(i,j,k) + M1(3,3)*dZdxi_ft(3)%values(i,j,k)) &
                                   & + MIN(valR(2),0.0_WP)* ( M2(3,1)*dZdxi_ft(1)%values(i,j,k) + M2(3,2)*dZdxi_ft(2)%values(i,j,k) + M2(3,3)*dZdxi_ft(3)%values(i,j,k)) &
                                   & + MIN(valR(3),0.0_WP)* ( M3(3,1)*dZdxi_ft(1)%values(i,j,k) + M3(3,2)*dZdxi_ft(2)%values(i,j,k) + M3(3,3)*dZdxi_ft(3)%values(i,j,k)))

          end do
       end do
    end do
#endif

    Mi(1)%values = (2.0*filtersize)**2.0 * ( temp(1)%values * dZdxi_ft(1)%values + &
                                             temp(2)%values * dZdxi_ft(2)%values + &
                                             temp(3)%values * dZdxi_ft(3)%values )

    ! Compute DynCoef = <Li(1)> / <Mi(1)>
    ie_gm = computeFieldAvg( Li(1) , spec_rank , nbcpus )
    qe_gm = computeFieldAvg( Mi(1) , spec_rank , nbcpus )
    DynCoef = ie_gm / qe_gm
    if (spec_rank.eq.0) print*,'DynCoef Sij^- (procedure 4): ',nd,DynCoef, 1.0/12.0

    dTidxi_M%values = DynCoef * real(filtersize)**2.0 * (phi(2)%values)

    TiM(1)%values = DynCoef * real(filtersize)**2.0 * SijmTj(1)%values
    TiM(2)%values = DynCoef * real(filtersize)**2.0 * SijmTj(2)%values
    TiM(3)%values = DynCoef * real(filtersize)**2.0 * SijmTj(3)%values

    ! Compute the optimal estimator of Model
    res = computeEONVar( spec_rank,(/phi(2)/),& !Les variables
                       & dTidxi_DNS,&                            !Le champ à estimer
                       & (/typeDisc/),&     !Les discrétisations par variable
                       & (/bin2/),&                 !le nombre d'intervalles par variable
                       & fieldOpt=temp(1))                    !champ optimal en sortie

    ! IE of GM
    ie_gm = computeMSE(dTidxi_DNS,temp(1),spec_rank,.true.)     ! true = adim par la variance

    ! QE of GM
    qe_gm = computeMSE(dTidxi_DNS,dTidxi_M,spec_rank,.true.)     ! true = adim par la variance

    if (spec_rank.EQ.0) then
       write(50,*) nd, ie_gm, qe_gm
    end if

    ! Dissipation term: !!! warning !!! Dissipation is defined as zdTidxi here

    ! DNS dissipation <- temp(2)
    temp(2)%values = Scalf%values*dTidxi_DNS%values
    ! model parameter -- already computed
    ! optimal estimator dissipation <- temp(1)
    res = computeEONVar( spec_rank,(/phidiss(2)/),& !Les variables
                       & temp(2),&                            !Le champ à estimer
                       & (/typeDisc/),&     !Les discrétisations par variable
                       & (/bin3/),&                 !le nombre d'intervalles par variable
                       & fieldOpt=temp(1))                    !champ optimal en sortie
    ! Model dissipation <- temp(3)
    temp(3)%values = Scalf%values*dTidxi_M%values
    ! IE et QE of GM for dissipation
    ie_gm = computeMSE(temp(2),temp(1),spec_rank,.true.)     ! true = adim par la variance
    qe_gm = computeMSE(temp(2),temp(3),spec_rank,.true.)     ! true = adim par la variance
    if (spec_rank.EQ.0) then
       write(51,*) nd, ie_gm, qe_gm
    end if
    ! global averaging of dissipation: ie_gm <- DNS dissipation and qe_gm <- model dissipation
    ie_gm = computeFieldAvg( temp(2) , spec_rank , nbcpus )
    qe_gm = computeFieldAvg( temp(3) , spec_rank , nbcpus )
    if (spec_rank.EQ.0) then
       write(52,*) nd, ie_gm, qe_gm
    end if

    temp(2)%values = TiM(1)%values*temp4(1)%values + TiM(2)%values*temp4(2)%values + TiM(3)%values*temp4(3)%values
    temp(2)%name = "S-2_Diss"
    if (.not. computeFieldPDF(nd,temp(2),bin1,spec_rank,n_iter)) then
         write (6,'()') '[ERROR] '
         return
    endif
    qe_gm = computeFieldAvg( temp(2) , spec_rank , nbcpus )
    if (spec_rank.EQ.0) then
       write(53,*) nd, ie_gm, qe_gm
    end if

    !!!!!!!!!!!!!!!!!
    !-> IE and QE for the regularized gradient model: dTidxi = delta^2 / 12 * ( d/dxi ( Sij^- * dZdxj ) + d/dxi ( Oij * dZdxj) )
    !   => Phi1 = d/dxi ( Sij^- * dZdxj ) + d/dxi ( Oij * dZdxj )
    !!!!!!!!!!!!!!!!!
    dTidxi_M%values = real(filtersize)**2.0 / 12.0 * (phi(2)%values + phi(3)%values)

    ! Compute the optimal estimator of Model
    phi(4)%values = phi(2)%values + phi(3)%values
    res = computeEONVar( spec_rank,(/phi(4)/),& !Les variables
                       & dTidxi_DNS,&                            !Le champ à estimer
                       & (/typeDisc/),&     !Les discrétisations par variable
                       & (/bin3/),&                 !le nombre d'intervalles par variable
                       & fieldOpt=temp(1))                    !champ optimal en sortie

    ! IE of GM
    ie_gm = computeMSE(dTidxi_DNS,temp(1),spec_rank,.true.)     ! true = adim par la variance

    ! QE of GM
    qe_gm = computeMSE(dTidxi_DNS,dTidxi_M,spec_rank,.true.)     ! true = adim par la variance

    if (spec_rank.EQ.0) then
       write(60,*) nd, ie_gm, qe_gm
    end if

    ! Dissipation term: !!! warning !!! Dissipation is defined as zdTidxi here

    ! DNS dissipation <- temp(2)
    temp(2)%values = Scalf%values*dTidxi_DNS%values
    ! model parameter -- already computed
    ! optimal estimator dissipation <- temp(1)
    phidiss(4)%values = phidiss(2)%values+phidiss(3)%values
    res = computeEONVar( spec_rank,(/phidiss(4)/),& !Les variables
                       & temp(2),&                            !Le champ à estimer
                       & (/typeDisc/),&     !Les discrétisations par variable
                       & (/bin3/),&                 !le nombre d'intervalles par variable
                       & fieldOpt=temp(1))                    !champ optimal en sortie
    ! Model dissipation <- temp(3)
    temp(3)%values = Scalf%values*dTidxi_M%values
    ! IE et QE of GM for dissipation
    ie_gm = computeMSE(temp(2),temp(1),spec_rank,.true.)     ! true = adim par la variance
    qe_gm = computeMSE(temp(2),temp(3),spec_rank,.true.)     ! true = adim par la variance
    if (spec_rank.EQ.0) then
       write(61,*) nd, ie_gm, qe_gm
    end if
    ! global averaging of dissipation: ie_gm <- DNS dissipation and qe_gm <- model dissipation
    ie_gm = computeFieldAvg( temp(2) , spec_rank , nbcpus )
    qe_gm = computeFieldAvg( temp(3) , spec_rank , nbcpus )
    if (spec_rank.EQ.0) then
       write(62,*) nd, ie_gm, qe_gm
    end if

    !!!!!!!!!!!!!!!!!
    !-> IE and QE for the regularized gradient model: dTidxi = delta^2 / 12 * ( d/dxi ( Sij^+ * dZdxj ) + d/dxi ( Oij * dZdxj) )
    !   => Phi1 = d/dxi ( Sij^+ * dZdxj ) + d/dxi ( Oij * dZdxj )
    !!!!!!!!!!!!!!!!!
    dTidxi_M%values = real(filtersize)**2.0 / 12.0 * (phi(1)%values + phi(3)%values)

    ! Compute the optimal estimator of Model
    phi(4)%values = phi(1)%values + phi(3)%values
    res = computeEONVar( spec_rank,(/phi(4)/),& !Les variables
                       & dTidxi_DNS,&                            !Le champ à estimer
                       & (/typeDisc/),&     !Les discrétisations par variable
                       & (/bin3/),&                 !le nombre d'intervalles par variable
                       & fieldOpt=temp(1))                    !champ optimal en sortie

    ! IE of GM
    ie_gm = computeMSE(dTidxi_DNS,temp(1),spec_rank,.true.)     ! true = adim par la variance

    ! QE of GM
    qe_gm = computeMSE(dTidxi_DNS,dTidxi_M,spec_rank,.true.)     ! true = adim par la variance

    if (spec_rank.EQ.0) then
       write(70,*) nd, ie_gm, qe_gm
    end if

    ! Dissipation term: !!! warning !!! Dissipation is defined as zdTidxi here

    ! DNS dissipation <- temp(2)
    temp(2)%values = Scalf%values*dTidxi_DNS%values
    ! model parameter -- already computed
    ! optimal estimator dissipation <- temp(1)
    phidiss(4)%values = phidiss(1)%values+phidiss(3)%values
    res = computeEONVar( spec_rank,(/phidiss(4)/),& !Les variables
                       & temp(2),&                            !Le champ à estimer
                       & (/typeDisc/),&     !Les discrétisations par variable
                       & (/bin3/),&                 !le nombre d'intervalles par variable
                       & fieldOpt=temp(1))                    !champ optimal en sortie
    ! Model dissipation <- temp(3)
    temp(3)%values = Scalf%values*dTidxi_M%values
    ! IE et QE of GM for dissipation
    ie_gm = computeMSE(temp(2),temp(1),spec_rank,.true.)     ! true = adim par la variance
    qe_gm = computeMSE(temp(2),temp(3),spec_rank,.true.)     ! true = adim par la variance
    if (spec_rank.EQ.0) then
       write(71,*) nd, ie_gm, qe_gm
    end if
    ! global averaging of dissipation: ie_gm <- DNS dissipation and qe_gm <- model dissipation
    ie_gm = computeFieldAvg( temp(2) , spec_rank , nbcpus )
    qe_gm = computeFieldAvg( temp(3) , spec_rank , nbcpus )
    if (spec_rank.EQ.0) then
       write(72,*) nd, ie_gm, qe_gm
    end if

    !!!!!!!!!!!!!!!!!
    !-> IE and QE for the regularized gradient model: dTidxi = delta^2 / 12 * ( d/dxi ( Sij^+ * dZdxj ) )
    !   => Phi1 = d/dxi ( Sij^+ * dZdxj )
    !!!!!!!!!!!!!!!!!
    dTidxi_M%values = real(filtersize)**2.0 / 12.0 * (phi(1)%values)

    ! Compute the optimal estimator of Model
    phi(4)%values = phi(1)%values
    res = computeEONVar( spec_rank,(/phi(4)/),& !Les variables
                       & dTidxi_DNS,&                            !Le champ à estimer
                       & (/typeDisc/),&     !Les discrétisations par variable
                       & (/bin3/),&                 !le nombre d'intervalles par variable
                       & fieldOpt=temp(1))                    !champ optimal en sortie

    ! IE of GM
    ie_gm = computeMSE(dTidxi_DNS,temp(1),spec_rank,.true.)     ! true = adim par la variance

    ! QE of GM
    qe_gm = computeMSE(dTidxi_DNS,dTidxi_M,spec_rank,.true.)     ! true = adim par la variance

    if (spec_rank.EQ.0) then
       write(90,*) nd, ie_gm, qe_gm
    end if

    ! Dissipation term: !!! warning !!! Dissipation is defined as zdTidxi here

    ! DNS dissipation <- temp(2)
    temp(2)%values = Scalf%values*dTidxi_DNS%values
    ! model parameter -- already computed
    ! optimal estimator dissipation <- temp(1)
    phidiss(4)%values = phidiss(1)%values
    res = computeEONVar( spec_rank,(/phidiss(4)/),& !Les variables
                       & temp(2),&                            !Le champ à estimer
                       & (/typeDisc/),&     !Les discrétisations par variable
                       & (/bin3/),&                 !le nombre d'intervalles par variable
                       & fieldOpt=temp(1))                    !champ optimal en sortie
    ! Model dissipation <- temp(3)
    temp(3)%values = Scalf%values*dTidxi_M%values
    ! IE et QE of GM for dissipation
    ie_gm = computeMSE(temp(2),temp(1),spec_rank,.true.)     ! true = adim par la variance
    qe_gm = computeMSE(temp(2),temp(3),spec_rank,.true.)     ! true = adim par la variance
    if (spec_rank.EQ.0) then
       write(91,*) nd, ie_gm, qe_gm
    end if
    ! global averaging of dissipation: ie_gm <- DNS dissipation and qe_gm <- model dissipation
    ie_gm = computeFieldAvg( temp(2) , spec_rank , nbcpus )
    qe_gm = computeFieldAvg( temp(3) , spec_rank , nbcpus )
    if (spec_rank.EQ.0) then
       write(92,*) nd, ie_gm, qe_gm
    end if

    !!!!!!!!!!!!!!!!!
    !-> IE and QE for the regularized gradient model: dTidxi = delta^2 / 12 * ( d/dxi ( Oij * dZdxj ) )
    !   => Phi1 = d/dxi ( Oij * dZdxj )
    !!!!!!!!!!!!!!!!!
    dTidxi_M%values = real(filtersize)**2.0 / 12.0 * (phi(3)%values)

    ! Compute the optimal estimator of Model
    phi(4)%values = phi(3)%values
    res = computeEONVar( spec_rank,(/phi(4)/),& !Les variables
                       & dTidxi_DNS,&                            !Le champ à estimer
                       & (/typeDisc/),&     !Les discrétisations par variable
                       & (/bin3/),&                 !le nombre d'intervalles par variable
                       & fieldOpt=temp(1))                    !champ optimal en sortie

    ! IE of GM
    ie_gm = computeMSE(dTidxi_DNS,temp(1),spec_rank,.true.)     ! true = adim par la variance

    ! QE of GM
    qe_gm = computeMSE(dTidxi_DNS,dTidxi_M,spec_rank,.true.)     ! true = adim par la variance

    if (spec_rank.EQ.0) then
       write(100,*) nd, ie_gm, qe_gm
    end if

    ! Dissipation term: !!! warning !!! Dissipation is defined as zdTidxi here

    ! DNS dissipation <- temp(2)
    temp(2)%values = Scalf%values*dTidxi_DNS%values
    ! model parameter -- already computed
    ! optimal estimator dissipation <- temp(1)
    phidiss(4)%values = phidiss(3)%values
    res = computeEONVar( spec_rank,(/phidiss(4)/),& !Les variables
                       & temp(2),&                            !Le champ à estimer
                       & (/typeDisc/),&     !Les discrétisations par variable
                       & (/bin3/),&                 !le nombre d'intervalles par variable
                       & fieldOpt=temp(1))                    !champ optimal en sortie
    ! Model dissipation <- temp(3)
    temp(3)%values = Scalf%values*dTidxi_M%values
    ! IE et QE of GM for dissipation
    ie_gm = computeMSE(temp(2),temp(1),spec_rank,.true.)     ! true = adim par la variance
    qe_gm = computeMSE(temp(2),temp(3),spec_rank,.true.)     ! true = adim par la variance
    if (spec_rank.EQ.0) then
       write(101,*) nd, ie_gm, qe_gm
    end if
    ! global averaging of dissipation: ie_gm <- DNS dissipation and qe_gm <- model dissipation
    ie_gm = computeFieldAvg( temp(2) , spec_rank , nbcpus )
    qe_gm = computeFieldAvg( temp(3) , spec_rank , nbcpus )
    if (spec_rank.EQ.0) then
       write(102,*) nd, ie_gm, qe_gm
    end if

!    !Calcul de l'EO
!    res = computeEONVar( spec_rank,(/phi(1),phi(2),phi(3)/),& !Les variables
!                       & temp(1),&                            !Le champ à estimer
!                       & (/typeDisc,typeDisc,typeDisc/),&     !Les discrétisations par variable
!                       & (/bin1,bin2,bin3/),&                 !le nombre d'intervalles par variable
!                       & fieldOpt=temp(2))                    !champ optimal en sortie
!    !Calcul MSE
!    mse = computeMSE(temp(1),temp(2),spec_rank,.true.)     ! true = adim par la variance de temp(1)

  enddo

  if (spec_rank.EQ.0) then
     close(14)
     close(11)
     close(12)
     close(20)
     close(21)
     close(22)
     close(30)
     close(31)
     close(32)
     close(40)
     close(41)
     close(42)
     close(50)
     close(51)
     close(52)
     close(60)
     close(61)
     close(62)
     close(70)
     close(71)
     close(72)
     close(80)
     close(81)
     close(82)
     close(90)
     close(91)
     close(92)
     close(100)
     close(101)
     close(102)
     close(110)
     close(111)
     close(112)
  endif
  res = deleteWorkArray(phi)
  res = deleteWorkArray(temp)
  res = deleteWorkArray(temp5)
  res = deleteWorkArray(temp2)
  res = deleteWorkArray(temp3)
  res = deleteWorkArray(temp4)

  success = .true.

 end function gradientmodOE

#endif

end module post_scalar
!> @}

