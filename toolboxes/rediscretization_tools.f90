!USEFORTEST toolbox
!USEFORTEST avgcond 
!USEFORTEST postprocess
!USEFORTEST io
!> @addtogroup toolbox
!! @{
!------------------------------------------------------------------------------
!
! MODULE: rediscretization
!
!! DESCRIPTION:
!> This module implements the interpolation routines for velocities.
!>
!> @author
!> Guillaume Balarac et Patrick BEGOU, LEGI
!>
!> @details
!> This module implements the interpolation routines for velocities.
!! Its aim is interpolate veocities from fourier space on the scalar resolution.
!! It provides the resulting velocity in the phisical space.
!!
!> @author
!> Patrick BEGOU, LEGI
!------------------------------------------------------------------------------
MODULE rediscretization_tools

USE datalayout
USE transforms_tools
USE mpilayout_tools
USE parallel_tools

  IMPLICIT NONE

INTERFACE transdiscretization
  MODULE PROCEDURE transdiscretization_dimension
  MODULE PROCEDURE transdiscretization_layout
END INTERFACE

  CONTAINS
!------------------------------------------------------------------------------
!> @author 
!> Patrick BEGOU, LEGI
!
!>
!> @details
!> This function interpolate velocities on scala resolution.
!> Optimization could be achieved with carrying ffts and interpolation in the same step (for each direction).
!>
!> @param [in] Uk  velocity to interpolate, in fourier space
!> @param [in] scalar scalar to use for velocity interpolation (contains dimensions informations) 
!> @param [in,out] Udisc velocity interpolated, in physical space
!> @return .TRUE. if initialization is successfull.
!------------------------------------------------------------------------------
  FUNCTION discretization (Uk, Scalar, Udisc, UdiscK) result(res)
  IMPLICIT NONE
    LOGICAL ::res
    TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Uk
    TYPE(REAL_DATA_LAYOUT),    INTENT(IN) :: Scalar
    TYPE(REAL_DATA_LAYOUT), INTENT(INOUT) :: Udisc
    TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT), OPTIONAL :: UdiscK

    TYPE(COMPLEX_DATA_LAYOUT) :: buffkz, buffky, buffkx
    INTEGER :: me, nbcpus

    res=.FALSE.
    me=getnbmycpu()
    nbcpus=getnbcpus()

    !DEBUG
    !IF(me.EQ.0) WRITE(6,'(6(a,i0),a)')'[DEBUG] discretization switching from [',Uk%nx,',',Uk%ny,',',Uk%nz,&
    !                         &'] to [',Udisc%nx/2+1,',',Udisc%ny,',',Udisc%nz,']'

    !Check that Uk is alongZ (default behavior)
    IF(Uk%curlayout.NE.AlongZ) THEN
       IF(me.EQ.0) WRITE(6,'(a)')'[ERROR] discretization: internal error: velocity not provided Along Z'
       RETURN
    ENDIF

    !Check the layout of the array for output storage
    IF (.NOT. sameLayout(Scalar,Udisc)) THEN
      IF(me.EQ.0) WRITE(6,'(a)') '[ERROR] discretization: output array as not the same layout than scalar'
      RETURN
    ENDIF

    !Allocate fisrt buffer alongZ with requested dimensions on Z
    IF (.NOT. initDataLayout("buffKz",buffkz,Uk%nx,Uk%ny,scalar%nz,Uk%Lx,UK%Ly,Uk%Lz,nbcpus,me,AlongZ)) THEN
      IF(me.EQ.0) WRITE(6,'(a)') '[ERROR] discretization: buffKz allocation fails'
      RETURN
    ENDIF

    !Proceed along Z dimension
    IF(buffKz%nz .GT. Uk%nz) THEN
       buffKz%values(:,:,1:Uk%nz/2+1)                   = Uk%values(:,:,1:Uk%nz/2+1)
       buffKz%values(:,:,buffKz%nz-Uk%nz/2+2:buffKz%nz) = Uk%values(:,:,Uk%nz/2+2:Uk%nz)
    ELSE
       buffKz%values(:,:,1:buffKz%nz/2+1)               = Uk%values(:,:,1:buffKz%nz/2+1)
       buffKz%values(:,:,buffKz%nz/2+2:buffKz%nz)       = Uk%values(:,:,Uk%nz-buffKz%nz/2+2:Uk%nz)
    ENDIF

    !May be we can run a backward FFT along Z here

    ! set buffer along Y for processing next direction
    IF(.NOT.setAlongY(me,buffKz)) THEN
      CALL deleteDataLayout(buffKz)
      RETURN
    ENDIF

    !Allocate second buffer alongY with requested dimensions on Y and Z
    IF (.NOT. initDataLayout("buffKy",buffky,Uk%nx,scalar%ny,scalar%nz,Uk%Lx,UK%Ly,Uk%Lz,nbcpus,me,AlongY)) THEN
      IF(me.EQ.0) WRITE(6,'(a)') '[ERROR] discretization: buffKy allocation fails'
      CALL deleteDataLayout(buffKz)
      RETURN
    ENDIF

    !Proceed along Y dimension
    IF(buffKy%ny .GT. buffKz%ny) THEN
       buffKy%values(:,1:buffKz%ny/2+1,:)                   = buffKz%values(:,1:buffKz%ny/2+1,:)
       buffKy%values(:,buffKy%ny-buffKz%ny/2+2:buffKy%ny,:) = buffKz%values(:,buffKz%ny/2+2:buffKz%ny,:)
    ELSE
       buffKy%values(:,1:buffKy%ny/2+1,:)               = buffKz%values(:,1:buffKy%ny/2+1,:)
       buffKy%values(:,buffKy%ny/2+2:buffKy%ny,:)       = buffKz%values(:,buffKz%ny-buffKy%ny/2+2:buffKz%ny,:)
    ENDIF
    CALL deleteDataLayout(buffKz)

    !May be we can run a backward FFT along Y here

    ! set buffer along X for processing next direction
    IF(.NOT.setAlongX(me,buffKy)) THEN
      CALL deleteDataLayout(buffKy)
      RETURN
    ENDIF

    !Allocate third buffer alongx with requested dimensions on X, Y and Z
    IF (.NOT. initDataLayout("buffKx",buffkx,scalar%nx/2+1,scalar%ny,scalar%nz,Uk%Lx,UK%Ly,Uk%Lz,nbcpus,me,AlongX)) THEN
      IF(me.EQ.0) WRITE(6,'(a)') '[ERROR] discretization: buffKx allocation fails'
      CALL deleteDataLayout(buffKy)
      RETURN
    ENDIF

    !Proceed along Y dimension
    buffKx%values(1:min(buffKx%nx,buffky%nx),:,:)=buffKy%values(1:min(buffKx%nx,buffky%nx),:,:)
    CALL deleteDataLayout(buffKy)
    
    
    !May be we can run a backward FFT along X here
    
    ! set buffer along Z for FFts
    
    ! If we run FFT along this function, these lines are unneeded.
    IF(.NOT.setAlongZ(me,buffKx)) THEN
      CALL deleteDataLayout(buffKx)
      RETURN
    ELSE
      IF(PRESENT(UdiscK)) THEN
        !Check the layout of the array for output storage
        IF (.NOT. sameLayout(buffkx,UdiscK)) THEN
          IF(me.EQ.0) WRITE(6,'(a)') '[ERROR] discretization: output spectral array as not the same layout &
                                     & than scalar in spectral'
          RETURN
        ENDIF
        UdiscK=buffkx
      ENDIF
      CALL btran(buffKx,Udisc,res)
      CALL deleteDataLayout(buffKx)
    ENDIF
    RETURN
  END FUNCTION discretization

!------------------------------------------------------------------------------
!> @author 
!> Patrick BEGOU, LEGI
!
!>
!> @details
!> This function interpolate a field in spectral space to a new grid 
!>
!> @param [in] Uk field to interpolate
!> @param [in] newNx new discretization along the first direction
!> @param [in] newNy new discretization along the second direction
!> @param [in] newNz new discretization along the third direction
!> @param [in,out] Udisc field interpolated, in physical space
!> @return .TRUE. if initialization is successfull.
!------------------------------------------------------------------------------
  FUNCTION transdiscretization_dimension (Uk,newNx,newNy,newNz,Udisc) result(res)
  IMPLICIT NONE
    LOGICAL ::res
    TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN)              :: Uk
    TYPE(REAL_DATA_LAYOUT), INTENT(INOUT)              :: Udisc
    INTEGER,INTENT(IN)                                 :: newNx,newNy,newNz

    include 'fftw3.f'
  
    TYPE(COMPLEX_DATA_LAYOUT) :: buffkz, buffky, buffkx
    INTEGER :: me, nbcpus
    
    !for the FFTs
    COMPLEX(WP), DIMENSION(:,:,:), ALLOCATABLE::buff
    INTEGER(KIND=8)                           :: plan
    INTEGER, DIMENSION(3)                     :: iarray, oarray
    INTEGER                                   :: ires,xx,yy,zz
    INTEGER                                   :: nbfft1d
  
    res=.FALSE.
    me=getnbmycpu()
    nbcpus=getnbcpus()
    
    
    !Check that Uk is alongZ (default behavior)
    IF(Uk%curlayout.NE.AlongZ) THEN
       IF(me.EQ.0) WRITE(6,'(a)')'[ERROR] discretization: internal error: velocity not provided Along Z'
       RETURN
    ENDIF
    
    !Check the layout of the array for output storage
    IF (Uk%Lx.NE.Udisc%Lx.OR.Uk%Ly.NE.Udisc%Ly.OR.Uk%Lz.NE.Udisc%Lz.OR.Udisc%nx.NE.newNx.OR.Udisc%ny.NE.newNy.OR.Udisc%nz.NE.newNz) THEN
      IF(me.EQ.0) WRITE(6,'(a)') '[ERROR] discretization: output array as wrong dimensions'
      RETURN
    ENDIF
    
    !Allocate fisrt buffer alongZ with requested dimensions on Z
    IF (.NOT. initDataLayout("buffKz",buffkz,Uk%nx,Uk%ny,newNz,Uk%Lx,UK%Ly,Uk%Lz,nbcpus,me,AlongZ)) THEN
      IF(me.EQ.0) WRITE(6,'(a)') '[ERROR] discretization: buffKz allocation fails'
      RETURN
    ENDIF
    
    !Proceed along Z dimension
    IF(buffKz%nz .GT. Uk%nz) THEN
       buffKz%values(:,:,1:Uk%nz/2+1)                   = Uk%values(:,:,1:Uk%nz/2+1)
       buffKz%values(:,:,buffKz%nz-Uk%nz/2+2:buffKz%nz) = Uk%values(:,:,Uk%nz/2+2:Uk%nz)
    ELSE
       buffKz%values(:,:,1:buffKz%nz/2+1)               = Uk%values(:,:,1:buffKz%nz/2+1)
       buffKz%values(:,:,buffKz%nz/2+2:buffKz%nz)       = Uk%values(:,:,Uk%nz-buffKz%nz/2+2:Uk%nz)
    ENDIF
    
    !-----------------------------------
    !run a backward FFT along Z here
    !-----------------------------------
    ! Use a temp array to reorganize datas with contiguous position along Z
    ALLOCATE (buff(buffKz%zmin:buffKz%zmax,buffKz%xmin:buffKz%xmax,buffKz%ymin:buffKz%ymax),stat=ires)
    IF (ires .NE.0) THEN
       WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": not enought memory for buff in transdiscretization!"
       RETURN
    ENDIF
    !put datas in buff
    DO zz=buffKz%zmin, buffKz%zmax
      DO yy=buffKz%ymin, buffKz%ymax
        DO xx=buffKz%xmin, buffKz%xmax
           buff(zz,xx,yy)=buffKz%values(xx,yy,zz)
        ENDDO  
      ENDDO  
    ENDDO  
    !proceed to the first backward FFT complex-complex "in place"
    ! Sizes of the local complex array
    iarray(1)=buffKz%nz
    iarray(2)=buffKz%xmax-buffKz%xmin+1
    iarray(3)=buffKz%ymax-buffKz%ymin+1
    ! number of ffts along Z to execute
    nbfft1d=(iarray(2)*iarray(3))
    ! build the plan for FFTs along X and run FFT real -> complex
    CALL dfftw_plan_many_dft(plan,1,iarray(1),nbfft1d,buff,iarray,1,iarray(1),&
                         & buff,iarray,1,iarray(1),FFTW_BACKWARD,FFTW_ESTIMATE)
    CALL dfftw_execute_dft(plan,buff,buff)
    CALL dfftw_destroy_plan(plan)
    !reorder datas in  output%values. Swap loops for cache optimization at read time.
    DO yy=buffKz%ymin, buffKz%ymax
      DO xx=buffKz%xmin, buffKz%xmax
        DO zz=buffKz%zmin, buffKz%zmax
           buffKz%values(xx,yy,zz)=buff(zz,xx,yy)
        ENDDO  
      ENDDO  
    ENDDO  
    DEALLOCATE (buff)


    ! set buffer along Y for processing next direction
    IF(.NOT.setAlongY(me,buffKz)) THEN
      CALL deleteDataLayout(buffKz)
      RETURN
    ENDIF

    !Allocate second buffer alongY with requested dimensions on Y and Z
    IF (.NOT. initDataLayout("buffKy",buffky,Uk%nx,newNy,newNz,Uk%Lx,UK%Ly,Uk%Lz,nbcpus,me,AlongY)) THEN
      IF(me.EQ.0) WRITE(6,'(a)') '[ERROR] discretization: buffKy allocation fails'
      CALL deleteDataLayout(buffKz)
      RETURN
    ENDIF
    
    !Proceed along Y dimension
    IF(buffKy%ny .GT. buffKz%ny) THEN
       buffKy%values(:,1:buffKz%ny/2+1,:)                   = buffKz%values(:,1:buffKz%ny/2+1,:)
       buffKy%values(:,buffKy%ny-buffKz%ny/2+2:buffKy%ny,:) = buffKz%values(:,buffKz%ny/2+2:buffKz%ny,:)
    ELSE
       buffKy%values(:,1:buffKy%ny/2+1,:)               = buffKz%values(:,1:buffKy%ny/2+1,:)
       buffKy%values(:,buffKy%ny/2+2:buffKy%ny,:)       = buffKz%values(:,buffKz%ny-buffKy%ny/2+2:buffKz%ny,:)
    ENDIF
    CALL deleteDataLayout(buffKz)
    
    !-----------------------------------
    ! run a backward FFT along Y here
    !-----------------------------------
    ALLOCATE (buff(buffKy%ymin:buffKy%ymax,buffKy%xmin:buffKy%xmax,buffKy%zmin:buffKy%zmax),stat=ires)
    IF (ires .NE.0) THEN
       WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": not enought memory for buff in ftran_wrap!"
       RETURN
    ENDIF
    !put datas in buff
    DO zz=buffKy%zmin, buffKy%zmax
      DO yy=buffKy%ymin, buffKy%ymax
        DO xx=buffKy%xmin, buffKy%xmax
           buff(yy,xx,zz)=buffKy%values(xx,yy,zz)
        ENDDO  
      ENDDO  
    ENDDO  
    !proceed to the second FFT complex-complex "in place"
    ! Sizes of the local complex array
    iarray(1)=buffKy%ny
    iarray(2)=buffKy%xmax-buffKy%xmin+1
    iarray(3)=buffKy%zmax-buffKy%zmin+1
    ! number of ffts along Y to execute
    nbfft1d=(iarray(2)*iarray(3))
    ! build the plan for FFTs along X and run FFT complex -> complex
    CALL dfftw_plan_many_dft(plan,1,iarray(1),nbfft1d,buff,iarray,1,iarray(1),&
                           & buff,iarray,1,iarray(1),FFTW_BACKWARD,FFTW_ESTIMATE)
    CALL dfftw_execute_dft(plan,buff,buff)
    CALL dfftw_destroy_plan(plan)
    !reorder datas in  output%values. Swap loops for cache optimization at read time.
    DO zz=buffKy%zmin, buffKy%zmax
      DO xx=buffKy%xmin, buffKy%xmax
        DO yy=buffKy%ymin, buffKy%ymax
           buffKy%values(xx,yy,zz)=buff(yy,xx,zz)
        ENDDO  
      ENDDO  
    ENDDO  
    DEALLOCATE (buff)
    

    ! set buffer along X for processing next direction
    IF(.NOT.setAlongX(me,buffKy)) THEN
      CALL deleteDataLayout(buffKy)
      RETURN
    ENDIF

    !Allocate third buffer alongx with requested dimensions on X, Y and Z
    IF (.NOT. initDataLayout("buffKx",buffkx,newNx/2+1,newNy,newNz,Uk%Lx,UK%Ly,Uk%Lz,nbcpus,me,AlongX)) THEN
      IF(me.EQ.0) WRITE(6,'(a)') '[ERROR] discretization: buffKx allocation fails'
      CALL deleteDataLayout(buffKy)
      RETURN
    ENDIF
    
    !Proceed along Y dimension
    buffKx%values(1:min(buffKx%nx,buffky%nx),:,:)=buffKy%values(1:min(buffKx%nx,buffky%nx),:,:)

    CALL deleteDataLayout(buffKy)
    
    !-----------------------------------
    !run a backward FFT along X here
    !-----------------------------------
    ! Sizes of the local real array
    iarray(1)=buffkx%nx
    iarray(2)=buffkx%ymax-buffkx%ymin+1
    iarray(3)=buffkx%zmax-buffkx%zmin+1
    ! Sizes of the local complex array
    oarray(1)=Udisc%nx
    oarray(2)=Udisc%ymax-Udisc%ymin+1
    oarray(3)=Udisc%zmax-Udisc%zmin+1
    ! number of ffts along X to execute
    nbfft1d=(iarray(2)*iarray(3))
    ! build the plan for FFTs along X and run FFT real -> complex
    CALL dfftw_plan_many_dft_c2r(plan,1,Udisc%nx,nbfft1d,buffkx%values,iarray,1,buffkx%nx,&
                               & Udisc%values,oarray,1,Udisc%nx,FFTW_ESTIMATE)
    CALL dfftw_execute_dft_c2r(plan,buffkx%values,Udisc%values)
    CALL dfftw_destroy_plan(plan)
    
    CALL deleteDataLayout(buffKx)
    res=.TRUE.
    RETURN
  END FUNCTION transdiscretization_dimension



!------------------------------------------------------------------------------
!> @author 
!> Patrick BEGOU, LEGI
!
!>
!> @details
!> This function interpolate velocities on scala resolution.
!> Optimization could be achieved with carrying ffts and interpolation in the same step (for each direction).
!>
!> @param [in] Uk  velocity to interpolate, in fourier space
!> @param [in] scalar scalar to use for velocity interpolation (contains dimensions informations) 
!> @param [in,out] Udisc velocity interpolated, in physical space
!> @return .TRUE. if initialization is successfull.
!------------------------------------------------------------------------------
  FUNCTION transdiscretization_layout (Uk, Scalar, Udisc) result(res)
  IMPLICIT NONE
    LOGICAL ::res
    TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN)              :: Uk
    TYPE(REAL_DATA_LAYOUT),    INTENT(IN)              :: Scalar
    TYPE(REAL_DATA_LAYOUT), INTENT(INOUT)              :: Udisc

    include 'fftw3.f'
  
    TYPE(COMPLEX_DATA_LAYOUT) :: buffkz, buffky, buffkx
    INTEGER :: me, nbcpus
    
    !for the FFTs
    COMPLEX(WP), DIMENSION(:,:,:), ALLOCATABLE::buff
    INTEGER(KIND=8)                           :: plan
    INTEGER, DIMENSION(3)                     :: iarray, oarray
    INTEGER                                   :: ires,xx,yy,zz
    INTEGER                                   :: nbfft1d
  
    res=.FALSE.
    me=getnbmycpu()
    nbcpus=getnbcpus()
    
    !DEBUG
    !IF(me.EQ.0) WRITE(6,'(6(a,i0),a)')'[DEBUG] discretization switching from [',Uk%nx,',',Uk%ny,',',Uk%nz,&
    !                         &'] to [',Udisc%nx/2+1,',',Udisc%ny,',',Udisc%nz,']'
    
    !Check that Uk is alongZ (default behavior)
    IF(Uk%curlayout.NE.AlongZ) THEN
       IF(me.EQ.0) WRITE(6,'(a)')'[ERROR] discretization: internal error: velocity not provided Along Z'
       RETURN
    ENDIF
    
    !Check the layout of the array for output storage
    IF (.NOT. sameLayout(Scalar,Udisc)) THEN
      IF(me.EQ.0) WRITE(6,'(a)') '[ERROR] discretization: output array as not the same layout than scalar'
      RETURN
    ENDIF
    
    !Allocate fisrt buffer alongZ with requested dimensions on Z
    IF (.NOT. initDataLayout("buffKz",buffkz,Uk%nx,Uk%ny,scalar%nz,Uk%Lx,UK%Ly,Uk%Lz,nbcpus,me,AlongZ)) THEN
      IF(me.EQ.0) WRITE(6,'(a)') '[ERROR] discretization: buffKz allocation fails'
      RETURN
    ENDIF
    
    !Proceed along Z dimension
    IF(buffKz%nz .GT. Uk%nz) THEN
       buffKz%values(:,:,1:Uk%nz/2+1)                   = Uk%values(:,:,1:Uk%nz/2+1)
       buffKz%values(:,:,buffKz%nz-Uk%nz/2+2:buffKz%nz) = Uk%values(:,:,Uk%nz/2+2:Uk%nz)
    ELSE
       buffKz%values(:,:,1:buffKz%nz/2+1)               = Uk%values(:,:,1:buffKz%nz/2+1)
       buffKz%values(:,:,buffKz%nz/2+2:buffKz%nz)       = Uk%values(:,:,Uk%nz-buffKz%nz/2+2:Uk%nz)
    ENDIF
    
    !-----------------------------------
    !run a backward FFT along Z here
    !-----------------------------------
    ! Use a temp array to reorganize datas with contiguous position along Z
    ALLOCATE (buff(buffKz%zmin:buffKz%zmax,buffKz%xmin:buffKz%xmax,buffKz%ymin:buffKz%ymax),stat=ires)
    IF (ires .NE.0) THEN
       WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": not enought memory for buff in transdiscretization!"
       RETURN
    ENDIF
    !put datas in buff
    DO zz=buffKz%zmin, buffKz%zmax
      DO yy=buffKz%ymin, buffKz%ymax
        DO xx=buffKz%xmin, buffKz%xmax
           buff(zz,xx,yy)=buffKz%values(xx,yy,zz)
        ENDDO  
      ENDDO  
    ENDDO  
    !proceed to the first backward FFT complex-complex "in place"
    ! Sizes of the local complex array
    iarray(1)=buffKz%nz
    iarray(2)=buffKz%xmax-buffKz%xmin+1
    iarray(3)=buffKz%ymax-buffKz%ymin+1
    ! number of ffts along Z to execute
    nbfft1d=(iarray(2)*iarray(3))
    ! build the plan for FFTs along X and run FFT real -> complex
    CALL dfftw_plan_many_dft(plan,1,iarray(1),nbfft1d,buff,iarray,1,iarray(1),&
                         & buff,iarray,1,iarray(1),FFTW_BACKWARD,FFTW_ESTIMATE)
    CALL dfftw_execute_dft(plan,buff,buff)
    CALL dfftw_destroy_plan(plan)
    !reorder datas in  output%values. Swap loops for cache optimization at read time.
    DO yy=buffKz%ymin, buffKz%ymax
      DO xx=buffKz%xmin, buffKz%xmax
        DO zz=buffKz%zmin, buffKz%zmax
           buffKz%values(xx,yy,zz)=buff(zz,xx,yy)
        ENDDO  
      ENDDO  
    ENDDO  
    DEALLOCATE (buff)


    ! set buffer along Y for processing next direction
    IF(.NOT.setAlongY(me,buffKz)) THEN
      CALL deleteDataLayout(buffKz)
      RETURN
    ENDIF

    !Allocate second buffer alongY with requested dimensions on Y and Z
    IF (.NOT. initDataLayout("buffKy",buffky,Uk%nx,scalar%ny,scalar%nz,Uk%Lx,UK%Ly,Uk%Lz,nbcpus,me,AlongY)) THEN
      IF(me.EQ.0) WRITE(6,'(a)') '[ERROR] discretization: buffKy allocation fails'
      CALL deleteDataLayout(buffKz)
      RETURN
    ENDIF
    
    !Proceed along Y dimension
    IF(buffKy%ny .GT. buffKz%ny) THEN
       buffKy%values(:,1:buffKz%ny/2+1,:)                   = buffKz%values(:,1:buffKz%ny/2+1,:)
       buffKy%values(:,buffKy%ny-buffKz%ny/2+2:buffKy%ny,:) = buffKz%values(:,buffKz%ny/2+2:buffKz%ny,:)
    ELSE
       buffKy%values(:,1:buffKy%ny/2+1,:)               = buffKz%values(:,1:buffKy%ny/2+1,:)
       buffKy%values(:,buffKy%ny/2+2:buffKy%ny,:)       = buffKz%values(:,buffKz%ny-buffKy%ny/2+2:buffKz%ny,:)
    ENDIF
    CALL deleteDataLayout(buffKz)
    
    !-----------------------------------
    ! run a backward FFT along Y here
    !-----------------------------------
    ALLOCATE (buff(buffKy%ymin:buffKy%ymax,buffKy%xmin:buffKy%xmax,buffKy%zmin:buffKy%zmax),stat=ires)
    IF (ires .NE.0) THEN
       WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": not enought memory for buff in ftran_wrap!"
       RETURN
    ENDIF
    !put datas in buff
    DO zz=buffKy%zmin, buffKy%zmax
      DO yy=buffKy%ymin, buffKy%ymax
        DO xx=buffKy%xmin, buffKy%xmax
           buff(yy,xx,zz)=buffKy%values(xx,yy,zz)
        ENDDO  
      ENDDO  
    ENDDO  
    !proceed to the second FFT complex-complex "in place"
    ! Sizes of the local complex array
    iarray(1)=buffKy%ny
    iarray(2)=buffKy%xmax-buffKy%xmin+1
    iarray(3)=buffKy%zmax-buffKy%zmin+1
    ! number of ffts along Y to execute
    nbfft1d=(iarray(2)*iarray(3))
    ! build the plan for FFTs along X and run FFT complex -> complex
    CALL dfftw_plan_many_dft(plan,1,iarray(1),nbfft1d,buff,iarray,1,iarray(1),&
                           & buff,iarray,1,iarray(1),FFTW_BACKWARD,FFTW_ESTIMATE)
    CALL dfftw_execute_dft(plan,buff,buff)
    CALL dfftw_destroy_plan(plan)
    !reorder datas in  output%values. Swap loops for cache optimization at read time.
    DO zz=buffKy%zmin, buffKy%zmax
      DO xx=buffKy%xmin, buffKy%xmax
        DO yy=buffKy%ymin, buffKy%ymax
           buffKy%values(xx,yy,zz)=buff(yy,xx,zz)
        ENDDO  
      ENDDO  
    ENDDO  
    DEALLOCATE (buff)
    

    ! set buffer along X for processing next direction
    IF(.NOT.setAlongX(me,buffKy)) THEN
      CALL deleteDataLayout(buffKy)
      RETURN
    ENDIF

    !Allocate third buffer alongx with requested dimensions on X, Y and Z
    IF (.NOT. initDataLayout("buffKx",buffkx,scalar%nx/2+1,scalar%ny,scalar%nz,Uk%Lx,UK%Ly,Uk%Lz,nbcpus,me,AlongX)) THEN
      IF(me.EQ.0) WRITE(6,'(a)') '[ERROR] discretization: buffKx allocation fails'
      CALL deleteDataLayout(buffKy)
      RETURN
    ENDIF
    
    !Proceed along Y dimension
    buffKx%values(1:min(buffKx%nx,buffky%nx),:,:)=buffKy%values(1:min(buffKx%nx,buffky%nx),:,:)

    CALL deleteDataLayout(buffKy)
    
    !-----------------------------------
    !run a backward FFT along X here
    !-----------------------------------
    ! Sizes of the local real array
    iarray(1)=buffkx%nx
    iarray(2)=buffkx%ymax-buffkx%ymin+1
    iarray(3)=buffkx%zmax-buffkx%zmin+1
    ! Sizes of the local complex array
    oarray(1)=Udisc%nx
    oarray(2)=Udisc%ymax-Udisc%ymin+1
    oarray(3)=Udisc%zmax-Udisc%zmin+1
    ! number of ffts along X to execute
    nbfft1d=(iarray(2)*iarray(3))
    ! build the plan for FFTs along X and run FFT real -> complex
    CALL dfftw_plan_many_dft_c2r(plan,1,Udisc%nx,nbfft1d,buffkx%values,iarray,1,buffkx%nx,&
                               & Udisc%values,oarray,1,Udisc%nx,FFTW_ESTIMATE)
    CALL dfftw_execute_dft_c2r(plan,buffkx%values,Udisc%values)
    CALL dfftw_destroy_plan(plan)
    
    CALL deleteDataLayout(buffKx)
    res=.TRUE.
    RETURN
  END FUNCTION transdiscretization_layout

!-------------------------------------------------------------------------------------------
!> transdiscretize field at avs format in order to be able to restart a stopped simulation.
!!    @param[in]    spec_rank       = rank inside the spectral communicator
!!    @param[in]    scalArray       = array containing all the scalar fields (scalars solved 
!!                                  with pseudo-spectral methods)
!!    @param[in]    scalArray_part  = array containing all the scalar fields (scalars solved 
!!                                  with pseudo-spectral methods)
!!    @param[in]    U               = velocity along X
!!    @param[in]    V               = velocity along Y
!!    @param[in]    W               = velocity along Z
!!    @param[in]    B               = magentic field (3D vector)
!!    @param[in]    imhd            = 1 if MHD is ON, 0 if it is OFF
!!    @param[in]    iter            = number of the current iteration
!!    @param[in]    iterMax         = total number of iterations
!!    @return       success         = .true. if everything is right !
!! @author  Antoine Vollant
!-------------------------------------------------------------------------------------------
 function transdiscretizeField(fieldIn,fieldOut,spec_rank,nbcpus,nx,ny,nz) result(success)


   !I/O data
   type(real_data_layout),intent(in)                      :: fieldIn
   type(real_data_layout),intent(inout)                   :: fieldOut
   integer,intent(in)                                     :: spec_rank
   integer,intent(in)                                     :: nbcpus
   integer,intent(in)                                     :: nx,ny,nz
   logical                                                :: success 

   !Local data
   type(complex_DATA_LAYOUT)      :: tempK
   logical                        :: res 

    success = .false.

    if (.not. initDataLayout("TempK", tempK,(fieldIn%nx/2)+1,fieldIn%ny,fieldIn%nz, &
               & fieldIn%Lx,fieldIn%Ly,fieldIn%Lz,nbcpus,spec_rank,alongZ)) then
      write(6,'(a)')'[ERROR] initDataLayout for TYPE (COMPLEX_DATA_LAYOUT) : not enought memory!'
      return
    endif
    if (.not. initDataLayout(trim(fieldIn%name)//'_Trans', fieldOut ,nx,ny,nz, &
               & fieldIn%Lx,fieldIn%Ly,fieldIn%Lz,nbcpus,spec_rank)) then
      write(6,'(a)')'[ERROR] initDataLayout for TYPE (REAL_DATA_LAYOUT) : not enought memory!'
      return
    endif
    call ftran(fieldIn,tempK,res)
    if (.not. res) return
    if (.not. transdiscretization_dimension(tempK,nx,ny,nz,fieldOut) ) then
      write(6,'(a)')'[ERROR] transdiscretization failed' 
      return
    end if
    call deletedatalayout(tempK) 

    success = .true.

 end function transdiscretizeField



END MODULE rediscretization_tools
!> @}
