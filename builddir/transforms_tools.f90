!USEFORTEST toolbox
!USEFORTEST avgcond 
!USEFORTEST postprocess
!USEFORTEST io
!> @addtogroup toolbox 
!! @{
!------------------------------------------------------------------------------
!
! MODULE: transform
!
!> @author
!> Patrick Begou, LEGI
!
! DESCRIPTION: 
!> The Fast Fourier Transforms
!>
!> @details
!> The aim of this module is to provide tools for 3D FFT (forward or bakward) on the 
!> distributed domain.
!> WARNING: data are normalized by the number of point when the backward FFT returns
!> in the physical space
!------------------------------------------------------------------------------
MODULE transforms_tools
include 'fftw3.f'


CONTAINS
  !------------------------------------------------------------------------------
  !> @author 
  !> Patrick Begou
  !
  !> \brief
  !> Forward Transform
  !
  !> @details
  !> The aim of this subroutine is to run a 3D FFT on the domain. The input domain should
  !> be provided "AlongX" (it is set AlongX in the other case) and the reukts in output
  !> is organized along Z. This version is very basic without optimizations. For
  !> complex/complex it uses the "in place" algorithm (same array for input and output)
  !> Plans are evaluated for each call. The only optimization is to provide data contiguous
  !> in the FFT direction in the 3 steps (3 time a 1D FFT)
  !> @param[in] input1 the real data in the physical space
  !> @param[in,out] output the complex data in the spectral space.
  !> @param[out] res set to .TRUE. if all is successfull.
  !------------------------------------------------------------------------------
  SUBROUTINE ftran(input1,output,res)

  USE parallel_tools
  USE datalayout

  IMPLICIT NONE
  TYPE(REAL_DATA_LAYOUT), INTENT(IN)       :: input1
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT) :: output
  LOGICAL, INTENT(OUT) :: res

  !> Warning, plan is a C pointer to store in a large integer (same size as a C pointer) 
  INTEGER(KIND=8) :: plan
  INTEGER :: me,nbfft1d,iret
  INTEGER, DIMENSION(3):: iarray, oarray
  !>local array to reshape the array with contiguous values in the FFT direction
  COMPLEX(WP), DIMENSION(:,:,:), ALLOCATABLE::buff
  INTEGER ::ires,xx,yy,zz
  TYPE(REAL_DATA_LAYOUT)    ::input

  res=.FALSE.
  me=getnbmycpu()
  input=input1

  ! Check the output data structure.
  IF (.NOT.(      output%nx .EQ. (input%nx/2)+1 & 
     &      .AND. output%ny .EQ. input%ny       &
     &      .AND. output%nz .EQ. input%nz       &
     &     ) &
     &) THEN
     WRITE(6,'(a,i0,a,6(i0,a))')&
   &  '[ERROR] ftran_wrap on ',me,': incompatible data sizes!'//   &
   &  ' Real part is [',input%nx,',',input%ny,',',input%nz,        &
   &  '] and Complex is [',output%nx,',',output%ny,',',output%nz,']' 
     RETURN
  ENDIF

  !Check we are along X for Input and Output for the REAL -> complex FFT
  IF (.NOT. setAlongX(me,input)) RETURN
  IF (.NOT. setAlongX(me,output)) RETURN

  ! clear the output 
  output%values=0.0_WP

  ! Sizes of the local real array
  iarray(1)=input%nx
  iarray(2)=input%ymax-input%ymin+1
  iarray(3)=input%zmax-input%zmin+1

  ! Sizes of the local complex array
  oarray(1)=output%nx
  oarray(2)=output%ymax-input%ymin+1
  oarray(3)=output%zmax-input%zmin+1

  ! number of ffts along X to execute
  nbfft1d=(iarray(2)*iarray(3))

  ! TRY TO SPEED-UP WITH OMP
!!  call dfftw_init_threads(iret)
!!  call dfftw_plan_with_nthreads(8)


  ! build the plan for FFTs along X and run FFT real -> complex
  CALL dfftw_plan_many_dft_r2c(plan,1,input%nx,nbfft1d,input%values,iarray,1,input%nx,&
                             & output%values,oarray,1,output%nx,FFTW_ESTIMATE)
  CALL dfftw_execute_dft_r2c(plan,input%values,output%values)
  CALL dfftw_destroy_plan(plan)

  ! Now proceed along Y
  IF (.NOT. setAlongy(me,output)) RETURN

  !first use a temp array to reorganize datas with contiguous position along Y
  ALLOCATE (buff(output%ymin:output%ymax,output%xmin:output%xmax,output%zmin:output%zmax),stat=ires)
  IF (ires .NE.0) THEN
     WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": not enought memory for buff in ftran_wrap!"
     RETURN
  ENDIF

  !put datas in buff
  DO zz=output%zmin, output%zmax
    DO yy=output%ymin, output%ymax
      DO xx=output%xmin, output%xmax
         buff(yy,xx,zz)=output%values(xx,yy,zz)
      ENDDO  
   ENDDO  
  ENDDO  


  !proceed to the second FFT complex-complex "in place"
  ! Sizes of the local complex array
  iarray(1)=output%ny
  iarray(2)=output%xmax-output%xmin+1
  iarray(3)=output%zmax-output%zmin+1

  ! number of ffts along X to execute
  nbfft1d=(iarray(2)*iarray(3))

  ! build the plan for FFTs along X and run FFT complex -> complex
  ! fft along Y
  CALL dfftw_plan_many_dft(plan,1,iarray(1),nbfft1d,buff,iarray,1,iarray(1),&
                         & buff,iarray,1,iarray(1),FFTW_FORWARD,FFTW_ESTIMATE)
  CALL dfftw_execute_dft(plan,buff,buff)
  CALL dfftw_destroy_plan(plan)

  !reorder datas in  output%values. Swap loops for cache optimization at read time.
  DO zz=output%zmin, output%zmax
    DO xx=output%xmin, output%xmax
      DO yy=output%ymin, output%ymax
         output%values(xx,yy,zz)=buff(yy,xx,zz)
      ENDDO  
    ENDDO  
  ENDDO  
  DEALLOCATE (buff)

  ! Now proceed along Z
  IF (.NOT. setAlongz(me,output)) RETURN

  !first use a temp array to reorganize datas with contiguous position along Z
  ALLOCATE (buff(output%zmin:output%zmax,output%xmin:output%xmax,output%ymin:output%ymax),stat=ires)
  IF (ires .NE.0) THEN
     WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": not enought memory for buff in ftran_wrap!"
     RETURN
  ENDIF

  !put datas in buff
  DO zz=output%zmin, output%zmax
    DO yy=output%ymin, output%ymax
      DO xx=output%xmin, output%xmax
         buff(zz,xx,yy)=output%values(xx,yy,zz)
      ENDDO  
    ENDDO  
  ENDDO  

  !proceed to the third FFT complex-complex "in place"
  ! Sizes of the local complex array
  iarray(1)=output%nz
  iarray(2)=output%xmax-output%xmin+1
  iarray(3)=output%ymax-output%ymin+1

  ! number of ffts along X to execute
  nbfft1d=(iarray(2)*iarray(3))

  ! build the plan for FFTs along X and run FFT real -> complex
  ! fft along Z
  CALL dfftw_plan_many_dft(plan,1,iarray(1),nbfft1d,buff,iarray,1,iarray(1),&
                         & buff,iarray,1,iarray(1),FFTW_FORWARD,FFTW_ESTIMATE)
  CALL dfftw_execute_dft(plan,buff,buff)
  CALL dfftw_destroy_plan(plan)

  !reorder datas in  output%values. Swap loops for cache optimization at read time.
  DO yy=output%ymin, output%ymax
    DO xx=output%xmin, output%xmax
      DO zz=output%zmin, output%zmax
         output%values(xx,yy,zz)=buff(zz,xx,yy)
      ENDDO  
    ENDDO  
  ENDDO  
  DEALLOCATE (buff)

  !normalize now by the number of points
  !should prevent any overflow at hight resolution.
  output%values=output%values/REAL(input%nx,WP)
  output%values=output%values/REAL(input%ny,WP)
  output%values=output%values/REAL(input%nz,WP)

  CALL deletedataLayout(input)
  res=.TRUE.
  RETURN

  END SUBROUTINE ftran


  !------------------------------------------------------------------------------
  !> @author
  !> Patrick Begou
  !
  !> \brief
  !> Backward Transform
  !
  !> @details
  !> The aim of this subroutine is to run a 3D FFT on the domain. The input domain should
  !> be provided "AlongZ" (it is set AlongZ if it is not) and the results in output
  !> is organized along X. This version is very basic without optimizations. For
  !> complex/complex it uses the "in place" algorithm (same array for input and output)
  !> Plans are evaluated for each call. The only optimization is to provide data contiguous
  !> in the FFT direction in the 3 steps (3 time a 1D FFT)
  !> WARNING: data are normalized by the number of points when the backward FFT returns
  !> in the pisical space
  !> @param[in] input1 the complex data in the spectral space.
  !> @param[in,out] output the real data in the physical space
  !> @param[out] res set to .TRUE. if all is successfull.
  !------------------------------------------------------------------------------
  SUBROUTINE btran(input1,output,res)

  USE parallel_tools
  USE datalayout

  IMPLICIT NONE
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) ::input1
  TYPE(REAL_DATA_LAYOUT), INTENT(INOUT) ::output
  LOGICAL, INTENT(OUT) :: res

  TYPE(COMPLEX_DATA_LAYOUT)    ::input
  INTEGER(KIND=8) :: plan
  INTEGER :: me,nbfft1d
  INTEGER, DIMENSION(3):: iarray, oarray
  COMPLEX(WP), DIMENSION(:,:,:), ALLOCATABLE::buff
  INTEGER  :: ires,xx,yy,zz

  res=.FALSE.
  me=getnbmycpu()

  input=input1


  ! Check the output data structure.
  IF (.NOT.(     (output%nx/2)+1 .EQ. input%nx  & 
     &      .AND. output%ny .EQ. input%ny       &
     &      .AND. output%nz .EQ. input%nz       &
     &     ) &
     &) THEN
     WRITE(6,'(a,i0,a,6(i0,a))')&
   &  '[ERROR] Btran_wrap on ',me,': incompatible data sizes!'//   &
   &  ' Complex part is [',input%nx,',',input%ny,',',input%nz,        &
   &  '] and Real is [',output%nx,',',output%ny,',',output%nz,']' 
     RETURN
  ENDIF

  !Check  we are along Z for Input and X for output (the last FFT Complex->Real is along X)
  IF (.NOT. setAlongZ(me,input)) RETURN
  IF (.NOT. setAlongX(me,output)) RETURN

  ! clear the output 
  output%values=0.0_WP

  ! First proceed along Z
  ! Use a temp array to reorganize datas with contiguous position along Z
  ALLOCATE (buff(input%zmin:input%zmax,input%xmin:input%xmax,input%ymin:input%ymax),stat=ires)
  IF (ires .NE.0) THEN
     WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": not enought memory for buff in Btran_wrap!"
     RETURN
  ENDIF

  !put datas in buff
  DO zz=input%zmin, input%zmax
    DO yy=input%ymin, input%ymax
      DO xx=input%xmin, input%xmax
         buff(zz,xx,yy)=input%values(xx,yy,zz)
      ENDDO  
    ENDDO  
  ENDDO  

  !proceed to the first backward FFT complex-complex "in place"
  ! Sizes of the local complex array
  iarray(1)=input%nz
  iarray(2)=input%xmax-input%xmin+1
  iarray(3)=input%ymax-input%ymin+1

  ! number of ffts along X to execute
  nbfft1d=(iarray(2)*iarray(3))

  ! build the plan for FFTs along X and run FFT real -> complex
  CALL dfftw_plan_many_dft(plan,1,iarray(1),nbfft1d,buff,iarray,1,iarray(1),&
                         & buff,iarray,1,iarray(1),FFTW_BACKWARD,FFTW_ESTIMATE)
  CALL dfftw_execute_dft(plan,buff,buff)
  CALL dfftw_destroy_plan(plan)

  !reorder datas in  output%values. Swap loops for cache optimization at read time.
  DO yy=input%ymin, input%ymax
    DO xx=input%xmin, input%xmax
      DO zz=input%zmin, input%zmax
         input%values(xx,yy,zz)=buff(zz,xx,yy)
      ENDDO  
    ENDDO  
  ENDDO  
  DEALLOCATE (buff)

  ! Now proceed along Y
  IF (.NOT. setAlongy(me,input)) RETURN

  !first use a temp array to reorganize datas with contiguous position along Y
  ALLOCATE (buff(input%ymin:input%ymax,input%xmin:input%xmax,input%zmin:input%zmax),stat=ires)
  IF (ires .NE.0) THEN
     WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": not enought memory for buff in ftran_wrap!"
     RETURN
  ENDIF

  !put datas in buff
  DO zz=input%zmin, input%zmax
    DO yy=input%ymin, input%ymax
      DO xx=input%xmin, input%xmax
         buff(yy,xx,zz)=input%values(xx,yy,zz)
      ENDDO  
    ENDDO  
  ENDDO  

  !proceed to the second FFT complex-complex "in place"
  ! Sizes of the local complex array
  iarray(1)=input%ny
  iarray(2)=input%xmax-input%xmin+1
  iarray(3)=input%zmax-input%zmin+1

  ! number of ffts along X to execute
  nbfft1d=(iarray(2)*iarray(3))

  ! build the plan for FFTs along X and run FFT complex -> complex
  CALL dfftw_plan_many_dft(plan,1,iarray(1),nbfft1d,buff,iarray,1,iarray(1),&
                         & buff,iarray,1,iarray(1),FFTW_BACKWARD,FFTW_ESTIMATE)
  CALL dfftw_execute_dft(plan,buff,buff)
  CALL dfftw_destroy_plan(plan)

  !reorder datas in  output%values. Swap loops for cache optimization at read time.
  DO zz=input%zmin, input%zmax
    DO xx=input%xmin, input%xmax
      DO yy=input%ymin, input%ymax
         input%values(xx,yy,zz)=buff(yy,xx,zz)
      ENDDO  
    ENDDO  
  ENDDO  
  DEALLOCATE (buff)

  ! Run the third FFT Complex->Real along X
  IF (.NOT. setAlongX(me,input)) RETURN

  ! Sizes of the local real array
  iarray(1)=input%nx
  iarray(2)=input%ymax-input%ymin+1
  iarray(3)=input%zmax-input%zmin+1

  ! Sizes of the local complex array
  oarray(1)=output%nx
  oarray(2)=output%ymax-output%ymin+1
  oarray(3)=output%zmax-output%zmin+1

  ! number of ffts along X to execute
  nbfft1d=(iarray(2)*iarray(3))

  ! build the plan for FFTs along X and run FFT real -> complex
  CALL dfftw_plan_many_dft_c2r(plan,1,output%nx,nbfft1d,input%values,iarray,1,input%nx,&
                             & output%values,oarray,1,output%nx,FFTW_ESTIMATE)
  CALL dfftw_execute_dft_c2r(plan,input%values,output%values)
  CALL dfftw_destroy_plan(plan)

  CALL deletedataLayout(input)

  res=.TRUE.
  RETURN

  END SUBROUTINE btran


  !------------------------------------------------------------------------------
  !> Forward Transform
  !> @author
  !> Patrick Begou and small modification by Jean-Baptiste Lagaert
  !
  !
  !> @details
  !> The aim of this subroutine is to run a 3D FFT on the domain. The input domain should
  !> be provided "AlongX" (it is set AlongX in the other case) and the reukts in output
  !> is organized along Z. This version is very basic without optimizations. For
  !> complex/complex it uses the "in place" algorithm (same array for input and output)
  !> Plans are evaluated for each call. The only optimization is to provide data contiguous
  !> in the FFT direction in the 3 steps (3 time a 1D FFT)
  !> @param[in] input1 the real data in the physical space
  !> @param[in,out] output the complex data in the spectral space.
  !> @param[out] res set to .TRUE. if all is successfull.
  !------------------------------------------------------------------------------
  SUBROUTINE ftran_vector(input1,output,res)

  USE parallel_tools
  USE datalayout

  IMPLICIT NONE
  TYPE(REAL_VECTOR_DATA_LAYOUT), INTENT(INOUT)    ::input1
  TYPE(COMPLEX_VECTOR_DATA_LAYOUT), INTENT(INOUT) ::output
  LOGICAL, INTENT(OUT) :: res

  !> Warning, plan is a C pointer to store in a large integer (same size as a C pointer)
  INTEGER(KIND=8) :: plan
  INTEGER :: me,nbfft1d
  INTEGER, DIMENSION(4):: iarray, oarray
  !>local array to reshape the array with contiguous values in the FFT direction
  COMPLEX(WP), DIMENSION(:,:,:,:), ALLOCATABLE::buff
  INTEGER ::ires,xx,yy,zz,comp
  !TYPE(REAL_VECTOR_DATA_LAYOUT)    ::input

  res=.FALSE.
  me=getnbmycpu()
  !input=input1

  ! Check the output data structure.
  IF (.NOT.(      output%nx .EQ. (input1%nx/2)+1 &
     &      .AND. output%ny .EQ. input1%ny       &
     &      .AND. output%nz .EQ. input1%nz       &
     &      .AND. output%nb_components .EQ. input1%nb_components       &
     &     ) &
     &) THEN
     WRITE(6,'(a,i0,a,8(i0,a))')&
   &  '[ERROR] ftran_wrap on ',me,': incompatible data sizes!'//    &
   &  ' Real part is [',input1%nx,',',input1%ny,',',input1%nz,      &
   &  ',',input1%nb_components ,                                    &
   &  '] and Complex is [',output%nx,',',output%ny,',',output%nz,   &
   &  ',',output%nb_components ,']'
     RETURN
  ENDIF

  !Check  we are along X for Input and Output for the REAL -> complex FFT
  IF (.NOT. setAlongX(me,input1)) RETURN
  IF (.NOT. setAlongX(me,output)) RETURN

  ! clear the output
  output%values=0.0_WP

  ! Sizes of the local real array
  iarray(1)=input1%nx
  iarray(2)=input1%ymax-input1%ymin+1
  iarray(3)=input1%zmax-input1%zmin+1
  iarray(4)=input1%nb_components

  ! Sizes of the local complex array
  oarray(1)=output%nx
  oarray(2)=output%ymax-input1%ymin+1
  oarray(3)=output%zmax-input1%zmin+1
  oarray(4)=output%nb_components

  ! number of ffts along X to execute
  nbfft1d=(iarray(2)*iarray(3)*iarray(4))

  ! build the plan for FFTs along X and run FFT real -> complex
  CALL dfftw_plan_many_dft_r2c(plan,1,input1%nx,nbfft1d,input1%values,iarray,1,input1%nx,&
                             & output%values,oarray,1,output%nx,FFTW_ESTIMATE)
  CALL dfftw_execute_dft_r2c(plan,input1%values,output%values)
  CALL dfftw_destroy_plan(plan)

  ! Now proceed along Y
  IF (.NOT. setAlongy(me,output)) RETURN

  ! first use a temp array to reorganize datas with contiguous position along Y
  ALLOCATE (buff(output%ymin:output%ymax,output%xmin:output%xmax,output%zmin:output%zmax,input1%nb_components),stat=ires)
  IF (ires .NE.0) THEN
     WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": not enought memory for buff in ftran_wrap!"
     RETURN
  ENDIF

  !put datas in buff
  DO comp = 1, output%nb_components
    DO zz=output%zmin, output%zmax
      DO yy=output%ymin, output%ymax
        DO xx=output%xmin, output%xmax
          buff(yy,xx,zz,comp)=output%values(xx,yy,zz,comp)
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  !proceed to the second FFT complex-complex "in place"
  ! Sizes of the local complex array
  iarray(1)=output%ny
  iarray(2)=output%xmax-output%xmin+1
  iarray(3)=output%zmax-output%zmin+1
  oarray(4)=output%nb_components

  ! number of ffts along Y to execute
  nbfft1d=(iarray(2)*iarray(3)*iarray(4))

  ! build the plan for FFTs along X and run FFT complex -> complex
  ! fft along Y
  CALL dfftw_plan_many_dft(plan,1,iarray(1),nbfft1d,buff,iarray,1,iarray(1),&
                         & buff,iarray,1,iarray(1),FFTW_FORWARD,FFTW_ESTIMATE)
  CALL dfftw_execute_dft(plan,buff,buff)
  CALL dfftw_destroy_plan(plan)

  !reorder datas in  output%values. Swap loops for cache optimization at read time.
  DO comp = 1, output%nb_components
    DO zz=output%zmin, output%zmax
      DO xx=output%xmin, output%xmax
        DO yy=output%ymin, output%ymax
          output%values(xx,yy,zz,comp)=buff(yy,xx,zz,comp)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  DEALLOCATE (buff)

  ! Now proceed along Z
  IF (.NOT. setAlongz(me,output)) RETURN

  !first use a temp array to reorganize datas with contiguous position along Z
  ALLOCATE (buff(output%zmin:output%zmax,output%xmin:output%xmax,output%ymin:output%ymax,input1%nb_components),stat=ires)
  IF (ires .NE.0) THEN
     WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": not enought memory for buff in ftran_wrap!"
     RETURN
  ENDIF

  !put datas in buff
  DO comp = 1, output%nb_components
    DO zz=output%zmin, output%zmax
      DO yy=output%ymin, output%ymax
        DO xx=output%xmin, output%xmax
           buff(zz,xx,yy,comp)=output%values(xx,yy,zz,comp)
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  !proceed to the third FFT complex-complex "in place"
  ! Sizes of the local complex array
  iarray(1)=output%nz
  iarray(2)=output%xmax-output%xmin+1
  iarray(3)=output%ymax-output%ymin+1
  iarray(4)=output%nb_components

  ! number of ffts along X to execute
  nbfft1d=(iarray(2)*iarray(3)*iarray(4))

  ! build the plan for FFTs along X and run FFT real -> complex
  ! fft along Z
  CALL dfftw_plan_many_dft(plan,1,iarray(1),nbfft1d,buff,iarray,1,iarray(1),&
                         & buff,iarray,1,iarray(1),FFTW_FORWARD,FFTW_ESTIMATE)
  CALL dfftw_execute_dft(plan,buff,buff)
  CALL dfftw_destroy_plan(plan)

  !reorder datas in  output%values. Swap loops for cache optimization at read time.
  DO comp = 1, output%nb_components
    DO yy=output%ymin, output%ymax
      DO xx=output%xmin, output%xmax
        DO zz=output%zmin, output%zmax
           output%values(xx,yy,zz,comp)=buff(zz,xx,yy,comp)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  DEALLOCATE (buff)

  !normalize now by the number of points
  !should prevent any overflow at hight resolution.
  output%values=output%values/REAL(input1%nx,WP)
  output%values=output%values/REAL(input1%ny,WP)
  output%values=output%values/REAL(input1%nz,WP)

  !CALL deletedataLayout(input)

  res=.TRUE.
  RETURN

  END SUBROUTINE ftran_vector


  !------------------------------------------------------------------------------
  !> Backward Transform
  !> @author
  !> Patrick Begou and small modification by Jean-Baptiste Lagaert
  !
  !
  !> @details
  !> The aim of this subroutine is to run a 3D FFT on the domain. The input domain should
  !> be provided "AlongZ" (it is set AlongZ if it is not) and the results in output
  !> is organized along X. This version is very basic without optimizations. For
  !> complex/complex it uses the "in place" algorithm (same array for input and output)
  !> Plans are evaluated for each call. The only optimization is to provide data contiguous
  !> in the FFT direction in the 3 steps (3 time a 1D FFT)
  !> WARNING: data are normalized by the number of points when the backward FFT returns
  !> in the physical space
  !> @param[in] input1 the complex data in the spectral space.
  !> @param[in,out] output the real data in the physical space
  !> @param[out] res set to .TRUE. if all is successfull.
  !------------------------------------------------------------------------------
  SUBROUTINE btran_vector(input1,output,res)

  USE parallel_tools
  USE datalayout

  IMPLICIT NONE
  TYPE(COMPLEX_VECTOR_DATA_LAYOUT), INTENT(INOUT)   ::input1
  TYPE(REAL_VECTOR_DATA_LAYOUT), INTENT(INOUT)      ::output
  LOGICAL, INTENT(OUT) :: res

  !TYPE(COMPLEX_VECTOR_DATA_LAYOUT)    ::input
  INTEGER(KIND=8) :: plan
  INTEGER :: me,nbfft1d
  INTEGER, DIMENSION(4):: iarray, oarray
  COMPLEX(WP), DIMENSION(:,:,:,:), ALLOCATABLE::buff
  INTEGER  :: ires,xx,yy,zz,comp

  res=.FALSE.
  me=getnbmycpu()

  !input=input1


  ! Check the output data structure.
  IF (.NOT.(     (output%nx/2)+1 .EQ. input1%nx  &
     &      .AND. output%ny .EQ. input1%ny       &
     &      .AND. output%nz .EQ. input1%nz       &
     &      .AND. output%nb_components .EQ. input1%nb_components       &
     &     ) &
     &) THEN
     WRITE(6,'(a,i0,a,8(i0,a))')&
   &  '[ERROR] Btran_wrap on ',me,': incompatible data sizes!'//    &
   &  ' Complex part is [',input1%nx,',',input1%ny,',',input1%nz,   &
   &  ',',input1%nb_components ,                                    &
   &  '] and Real is [',output%nx,',',output%ny,',',output%nz,      &
   &  ',',output%nb_components ,']'
     RETURN
  ENDIF

  !Check  we are along Z for Input and X for output (the last FFT Complex->Real is along X)
  IF (.NOT. setAlongZ(me,input1)) RETURN
  IF (.NOT. setAlongX(me,output)) RETURN

  ! clear the output
  output%values=0.0_WP

  ! First proceed along Z
  ! Use a temp array to reorganize datas with contiguous position along Z
  ALLOCATE (buff(input1%zmin:input1%zmax,input1%xmin:input1%xmax,input1%ymin:input1%ymax,input1%nb_components),stat=ires)
  IF (ires .NE.0) THEN
     WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": not enought memory for buff in Btran_wrap!"
     RETURN
  ENDIF

  !put datas in buff
  DO comp = 1, output%nb_components
    DO zz=input1%zmin, input1%zmax
      DO yy=input1%ymin, input1%ymax
        DO xx=input1%xmin, input1%xmax
           buff(zz,xx,yy,comp)=input1%values(xx,yy,zz,comp)
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  !proceed to the first backward FFT complex-complex "in place"
  ! Sizes of the local complex array
  iarray(1)=input1%nz
  iarray(2)=input1%xmax-input1%xmin+1
  iarray(3)=input1%ymax-input1%ymin+1
  iarray(4)=input1%nb_components

  ! number of ffts along X to execute
  nbfft1d=(iarray(2)*iarray(3)*iarray(4))

  ! build the plan for FFTs along X and run FFT real -> complex
  CALL dfftw_plan_many_dft(plan,1,iarray(1),nbfft1d,buff,iarray,1,iarray(1),&
                         & buff,iarray,1,iarray(1),FFTW_BACKWARD,FFTW_ESTIMATE)
  CALL dfftw_execute_dft(plan,buff,buff)
  CALL dfftw_destroy_plan(plan)

  !reorder datas in  output%values. Swap loops for cache optimization at read time.
  DO comp = 1, input1%nb_components
    DO yy=input1%ymin, input1%ymax
      DO xx=input1%xmin, input1%xmax
        DO zz=input1%zmin, input1%zmax
           input1%values(xx,yy,zz,comp)=buff(zz,xx,yy,comp)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  DEALLOCATE (buff)

  ! Now proceed along Y
  IF (.NOT. setAlongy(me,input1)) RETURN

  !first use a temp array to reorganize datas with contiguous position along Y
  ALLOCATE (buff(input1%ymin:input1%ymax,input1%xmin:input1%xmax,input1%zmin:input1%zmax,input1%nb_components),stat=ires)
  IF (ires .NE.0) THEN
     WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": not enought memory for buff in ftran_wrap!"
     RETURN
  ENDIF

  !put datas in buff
  DO comp = 1, input1%nb_components
    DO zz=input1%zmin, input1%zmax
      DO yy=input1%ymin, input1%ymax
        DO xx=input1%xmin, input1%xmax
           buff(yy,xx,zz,comp)=input1%values(xx,yy,zz,comp)
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  !proceed to the second FFT complex-complex "in place"
  ! Sizes of the local complex array
  iarray(1)=input1%ny
  iarray(2)=input1%xmax-input1%xmin+1
  iarray(3)=input1%zmax-input1%zmin+1
  iarray(4)=input1%nb_components

  ! number of ffts along X to execute
  nbfft1d=(iarray(2)*iarray(3)*iarray(4))

  ! build the plan for FFTs along X and run FFT complex -> complex
  CALL dfftw_plan_many_dft(plan,1,iarray(1),nbfft1d,buff,iarray,1,iarray(1),&
                         & buff,iarray,1,iarray(1),FFTW_BACKWARD,FFTW_ESTIMATE)
  CALL dfftw_execute_dft(plan,buff,buff)
  CALL dfftw_destroy_plan(plan)

  !reorder datas in  output%values. Swap loops for cache optimization at read time.
  DO comp = 1, input1%nb_components
    DO zz=input1%zmin, input1%zmax
      DO xx=input1%xmin, input1%xmax
        DO yy=input1%ymin, input1%ymax
           input1%values(xx,yy,zz,comp)=buff(yy,xx,zz,comp)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  DEALLOCATE (buff)

  ! Run the third FFT Complex->Real along X
  IF (.NOT. setAlongX(me,input1)) RETURN

  ! Sizes of the local real array
  iarray(1)=input1%nx
  iarray(2)=input1%ymax-input1%ymin+1
  iarray(3)=input1%zmax-input1%zmin+1
  iarray(4)=input1%nb_components

  ! Sizes of the local complex array
  oarray(1)=output%nx
  oarray(2)=output%ymax-output%ymin+1
  oarray(3)=output%zmax-output%zmin+1
  iarray(4)=output%nb_components

  ! number of ffts along X to execute
  nbfft1d=(iarray(2)*iarray(3)*iarray(4))

  ! build the plan for FFTs along X and run FFT real -> complex
  CALL dfftw_plan_many_dft_c2r(plan,1,output%nx,nbfft1d,input1%values,iarray,1,input1%nx,&
                             & output%values,oarray,1,output%nx,FFTW_ESTIMATE)
  CALL dfftw_execute_dft_c2r(plan,input1%values,output%values)
  CALL dfftw_destroy_plan(plan)

  !CALL deletedataLayout(input)

  res=.TRUE.
  RETURN

  END SUBROUTINE btran_vector


  !------------------------------------------------------------------------------
  !> Compute backward transform along Z only. The purpose is to obtain value in
  !! spectral space along X and Y only.
  !! @author
  !! Jean-Baptiste Lagaert and Patrick Begou, LEGI
  !
  !! @details
  !! It is useful, for instance, for some plane-jet post-process like 2D
  !! spectrum along X and Y in function of Z.
  !> The input domain should be provided "AlongZ" (it is set AlongZ if it is not)
  !! and the results in output is still organized along Z as it use suppose to
  !! use "like that" wich allows to avoid communications.
  !> This version is very basic without optimizations. For complex/complex it uses the
  !! "in place" algorithm (same array for input and output). Plans are evaluated for 
  !! each call. The only optimization is to provide data contiguous in the FFT 
  !! direction in the 3 steps (3 time a 1D FFT)
  !! WARNING: data are normalized by the number of points when the backward FFT returns
  !! in the physical space
  !! @param[in] input1 the complex data in the spectral space.
  !! @param[in,out] output the complex data in the physical space along Z and
  !!                    spectral spaces along X and Y.
  !! @param[out] res set to .TRUE. if all is successfull.
  !------------------------------------------------------------------------------
  SUBROUTINE btran_alongZ(input1,output,res)

  USE parallel_tools
  USE datalayout

  IMPLICIT NONE
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) ::input1
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT) ::output
  LOGICAL, INTENT(OUT) :: res

  INTEGER(KIND=8) :: plan
  INTEGER :: me,nbfft1d
  INTEGER, DIMENSION(3):: iarray
  COMPLEX(WP), DIMENSION(:,:,:), ALLOCATABLE::buff
  INTEGER  :: ires,xx,yy,zz

  res=.FALSE.
  me=getnbmycpu()

  output=input1


  ! Check the output data structure.
  IF (.NOT. samelayout(input1, output)) then
     WRITE(6,'(a,i0,a,6(i0,a))')&
   &  '[ERROR] Btran_wrap on ',me,': incompatible data sizes!'//   &
   &  ' Complex input part is [',input1%nx,',',input1%ny,',',input1%nz,        &
   &  '] and output is [',output%nx,',',output%ny,',',output%nz,']' 
     RETURN
  ENDIF

  !Check  we are along Z for Input and X for output (the last FFT Complex->Real is along X)
  IF (.NOT. setAlongZ(me,output)) RETURN

  ! clear the output
  output%values=0.0_WP

  ! First proceed along Z
  ! Use a temp array to reorganize datas with contiguous position along Z
  ALLOCATE (buff(output%zmin:output%zmax,output%xmin:output%xmax,output%ymin:output%ymax),stat=ires)
  IF (ires .NE.0) THEN
     WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": not enought memory for buff in Btran_wrap!"
     RETURN
  ENDIF

  !put datas in buff
  DO zz=output%zmin, output%zmax
    DO yy=output%ymin, output%ymax
      DO xx=output%xmin, output%xmax
         buff(zz,xx,yy)=input1%values(xx,yy,zz)
      ENDDO
    ENDDO
  ENDDO

  !proceed to the first backward FFT complex-complex "in place"
  ! Sizes of the local complex array
  iarray(1)=input1%nz
  iarray(2)=input1%xmax-input1%xmin+1
  iarray(3)=input1%ymax-input1%ymin+1

  ! number of ffts along X to execute
  nbfft1d=(iarray(2)*iarray(3))

  ! build the plan for FFTs along X and run FFT real -> complex
  CALL dfftw_plan_many_dft(plan,1,iarray(1),nbfft1d,buff,iarray,1,iarray(1),&
                         & buff,iarray,1,iarray(1),FFTW_BACKWARD,FFTW_ESTIMATE)
  CALL dfftw_execute_dft(plan,buff,buff)
  CALL dfftw_destroy_plan(plan)

  !reorder datas in  output%values. Swap loops for cache optimization at read time.
  DO yy=output%ymin, output%ymax
    DO xx=output%xmin, output%xmax
      DO zz=output%zmin, output%zmax
         output%values(xx,yy,zz)=buff(zz,xx,yy)
      ENDDO
    ENDDO
  ENDDO
  DEALLOCATE (buff)

  res=.TRUE.
  RETURN

  END SUBROUTINE btran_alongZ
  !------------------------------------------------------------------------------

  !! @author
  !! Luca Marradi, LIPhy 2014
  !
  !! @details
  !! Convolution product is computed 
  !> The input domain should be provided "AlongZ" (it is set AlongZ if it is
  !not)
  !returns
  !! @param[in]    G the propagator in the spectral space
  !! @param[in]    sigmak_ij in the spectral space
  !! @param[out]   output store the convolution in the spectral space
  !! @param[out]   res set to .TRUE. if all is successfull.
  !------------------------------------------------------------------------------

 SUBROUTINE convolution(convol,sigmak,G,res)
 
    USE parallel_tools
    USE datalayout
   
    TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT) :: G,sigmak
    TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT) :: convol
  
    INTEGER ::xx,yy,zz
    LOGICAL, INTENT(OUT) :: res

    res=.FALSE.

    DO zz=convol%zmin, convol%zmax
     DO yy=convol%ymin, convol%ymax
      DO xx=convol%xmin, convol%xmax
         convol%values(xx,yy,zz)=G%values(xx,yy,zz) * sigmak%values(xx,yy,zz)
      ENDDO
     ENDDO
    ENDDO

   res=.TRUE.
   RETURN 

  END SUBROUTINE convolution

END MODULE transforms_tools
!> @}
