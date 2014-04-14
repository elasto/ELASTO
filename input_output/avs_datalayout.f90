!USEFORTEST io
!> @addtogroup output
!! @{

!------------------------------------------------------------------------------
!
! MODULE:  avs
!
!> @author
!> Patrick Begou, LEGI
!
! DESCRIPTION: 
!> The aim of this module is to provide utilties for creating files for AVS
!> vizexpress software.
!------------------------------------------------------------------------------
MODULE avs

    use precision_tools
    use avs_real
    use avs_cmplx

    implicit none

    private buildFldFile

    interface ReadFromFile
        module procedure real_ReadFromFile
        module procedure cmplx_ReadFromFile 
    end interface ReadFromFile

    interface WriteToFile
       module procedure real_WriteToFile
       module procedure cmplx_WriteToFile
    end interface WriteToFile

contains

!------------------------------------------------------------------------------
!> @author 
!> Patrick Begou
!
!> @details
!> This subroutine output files tailored for use with AVS Vizexpress software.
!> It generates binary data files, asscii coordinate files and a field file
!> describing the datas. It can be called for one to four REAL_DATA_LAYOUT items
!> but all must have the same layout (same global and local grid size).
!
!> @param[in] name the basename of the file to create. We'll get name.fld, name.coordx,
!> name.coordy...... name_//U%name.data.... etc
!> @param[in] me the rank number of the process.
!> @param[in] lx the length of the global domain along X.
!> @param[in] ly the length of the global domain along Y.
!> @param[in] lz the length of the global domain along Z.
!> @param[in] U is a variable of type REAL_DATA_LAYOUT to be dumped in avs files.
!> it could be the U velocity but every thing else!
!> @param[in] V is an optional variable of type REAL_DATA_LAYOUT to be dumped in avs files. It
!> should have the same layout than U.
!> @param[in] W is an optional variable of type REAL_DATA_LAYOUT to be dumped in avs files. It
!> should have the same layout than U.
!> @param[in] P is an optional variable of type REAL_DATA_LAYOUT to be dumped in avs files. It
!> should have the same layout than U.
!> @return true if files are successfuly created
!------------------------------------------------------------------------------
FUNCTION dump_UVWP_ForAVS(name,me,lx,ly,lz,U,V,W,P,usename)
  USE datalayout
  USE parallel_tools

  CHARACTER(LEN=*), INTENT(IN)                :: name
  INTEGER, INTENT(IN)                         :: me
  TYPE(REAL_DATA_LAYOUT),INTENT(IN)           :: U
  TYPE(REAL_DATA_LAYOUT),INTENT(IN), OPTIONAL :: V
  TYPE(REAL_DATA_LAYOUT),INTENT(IN), OPTIONAL :: W
  TYPE(REAL_DATA_LAYOUT),INTENT(IN), OPTIONAL :: P
  LOGICAL,INTENT(IN), OPTIONAL :: usename
  REAL(WP), INTENT(IN)::lx,ly,lz
  LOGICAL :: dump_UVWP_ForAVS


  CHARACTER(len=str_medium) ::labeltab(4)
  INTEGER:: nchamp
  LOGICAL:: withname

  dump_UVWP_ForAVS=.FALSE.
  !IF(me.EQ.0) WRITE(6,'(a)')'[INFO] Writing AVS data files'
  if (present(usename))then
    if (usename) then
      withname = .true.
    else
      withname = .false.
    endif
  else
    withname = .false.
  endif

  ! How many data values to dump? Are they compatibles ?
  IF (PRESENT(P)) THEN
    nchamp=4
    IF (.NOT. samelayout(U,V,W,P)) THEN
         WRITE(6,'(a,i0,a)') &
    &   "ERROR dump_UVWP_ForAVS on process ",me, &
    &   ": the 4 data fields have not the same layout!"
         RETURN
     ENDIF
     labeltab(1)=TRIM(ADJUSTL(U%name))
     labeltab(2)=TRIM(ADJUSTL(V%name))
     labeltab(3)=TRIM(ADJUSTL(W%name))
     labeltab(4)=TRIM(ADJUSTL(P%name))
     if (withname) then
       IF(.NOT. real_WriteToFile(me,TRIM(labeltab(1))//'.data',U)) RETURN
       IF(.NOT. real_WriteToFile(me,TRIM(labeltab(2))//'.data',V)) RETURN
       IF(.NOT. real_WriteToFile(me,TRIM(labeltab(3))//'.data',W)) RETURN
       IF(.NOT. real_WriteToFile(me,TRIM(labeltab(4))//'.data',P)) RETURN
     else
       IF(.NOT. real_WriteToFile(me,TRIM(name)//'_'//TRIM(labeltab(1))//'.data',U)) RETURN
       IF(.NOT. real_WriteToFile(me,TRIM(name)//'_'//TRIM(labeltab(2))//'.data',V)) RETURN
       IF(.NOT. real_WriteToFile(me,TRIM(name)//'_'//TRIM(labeltab(3))//'.data',W)) RETURN
       IF(.NOT. real_WriteToFile(me,TRIM(name)//'_'//TRIM(labeltab(4))//'.data',P)) RETURN
     endif
  ELSE IF(PRESENT(W)) THEN
     nchamp=3
     IF (.NOT. samelayout(U,V,W)) THEN
         WRITE(6,'(a,i0,a)') &
        &   "ERROR dump_UVWP_ForAVS on process ",me, &
        &   ": the 3 data fields have not the same layout!"
         RETURN
     ENDIF
     labeltab(1)=TRIM(ADJUSTL(U%name))
     labeltab(2)=TRIM(ADJUSTL(V%name))
     labeltab(3)=TRIM(ADJUSTL(W%name))
     if (withname) then
       IF(.NOT. real_WriteToFile(me,TRIM(labeltab(1))//'.data',U)) RETURN
       IF(.NOT. real_WriteToFile(me,TRIM(labeltab(2))//'.data',V)) RETURN
       IF(.NOT. real_WriteToFile(me,TRIM(labeltab(3))//'.data',W)) RETURN
     else
       IF(.NOT. real_WriteToFile(me,TRIM(name)//'_'//TRIM(labeltab(1))//'.data',U)) RETURN
       IF(.NOT. real_WriteToFile(me,TRIM(name)//'_'//TRIM(labeltab(2))//'.data',V)) RETURN
       IF(.NOT. real_WriteToFile(me,TRIM(name)//'_'//TRIM(labeltab(3))//'.data',W)) RETURN
     endif
  ELSE IF(PRESENT(V)) THEN
     nchamp=2
     IF (.NOT. samelayout(U,V)) THEN
         WRITE(6,'(a,i0,a)') &
     &   "ERROR dump_UVWP_ForAVS on process ",me, &
     &   ": the 2 data fields have not the same layout!"
         RETURN
     ENDIF
     labeltab(1)=TRIM(ADJUSTL(U%name))
     labeltab(2)=TRIM(ADJUSTL(V%name))
     if (withname) then
       IF(.NOT. real_WriteToFile(me,TRIM(labeltab(1))//'.data',U)) RETURN
       IF(.NOT. real_WriteToFile(me,TRIM(labeltab(2))//'.data',V)) RETURN
     else
       IF(.NOT. real_WriteToFile(me,TRIM(name)//'_'//TRIM(labeltab(1))//'.data',U)) RETURN
       IF(.NOT. real_WriteToFile(me,TRIM(name)//'_'//TRIM(labeltab(2))//'.data',V)) RETURN
     endif
  ELSE
     nchamp=1
     labeltab(1)=TRIM(ADJUSTL(U%name))
     if (withname) then
       IF(.NOT. real_WriteToFile(me,TRIM(labeltab(1))//'.data',U)) RETURN
     else
       IF(.NOT. real_WriteToFile(me,TRIM(name)//'.data',U)) RETURN
     endif
  ENDIF

  ! creates ASCII coordinates files
  IF (me.EQ.0) THEN  
    IF (.NOT. buildGridFile(name,U%nx,U%ny,U%nz,lx,ly,lz)) RETURN
    IF (withname) THEN
      IF (.NOT. buildFldFile(name,labeltab,nchamp,U%nx,U%ny,U%nz,withname)) RETURN
    ELSE
      IF (.NOT. buildFldFile(name,labeltab,nchamp,U%nx,U%ny,U%nz)) RETURN
    ENDIF
  ENDIF

  dump_UVWP_ForAVS=.TRUE.
  RETURN

END FUNCTION dump_UVWP_ForAVS

!------------------------------------------------------------------------------
!> @author 
!> Luca MARRADI LIPhy 2014
!
!> @details
!> This subroutine output files tailored for use with AVS Vizexpress software.
!> It generates binary data files, asscii coordinate files and a field file
!> describing the datas. It can be called for one to four REAL_DATA_LAYOUT items
!> but all must have the same layout (same global and local grid size).
!
!> @param[in] name the basename of the file to create. We'll get name.fld, name.coordx,
!> name.coordy...... name_//U%name.data.... etc
!> @param[in] me the rank number of the process.
!> @param[in] lx the length of the global domain along X.
!> @param[in] ly the length of the global domain along Y.
!> @param[in] lz the length of the global domain along Z.
!> @param[in] sigma is a variable of type REAL_DATA_LAYOUT to be dumped in avs files.
!> @return true if files are successfuly created
!------------------------------------------------------------------------------
FUNCTION dump_UVWP_ForAVS_luca(name,me,lx,ly,lz,sigma,usename)
  USE datalayout
  USE parallel_tools

  IMPLICIT NONE

  !INPUT/OUTPUT DATA
  CHARACTER(LEN=*), INTENT(IN)                :: name
  INTEGER, INTENT(IN)                         :: me
  TYPE(REAL_DATA_LAYOUT),INTENT(IN)           :: sigma
  LOGICAL,INTENT(IN), OPTIONAL :: usename
  REAL(WP), INTENT(IN)::lx,ly,lz
  LOGICAL :: dump_UVWP_ForAVS_luca

  !== LOCAL DATA ==
  CHARACTER(len=str_medium) ::labeltab(4)
  INTEGER:: nchamp
  LOGICAL:: withname

  dump_UVWP_ForAVS_luca=.FALSE.
  !IF(me.EQ.0) WRITE(6,'(a)')'[INFO] Writing AVS data files'
  IF (present(usename))THEN
    if (usename) then
      withname = .true.
    else
      withname = .false.
    endif
  else
    withname = .false.
  ENDIF

     nchamp=1
     labeltab(1)=TRIM(ADJUSTL(sigma%name))
     IF (withname) THEN
       IF(.NOT. real_WriteToFile(me,TRIM(labeltab(1))//'.data',sigma)) RETURN
     else
       IF(.NOT. real_WriteToFile(me,TRIM(name)//'.data',sigma)) RETURN
     ENDIF

  !== CREATES ASCII COORDINATES FILE ==
  IF (me.EQ.0) THEN  
    IF (.NOT. buildGridFile(name,sigma%nx,sigma%ny,sigma%nz,lx,ly,lz)) RETURN
    IF (withname) THEN
      IF (.NOT. buildFldFile(name,labeltab,nchamp,sigma%nx,sigma%ny,sigma%nz,withname)) RETURN
    ELSE
      IF (.NOT. buildFldFile(name,labeltab,nchamp,sigma%nx,sigma%ny,sigma%nz)) RETURN
    ENDIF
  ENDIF

  dump_UVWP_ForAVS_luca=.TRUE.
  RETURN

END FUNCTION dump_UVWP_ForAVS_luca

!------------------------------------------------------------------------------
!> @author 
!> Patrick Begou
!
!> @details
!> Create *.fld file
!
!> @return true if files are successfuly created
!------------------------------------------------------------------------------
FUNCTION buildFldFile(name,labeltab,nchamp,nx,ny,nz,usename)

  USE fileio

  CHARACTER(LEN=*), INTENT(IN) :: name
  INTEGER, INTENT(IN)          :: nchamp
  CHARACTER(LEN=*), INTENT(IN) :: labeltab(nchamp)
  INTEGER, INTENT(IN)          :: nx,ny,nz
  LOGICAL, INTENT(IN), OPTIONAL:: usename
  LOGICAL ::buildFldFile

  CHARACTER(LEN=str_long)::filename,buff
  INTEGER :: fich,i
  LOGICAL :: withname
  
  buildFldFile=.FALSE.
  if (present(usename))then
    if (usename) then
      withname = .true.
    else
      withname = .false.
    endif
  else
    withname = .false.
  endif

  fich=iopen()
  filename=TRIM(ADJUSTL(name))//'.fld'
  OPEN(UNIT=fich, FILE=filename, FORM="FORMATTED", ERR=100)
  REWIND(fich)
  
  WRITE(fich,'(a)',err=101)'# AVS field'
  WRITE(fich,'(a)',err=101)'# Generated by Codescalar'
  WRITE(fich,'(2x,a)',err=101)   'ndim   =  3'           ! Dim dataset
  WRITE(fich,'(2x,a)',err=101)   'nspace =  3'           ! Dim computational space
  WRITE(fich,'(2x,a,i0)',err=101)   'veclen =  ',nchamp  ! number of scalar fields
  WRITE(fich,'(2x,a,i3)',err=101)'dim1   = ',nx          ! first dimension
  WRITE(fich,'(2x,a,i3)',err=101)'dim2   = ',ny          ! second dimension
  WRITE(fich,'(2x,a,i3)',err=101)'dim3   = ',nz          ! third dimension
  WRITE(fich,'(2x,a)',err=101)   'data   = double'       ! type of data
  WRITE(fich,'(2x,a)',err=101)   'field  = rectilinear'  ! type of field
  WRITE(fich,'(a)',err=101)'#'
  
  ! Labels
  buff='label ='
  DO i=1,nchamp
      buff=TRIM(buff)//" "//TRIM(labeltab(i))
  ENDDO
  WRITE(fich,'(a)',err=101) TRIM(buff)

  ! Data files
  IF (withname) THEN
    DO i=1,nchamp
      WRITE(fich,'(a,i0,a)',err=101) 'variable ',i, &
      &' file = '//TRIM(labeltab(i))//'.data type=binary skip=12'
    ENDDO
  ELSE
    if ( nchamp .eq. 1 ) then
      WRITE(fich,'(a)',err=101) 'variable 1 &
      & file = '//TRIM(name)//'.data type=binary skip=12'
    else 
      DO i=1,nchamp
        WRITE(fich,'(a,i0,a)',err=101) 'variable ',i, &
        &' file = '//TRIM(name)//'_'//TRIM(labeltab(i))//'.data type=binary skip=12'
      ENDDO
    ENDIF
  ENDIF

  ! Coordinates  
  WRITE(fich,'(a)',err=101)'coord 1 file =  '//TRIM(name)//'.coordX filetype = ascii  skip = 0'
  WRITE(fich,'(a)',err=101)'coord 2 file =  '//TRIM(name)//'.coordY filetype = ascii  skip = 0'
  WRITE(fich,'(a)',err=101)'coord 3 file =  '//TRIM(name)//'.coordZ filetype = ascii  skip = 0'
  
  CLOSE(fich)
    
  fich=iclose(fich)
  
  buildFldFile=.TRUE.
  RETURN
100 WRITE(6,'(a)') 'Unable to open file "'//TRIM(filename)//'"'
  RETURN
101 WRITE(6,'(a)') 'Write error on file "'//TRIM(filename)//'"'
  RETURN

END FUNCTION buildFldFile


!------------------------------------------------------------------------------
!> @author 
!> Patrick Begou
!
!> @details
!> This subroutine output ascii coordinates files tailored for use with AVS Vizexpress software.
!
!> @param[in] name the basename of the file to create. We'll create
!> name.coordX, name.coordY and name.coordZ ascii files.
!> @param[in] nx the size of the global domain along X.
!> @param[in] ny the size of the global domain along Y.
!> @param[in] nz the size of the global domain along Z.
!> @param[in] lx the length of the global domain along X.
!> @param[in] ly the length of the global domain along Y.
!> @param[in] lz the length of the global domain along Z.
!> @return true if files are successfuly created
!------------------------------------------------------------------------------
FUNCTION buildGridFile(name,nx,ny,nz,lx,ly,lz)

  USE fileio

  CHARACTER(LEN=*), INTENT(IN) :: name
  INTEGER, INTENT(IN)          :: nx,ny,nz
  REAL(WP), INTENT(IN)         :: lx,ly,lz
  LOGICAL ::buildGridFile
  
  CHARACTER(LEN=str_long)::filename
  INTEGER :: fich,i
  
  buildGridFile=.FALSE.
  fich=iopen()

  !X coordinates
  filename=TRIM(ADJUSTL(name))//'.coordX'
  OPEN(UNIT=fich, FILE=filename, FORM="FORMATTED", ERR=100)
  REWIND(fich)
  DO i=0, nx-1
     WRITE(fich,'(f12.6)',err=110) i*lx/(nx)
  ENDDO  
  CLOSE(fich)

  !Y coordinates
  filename=TRIM(ADJUSTL(name))//'.coordY'
  OPEN(UNIT=fich, FILE=filename, FORM="FORMATTED", ERR=100)
  REWIND(fich)
  DO i=0, ny-1
     WRITE(fich,'(f12.6)',err=110) i*ly/(ny)
  ENDDO  
  CLOSE(fich)

  !Z coordinates
  filename=TRIM(ADJUSTL(name))//'.coordZ'
  OPEN(UNIT=fich, FILE=filename, FORM="FORMATTED", ERR=100)
  REWIND(fich)
  DO i=0, nz-1
     WRITE(fich,'(f12.6)',err=110) i*lz/(nz)
  ENDDO  
  CLOSE(fich)
    
  fich=iclose(fich)
  buildGridFile=.TRUE.
  RETURN

100 WRITE(6,'(a)') 'Unable to open file "'//TRIM(filename)//'"'
  RETURN
110 WRITE(6,'(a)') 'Write error on file "'//TRIM(filename)//'"'
  RETURN
  

END FUNCTION buildGridFile

END MODULE avs
!> @}
