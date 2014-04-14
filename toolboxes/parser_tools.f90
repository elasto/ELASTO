!USEFORTEST toolbox
!USEFORTEST avgcond 
!USEFORTEST postprocess
!USEFORTEST advec
!USEFORTEST io
!USEFORTEST topo
!------------------------------------------------------------------------------
!
! MODULE: Parser
!
!> @author
!> Guillaume Balarac, LEGI
!
! DESCRIPTION:
!> The aim of this module is to store in a large array the parameters of
!! the simulation provided in an imput file.
!! It provides access functions to
!! the informations for the code with data conversion (all is stored as strings)
!! Provides:
!!   SUBROUTINE parser_parsefile()
!!   SUBROUTINE parser_show()
!!   SUBROUTINE parser_init()
!!   SUBROUTINE parser_read()
!
!------------------------------------------------------------------------------
MODULE parser_tools

  USE precision_tools
  IMPLICIT NONE

  ! Input values storage
  INTEGER,PRIVATE, PARAMETER :: tag_length=128
  INTEGER,PRIVATE, PARAMETER :: line_length=4960
  INTEGER, PRIVATE           :: nfields=0
  TYPE entry_type
     CHARACTER(tag_length)     :: tag
     CHARACTER(line_length)    :: value
     TYPE(entry_type), POINTER :: next
  END TYPE entry_type

  TYPE(entry_type), POINTER, PRIVATE :: entries=>NULL()
  TYPE(entry_type), POINTER, PRIVATE :: current=>NULL()

   !---------------------------------------------------------------------------
   !> @author
   !> Guillaume Balarac, LEGI.
   !
   ! DESCRIPTION:
   !> Public interface for datas access in the entries array built from input file.
   !
   !> @return none
   !---------------------------------------------------------------------------
  INTERFACE parser_read
     MODULE PROCEDURE parser_readlogical
     MODULE PROCEDURE parser_readint
     MODULE PROCEDURE parser_readintarray
     MODULE PROCEDURE parser_readfloat
     MODULE PROCEDURE parser_readfloatarray
     MODULE PROCEDURE parser_readfloatarray2D
     MODULE PROCEDURE parser_readchar
     MODULE PROCEDURE parser_readchararray
  END INTERFACE

    PUBLIC     parser_parsefile
    PUBLIC     parser_newentry
    PUBLIC     parser_is_defined
    PUBLIC     parser_read
    PUBLIC     parser_reset
    PUBLIC     parser_show
    PUBLIC     parser_getsize
    PUBLIC     parseFilter

    PRIVATE    parser_readlogical
    PRIVATE    parser_readint
    PRIVATE    parser_readintarray
    PRIVATE    parser_readfloat
    PRIVATE    parser_readfloatarray
    PRIVATE    parser_readfloatarray2D
    PRIVATE    parser_readchar
    PRIVATE    parser_readchararray
    PRIVATE    parser_fieldfortag
    PRIVATE    parser_free_entry
    PRIVATE    parser_showAlltag

CONTAINS

   !---------------------------------------------------------------------------
   !> @author
   !> Patrick BEGOU, LEGI.
   !
   ! DESCRIPTION:
   !> Free the list of entries.
   !> Set to private because should only be called internaly
   !
   !> @param [in,out] rec an associated pointer to an entry type.
   !> @return none
   !---------------------------------------------------------------------------
   RECURSIVE SUBROUTINE parser_free_entry(rec)
      IMPLICIT NONE
      TYPE(entry_type), POINTER, INTENT(INOUT) :: rec
     
      IF (associated(rec%next)) CALL parser_free_entry(rec%next)
      DEALLOCATE(rec)
      NULLIFY(rec)
     
      RETURN
   END SUBROUTINE parser_free_entry

   !---------------------------------------------------------------------------
   !> @author
   !> Guillaume Balarac, LEGI.
   !
   ! DESCRIPTION:
   !> Initialize the parser reseting all allocated datas.
   !> Set to private because should only be called when reading a new file
   !> by parser_parsefile
   !
   !> @return none
   !---------------------------------------------------------------------------
   SUBROUTINE parser_reset()
      IMPLICIT NONE
     
      nfields = 0
      IF (associated(entries)) THEN
         CALL parser_free_entry(entries)
         NULLIFY(entries)
      END IF
      entries=>null()
     
      RETURN
   END SUBROUTINE parser_reset

   !---------------------------------------------------------------------------
   !> @author
   !> Patrick BEGOU, LEGI.
   !
   ! DESCRIPTION:
   !> Print the key / values stored in the array.
   !> This subroutine is mainly for debuging, showing the content of the
   !> array of parameters.
   !
   !> @param [in] mytag the parameter tag
   !> @return none
   !---------------------------------------------------------------------------

   SUBROUTINE parser_showtag(mytag)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) ::mytag
      TYPE(entry_type), POINTER    ::myfield
      LOGICAL :: isdef 
      
      CALL parser_fieldfortag(mytag,myfield,isdef)
      IF(isdef) &
      & WRITE(6,'(a30,1x,":",a20)') TRIM(myfield%tag),TRIM(myfield%value)
   END SUBROUTINE parser_showtag

   !---------------------------------------------------------------------------
   !> @author
   !> Patrick BEGOU, LEGI.
   !
   ! DESCRIPTION:
   !> Print recursively the key / values stored in the structure.
   !
   !> @param [in] myfield the pointer to the structure list
   !> @return none
   !---------------------------------------------------------------------------  

   RECURSIVE SUBROUTINE parser_showAlltag(myfield)
      IMPLICIT NONE
      TYPE(entry_type), POINTER, INTENT(IN)    ::myfield
      
      IF(.NOT. ASSOCIATED(myfield)) RETURN
      WRITE(6,'("|",a40,1x,":",a75,"|")') TRIM(myfield%tag),TRIM(myfield%value)
      CALL parser_showAlltag(myfield%next)
   END SUBROUTINE parser_showAlltag

   !---------------------------------------------------------------------------  
   !> @author 
   !> Patrick BEGOU, LEGI.
   !
   ! DESCRIPTION: 
   !> Print the key / values stored in the array.
   !> This subroutine is mainly for debuging, showing the content of the
   !> array of parameters. 
   !
   !> @return none
   !---------------------------------------------------------------------------  

   SUBROUTINE parser_show()
      IMPLICIT NONE
      WRITE(6,'(119("="))')
      CALL parser_showAlltag(entries)
      WRITE(6,'(119("="))')
   END SUBROUTINE parser_show

   !---------------------------------------------------------------------------  
   !> @author 
   !> Guillaume Balarac, Patrick Begou LEGI.
   !
   ! DESCRIPTION: 
   !> Read and parse the input file.
   !> This subroutine processes the input file and stores key / values in a large
   !> array of strings. 
   !> 
   ! 
   !> @param [in]  input  The name of the input file
   !> @return .TRUE. if no error, else retruns .FALSE. 
   !---------------------------------------------------------------------------  
   LOGICAL FUNCTION parser_parsefile(input,OK)

     USE fileio
     
     IMPLICIT NONE
     CHARACTER(LEN=*), INTENT(IN)   :: input
     LOGICAL, INTENT(OUT), OPTIONAL ::OK
     
     INTEGER :: iunit,ierr,limiter,nlines,i,j,ntags,comment
     INTEGER, DIMENSION(:), ALLOCATABLE :: limit,line
     CHARACTER(LEN=line_length) :: buffer
     CHARACTER(LEN=line_length), DIMENSION(:), ALLOCATABLE :: file
     CHARACTER(LEN=line_length) :: value
     CHARACTER(LEN=tag_length)  :: tag
     
     parser_parsefile=.FALSE.
     
     !CALL parser_reset()
     
     ! Open the file
     ierr = 0
     iunit = iopen()
     OPEN (iunit,FILE=input,FORM='formatted',STATUS='old',IOSTAT=ierr)
     IF (ierr .NE. 0) THEN
         WRITE(6,'(a)')"Parser : unable to find input file ["//TRIM(input)//"]"
     RETURN
     ENDIF
     
     ! Count the number of lines in the file
     ierr = 0
     nlines = 0
     DO WHILE (ierr .EQ. 0)
        READ(iunit,'(a)',IOSTAT=ierr) buffer
        nlines = nlines + 1
     END DO
     REWIND(iunit)
     
     ! Allocate to the right size
     ALLOCATE(file(nlines+1),limit(nlines+1),line(nlines+1))
     
     ! Read everything in the buffer and put it in the file array
     ierr = 0
     nlines = 0
     DO WHILE (ierr .EQ. 0)
        READ(iunit,'(a)',IOSTAT=ierr) buffer
        IF (ierr.NE.0) EXIT 
        ! Remove the tabs
        DO j=1,line_length
           IF (ichar(buffer(j:j)).EQ.9) buffer(j:j)=' '
        END DO
        ! Find comments
        comment = SCAN(buffer,'!#%')
        ! Remove them
        IF (comment.NE.0) buffer(comment:) = ''
        ! Trim
        buffer = ADJUSTL(buffer)
        ! Add line
        IF (LEN_TRIM(buffer).NE.0) THEN
           nlines = nlines + 1
           file(nlines) = TRIM(buffer)
        END IF
     ENDDO 
     
     ! Close de file
     CLOSE(iunit)
     ierr = iclose(iunit)
     
     ! Get the tags
     ntags = 0
     DO i=1,nlines
        limiter = index(file(i),':')
        IF (limiter.NE.0) THEN
           ntags = ntags + 1
           line(ntags) = i
           line(ntags+1) = nlines+1
           limit(ntags) = limiter
        END IF
     END DO
     
     ! Read everything now
     DO i=1,ntags
        buffer = ''
        DO j=line(i),line(i+1)-1
           IF (j==line(i)) THEN
              buffer = trim(buffer) // trim(file(j))
           ELSE
              buffer = trim(buffer) // ' ' // trim(file(j))
           END IF
        END DO
        READ(buffer(1:limit(i)-1),'(a)') tag
        READ(buffer(limit(i)+1:),'(a)') value
        IF (len_trim(value).NE.0) THEN
           value = adjustl(value)
           IF(.NOT. parser_newentry(TRIM(tag),TRIM(value))) THEN
             IF (PRESENT(OK)) THEN
                OK=.FALSE.
                RETURN
             ELSE
                STOP 'FAILED'
             ENDIF
           ENDIF
        END IF
     END DO
     
     DEALLOCATE(file,limit,line)
     parser_parsefile=.TRUE.
   
     RETURN
   END FUNCTION parser_parsefile


   !---------------------------------------------------------------------------
   !> @author
   !> Guillaume Balarac LEGI.
   !
   ! DESCRIPTION:
   !> Add an additionnal parameter (Tag + Value) in the parameters structure.
   !> If a value exist for this key word, the new value replace the previous one.
   !> Internal use only
   !
   !> @param [in] mytag the parameter tag
   !> @param [in] myvalue the parameter value
   !---------------------------------------------------------------------------

   FUNCTION parser_newentry(mytag,myvalue)
    IMPLICIT NONE
    CHARACTER(*),INTENT(IN) :: mytag
    CHARACTER(*),INTENT(IN) :: myvalue
    LOGICAL :: parser_newentry

    LOGICAL :: isdef
    INTEGER ::ires
    TYPE(entry_type), POINTER    ::ifield

    parser_newentry=.FALSE.
    CALL parser_fieldfortag(mytag,ifield,isdef)
    IF (.NOT. isdef) THEN
       IF(.NOT. ASSOCIATED(entries)) THEN
          ALLOCATE(entries,stat=ires)
          IF (ires .NE.0) THEN
              WRITE(6,'(a)')"[ERROR] Parser: not enought memory for tag=[" &
	      &//TRIM(mytag)//"] and value= ["//TRIM(myvalue)//"]"
              RETURN
          ELSE
             entries%tag=mytag
             entries%value=myvalue
             entries%next=>null()
             current=>entries
          ENDIF
       ELSE
          ALLOCATE(ifield,stat=ires)
          IF (ires .NE.0) THEN
              WRITE(6,'(a)')"[ERROR] Parser: not enought memory for tag=[" &
	      &//TRIM(mytag)//"] and value= ["//TRIM(myvalue)//"]"
              RETURN
          ELSE
              ifield%tag=mytag
              ifield%value=myvalue
              ifield%next=>null()
              current%next=>ifield
              current=>ifield
          ENDIF
       ENDIF
    ELSE
       ifield%value=myvalue
    END IF
    parser_newentry=.TRUE.
    
    RETURN
   END FUNCTION parser_newentry
  
   !---------------------------------------------------------------------------  
   !> @author 
   !> Guillaume Balarac LEGI.
   !
   ! DESCRIPTION: 
   !> Search in the current parameters structure the occurence for the tag mytag.
   !> If found returns it's position (index in the array) in myfield argument and set
   !> the optional isdef argument to .TRUE. If the isdef argument is not provided,
   !> the returned value in myfield is undefine when the tag is not found
   !> Internal use only
   ! 
   !> @param [in] mytag The key word to search for 
   !> @param [out] myfield The position of the tag/value entry if it exists. Nndefined otherwise.
   !> @param [out] isdef .TRUE. if mytag exist, else .FALSE. 
   !---------------------------------------------------------------------------  
  SUBROUTINE parser_fieldfortag(mytag,myfield,isdef)
    IMPLICIT NONE
    CHARACTER(*),INTENT(IN)               :: mytag
    TYPE(entry_type), POINTER,INTENT(OUT) :: myfield
    LOGICAL,OPTIONAL,INTENT(OUT)          :: isdef

    LOGICAL :: OK

    OK=.FALSE.

    IF(ASSOCIATED(entries)) THEN
      myfield=>entries
      ! First element
      IF(myfield%tag .EQ. mytag) OK=.TRUE.

      DO WHILE(.NOT. OK .AND. ASSOCIATED(myfield%next))
         myfield=>myfield%next
         IF(myfield%tag .EQ. mytag) OK=.TRUE.
      END DO
    ENDIF
    IF(.NOT. OK) myfield=>null()
    IF(PRESENT(isdef)) isdef=OK

    RETURN
  END SUBROUTINE parser_fieldfortag

   !---------------------------------------------------------------------------  
   !> @author 
   !> Guillaume Balarac LEGI.
   !
   ! DESCRIPTION:
   !> Search in the current parameters array for the occurence for the key word provided in mytag.
   !> If found returns .TRUE. in the isdef argument. Otherwise sets isdef argument to .FALSE.
   !> Internal use only
   !
   !> @param [in] mytag the key word to search for
   !> @return .TRUE. if mytag exist, else .FALSE.
   !---------------------------------------------------------------------------
   ! Check whether the field is defined -------------------------------------
  LOGICAL FUNCTION parser_is_defined(mytag)
    IMPLICIT NONE

    CHARACTER(*),INTENT(IN) :: mytag

    TYPE(entry_type), POINTER :: ifield
    LOGICAL :: isdef

    CALL parser_fieldfortag(mytag,ifield,isdef)
    parser_is_defined=isdef

    RETURN
  END FUNCTION parser_is_defined
  
   !---------------------------------------------------------------------------
   !> @author
   !> Guillaume Balarac, LEGI.
   !
   ! DESCRIPTION:
   !> Private implementation of the generic parser_read() interface
   !> for logical datas.
   !
   !> @param [in]  mytag  The key word of this parameter
   !> @param [out]  value  The value of this parameter read in the input file
   !> @return none
   !---------------------------------------------------------------------------
   SUBROUTINE parser_readlogical(mytag,value)
     IMPLICIT NONE

     CHARACTER(*),INTENT(IN) :: mytag
     LOGICAL,INTENT(OUT) :: value

     TYPE(entry_type), POINTER :: ifield
     LOGICAL :: isdef

     CALL parser_fieldfortag(mytag,ifield,isdef)
     IF (.NOT.isdef) THEN
        PRINT*,'Parser : '// mytag //' not defined'
        STOP
     ELSE
        READ(ifield%value,*) value
     END IF

     RETURN
   END SUBROUTINE parser_readlogical

   !---------------------------------------------------------------------------  
   !> @author 
   !> Guillaume Balarac, LEGI.
   !
   ! DESCRIPTION:
   !> Private implementation of the generic parser_read() interface
   !> for INTEGER datas.
   !
   !> @param [in]  mytag  The key word of this parameter
   !> @param [out]  value  The value of this parameter read in the input file
   !> @return none
   !---------------------------------------------------------------------------
  SUBROUTINE parser_readint(mytag,value)
    IMPLICIT NONE

    CHARACTER(*),INTENT(IN) :: mytag
    INTEGER,INTENT(OUT) :: value

    TYPE(entry_type), POINTER :: ifield
    LOGICAL :: isdef

    CALL parser_fieldfortag(mytag,ifield,isdef)
    IF (.NOT.isdef) THEN
       PRINT*,'Parser : '// mytag //' not defined'
       STOP
    ELSE
       READ(ifield%value,*) value
    END IF
    
    RETURN
  END SUBROUTINE parser_readint
  
  
   !---------------------------------------------------------------------------  
   !> @author 
   !> Guillaume Balarac, LEGI.
   !
   ! DESCRIPTION: 
   !> Private implementation of the generic parser_read() interface 
   !> for REAL datas. 
   !
   !> @param [in]  mytag  The key word of this parameter
   !> @param [out]  value  The value of this parameter read in the input file
   !> @return none
   !---------------------------------------------------------------------------  
   ! Read floats ---------------------------------------------
  SUBROUTINE parser_readfloat(mytag,value)
    IMPLICIT NONE
    
    CHARACTER(*),INTENT(IN) :: mytag
    REAL(WP),INTENT(OUT) :: value

    TYPE(entry_type), POINTER :: ifield
    LOGICAL :: isdef
    
    CALL parser_fieldfortag(mytag,ifield,isdef)
    IF (.NOT.isdef) THEN
       PRINT*,'Parser : '// mytag //' not defined'
       STOP
    ELSE
       READ(ifield%value,*) value
    END IF
    
    RETURN
  END SUBROUTINE parser_readfloat

  
   !---------------------------------------------------------------------------  
   !> @author 
   !> Guillaume Balarac, LEGI.
   !
   ! DESCRIPTION: 
   !> Private implementation of the generic parser_read() interface 
   !> for STRING of character. 
   !
   !> @param [in]  mytag  The key word of this parameter
   !> @param [out]  value  The value of this parameter read in the input file
   !> @return none
   !---------------------------------------------------------------------------  
   ! Read characters ---------------------------------------------
  SUBROUTINE parser_readchar(mytag,value)
    IMPLICIT NONE
    
    CHARACTER(*),INTENT(IN) :: mytag
    CHARACTER(LEN=*),INTENT(OUT) :: value

    TYPE(entry_type), POINTER :: ifield
    LOGICAL :: isdef
    
    CALL parser_fieldfortag(mytag,ifield,isdef)
    IF (.NOT.isdef) THEN
       PRINT*,'Parser : '// mytag //' not defined'
       STOP
    ELSE
       READ(ifield%value,'(a)') value
    END IF
    
    RETURN
  END SUBROUTINE parser_readchar
  
  
   !---------------------------------------------------------------------------  
   !> @author 
   !> Guillaume Balarac, LEGI.
   !
   ! DESCRIPTION: 
   !> Count the number of elements in the value field of the entry corresponding
   !> to the key word.
   !
   !> @param [in]  mytag  The key word of this parameter
   !> @param [out]  numb the number of fields
   !> @return none
   !---------------------------------------------------------------------------  
  ! Count size of the arrays ----------------------------------------
  SUBROUTINE parser_getsize(mytag,numb)
    IMPLICIT NONE
    
    CHARACTER(*),INTENT(IN) :: mytag
    INTEGER,INTENT(OUT) :: numb

    TYPE(entry_type), POINTER :: ifield
    LOGICAL :: isdef
    INTEGER :: i
    INTEGER, DIMENSION(line_length) :: counter
    
    ! Read it now                                                                                    
    CALL parser_fieldfortag(mytag,ifield,isdef)
    IF (.NOT. isdef) THEN
       PRINT*,'Parser : '// mytag //' not defined'
       STOP
    END IF
    
    ! Count the number of entries                                                                    
    counter = 0
    DO i=1,len_trim(ifield%value)
       IF (ifield%value(i:i).EQ.' ') counter(i)=1
    END DO
    DO i=1+1,len_trim(ifield%value)
       IF (counter(i).EQ.1 .AND. counter(i-1).EQ.1) counter(i-1)=0
    END DO
    numb = sum(counter)+1
    
    RETURN
  END SUBROUTINE parser_getsize
  
  
   !---------------------------------------------------------------------------  
   !> @author 
   !> Guillaume Balarac, Patrick Begou LEGI.
   !
   ! DESCRIPTION: 
   !> Read a list of integer values in the value field of the entry corresponding
   !> to the key word. The array size should correspond exactly to the integer items number.
   !
   !> @param [in]  mytag  The key word of this parameter
   !> @param [out]  value an integer array with the read values.
   !> @return none
   !---------------------------------------------------------------------------  
    ! Read integer arrays ---------------------------------------------   
  SUBROUTINE parser_readintarray(mytag,value)
    IMPLICIT NONE

    CHARACTER(*),INTENT(IN)            :: mytag
    INTEGER,DIMENSION(:),INTENT(INOUT) :: value
    
    INTEGER :: nbval
    TYPE(entry_type), POINTER :: ifield
    LOGICAL :: isdef
    
    ! Read them
    CALL parser_fieldfortag(mytag,ifield,isdef)
    IF (.NOT. isdef) THEN
       PRINT*,'Parser : '// mytag //' not defined'
       STOP
    END IF
    
    CALL parser_getsize(mytag,nbval)
    IF (nbval .NE. SIZE(value)) THEN
       WRITE(6,'(2(a,i0))')'[ERROR] Parser: Try to get ',SIZE(value), &
                          &' values and item ['//mytag//'] has ',nbval
       STOP 'FAILED'
    END IF
    
    READ(ifield%value,*,err=100) value

    RETURN
    
100 WRITE(6,'(a,i0,a)')'[ERROR] Parser: Failed to read ',SIZE(value), &
        &' integer values in item ['//mytag//'] which is ['//TRIM(ifield%value)//']'
    STOP 'FAILED'

  END SUBROUTINE parser_readintarray

  
   !---------------------------------------------------------------------------  
   !> @author 
   !> Guillaume Balarac, LEGI.
   !
   ! DESCRIPTION: 
   !> Read a list of real values in the value field of the entry corresponding
   !> to the key word. 
   !> WARNING this routine is not secure, it suppose the array provided has the right size!
   !
   !> @param [in]  mytag  The key word of this parameter
   !> @param [out]  value a real array with the read values.
   !> @return none
   !---------------------------------------------------------------------------  
  ! Read float arrays --------------------------------------------- 
  SUBROUTINE parser_readfloatarray(mytag,value)
    IMPLICIT NONE

    CHARACTER(*),INTENT(IN) :: mytag
    REAL(WP),DIMENSION(:),INTENT(INOUT) :: value

    TYPE(entry_type), POINTER  :: ifield
    INTEGER :: nbval
    LOGICAL :: isdef
    
    ! Read them
    CALL parser_fieldfortag(mytag,ifield,isdef)
    IF (.NOT. isdef) THEN
       PRINT*,'Parser : '// mytag //' not defined'
       STOP
    END IF   
                                                                   
    CALL parser_getsize(mytag,nbval)
    IF (nbval .NE. SIZE(value)) THEN
       WRITE(6,'(2(a,i0))')'[ERROR] Parser: Try to get ',SIZE(value), &
                          &' values and item ['//mytag//'] has ',nbval
       STOP 'FAILED'
    END IF
    
    READ(ifield%value,*,err=100) value

    RETURN
    
100 WRITE(6,'(a,i0,a)')'[ERROR] Parser: Failed to read ',SIZE(value), &
        &' real values in item ['//mytag//'] which is ['//TRIM(ifield%value)//']'
    STOP 'FAILED'
    RETURN
  END SUBROUTINE parser_readfloatarray

  
   !---------------------------------------------------------------------------  
   !> @author 
   !> Guillaume Balarac, LEGI.
   !
   ! DESCRIPTION: 
   !> Read a list of real values in the value field of the entry corresponding
   !> to the key word. 
   !> WARNING this routine is not secure, it suppose the array provided has the right size!
   !
   !> @param [in]  mytag  The key word of this parameter
   !> @param [out]  value a 2D real array with the read values.
   !> @return none
   !---------------------------------------------------------------------------  
  SUBROUTINE parser_readfloatarray2D(mytag,value)
    IMPLICIT NONE
    
    CHARACTER(*),INTENT(IN) :: mytag
    REAL(WP),DIMENSION(:,:),INTENT(OUT) :: value

    TYPE(entry_type), POINTER :: ifield
    LOGICAL :: isdef
    
    ! Read them
    CALL parser_fieldfortag(mytag,ifield,isdef)
    IF (.NOT. isdef) THEN
       PRINT*,'Parser : '// mytag //' not defined'
       STOP
    END IF
    READ(ifield%value,*) value
    
    RETURN
  END SUBROUTINE parser_readfloatarray2D
  
   !---------------------------------------------------------------------------  
   !> @author 
   !> Guillaume Balarac, LEGI.
   !
   ! DESCRIPTION: 
   !> Read a list of string values in the value field of the entry corresponding
   !> to the key word. 
   !> WARNING this routine is not secure, it suppose the array provided has the right size!
   !
   !> @param [in]  mytag  The key word of this parameter
   !> @param [out]  value a character(*) array with the read values.
   !> @return none
   !---------------------------------------------------------------------------  
  SUBROUTINE parser_readchararray(mytag,value)
    IMPLICIT NONE

    CHARACTER(*),INTENT(IN) :: mytag
    CHARACTER(*),DIMENSION(:),INTENT(OUT) :: value

    TYPE(entry_type), POINTER :: ifield
    LOGICAL :: isdef
    
    ! Read them
    CALL parser_fieldfortag(mytag,ifield,isdef)
    IF (.NOT. isdef) THEN
       PRINT*,'Parser : '// mytag //' not defined'
       STOP
    END IF                                                                   
    READ(ifield%value,*) value

    RETURN
  END SUBROUTINE parser_readchararray


  !------------------------------------------------------------------------------
  !> @author 
  !> Patrick BEGOU, LEGI
  !
  !>
  !> @details
  !> This function returns nx, ny,nz (3 positive integers) for the current tag
  !> and returns .TRUE. if all is OK, .FALSE. otherwise.
  !>@param  [in]  tag the key word to get the values from
  !>@param  [out] nx first integer value for the key word
  !>@param  [out] ny second integer value for the key word
  !>@param  [out] nz third integer value for the key word
  !>@param  [in]  default some default value for nx, ny and nz if the tag is not found.
  !>@return .TRUE. if 3 positive integer values found, .FALSE. otherwise
  !------------------------------------------------------------------------------
  LOGICAL FUNCTION GetNxNyNz(tag,nx,ny,nz, default)

  IMPLICIT NONE
  CHARACTER(len=*), INTENT(IN)  ::tag
  INTEGER, INTENT(OUT)          ::nx,ny,nz
  INTEGER, INTENT(IN), OPTIONAL :: default
  
  INTEGER:: n, nxyz(3)
  
    GetNxNyNz=.FALSE.
  
    IF (.NOT. parser_is_defined(tag)) then
        nxyz = default
        return
    ENDIF
        
    
    CALL parser_getsize(tag,n)
    SELECT CASE (n)
    CASE(1)
         CALL parser_read(tag,nx)
         ny=nx
         nz=nx
    CASE(3)
         CALL parser_read(tag,nxyz)
         nx=nxyz(1)
         ny=nxyz(2)
         nz=nxyz(3)
    CASE DEFAULT
        WRITE(6,'(a,i0)') &
        &'[ERROR] GetNxNyNz: waiting for 1 or 3 values for domain size and got ',n
        CALL parser_showtag(tag)
        RETURN 
    END SELECT
    IF (nx.LE.1 .OR. ny.LE.1 .OR. nz.LE.1) THEN
        WRITE(6,'(a,3(1x,i0))')'[ERROR] GetNxNyNz: values should be greater than 1: ',nx,ny,nz
        CALL parser_showtag(tag)
        RETURN
    ENDIF
    
    GetNxNyNz=.TRUE.
  END FUNCTION GetNxNyNz


  !------------------------------------------------------------------------------
  !> @author 
  !> Patrick BEGOU, LEGI
  !
  !>
  !> @details
  !> This function returns lx, ly,lz (3 positive reals) for the current tag
  !> and returns .TRUE. if all is OK, .FALSE. otherwise.
  !>@param  [in]  tag the key word to get the values from
  !>@param  [out] lx first float value for the key word
  !>@param  [out] ly second float value for the key word
  !>@param  [out] lz third float value for the key word
  !>@return .TRUE. if 3 positive float values found, .FALSE. otherwise
  !------------------------------------------------------------------------------
  LOGICAL FUNCTION GetLxLyLz(tag,lx,ly,lz)

  IMPLICIT NONE
  CHARACTER(len=*), INTENT(IN) ::tag
  REAL(WP), INTENT(OUT)        ::lx,ly,lz
  
  INTEGER  :: n
  REAL(WP) ::lxyz(3)
  
    GetLxLyLz=.FALSE.
  
    IF (.NOT. parser_is_defined(tag)) RETURN
    
    CALL parser_getsize(tag,n)
    SELECT CASE (n)
    CASE(1)
         CALL parser_read(tag,lx)
         ly=lx
         lz=lx
    CASE(3)
         CALL parser_read(tag,lxyz)
         lx=lxyz(1)
         ly=lxyz(2)
         lz=lxyz(3)
    CASE DEFAULT
        WRITE(6,'(a,i0)') &
        &'[ERROR] GetLxLyLz: waiting for 1 or 3 values for domain lenght and got ',n
        CALL parser_showtag(tag)
        RETURN 
    END SELECT
    IF (lx.LE.0.0 .OR. ly.LE.0.0 .OR. lz.LE.0.0) THEN
        WRITE(6,'(a,3(1x,f10.6))')'[ERROR] GetLxLyLz: values should be strictly positive: ',lx,ly,lz
        CALL parser_showtag(tag)
        RETURN
    ENDIF
    
    GetLxLyLz=.TRUE.
  END FUNCTION GetLxLyLz


!---------------------------------------------------------------------------
!> @details 
!>    This routine check the type of filetring choosen for models which need 
!> @author = Antoine Vollant, LEGI
!> @param[inout]     valFilter contains the number of filter to use 
!> @param[in]        spec_rank is the number of processus in cpu pool 
!> @param[in]        if the routine check the filter for scalars ,sca is the
!>                   number of scalar (optional) 
!---------------------------------------------------------------------------
 function parseFilter(valFilter,spec_rank,sca,other) result(res)

   !I/O
   integer,intent(inout)                         :: valFilter
   integer,intent(in)                            :: spec_rank
   integer,intent(in),optional                   :: sca
   character(len=*),intent(in),optional          :: other 
   logical                                       :: res

   !Local data
   character(str_long)           :: request
   character(str_long)           :: namefiltervel 

   res =.true.

   if (present(other) .and. present(sca)) then
     write(request,'(a,i0)') trim(adjustl(other)),sca
   elseif ((present(sca)).and.(.not. present(other))) then
     write(request,'(a,i0)') 'Filter for scalar model ',sca
   elseif (.not.present(sca) .and. present(other)) then
     write(request,'(a)') trim(adjustl(other))
   else
     write(request,'(a)') 'Filter for velocity model'
   endif
   if (parser_is_defined(request)) then
     call parser_read(request,namefiltervel)
     if ( (namefiltervel .eq. 'Cutoff') .or. (namefiltervel .eq. 'cutoff')  ) then
       valFilter = 1
     elseif ( (namefiltervel .eq. 'Box') .or. (namefiltervel .eq. 'box') ) then
       valFilter = 2
     elseif ( (namefiltervel .eq. 'Gaussian') .or. (namefiltervel .eq. 'gaussian') ) then
       valFilter = 3
     elseif ( (namefiltervel .eq. 'Cutoff+Box') .or. (namefiltervel .eq. 'cutoff+box') ) then 
       valFilter = 4
     else
       if (spec_rank .eq. 0) then
        if (  present (sca) ) then
          write(6,'(a,1x,i0,1x,a,a)')'[ERROR] you made a spelling mistake for&
                               & the "Filter for scalar model',sca,'". You said ',namefiltervel
        else
          write(6,'(a,1x,a)')'[ERROR] you made a spelling mistake for&
                               & the "Filter for velocity" You said ',namefiltervel
        endif
         write(6,'(a)')'[Info]'
         write(6,'(a)')'[Info]  Please choose between:'
         write(6,'(a)')'[Info]  -> Box'
         write(6,'(a)')'[Info]  -> Cutoff'
         write(6,'(a)')'[Info]  -> Gaussian'
         write(6,'(a)')'[Info]  -> Cutoff+Box'
         write(6,'(a)')'[Info]'
       endif
       res = .false.
     endif
   else
     if (present(sca)) then
       if (spec_rank .eq. 0 ) write(6,'(a,1x,i0,a)')'[ERROR] You did not set &
                              & the "Filter for scalar model',sca,'"'
     else
       if (spec_rank .eq. 0 ) write(6,'(a,i0)')'[ERROR] You did not set the &
                              &"Filter for velocity model'
     endif
     if (spec_rank .eq. 0 ) then
       write(6,'(a)')'[Info]'
       write(6,'(a)')'[Info]  Please choose between:'
       write(6,'(a)')'[Info]  -> Box'
       write(6,'(a)')'[Info]  -> Cutoff'
       write(6,'(a)')'[Info]  -> Gaussian'
       write(6,'(a)')'[Info]  -> Cuttoff+Box'
       write(6,'(a)')'[Info]'
     endif
     valFilter = 0
     res = .false.
   endif

 end function 


END MODULE parser_tools

