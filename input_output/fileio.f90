!USEFORTEST toolbox
!USEFORTEST avgcond
!USEFORTEST postprocess
!USEFORTEST advec
!USEFORTEST io
!USEFORTEST topo
MODULE fileio
  IMPLICIT NONE

  INTEGER, PRIVATE, PARAMETER           :: firstfile=10
  INTEGER, PRIVATE, PARAMETER           :: lastfile=128
  INTEGER, PRIVATE, DIMENSION(firstfile:lastfile) :: iunits
  LOGICAL, PRIVATE                      :: initialized=.FALSE.

  CONTAINS

  ! ====================== !
  ! File index management: !
  !   - open a file        !
  !   - add it to the list !
  ! ====================== !
  INTEGER FUNCTION iopen()
    IMPLICIT NONE
    INTEGER ::i

    IF (.NOT. initialized) THEN
       iunits=0
       initialized=.TRUE.
    ENDIF
    i=firstfile
    DO WHILE (i.LE.lastfile)
       IF (iunits(i).EQ.0) EXIT
       i=i+1
    ENDDO

    IF (i.GT. lastfile) STOP "iopen: maximum units number reached"
    iunits(i)=1
    iopen=i
    RETURN
  END FUNCTION iopen

  ! ======================= !
  ! File index management:  !
  !   - close a file        !
  !   - remove it from list !
  ! ======================= !
  INTEGER FUNCTION iclose(iu)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: iu

    IF (.NOT. initialized) THEN
       STOP 'ERROR iclose(): trying a close before opening any file!'
    ENDIF
    IF (iu .ge.firstfile .and. iu .le. lastfile) then
       IF (iunits(iu).EQ.1) THEN
          iunits(iu)= 0
       ELSE
          WRITE(6,'(a,i0,a)')"WARNING iclose(): File n. ",iu,' no registred!'
       ENDIF
       iclose = iu
    ELSE
       iclose = -1
    ENDIF
    RETURN
  END FUNCTION iclose

END MODULE fileio
