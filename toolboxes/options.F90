!USEFORTEST toolbox
!USEFORTEST avgcond 
!USEFORTEST postprocess
!USEFORTEST advec
!USEFORTEST io
!USEFORTEST topo
      MODULE options
      
      IMPLICIT NONE

      ! Attributs
      INTEGER, PRIVATE :: myNbCpus1=-1
      CHARACTER(LEN=128), PRIVATE:: myInputFile='input'
      CHARACTER(LEN=128), PRIVATE:: myOutputFile='out'
      CHARACTER(LEN=128), PRIVATE:: prog=''
      LOGICAL, PRIVATE:: debugRandom=.FALSE.
      
      !Methodes publiques
      PUBLIC get_options
      PUBLIC GetNbCpus1
      PUBLIC GetInputFile
      PUBLIC GetOutputFile
      PUBLIC GetDebugRandom
      
      !Methodes publiques
      PRIVATE usage
      
      CONTAINS
      !-----------------------------------
      FUNCTION GetDebugRandom()
      !-----------------------------------
      IMPLICIT NONE
      LOGICAL:: GetDebugRandom
      GetDebugRandom=debugRandom
      RETURN
      END FUNCTION GetDebugRandom
      !-----------------------------------
      FUNCTION GetNbCpus1()
      !-----------------------------------
      IMPLICIT NONE
      INTEGER:: GetNbCpus1
      GetNbCpus1=myNbCpus1
      RETURN
      END FUNCTION GetNbCpus1

      !-----------------------------------
      FUNCTION GetInputFile()
      !-----------------------------------
      IMPLICIT NONE
      CHARACTER(LEN=128) GetInputFile
      GetInputFile=myInputFile
      RETURN
      END FUNCTION GetInputFile
      
      !-----------------------------------
      FUNCTION GetOutputFile()
      !-----------------------------------
      IMPLICIT NONE
      CHARACTER(LEN=128) GetOutputFile
      GetOutputFile=myOutputFile
      RETURN
      END FUNCTION GetOutputFile
      
      
      !-----------------------------------
      SUBROUTINE usage()
      !-----------------------------------
      IMPLICIT NONE
      
      WRITE(6,'(a)') ""  
      WRITE(6,'(a)') "Usage:"
      WRITE(6,'(a)') trim(prog)// " [-h] "
      WRITE(6,'(a)') trim(prog)// " [-i infile] [-o ofile] [-nx nbr]"
      WRITE(6,'(a)') "  -h         show this help information"
      WRITE(6,'(a)') "  -i  ifile  provides the name of the input file. "
      WRITE(6,'(a)') "             Default is: input. "
      WRITE(6,'(a)') "  -o  ofile  provides the base name for the output files. "
      WRITE(6,'(a)') "             Default is: out. "
      WRITE(6,'(a)') "             (To be implemented) "
      WRITE(6,'(a)') "  -nx nbr    sets the processus numbe to use in the X direction"
      WRITE(6,'(a)') "             for layout along Y or along Z."   
      WRITE(6,'(a)') "             It deactivates the distribution algorithm"  
      WRITE(6,'(a)') "  -random    sets the random number generator in debug mode."
      WRITE(6,'(a)') "             This allow to launch successive run with exactly the"   
      WRITE(6,'(a)') "             same output (exepted iteration elapsed time"  
      WRITE(6,'(a)') "             Usefull onfly for debug and comparison between two versions."
      WRITE(6,'(a)') ""  
      RETURN
      END SUBROUTINE

!-----------------------------------
!> The function get_options read the options  
!> associated with scaleExe file. If wrong options
!> are inserted the available full list of them is  
!> printed in the output file.
!----------------------------------


      !-----------------------------------
      FUNCTION get_options(myrank)
      !-----------------------------------
      IMPLICIT NONE
      LOGICAL          :: get_options
      INTEGER, INTENT(IN), OPTIONAL ::myrank
      INTEGER(kind=4)  :: i
      CHARACTER(LEN=80):: name
      INTEGER(kind=4)  :: iargc,rank,res
      
      get_options=.TRUE.
      rank=0
      IF (PRESENT(myrank)) rank=myrank
      
      ! le nom de cet executable tappe par l'utilisateur
      CALL getarg(0, prog)

      ! recherche des options
      i=1
      DO WHILE (i .LE. iargc() .AND. get_options)
       CALL getarg(i, name)
       IF (trim(name).EQ."-i") THEN
          i=i+1
          IF (i .GT.iargc()) THEN
             IF(rank .EQ. 0) CALL usage()
             get_options=.FALSE.
          ELSE
             CALL getarg(i,myInputFile )
             i=i+1
          ENDIF
       ELSE IF (trim(name).EQ."-o") THEN
          i=i+1
          IF (i .GT.iargc()) THEN
             IF(rank .EQ. 0) CALL usage()
             get_options=.FALSE.
          ELSE
             CALL getarg(i,myOutputFile )
             i=i+1
          ENDIF
       ELSE IF (trim(name).EQ."-nx") THEN
          i=i+1
          IF (i .GT.iargc()) THEN
             IF(rank .EQ. 0) CALL usage()
             get_options=.FALSE.
          ELSE
             CALL getarg(i,name)
             READ(name,*,iostat=res)myNbCpus1
             IF (res .EQ. 0) THEN 
               IF (myNbCpus1.LE. 0) THEN
                  get_options=.FALSE.
                  IF(rank .EQ. 0) CALL usage()
               ENDIF
             ELSE
               get_options=.FALSE.
             ENDIF
             i=i+1
          ENDIF
       ELSE IF (trim(name).EQ."-random") THEN
           IF (debugRandom) THEN
              IF(rank .EQ. 0) CALL usage()
              get_options=.FALSE.
           ELSE
              debugRandom=.TRUE. 
              i=i+1
           ENDIF
       ELSE IF (trim(name).EQ."-h") THEN
          IF(rank .EQ. 0) CALL usage()
           get_options=.FALSE.
      ELSE
           IF(rank .EQ. 0) CALL usage()
          get_options=.FALSE.
          RETURN
        ENDIF
      ENDDO
      
      RETURN
      END FUNCTION
      
      END MODULE options
