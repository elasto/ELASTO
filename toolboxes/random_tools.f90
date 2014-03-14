!USEFORTEST toolbox
!USEFORTEST avgcond 
!> @addtogroup toolbox 
!! @{
!------------------------------------------------------------------------------
!
! MODULE: random
!
!> @author
!> Guillaume Balarac, LEGI
!
! DESCRIPTION: 
!> The aim of this module is to provide global functionalities for random numbers
!> generation
!------------------------------------------------------------------------------
MODULE random_tools

  USE precision_tools
  IMPLICIT NONE
  
  LOGICAL, PRIVATE::initialized=.FALSE.
  
CONTAINS
  
  !------------------------------------------------------------------------------
  !> @author 
  !> Guillaume Balarac, LEGI
  !
  !> \brief
  !> Returns a normally distributed pseudo-random number
  !
  !> @details
  !>  The function random_normal() returns a normally distributed       
  !>  pseudo-random number with zero mean and unit variance.            
  !>  The algorithm uses the ratio of uniforms method of A.J. Kinderman 
  !>  and J.F. Monahan augmented with quadratic bounding curves.        
  !>  Adapted from the following Fortran 77 code:                         
  !>      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.                 
  !>      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE, 
  !>      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.                  
  !
  !> @param[in] sd standard deviation
  !> @param[in] m mean
  !> @return a normally distributed pseudo-random number with zero mean and unit variance
  ! ------------------------------------------------------------------ !
  REAL(WP) FUNCTION random_normal(m,sd)
    IMPLICIT NONE
    
    REAL(WP), INTENT(IN), OPTIONAL :: m
    REAL(WP), INTENT(IN), OPTIONAL :: sd
    REAL(WP), PARAMETER :: s = 0.449871_WP
    REAL(WP), PARAMETER :: t = -0.386595_WP
    REAL(WP), PARAMETER :: a = 0.19600_WP
    REAL(WP), PARAMETER :: b = 0.25472_WP
    REAL(WP), PARAMETER :: r1 = 0.27597_WP
    REAL(WP), PARAMETER :: r2 = 0.27846_WP
    REAL(WP) :: u,v,x,y,q
    
    IF(.NOT. initialized) THEN
       WRITE(6,'(a)')'[WARNING] random_normal: Call to random_normal without a previous call to random_init.'
    ENDIF
    
    ! Generate P = (u,v) uniform in rectangle enclosing acceptance region
    DO    
       CALL RANDOM_NUMBER(u)
       CALL RANDOM_NUMBER(v)
       v=1.7156_WP*(v-0.5_WP)
       ! Evaluate the quadratic form
       x=u-s
       y=abs(v)-t
       q=x**2+y*(a*y-b*x)
       ! Accept P if inside inner ellipse
       if (q<r1) exit
       ! Reject P if outside outer ellipse
       if (q>r2) cycle
       ! Reject P if outside acceptance region
       if (v**2<-4.0_WP*log(u)*u**2) exit
    END DO
    
    ! Return ratio of P's coordinates as the normal deviate
    random_normal = v/u
    
    ! Modify to give correct mean and standard deviation
    IF (PRESENT(sd)) random_normal = random_normal*sd
    IF (PRESENT(m))  random_normal = random_normal+m
    
    RETURN
  END FUNCTION random_normal
  
  !------------------------------------------------------------------------------
  !> @author 
  !> Guillaume Balarac, LEGI
  !
  !> \brief
  !> Returns a log-normal distributed pseudo-random number
  !
  !> @details
  !> If X has a lognormal distribution, then log(X) is normally distributed.      
  !> Here the logarithm is the natural logarithm, that is to base e, sometimes    
  !> denoted as ln.  To generate random variates from this distribution, generate 
  !> a random deviate from the normal distribution with mean and variance equal   
  !> to the mean and variance of the logarithms of X, then take its exponential.  
  !>                                                                              
  !> Relationship between the mean & variance of log(X) and the mean & variance   
  !> of X, when X has a lognormal distribution.                                   
  !> Let m = mean of log(X), and s^2 = variance of log(X)                         
  !> Then                                                                         
  !> mean of X     = exp(m + 0.5s^2)                                              
  !> variance of X = (mean(X))^2.[exp(s^2) - 1]                                   
  !>                                                                              
  !> In the reverse direction                                                     
  !> variance of log(X) = log[1 + var(X)/(mean(X))^2]                             
  !> mean of log(X)     = log(mean(X) - 0.5var(log(X))                            
  !> @param[in] sd standard deviation
  !> @param[in] m mean
  !> @return a log-normal distributed pseudo-random number with zero mean and unit variance
  ! ---------------------------------------------------------------------------- !
  REAL(WP) FUNCTION random_lognormal(m,sd)
    IMPLICIT NONE
    
    REAL(WP), INTENT(IN) :: m
    REAL(WP), INTENT(IN) :: sd
    REAL(WP) :: x,mlog,sdlog
    
    sdlog = SQRT(LOG(1.0_WP+(sd/m)**2))
    mlog  = LOG(m)-0.5_WP*sdlog**2
    x     = random_normal(mlog,sdlog)
    random_lognormal = EXP(x)
    
    RETURN
  END FUNCTION random_lognormal
    
  !------------------------------------------------------------------------------
  !> @author 
  !> Guillaume Balarac, LEGI
  !
  !> \brief
  !> Initialization of the random number generator based on parallel partitioning
  !
  !> @param[in] rank the rank number of the current MPI process
  !------------------------------------------------------------------------------
  SUBROUTINE random_init(rank,debug)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: rank
    LOGICAL, OPTIONAL, INTENT(IN) :: debug
    
    INTEGER :: k
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed
    LOGICAL ::RunInDebug
    
    CALL RANDOM_SEED(size=k)
    ALLOCATE(seed(k))
    
    IF (PRESENT(debug)) THEN
       RunInDebug=debug
    ELSE 
       RunInDebug=.FALSE.
    ENDIF

    IF(RunInDebug) THEN
      seed=1
    ELSE
      CALL SYSTEM_CLOCK(COUNT=seed(1))
      seed(1:k) = seed(1:k)*(rank+1)
    ENDIF

    CALL RANDOM_SEED(PUT=seed)
    DEALLOCATE(seed)
    initialized=.TRUE.
    RETURN
  END SUBROUTINE random_init

!------------------------------------------------------------------------------
!> @author 
!> From "Numerical recepite", chapter 7 (Random numbers)
!
!> \brief
!> random number generator called ran3 in chapter 7 of "Numerical recepite" .
!> It returns a number distributed uniformely between [0,1).
!!
!> @param[in] rank of the the current MPI process
!------------------------------------------------------------------------------

FUNCTION ran3(idum)

     USE precision_tools

     INTEGER, INTENT(INOUT)      :: idum
     !INTEGER      :: MBIG,MSEED,MZ
     REAL(WP)     :: MBIG,MSEED,MZ
     REAL(WP)     :: ran3,FAC
     !PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
     PARAMETER(MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
     INTEGER      :: i,iff,ii,inext,inextp,k
     !INTEGER      :: mj,mk,ma(55)

     REAL(WP)     :: mj,mk,ma(55)
     SAVE iff,inext,inextp,ma
     DATA iff /0/

     IF(idum.lt.0.or.iff.eq.0) THEN
      iff=1
      mj=abs(MSEED-abs(idum))
      mj=mod(mj,MBIG)
      ma(55)=mj
      mk=1

      do i=1,54
        ii=mod(21*i,55)
        ma(ii)=mk
        mk=mj-mk
        if(mk.lt.MZ) mk=mk+MBIG
        mj=ma(ii)
      enddo

        do k=1,4
         do i =1,55
          ma(i)=ma(i)-ma(1+mod(i+30,55))
          if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
         enddo
        enddo
     
        inext=0
        inextp=31
        idum=1
      ENDIF

       inext=inext+1
       if(inext.eq.56) inext=1
       inextp=inextp+1
       if(inextp.eq.56)inextp=1
       mj=ma(inext)-ma(inextp)
       if(mj.lt.MZ)mj=mj+MBIG
       ma(inext)=mj
       ran3=mj*FAC

       RETURN
END FUNCTION ran3

END MODULE random_tools 
!> @}
