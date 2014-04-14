      FUNCTION ran3(idum)
!Returns a uniform random deviate betweena [0,1)
      
       USE precision_tools

       INTEGER      :: idum
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
       return
      END FUNCTION 

