!> @addtogroup post_process
!! @{

!------------------------------------------------------------------------------

!
! MODULE: post_mhd
!
!
! DESCRIPTION: 
!>  This module provide a library of post-process routines, specifics to MHD simulations. This library
!! provide all that you need to build your own post-process.
!
!> @details 
!! 
!> @author
!! Mouloud KESSAR, LEGI
!
!------------------------------------------------------------------------------

module post_mhd_helic_lib

    use datalayout
    use toolbox 
    use stat_tools 
    use differential_tools 
    USE wavenumber_tools 
    use transforms_tools
    use subgrid_tools
    use differential_tools 
    use filtering_tools 
    use mhd_model
    use velmodels
    implicit none
    private

     ! ==== Public procedures ==== 
!     public          :: eqmhd_lib
!     public          :: EmhdGS_lib
!     public          :: aprioriMHDEm
    public          :: aprioriUDivtb !aprioriEQDmMHD
    public          :: aprioBrotDivtub
    public          :: aprioJDivtub
    public          :: aprioriBdivTb
    public          :: aprioriBdivTu
contains


!!!!!!Hcurr!!!

subroutine aprioBrotDivtub(Uin,Vin,Win,Bxin,Byin,Bzin,Ukin,Vkin,Wkin,&
                    & Bxkin,Bykin,Bzkin,Ufin,Vfin,Wfin,Bfxin,Bfyin,Bfzin,& 
                    & Wave,spec_rank,delta,filter,nbcpus,meanbrotdivTubexact, meanbrotdivTubSmag)

  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Uin,Vin,Win
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Bxin,Byin,Bzin 
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Ufin,Vfin,Wfin,Bfxin,Bfyin,Bfzin
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Ukin,Vkin,Wkin
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Bxkin,Bykin,Bzkin
  TYPE(WaveNumbers), INTENT(IN) :: Wave
  
  TYPE(COMPLEX_DATA_LAYOUT) :: div1,div2,div3
  TYPE(REAL_DATA_LAYOUT) :: div1R,div2R,div3R

  TYPE(COMPLEX_DATA_LAYOUT) :: rotdiv1,rotdiv2,rotdiv3
  TYPE(REAL_DATA_LAYOUT) :: rotdiv1R,rotdiv2R,rotdiv3R

  INTEGER :: spec_rank,nbcpus,filter
  TYPE(REAL_DATA_LAYOUT) :: T12,T13,T23
  TYPE(REAL_DATA_LAYOUT) :: Wtab
  logical :: res
  REAL(WP) :: delta
  REAL(WP) :: meanbrotdivTubexact, meanbrotdivTubSmag

  IF ((.NOT. copyStructOnly(Bxin,Wtab)).OR. &
     &(.NOT. copyStructOnly(Bxin,T12)) .OR. &
     &(.NOT. copyStructOnly(Bxin,T13)) .OR. &
     &(.NOT. copyStructOnly(Bxin,T23))  ) RETURN 

  IF ((.NOT. copyStructOnly(Bxkin,div1)) .OR. &
     &(.NOT. copyStructOnly(Bxkin,div2)) .OR. &
     &(.NOT. copyStructOnly(Bxkin,div3))   ) RETURN 


  IF ((.NOT. copyStructOnly(Bxin,div1R)) .OR. &
     &(.NOT. copyStructOnly(Bxin,div2R)) .OR. &
     &(.NOT. copyStructOnly(Bxin,div3R))  ) RETURN 

  IF ((.NOT. copyStructOnly(Bxkin,rotdiv1)) .OR. &
     &(.NOT. copyStructOnly(Bxkin,rotdiv2)) .OR. &
     &(.NOT. copyStructOnly(Bxkin,rotdiv3))   ) RETURN 


  IF ((.NOT. copyStructOnly(Bxin,rotdiv1R)) .OR. &
     &(.NOT. copyStructOnly(Bxin,rotdiv2R)) .OR. &
     &(.NOT. copyStructOnly(Bxin,rotdiv3R))  ) RETURN 

     !!!!div Exact
     call computeT_ijVecAVec_B( Ukin,Uin,Vin,Win,Ufin,Vfin,Wfin,&
                              & Bxin,Byin,Bzin,Bfxin,Bfyin,Bfzin,&
                              & wave,T12,T13,T23,filter,res,delta)

     call ComputeDivTijVec1Vec2(div1,div2,div3,T12,T13,T23,wave,res)

     call computeCurlKmain(wave,div1,div2,div3,rotdiv1,rotdiv2,rotdiv3)

           CALL btran(rotdiv1,rotdiv1R,res)
           IF (.NOT.res) RETURN
           CALL btran(rotdiv2,rotdiv2R,res)
           IF (.NOT.res) RETURN
           CALL btran(rotdiv3,rotdiv3R,res)
           IF (.NOT.res) RETURN


     Wtab%values = rotdiv1R%values*Bfxin%values + rotdiv2R%values*Bfyin%values + rotdiv3R%values*Bfzin%values
     meanbrotdivTubexact  = computeFieldAvg(Wtab, spec_rank,nbcpus)

     !!!!! div smag
     call smagorinskyMHDTij(T12,T13,T23,Bxkin,Bykin,Bzkin,&
                           &wave,spec_rank,res,delta)

     call ComputeDivTijVec1Vec2(div1,div2,div3,T12,T13,T23,wave,res)


     call computeCurlKmain(wave,div1,div2,div3,rotdiv1,rotdiv2,rotdiv3)

           CALL btran(rotdiv1,rotdiv1R,res)
           IF (.NOT.res) RETURN
           CALL btran(rotdiv2,rotdiv2R,res)
           IF (.NOT.res) RETURN
           CALL btran(rotdiv3,rotdiv3R,res)
           IF (.NOT.res) RETURN


!      Wtab%values = div1R%values*Bfxin%values + div2R%values*Bfyin%values + div3R%values*Bfzin%values
     Wtab%values = rotdiv1R%values*Bfxin%values + rotdiv2R%values*Bfyin%values + rotdiv3R%values*Bfzin%values

     meanbrotdivTubSmag  = computeFieldAvg(Wtab, spec_rank,nbcpus)

    CALL deleteDataLayout(T12)
    CALL deleteDataLayout(T13)
    CALL deleteDataLayout(T23)

    CALL deleteDataLayout(Wtab)

    CALL deleteDataLayout(div3)
    CALL deleteDataLayout(div2)
    CALL deleteDataLayout(div1)

    CALL deleteDataLayout(div3R)
    CALL deleteDataLayout(div2R)
    CALL deleteDataLayout(div1R)

    CALL deleteDataLayout(rotdiv3)
    CALL deleteDataLayout(rotdiv2)
    CALL deleteDataLayout(rotdiv1)

    CALL deleteDataLayout(rotdiv3R)
    CALL deleteDataLayout(rotdiv2R)
    CALL deleteDataLayout(rotdiv1R)

end subroutine aprioBrotDivtub

subroutine aprioJDivtub(Uin,Vin,Win,Bxin,Byin,Bzin,Ukin,Vkin,Wkin,&
                    & Bxkin,Bykin,Bzkin,Ufin,Vfin,Wfin,Bfxin,Bfyin,Bfzin,Jxf,Jyf,Jzf,& 
                    & Wave,spec_rank,delta,filter,nbcpus,meanJdivTubexact, meanJdivTubSmag)

  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Uin,Vin,Win
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Bxin,Byin,Bzin 
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Ufin,Vfin,Wfin,Bfxin,Bfyin,Bfzin
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Ukin,Vkin,Wkin
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Bxkin,Bykin,Bzkin
  TYPE(WaveNumbers), INTENT(IN) :: Wave
  
  TYPE(COMPLEX_DATA_LAYOUT) :: div1,div2,div3
  TYPE(REAL_DATA_LAYOUT) :: div1R,div2R,div3R

!   TYPE(COMPLEX_DATA_LAYOUT) :: Jxk,Jyk,Jzk
  TYPE(REAL_DATA_LAYOUT) :: Jxf,Jyf,Jzf

  INTEGER :: spec_rank,nbcpus,filter
  TYPE(REAL_DATA_LAYOUT) :: T12,T13,T23
  TYPE(REAL_DATA_LAYOUT) :: Wtab
  logical :: res
  REAL(WP) :: delta
  REAL(WP) :: meanJdivTubexact, meanJdivTubSmag

  IF ((.NOT. copyStructOnly(Bxin,Wtab)).OR. &
     &(.NOT. copyStructOnly(Bxin,T12)) .OR. &
     &(.NOT. copyStructOnly(Bxin,T13)) .OR. &
     &(.NOT. copyStructOnly(Bxin,T23))  ) RETURN 

  IF ((.NOT. copyStructOnly(Bxkin,div1)) .OR. &
     &(.NOT. copyStructOnly(Bxkin,div2)) .OR. &
     &(.NOT. copyStructOnly(Bxkin,div3))   ) RETURN 


  IF ((.NOT. copyStructOnly(Bxin,div1R)) .OR. &
     &(.NOT. copyStructOnly(Bxin,div2R)) .OR. &
     &(.NOT. copyStructOnly(Bxin,div3R))  ) RETURN 


     !!!!div Exact
     call computeT_ijVecAVec_B( Ukin,Uin,Vin,Win,Ufin,Vfin,Wfin,&
                              & Bxin,Byin,Bzin,Bfxin,Bfyin,Bfzin,&
                              & wave,T12,T13,T23,filter,res,delta)

     call ComputeDivTijVec1Vec2(div1,div2,div3,T12,T13,T23,wave,res)

           CALL btran(div1,div1R,res)
           IF (.NOT.res) RETURN
           CALL btran(div2,div2R,res)
           IF (.NOT.res) RETURN
           CALL btran(div3,div3R,res)
           IF (.NOT.res) RETURN
     Wtab%values = div1R%values*Jxf%values + div2R%values*Jyf%values + div3R%values*Jzf%values
     meanJdivTubexact  = computeFieldAvg(Wtab, spec_rank,nbcpus)

     !!!!! div smag
     call smagorinskyMHDTij(T12,T13,T23,Bxkin,Bykin,Bzkin,&
                           &wave,spec_rank,res,delta)

     call ComputeDivTijVec1Vec2(div1,div2,div3,T12,T13,T23,wave,res)
           CALL btran(div1,div1R,res)
           IF (.NOT.res) RETURN
           CALL btran(div2,div2R,res)
           IF (.NOT.res) RETURN
           CALL btran(div3,div3R,res)
           IF (.NOT.res) RETURN

!      Wtab%values = div1R%values*Bfxin%values + div2R%values*Bfyin%values + div3R%values*Bfzin%values
     Wtab%values = div1R%values*Jxf%values + div2R%values*Jyf%values + div3R%values*Jzf%values

     meanJdivTubSmag  = computeFieldAvg(Wtab, spec_rank,nbcpus)

    CALL deleteDataLayout(T12)
    CALL deleteDataLayout(T13)
    CALL deleteDataLayout(T23)

    CALL deleteDataLayout(Wtab)

    CALL deleteDataLayout(div3)
    CALL deleteDataLayout(div2)
    CALL deleteDataLayout(div1)

    CALL deleteDataLayout(div3R)
    CALL deleteDataLayout(div2R)
    CALL deleteDataLayout(div1R)

!     CALL deleteDataLayout(rotdiv3)
!     CALL deleteDataLayout(rotdiv2)
!     CALL deleteDataLayout(rotdiv1)
! 
!     CALL deleteDataLayout(rotdiv3R)
!     CALL deleteDataLayout(rotdiv2R)
!     CALL deleteDataLayout(rotdiv1R)

end subroutine aprioJDivtub


!!!!Hcross!!!

subroutine aprioUDivtub(Uin,Vin,Win,Bxin,Byin,Bzin,Ukin,Vkin,Wkin,&
                    & Bxkin,Bykin,Bzkin,Ufin,Vfin,Wfin,Bfxin,Bfyin,Bfzin,& 
                    & Wave,spec_rank,delta,filter,nbcpus,meanUdivTubexact, meanUdivTubSmag)

  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Uin,Vin,Win
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Bxin,Byin,Bzin 
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Ufin,Vfin,Wfin,Bfxin,Bfyin,Bfzin
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Ukin,Vkin,Wkin
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Bxkin,Bykin,Bzkin
  TYPE(WaveNumbers), INTENT(IN) :: Wave
  
  TYPE(COMPLEX_DATA_LAYOUT) :: div1,div2,div3
  TYPE(REAL_DATA_LAYOUT) :: div1R,div2R,div3R

!   TYPE(COMPLEX_DATA_LAYOUT) :: Jxk,Jyk,Jzk
!   TYPE(REAL_DATA_LAYOUT) :: Jxf,Jyf,Jzf

  INTEGER :: spec_rank,nbcpus,filter
  TYPE(REAL_DATA_LAYOUT) :: T12,T13,T23
  TYPE(REAL_DATA_LAYOUT) :: Wtab
  logical :: res
  REAL(WP) :: delta
  REAL(WP) :: meanUdivTubexact, meanUdivTubSmag

  IF ((.NOT. copyStructOnly(Bxin,Wtab)).OR. &
     &(.NOT. copyStructOnly(Bxin,T12)) .OR. &
     &(.NOT. copyStructOnly(Bxin,T13)) .OR. &
     &(.NOT. copyStructOnly(Bxin,T23))  ) RETURN 

  IF ((.NOT. copyStructOnly(Bxkin,div1)) .OR. &
     &(.NOT. copyStructOnly(Bxkin,div2)) .OR. &
     &(.NOT. copyStructOnly(Bxkin,div3))   ) RETURN 


  IF ((.NOT. copyStructOnly(Bxin,div1R)) .OR. &
     &(.NOT. copyStructOnly(Bxin,div2R)) .OR. &
     &(.NOT. copyStructOnly(Bxin,div3R))  ) RETURN 


     !!!!div Exact
     call computeT_ijVecAVec_B( Ukin,Uin,Vin,Win,Ufin,Vfin,Wfin,&
                              & Bxin,Byin,Bzin,Bfxin,Bfyin,Bfzin,&
                              & wave,T12,T13,T23,filter,res,delta)

     call ComputeDivTijVec1Vec2(div1,div2,div3,T12,T13,T23,wave,res)

           CALL btran(div1,div1R,res)
           IF (.NOT.res) RETURN
           CALL btran(div2,div2R,res)
           IF (.NOT.res) RETURN
           CALL btran(div3,div3R,res)
           IF (.NOT.res) RETURN
     Wtab%values = -div1R%values*Ufin%values -div2R%values*Vfin%values -div3R%values*Wfin%values
     meanUdivTubexact  = computeFieldAvg(Wtab, spec_rank,nbcpus)

     !!!!! div smag
     call smagorinskyMHDTij(T12,T13,T23,Bxkin,Bykin,Bzkin,&
                           &wave,spec_rank,res,delta)

     call ComputeDivTijVec1Vec2(div1,div2,div3,T12,T13,T23,wave,res)
           CALL btran(div1,div1R,res)
           IF (.NOT.res) RETURN
           CALL btran(div2,div2R,res)
           IF (.NOT.res) RETURN
           CALL btran(div3,div3R,res)
           IF (.NOT.res) RETURN

!      Wtab%values = div1R%values*Bfxin%values + div2R%values*Bfyin%values + div3R%values*Bfzin%values
     Wtab%values = -div1R%values*Ufin%values -div2R%values*Vfin%values -div3R%values*Wfin%values
     meanUdivTubSmag  = computeFieldAvg(Wtab, spec_rank,nbcpus)

    CALL deleteDataLayout(T12)
    CALL deleteDataLayout(T13)
    CALL deleteDataLayout(T23)

    CALL deleteDataLayout(Wtab)

    CALL deleteDataLayout(div3)
    CALL deleteDataLayout(div2)
    CALL deleteDataLayout(div1)

    CALL deleteDataLayout(div3R)
    CALL deleteDataLayout(div2R)
    CALL deleteDataLayout(div1R)



end subroutine aprioUDivtub

subroutine aprioriBDivtu(Uin,Vin,Win,Bxin,Byin,Bzin,Ukin,Vkin,Wkin,&
                    & Uinf,Vinf,Winf,Bxinf,Byinf,Bzinf,&
                    & Bxkin,Bykin,Bzkin,Wave,spec_rank,filter,delta,&
                    & meanBdivTuexact,meanBdivTuSmag,nbcpus)

  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Uin,Vin,Win
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Bxin,Byin,Bzin 
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Uinf,Vinf,Winf
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Bxinf,Byinf,Bzinf 
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Ukin,Vkin,Wkin
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Bxkin,Bykin,Bzkin
  TYPE(WaveNumbers), INTENT(IN) :: Wave
  
  TYPE(COMPLEX_DATA_LAYOUT) :: div1,div2,div3
  TYPE(REAL_DATA_LAYOUT) :: div1R,div2R,div3R

  INTEGER :: spec_rank,nbcpus,numVel,type_avg,filter
  TYPE(REAL_DATA_LAYOUT) :: T12,T13,T23
  TYPE(REAL_DATA_LAYOUT) :: T11,T22,T33
  TYPE(REAL_DATA_LAYOUT) :: S12,S13,S23
  TYPE(REAL_DATA_LAYOUT) :: S11,S22,S33
  TYPE(REAL_DATA_LAYOUT) :: Wtab
  logical :: res
  REAL(WP) :: delta!, meanEqdmTUuuSimil,meanEqdmTUubSimil
  REAL(WP) :: meanBdivTuexact,meanBdivTuSmag
!   REAL(WP) :: meanBdivTubGrad, meanEqdmTUubGrad

  IF ((.NOT. copyStructOnly(Bxin,Wtab)) .OR. &
     &(.NOT. copyStructOnly(Bxin,T12))  .OR. &
     &(.NOT. copyStructOnly(Bxin,T13))  .OR. &
     &(.NOT. copyStructOnly(Bxin,T23))  .OR. &
     &(.NOT. copyStructOnly(Bxin,T11))  .OR. &
     &(.NOT. copyStructOnly(Bxin,T22))  .OR. &
     &(.NOT. copyStructOnly(Bxin,T33))  ) RETURN 
 
  IF ((.NOT. copyStructOnly(Bxkin,div1)) .OR. &
     &(.NOT. copyStructOnly(Bxkin,div2)) .OR. &
     &(.NOT. copyStructOnly(Bxkin,div3))   ) RETURN 

  IF ((.NOT. copyStructOnly(Bxin,div1R)) .OR. &
     &(.NOT. copyStructOnly(Bxin,div2R)) .OR. &
     &(.NOT. copyStructOnly(Bxin,div3R))   ) RETURN 

     numVel=1
     type_avg=0
!      call computeStressTensor(Ukin,VKin,WKin,wave,S11,S12,S13,S22,S23,S33,res)

     call smagorinskyVel(div1,div2,div3,Uin,Vin,Win,Ukin,Vkin,Wkin,wave,numVel,spec_rank,res, &
                          &  filter,type_avg)

           CALL btran(div1,div1R,res)
           IF (.NOT.res) RETURN
           CALL btran(div2,div2R,res)
           IF (.NOT.res) RETURN
           CALL btran(div3,div3R,res)
           IF (.NOT.res) RETURN
    Wtab%values = -div1R%values*Bxinf%values -div2R%values*Byinf%values -div3R%values*Bzinf%values
    meanBdivTuSmag  = computeFieldAvg(Wtab, spec_rank,nbcpus)

    call computeT_ijVecA(Uin,Vin,Win,Uinf,Vinf,Winf,wave,T11,T12,T13, &
            & T22,T23,T33,filter,res,delta)

    call ComputeDivSymTijVec1Vec2(div1,div2,div3,T12,T13,T23,T11,T22,T33,wave,res)
           CALL btran(div1,div1R,res)
           IF (.NOT.res) RETURN
           CALL btran(div2,div2R,res)
           IF (.NOT.res) RETURN
           CALL btran(div3,div3R,res)
    Wtab%values = -div1R%values*Bxinf%values -div2R%values*Byinf%values -div3R%values*Bzinf%values
    meanBdivTuexact  = computeFieldAvg(Wtab, spec_rank,nbcpus)

    CALL deleteDataLayout(T12)
    CALL deleteDataLayout(T13)
    CALL deleteDataLayout(T23)
    CALL deleteDataLayout(T11)
    CALL deleteDataLayout(T22)
    CALL deleteDataLayout(T33)

    CALL deleteDataLayout(Wtab)

    CALL deleteDataLayout(div3)
    CALL deleteDataLayout(div2)
    CALL deleteDataLayout(div1)

end subroutine aprioriBdivTu

subroutine aprioriBDivtb(Bxin,Byin,Bzin,Bxinf,Byinf,Bzinf,&
                    & Bxkin,Bykin,Bzkin,Wave,spec_rank,filter,delta,&
                    & meanBdivTbexact,meanBdivTbSmag,nbcpus)

!   TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Uin,Vin,Win
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Bxin,Byin,Bzin 
!   TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Ukin,Vkin,Wkin
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Bxkin,Bykin,Bzkin
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Bxinf,Byinf,Bzinf 
  TYPE(COMPLEX_DATA_LAYOUT) :: Bxkinf,Bykinf,Bzkinf
  TYPE(WaveNumbers), INTENT(IN) :: Wave
  TYPE(COMPLEX_DATA_LAYOUT) :: div1,div2,div3
  TYPE(REAL_DATA_LAYOUT) :: div1R,div2R,div3R

  INTEGER :: spec_rank,nbcpus,numVel,type_avg,filter
  TYPE(REAL_DATA_LAYOUT) :: T12,T13,T23
  TYPE(REAL_DATA_LAYOUT) :: T11,T22,T33
  TYPE(REAL_DATA_LAYOUT) :: S12,S13,S23
  TYPE(REAL_DATA_LAYOUT) :: S11,S22,S33
  TYPE(REAL_DATA_LAYOUT) :: Wtab
  logical :: res
  REAL(WP) :: delta!, meanEqdmTUuuSimil,meanEqdmTUubSimil
  REAL(WP) :: meanBdivTbexact,meanBdivTbSmag
!   REAL(WP) :: meanBdivTubGrad, meanEqdmTUubGrad

  IF ((.NOT. copyStructOnly(Bxin,Wtab)) .OR. &
     &(.NOT. copyStructOnly(Bxin,T12))  .OR. &
     &(.NOT. copyStructOnly(Bxin,T13))  .OR. &
     &(.NOT. copyStructOnly(Bxin,T23))  .OR. &
     &(.NOT. copyStructOnly(Bxin,T11))  .OR. &
     &(.NOT. copyStructOnly(Bxin,T22))  .OR. &
     &(.NOT. copyStructOnly(Bxin,T33))  ) RETURN 
 
  IF ((.NOT. copyStructOnly(Bxkin,div1)) .OR. &
     &(.NOT. copyStructOnly(Bxkin,div2)) .OR. &
     &(.NOT. copyStructOnly(Bxkin,div3))   ) RETURN 

  IF ((.NOT. copyStructOnly(Bxin,div1R)) .OR. &
     &(.NOT. copyStructOnly(Bxin,div2R)) .OR. &
     &(.NOT. copyStructOnly(Bxin,div3R))   ) RETURN 

     numVel=1
     type_avg=0

     call smagorinskyVel(div1,div2,div3,Bxinf,Byinf,Bzinf,Bxkinf,Bykinf,Bzkinf,wave,numVel,spec_rank,res, &
                          &  filter,type_avg)


    Wtab%values = -div1R%values*Bxinf%values -div2R%values*Byinf%values -div3R%values*Bzinf%values
    meanBdivTbSmag  = computeFieldAvg(Wtab, spec_rank,nbcpus)

    call computeT_ijVecA(Bxin,Byin,Bzin,Bxinf,Byinf,Bzinf,wave,T11,T12,T13, &
            & T22,T23,T33,filter,res,delta)

    call ComputeDivSymTijVec1Vec2(div1,div2,div3,T12,T13,T23,T11,T22,T33,wave,res)
    Wtab%values = -div1R%values*Bxinf%values -div2R%values*Byinf%values -div3R%values*Bzinf%values
    meanBdivTbexact  = computeFieldAvg(Wtab, spec_rank,nbcpus)

    CALL deleteDataLayout(T12)
    CALL deleteDataLayout(T13)
    CALL deleteDataLayout(T23)
    CALL deleteDataLayout(T11)
    CALL deleteDataLayout(T22)
    CALL deleteDataLayout(T33)

    CALL deleteDataLayout(Wtab)


    CALL deleteDataLayout(div3)
    CALL deleteDataLayout(div2)
    CALL deleteDataLayout(div1)

end subroutine aprioriBdivTb

subroutine aprioriUDivtb(Bxin,Byin,Bzin,Bxinf,Byinf,Bzinf,&
                    & Bxkin,Bykin,Bzkin,Wave,spec_rank,filter,delta,meanUdivTbexact,&
                    & meanUdivTbSmag,nbcpus)

  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Bxin,Byin,Bzin 
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Bxkin,Bykin,Bzkin
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Bxinf,Byinf,Bzinf 
  TYPE(WaveNumbers), INTENT(IN) :: Wave
  
  TYPE(COMPLEX_DATA_LAYOUT) :: div1,div2,div3
  TYPE(REAL_DATA_LAYOUT) :: div1R,div2R,div3R

  INTEGER :: spec_rank,nbcpus,numVel,type_avg,filter
  TYPE(REAL_DATA_LAYOUT) :: T12,T13,T23
  TYPE(REAL_DATA_LAYOUT) :: T11,T22,T33
  TYPE(REAL_DATA_LAYOUT) :: S12,S13,S23
  TYPE(REAL_DATA_LAYOUT) :: S11,S22,S33
  TYPE(REAL_DATA_LAYOUT) :: Wtab
  logical :: res
  REAL(WP) :: delta
  REAL(WP) :: meanUdivTbexact,meanUdivTbSmag

  IF ((.NOT. copyStructOnly(Bxin,Wtab)) .OR. &
     &(.NOT. copyStructOnly(Bxin,T12))  .OR. &
     &(.NOT. copyStructOnly(Bxin,T13))  .OR. &
     &(.NOT. copyStructOnly(Bxin,T23))  .OR. &
     &(.NOT. copyStructOnly(Bxin,T11))  .OR. &
     &(.NOT. copyStructOnly(Bxin,T22))  .OR. &
     &(.NOT. copyStructOnly(Bxin,T33))  ) RETURN 
 
  IF ((.NOT. copyStructOnly(Bxkin,div1)) .OR. &
     &(.NOT. copyStructOnly(Bxkin,div2)) .OR. &
     &(.NOT. copyStructOnly(Bxkin,div3))   ) RETURN 

  IF ((.NOT. copyStructOnly(Bxin,div1R)) .OR. &
     &(.NOT. copyStructOnly(Bxin,div2R)) .OR. &
     &(.NOT. copyStructOnly(Bxin,div3R))   ) RETURN 

     numVel=1
     type_avg=0

     call smagorinskyVel(div1,div2,div3,Bxin,Byin,Bzin,Bxkin,Bykin,Bzkin,wave,numVel,spec_rank,res, &
                          &  filter,type_avg)


    Wtab%values = -div1R%values*Bxinf%values -div2R%values*Byinf%values -div3R%values*Bzinf%values
    meanUdivTbSmag  = computeFieldAvg(Wtab, spec_rank,nbcpus)

    call computeT_ijVecA(Bxin,Byin,Bzin,Bxinf,Byinf,Bzinf,wave,T11,T12,T13,T22,T23,T33,filter,res,delta)

    call ComputeDivSymTijVec1Vec2(div1,div2,div3,T12,T13,T23,T11,T22,T33,wave,res)
    Wtab%values = -div1R%values*Bxinf%values -div2R%values*Byinf%values -div3R%values*Bzinf%values
    meanUdivTbexact = computeFieldAvg(Wtab, spec_rank,nbcpus)

    CALL deleteDataLayout(T12)
    CALL deleteDataLayout(T13)
    CALL deleteDataLayout(T23)
    CALL deleteDataLayout(T11)
    CALL deleteDataLayout(T22)
    CALL deleteDataLayout(T33)

    CALL deleteDataLayout(Wtab)


    CALL deleteDataLayout(div3)
    CALL deleteDataLayout(div2)
    CALL deleteDataLayout(div1)

end subroutine aprioriUdivTb

end module post_mhd_helic_lib
!> @}
