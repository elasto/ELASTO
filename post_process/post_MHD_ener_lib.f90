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

module post_mhd_ener_lib

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
    public          :: eqmhd_lib
    public          :: EmhdGS_lib
    public          :: aprioriMHDEm
    public          :: aprioriEQDmMHD
    public          :: aprioBDivtub
contains
!======== EQMHD==========
!------------------------------------------------------------------------------
!> Compute and save mean values of all quantities in magnetic energie equation.
!! @author Mouloud KESSAR
!!    @param[in]    spec_rank     = mpi rank inside spectral communicator
!!    @param[in]    Uin,Vin,Win   = Velocity fields ins physical space
!!    @param[in]    Bin           = magnetic Field
!!    @param[in]    Wave          = wavenumber
!!    @param[in]    nbcpus        = numbers of cpus
!! @details
!!    
!------------------------------------------------------------------------------
subroutine eqmhd_lib(Uin,Vin,Win,Bxin,Byin,Bzin,Ukin,Vkin,Wkin,&
               & Bxink,Byink,Bzink,Wave,nbcpus,spec_rank,&
               & meanadvecMHD, meandiffvisMHD, meandissvisMHD, meanlorrMHD, meanforMHD ,somme,etain,meanEM)


 implicit none
 integer, intent(in) :: nbcpus,spec_rank
 TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Uin,Vin,Win
 TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Bxin,Byin,Bzin
 TYPE(WaveNumbers), INTENT(IN) :: Wave
 TYPE(REAL_DATA_LAYOUT) :: tab1, tab2, tab3
 TYPE(REAL_DATA_LAYOUT) :: dfdx, dfdy, dfdz
 TYPE(REAL_DATA_LAYOUT) :: advecMHD,diffvisMHD,dissvisMHD,lorrMHD,forMHD,EnMHD
 REAL(WP) :: meanadvecMHD, meandiffvisMHD, meandissvisMHD, meanlorrMHD, meanforMHD ,somme,meanEM
 REAL(WP) :: etain
 TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Ukin,Vkin,Wkin
 TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Bxink,Byink,Bzink
 TYPE(COMPLEX_DATA_LAYOUT) :: tabk1,tabk2,tabk3,Wtabk
 logical :: res


     IF ((.NOT. copyStructOnly(Uin,advecMHD)) .OR. &
        &(.NOT. copyStructOnly(Uin,diffvisMHD)) .OR. &
        &(.NOT. copyStructOnly(Uin,dissvisMHD)) .OR. &
        &(.NOT. copyStructOnly(Uin,lorrMHD)) .OR. &
        &(.NOT. copyStructOnly(Uin,forMHD)) .OR. &
        &(.NOT. copyStructOnly(Uin,EnMHD)) .OR. &
        &(.NOT. copyStructOnly(Uin,dfdx)) .OR. &
        &(.NOT. copyStructOnly(Uin,dfdy)) .OR. &
        &(.NOT. copyStructOnly(Uin,dfdz)) .OR. &
        &(.NOT. copyStructOnly(Uin,tab1)) .OR. &
        &(.NOT. copyStructOnly(Uin,tab2)) .OR. &
        &(.NOT. copyStructOnly(Uin,tab3)) ) RETURN
 

  
     IF ( (.NOT. copyStructOnly(Ukin,tabk1 ))) RETURN
     IF ( (.NOT. copyStructOnly(Ukin,tabk2 ))) RETURN
     IF ( (.NOT. copyStructOnly(Ukin,tabk3 ))) RETURN
     IF ( (.NOT. copyStructOnly(Ukin,Wtabk ))) RETURN

  !!!Advection
  EnMHD%values = 0.5_WP*(Bxin%values*Bxin%values + Byin%values*Byin%values + Bzin%values*Bzin%values)
!   call computePhysicalGradient(Wave,EnMHD,dfxdx,dfxdy,dfxdz)

           CALL ftran(EnMHD,Wtabk,res)
           call computeGradientKmain(wave,Wtabk,tabk1,tabk2,tabk3)

           CALL btran(tabk1,dfdx,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk2,dfdy,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk3,dfdz,res)
           IF (.NOT.res) RETURN

  advecMHD%values = -Uin%values*dfdx%values - Vin%values*dfdy%values - Win%values*dfdz%values
  
  
  !!!Diffusion visqueuse
!   tab1%values = dfxdx%values
!   tab2%values = dfxdy%values
!   tab3%values = dfxdz%values
!   call computePhysicalGradient(Wave,tab1,dfxdx,dfxdy,dfxdz)

!   call computePhysicalGradient(Wave,tab2,dfydx,dfydy,dfydz)
!   call computePhysicalGradient(Wave,tab3,dfzdx,dfzdy,dfzdz)
!   diffvisMHD%values = etain*( dfxdx%values + dfydy%values + dfzdz%values )
    !!!Diffusion visqueuse
    tab1%values = dfdx%values
    tab2%values = dfdy%values
    tab3%values = dfdz%values 
!     call computePhysicalGradient(Wave,tab1,dfxdx,dfxdy,dfxdz)
!     call computePhysicalGradient(Wave,tab2,dfydx,dfydy,dfydz)
!     call computePhysicalGradient(Wave,tab3,dfzdx,dfzdy,dfzdz)
           CALL ftran(tab1,Wtabk,res)
           call computeGradientKmain(wave,Wtabk,tabk1,tabk2,tabk3)

           CALL btran(tabk1,dfdx,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk2,dfdy,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk3,dfdz,res)
           IF (.NOT.res) RETURN
       diffvisMHD%values = etain*dfdx%values

           CALL ftran(tab2,Wtabk,res)
           call computeGradientKmain(wave,Wtabk,tabk1,tabk2,tabk3)

           CALL btran(tabk1,dfdx,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk2,dfdy,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk3,dfdz,res)
           IF (.NOT.res) RETURN
       diffvisMHD%values = diffvisMHD%values + etain*dfdy%values

           CALL ftran(tab3,Wtabk,res)
           call computeGradientKmain(wave,Wtabk,tabk1,tabk2,tabk3)

           CALL btran(tabk1,dfdx,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk2,dfdy,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk3,dfdz,res)
           IF (.NOT.res) RETURN

       diffvisMHD%values = diffvisMHD%values + etain*dfdz%values  
  
  !!!Dissipation visqueuse
           call computeGradientKmain(wave,Bxink,tabk1,tabk2,tabk3)
           CALL btran(tabk1,dfdx,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk2,dfdy,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk3,dfdz,res)
           IF (.NOT.res) RETURN
  dissvisMHD%values = -2.0_WP*etain*( dfdx%values*dfdx%values + dfdy%values*dfdy%values + dfdz%values*dfdz%values) 

           call computeGradientKmain(wave,Byink,tabk1,tabk2,tabk3)
           CALL btran(tabk1,dfdx,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk2,dfdy,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk3,dfdz,res)
           IF (.NOT.res) RETURN
  dissvisMHD%values =  dissvisMHD%values -2.0_WP*etain*( dfdx%values*dfdx%values + dfdy%values*dfdy%values + dfdz%values*dfdz%values) 

           call computeGradientKmain(wave,Bzink,tabk1,tabk2,tabk3)

           CALL btran(tabk1,dfdx,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk2,dfdy,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk3,dfdz,res)
           IF (.NOT.res) RETURN
  dissvisMHD%values =   dissvisMHD%values -2.0_WP*etain*( dfdx%values*dfdx%values + dfdy%values*dfdy%values + dfdz%values*dfdz%values) 
  

  !!Force de lorrentz
  tab1%values = Bxin%values*Bxin%values
  tab2%values = Bxin%values*Byin%values
  tab3%values = Bxin%values*Bzin%values

           CALL ftran(tab1,Wtabk,res)
           call computeGradientKmain(wave,Wtabk,tabk1,tabk2,tabk3)
           CALL btran(tabk1,dfdx,res)
           IF (.NOT.res) RETURN

           CALL ftran(tab2,Wtabk,res)
           call computeGradientKmain(wave,Wtabk,tabk1,tabk2,tabk3)
           CALL btran(tabk2,dfdy,res)
           IF (.NOT.res) RETURN

           CALL ftran(tab3,Wtabk,res)
           call computeGradientKmain(wave,Wtabk,tabk1,tabk2,tabk3)
           CALL btran(tabk3,dfdz,res)
           IF (.NOT.res) RETURN




!   call computePhysicalGradient(Wave,tab1,dfxdx,dfxdy,dfxdz)
!   call computePhysicalGradient(Wave,tab2,dfydx,dfydy,dfydz)
!   call computePhysicalGradient(Wave,tab3,dfzdx,dfzdy,dfzdz)

  lorrMHD%values = -2.*Uin%values*( dfdx%values + dfdy%values + dfdz%values )
  tab1%values = Byin%values*Bxin%values
  tab2%values = Byin%values*Byin%values
  tab3%values = Byin%values*Bzin%values

           CALL ftran(tab1,Wtabk,res)
           call computeGradientKmain(wave,Wtabk,tabk1,tabk2,tabk3)
           CALL btran(tabk1,dfdx,res)
           IF (.NOT.res) RETURN

           CALL ftran(tab2,Wtabk,res)
           call computeGradientKmain(wave,Wtabk,tabk1,tabk2,tabk3)
           CALL btran(tabk2,dfdy,res)
           IF (.NOT.res) RETURN

           CALL ftran(tab3,Wtabk,res)
           call computeGradientKmain(wave,Wtabk,tabk1,tabk2,tabk3)
           CALL btran(tabk3,dfdz,res)
           IF (.NOT.res) RETURN

 
!   call computePhysicalGradient(Wave,tab1,dfxdx,dfxdy,dfxdz)
!   call computePhysicalGradient(Wave,tab2,dfydx,dfydy,dfydz)
!   call computePhysicalGradient(Wave,tab3,dfzdx,dfzdy,dfzdz)

  lorrMHD%values = lorrMHD%values - 2.*Vin%values*( dfdx%values + dfdy%values + dfdz%values )

  tab1%values = Bzin%values*Bxin%values
  tab2%values = Bzin%values*Byin%values
  tab3%values = Bzin%values*Bzin%values
           CALL ftran(tab1,Wtabk,res)
           call computeGradientKmain(wave,Wtabk,tabk1,tabk2,tabk3)
           CALL btran(tabk1,dfdx,res)
           IF (.NOT.res) RETURN

           CALL ftran(tab2,Wtabk,res)
           call computeGradientKmain(wave,Wtabk,tabk1,tabk2,tabk3)
           CALL btran(tabk2,dfdy,res)
           IF (.NOT.res) RETURN

           CALL ftran(tab3,Wtabk,res)
           call computeGradientKmain(wave,Wtabk,tabk1,tabk2,tabk3)
           CALL btran(tabk3,dfdz,res)
           IF (.NOT.res) RETURN
! 
!   call computePhysicalGradient(Wave,tab1,dfxdx,dfxdy,dfxdz)
!   call computePhysicalGradient(Wave,tab2,dfydx,dfydy,dfydz)
!   call computePhysicalGradient(Wave,tab3,dfzdx,dfzdy,dfzdz)

  lorrMHD%values = lorrMHD%values - 2.0_WP*Win%values*( dfdx%values + dfdy%values + dfdz%values )
 
  
  !!!Force
  tab1%values = Uin%values*Bxin%values*Bxin%values
  tab2%values = Uin%values*Bxin%values*Byin%values
  tab3%values = Uin%values*Bxin%values*Bzin%values

!   call computePhysicalGradient(Wave,tab1,dfxdx,dfxdy,dfxdz)
!   call computePhysicalGradient(Wave,tab2,dfydx,dfydy,dfydz)
!   call computePhysicalGradient(Wave,tab3,dfzdx,dfzdy,dfzdz)
           CALL ftran(tab1,Wtabk,res)
           call computeGradientKmain(wave,Wtabk,tabk1,tabk2,tabk3)
           CALL btran(tabk1,dfdx,res)
           IF (.NOT.res) RETURN

           CALL ftran(tab2,Wtabk,res)
           call computeGradientKmain(wave,Wtabk,tabk1,tabk2,tabk3)
           CALL btran(tabk2,dfdy,res)
           IF (.NOT.res) RETURN

           CALL ftran(tab3,Wtabk,res)
           call computeGradientKmain(wave,Wtabk,tabk1,tabk2,tabk3)
           CALL btran(tabk3,dfdz,res)
           IF (.NOT.res) RETURN

  forMHD%values = 2.0_WP*( dfdx%values + dfdy%values + dfdz%values )

  tab1%values = Vin%values*Byin%values*Bxin%values
  tab2%values = Vin%values*Byin%values*Byin%values
  tab3%values = Vin%values*Byin%values*Bzin%values
           CALL ftran(tab1,Wtabk,res)
           call computeGradientKmain(wave,Wtabk,tabk1,tabk2,tabk3)
           CALL btran(tabk1,dfdx,res)
           IF (.NOT.res) RETURN

           CALL ftran(tab2,Wtabk,res)
           call computeGradientKmain(wave,Wtabk,tabk1,tabk2,tabk3)
           CALL btran(tabk2,dfdy,res)
           IF (.NOT.res) RETURN

           CALL ftran(tab3,Wtabk,res)
           call computeGradientKmain(wave,Wtabk,tabk1,tabk2,tabk3)
           CALL btran(tabk3,dfdz,res)
           IF (.NOT.res) RETURN
! 
!   call computePhysicalGradient(Wave,tab1,dfxdx,dfxdy,dfxdz)
!   call computePhysicalGradient(Wave,tab2,dfydx,dfydy,dfydz)
!   call computePhysicalGradient(Wave,tab3,dfzdx,dfzdy,dfzdz)

  forMHD%values = forMHD%values + 2.0_WP*( dfdx%values + dfdy%values + dfdz%values )

  tab1%values = Win%values*Bzin%values*Bxin%values
  tab2%values = Win%values*Bzin%values*Byin%values
  tab3%values = Win%values*Bzin%values*Bzin%values
           CALL ftran(tab1,Wtabk,res)
           call computeGradientKmain(wave,Wtabk,tabk1,tabk2,tabk3)
           CALL btran(tabk1,dfdx,res)
           IF (.NOT.res) RETURN

           CALL ftran(tab2,Wtabk,res)
           call computeGradientKmain(wave,Wtabk,tabk1,tabk2,tabk3)
           CALL btran(tabk2,dfdy,res)
           IF (.NOT.res) RETURN

           CALL ftran(tab3,Wtabk,res)
           call computeGradientKmain(wave,Wtabk,tabk1,tabk2,tabk3)
           CALL btran(tabk3,dfdz,res)
           IF (.NOT.res) RETURN
!   call computePhysicalGradient(Wave,tab1,dfxdx,dfxdy,dfxdz)
!   call computePhysicalGradient(Wave,tab2,dfydx,dfydy,dfydz)
!   call computePhysicalGradient(Wave,tab3,dfzdx,dfzdy,dfzdz)
  forMHD%values = forMHD%values + 2.0_WP*(dfdx%values + dfdy%values + dfdz%values )
 
 
 
!   if (rank .eq. 0) print*,'Mean calculation'

    meanadvecMHD   = computeFieldAvg(advecMHD,spec_rank,nbcpus) 
    meandiffvisMHD = computeFieldAvg(diffvisMHD,spec_rank,nbcpus) 
    meandissvisMHD = computeFieldAvg(dissvisMHD,spec_rank,nbcpus) 
    meanlorrMHD    = computeFieldAvg(lorrMHD,spec_rank,nbcpus)
    meanforMHD     = computeFieldAvg(forMHD,spec_rank,nbcpus)
    meanEM         = computeFieldAvg(EnMHD,spec_rank,nbcpus)
    somme = meanadvecMHD + meandiffvisMHD + meandissvisMHD + meanlorrMHD + meanforMHD


  CALL deleteDataLayout(advecMHD)
  CALL deleteDataLayout(diffvisMHD)
  CALL deleteDataLayout(dissvisMHD)
  CALL deleteDataLayout(lorrMHD)
  CALL deleteDataLayout(forMHD)
  CALL deleteDataLayout(EnMHD)
  CALL deleteDataLayout(dfdx)
  CALL deleteDataLayout(dfdy)
  CALL deleteDataLayout(dfdz)
  CALL deleteDataLayout(tab1)
  CALL deleteDataLayout(tab2)
  CALL deleteDataLayout(tab3)

end subroutine eqmhd_lib


!======= EmhdGS =========
!------------------------------------------------------------------------------
!> Compute and save mean values of all quantities in
!! @author Mouloud KESSAR
!!    @param[in]    spec_rank     = mpi rank inside spectral communicator
!!    @param[in]    spec_rank     = mpi rank inside spectral communicator
!!    @param[in]    Uin,Vin,Win   = Velocity fields ins physical space
!!    @param[in]    Bin           = magnetic Field
!!    @param[in]    Wave          = wavenumber
!!    @param[in]    nbcpus        = numbers of cpus
!!    @param[in]    spec_rank     = mpi rank inside spectral communicator
!!    @param[in]    Ufin,Vfin,Wfin   = filtered Velocity fields ins physical space
!!    @param[in]    Bfxin, Bfyin, Bfzin    = filtered magnetic Field
!!    @param[in]    delta     = filtersize
!!    @param[in]    ndelta    = filter number
!!    @param[in]    filter    = filter type
!! @details
!!
!------------------------------------------------------------------------------

subroutine EmhdGS_lib( Ufin,Vfin,Wfin,Bfxin,Bfyin,Bfzin,&
                 & Ukin,Vkin,Wkin, &
                 & Uin,Vin,Win,Bxin,Byin,Bzin,Bfxkin,Bfykin,Bfzkin,&
                 & Wave,delta,ndelta,filter,nbcpus,spec_rank,etain,&
                 & meanDiffUBSGSEMHD,meanDissUBSGSEMHD,meanadvecGSEMHD,&
                 & meanDifGSEMHD, meanDissGSEMHD, TransGSEmEk, meanProdEMHD ,sommeGS,meanEmGS)

  IMPLICIT NONE
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Bfxin,Bfyin,Bfzin,Ufin,Vfin,Wfin
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Uin,Vin,Win
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Bxin,Byin,Bzin
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Bfxkin,Bfykin,Bfzkin
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) ::Ukin,Vkin,Wkin
!   TYPE(REAL_DATA_LAYOUT) :: DissGSEMHD, AdvecGSEMHD, DifGSEMHD
  TYPE(REAL_DATA_LAYOUT) :: DiffUBSGSEMHD, DissUBSGSEMHD!,ProdEMhd, DissSGSEMHD
  REAL(WP) :: delta,etain
  REAL(WP) :: meanDissGSEMHD, meanAdvecGSEMHD, meanDifGSEMHD,meanDiffUBSGSEMHD
  REAL(WP) :: meanDissUBSGSEMHD,meanProdEMHD , meanDissSGSEMHD 
  REAL(WP) :: TransGSEmEk,sommeGS,meanEmGS
  INTEGER, INTENT(IN) :: nbcpus, spec_rank, filter,ndelta
  TYPE(REAL_DATA_LAYOUT) :: tab1, tab2, tab3
  TYPE(COMPLEX_DATA_LAYOUT) :: tabk1, tabk2, tabk3
  TYPE(COMPLEX_DATA_LAYOUT) :: Wtabk
  type(waveNumbers), intent(in)        :: wave
  TYPE(REAL_DATA_LAYOUT) :: T12,T13,T23
  TYPE(REAL_DATA_LAYOUT) :: J12,J13,J23
  TYPE(REAL_DATA_LAYOUT) :: dfdx, dfdy, dfdz
  logical :: res


  
  call eqmhd_lib(Ufin,Vfin,Wfin,Bfxin,Bfyin,Bfzin,Ukin,Vkin,Wkin,&
               & Bfxkin,Bfykin,Bfzkin,Wave,nbcpus,spec_rank,&
               & meanadvecGSEMHD, meanDifGSEMHD, meanDissGSEMHD, TransGSEmEk, meanProdEMHD ,sommeGS,etain,meanEmGS)

  IF ((.NOT. copyStructOnly(Uin,DiffUBSGSEMHD)) .OR. &
     &(.NOT. copyStructOnly(Uin,DissUBSGSEMHD)) ) RETURN

     IF ((.NOT. copyStructOnly(Uin,J12)) .OR. &
        &(.NOT. copyStructOnly(Uin,J13)) .OR. &
        &(.NOT. copyStructOnly(Uin,J23))  ) RETURN
     IF ((.NOT. copyStructOnly(Uin,T12)) .OR. &
        &(.NOT. copyStructOnly(Uin,T13)) .OR. &
        &(.NOT. copyStructOnly(Uin,T23))  ) RETURN


  call computeJijTensor(Bfxkin,Bfykin,Bfzkin,wave,J12,J13,J23, &
             & res)

  call computeT_ijVecAVec_B(Ukin,Uin,Vin,Win,Ufin,Vfin,Wfin,&
                          & Bxin,Byin,Bzin,Bfxin,Bfyin,Bfzin,&
                          & wave,T12,T13,T23,filter,res,delta)

     IF ( (.NOT. copyStructOnly(Ukin,tabk1 ))) RETURN
     IF ( (.NOT. copyStructOnly(Ukin,tabk2 ))) RETURN
     IF ( (.NOT. copyStructOnly(Ukin,tabk3 ))) RETURN
     IF ( (.NOT. copyStructOnly(Ukin,Wtabk ))) RETURN

!      IF ( (.NOT. copyStructOnly(Uin,tab1 ))) RETURN
     IF ( (.NOT. copyStructOnly(Uin,tab2 ))) RETURN
     IF ( (.NOT. copyStructOnly(Uin,tab3)) .OR. &
        & (.NOT. copyStructOnly(Uin,dfdx)) .OR. &
        & (.NOT. copyStructOnly(Uin,dfdy)) .OR. &
        & (.NOT. copyStructOnly(Uin,dfdz)) .OR. &
        & (.NOT. copyStructOnly(Uin,tab1)) ) RETURN

     tab1%values = 2.0_WP*( Byin%values*T12%values + Bzin%values*T13%values )
     tab2%values = 2.0_WP*(-Bxin%values*T12%values + Bzin%values*T23%values )
     tab3%values = 2.0_WP*(-Bxin%values*T13%values - Byin%values*T23%values )

           CALL ftran(tab1,Wtabk,res)
           call computeGradientKmain(wave,Wtabk,tabk1,tabk2,tabk3)
           CALL btran(tabk1,dfdx,res)
           IF (.NOT.res) RETURN
       
           CALL ftran(tab2,Wtabk,res)
           call computeGradientKmain(wave,Wtabk,tabk1,tabk2,tabk3)
           CALL btran(tabk2,dfdy,res)
           IF (.NOT.res) RETURN

           CALL ftran(tab3,Wtabk,res)
           call computeGradientKmain(wave,Wtabk,tabk1,tabk2,tabk3)
           CALL btran(tabk3,dfdz,res)
           IF (.NOT.res) RETURN


  DiffUBSGSEMHD%values = -dfdx%values - dfdy%values - dfdz%values
  DissUBSGSEMHD%values = 4.0_WP*(T13%values*J13%values + T12%values*J12%values + T23%values*J23%values)

  meanDiffUBSGSEMHD = computeFieldAvg(DiffUBSGSEMHD, spec_rank,nbcpus)
  meanDissUBSGSEMHD = computeFieldAvg(DissUBSGSEMHD, spec_rank,nbcpus)

  CALL deleteDataLayout(DiffUBSGSEMHD) 
  CALL deleteDataLayout(DissUBSGSEMHD)

  CALL deleteDataLayout(J12)
  CALL deleteDataLayout(J13)
  CALL deleteDataLayout(J23)

  CALL deleteDataLayout(T12)
  CALL deleteDataLayout(T13)
  CALL deleteDataLayout(T23)

  CALL deleteDataLayout(tab1)
  CALL deleteDataLayout(tab2)
  CALL deleteDataLayout(tab3)

  CALL deleteDataLayout(Wtabk)
  CALL deleteDataLayout(tabk3)
  CALL deleteDataLayout(tabk2)
  CALL deleteDataLayout(tabk1)
  CALL deleteDataLayout(dfdx)
  CALL deleteDataLayout(dfdy)
  CALL deleteDataLayout(dfdz)

end subroutine EmhdGS_lib

subroutine aprioriMHDEm(Uin,Vin,Win,Bxin,Byin,Bzin,Ukin,Vkin,Wkin,&
                    & Bxkin,Bykin,Bzkin,Wave,spec_rank,delta,filter,meanEmTBbmodel,nbcpus)

  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Uin,Vin,Win
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Bxin,Byin,Bzin 
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Ukin,Vkin,Wkin
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Bxkin,Bykin,Bzkin
  TYPE(WaveNumbers), INTENT(IN) :: Wave

  INTEGER :: spec_rank,nbcpus,filter
  TYPE(REAL_DATA_LAYOUT) :: T12,T13,T23
  TYPE(REAL_DATA_LAYOUT) :: J12,J13,J23
  TYPE(REAL_DATA_LAYOUT) :: EmTBbmodel
  logical :: res
  REAL(WP) :: delta,meanEmTBbmodel

  IF ((.NOT. copyStructOnly(Bxin,J12)) .OR. &
     &(.NOT. copyStructOnly(Bxin,J13)) .OR. &
     &(.NOT. copyStructOnly(Bxin,J23)) .OR. &
     &(.NOT. copyStructOnly(Bxin,EmTBbmodel)) .OR. &
     &(.NOT. copyStructOnly(Bxin,t12)) .OR. &
     &(.NOT. copyStructOnly(Bxin,t13)) .OR. &
     &(.NOT. copyStructOnly(Bxin,t23))    ) RETURN 



     call smagorinskyMHDTij(T12,T13,T23,Bxkin,Bykin,Bzkin,&
                           &wave,spec_rank,res, delta)
  
     call computeJijTensor(BxKin,ByKin,BzKin,wave,J12,J13,J23, res)
!!!!! facteur 4, car on a un deux dans l'expression et le tenseur est antisymétrique
     EmTBbmodel%values = 4.0_WP*(J12%values*T12%values &
                               &+J13%values*T13%values &
                               &+J23%values*T23%values )
    meanEmTBbmodel  = computeFieldAvg(EmTBbmodel, spec_rank,nbcpus)

    CALL deleteDataLayout(J12)
    CALL deleteDataLayout(J13)
    CALL deleteDataLayout(J23)
    CALL deleteDataLayout(T12)
    CALL deleteDataLayout(T13)
    CALL deleteDataLayout(T23)
    CALL deleteDataLayout(EmTBbmodel)



end subroutine aprioriMHDEm

subroutine aprioriEQDmMHD(Uin,Vin,Win,Bxin,Byin,Bzin,Ukin,Vkin,Wkin,&
                    & Bxkin,Bykin,Bzkin,Wave,spec_rank,filter,delta,meanEqdmTUuuSmag,&
                    & meanEqdmTUubSmag,meanEqdmTUuuSimil,meanEqdmTUubSimil,meanEqdmTUuuGrad, meanEqdmTUubGrad,nbcpus)

  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Uin,Vin,Win
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Bxin,Byin,Bzin 
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Ukin,Vkin,Wkin
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Bxkin,Bykin,Bzkin
  TYPE(WaveNumbers), INTENT(IN) :: Wave
  
  TYPE(COMPLEX_DATA_LAYOUT) :: div1,div2,div3
  INTEGER :: spec_rank,nbcpus,numVel,type_avg,filter
  TYPE(REAL_DATA_LAYOUT) :: T12,T13,T23
  TYPE(REAL_DATA_LAYOUT) :: T11,T22,T33
  TYPE(REAL_DATA_LAYOUT) :: S12,S13,S23
  TYPE(REAL_DATA_LAYOUT) :: S11,S22,S33
  TYPE(REAL_DATA_LAYOUT) :: Wtab
  logical :: res
  REAL(WP) :: delta,meanEqdmTUuuSmag, meanEqdmTUubSmag,meanEqdmTUuuSimil,meanEqdmTUubSimil
  REAL(WP) :: meanEqdmTUuuGrad, meanEqdmTUubGrad

  IF ((.NOT. copyStructOnly(Bxin,S12)) .OR. &
     &(.NOT. copyStructOnly(Bxin,S13)) .OR. &
     &(.NOT. copyStructOnly(Bxin,S23)) .OR. &
     &(.NOT. copyStructOnly(Bxin,S11)) .OR. &
     &(.NOT. copyStructOnly(Bxin,S22)) .OR. &
     &(.NOT. copyStructOnly(Bxin,S33)) .OR. &
     &(.NOT. copyStructOnly(Bxin,Wtab)) .OR. &
     &(.NOT. copyStructOnly(Bxin,T12)) .OR. &
     &(.NOT. copyStructOnly(Bxin,T13)) .OR. &
     &(.NOT. copyStructOnly(Bxin,T23))  .OR. &
     &(.NOT. copyStructOnly(Bxin,T11)) .OR. &
     &(.NOT. copyStructOnly(Bxin,T22)) .OR. &
     &(.NOT. copyStructOnly(Bxin,T33))   ) RETURN 

  IF ((.NOT. copyStructOnly(Bxkin,div1)) .OR. &
     &(.NOT. copyStructOnly(Bxkin,div2)) .OR. &
     &(.NOT. copyStructOnly(Bxkin,div3))   ) RETURN 

     numVel=1
     type_avg=0
     call computeStressTensor(Ukin,VKin,WKin,wave,S11,S12,S13,S22,S23,S33,res)

     call smagorinskyVel(div1,div2,div3,Uin,Vin,Win,Ukin,Vkin,Wkin,wave,numVel,spec_rank,res, &
                          &  filter,type_avg,T11,T12,T13,T22,T23,T33,delta)


     Wtab%values = 4.0_WP*(   S12%values*T12%values &
                          & + S13%values*T13%values &
                          & + S23%values*T23%values) &
               & + 2.0_WP*(   S11%values*T11%values &
                          & + S22%values*T22%values &
                          & + S23%values*T33%values )

    meanEqdmTUuuSmag  = computeFieldAvg(Wtab, spec_rank,nbcpus)


     call smagorinskyVel(div1,div2,div3,Bxin,Byin,Bzin,Bxkin,Bykin,Bzkin,wave,numVel,spec_rank,res, &
                          &  filter,type_avg,T11,T12,T13,T22,T23,T33,delta)


     Wtab%values = - 4.0_WP*(  S12%values*T12%values &
                           & + S13%values*T13%values &
                           & + S23%values*T23%values ) &
                 & - 2.0_WP*(  S11%values*T11%values &
                           & + S22%values*T22%values &
                           & + S23%values*T33%values )

    meanEqdmTUubSmag = computeFieldAvg(Wtab, spec_rank,nbcpus)

    !!!!!!similarity

    call SimilVel(t11,t12,t13,t22,t23,t33,Uin,Vin,Win,Ukin,Vkin,Wkin,Wave,delta,filter,spec_rank,nbcpus)

     Wtab%values = 4.0_WP*(   S12%values*T12%values &
                          & + S13%values*T13%values &
                          & + S23%values*T23%values) &
               & + 2.0_WP*(   S11%values*T11%values &
                          & + S22%values*T22%values &
                          & + S23%values*T33%values )

    meanEqdmTUuuSimil  = computeFieldAvg(Wtab, spec_rank,nbcpus)

    call SimilVel(t11,t12,t13,t22,t23,t33,Bxin,Byin,Bzin,Bxkin,Bykin,Bzkin,Wave,delta,filter,spec_rank,nbcpus)

     Wtab%values = - 4.0_WP*(   S12%values*T12%values &
                            & + S13%values*T13%values &
                            & + S23%values*T23%values) &
                 & - 2.0_WP*(   S11%values*T11%values &
                            & + S22%values*T22%values &
                            & + S23%values*T33%values )

    meanEqdmTUubSimil  = computeFieldAvg(Wtab, spec_rank,nbcpus)

    !!!!!!Gradient

    call gradientVel(t11,t12,t13,t22,t23,t33,Uin,Vin,Win,Wave,delta,spec_rank,nbcpus)
     Wtab%values = 4.0_WP*(   S12%values*T12%values &
                          & + S13%values*T13%values &
                          & + S23%values*T23%values) &
               & + 2.0_WP*(   S11%values*T11%values &
                          & + S22%values*T22%values &
                          & + S23%values*T33%values )

    meanEqdmTUuuGrad  = computeFieldAvg(Wtab, spec_rank,nbcpus)

    call gradientVel(t11,t12,t13,t22,t23,t33,Bxin,Byin,Bzin,Wave,delta,spec_rank,nbcpus)
     Wtab%values = - 4.0_WP*(   S12%values*T12%values &
                            & + S13%values*T13%values &
                            & + S23%values*T23%values) &
                 & - 2.0_WP*(   S11%values*T11%values &
                            & + S22%values*T22%values &
                            & + S23%values*T33%values )

    meanEqdmTUubGrad  = computeFieldAvg(Wtab, spec_rank,nbcpus)

    CALL deleteDataLayout(S11)
    CALL deleteDataLayout(S22)
    CALL deleteDataLayout(S33)
    CALL deleteDataLayout(S12)
    CALL deleteDataLayout(S13)
    CALL deleteDataLayout(S23)

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

end subroutine aprioriEQDmMHD

subroutine aprioBDivtub(Uin,Vin,Win,Bxin,Byin,Bzin,Ukin,Vkin,Wkin,&
                    & Bxkin,Bykin,Bzkin,Ufin,Vfin,Wfin,Bfxin,Bfyin,Bfzin,& 
                    & Wave,spec_rank,delta,filter,nbcpus,meanbdivTubexact, meanbdivTubSmag)

  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Uin,Vin,Win
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Bxin,Byin,Bzin 
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Ufin,Vfin,Wfin,Bfxin,Bfyin,Bfzin
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Ukin,Vkin,Wkin
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Bxkin,Bykin,Bzkin
  TYPE(WaveNumbers), INTENT(IN) :: Wave
  
  TYPE(COMPLEX_DATA_LAYOUT) :: div1,div2,div3
  TYPE(REAL_DATA_LAYOUT) :: div1R,div2R,div3R

  INTEGER :: spec_rank,nbcpus,filter
  TYPE(REAL_DATA_LAYOUT) :: T12,T13,T23
  TYPE(REAL_DATA_LAYOUT) :: Wtab
  logical :: res
  REAL(WP) :: delta
  REAL(WP) :: meanbdivTubexact, meanbdivTubSmag

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

     Wtab%values = div1R%values*Bfxin%values + div2R%values*Bfyin%values + div3R%values*Bfzin%values
     meanbdivTubexact  = computeFieldAvg(Wtab, spec_rank,nbcpus)

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

     Wtab%values = div1R%values*Bfxin%values + div2R%values*Bfyin%values + div3R%values*Bfzin%values
     meanbdivTubSmag  = computeFieldAvg(Wtab, spec_rank,nbcpus)

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

end subroutine aprioBDivtub

end module post_mhd_ener_lib
!> @}
