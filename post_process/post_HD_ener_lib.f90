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

module post_hd_ener_lib

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
    use physical_values_tools
    implicit none
    private

     ! ==== Public procedures ==== 
    public          :: EQDMGSHD
    public          :: eqdm
    public          :: aprioriEQDHD
    ! ======= Mean terms ========
 contains
!======== EQDM==========
!------------------------------------------------------------------------------
!> Compute and save mean values of all quantities in kinetic energie equation.
!! @author Mouloud KESSAR
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in]    Uin,Vin,Win     = Velocity fields ins physical space
!!    @param[in]    Bin             = magnetic Field
!!    @param[in]    Wave            = wavenumber
!! @details
!!    
!------------------------------------------------------------------------------

subroutine eqdm(Uin,Vin,Win,Ukin,VKin,WKin,Wave,nbcpus,spec_rank,mu,meanadvecQDM, meandiffvisQDM, meandissvisQDM, somme, meanEk)


   implicit none

   integer, intent(in)    :: nbcpus,spec_rank
   TYPE(REAL_DATA_LAYOUT), INTENT(IN)          :: Uin,Vin,Win
   TYPE(WaveNumbers), INTENT(IN) :: Wave
   TYPE(REAL_DATA_LAYOUT) :: tab1, tab2, tab3
   TYPE(REAL_DATA_LAYOUT) :: dfdx, dfdy, dfdz
   TYPE(REAL_DATA_LAYOUT) :: advecQDM,diffvisQDM,dissvisQDM,EnQDM
   REAL(WP) :: meanadvecQDM, meandiffvisQDM, meandissvisQDM, somme, meanEk
   REAL(WP) :: mu
   TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Ukin,Vkin,Wkin
   logical :: res
   TYPE(COMPLEX_DATA_LAYOUT) :: tabk1,tabk2,tabk3,Wtabk

     IF ((.NOT. copyStructOnly(Uin,advecQDM)) .OR. &
        &(.NOT. copyStructOnly(Uin,diffvisQDM)) .OR. &
        &(.NOT. copyStructOnly(Uin,dissvisQDM)) .OR. &
        &(.NOT. copyStructOnly(Uin,EnQDM)) .OR. &
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

     !!Energie!!
    EnQDM%values = 0.5_WP*( Uin%values*Uin%values + Vin%values*Vin%values + Win%values*Win%values)

    !!!Advection
!     call computePhysicalGradient(Wave,EnQDM,dfxdx,dfxdy,dfxdz)
           CALL ftran(EnQDM,Wtabk,res)
           call computeGradientKmain(wave,Wtabk,tabk1,tabk2,tabk3)

           CALL btran(tabk1,dfdx,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk2,dfdy,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk3,dfdz,res)
           IF (.NOT.res) RETURN
    
    advecQDM%values = -Uin%values*dfdx%values -Vin%values*dfdy%values -Win%values*dfdz%values
 
  
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
       diffvisQDM%values = mu*dfdx%values

           CALL ftran(tab2,Wtabk,res)
           call computeGradientKmain(wave,Wtabk,tabk1,tabk2,tabk3)

           CALL btran(tabk1,dfdx,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk2,dfdy,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk3,dfdz,res)
           IF (.NOT.res) RETURN
       diffvisQDM%values = diffvisQDM%values + mu*dfdy%values

           CALL ftran(tab3,Wtabk,res)
           call computeGradientKmain(wave,Wtabk,tabk1,tabk2,tabk3)

           CALL btran(tabk1,dfdx,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk2,dfdy,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk3,dfdz,res)
           IF (.NOT.res) RETURN

       diffvisQDM%values = diffvisQDM%values + mu*dfdz%values
  
    !!!Dissipation visqueuse
!     call computePhysicalGradient(Wave,Uin,dfxdx,dfxdy,dfxdz)
!     call computePhysicalGradient(Wave,Vin,dfydx,dfydy,dfydz)
!     call computePhysicalGradient(Wave,Win,dfzdx,dfzdy,dfzdz)
           call computeGradientKmain(wave,Ukin,tabk1,tabk2,tabk3)

           CALL btran(tabk1,dfdx,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk2,dfdy,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk3,dfdz,res)
           IF (.NOT.res) RETURN
dissvisQDM%values = -2*mu*(dfdx%values*dfdx%values + dfdy%values*dfdy%values + dfdz%values*dfdz%values )
  
           call computeGradientKmain(wave,Vkin,tabk1,tabk2,tabk3)
       
           CALL btran(tabk1,dfdx,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk2,dfdy,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk3,dfdz,res)
           IF (.NOT.res) RETURN
dissvisQDM%values = dissvisQDM%values -2*mu*(dfdx%values*dfdx%values + dfdy%values*dfdy%values + dfdz%values*dfdz%values )

           call computeGradientKmain(wave,Wkin,tabk1,tabk2,tabk3)

           CALL btran(tabk1,dfdx,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk2,dfdy,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk3,dfdz,res)
           IF (.NOT.res) RETURN
dissvisQDM%values = dissvisQDM%values -2*mu*(dfdx%values*dfdx%values + dfdy%values*dfdy%values + dfdz%values*dfdz%values )

 
  
!     if (rank .eq. 0) print*,'Mean calculation'
    meanadvecQDM   = computeFieldAvg(advecQDM,spec_rank,nbcpus) 
    meandiffvisQDM = computeFieldAvg(diffvisQDM,spec_rank,nbcpus) 
    meandissvisQDM = computeFieldAvg(dissvisQDM,spec_rank,nbcpus) 
    meanEk = computeFieldAvg(EnQDM,spec_rank,nbcpus) 

  CALL deleteDataLayout(advecQDM)
  CALL deleteDataLayout(diffvisQDM)
  CALL deleteDataLayout(dissvisQDM)
  CALL deleteDataLayout(EnQDM)
  CALL deleteDataLayout(dfdx)
  CALL deleteDataLayout(dfdy)
  CALL deleteDataLayout(dfdz)
  CALL deleteDataLayout(tab1)
  CALL deleteDataLayout(tab2)
  CALL deleteDataLayout(tab3)

  CALL deleteDataLayout(tabk1)
  CALL deleteDataLayout(tabk2)
  CALL deleteDataLayout(tabk3)

  CALL deleteDataLayout(Wtabk)

end subroutine eqdm
 


!======= EQDMGS =========
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


subroutine EQDMGSHD(Ufin,Vfin,Wfin,Uin,Vin,Win,Ufkin,Vfkin,Wfkin,Ukin,VKin,WKin,Wave,delta,ndelta,filter,nbcpus,spec_rank,mu,&
                   & meanDissGSEQDM, meanAdvecGSEQDM, meanDifGSEQDM,&
                   & meanDiffUSGSEqdm, meanDissUSGSEqdm, meanDissSGSEQDM, meanEkGS)

  IMPLICIT NONE
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Ufin,Vfin,Wfin
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Uin,Vin,Win
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Ukin,Vkin,Wkin
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Ufkin,Vfkin,Wfkin
  TYPE(REAL_DATA_LAYOUT) :: DissGSEQDM, AdvecGSEQDM, DifGSEQDM, EkGS
  TYPE(REAL_DATA_LAYOUT) :: DiffUSGSEqdm, DissUSGSEqdm, DissSGSEQDM
  REAL(WP) :: delta
  REAL(WP) :: meanDissGSEQDM, meanAdvecGSEQDM, meanDifGSEQDM
  REAL(WP) :: meanDiffUSGSEqdm, meanDissUSGSEqdm, meanDissSGSEQDM 
  REAL(WP) :: meanEkGS
  INTEGER, INTENT(IN) :: nbcpus, spec_rank, filter,ndelta
  TYPE(REAL_DATA_LAYOUT) :: dfdx, dfdy, dfdz
!   TYPE(REAL_DATA_LAYOUT) :: dfydx, dfydy, dfydz
!   TYPE(REAL_DATA_LAYOUT) :: dfzdx, dfzdy, dfzdz
  TYPE(REAL_DATA_LAYOUT) :: tab1, tab2, tab3
  type(waveNumbers), intent(in)        :: wave
  TYPE(REAL_DATA_LAYOUT) :: T12,T13,T23
  TYPE(REAL_DATA_LAYOUT) :: T11,T22,T33
  TYPE(REAL_DATA_LAYOUT) :: S12,S13,S23
  TYPE(REAL_DATA_LAYOUT) :: S11,S22,S33
  REAL(WP) :: mu
  logical :: res
  TYPE(COMPLEX_DATA_LAYOUT) :: tabk1,tabk2,tabk3,Wtabk

     IF ((.NOT. copyStructOnly(Uin,dfdx)) .OR. &
        &(.NOT. copyStructOnly(Uin,dfdy)) .OR. &
        &(.NOT. copyStructOnly(Uin,dfdz)) .OR. &
        &(.NOT. copyStructOnly(Uin,tab1)) .OR. &
        &(.NOT. copyStructOnly(Uin,tab2)) .OR. &
        &(.NOT. copyStructOnly(Uin,tab3)) ) RETURN

     IF ( (.NOT. copyStructOnly(Uin,DissGSEQDM ))  .OR. &
         &(.NOT. copyStructOnly(Uin,DissSGSEQDM )) .OR. &
         &(.NOT. copyStructOnly(Uin,AdvecGSEQDM))  .OR. &
         &(.NOT. copyStructOnly(Uin,EkGS))  .OR. &
         &(.NOT. copyStructOnly(Uin,DifGSEQDM)) ) RETURN

     IF ( (.NOT. copyStructOnly(Ukin,tabk1 ))) RETURN
     IF ( (.NOT. copyStructOnly(Ukin,tabk2 ))) RETURN
     IF ( (.NOT. copyStructOnly(Ukin,tabk3 ))) RETURN
     IF ( (.NOT. copyStructOnly(Ukin,Wtabk ))) RETURN


  EkGS%values = 0.5_WP*( Ufin%values**2 + Vfin%values**2 + Wfin%values**2  )

           call computeGradientKmain(wave,Ufkin,tabk1,tabk2,tabk3)

           CALL btran(tabk1,dfdx,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk2,dfdy,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk3,dfdz,res)
           IF (.NOT.res) RETURN
  DissGSEQDM%values = -2.0_WP*mu*( dfdx%values**2 + dfdy%values**2 + dfdz%values**2 )

           call computeGradientKmain(wave,Vfkin,tabk1,tabk2,tabk3)
       
           CALL btran(tabk1,dfdx,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk2,dfdy,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk3,dfdz,res)
           IF (.NOT.res) RETURN
  DissGSEQDM%values = DissGSEQDM%values -2.0_WP*mu*( dfdx%values**2 + dfdy%values**2 + dfdz%values**2 )

           call computeGradientKmain(wave,Wfkin,tabk1,tabk2,tabk3)

           CALL btran(tabk1,dfdx,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk2,dfdy,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk3,dfdz,res)
           IF (.NOT.res) RETURN


  DissGSEQDM%values = DissGSEQDM%values -2.0_WP*mu*( dfdx%values**2 + dfdy%values**2 + dfdz%values**2 )


 

           call computeGradientKmain(wave,Ukin,tabk1,tabk2,tabk3)

           CALL btran(tabk1,dfdx,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk2,dfdy,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk3,dfdz,res)
           IF (.NOT.res) RETURN
  tab1%values  = -2.0_WP*mu*( dfdx%values**2 + dfdy%values**2 + dfdz%values**2 )

           call computeGradientKmain(wave,Vkin,tabk1,tabk2,tabk3)
       
           CALL btran(tabk1,dfdx,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk2,dfdy,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk3,dfdz,res)
           IF (.NOT.res) RETURN
  tab1%values  = tab1%values -2.0_WP*mu*( dfdx%values**2 + dfdy%values**2 + dfdz%values**2 )

           call computeGradientKmain(wave,Wkin,tabk1,tabk2,tabk3)

           CALL btran(tabk1,dfdx,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk2,dfdy,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk3,dfdz,res)
           IF (.NOT.res) RETURN

  tab1%values  = tab1%values -2.0_WP*mu*( dfdx%values**2 + dfdy%values**2 + dfdz%values**2 )


!   CALL computePhysicalFilter(wave,tab1,DissSGSEQDM,delta,filter)
   CALL ftran(tab1,tabk1,res)
   IF (.NOT.res) RETURN

   call computeFilter(Wave, delta, tabk1, tabk2, filter)

   CALL btran(tabk2,DissSGSEQDM,res)
   IF (.NOT.res) RETURN


  DissSGSEQDM%values= DissSGSEQDM%values - DissGSEQDM%values 



  tab1%values = Ufin%values*Ufin%values + Vfin%values*Vfin%values + Wfin%values*Wfin%values
!   call computePhysicalGradient(wave,tab1,dfxdx,dfxdy,dfxdz)
           CALL ftran(tab1,Wtabk,res)
           call computeGradientKmain(wave,Wtabk,tabk1,tabk2,tabk3)

           CALL btran(tabk1,dfdx,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk2,dfdy,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk3,dfdz,res)
           IF (.NOT.res) RETURN
    
  AdvecGSEQDM%values = Ufin%values*dfdx%values + Vfin%values*dfdy%values + Wfin%values*dfdz%values

  tab1%values = dfdx%values
  tab2%values = dfdy%values
  tab3%values = dfdz%values


           CALL ftran(tab1,Wtabk,res)
           call computeGradientKmain(wave,Wtabk,tabk1,tabk2,tabk3)

           CALL btran(tabk1,dfdx,res)
           IF (.NOT.res) RETURN


  DifGSEQDM%values = mu*dfdx%values
       
           CALL ftran(tab2,Wtabk,res)
           call computeGradientKmain(wave,Wtabk,tabk1,tabk2,tabk3)



           CALL btran(tabk2,dfdy,res)
           IF (.NOT.res) RETURN


  DifGSEQDM%values = DifGSEQDM%values + mu*dfdy%values


           CALL ftran(tab3,Wtabk,res)
           call computeGradientKmain(wave,Wtabk,tabk1,tabk2,tabk3)


           CALL btran(tabk3,dfdz,res)
           IF (.NOT.res) RETURN

  DifGSEQDM%values = DifGSEQDM%values + mu*dfdz%values

  IF ((.NOT. copyStructOnly(Uin,DiffUSGSEqdm)) .OR. &
     &(.NOT. copyStructOnly(Uin,DissUSGSEqdm))  ) RETURN

     IF ((.NOT. copyStructOnly(Uin,S12)) .OR. &
        &(.NOT. copyStructOnly(Uin,S13)) .OR. &
        &(.NOT. copyStructOnly(Uin,S23)) .OR. &
        &(.NOT. copyStructOnly(Uin,S11)) .OR. &
        &(.NOT. copyStructOnly(Uin,S22)) .OR. &
        &(.NOT. copyStructOnly(Uin,S33)) ) RETURN
     IF ((.NOT. copyStructOnly(Uin,T12)) .OR. &
        &(.NOT. copyStructOnly(Uin,T13)) .OR. &
        &(.NOT. copyStructOnly(Uin,T23)) .OR. &
        &(.NOT. copyStructOnly(Uin,T11)) .OR. &
        &(.NOT. copyStructOnly(Uin,T22)) .OR. &
        &(.NOT. copyStructOnly(Uin,T33)) ) RETURN

  call computeStressTensor(Ufkin,Vfkin,Wfkin,wave,S11,S12,S13,S22,S23,S33,res)

  !call computeT_ijVecA(Ukin,VKin,WKin,Uin,Vin,Win,Ufin,Vfin,Wfin,wave,T11,T12,T13, &
  !          & T22,T23,T33,filter,res,delta)
  call computeT_ijVecA(Uin,Vin,Win,Ufin,Vfin,Wfin,wave,T11,T12,T13, &
            & T22,T23,T33,filter,res,delta)

  DiffUSGSEqdm%values = 0.0_WP

  tab1%values = 2*(Uin%values*T11%values + Vin%values*T12%values + Win%values*T13%values )
  tab2%values = 2*(Uin%values*T12%values + Vin%values*T22%values + Win%values*T23%values )
  tab3%values = 2*(Uin%values*T13%values + Vin%values*T23%values + Win%values*T33%values )

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

  DiffUSGSEqdm%values = -dfdx%values - dfdy%values - dfdz%values

  DissUSGSEqdm%values = 2.0_WP*(T11%values*S11%values + T22%values*S22%values + T33%values*S33%values) & 
             & + 4.0_WP*(T13%values*S13%values + T12%values*S12%values + T23%values*S23%values)

  CALL deleteDataLayout(dfdx)
  CALL deleteDataLayout(dfdy)
  CALL deleteDataLayout(dfdz)
  CALL deleteDataLayout(tab1)
  CALL deleteDataLayout(tab2)
  CALL deleteDataLayout(tab3)

  CALL deleteDataLayout(tabk3)
  CALL deleteDataLayout(tabk2)
  CALL deleteDataLayout(tabk1)

  CALL deleteDataLayout(Wtabk)

  CALL deleteDataLayout(S12)
  CALL deleteDataLayout(S13)
  CALL deleteDataLayout(S23)
  CALL deleteDataLayout(S11)
  CALL deleteDataLayout(S22)
  CALL deleteDataLayout(S33)

  CALL deleteDataLayout(T12)
  CALL deleteDataLayout(T13)
  CALL deleteDataLayout(T23)
  CALL deleteDataLayout(T11)
  CALL deleteDataLayout(T22)
  CALL deleteDataLayout(T33)

  meanDifGSEQDM    = computeFieldAvg( DifGSEQDM,  spec_rank,nbcpus)
  meanAdvecGSEQDM  = computeFieldAvg( AdvecGSEQDM, spec_rank,nbcpus)
  meanDissGSEQDM   = computeFieldAvg( DissGSEQDM, spec_rank,nbcpus)
  meanDiffUSGSEqdm = computeFieldAvg( DiffUSGSEqdm, spec_rank,nbcpus)
  meanDissUSGSEqdm = computeFieldAvg( DissUSGSEqdm, spec_rank,nbcpus)
  meanDissSGSEQDM  = computeFieldAvg( DissSGSEQDM, spec_rank,nbcpus)
  meanEkGS         = computeFieldAvg( EkGS, spec_rank,nbcpus)

  CALL deleteDataLayout(DifGSEQDM)
  CALL deleteDataLayout(AdvecGSEQDM)
  CALL deleteDataLayout(DissGSEQDM)
  CALL deleteDataLayout(DiffUSGSEqdm)
  CALL deleteDataLayout(DissUSGSEqdm)
  CALL deleteDataLayout(EkGS)
  CALL deleteDataLayout(DissSGSEQDM)

end subroutine EQDMGSHD


subroutine aprioriEQDHD(Uin,Vin,Win,Ukin,Vkin,Wkin,&
                    & Wave,spec_rank,delta,meanEqdmTUuu,nbcpus,filterType,type_avg)

  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Uin,Vin,Win
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Ukin,Vkin,Wkin
  TYPE(WaveNumbers), INTENT(IN) :: Wave
  
  TYPE(COMPLEX_DATA_LAYOUT) :: div1,div2,div3
  INTEGER :: spec_rank,nbcpus,numVel,filterType,type_avg
  TYPE(REAL_DATA_LAYOUT) :: T12,T13,T23
  TYPE(REAL_DATA_LAYOUT) :: T11,T22,T33
  TYPE(REAL_DATA_LAYOUT) :: S12,S13,S23
  TYPE(REAL_DATA_LAYOUT) :: S11,S22,S33
  TYPE(REAL_DATA_LAYOUT) :: Wtab
  logical :: res
  REAL(WP) :: delta,meanEqdmTUuu

  IF ((.NOT. copyStructOnly(Uin,S12)) .OR. &
     &(.NOT. copyStructOnly(Uin,S13)) .OR. &
     &(.NOT. copyStructOnly(Uin,S23)) .OR. &
     &(.NOT. copyStructOnly(Uin,S11)) .OR. &
     &(.NOT. copyStructOnly(Uin,S22)) .OR. &
     &(.NOT. copyStructOnly(Uin,S33)) .OR. &
     &(.NOT. copyStructOnly(Uin,Wtab)) .OR. &
     &(.NOT. copyStructOnly(Uin,T12)) .OR. &
     &(.NOT. copyStructOnly(Uin,T13)) .OR. &
     &(.NOT. copyStructOnly(Uin,T23))  .OR. &
     &(.NOT. copyStructOnly(Uin,T11)) .OR. &
     &(.NOT. copyStructOnly(Uin,T22)) .OR. &
     &(.NOT. copyStructOnly(Uin,T33))   ) RETURN 

  IF ((.NOT. copyStructOnly(Ukin,div1)) .OR. &
     &(.NOT. copyStructOnly(Ukin,div2)) .OR. &
     &(.NOT. copyStructOnly(Ukin,div3))   ) RETURN 

     numVel=1
     type_avg=0

     call smagorinskyVel(div1,div2,div3,Uin,Vin,Win,Ukin,Vkin,Wkin,wave,numVel,spec_rank,res, &
                          &  filterType,type_avg,T11,T12,T13,T22,T23,T33,delta)
    CALL deleteDataLayout(div3)
    CALL deleteDataLayout(div2)
    CALL deleteDataLayout(div1)

     call computeStressTensor(Ukin,VKin,WKin,wave,S11,S12,S13,S22,S23,S33,res)

     Wtab%values = 4.0_WP*(   S12%values*T12%values &
                          & + S13%values*T13%values &
                          & + S23%values*T23%values)&
               & + 2.0_WP*(   S11%values*T11%values &
                          & + S22%values*T22%values &
                          & + S23%values*T33%values )

    meanEqdmTUuu  = computeFieldAvg(Wtab, spec_rank,nbcpus)


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



end subroutine aprioriEQDHD



end module post_hd_ener_lib
!> @}
