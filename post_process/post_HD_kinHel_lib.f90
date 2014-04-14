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

module post_hd_kinHel_lib

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
    public          :: EqHkinGS
    public          :: EqHkin
!     public          :: aprioriEQDHD
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

subroutine EqHkin(Uin,Vin,Win,Ukin,VKin,WKin,Wave,Vortx,Vorty,Vortz,Vortxk,Vortyk,Vortzk,&
                 &nbcpus,spec_rank,mu,meanadvecHkin, meandiffvisHkin, meandissvisHkin, meanHkin)


   implicit none

   integer, intent(in)    :: nbcpus,spec_rank
   TYPE(REAL_DATA_LAYOUT), INTENT(IN)          :: Uin,Vin,Win
   TYPE(REAL_DATA_LAYOUT), INTENT(IN)          :: Vortx,Vorty,Vortz
   TYPE(WaveNumbers), INTENT(IN) :: Wave
   TYPE(REAL_DATA_LAYOUT) :: tab1, tab2, tab3
   TYPE(REAL_DATA_LAYOUT) :: dfdx, dfdy, dfdz
   TYPE(REAL_DATA_LAYOUT) :: df2dx, df2dy, df2dz
   TYPE(REAL_DATA_LAYOUT) :: advecHkin,diffvisHkin,dissvisHkin,Hkin
   REAL(WP) :: meanadvecHkin, meandiffvisHkin, meandissvisHkin, meanHkin
   REAL(WP) :: mu
   TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Ukin,Vkin,Wkin
   TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Vortxk,Vortyk,Vortzk
   logical :: res
   TYPE(COMPLEX_DATA_LAYOUT) :: tabk1,tabk2,tabk3,Wtabk

     IF ((.NOT. copyStructOnly(Uin,advecHkin)) .OR. &
        &(.NOT. copyStructOnly(Uin,diffvisHkin)) .OR. &
        &(.NOT. copyStructOnly(Uin,dissvisHkin)) .OR. &
        &(.NOT. copyStructOnly(Uin,Hkin)) .OR. &
        &(.NOT. copyStructOnly(Uin,dfdx)) .OR. &
        &(.NOT. copyStructOnly(Uin,dfdy)) .OR. &
        &(.NOT. copyStructOnly(Uin,dfdz)) .OR. &
        &(.NOT. copyStructOnly(Uin,df2dx)) .OR. &
        &(.NOT. copyStructOnly(Uin,df2dy)) .OR. &
        &(.NOT. copyStructOnly(Uin,df2dz)) .OR. &
        &(.NOT. copyStructOnly(Uin,tab1)) .OR. &
        &(.NOT. copyStructOnly(Uin,tab2)) .OR. &
        &(.NOT. copyStructOnly(Uin,tab3)) ) RETURN


  
     IF ( (.NOT. copyStructOnly(Ukin,tabk1 ))) RETURN
     IF ( (.NOT. copyStructOnly(Ukin,tabk2 ))) RETURN
     IF ( (.NOT. copyStructOnly(Ukin,tabk3 ))) RETURN
     IF ( (.NOT. copyStructOnly(Ukin,Wtabk ))) RETURN

     !!Helicitée cinétique!!
    Hkin%values = 0.5_WP*( Uin%values*Vortx%values + Vin%values*Vorty%values + Win%values*Vortz%values)

    !!!Advection
!     call computePhysicalGradient(Wave,EnQDM,dfxdx,dfxdy,dfxdz)
           CALL ftran(Hkin,Wtabk,res)
           call computeGradientKmain(wave,Wtabk,tabk1,tabk2,tabk3)

           CALL btran(tabk1,dfdx,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk2,dfdy,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk3,dfdz,res)
           IF (.NOT.res) RETURN
    
    advecHkin%values = -Uin%values*dfdx%values -Vin%values*dfdy%values -Win%values*dfdz%values
 
  
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
       diffvisHkin%values = mu*dfdx%values

           CALL ftran(tab2,Wtabk,res)
           call computeGradientKmain(wave,Wtabk,tabk1,tabk2,tabk3)

           CALL btran(tabk1,dfdx,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk2,dfdy,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk3,dfdz,res)
           IF (.NOT.res) RETURN
       diffvisHkin%values = diffvisHkin%values + mu*dfdy%values

           CALL ftran(tab3,Wtabk,res)
           call computeGradientKmain(wave,Wtabk,tabk1,tabk2,tabk3)

           CALL btran(tabk1,dfdx,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk2,dfdy,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk3,dfdz,res)
           IF (.NOT.res) RETURN

       diffvisHkin%values = diffvisHkin%values + mu*dfdz%values
  
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

           call computeGradientKmain(wave,VortxK,tabk1,tabk2,tabk3)

           CALL btran(tabk1,df2dx,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk2,df2dy,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk3,df2dz,res)
           IF (.NOT.res) RETURN

dissvisHkin%values = -2*mu*(dfdx%values*df2dx%values + dfdy%values*df2dy%values + dfdz%values*df2dz%values )
  
           call computeGradientKmain(wave,Vkin,tabk1,tabk2,tabk3)
       
           CALL btran(tabk1,dfdx,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk2,dfdy,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk3,dfdz,res)
           IF (.NOT.res) RETURN

           call computeGradientKmain(wave,VortyK,tabk1,tabk2,tabk3)

           CALL btran(tabk1,df2dx,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk2,df2dy,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk3,df2dz,res)
           IF (.NOT.res) RETURN

dissvisHkin%values = dissvisHkin%values -2*mu*(dfdx%values*df2dx%values + dfdy%values*df2dy%values + dfdz%values*df2dz%values )

           call computeGradientKmain(wave,Wkin,tabk1,tabk2,tabk3)

           CALL btran(tabk1,dfdx,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk2,dfdy,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk3,dfdz,res)
           IF (.NOT.res) RETURN
           call computeGradientKmain(wave,VortzK,tabk1,tabk2,tabk3)

           CALL btran(tabk1,df2dx,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk2,df2dy,res)
           IF (.NOT.res) RETURN
           CALL btran(tabk3,df2dz,res)
           IF (.NOT.res) RETURN

dissvisHkin%values = dissvisHkin%values -2*mu*(dfdx%values*df2dx%values + dfdy%values*df2dy%values + dfdz%values*df2dz%values )

 
  
!     if (rank .eq. 0) print*,'Mean calculation'
    meanadvecHkin   = computeFieldAvg(advecHkin,spec_rank,nbcpus) 
    meandiffvisHkin = computeFieldAvg(diffvisHkin,spec_rank,nbcpus) 
    meandissvisHkin = computeFieldAvg(dissvisHkin,spec_rank,nbcpus) 
    meanHkin = computeFieldAvg(Hkin,spec_rank,nbcpus) 


  CALL deleteDataLayout(advecHkin)
  CALL deleteDataLayout(diffvisHkin)
  CALL deleteDataLayout(dissvisHkin)
  CALL deleteDataLayout(Hkin)
  CALL deleteDataLayout(dfdx)
  CALL deleteDataLayout(dfdy)
  CALL deleteDataLayout(dfdz)
  CALL deleteDataLayout(tab1)
  CALL deleteDataLayout(tab2)
  CALL deleteDataLayout(tab3)
  CALL deleteDataLayout(df2dx)
  CALL deleteDataLayout(df2dy)
  CALL deleteDataLayout(df2dz)

  CALL deleteDataLayout(tabk1)
  CALL deleteDataLayout(tabk2)
  CALL deleteDataLayout(tabk3)

  CALL deleteDataLayout(Wtabk)

end subroutine EqHkin
 


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


subroutine EqHkinGS(Ufin,Vfin,Wfin,Uin,Vin,Win,Ufkin,Vfkin,Wfkin,Ukin,VKin,WKin,&
                   & Vortx,Vorty,Vortz,Vortxf,Vortyf,Vortzf,&
                   & Vortxk,Vortyk,Vortzk,Vortxfk,Vortyfk,Vortzfk,&
                   & Wave,delta,filter,nbcpus,spec_rank,mu,&
                   & meanDissGSEqHkin, meanAdvecGSEqHkin, meanDifGSEqHkin,&
                   & meanDissTUWUijSGSHkin, meanDiffTUWSGSHkin,meanDiffTUSGSHkin, meanDissTUSwSGSHkin, meanHkinGS)

  IMPLICIT NONE
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Ufin,Vfin,Wfin
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Uin,Vin,Win
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Ukin,Vkin,Wkin
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Ufkin,Vfkin,Wfkin
  TYPE(REAL_DATA_LAYOUT) :: DissTUWUijSGSHkin, DiffTUWSGSHkin
  TYPE(REAL_DATA_LAYOUT) :: DiffTUSGSHkin, DissTUSwSGSHkin
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Vortx,Vorty,Vortz,Vortxf,Vortyf,Vortzf
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Vortxk,Vortyk,Vortzk,Vortxfk,Vortyfk,Vortzfk
  REAL(WP) :: delta
  REAL(WP) :: meanDissGSEqHkin, meanAdvecGSEqHkin, meanDifGSEqHkin
  REAL(WP) :: meanHkinGS
  REAL(WP) :: meanDissTUWUijSGSHkin, meanDiffTUWSGSHkin
  REAL(WP) :: meanDiffTUSGSHkin, meanDissTUSwSGSHkin
  INTEGER, INTENT(IN) :: nbcpus, spec_rank, filter
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

  call EqHkin(Ufin,Vfin,Wfin,Ufkin,Vfkin,Wfkin,Wave,Vortxf,Vortyf,Vortzf,Vortxfk,Vortyfk,Vortzfk,&
                 &nbcpus,spec_rank,mu,meanAdvecGSEqHkin, meanDifGSEqHkin, meanDissGSEqHkin, meanHkinGS)

     IF ((.NOT. copyStructOnly(Uin,dfdx)) .OR. &
        &(.NOT. copyStructOnly(Uin,dfdy)) .OR. &
        &(.NOT. copyStructOnly(Uin,dfdz)) .OR. &
        &(.NOT. copyStructOnly(Uin,tab1)) .OR. &
        &(.NOT. copyStructOnly(Uin,tab2)) .OR. &
        &(.NOT. copyStructOnly(Uin,tab3)) ) RETURN

     IF ( (.NOT. copyStructOnly(Ukin,tabk1 ))) RETURN
     IF ( (.NOT. copyStructOnly(Ukin,tabk2 ))) RETURN
     IF ( (.NOT. copyStructOnly(Ukin,tabk3 ))) RETURN
     IF ( (.NOT. copyStructOnly(Ukin,Wtabk ))) RETURN

  IF ((.NOT. copyStructOnly(Uin,DiffTUSGSHkin)) .OR. &
     &(.NOT. copyStructOnly(Uin,DissTUSwSGSHkin))  ) RETURN


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

  call computeStressTensor(Vortxfk,Vortyfk,Vortzfk,wave,S11,S12,S13,S22,S23,S33,res)

!  call computeT_ijVecA(Vortxk,Vortyk,Vortzk,Vortx,Vorty,Vortz,Vortxf,Vortyf,Vortzf,wave,T11,T12,T13, &
!            & T22,T23,T33,filter,res,delta)
  call computeT_ijVecA(Vortx,Vorty,Vortz,Vortxf,Vortyf,Vortzf,wave,T11,T12,T13, &
            & T22,T23,T33,filter,res,delta)


  tab1%values = 2.0_WP*(Uin%values*T11%values + Vin%values*T12%values + Win%values*T13%values )
  tab2%values = 2.0_WP*(Uin%values*T12%values + Vin%values*T22%values + Win%values*T23%values )
  tab3%values = 2.0_WP*(Uin%values*T13%values + Vin%values*T23%values + Win%values*T33%values )

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

  DiffTUSGSHkin%values = -dfdx%values - dfdy%values - dfdz%values

  DissTUSwSGSHkin%values = 2.0_WP*(T11%values*S11%values + T22%values*S22%values + T33%values*S33%values) & 
             & + 4.0_WP*(T13%values*S13%values + T12%values*S12%values + T23%values*S23%values)


  CALL deleteDataLayout(S11)
  CALL deleteDataLayout(S22)
  CALL deleteDataLayout(S33)

  CALL deleteDataLayout(T11)
  CALL deleteDataLayout(T22)
  CALL deleteDataLayout(T33)

  IF ((.NOT. copyStructOnly(Uin,DiffTUWSGSHkin)) .OR. &
     &(.NOT. copyStructOnly(Uin,DissTUWUijSGSHkin))  ) RETURN

  call computeJijTensor(Ufkin,Vfkin,Wfkin,wave,S12,S13,S23, &
             & res)

  call computeT_ijVecAVec_B(Ukin,Uin,Vin,Win,Ufin,Vfin,Wfin,&
                          & Vortx,Vorty,Vortz,Vortxf,Vortyf,Vortzf,&
                          & wave,T12,T13,T23,filter,res,delta)
  tab1%values = 2.0_WP*( Vin%values*T12%values + Win%values*T13%values )
  tab2%values = 2.0_WP*(-Uin%values*T12%values + Win%values*T23%values )
  tab3%values = 2.0_WP*(-Uin%values*T13%values - Vin%values*T23%values )

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


  DiffTUWSGSHkin%values= -dfdx%values - dfdy%values - dfdz%values
  DissTUWUijSGSHkin%values= 4.0_WP*(T13%values*S13%values + T12%values*S12%values + T23%values*S23%values)

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

  CALL deleteDataLayout(T12)
  CALL deleteDataLayout(T13)
  CALL deleteDataLayout(T23)

  meanDiffTUSGSHkin = computeFieldAvg( DiffTUSGSHkin, spec_rank,nbcpus)
  meanDissTUSwSGSHkin = computeFieldAvg( DissTUSwSGSHkin, spec_rank,nbcpus)

  meanDiffTUWSGSHkin = computeFieldAvg( DiffTUWSGSHkin, spec_rank,nbcpus)
  meanDissTUWUijSGSHkin = computeFieldAvg( DissTUWUijSGSHkin, spec_rank,nbcpus)

  CALL deleteDataLayout(DiffTUWSGSHkin)
  CALL deleteDataLayout(DissTUWUijSGSHkin)
  CALL deleteDataLayout(DiffTUSGSHkin)
  CALL deleteDataLayout(DissTUSwSGSHkin)

end subroutine EqHkinGS


end module post_hd_kinHel_lib
!> @}
