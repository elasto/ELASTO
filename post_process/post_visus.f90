!> @addtogroup post_process
!! @{

!------------------------------------------------------------------------------

!
! MODULE: post_hd_out
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

module post_visus

    USE avs
    USE datalayout
    USE toolbox 
    USE stat_tools 
    USE differential_tools 
    USE data
    USE post_hd_out
    USE physical_values_tools
    USE transforms_tools
    USE post_lib
    implicit none
    private

     ! ==== Public procedures ==== 
    public          :: Helicity_Fieldsave
    !public          :: dumpSGSFieldsForScalar 
    !public          :: dumpUScalaireRotdTauijdxi

   contains
!======= compute all the grid and subgrid terms =========
!------------------------------------------------------------------------------
!> Compute and save mean values of all quantities in GS and SGS énergie equations
!! @author Mouloud KESSAR
!!    @param[in]    spec_rank     = mpi rank inside spectral communicator
!!    @param[in]    Uin,Vin,Win   = Velocity fields ins physical space
!!    @param[in]    Bin           = magnetic Field
!!    @param[in]    Wave          = wavenumber
!!    @param[in]    nbcpus        = numbers of cpus
!! @details
!!    
!------------------------------------------------------------------------------
subroutine Helicity_Fieldsave(Uin,Vin,Win,&
                         & nbcpus,spec_rank)

    implicit none

    TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Uin,Vin,Win
    TYPE(REAL_DATA_LAYOUT) :: Vortx,Vorty,Vortz
    TYPE(COMPLEX_DATA_LAYOUT) :: Vortxk,Vortyk,Vortzk
    TYPE(REAL_DATA_LAYOUT) :: Helicity

    INTEGER  :: nbcpus,spec_rank
    logical :: res

  call post_HD_init(Uin,Vin,Win,&
                         & nbcpus,spec_rank)
 
     IF ( (.NOT. copyStructOnly(Uin,Vortx ))) RETURN
     IF ( (.NOT. copyStructOnly(Uin,Vorty ))) RETURN
     IF ( (.NOT. copyStructOnly(Uin,Vortz ))) RETURN
     IF ( (.NOT. copyStructOnly(Uin,Helicity ))) RETURN



    call computeVorticity(Uk,Vk,Wk,VelWN,Vortx,Vorty,Vortz,res)

    Helicity%values=0.5_WP*( Uin%values*Vortx%values &
                         & + Vin%values*Vorty%values &
                         & + Win%values*Vortz%values )
!      labeltab(1)=TRIM(ADJUSTL(U%name))
!      if (withname) then
       IF(.NOT. real_WriteToFile(spec_rank,'helicitykin.data',Helicity)) RETURN
     res=computeFieldPDF(1,Helicity,1000,spec_rank,1)

     IF ( (.NOT. copyStructOnly(Uk,Vortxk ))) RETURN
     IF ( (.NOT. copyStructOnly(Uk,Vortyk ))) RETURN
     IF ( (.NOT. copyStructOnly(Uk,Vortzk ))) RETURN
          CALL ftran(Vortx,Vortxk,res)
          IF (.NOT.res) RETURN
          CALL ftran(Vorty,Vortyk,res)
          IF (.NOT.res) RETURN
          CALL ftran(Vortz,Vortzk,res)
          IF (.NOT.res) RETURN

  call compute_Helicity_spectrum(Uk,Vk,Wk,&
                                  &Vortxk,Vortyk,Vortzk,VelWN, spec_rank) 

  CALL deleteDataLayout(Vortx)
  CALL deleteDataLayout(Vorty)
  CALL deleteDataLayout(Vortz)
  CALL deleteDataLayout(Vortxk)
  CALL deleteDataLayout(Vortyk)
  CALL deleteDataLayout(Vortzk)

end subroutine Helicity_Fieldsave

end module post_visus
!> @}
