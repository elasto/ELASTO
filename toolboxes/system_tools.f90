!> @addtogroup toolbox 
!! @{
!------------------------------------------------------------------------------
!
! MODULE: system_tools 
!
!! DESCRIPTION:
!> This module implements dealiasing methods.
!> @author
!> Guillaume Balarac et Patrick BEGOU, LEGI
!>
!> @details
!! This module implements 3 choices for dealiasing: none, isotropic and anisotropic.
!! The last one is usefull when the number point in the 3 dimensions are very heterogeneous
!! to avoid loosing to much informations.
!! \image html dealias.png "Influence of dealiasing method on scalar min/max behavior."
!> @author
!> Guillaume Balarac et Patrick BEGOU, LEGI
!------------------------------------------------------------------------------
module system_tools 

    implicit none

    public  

contains

!------------------------------------------------------------------------------
!> Show memory used by the simulation - to trak memory leak
!!    @param[in]    spec_rank   = mpi rank inside the spectral communicator
!! @autor Jean-Baptiste Lagaert
!------------------------------------------------------------------------------
subroutine show_memory(spec_rank)

    integer, intent(in)                                 :: spec_rank

    ! Todo : detecter la machine (UNIX-like ou vargas/babel et switcher sur la
    ! bonne ligne)

    if (spec_rank ==0) then
        write(*,'(a)') '   [MEM] Memory used (% and absolute value), and cpu and user'
        call system('ps ww -A -o %mem,vsz,user,command | grep elasto.exe | grep -v grep | head -2 ')
        ! For vargas and babel, activate the following line
        !call system('ps -AkF " %z : < %u > %a" | grep scaleExe | grep -v grep | head -1 ')
    end if

end subroutine show_memory

end module system_tools 
!> @}
