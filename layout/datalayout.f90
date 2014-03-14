!USEFORTEST toolbox
!USEFORTEST avgcond 
!USEFORTEST postprocess
!USEFORTEST advec
!USEFORTEST io
!USEFORTEST topo
!> @addtogroup layout
!! @{
!------------------------------------------------------------------------------
!
! MODULE: datalayout
!
!> @author
!> Patrick Begou, LEGI
!
! DESCRIPTION:
!> The aim of this module is to provide global functionalities and parameters
!> to the underliying datalayout architecture. This module provides generic
!> access to the functionalities. The hierarchical organisation is shown in the 
!> figure. The generic implementation is located in the file implementdatalayout.Fortran.
!> It is processed to generate realdatalayout.F90 (module realdatalayout for REAL datas) and
!> cmplxdatalayout.F90 (module cmplxdatalayout for COMPLEX data).
!>\image html datalayout.png "Hierarchical organisation of the datalayout modules and files"
!------------------------------------------------------------------------------

MODULE datalayout
USE cmplxdatalayout
USE realdatalayout
USE cmplxVectordatalayout
USE realVectordatalayout

INTERFACE initDataLayout
   module procedure real_initDataLayout
   module procedure cmplx_initDataLayout
   module procedure real_vector_initDataLayout
   module procedure cmplx_vector_initDataLayout
END INTERFACE initDataLayout

INTERFACE deleteDataLayout
   module procedure real_deleteDataLayout
   module procedure cmplx_deleteDataLayout
   module procedure real_vector_deleteDataLayout
   module procedure cmplx_vector_deleteDataLayout
END INTERFACE deleteDataLayout

INTERFACE showLayout
   module procedure real_showLayout
   module procedure cmplx_showLayout
   module procedure real_vector_showLayout
   module procedure cmplx_vector_showLayout
END INTERFACE showLayout

INTERFACE sameLayout
   module procedure real_sameLayout
   module procedure cmplx_sameLayout
   module procedure real_vector_sameLayout
   module procedure cmplx_vector_sameLayout
END INTERFACE sameLayout

INTERFACE copyStructOnly
   module procedure real_copyStructOnly
   module procedure cmplx_copyStructOnly
   module procedure real_vector_copyStructOnly
   module procedure cmplx_vector_copyStructOnly
END INTERFACE copyStructOnly

INTERFACE permute
  module procedure real_permute
  module procedure cmplx_permute
END INTERFACE permute

INTERFACE foundProcess
  module procedure real_foundProcess
  module procedure cmplx_foundProcess
END INTERFACE foundProcess

END MODULE datalayout
!> @}
