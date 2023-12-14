!
!     ******************************************************************
!
      MODULE SOLV_VAR
!
!     ******************************************************************
!     *                                                                *
!     *   DYNAMICALLY ALLOCATED VARIABLES FOR THE FLOW                 *
!     *                                                                *
!     ******************************************************************
!

      REAL, DIMENSION(:,:,:)    , ALLOCATABLE :: DW,FW,VW
      REAL, DIMENSION(:,:,:)    , ALLOCATABLE :: W0

      END MODULE SOLV_VAR
