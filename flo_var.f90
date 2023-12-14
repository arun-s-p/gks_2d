!
!     ******************************************************************
!
      MODULE FLO_VAR
!
!     ******************************************************************
!     *                                                                *
!     *   DYNAMICALLY ALLOCATED VARIABLES FOR THE FLOW SOLUTION        *
!     *                                                                *
!     ******************************************************************
!

      REAL, DIMENSION(:,:,:)    , ALLOCATABLE :: W
      REAL, DIMENSION(:,:,:,:)  , ALLOCATABLE :: WX
      REAL, DIMENSION(:,:)      , ALLOCATABLE :: RLV

      REAL	:: RHOM,UM

      END MODULE FLO_VAR
