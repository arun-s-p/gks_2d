!
!     ******************************************************************
!
      MODULE MESH_VAR
!
!     ******************************************************************
!     *                                                                *
!     *   MESH VARIABLES                                               *
!     *                                                                *
!     ******************************************************************
!
      REAL,   DIMENSION(:),        ALLOCATABLE :: XX,YY
      REAL,   DIMENSION(:,:,:),    ALLOCATABLE :: X
      REAL,   DIMENSION(:,:),      ALLOCATABLE :: VOL

      REAL :: RLEN,DX,DY

      END MODULE MESH_VAR
