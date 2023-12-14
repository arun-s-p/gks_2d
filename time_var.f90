!
!     ******************************************************************
!
      MODULE TIME_VAR
!
!     ******************************************************************
!     *                                                                *
!     *   DYNAMICALLY ALLOCATED VARIABLES FOR THE FLOW SOLUTION        *
!     *                                                                *
!     ******************************************************************
!

      REAL :: TIME,DTMIN,viztime,vtime
      REAL,   DIMENSION(:,:),      ALLOCATABLE :: DTL

      INTEGER :: CYC,NCYC,counter


      END MODULE TIME_VAR
