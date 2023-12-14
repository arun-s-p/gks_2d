!
!     ******************************************************************
!
      SUBROUTINE ALLOC_FLO
!
!     ******************************************************************
!     *                                                                *
!     *   ALLOCATE VARIABLES FOR THE FLOW SOLUTION                     *
!     *                                                                *
!     ******************************************************************
!
      USE DIMS                 ,ONLY : IE,JE
      USE FLO_VAR
      USE SOLV_VAR

      IMPLICIT NONE

      ALLOCATE(W(4,IE,JE))
      ALLOCATE(WX(2,4,IE,JE))
      ALLOCATE(RLV(IE,JE))

      ALLOCATE(W0(4,IE,JE))
      ALLOCATE(DW(4,IE,JE))
      ALLOCATE(FW(4,IE,JE))
      ALLOCATE(VW(4,IE,JE))
      END SUBROUTINE ALLOC_FLO
