
!
!     ******************************************************************
!
      SUBROUTINE ALLOC_MESH

!     ******************************************************************
!


      USE DIMS        ,ONLY : IE,JE,IB,JB
      USE MESH_VAR    
      USE TIME_VAR, ONLY: DTL

      IMPLICIT NONE


      ALLOCATE(XX(IB))
      ALLOCATE(YY(JB))
      ALLOCATE(X(2,IB,JB))
      ALLOCATE(VOL(IE,JE))
      ALLOCATE(DTL(IE,JE))

      END SUBROUTINE ALLOC_MESH
