!
!     ******************************************************************
!
      SUBROUTINE BCONDS

!     ******************************************************************
!

!
! --- ******************************************************************
!     *                                                                *
!     * SET BOUNDARY CONDITIONS VIA HALO CELLS                         *
!     *                                                                *
! --- ******************************************************************
!

      USE DIMS              ,ONLY : IL,JL,IE,JE
      USE FLO_VAR           ,ONLY : W

! --- ******************************************************************
 
      IMPLICIT NONE
 
! --- ******************************************************************


! --- ******************************************************************
!
!          LOCAL VARIABLES
!
! --- ******************************************************************

      INTEGER :: I,J,N
      INTEGER :: JHALO,JINT


! --- ******************************************************************
!     *      SET SLIP CONDITIONS AT THE BOTTOM AND AT THE TOP          *
! --- ******************************************************************

      DO I = 3, IL -1

         JHALO = 2; JINT = 3
         W(1:4,I,JHALO) =     W(1:4,I,JINT)
         W(3,I,JHALO)   =    -W(3,I,JINT)

         JHALO = 1; JINT = 4
         W(1:4,I,JHALO) =     W(1:4,I,JINT)
         W(3,I,JHALO)   =    -W(3,I,JINT)


         JHALO = JL; JINT = JL-1
         W(1:4,I,JHALO) =     W(1:4,I,JINT)
         W(3,I,JHALO)   =    -W(3,I,JINT)

         JHALO = JE; JINT = JL-2
         W(1:4,I,JHALO) =     W(1:4,I,JINT)
         W(3,I,JHALO)   =    -W(3,I,JINT)
      ENDDO

! --- ******************************************************************
!     *        SET PERIODIC BOUNDARY CONDITIONS ON LEFT AND RIGHT      *
! --- ******************************************************************


      DO J = 3,JL-1
! ---  LEFT
         W(1:4,2,J) = W(1:4,IL-1,J)
         W(1:4,1,J) = W(1:4,IL-2,J)

! ---  RIGHT EXIT
         W(1:4,IL,J) = W(1:4,3,J)
         W(1:4,IE,J) = W(1:4,4,J)
      ENDDO
!!$
!!$      DO J = 1,JE
!!$         DO I = 1,IE
!!$            WRITE(31,12) I,J,W(1:4,I,J)
!!$         ENDDO
!!$      ENDDO
!!$
!!$12    FORMAT(1X,2I6,4E18.10)
      END SUBROUTINE BCONDS
      
