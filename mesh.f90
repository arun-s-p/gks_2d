!
!     ******************************************************************
!
      SUBROUTINE MESH

!     ******************************************************************
!

!
! --- ******************************************************************
!     *                                                                *
!     * SET LOCATIONS FOR VERTICES AND SET VOLUMES OF CELLS            *
!     *                                                                *
! --- ******************************************************************
!

      USE DIMS
      USE MESH_VAR
      USE FLO_PARAM

! --- ******************************************************************
 
      IMPLICIT NONE
 
! --- ******************************************************************


! --- ******************************************************************
!
!          LOCAL VARIABLES
!
! --- ******************************************************************


      REAL :: ETAXL,ETAXR
      REAL :: ETAY,pi
!      REAL :: DX,DY

      INTEGER :: I,J
      INTEGER :: IMSH

! --- ******************************************************************

      CALL ALLOC_MESH

! --- ******************************************************************
!
!       READ MESH FROM EXISTING DATA FILE
!
! --- ******************************************************************
      pi = 4.*ATAN(1.)
      DX   = 20.*VTHICK0/(FLOAT(NX) - 1.);
      XX(3) = 0.!-10.*VTHICK0
      DO I = 4,IL
         XX(I) = XX(I-1) + DX
      ENDDO

      DY   = 20.*VTHICK0/(FLOAT(NY) - 1.);
      YY(3) = 0.!-10.*VTHICK0
      DO J = 4,JL
         YY(J) = YY(J-1) + DY
      ENDDO

! --- CREATE GHOST CELL VERTICES BY REFLECTING THEM ABOUT THE LAST VERTEX
      XX(2)    = 2.*XX(3)  - XX(4)
      XX(1)    = 2.*XX(3)  - XX(5)

      XX(IE)   = 2.*XX(IL) - XX(IL-1)
      XX(IB)   = 2.*XX(IL) - XX(IL-2)

      YY(2)    = 2.*YY(3)  - YY(4)
      YY(1)    = 2.*YY(3)  - YY(5)

      YY(JE)   = 2.*YY(JL) - YY(JL-1)
      YY(JB)   = 2.*YY(JL) - YY(JL-2)

      DO I=1, IB
         WRITE(51,*)  I, XX(I)
      ENDDO

      DO J=1, JB
         WRITE(52,*)  J, YY(J)
      ENDDO

! --- COMPUTE VOLUMES OF CELLS

      DO J = 1, JE
         DO I = 1,IE
            VOL(I,J) = (XX(I+1) -XX(I))*(YY(J+1) -YY(J))
         ENDDO
      ENDDO

! --- DEFINE ONE 2D MESH VARIABLE CONTAINING THE VERTEX LOCATIONS

      DO J = 1,JB
         DO I = 1,IB
            X(1,I,J) = XX(I)
            X(2,I,J) = YY(J)
         ENDDO
      ENDDO

      END SUBROUTINE MESH
