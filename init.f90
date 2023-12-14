!
!     ******************************************************************
!
      SUBROUTINE INIT

!     ******************************************************************
!

      USE DIMS
      USE FLO_VAR
      USE FLO_PARAM
      USE MESH_VAR       ,ONLY: RLEN,X
      USE IODEFS         ,ONLY: IREAD,IOUT,FNAME0
      USE TIME_VAR       ,ONLY: TIME,CYC

!     ******************************************************************
!
      IMPLICIT NONE
!
!     ******************************************************************


!     ******************************************************************
!
!          LOCAL VARIABLES
!
!     ******************************************************************

      REAL , DIMENSION(4) :: DWRAND
      REAL                :: EPS,PI
      REAL                :: XC,YC,ETAX,ETAY,EE
      INTEGER :: I,J

!     ******************************************************************

      PI = 4.*ATAN(1.)
      IOUT = 17
      
      CALL ALLOC_FLO


      EPS = 0.02
         
      IF(IREAD==1) THEN !READ FROM PREVIOUS SOLUTION FOR RESTARTING
         OPEN (IOUT,FILE=FNAME0,ACCESS='SEQUENTIAL')         
         DO J = 1,JE
            DO I = 1,IE
               READ(IOUT,'(4(1X,E18.10))') W(1,I,J),W(2,I,J),W(3,I,J),W(4,I,J)
               RLV(I,J) = RMU0
            ENDDO
         ENDDO
         READ(IOUT,*) TIME
         CLOSE(IOUT)
      ENDIF

      IF(IREAD==0) THEN
         DO J = 1,JE
            DO I = 1,IE
               XC = .5*(X(1,I,J) + X(1,I+1,J))
               YC = .5*(X(2,I,J) + X(2,I,J+1))-10.*VTHICK0
               ETAY = 2.*YC/RLEN
               ETAX = 2.*PI*XC
               EE   = EXP(-ETAY*ETAY)

               W(1,I,J) = RHO0
               W(2,I,J) = 0.5*U0*TANH(2.*YC/VTHICK0) &
		-EE*.05*YC*20.*VTHICK0/(2.*PI*10.)*COS(4.*PI*XC/(20.*VTHICK0))*EXP(-YC**2/10.) &
		-EE*.025*YC*20.*VTHICK0/(PI*10.)*COS(2.*PI*XC/(20.*VTHICK0))*EXP(-YC**2/10.)
               W(3,I,J) = EE*.05*SIN(4.*PI*XC/(20.*VTHICK0))*EXP(-YC**2/10.) &
		-EE*.025*SIN(2.*PI*XC/(20.*VTHICK0))*EXP(-YC**2/10.)
!               W(2,I,J) = -0.5*U0*TANH(2.*YC/VTHICK0) &
!                -.05*YC*20.*VTHICK0/(2.*PI*10.)*COS(4.*PI*XC/(20.*VTHICK0))*EXP(-YC**2/10.) &
!                -.025*YC*20.*VTHICK0/(PI*10.)*COS(2.*PI*XC/(20.*VTHICK0))*EXP(-YC**2/10.)
!               W(3,I,J) = .05*SIN(4.*PI*XC/(20.*VTHICK0))*EXP(-YC**2/10.) &
!                -.025*SIN(2.*PI*XC/(20.*VTHICK0))*EXP(-YC**2/10.)
               W(4,I,J) = P0

               RLV(I,J) = RMU0
            END DO
         END DO
      ENDIF
      DO J = 1,JE
         DO I = 1,IE
            WRITE(30,*) I,J,W(1:4,I,J)
         ENDDO
      ENDDO

      CALL BCONDS
      CALL OUTPUT 
      END SUBROUTINE INIT
