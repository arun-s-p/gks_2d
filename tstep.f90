!
!     ******************************************************************
!
      SUBROUTINE TSTEP

!     ******************************************************************
!

!
! --- ******************************************************************
!     *                                                                *
!     * COMPUTE TIME STEP IN EACH CELL                                 *
!     *                                                                *
! --- ******************************************************************
!

      USE DIMS                        ,ONLY : IL,JL,IE,JE
      USE MESH_VAR                    ,ONLY : XX,YY,X,VOL
      USE TIME_VAR                   ! ,ONLY : DTL,DTMIN
      USE FLO_VAR                     ,ONLY : W,RLV
      USE FLO_PARAM                   ,ONLY : GAMMA,PRN,KVIS,CFL

! --- ******************************************************************
 
      IMPLICIT NONE
 
! --- ******************************************************************


! --- ******************************************************************
!
!          LOCAL VARIABLES
!
! --- ******************************************************************


      REAL :: DX,DY,DS
      REAL :: A
      REAL :: PRG
      REAL :: PA,RA
      REAL :: CS,QS,safe

      INTEGER :: I,J,ramp,order
! --- ******************************************************************


! --- INITIALIZE TIME STEP IN ALL CELLS

      DO J = 1,JE
         DO I = 1,IE
            DTL(I,J) = 0.
         ENDDO
      ENDDO

 
! --- SPECTRAL RADIUS IN I DIRECTION

      DO J = 3, JL -1
         DO I = 3, IL
            PA = W(4,I-1,J) + W(4,I,J)
            RA = W(1,I-1,J) + W(1,I,J)
            DS = YY(J+1) - YY(J)
!            DS = X(2,I,J+1) - X(2,I,J)

            CS = SQRT(GAMMA*PA/RA)*DS
            QS = .5*( W(2,I-1,J) + W(2,I,J) ) * DS
            A  = ABS(QS) + CS

            DTL(I-1,J) = DTL(I-1,J) + A
            DTL(I  ,J) = DTL(I  ,J) + A
         ENDDO
      ENDDO

! --- SPECTRAL RADIUS IN J DIRECTION
      DO J = 3, JL
         DO I = 3, IL -1
            PA = W(4,I,J-1)   + W(4,I,J)
            RA = W(1,I,J-1) + W(1,I,J)
            DS = XX(I+1) - XX(I)
!            DS = X(1,I+1,J) - X(1,I,J)

            CS = SQRT(GAMMA*PA/RA)*DS
            QS = .5*( W(3,I,J-1) + W(3,I,J)) * DS
            A  = ABS(QS) + CS

            DTL(I,J-1) = DTL(I,J-1) + A
            DTL(I,J  ) = DTL(I,J  ) + A
         ENDDO
      ENDDO

! --- IN CASE OF VISCOUS FLOW TAKE INTO ACCOUNT VISCOSITY
      
      IF(KVIS > 0 ) THEN
!     IF(KVIS > 1 ) THEN

         PRG = GAMMA/PRN

! --- SPECTRAL RADIUS IN I DIRECTION

         DO J = 3, JL -1
            DO I = 3, IL
               DTL(I-1,J) = DTL(I-1,J) + &
                            2.*PRG*RLV(I-1,J)/(W(1,I-1,J)*(XX(I)-XX(I-1)))
               DTL(I  ,J) = DTL(I  ,J) + &
                            2.*PRG*RLV(I  ,J)/(W(1,I  ,J)*(XX(I+1)-XX(I)))
            ENDDO
         ENDDO
         
! --- SPECTRAL RADIUS IN J DIRECTION
         DO J = 3, JL
            DO I = 3, IL -1
               DTL(I,J-1) = DTL(I,J-1) + &
                            2.*PRG*RLV(I,J-1)/(W(1,I,J-1)*(YY(J)-YY(J-1)))
               DTL(I,J  ) = DTL(I  ,J) + &
                            2.*PRG*RLV(I,J  )/(W(1,I,J  )*(YY(J+1)-YY(J)))
            ENDDO
         ENDDO
      ENDIF

! --- DIVIDE BY RESPECTIVE VOLUMES

      DO J = 3, JL - 1
         DO I = 3, IL -1
            DTL(I,J) = 4.*VOL(I,J)/DTL(I,J)
         ENDDO
      ENDDO



! --- SET TIME STEP IN HALO CELLS
      DO J = 3, JL
         DTL(1,J )  = DTL(4,J)
         DTL(2,J )  = DTL(3,J)
         DTL(IL,J)  = DTL(IL-1,J)
         DTL(IE,J)  = DTL(IL-2,J)
      ENDDO

      DO I = 3, IL
         DTL(I,1 )  = DTL(I,4)
         DTL(I,2 )  = DTL(I,3)
         DTL(I,JL)  = DTL(I,JL-1)
         DTL(I,JE)  = DTL(I,JL-2)
      ENDDO

      DTL   = CFL*DTL
! --- FIND MINIMUM TIME STEP
      DTMIN = DTL(3,3)
      DO J = 3,JL -1
         DO I = 3, IL -1
            IF(DTL(I,J) < DTMIN) DTMIN = DTL(I,J)
         ENDDO
      ENDDO

! --- CONSTANT TIME STEP
      DO J = 1,JE
         DO I = 1,IE
            DTL(I,J) = DTMIN
         ENDDO
      ENDDO

 !print *, DTMIN
 

      END SUBROUTINE TSTEP
