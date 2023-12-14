!
!     ******************************************************************
!
      SUBROUTINE DERIVS_BGK

!     ******************************************************************
!

!
! --- ******************************************************************
!     *                                                                *
!     * COMPUTE RESIDUALS IN EACH CELL USING NS EQUATIONS              *
!     *                                                                *
! --- ******************************************************************
!

      USE DIMS              ,ONLY : IL,JL
      USE MESH_VAR          ,ONLY : X,VOL
      USE FLO_VAR           ,ONLY : W,WX
      USE SOLV_VAR          ,ONLY : DW,FW,VW
      USE FLO_PARAM         !,ONLY : GAMMA,PRN,RMU0
      USE TIME_VAR          ,ONLY : DTMIN
! --- ******************************************************************
 
      IMPLICIT NONE
 
! --- ******************************************************************


! --- ******************************************************************
!
!          LOCAL VARIABLES
!
! --- ******************************************************************

      REAL , PARAMETER :: EPS = 1.E-7 
      REAL             :: DWL,DWR
      REAL             :: DSL,DSR,DX,DY
      REAL	       :: XY
!      REAL             :: DX,XY

      REAL , DIMENSION(4) :: FC,FD,FV,FCP,FDP,FVP
      REAL , DIMENSION(4) :: WL,WR,WLP,WRP,W1P,W2P
      REAL , DIMENSION(4) :: DWT0
      REAL                :: MU
      INTEGER :: I,J,N
      INTEGER :: IP1,JP1
      REAL, DIMENSION(2) :: PLW, PRW, betaL, betaR, alphaL, alphaR
      REAL, DIMENSION(2) :: omegaL, omegaR
      REAL :: SMV
      
      SMV = 1E-6

      
! --- ******************************************************************

      DW = 0.
      FW = 0.
      VW = 0.
      WX = 0.
      !epsw = 1.0E-5

! --- FIND LIMITED FLOW SLOPES IN I & J DIRECTIONS IN EACH CELL
      DWT0 = 0.

!------------------------------------------------------------
!
! --- FIND FLUX CONTRIBUTIONS IN I DIRECTION
!
!------------------------------------------------------------

      DO J = 2,JL -1
         DO I = 2,IL -1
! ---   LENGTH OF FACE I + 1/2
            DY   = X(2,I,J+1) - X(2,I,J)
            DSL  = .5*(X(1,I+1,J) - X(1,I,J))
            DSR  = .5*(X(1,I+2,J) - X(1,I+1,J))

!----------WENO interpolation
            DO N = 1,4

! RECONSTRUCT WL AND WR FOR N=1,4 FROM W - XDIRECTION I INDEX

                PLW(1) = 0.5*W(N,I,J)+0.5*W(N,I+1,J)
                PLW(2) = -0.5*W(N,I-1,J)+1.5*W(N,I,J)

                PRW(1) = 1.5*W(N,I+1,J)-0.5*W(N,I+2,J)
                PRW(2) = 0.5*W(N,I,J)+0.5*W(N,I+1,J)

                betaL(1) = (W(N,I+1,J)-W(N,I,J))**2
                betaL(2) = (W(N,I,J)-W(N,I-1,J))**2
            
                betaR(1) = (W(N,I+2,J)-W(N,I+1,J))**2
                betaR(2) = (W(N,I+1,J)-W(N,I,J))**2

                alphaL(1) = 2.0/(3.0*(SMV+betaL(1))**2)
                alphaL(2) = 1.0/(3.0*(SMV+betaL(2))**2)

                alphaR(1) = 1.0/(3.0*(SMV+betaR(1))**2)
                alphaR(2) = 2.0/(3.0*(SMV+betaR(2))**2)

                omegaL(1) = alphaL(1)/(alphaL(1)+alphaL(2))
                omegaL(2) = alphaL(2)/(alphaL(1)+alphaL(2))

                omegaR(1) = alphaR(1)/(alphaR(1)+alphaR(2))
                omegaR(2) = alphaR(2)/(alphaR(1)+alphaR(2))

                WL(N) = omegaL(1)*PLW(1)+omegaL(2)*PLW(2)
                WR(N) = omegaR(1)*PRW(1)+omegaR(2)*PRW(2)
            ENDDO
            W1P  = W(1:4,I,J)
            W2P  = W(1:4,I+1,J)

! --- NO NEED TO ROTATE VARIABLES AS IT IS THE I FACE

! --- COMPUTE CONVECTIVE FLUXES
      CALL BGKFLUX(WL,WR,W1P,W2P,DSL,DSR,DTMIN,RMU0,FC,GAMMA,PRN)
            FC = FC*DY
!!$            FD = FD*DY
!!$
!!$! --- COMPUTE VISCOUS FLUXES
!!$!            MU = 0.5*(RLV(I,J)+ RLV(I+1,J))
!!$            CALL NSFLUX_2D(WL,WR,W1P,W2P,DWT0,DSL,DSR,FV,RMU0,GAMMA,PRN)
!!$            FV = FV*DY

            DO N = 1,4
! --- ACCUMULATE COMPLETE CONVECTIVE FLUX
               DW(N,I,J)   = DW(N,I,J)   + FC(N)
               DW(N,I+1,J) = DW(N,I+1,J) - FC(N)
!!$
!!$! --- ACCUMULATE JUST THE ARTIFICIAL DISSIPATION PORTION
!!$               FW(N,I,J)   = FW(N,I,J)   + FD(N)
!!$               FW(N,I+1,J) = FW(N,I+1,J) - FD(N)
!!$
!!$! --- ACCUMULATE THE VISCOUS DISSIPATION
!!$               VW(N,I,J)   = VW(N,I,J)   + FV(N)
!!$               VW(N,I+1,J) = VW(N,I+1,J) - FV(N)
            ENDDO
         ENDDO
      ENDDO


!------------------------------------------------------------
!
! --- FIND FLUX CONTRIBUTIONS IN J DIRECTION
!
!------------------------------------------------------------

      DO J = 2,JL -1
         DO I = 2,IL -1
! ---   LENGTH OF FACE J + 1/2
            DX   = X(1,I+1,J) - X(1,I,J)

            DSL  = .5*(X(2,I,J+1) - X(2,I,J))
            DSR  = .5*(X(2,I,J+2) - X(2,I,J+1))

!----------WENO interpolation
            DO N = 1,4

! RECONSTRUCT WL AND WR FOR N=1,4 FROM W - YDIRECTION J INDEX

               PLW(1) = 0.5*W(N,I,J)+0.5*W(N,I,J+1)
               PLW(2) = -0.5*W(N,I,J-1)+1.5*W(N,I,J)

               PRW(1) = 1.5*W(N,I,J+1)-0.5*W(N,I,J+2)
               PRW(2) = 0.5*W(N,I,J)+0.5*W(N,I,J+1)

               betaL(1) = (W(N,I,J+1)-W(N,I,J))**2
               betaL(2) = (W(N,I,J)-W(N,I,J-1))**2
            
               betaR(1) = (W(N,I,J+2)-W(N,I,J+1))**2
               betaR(2) = (W(N,I,J+1)-W(N,I,J))**2

               alphaL(1) = 2.0/(3.0*(SMV+betaL(1))**2)
               alphaL(2) = 1.0/(3.0*(SMV+betaL(2))**2)

               alphaR(1) = 1.0/(3.0*(SMV+betaR(1))**2)
               alphaR(2) = 2.0/(3.0*(SMV+betaR(2))**2)

               omegaL(1) = alphaL(1)/(alphaL(1)+alphaL(2))
               omegaL(2) = alphaL(2)/(alphaL(1)+alphaL(2))

               omegaR(1) = alphaR(1)/(alphaR(1)+alphaR(2))
               omegaR(2) = alphaR(2)/(alphaR(1)+alphaR(2))

               WL(N) = omegaL(1)*PLW(1)+omegaL(2)*PLW(2)
               WR(N) = omegaR(1)*PRW(1)+omegaR(2)*PRW(2)

            ENDDO
		IF((J == 2).OR.(J==JL-1)) THEN
               	WL(1:4) = .5*(W(1:4,I,J) + W(1:4,I,J+1))
               	WR      =     WL
            	ENDIF
! -- ROTATE THE VARIABLES
            WLP(1) = WL(1)
            WLP(2) = WL(3)
            WLP(3) =-WL(2)
            WLP(4) = WL(4)

            WRP(1) = WR(1)
            WRP(2) = WR(3)
            WRP(3) =-WR(2)
            WRP(4) = WR(4)

            W1P(1) = W(1,I,J)
            W1P(2) = W(3,I,J)
            W1P(3) =-W(2,I,J)
            W1P(4) = W(4,I,J)

            W2P(1) = W(1,I,J+1)
            W2P(2) = W(3,I,J+1)
            W2P(3) =-W(2,I,J+1)
            W2P(4) = W(4,I,J+1)

! --- COMPUTE CONVECTIVE FLUXES
      CALL BGKFLUX(WLP,WRP,W1P,W2P,DSL,DSR,DTMIN,RMU0,FCP,GAMMA,PRN)
! --- ROTATE THE FLUXES BACK
            FC(1)  = FCP(1)*DX
            FC(2)  =-FCP(3)*DX
            FC(3)  = FCP(2)*DX
            FC(4)  = FCP(4)*DX
!!$
!!$            FD(1)  = FDP(1)*DX
!!$            FD(2)  =-FDP(3)*DX
!!$            FD(3)  = FDP(2)*DX
!!$            FD(4)  = FDP(4)*DX
!!$! --- COMPUTE VISCOUS FLUXES
!!$!            MU = 0.5*(RLV(I,J)+ RLV(I,J+1))
!!$            CALL NSFLUX_2D(WLP,WRP,W1P,W2P,DWT0,DSL,DSR,FVP,RMU0,GAMMA,PRN)
!!$            FV(1)  = FVP(1)*DX
!!$            FV(2)  =-FVP(3)*DX
!!$            FV(3)  = FVP(2)*DX
!!$            FV(4)  = FVP(4)*DX

            DO N = 1,4
! --- ACCUMULATE COMPLETE CONVECTIVE FLUX
               DW(N,I,J)   = DW(N,I,J)   + FC(N)
               DW(N,I,J+1) = DW(N,I,J+1) - FC(N)
!!$               
!!$! --- ACCUMULATE JUST THE ARTIFICIAL DISSIPATION PORTION
!!$               FW(N,I,J)   = FW(N,I,J)   + FD(N)
!!$               FW(N,I,J+1) = FW(N,I,J+1) - FD(N)
!!$
!!$! --- ACCUMULATE THE VISCOUS DISSIPATION
!!$               VW(N,I,J)   = VW(N,I,J)   + FV(N)
!!$               VW(N,I,J+1) = VW(N,I,J+1) - FV(N)
            ENDDO
         ENDDO
      ENDDO

!!$
!!$!---------------------------------------------------------------------!
!!$!                                                                     !
!!$! --- COMPARE DIFFUSIVE AND VISCOUS FLUXES AND ADD GREATER OF THE TWO !
!!$!                                                                     !
!!$!---------------------------------------------------------------------!
!!$
!!$      DO J = 3,JL-1
!!$         DO I = 3,IL-1
!!$            DO N = 2,4
!!$               IF( ABS(VW(N,I,J)) > ABS(FW(N,I,J))) THEN
!!$                  DW(N,I,J) = DW(N,I,J) - FW(N,I,J) + VW(N,I,J)
!!$               ENDIF
!!$            ENDDO
!!$         ENDDO
!!$      ENDDO

      END SUBROUTINE DERIVS_BGK
