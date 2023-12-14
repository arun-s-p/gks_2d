!
!     ******************************************************************
!
      SUBROUTINE UPDATE

!     ******************************************************************
!

!
! --- ******************************************************************
!     *                                                                *
!     * UPDATE FLOW SOLUTION VIA SOME RK TYPE SCHEME                   *
!     *                                                                *
! --- ******************************************************************
!

      USE DIMS              ,ONLY : IL,JL,NX
      USE MESH_VAR          ,ONLY : VOL,RLEN,XX,YY,X
      USE TIME_VAR          ,ONLY : DTL,TIME,DTMIN,CYC
      USE FLO_VAR           ,ONLY : W
      USE SOLV_VAR          ,ONLY : W0,DW
      USE FLO_PARAM         ,ONLY : GAMMA,SELECTED_SCHEME,BGK_SCHEME,UPS_SCHEME, &
                                    U0

! --- ******************************************************************
 
      IMPLICIT NONE
 
! --- ******************************************************************


! --- ******************************************************************
!
!          LOCAL VARIABLES
!
! --- ******************************************************************
      REAL    :: DU,DV,DP
      REAL    :: DT
      REAL    :: EL,RTMAX,RTRMS
      REAL , DIMENSION(4) :: WC
      REAL    :: GM1
      REAL    :: EDDY_T

!      REAL    :: DU,DV
      REAL    :: DX,DY,DXG
      REAL    :: MEANOMEGA, MAXOMEGA, OMEGA, RELVOR
      
      INTEGER :: I,J,N
      INTEGER :: IMAX,JMAX
! --- ******************************************************************

      GM1  = GAMMA - 1.D0
      W0 = W
      CALL DERIVS_BGK

      RTMAX  = 0.
      RTRMS  = 0.
      
      IMAX = 3; JMAX = 3
      DO J = 3, JL - 1
         DO I = 3, IL - 1
            EL   = W0(2,I,J)
            DT   = DTL(I,J)/VOL(I,J)

! --- COMPUTE CONSERVATIVE VARIABLES FROM PRIMITIVE VARIABLES

            WC(1) = W0(1,I,J)
            WC(2) = W0(1,I,J)*W0(2,I,J)
            WC(3) = W0(1,I,J)*W0(3,I,J)
            WC(4) = W0(4,I,J)/GM1 + .5*W0(1,I,J)*(W0(2,I,J)**2 + W0(3,I,J)**2)

! --- UPDATE CONSERVATIVE VARIABLES
            DO N = 1,4
               WC(N) = WC(N) - DT*DW(N,I,J)
            ENDDO
!-----GRAVITY
!		WC(3) = WC(3) - DTL(I,J)*W0(1,I,J)*9.81
!		WC(4) = WC(4) - DTL(I,J)*W0(1,I,J)*W0(3,I,J)*9.81

! --- RECOMPUTE PRIMITIVE VARIABLES FROM CONSERVATIVE VARIABLES
            W(1,I,J)  = WC(1)
            W(2,I,J)  = WC(2)/WC(1)
            W(3,I,J)  = WC(3)/WC(1)
            W(4,I,J)  = GM1*(WC(4) - .5*(WC(2)**2 + WC(3)**2)/WC(1))

            RTRMS = RTRMS + ( W(2,I,J) - EL )**2
            IF(ABS(W(2,I,J)-EL) > RTMAX) THEN
               RTMAX = MAX(RTMAX, ABS(W(2,I,J) - EL ) )
               IMAX = I
               JMAX = J
            ENDIF
         ENDDO
      ENDDO
   

      CALL BCONDS

      RTRMS    = SQRT(RTRMS   /((IL-3)*(JL-3)))

      IF(MOD(CYC,50) == 0) THEN

!  ----------- VORTICITY THICKNESS ------------

         DXG   = 2./(FLOAT(NX) - 1.);
         MAXOMEGA  = 0.
         
         DO J = 3,JL-1
            MEANOMEGA = 0.
            DO I = 3,IL-1
!   ----------- COMPUTE VORTICITY ---------- !
               DX = X(1,I+1,J) - X(1,I-1,J)
               DY = X(2,I,J+1) - X(2,I,J-1)
               
               DU = W(2,I,J+1) - W(2,I,J-1)
               DV = W(3,I+1,J) - W(3,I-1,J)
               OMEGA = DV/DX - DU/DY
               
               MEANOMEGA = MEANOMEGA + .5*OMEGA*DXG
               
            END DO
            MAXOMEGA = MAX(MAXOMEGA,ABS(MEANOMEGA))
         END DO
         
         EDDY_T = TIME/(RLEN/U0)
         RELVOR = 28.D0/MAXOMEGA
         !WRITE(33,*)  CYC, EDDY_T, TIME, DTMIN,RELVOR !, RTRMS,RTMAX
         !CALL FLUSH(33)
      ENDIF

      END SUBROUTINE UPDATE
