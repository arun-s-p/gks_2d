!
!     ******************************************************************
!
      MODULE FLO_PARAM
!
!     ******************************************************************
!     *                                                                *
!     *   FLOW PARAMETERS                                              *
!     *                                                                *
!     ******************************************************************
!
      REAL      :: GAMMA,RM,RHO0,P0,C0,U0,V0,HH0,EI0,PTOT,CVTTOT
      REAL      :: RE,PRN,PRT,T0,RMU0
      REAL      :: CFL,PWENO
      INTEGER   :: KVIS
      INTEGER   :: SELECTED_SCHEME,WENO_FLAG
      INTEGER, PARAMETER :: UPS_SCHEME = 0,BGK_SCHEME = 1,CUSP_SCHEME = 2
      REAL , DIMENSION(4) :: gammaw
      REAL		:: MTHICK0,VTHICK0
      END MODULE FLO_PARAM
