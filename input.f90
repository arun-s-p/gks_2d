!
!     ******************************************************************
!
      SUBROUTINE INPUT

!     ******************************************************************
!

      USE DIMS
      USE FLO_VAR
      USE MESH_VAR           ,ONLY: RLEN
      USE FLO_PARAM
      USE TIME_VAR           ,ONLY: NCYC
      USE IODEFS             ,ONLY : IREAD,FNAME0,FNAME1,FNAME2

      IMPLICIT NONE

!     ******************************************** 
!     READ GRID LIMITS
!     ******************************************** 

      NX = 128
      NY = 151

      IL = NX + 2
      JL = NY + 2

      IE = NX + 3
      IB = NX + 4

      JE = NY + 3
      JB = NY + 4
!     ******************************************** 
!     CHOICE OF SCHEME
!     ******************************************** 
!!$      SELECTED_SCHEME = UPS_SCHEME
!!$      SELECTED_SCHEME = CUSP_SCHEME
      SELECTED_SCHEME = BGK_SCHEME

!     ******************************************** 
!     READ FLOW PARAMETERS
!     ******************************************** 

      KVIS  = 1    ! CHANGE THIS TO READ STATEMENT      
      GAMMA = 1.667 !1.4  ! CHANGE THIS TO READ STATEMENT      
      PRN   = 1.0  ! CHANGE THIS TO READ STATEMENT    


      RHO0  = 1  ! CHANGE THIS TO READ STATEMENT
      RMU0  =  2.2360612501607223E-004!1.0E-5
      P0     = 10.
      C0     = SQRT(GAMMA*P0/RHO0)


      RM    = 0.2  ! CHANGE THIS TO READ STATEMENT
      U0     = 2.*RM*C0
      V0     = 0.

!     REYNOLDS NUMBER
      RE    = 200.
      RLEN  = 1./14.  ! CHANGE THIS TO READ STATEMENT
      VTHICK0 = 2.*RE*RMU0/(U0*RHO0)
      MTHICK0 = 4.*VTHICK0
      
!     ******************************************** 
!     READ RUN PARAMETERS
!     ******************************************** 

!      NCYC =31000      ! CHANGE THIS TO READ STATEMENT
      NCYC = 1000
      CFL  = .6         ! CHANGE THIS TO READ STATEMENT

!     ******************************************** 
!     READ FILE NAMES
!     ******************************************** 
      IREAD = 0                          ! IREAD = 1 FOR RESTART
      WRITE(FNAME0,"('MLAYER.rst')")     !RESTART FILE
      WRITE(FNAME1,"('MLAYER.dat')")     !TECPLOT OUTPUT
      WRITE(FNAME2,"('growth_rate.dat')") !VARIABLES IN SIMILARITY COORDS

      END SUBROUTINE INPUT
