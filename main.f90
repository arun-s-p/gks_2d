      PROGRAM MIXINGLAYER

!
! --- ******************************************************************
!     *                                                                *
!     * TRIES TO SIMULATE A MIXING LAYER VIA 'MILES'                   *
!     *                                                                *
! --- ******************************************************************
!

      USE DIMS
      USE TIME_VAR
      USE IODEFS
      USE MESH_VAR    ,ONLY :RLEN
      USE FLO_PARAM   !,ONLY :U0
      USE FLO_VAR
      USE MESH_VAR
! --- ******************************************************************

      IMPLICIT NONE

      INTEGER   :: NSTART, ct, jt
      CHARACTER*20 :: FNAME
      DOUBLE PRECISION 	:: mom,tstop,tscale
! --- ******************************************************************

      CALL INPUT
      CALL MESH


      CALL INIT
	PRINT *,"U0 = ",U0,"Vt = ",VTHICK0,"RE = ",RE,"M = ",RM

      counter=0
      viztime=0.
      vtime=1.
      tstop=20.

!CALL OUTPUT
      counter=counter+1
      IF (IREAD==0) THEN
         TIME = 0.
         NSTART = 100
! ---- DO NSTART CYCLES TO GET A GOOD ESTIMATE OF INITIAL PRESSURE
         CALL TSTEP
         DTL   = DTMIN/FLOAT(NSTART)
	 PRINT *, DTMIN
         DTMIN = DTMIN/FLOAT(NSTART)
         DO CYC = 1,NSTART
            CALL UPDATE
            TIME = TIME + DTMIN
         ENDDO
      ENDIF


viztime = viztime+vtime
! ---- ACTUAL TIME STEPPING STARTS HERE
      DO CYC = 1,NCYC
         CALL TSTEP
         CALL UPDATE
!	 PRINT *, TIME
         TIME = TIME + DTMIN
         IF(MOD(CYC,50) == 0) PRINT *, CYC, DTMIN, TIME
!         IF(TIME/(RLEN/U0) >= 100) THEN
!            CALL OUTPUT
!            STOP
!         ENDIF
	tscale=TIME*0.5*U0/VTHICK0
	IF (tscale>viztime) THEN
		CALL OUTPUT
		PRINT *,'----------writing---------- t = ', TIME
		counter=counter+1
		viztime=viztime+vtime
	ENDIF
	IF (tscale>tstop) STOP
      ENDDO

   
	
      CALL OUTPUT

     END PROGRAM MIXINGLAYER
 
