!
!     ******************************************************************
!
      SUBROUTINE OUTPUT

!     ******************************************************************
!

      USE DIMS
      USE FLO_VAR
      USE MESH_VAR
      USE FLO_PARAM
      USE IODEFS
      USE TIME_VAR  !,ONLY: TIME,CYC

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

      REAL :: XC,YC
      REAL :: DU,DV
      REAL :: OMEGA
      REAL :: EDDY_T,UTHICK,MT_T
      REAL :: URA,RRA,UFAVRE,MTHICK
      REAL :: VTHICK,DUDY
      REAL, DIMENSION(JL) :: UAVG
      INTEGER :: I,J,N,IR

      CHARACTER(LEN=30) 	:: COMMAND,jobname = 'kh'

!     ******************************************************************

      EDDY_T = TIME/(RLEN/U0)

      IFLO = 18
      ISIM = 28
      WRITE(FNAME1,'(A,I4.4,A)') TRIM(jobname), counter, '.tec'
      OPEN  (IOUT,FILE=FNAME0,ACCESS='SEQUENTIAL')
      OPEN  (IFLO,FILE=FNAME1,status='unknown')
      
! --------------------------------
!     DUMP SOLUTION FOR RESTARTING
! --------------------------------


      DO J = 1,JE
         DO I = 1,IE
            N = INT(EDDY_T) + RE
            WRITE(IOUT,'(4(1X,E18.10))') W(1,I,J),W(2,I,J),W(3,I,J),W(4,I,J)
         ENDDO
      ENDDO
      WRITE(IOUT,*) TIME
      CLOSE (IOUT)

! -------------------------------------
!    WRITE SOLUTION IN TECPLOT FORMAT
! -------------------------------------
      WRITE (IFLO,'(A)') 'TITLE = BLAYER NEW SCHEME'

      WRITE (IFLO,'(A)')                               &
      'VARIABLES = "X" "Y" "RHO" "U" "V" "PRESSURE" "VORTICITY"'

      WRITE ( IFLO, '(A,I6,A,I6,A)' ) 'Zone I=',NX-1, ', J=', NY-1, &
             ', F=POINT'

       DO J = 3,JL-1
         DO I = 3,IL-1
            XC   = .5*(XX(I) + XX(I+1))
            YC   = .5*(YY(J) + YY(J+1))
!   ----------- COMPUTE VORTICITY ---------- !
            DX = X(1,I+1,J) - X(1,I-1,J)
            DY = X(2,I,J+1) - X(2,I,J-1)

            DU = W(2,I,J+1) - W(2,I,J-1)
            DV = W(3,I+1,J) - W(3,I-1,J)
            OMEGA = DV/DX - DU/DY

            WRITE(IFLO,12) XC,YC,W(1,I,J),W(2,I,J),W(3,I,J),W(4,I,J),OMEGA
         END DO
      END DO
      CLOSE (IFLO)

   12 FORMAT(1X,9F18.10)
   13 FORMAT(1X,2I5,1X,9F18.10)

!      WRITE(COMMAND,'(A,A)') 'preplot ',TRIM(FNAME1)
!	CALL SYSTEM(COMMAND)
!     WRITE(COMMAND,'(A,A)') 'rm ',TRIM(FNAME1)
!	CALL SYSTEM(COMMAND)

!----write momentum thickness----
	IR=IE/2
	MT_T = TIME*0.5*U0/VTHICK0
	!IF (MT_T .gt. 600.) THEN
	OPEN  (ISIM,FILE=FNAME2,ACCESS='APPEND')
	MTHICK = 0.
	
      DO J = 3,JL-1
	URA = 0.
	RRA = 0.
	DO I = 3,IL-1
	   URA = URA+(W(1,I,J)*W(2,I,J))
	   RRA = RRA+W(1,I,J)
	END DO
	URA = URA/(IL-3)
	RRA = RRA/(IL-3)
	UFAVRE = URA/RRA
	UTHICK = (URA/U0)**2
         MTHICK=MTHICK+RRA*(.25-UTHICK)*DY
	UAVG(J) = UFAVRE
      ENDDO

      VTHICK = 0.
      DO J = 4,JL-2
	DUDY = ABS((UAVG(J+1)-UAVG(J-1))/DY)
	IF (DUDY .gt. VTHICK) VTHICK = DUDY
      ENDDO
	MTHICK = 2.*MTHICK/MTHICK0

	VTHICK = U0/VTHICK
!print *, vthick
	IF (CYC .eq. 0) RE = VTHICK0*U0*RHO0/(2.*RMU0)
        VTHICK = VTHICK/VTHICK0
	 WRITE(ISIM,'(4(1X,E18.10))') MT_T,MTHICK,VTHICK,URA
	CLOSE (ISIM)
	!END IF
      END SUBROUTINE OUTPUT
