!-------------------------------------------------------------------------------------------------------!
!													!
!  Add the following lines in the output subroutine to convert to binary and thereupon remove the ASCII	!
!													!
!      	WRITE(COMMAND,'(A,A)') 'preplot ',TRIM(FNAME1)							!
!	CALL SYSTEM(COMMAND)										!
!     	WRITE(COMMAND,'(A,A)') 'rm ',TRIM(FNAME1)							!
!	CALL SYSTEM(COMMAND)										!
!-------------------------------------------------------------------------------------------------------!




      PROGRAM OUTPUT

	INTEGER 		:: I
      CHARACTER(LEN=30) 	:: jobname = 'kh'
      CHARACTER(LEN=30) 	:: FNAME1,COMMAND

	DO I=0,50
      WRITE(FNAME1,'(A,I4.4,A)') TRIM(jobname), I, '.tec'

      WRITE(COMMAND,'(A,A)') 'preplot ',TRIM(FNAME1)
	CALL SYSTEM(COMMAND)
     WRITE(COMMAND,'(A,A)') 'rm ',TRIM(FNAME1)
	CALL SYSTEM(COMMAND)

	END DO
   12 FORMAT(1X,9F18.10)
   13 FORMAT(1X,2I5,1X,9F18.10)


      END PROGRAM OUTPUT
