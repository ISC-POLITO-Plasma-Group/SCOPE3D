       INTEGER FUNCTION LOCALDIM(GDIM,NP,ME)

!   GDIM = GLOBAL DIMENSION OF AN ARRAY TO BE SHEARED AMONG PROCESSORS
!   NP   = NUMBER OF PROCESSOR
!   ME   = INDEX OF THE CALLING PROCESSOR (starting from 0)
!  
!   THIS FUNCTION RETURN THE NUMBER OF ELEMENTS OF THE ARRAY STORED
!   IN THE LOCAL MEMORY OF THE PROCESSOR "ME"

       IMPLICIT NONE
       INTEGER GDIM,NP,ME, R,Q

       Q = INT(GDIM/NP)
       R = MOD(GDIM,NP)

       IF((ME+1).LE.R) THEN
         LOCALDIM = Q+1
       ELSE
         LOCALDIM = Q
       END IF
 
       RETURN
       END