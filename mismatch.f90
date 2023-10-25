PROGRAM MAIN
  
  IMPLICIT NONE

  INTEGER,PARAMETER        :: dp = selected_real_kind(p=15,r=300), repeats=12
  REAL(KIND=dp),PARAMETER  :: pi = 4.0_dp*ATAN(1.0_dp)
  REAL, PARAMETER          :: mm = pi/3

  REAL                     :: qH2O_O, qH2O_H, qN, qNH, idipolez, mdipolez, angle
  REAL, DIMENSION(3)       :: Cpi, Cpm, Hpi, Hpm, CHi, CHm
  CHARACTER(LEN=500)       :: initial, minimised, dump
  INTEGER                  :: ifnum, mfnum, atoms, i, j, k, n
  
  qH2O_O = -1.1794
  qN     = -0.24
  qH2O_H = 0.5897
  qNH    = 0.31

  CALL GET_COMMAND_ARGUMENT(1, initial) 
  CALL GET_COMMAND_ARGUMENT(2, minimised) 

  ifnum=77
  mfnum=53

  idipolez = 0
  mdipolez = 0
  atoms = 384

  OPEN(ifnum, file = TRIM(initial), status = 'old')  
  READ(ifnum,*) ; READ(ifnum, *) dump

  OPEN(mfnum, file = TRIM(minimised), status = 'old')  
  READ(mfnum,*) ; READ(mfnum, *) dump

  n = 0

  DO i = 1, repeats
     READ(ifnum, *) dump, dump, dump, Cpi(1), Cpi(2), Cpi(3)
     idipolez = idipolez + Cpi(3)*qN
     READ(mfnum, *) dump, dump, dump, Cpm(1), Cpm(2), Cpm(3)
     mdipolez = mdipolez + Cpm(3)*qN
     n = n + 1
     DO j = 1, 4
        READ(ifnum, *) dump, dump, dump, Hpi(1), Hpi(2), Hpi(3)
        IF (Hpi(3) .GT. 25) THEN
           idipolez = idipolez + (Hpi(3)-25.5)*qNH
        ELSE
           idipolez = idipolez + Hpi(3)*qNH
        END IF
        READ(mfnum, *) dump, dump, dump, Hpm(1), Hpm(2), Hpm(3)
        IF (Hpm(3) .GT. 25) THEN
           mdipolez = mdipolez + (Hpm(3)-25.5)*qNH
        ELSE
           mdipolez = mdipolez + Hpm(3)*qNH
        END IF
        n = n + 1
        CHi = Hpi - Cpi  ; CHi = CHi/SQRT(DOT_PRODUCT(CHi, CHi))
        CHm = Hpm - Cpm  ; CHm = CHm/SQRT(DOT_PRODUCT(CHm, CHm))
        angle = ACOS(DOT_PRODUCT(CHi, CHm))
        IF (angle .GT. mm)  WRITE(*, *) TRIM(minimised), " ammonium ", n, " rotated by ", angle*180/pi
     END DO
  END DO

  DO i = repeats+1, atoms
     READ(ifnum, *) dump, dump, dump, Cpi(1), Cpi(2), Cpi(3)
     READ(mfnum, *) dump, dump, dump, Cpm(1), Cpm(2), Cpm(3)
     n = n + 1 
     DO j = 1, 2
        READ(ifnum, *) dump, dump, dump, Hpi(1), Hpi(2), Hpi(3)
        IF (Hpi(3) .GT. 25) THEN
           idipolez = idipolez + (Hpi(3)-25.5)*qH2O_H
        ELSE
           idipolez = idipolez + Hpi(3)*qH2O_H
        END IF
        READ(mfnum, *) dump, dump, dump, Hpm(1), Hpm(2), Hpm(3)
        IF (Hpm(3) .GT. 25) THEN
           mdipolez = mdipolez + (Hpm(3)-25.5)*qH2O_H
        ELSE
           mdipolez = mdipolez + Hpm(3)*qH2O_H
        END IF
        n = n + 1
        CHi = Hpi - Cpi  ; CHi = CHi/SQRT(DOT_PRODUCT(CHi, CHi))
        CHm = Hpm - Cpm  ; CHm = CHm/SQRT(DOT_PRODUCT(CHm, CHm))
        angle = ACOS(DOT_PRODUCT(CHi, CHm))
        IF (angle .GT. mm) WRITE(*, *) TRIM(minimised), " water ", n, " rotated by ", angle*180/pi
     END DO
     READ(ifnum, *) dump, dump, dump, Hpi(1), Hpi(2), Hpi(3)
     IF (Hpi(3) .GT. 25) THEN
        idipolez = idipolez + (Hpi(3)-25.5)*qH2O_O
     ELSE
        idipolez = idipolez + Hpi(3)*qH2O_O
     END IF
     READ(mfnum, *) dump, dump, dump, Hpm(1), Hpm(2), Hpm(3)
     IF (Hpm(3) .GT. 25) THEN
        mdipolez = mdipolez + (Hpm(3)-25.5)*qH2O_O
     ELSE
        mdipolez = mdipolez + Hpm(3)*qH2O_O
     END IF
     n = n + 1
     !CHi = Hpi - Cpi  ; CHi = CHi/SQRT(DOT_PRODUCT(CHi, CHi))
     !CHm = Hpm - Cpm  ; CHm = CHm/SQRT(DOT_PRODUCT(CHm, CHm))
     !angle = ACOS(DOT_PRODUCT(CHi, CHm))
     !IF (angle .GT. mm) WRITE(*, *) TRIM(minimised), " Msite on ", n, " rotated by ", angle*180/pi
  END DO

  READ(ifnum, *) dump, dump, dump 

  WRITE(*,*) TRIM(initial), idipolez
  WRITE(*,*) TRIM(minimised), mdipolez



END PROGRAM MAIN
