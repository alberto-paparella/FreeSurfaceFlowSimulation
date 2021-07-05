SUBROUTINE DataOutput(timestep)
    !------------------------------------------------------------!
    USE VarDef_mod
    IMPLICIT NONE
    !------------------------------------------------------------!
    INTEGER, INTENT(IN) :: timestep
    !
    INTEGER             :: i, DataUnit
    INTEGER             :: j
    REAL                :: ub
    REAL                :: vb
    CHARACTER(LEN=10)   :: citer
    CHARACTER(LEN=200)  :: IOFileName
    !------------------------------------------------------------!
    !
    WRITE(citer,'(I4.4)') timestep                        ! convert iteration number to string
    IOFileName = TRIM(TestName)//'-'//TRIM(citer)//'.dat' ! name of output file
    DataUnit   = 100                                      ! unit for output file
    !
    OPEN(UNIT=DataUnit, FILE=TRIM(IOFilename), STATUS='UNKNOWN', ACTION='WRITE')
    !
! WARNING - At the moment, we do not consider this case (matlab output)
#ifdef MATLAB_OUT
    ! Header
    WRITE(DataUnit,*) IMAX
    ! Coordinates
    DO i = 1, IMAX
      WRITE(DataUnit,*) xb(i)
    ENDDO  
      DO j = 1, JMAX
      WRITE(DataUnit,*) yb(j)
    ENDDO  
    ! Pressure
    DO j = 1, JMAX
        DO i= 1,IMAX     
            WRITE(DataUnit,*) eta(i,j)
        ENDDO
    ENDDO
    
#else
    ! Current time 
    WRITE(DataUnit,*) 'TITLE = "CURRENT TIME ', time, ' "'   
    ! Variables
    ! WRITE(DataUnit,*) 'VARIABLES = "x" "eta" "u" '
    WRITE(DataUnit,*) 'VARIABLES = "x" "y" "eta"  '  ! 2D
    ! Header
    ! WRITE(DataUnit,*) 'ZONE T="Only Zone", I=', IMAX, ' F=POINT'
    WRITE(DataUnit,*) 'ZONE T="Only Zone", I=', IMAX, ' J=', JMAX, ' F=POINT' ! 2D
    !
    DO i = 1, IMAX
      DO j = 1, JMAX    ! 2D
          ! TODO Fix ub fot 2D versione: how to update ub?
          ! ub = 0.5 * ( u(i) + u(i+1) )   ! interpolate velocity at barycenters
         ! ub = 0.5 * ( u(i,j) + u(i+1, j+1) )   ! 2D: interpolate velocity at barycenters
          !vb = 0.0  ! 2D: TODO interpolate velocity at barycenters  for v
          WRITE(DataUnit,*) xb(i), yb(j), eta(i,j)
      ENDDO ! 2D
    ENDDO  
#endif    
    !
    CLOSE(DataUnit)
    !
END SUBROUTINE DataOutput  
  