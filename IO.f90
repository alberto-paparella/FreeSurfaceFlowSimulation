SUBROUTINE DataOutput(timestep)
    !------------------------------------------------------------!
    USE VarDef_mod
    IMPLICIT NONE
    !------------------------------------------------------------!
    INTEGER, INTENT(IN) :: timestep
    !
    INTEGER             :: i, DataUnit
    REAL                :: ub
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
#ifdef MATLAB_OUT
    ! Header
    WRITE(DataUnit,*) IMAX
    ! Coordinates
    DO i = 1, IMAX
      WRITE(DataUnit,*) xb(i)
    ENDDO  
    ! Pressure
    DO i = 1, IMAX
      WRITE(DataUnit,*) eta(i)
    ENDDO
    ! Velocity (interpolation at barycenters)
    DO i = 1, IMAX
      ub = 0.5 * ( u(i) + u(i+1) )  
      WRITE(DataUnit,*) ub
    ENDDO
#else
    ! Current time 
    WRITE(DataUnit,*) 'TITLE = "CURRENT TIME ', time, ' "'   
    ! Variables
    WRITE(DataUnit,*) 'VARIABLES = "x" "eta" "u" '
    ! Header
    WRITE(DataUnit,*) 'ZONE T="Only Zone", I=', IMAX, ' F=POINT'
    !
    DO i = 1, IMAX
      ub = 0.5 * ( u(i) + u(i+1) )   ! interpolate velocity at barycenters
      WRITE(DataUnit,*) xb(i), eta(i), ub
    ENDDO  
#endif    
    !
    CLOSE(DataUnit)
    !
END SUBROUTINE DataOutput  
  