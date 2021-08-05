!======================================================================================================!
! Progetto II
! Equazione dei flussi a superficie libera su griglie strutturate
!======================================================================================================!
! UniversitÓ degli studi di Ferrara - Dipartimento di Matematica e Informatica
! Corso di Algoritmi per il Calcolo Parallelo
! Prof. Walter Boscheri
! A.A. 2020/2021
!======================================================================================================!
! Studenti.
! Alberto Paparella - 144261
! Elena Zerbin      - 145302
!======================================================================================================!
! IO.f90
! Usefull subroutines to print results (can be used for plotting)
!======================================================================================================!   
!#define PARALLEL  
    SUBROUTINE DataOutput(timestep,istart, iend, myrank,cmyrank)
    !==================================================================================================!
    USE VarDef_mod
    IMPLICIT NONE
    !==================================================================================================!
    INTEGER, INTENT(IN) :: timestep, istart, iend, myrank
    !==================================================================================================!
    INTEGER             :: i, j, DataUnit
    REAL                :: ub, vb
    CHARACTER(LEN=10)   :: citer,cmyrank
    CHARACTER(LEN=200)  :: IOFileName
    !==================================================================================================!
    WRITE(citer,'(I8.8)') timestep 
#ifdef PARALLEL  
    WRITE(citer,'(I8.8)') timestep
    WRITE(cmyrank,'(I4.4)') myrank                        ! convert iteration number to string
    IOFileName = TRIM(TestName)//'-'//TRIM(citer)//'-'//TRIM(cmyrank)//'.dat' ! name of output file
    DataUnit   = 100+myrank                               ! unit for output file
    !
#endif
    
    ! Convert iteration number to string
    IOFileName = TRIM(TestName)//'-'//TRIM(citer)//'.dat' ! Name of output file
    DataUnit   = 100                                      ! Unit for output file
    !==================================================================================================!
    OPEN(UNIT=DataUnit, FILE=TRIM(IOFilename), STATUS='UNKNOWN', ACTION='WRITE')
    !==================================================================================================! 
#ifdef MATLAB_OUT
    ! Header
    WRITE(DataUnit,*) IMAX
    WRITE(DataUnit,*) JMAX
    ! Coordinates
    ! Note: these are 2 vectors, they will be the coordinates of a matrix (eta)
  ! #ifdef PARALLEL 
        WRITE(DataUnit,*) IMAX
        WRITE(DataUnit,*) JMAX
        ! Coordinates
        DO i = istart, iend
          WRITE(DataUnit,*) x(i)
        ENDDO  
      !
         DO i = istart, iend
            DO j= istart, iend
                WRITE(DataUnit,*) eta(i,j)
            ENDDO
         ENDDO
     
           ! Velocity (interpolation at barycenters)
        ! Velocity on the x axys
        DO i = istart, iend
            DO j= istart, iend
                ub = 0.5 * ( u(i,j) + u(i+1,j) )  
                WRITE(DataUnit,*) ub
            ENDDO
        ENDDO
        ! Velocity on the y axys
        DO i = istart, iend
            DO j= istart, iend
                vb = 0.5 * ( v(i,j) + v(i,j+1) )  
                WRITE(DataUnit,*) vb
            ENDDO
        ENDDO
    !#endif       
    DO i = 1, IMAX
        WRITE(DataUnit,*) xb(i)
    ENDDO  
    DO j = 1, JMAX
        WRITE(DataUnit,*) yb(j)
    ENDDO  
    ! Pressure
    ! For each x eta coordinate, i print JMAX values for y eta coordinate
    DO i = 1, IMAX
        DO j= 1, JMAX     
            WRITE(DataUnit,*) eta(i,j)
        ENDDO
    ENDDO    
    ! Velocity (interpolation at barycenters)
    ! Velocity on the x axys
    DO i = 1, IMAX
        DO j = 1, JMAX
            ub = 0.5 * ( u(i,j) + u(i+1,j) )  
            WRITE(DataUnit,*) ub
        ENDDO
    ENDDO
    ! Velocity on the y axys
    DO i = 1, IMAX
        DO j = 1, JMAX
            vb = 0.5 * ( v(i,j) + v(i,j+1) )  
            WRITE(DataUnit,*) vb
        ENDDO
    ENDDO
#else
    ! Current time 
    WRITE(DataUnit,*) 'TITLE = "CURRENT TIME ', time, ' "'   
    ! Variables
    WRITE(DataUnit,*) 'VARIABLES = "x" "y" "eta" "u" "v" '
    ! Header
    WRITE(DataUnit,*) 'ZONE T="Only Zone", I=', IMAX, ' J=', JMAX, ' F=POINT'
    !==================================================================================================!
    DO i = 1, IMAX
      DO j = 1, JMAX
          ub = 0.5 * ( u(i,j) + u(i+1,j) )  ! Interpolate x velocity at barycenters
          vb = 0.5 * ( u(i,j) + u(i,j+1) )  ! Interpolate y velocity at barycenters
          WRITE(DataUnit,*) xb(i), yb(j), eta(i,j), ub, vb
      ENDDO
    ENDDO  
#endif    
    !==================================================================================================!
    CLOSE(DataUnit)
!======================================================================================================!
    END SUBROUTINE DataOutput  
 

  