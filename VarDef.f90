!======================================================================================================!
! Progetto II
! Equazione dei flussi a superficie libera su griglie strutturate
!======================================================================================================!
! Università degli studi di Ferrara - Dipartimento di Matematica e Informatica
! Corso di Algoritmi per il Calcolo Parallelo
! Prof. Walter Boscheri
! A.A. 2020/2021
!======================================================================================================!
! Studenti.
! Alberto Paparella - 144261
! Elena Zerbin      - 145302
!======================================================================================================!
! VarDef.f90
! In this file we defined the variables we are going to use
!======================================================================================================!
!#define PARALLEL
MODULE VarDef_mod
    IMPLICIT NONE
    PUBLIC
    !==================================================================================================!
    ! Geometry
    !==================================================================================================!
    ! Note: we decided to separate all properties concerning the
    ! problem on the x axys from the ones on the y axys.
    ! The meaning of that is to give the user the possibility to
    ! calculate the results on non-sqared grids.
    !==================================================================================================!
    INTEGER             :: IMAX     ! Total number of cells on the x axys
    INTEGER             :: JMAX     ! Total number of cells on the y axys
    REAL                :: xL, xR   ! Domain boundaries for axis x
    REAL                :: yD, yU   ! Domain boundaries for axis y
    real                :: s        ! Initial condition
    REAL, ALLOCATABLE   :: x(:)     ! Vertex x coords (where u is defined)  
    REAL, ALLOCATABLE   :: y(:)     ! Vertex y coords (where v is defined)
    ! TODO Why not use a single matrix for barycenter coords?
    ! REAL, ALLOCATABLE :: xyb(:,:) ! Barycenter coords (where eta is defined)
    REAL, ALLOCATABLE   :: xb(:)    ! Barycenter x coords (where eta x is defined)
    REAL, ALLOCATABLE   :: yb(:)    ! Barycenter y coords (where eta y is defined)    
    REAL                :: dx, dx2  ! Mesh spacing on the x axis
    REAL                :: dy, dy2  ! Mesh spacing on the y axis
    !==================================================================================================!
    ! Discretization
    !==================================================================================================!
    REAL                :: CFL          ! CFL number for time step
    REAL, PARAMETER     :: G=9.81       ! Gravity acceleration    
    REAL, ALLOCATABLE   :: u  (:, :)    ! Velocity on the x axys
    REAL, ALLOCATABLE   :: v  (:, :)    ! Velocity on the y axis
    REAL, ALLOCATABLE   :: Fu (:, :)    ! Convective and viscous terms operator
    REAL, ALLOCATABLE   :: Fv (:, :)    ! Convective and viscous terms operator for y axis
    REAL, ALLOCATABLE   :: eta(:, :)    ! Pressure
    REAL, ALLOCATABLE   :: Hu (:, :)    ! Total wather depth for u mesh
    REAL, ALLOCATABLE   :: Hv (:, :)    ! Total wather depth for v mesh
    REAL, ALLOCATABLE   :: bu (:, :)    ! Bottom elevation for u mesh
    REAL, ALLOCATABLE   :: bv (:, :)    ! Bottom elevation for v mesh
    REAL                :: nu           ! kinematic viscosity coefficient
    REAL, ALLOCATABLE   :: rhs(:, :)    ! rhs of the pressure system 
    !==================================================================================================!
    !variables for parallelization
    !==================================================================================================!
#ifdef PARALLEL
    INTEGER             :: LCPU, RCPU, MsgLength, nMsg 
    REAL                :: send_messageL, send_messageR, recv_messageL, recv_messageR
    INTEGER             :: send_request(2), recv_request(2) 
    !INTEGER             :: send_status_list(MPI_STATUS_SIZE,2),recv_status_list(MPI_STATUS_SIZE,2)   
#endif     
    !------------------------------------------------------------!
  
    !==================================================================================================!
    ! Concerning time
    !==================================================================================================!
    REAL                :: time             ! Time iterator
    REAL                :: dt, dt2, dt_fix  ! Time steps
    REAL                :: tend             ! Ending time of the computation
    INTEGER, PARAMETER  :: NMAX=1e6         ! Maximum number of iterations
    !==================================================================================================!
    ! Output
    !==================================================================================================!
    CHARACTER(LEN=200)  :: TestName     ! Test name (also name of the output file)
    real                :: tio, dtio    ! Time value and step for output file
    !==================================================================================================!  
    CONTAINS
!======================================================================================================!
! Routine for allocating all variables
!======================================================================================================! 
SUBROUTINE Allocate_Var
    IMPLICIT NONE
    !==================================================================================================!
    ALLOCATE( x( IMAX+1 ), xb( IMAX ) )
    ALLOCATE( y( JMAX+1 ), yb( JMAX ) )
    x  = 0.
    xb = 0.
    y  = 0.
    yb = 0.
    !==================================================================================================!
    ALLOCATE( u(   IMAX+1, JMAX   ), Fu( IMAX+1, JMAX   ) )
    ALLOCATE( v(   IMAX,   JMAX+1 ), Fv( IMAX,   JMAX+1 ) )
    ALLOCATE( Hu(  IMAX+1, JMAX   ), bu( IMAX+1, JMAX   ) )
    ALLOCATE( Hv(  IMAX,   JMAX+1 ), bv( IMAX,   JMAX+1 ) )
#IFNDEF PARALLEL    
    ALLOCATE( eta( IMAX, JMAX )                           )    
#ENDIF    
    
  
    u   = 0.
    v   = 0.
    Fu  = 0.
    Fv  = 0.
    Hu  = 0.
    Hv  = 0.
    bu  = 0.
    bv  = 0.
    eta = 0.
    !==================================================================================================!
END SUBROUTINE Allocate_Var
!======================================================================================================!
! Routine for deallocating all variables
!======================================================================================================!
SUBROUTINE Deallocate_Var
    IMPLICIT NONE
    DEALLOCATE( x,  xb )
    DEALLOCATE( y,  yb )
    DEALLOCATE( u,  Fu )
    DEALLOCATE( v,  Fv )
    DEALLOCATE( Hu, bu )
    DEALLOCATE( Hv, bv )
    DEALLOCATE( eta    )
    DEALLOCATE( rhs)
END SUBROUTINE Deallocate_Var
!======================================================================================================!
END MODULE VarDef_mod   ! End of module VarDef