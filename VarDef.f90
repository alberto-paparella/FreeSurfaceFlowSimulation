MODULE VarDef_mod
    IMPLICIT NONE
    PUBLIC
    
    ! Geometry
    INTEGER             :: IMAX     ! total number of cells
    INTEGER             :: JMAX     ! 2D: Max index for y coordinates
    REAL                :: xL, xR   ! domain boundaries for axis x
    REAL                :: yD,yU    ! 2D: domain boundaries for axis y
    real                :: s        ! 2D: initial condition
    REAL, ALLOCATABLE   :: x(:)     ! vertex coords     (where u is defined)
    REAL, ALLOCATABLE   :: y(:)     ! 2D: vertex y coords (where v is defined)
    ! TODO Why not use a single matrix for barycenter coords?
    REAL, ALLOCATABLE   :: xb(:)    ! barycenter coords (where eta is defined)
    REAL, ALLOCATABLE   :: yb(:)    ! 2D: barycenter y coords
    !REAL, ALLOCATABLE   :: xyb(:,:) ! 2D: barycenter coords
    REAL                :: dx, dx2  ! mesh spacing
    REAL                :: dy, dy2  ! 2D: mesh spacing on the y axis
    
    ! Discretization
    REAL                :: CFL      ! CFL number for time step
    REAL, PARAMETER     :: G=9.81   ! gravity acceleration
    
    REAL, ALLOCATABLE   :: u(:, :)      ! velocity
    REAL, ALLOCATABLE   :: v(:, :)      ! 2D: velocity on the y axis
    REAL, ALLOCATABLE   :: Fu(:, :)     ! convective and viscous terms operator
    REAL, ALLOCATABLE   :: Fv(:, :)     ! 2D: convective and viscous terms operator for y axis
    REAL, ALLOCATABLE   :: eta(:, :)    ! pressure
    REAL, ALLOCATABLE   :: Hu(:, :)     ! total wather depth for u mesh
    REAL, ALLOCATABLE   :: Hv(:, :)     ! total wather depth for v mesh
    REAL, ALLOCATABLE   :: bu(:, :)     ! bottom elevation for u mesh
    REAL, ALLOCATABLE   :: bv(:, :)     ! bottom elevation for v mesh
    REAL                :: nu           ! kinematic viscosity coefficient
    REAL, ALLOCATABLE   :: rhs(:, :)    ! rhs of the pressure system 
    
    REAL                :: time
    REAL                :: dt, dt2, dt_fix
    REAL                :: tend
    INTEGER, PARAMETER  :: NMAX=1e6
    
    ! Output
    CHARACTER(LEN=200)  :: TestName
    real                :: tio, dtio
    !--------------------------------------------------------------------------------
    
    CONTAINS
    
SUBROUTINE Allocate_Var
    IMPLICIT NONE
    
    ALLOCATE( x(IMAX+1), xb(IMAX) )
    ALLOCATE( y(JMAX+1), yb(JMAX) )
    x  = 0.
    xb = 0.
    y  = 0.
    yb = 0.
    
    ALLOCATE( u(   IMAX+1, JMAX   ), Fu( IMAX+1, JMAX   ) )
    ALLOCATE( v(   IMAX,   JMAX+1 ), Fv( IMAX,   JMAX+1 ) )
    ALLOCATE( Hu(  IMAX+1, JMAX   ), bu( IMAX+1, JMAX   ) )
    ALLOCATE( Hv(  IMAX,   JMAX+1 ), bv( IMAX,   JMAX+1 ) )
    ALLOCATE( eta( IMAX, JMAX )                           )
    !ALLOCATE (rhs(IMAX,JMAX)                              )
    !ALLOCATE (chs(IMAX,JMAX)                              )
    !ALLOCATE( eta( IMAX, JMAX ) ,rhs(IMAX,JMAX)           )
   ! ALLOCATE( eta( IMAX, JMAX ) ,chs(IMAX,JMAX)           )
    u   = 0.
    v   = 0.
    Fu  = 0.
    Fv  = 0.
    Hu  = 0.
    Hv  = 0.
    bu  = 0.
    bv  = 0.
    eta = 0.
    rhs = 0.
    
END SUBROUTINE Allocate_Var

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
    
    
END MODULE VarDef_mod