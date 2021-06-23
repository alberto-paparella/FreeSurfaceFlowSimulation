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
    REAL, ALLOCATABLE   :: y(:)     ! 2D: vertex y coords
    REAL, ALLOCATABLE   :: xb(:)    ! barycenter coords (where eta is defined)
    REAL, ALLOCATABLE   :: yb(:)    ! 2D: barycenter y coords
    REAL                :: dx, dx2  ! mesh spacing
    REAL                :: dy, dy2  ! 2D: mesh spacing on the y axis
    
    ! Discretization
    REAL                :: CFL      ! CFL number for time step
    REAL, PARAMETER     :: G=9.81   ! gravity acceleration
    
    REAL, ALLOCATABLE   :: u(:, :)     ! velocity
    REAL, ALLOCATABLE   :: v(:, :)     ! 2D, TODO who exactly is this? (re-watch lesson)
    REAL, ALLOCATABLE   :: Fu(:, :)    ! convective and viscous terms operator
    REAL, ALLOCATABLE   :: Fv(:, :)    ! 2D, TODO who exactly is this? (re-watch lesson)
    REAL, ALLOCATABLE   :: eta(:, :)   ! pressure
    REAL, ALLOCATABLE   :: H(:, :)     ! total wather depth
    REAL, ALLOCATABLE   :: b(:, :)     ! bottom elevation
    !REAL, ALLOCATABLE   :: rhs(:)  ! rhs of the pressure system
    REAL                :: nu       ! kinematic viscosity coefficient
    
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
    
    ALLOCATE( x(IMAX+1), xb(IMAX+1) )
    ALLOCATE( y(JMAX+1), yb(JMAX+1) )   ! 2D
    x  = 0.
    xb = 0.
    y  = 0. ! 2D
    yb = 0. ! 2D
    
    ! ALLOCATE( u(IMAX+1), Fu(IMAX+1) )
    ALLOCATE( u(IMAX+1, JMAX+1), Fu(IMAX+1, JMAX+1) )   ! 2D
    ALLOCATE( v(IMAX+1, JMAX+1), Fv(IMAX+1, JMAX+1) )   ! 2D
    ! ALLOCATE( H(IMAX+1), b(IMAX+1)  )
    ALLOCATE( H(IMAX+1, JMAX+1), b(IMAX+1, JMAX+1)  )   ! 2D
    ! ALLOCATE( eta(IMAX)             )
    ALLOCATE( eta(IMAX, JMAX)                       )   ! 2D
    !ALLOCATE( eta(IMAX), rhs(IMAX)   )
    u   = 0.
    Fu  = 0.
    H   = 0.
    b   = 0.
    eta = 0.
    !rhs = 0.
    
END SUBROUTINE Allocate_Var

SUBROUTINE Deallocate_Var
    IMPLICIT NONE
    DEALLOCATE( x, xb    )
    DEALLOCATE( y, yb    )
    DEALLOCATE( u, Fu    )
    DEALLOCATE( v, Fv    ) ! 2D
    DEALLOCATE( H, b     )
    DEALLOCATE( eta      )
    !DEALLOCATE( eta, rhs )
END SUBROUTINE Deallocate_Var
    
    
    
END MODULE VarDef_mod