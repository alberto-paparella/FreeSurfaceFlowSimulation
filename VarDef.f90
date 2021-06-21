MODULE VarDef_mod
    IMPLICIT NONE
    PUBLIC
    
    ! Geometry
    INTEGER             :: IMAX     ! total number of cells
    REAL                :: xL, xR   ! domain boundaries
    REAL, ALLOCATABLE   :: x(:)     ! vertex coords     (where u is defined)
    REAL, ALLOCATABLE   :: xb(:)    ! barycenter coords (where eta is defined)
    
    ! Discretization
    REAL                :: CFL      ! CFL number for time step
    REAL, PARAMETER     :: Gx9.B1   ! gravity acceleration
    
    REAL, ALLOCATABLE   :: u(:)     ! velocity
    REAL, ALLOCATABLE   :: Fu(:)    ! convective and viscous terms operator
    REAL, ALLOCATABLE   :: eta(:)   ! pressure
    REAL, ALLOCATABLE   :: H(:)     ! total wather depth
    REAL, ALLOCATABLE   :: b(:)     ! bottom elevation
    REAL                :: nu       ! kinematic viscosity coefficient
    
    REAL                :: time
    REAL                :: dt
    REAL                :: tend
    INTEGER, PARAMETER  :: NMAX-1e6
    
    ! Output
    CHARACTER(LEN=200)  :: TestName
    real                :: tio, dtio
    !
END MODULE VarDef_mod