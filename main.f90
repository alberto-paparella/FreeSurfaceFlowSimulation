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
! main.f90
! This is the main file
!======================================================================================================!
PROGRAM FS2D
    !==================================================================================================!
    USE VarDef_mod  ! In this module we defined the variables we are going to use
    IMPLICIT NONE
    !==================================================================================================!
    INTEGER            :: i, j, n
    REAL               :: umax, au, av, s0, ct, cs , a , b
    REAL, ALLOCATABLE  :: amat(:,:), amatj(:,:), bmat(:,:), cmat(:,:), cmatj(:,:)
    !==================================================================================================!
    WRITE(*,'(a)') ' | ================================================================== | '
    WRITE(*,'(a)') ' |  Equazione dei flussi a superficie libera su griglie strutturate.  | '
    WRITE(*,'(a)') ' | ================================================================== | '
    WRITE(*,'(a)') ' | '
    WRITE(*,'(a)') ' | Soluzione con il metodo del gradiente coniugato. '
    !==================================================================================================!
    ! SETTINGS OF THE COMPUTATION
    !==================================================================================================!
    TestName = 'Gaussian_test'
    ! Concerning x axys
    IMAX   = 500    ! Max index for x coordinates
    xL     = -0.5   ! Left boundary for x coordinates
    xR     = 0.5    ! Right boundary for y coordinates
    ! Concerning y axys
    JMAX   = 500    ! Max index for y coordinates
    yD     = -0.5   ! Left boundary for y coordinates
    yU     = 0.5    ! Right boundary for y coordinates
    ! Initial Conditions
    s      = 0.1    ! Initial condition
    ! Concerning time
    time   = 0.0    ! Starting of the computation time
    tend   = 0.1    ! Ending of the computation time
    tio    = 0.0    ! Time iterator for output file
    dtio   = 5e-2   ! Time step for output file
    CFL    = 0.9    ! CFL number for time step
    dt_fix = 2e-2   ! Fixed time step for computation
    nu     = 0.e-2  ! Kinematic viscosity coefficient
    !==================================================================================================!
    CALL Allocate_Var  ! Allocate all variables
    !==================================================================================================!
    ! Matrixes for assembling the tridiagonal linear system
    !==================================================================================================!
    ALLOCATE( amat  ( IMAX, JMAX ) )    ! Matrix a for i system
    ALLOCATE( amatj ( IMAX, JMAX ) )    ! Matrix a for j system
    ALLOCATE( bmat  ( IMAX, JMAX ) )    ! Matrix b for both systems
    ALLOCATE( cmat  ( IMAX, JMAX ) )    ! Matrix c for i system
    ALLOCATE( cmatj ( IMAX, JMAX ) )    ! Matrix c for j system
    ALLOCATE( rhs   ( IMAX, JMAX ) )    ! Matrix rhs 
    !==================================================================================================!
    ! 1) Computational domain 
    !==================================================================================================!
    WRITE(*,'(a)') ' | Building computational domain... '
    ! Domain on the x axys
    dx    = (xR - xL) / REAL(IMAX)
    dx2   = dx**2
    x(1)  = xL
    ! Domain on the y axys
    dy    = (yU - yD) / REAL(IMAX)
    dy2   = dy**2
    y(1)  = yD
    !==================================================================================================!
    ! We decided to use two separates dimensions, IMAX and JMAX,
    ! for x and y coords, to offer more flexibility
    !==================================================================================================!
    DO i = 1, IMAX
        x(i+1) = x(i) + dx
        xb(i)  = x(i) + dx/2.
    ENDDO
    DO j = 1, JMAX
        y(j+1) = y(j) + dy
        yb(j)  = y(j) + dy/2.
    ENDDO
    !==================================================================================================!
    ! 2) Initial condition
    !==================================================================================================!
    WRITE(*,'(a)') ' | Assigning initial condition... '
    !==================================================================================================!
    ! 2.1) Free surface elevation (barycenters)
    !==================================================================================================!
    DO i = 1, IMAX
        DO j = 1, JMAX
            !===========================================================================================!
            ! Gaussian profile
            !===========================================================================================!
            ! Note: we use xb and yb instead of x and y because x and y represents the vertex
            ! coords for u and v, while xb and yb represents the barycenter coords for eta
            !===========================================================================================!
            eta(i,j) = 1.0 + EXP( (-1.0 / (2 * (s**2) ) ) * ( xb(i)**2 + yb(j)**2 ) )
        ENDDO
    ENDDO
    !==================================================================================================!
    ! 2.1) Velocity, bottom and total water depth (interfaces)
    !==================================================================================================!
    ! Note: u dimensions (IMAX + 1) * JMAX, while v dimensions are IMAX * (JMAX + 1)
    !==================================================================================================!
    DO i = 1, IMAX + 1
        DO j = 1, JMAX
            u(i, j)  = 0.0  ! Fluid at rest
            bu(i, j) = 0.0  ! Zero bottom elevation
        ENDDO
    ENDDO
    DO i = 1, IMAX
        DO j = 1, JMAX + 1
            v(i, j)  = 0.0  ! Fluid at rest
            bv(i, j) = 0.0  ! Zero bottom elevation
        ENDDO
    ENDDO
    !==================================================================================================!
    ! Total water depth
    !==================================================================================================!
    DO j = 1, JMAX
        ! Initialize first and last columns for Hu
        Hu( 1,      j ) = MAX( 0.0, bu( 1,      j ) + eta( 1,    j ) )
        Hu( IMAX+1, j ) = MAX( 0.0, bu( IMAX+1, j ) + eta( IMAX, j ) )
    ENDDO
    DO i = 1, IMAX
        ! Initialize first and last rows for Hv
        Hv( i, 1      ) = MAX( 0.0, bv( i, 1      ) + eta( i, 1    ) )
        Hv( i, JMAX+1 ) = MAX( 0.0, bv( i, JMAX+1 ) + eta( i, JMAX ) )
    ENDDO
    DO i = 2, IMAX
        DO j = 1, JMAX
            Hu( i , j ) = MAXVAL( (/ 0.0, bu( i, j ) + eta( i - 1, j ), bu( i, j )+ eta( i, j ) /) )
        ENDDO
    ENDDO
    DO i = 1, IMAX
        DO j = 2, JMAX
            Hv( i , j ) = MAXVAL( (/ 0.0, bv( i, j ) + eta( i, j - 1 ), bv( i, j )+ eta( i, j ) /) )
        ENDDO
    ENDDO
    !==================================================================================================!
    ! Plot initial condition
    !==================================================================================================!
    CALL DataOutput(0)
    tio = tio + dtio
    !==================================================================================================!
    ! 3) Computation: main loop in time
    !==================================================================================================!
    WRITE(*,'(a)') ' | '
    WRITE(*,'(a)') ' | START of the COMPUTATION '
    !==================================================================================================!
    DO n = 1, NMAX
        !==============================================================================================!
        ! 3.1) Compute the time step
        !==============================================================================================!
        IF( time.GE.tend ) THEN
            EXIT    
        ENDIF
        !==============================================================================================!
        umax = MAXVAL( ABS(u) )             ! Semi-implicit time step
        dt   = MIN( dt_fix, CFL / ( ( umax + 1e-14 ) / dx + 2. * nu / dx2 ) )
        dt2   = dt**2
        IF( ( time + dt ).GT.tend ) THEN
            dt  = tend - time   
            tio = tend
        ENDIF    
        IF( ( time + dt ).GT.tio ) THEN
            dt = tio - time
        ENDIF
        !==============================================================================================!
        ! 3.2) Compute the operator Fu
        !==============================================================================================!        
        ! Fu = u    ! Case without considering Fu
        ! BC: no-slip wall
        ! First row and last row initialized to zeros
        Fu( 1      , : ) = Fu(1,:)
        Fu( IMAX+1 , : ) = Fu(IMAX+1,:)
        DO i = 2, IMAX
            DO j = 1, JMAX
                au = ABS( u(i,j) )
                ! Explicit upwind
                ! In our case, nu = 0
                !Fu(i,j) = ( 1. - dt * ( au/dx + 2.*nu/dx2 ) ) * u(i,j)   &     ! q^N
                !        + dt * ( nu/dx2 + (au-u(i,j))/(2.*dx) ) * u(i+1,j) &   ! a^+(...)
                !        + dt * ( nu/dx2 + (au+u(i,j))/(2.*dx) ) * u(i-1,j)     ! a^-(...)
                Fu(i,j) = ( 1. - dt * ( au/dx ) ) * u(i,j)   &
                       + dt * ( (au-u(i,j))/(2.*dx) ) * u(i+1,j) &
                       + dt * ( (au+u(i,j))/(2.*dx) ) * u(i-1,j)
            ENDDO
        ENDDO
        !==============================================================================================!
        ! 3.3) Compute the operator Fv
        !==============================================================================================!
        ! Fv = v    ! Case without considering Fv
        ! BC: no-slip wall
        ! First column and last column initialized to zeros
        Fv( : , 1      ) = Fv(:,1)
        Fv( : , JMAX+1 ) = Fv(:,JMAX+1)
        DO i = 1, IMAX
            DO j = 2, JMAX
                av = ABS( v(i,j) )
                ! Explicit upwind
                Fv(i,j) = ( 1. - dt * ( av/dx + 2.*nu/dy2 ) ) * v(i,j)     &
                        + dt * ( nu/dy2 + (av-v(i,j))/(2.*dy) ) * v(i,j+1) &
                        + dt * ( nu/dy2 + (av+v(i,j))/(2.*dy) ) * v(i,j-1)
            ENDDO
        ENDDO
        !==============================================================================================!
        ! 3.4) Solve the free surface equation
        !==============================================================================================!
        ! CONJUGATE GRADIENT METHOD
        !==============================================================================================!
        DO i = 1, IMAX
            DO j = 1, JMAX
                rhs(i,j) = eta(i,j) - dt/dx * ( Hu(i+1,j)*Fu(i+1,j) - Hu(i,j)*Fu(i,j)) &
                                    - dt/dy * ( Hv(i,j+1)*Fv(i,j+1) - Hv(i,j)*Fv(i,j))
            ENDDO
        ENDDO
        ! Here, we suppose that IMAX and JMAX are the same, at least for the moment
        CALL CG(IMAX,eta,rhs)
        !==============================================================================================!
        ! 3.5) Update the velocity (momentum equation)
        !==============================================================================================!
        ct = g*dt/dx    ! Temporary coefficient for x
        cs = g*dt/dy    ! Temporary coefficient for y
        !==============================================================================================!
        u(1,:)      = Fu(1,:) 
        u(IMAX+1,:) = Fu(IMAX+1,:)
        DO i = 2, IMAX
            DO j = 1, JMAX
                u(i,j) = Fu(i,j) - ct * ( eta(i,j) - eta(i-1,j) )
            ENDDO
        ENDDO
        !==============================================================================================!
        v(:,1)      = Fv(:,1)
        v(:,JMAX+1) = Fv(:,JMAX+1)
        DO i = 1, IMAX
            DO j = 2, JMAX
                v(i,j) = Fv(i,j) - cs * ( eta(i,j) - eta(i,j-1) )
            ENDDO
        ENDDO
        !==============================================================================================!
        ! 3.6) Update total water depth
        !==============================================================================================!
        DO j = 1, JMAX
            ! Initialize first and last columns for Hu
            Hu( 1,      j ) = MAX( 0.0, bu( 1,      j ) + eta( 1,    j ) )
            Hu( IMAX+1, j ) = MAX( 0.0, bu( IMAX+1, j ) + eta( IMAX, j ) )
        ENDDO
        DO i = 1, IMAX
            ! Initialize first and last rows for Hv
            Hv( i, 1      ) = MAX( 0.0, bv( i, 1      ) + eta( i, 1    ) )
            Hv( i, JMAX+1 ) = MAX( 0.0, bv( i, JMAX+1 ) + eta( i, JMAX ) )
        ENDDO
        DO i = 2, IMAX
            DO j = 1, JMAX
                Hu( i , j ) = MAXVAL( (/ 0.0, bu( i, j ) + eta( i - 1, j ), bu( i, j )+ eta( i, j ) /) )
            ENDDO
        ENDDO
        DO i = 1, IMAX
            DO j = 2, JMAX
                Hv( i , j ) = MAXVAL( (/ 0.0, bv( i, j ) + eta( i, j - 1 ), bv( i, j )+ eta( i, j ) /) )
            ENDDO
        ENDDO
        !==========================================================================================!
        time = time + dt  ! Update time
        !==========================================================================================!
        ! 3.7) Eventually plot the results
        !CALL DataOutput(n)  ! DEBUG
        !STOP                ! DEBUG
        IF(ABS(time-tio).LT.1e-12) THEN
            WRITE(*,'(a,f15.7)') ' |   plotting data output at time ', time
            CALL DataOutput(n)
            tio = tio + dtio
        ENDIF  
    !==============================================================================================!
    ENDDO ! n cycle
    !==============================================================================================!
    WRITE(*,'(a)') ' | END of the COMPUTATION '
    !==============================================================================================!
    ! Empty memory
    !==============================================================================================!
    CALL Deallocate_Var
    !==============================================================================================!
    DEALLOCATE( amat,bmat,cmat, amatj,cmatj )
    !==============================================================================================!
    WRITE(*,'(a)') ' | '
    WRITE(*,'(a)') ' |         Finalization was successful. Bye :-)           | '
    WRITE(*,'(a)') ' | ====================================================== | '
    !==============================================================================================!    
    PAUSE
!======================================================================================================!
END PROGRAM FS2D
!======================================================================================================!
SUBROUTINE matop2D(Ap,p,N)
    !==================================================================================================!
    USE VarDef_mod
    IMPLICIT NONE
    !==================================================================================================!
    INTEGER  :: N 
    REAL     :: Ap(N,N), p(N,N)
    !==================================================================================================!
    INTEGER  :: i, j
    REAL     :: cx, cy
    !==================================================================================================!
    cx = g*dt**2/dx2    ! Temporary coefficient
    cy = g*dt**2/dy2
    !==================================================================================================!	
	! 1) Diagonal term
    !==================================================================================================!	
	Ap = p
    !==================================================================================================!	
	! 2) Fluxes in x-direction
    !==================================================================================================!	
	DO j = 1, N
		DO i = 1, N
			IF(i.EQ.1) THEN
				Ap(i,j) = Ap(i,j) - cx * ( Hu(i+1,j)*(p(i+1,j)-p(i,j)) - 0.0 )
			ELSEIF(i.EQ.N) THEN  
				Ap(i,j) = Ap(i,j) - cx * ( 0.0 - Hu(i,j)*(p(i,j)-p(i-1,j)) )
			ELSE  
				Ap(i,j) = Ap(i,j) - cx * ( Hu(i+1,j)*(p(i+1,j)-p(i,j)) - Hu(i,j)*(p(i,j)-p(i-1,j)) )
			ENDIF 
		ENDDO
    ENDDO
    !==================================================================================================!	
    ! 3) Fluxes in y-direction
    !==================================================================================================!	
	DO j = 1, N
		DO i = 1, N
			IF(j.EQ.1) THEN
				Ap(i,j) = Ap(i,j) - cy * ( Hv(i,j+1)*(p(i,j+1)-p(i,j)) - 0.0 )
			ELSEIF(j.EQ.N) THEN  
				Ap(i,j) = Ap(i,j) - cy * ( 0.0 - Hv(i,j)*(p(i,j)-p(i,j-1)) )
			ELSE  
				Ap(i,j) = Ap(i,j) - cy * ( Hv(i,j+1)*(p(i,j+1)-p(i,j)) - Hv(i,j)*(p(i,j)-p(i,j-1)) )
			ENDIF 
		ENDDO
    ENDDO    
!======================================================================================================!
    END SUBROUTINE matop2D
    
    

!======================================================================================================!
!PARALLELIZATION 
!======================================================================================================!
    
#define PARALLEL
  
PROGRAM FS2D_MPI
    USE VarDef_MPI  ! In this module we defined the variables we are going to use
    !------------------------------------------------------------!
    IMPLICIT NONE
#ifdef PARALLEL     
    INCLUDE 'mpif.h'
#endif     
    !------------------------------------------------------------!
    INTEGER            :: i, j, n
    REAL               :: umax, au, av, s0, ct, cs , a , b
    REAL, ALLOCATABLE  :: amat(:,:), amatj(:,:), bmat(:,:), cmat(:,:), cmatj(:,:)
    !==================================================================================================!
    !
    !INTEGER, PARAMETER  :: IMAX=500  ! total number of cells 
    !REAL                :: xL, xR      ! left and right domain definition
    !REAL                :: xD          ! location of discontinuity
    !REAL                :: dx, dx2     ! mesh size and square of mesh size
    !REAL, ALLOCATABLE   :: x(:)        ! vertex coords
    !!
    !REAL                :: CFL         ! Courant-Friedrichs-Lewy number for stability condition (CFL<1) 
    !REAL                :: time        ! current time
    !REAL                :: dt          ! time step
    !REAL                :: tend        ! final time
    !INTEGER, PARAMETER  :: NMAX = 1e6  ! maximum number of time steps 
    !REAL, ALLOCATABLE   :: T(:), T1(:) ! temperature 
    !!
    !CHARACTER(LEN=200)  :: TestName    ! name of test problem
    !REAL                :: tio         ! output time and time step 
    !REAL                :: dtio        ! output time step
    !!
    !REAL                :: d           ! stability parameter  
    !REAL                :: kappa       ! free-surface conduction coefficient
    !REAL                :: TL, TR      ! boundary conditions
!#ifdef PARALLEL
!    INTEGER             :: LCPU, RCPU, MsgLength, nMsg 
!    REAL                :: send_messageL, send_messageR, recv_messageL, recv_messageR
!    INTEGER             :: send_request(2), recv_request(2) 
!    INTEGER             :: send_status_list(MPI_STATUS_SIZE,2),recv_status_list(MPI_STATUS_SIZE,2)   
!#endif     
!    !     
!    TYPE tMPI
!    INTEGER :: myrank, nCPU, iErr 
!    INTEGER :: AUTO_REAL
!    INTEGER :: nElem,istart,iend 
!    END TYPE tMPI 
!    TYPE(tMPI) :: MPI
!    REAL       :: realtest
!    REAL       :: WCT1, WCT2           ! Wall clock times 
!    !
!    !------------------------------------------------------------!
!    !
!    WCT1 = 0. 
!    WCT2 = 0. 
!    !
#ifdef PARALLEL
    CALL MPI_INIT(MPI%iErr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, MPI%myrank, MPI%iErr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, MPI%nCPU,   MPI%iErr)
    SELECT CASE(KIND(realtest))
    CASE(4)
      IF(MPI%myrank.EQ.0) THEN
          PRINT *, ' Single precision used for real. '
          PRINT *, ' Setting MPI%AUTO_REAL to MPI_REAL. '
      ENDIF
      MPI%AUTO_REAL = MPI_REAL
    CASE(8)
      IF(MPI%myrank.EQ.0) THEN
          PRINT *,  '  Double precision used for real.'
          PRINT *,  '  Setting MPI%AUTO_REAL to MPI_DOUBLE_PRECISION'
      ENDIF 
      MPI%AUTO_REAL = MPI_DOUBLE_PRECISION
    END SELECT 
    IF(MOD(IMAX,MPI%nCPU).NE.0) THEN
         IF(MPI%myrank.EQ.0) THEN
             PRINT *, ' Error. The number of mesh points must be a multiple of the number of CPUs. '
             PRINT *, ' IMAX = ', IMAX, ' nCPU = ', MPI%nCPU
         ENDIF  
         CALL MPI_FINALIZE(MPI%iErr) 
         STOP 
    ENDIF 
    MPI%nElem  = IMAX/MPI%nCPU 
    MPI%istart = 1 + MPI%myrank*MPI%nElem 
    MPI%iend  = MPI%istart + MPI%nElem - 1 
    ! 
#else
    MPI%myrank = 0              ! If the code is not compiled in parallel, then there is only one CPU ...
    MPI%nCPU   = 1              ! ...
    MPI%AUTO_REAL = -1          ! ... and the data type MPI_AUTO_REAL is not needed 
    MPI%istart = 1
    MPI%iend   = IMAX     
#endif    
    !
    IF(MPI%myrank.EQ.0) THEN
      WRITE(*,'(a)') ' '
      WRITE(*,'(a)') ' ====================================================== '
      WRITE(*,'(a)') '              Finite difference code for 2D free-surface flow equation       '
      WRITE(*,'(a,i10)') '          Total number of CPUs used : ', MPI%nCPU       
      WRITE(*,'(a)') ' ====================================================== '
      WRITE(*,'(a)') ' '
    ENDIF  
    !
    
    !==================================================================================================!
    ! SETTINGS OF THE COMPUTATION
    !==================================================================================================!
    TestName = 'Gaussian_test'
    ! Concerning x axys
    IMAX   = 500    ! Max index for x coordinates
    xL     = -0.5   ! Left boundary for x coordinates
    xR     = 0.5    ! Right boundary for y coordinates
    ! Concerning y axys
    JMAX   = 500    ! Max index for y coordinates
    yD     = -0.5   ! Left boundary for y coordinates
    yU     = 0.5    ! Right boundary for y coordinates
    ! Initial Conditions
    s      = 0.1    ! Initial condition
    ! Concerning time
    time   = 0.0    ! Starting of the computation time
    tend   = 0.1    ! Ending of the computation time
    tio    = 0.0    ! Time iterator for output file
    dtio   = 5e-2   ! Time step for output file
    CFL    = 0.9    ! CFL number for time step
    dt_fix = 2e-2   ! Fixed time step for computation
    nu     = 0.e-2  ! Kinematic viscosity coefficient
    !==================================================================================================!
    CALL Allocate_Var  ! Allocate all variables
    !==================================================================================================!
    ! Matrixes for assembling the tridiagonal linear system
    !==================================================================================================!
    ALLOCATE( amat  ( IMAX, JMAX ) )    ! Matrix a for i system
    ALLOCATE( amatj ( IMAX, JMAX ) )    ! Matrix a for j system
    ALLOCATE( bmat  ( IMAX, JMAX ) )    ! Matrix b for both systems
    ALLOCATE( cmat  ( IMAX, JMAX ) )    ! Matrix c for i system
    ALLOCATE( cmatj ( IMAX, JMAX ) )    ! Matrix c for j system
    ALLOCATE( rhs   ( IMAX, JMAX ) )    ! Matrix rhs 
    !==================================================================================================!
    ! 1) Computational domain 
    !==================================================================================================!
    WRITE(*,'(a)') ' | Building computational domain... '
    ! Domain on the x axys
    dx    = (xR - xL) / REAL(IMAX)
    dx2   = dx**2
    x(1)  = xL
    ! Domain on the y axys
    dy    = (yU - yD) / REAL(IMAX)
    dy2   = dy**2
    y(1)  = yD
    !==================================================================================================!
    ! We decided to use two separates dimensions, IMAX and JMAX,
    ! for x and y coords, to offer more flexibility
    !==================================================================================================!
    DO i = 1, IMAX
        x(i+1) = x(i) + dx
        xb(i)  = x(i) + dx/2.
    ENDDO
    DO j = 1, JMAX
        y(j+1) = y(j) + dy
        yb(j)  = y(j) + dy/2.
    ENDDO
    !!======================================== 
    !! SETTINGS OF THE COMPUTATION
    !!========================================
    !!
    !TestName = 'Free-surface 2D_explicit'
    !!
    !xL     = -1.0
    !xR     = 1.0
    !xD     = 0.0  ! location of the discontinuity
    !!
    !time   = 0.0
    !tend   = 0.02    
    !dtio   = 2e-2
    !tio    = dtio
    !CFL    = 0.9   ! for stability condition CFL<1
    !d      = 0.45  ! for stability condition d<0.5   
    !!
    !kappa  = 1     ! heat conduction coefficient 
    !!
    !! Boundary conditions
    !TL = 100
    !TR = 50
    !
    !========================================
    !
    ! Allocate variables
    !ALLOCATE( x(IMAX)  )
    !ALLOCATE( T1(IMAX) )
    !ALLOCATE( T(MPI%istart-1:MPI%iend+1) )
    !
    !---------------------------------------- 
    !---------------------------------------- 
    !!
    !! 1) Computational domain 
    !!
    !IF(MPI%myrank.EQ.0) THEN
    !  WRITE(*,'(a)') '  Building computational domain... '
    !ENDIF  
    !!    
    !dx    = (xR - xL) / REAL(IMAX)
    !dx2   = dx**2
    !x(1)  = xL
    !!
    !DO i = 1, IMAX
    !  x(i) = xL + 0.5*dx + (i-1)*dx
    !ENDDO  
    !
    !---------------------------------------- 
    !---------------------------------------- 
    !
     !==================================================================================================!
    ! 2) Initial condition
    !==================================================================================================!
    WRITE(*,'(a)') ' | Assigning initial condition... '
    !==================================================================================================!
    ! 2.1) Free surface elevation (barycenters)
    !==================================================================================================!
    DO i = 1, IMAX
        DO j = 1, JMAX
            !===========================================================================================!
            ! Gaussian profile
            !===========================================================================================!
            ! Note: we use xb and yb instead of x and y because x and y represents the vertex
            ! coords for u and v, while xb and yb represents the barycenter coords for eta
            !===========================================================================================!
            eta(i,j) = 1.0 + EXP( (-1.0 / (2 * (s**2) ) ) * ( xb(i)**2 + yb(j)**2 ) )
        ENDDO
    ENDDO
    !==================================================================================================!
    ! 2.1) Velocity, bottom and total water depth (interfaces)
    !==================================================================================================!
    ! Note: u dimensions (IMAX + 1) * JMAX, while v dimensions are IMAX * (JMAX + 1)
    !==================================================================================================!
    DO i = 1, IMAX + 1
        DO j = 1, JMAX
            u(i, j)  = 0.0  ! Fluid at rest
            bu(i, j) = 0.0  ! Zero bottom elevation
        ENDDO
    ENDDO
    DO i = 1, IMAX
        DO j = 1, JMAX + 1
            v(i, j)  = 0.0  ! Fluid at rest
            bv(i, j) = 0.0  ! Zero bottom elevation
        ENDDO
    ENDDO
    !==================================================================================================!
    ! Total water depth
    !==================================================================================================!
    DO j = 1, JMAX
        ! Initialize first and last columns for Hu
        Hu( 1,      j ) = MAX( 0.0, bu( 1,      j ) + eta( 1,    j ) )
        Hu( IMAX+1, j ) = MAX( 0.0, bu( IMAX+1, j ) + eta( IMAX, j ) )
    ENDDO
    DO i = 1, IMAX
        ! Initialize first and last rows for Hv
        Hv( i, 1      ) = MAX( 0.0, bv( i, 1      ) + eta( i, 1    ) )
        Hv( i, JMAX+1 ) = MAX( 0.0, bv( i, JMAX+1 ) + eta( i, JMAX ) )
    ENDDO
    DO i = 2, IMAX
        DO j = 1, JMAX
            Hu( i , j ) = MAXVAL( (/ 0.0, bu( i, j ) + eta( i - 1, j ), bu( i, j )+ eta( i, j ) /) )
        ENDDO
    ENDDO
    DO i = 1, IMAX
        DO j = 2, JMAX
            Hv( i , j ) = MAXVAL( (/ 0.0, bv( i, j ) + eta( i, j - 1 ), bv( i, j )+ eta( i, j ) /) )
        ENDDO
    ENDDO
    
    !! 2) Initial condition
    !!
    !IF(MPI%myrank.EQ.0) THEN
    !  WRITE(*,'(a)') '  Assigning initial condition... '
    !ENDIF  
    !!
    !DO i = MPI%istart, MPI%iend
    !  IF(x(i).LT.xD) THEN
    !    T(i) = TL
    !  ELSE
    !    T(i) = TR
    !  ENDIF  
    !ENDDO  
    !!    
    !---------------------------------------- 
    !---------------------------------------- 
    !
    ! 3) Computation: main loop in time
    !
    IF(MPI%myrank.EQ.0) THEN
      WRITE(*,'(a)') '  '
      WRITE(*,'(a)') '  Explicit finite difference solver. '
      WRITE(*,'(a)') '  '
    ENDIF  
    !
#ifdef PARALLEL
    WCT1 = MPI_WTime() 
#else
    CALL CPU_TIME(WCT1) 
#endif    
    !
    DO n = 1, NMAX
      !
      IF(MPI%myrank.EQ.0) THEN
        IF(MOD(n,100).EQ.0) THEN
          WRITE(*,'(a,i10)') '  Current timestep = ', n 
        ENDIF
      ENDIF
      !
      ! 3.1) Compute the time step
      !
      IF(time.GE.tend) THEN
        EXIT    
      ENDIF
      !
      dt = d * dx2 / kappa
      IF((time+dt).GT.tend) THEN
        dt  = tend - time   
        tio = tend
      ENDIF    
      IF((time+dt).GT.tio) THEN
        dt = tio - time
      ENDIF
      !
      ! 3.2) Numerical scheme: FTCS
      !
#ifdef PARALLEL
      !
      ! The only real MPI part is here: exchange of the boundary values between CPUs 
      !
      MsgLength = 1 
      RCPU = MPI%myrank + 1 
      IF(RCPU.GT.MPI%nCPU-1) THEN
          RCPU = 0
      ENDIF
      LCPU = MPI%myrank - 1 
      IF(LCPU.LT.0) THEN 
          LCPU = MPI%nCPU - 1 
      ENDIF 
      ! send your leftmost state to your left neighbor CPU 
      send_messageL = T(MPI%istart) 
      CALL MPI_ISEND(send_messageL, MsgLength, MPI%AUTO_REAL, LCPU, 1, MPI_COMM_WORLD, send_request(1), MPI%iErr)
      ! send your rightmost state to your right neighbor CPU => post a send request  
      send_messageR = T(MPI%iend)
      CALL MPI_ISEND(send_messageR, MsgLength, MPI%AUTO_REAL, RCPU, 2, MPI_COMM_WORLD, send_request(2), MPI%iErr)
      ! obviously, the right neighbor must expect this message, so tell him to go regularly to the post office...  
      CALL MPI_IRECV(recv_messageL, MsgLength, MPI%AUTO_REAL, LCPU, 2, MPI_COMM_WORLD, recv_request(1), MPI%iErr)            
      ! obviously, the left neighbor must expect this message, so tell him to go regularly to the post office...  
      CALL MPI_IRECV(recv_messageR, MsgLength, MPI%AUTO_REAL, RCPU, 1, MPI_COMM_WORLD, recv_request(2), MPI%iErr)               
      !
      ! We do some useful work here, e.g. update all internal elements of the domain, since we do non-blocking communication ! 
      !
      ! The main finite difference scheme is IDENTICAL to the serial code ! 
      ! This is the simplicify and the beauty of MPI :-) 
      !   
      DO i = MPI%istart+1, MPI%iend-1  
          T1(i)= T(i) + kappa*dt/(dx2) * ( T(i-1) - 2*T(i) + T(i+1) ) 
      ENDDO   
      !
      ! Wait until all communication has finished
      ! 
      nMsg = 2 
      CALL MPI_WAITALL(nMsg,send_request,send_status_list,MPI%ierr)
      CALL MPI_WAITALL(nMsg,recv_request,recv_status_list,MPI%ierr)
      !
      ! Now we are sure all the messages have been sent around, so let's put the data of the messages where it actually belongs... 
      !
      T(MPI%istart-1) = recv_messageL 
      T(MPI%iend+1)   = recv_messageR
      !
      ! Update the MPI boundary cell, now that we have all data
      !
      DO i = MPI%istart, MPI%iend, (MPI%iend-MPI%iStart)
        IF(i.EQ.1) THEN
          T1(i) = TL  ! left BC 
        ELSEIF(i.EQ.IMAX) THEN
          T1(i) = TR  ! right BC
        ELSE  
          T1(i)= T(i) + kappa*dt/(dx2) * ( T(i-1) - 2*T(i) + T(i+1) ) 
        ENDIF  
      ENDDO
      !
#else      
      DO i = 2, (IMAX-1)
        T1(i)= T(i) + kappa*dt/(dx2) * ( T(i-1) - 2*T(i) + T(i+1) ) 
      ENDDO
      ! BCs
      T1(1)    = TL   
      T1(IMAX) = TR   
#endif      
      !
      ! update time
      time = time + dt  
      ! overwrite current solution
      T(MPI%istart:MPI%iend) = T1(MPI%istart:MPI%iend)         
      !
      ! 3.3) Eventually plot the results
      IF(ABS(time-tio).LT.1e-12) THEN
        IF(MPI%myrank.EQ.0) THEN
          WRITE(*,'(a,f15.7)') '   => plotting data output at time ', time
        ENDIF  
#ifdef JOINT_OUTPUT
#ifdef PARALLEL  
        ! Collect all the data from each MPI subdomain 
        CALL MPI_ALLGATHER(T(MPI%istart:MPI%iend),MPI%nElem,MPI%AUTO_REAL,T1,MPI%nElem,MPI%AUTO_REAL,MPI_COMM_WORLD,MPI%iErr)
        IF(MPI%myrank.EQ.0) THEN
            CALL DataOutput(n,TestName,MPI%myrank,1,IMAX,IMAX,x,T1)
        ENDIF
#endif       
#else   
        CALL DataOutput(n,TestName,MPI%myrank,MPI%istart,MPI%iend,IMAX,x,T(MPI%istart:MPI%iend))
#endif 
        tio = tio + dtio
      ENDIF  
      !
    ENDDO !n
    !
#ifdef PARALLEL
    WCT2 = MPI_WTime() 
    CALL MPI_FINALIZE(MPI%iErr) 
#else
    CALL CPU_TIME(WCT2) 
#endif    
    !---------------------------------------- 
    !---------------------------------------- 
    !
    ! Empty memory
    !
    DEALLOCATE( x, T, T1 )
    !
    IF(MPI%myrank.EQ.0) THEN
      WRITE(*,'(a)') ' '
      WRITE(*,'(a,f15.7)'), '         Total wallclock time needed : ', WCT2-WCT1
      WRITE(*,'(a)') '         Finalization was successful. Bye :-)           '
      WRITE(*,'(a)') ' ====================================================== '
    ENDIF  
    !
    END PROGRAM FS2D_MPI

  
    
    #define noIMP_SOLVER    ! preprocessor for implicit solver
  
!PROGRAM Heat1D
!    !------------------------------------------------------------!
!    IMPLICIT NONE
!    !------------------------------------------------------------!
!    INTEGER             :: i, n
!    !
!    INTEGER             :: IMAX        ! total number of cells 
!    REAL                :: xL, xR      ! left and right domain definition
!    REAL                :: xD          ! location of discontinuity
!    REAL                :: dx, dx2     ! mesh size and square of mesh size
!    REAL, ALLOCATABLE   :: x(:)        ! vertex coords
!    !
!    REAL                :: CFL         ! Courant-Friedrichs-Lewy number for stability condition (CFL<1) 
!    REAL                :: time        ! current time
!    REAL                :: dt          ! time step
!    REAL                :: dt_fix      ! user defined time step
!    REAL                :: tend        ! final time
!    INTEGER, PARAMETER  :: NMAX = 1e6  ! maximum number of time steps 
!    REAL, ALLOCATABLE   :: T(:), T1(:) ! temperature 
!    !
!    CHARACTER(LEN=200)  :: TestName    ! name of test problem
!    REAL                :: tio         ! output time and time step 
!    REAL                :: dtio        ! output time step
!    !
!    REAL                :: d           ! stability parameter  
!    REAL                :: kappa       ! heat conduction coefficient
!    REAL                :: TL, TR      ! boundary conditions
!    !
!#ifdef IMP_SOLVER
!    REAL, ALLOCATABLE   :: av(:), bv(:), cv(:), rhs(:)
!#endif
!    !------------------------------------------------------------!
!    !
!    WRITE(*,'(a)') ' | ====================================================== | '
!    WRITE(*,'(a)') ' |      Finite difference code for 1D heat equation       | '
!    WRITE(*,'(a)') ' | ====================================================== | '
!    WRITE(*,'(a)') ' | '
!    !
!    !======================================== 
!    ! SETTINGS OF THE COMPUTATION
!    !========================================
!    !
!#ifdef IMP_SOLVER
!    TestName = 'Heat1D_implicit'
!    d        = 2.0  ! for stability condition ANY value of d is OK :-)
!#else
!    TestName = 'Heat1D_explicit'
!    d        = 0.45  ! for stability condition d<0.5
!#endif    
!    !
!    IMAX   = 500
!    xL     = -1.0
!    xR     = 1.0
!    xD     = 0.0  ! location of the discontinuity
!    !
!    time   = 0.0
!    tend   = 0.01
!    tio    = 0.0
!    dtio   = 1e-2
!    CFL    = 0.9
!    dt_fix = 1e-2
!    !
!    kappa  = 1     ! heat conduction coefficient 
!    !
!    ! Boundary condition
!    TL = 100
!    TR = 50
!    !
!    !========================================
!    !
!    ! Allocate variables
!    ALLOCATE( x(IMAX)  )
!    ALLOCATE( T(IMAX)  )
!    ALLOCATE( T1(IMAX) )
!    !
!#ifdef IMP_SOLVER
!    ! Allocate also arrays for Thomas algorithm
!    ALLOCATE( av(IMAX)  )  
!    ALLOCATE( bv(IMAX)  )
!    ALLOCATE( cv(IMAX)  )
!    ALLOCATE( rhs(IMAX) )
!#endif
!    !
!    !---------------------------------------- 
!    !---------------------------------------- 
!    !
!    ! 1) Computational domain 
!    !
!    WRITE(*,'(a)') ' | Building computational domain... '
!    !    
!    dx    = (xR - xL) / REAL(IMAX-1)
!    dx2   = dx**2
!    x(1)  = xL
!    !
!    DO i = 1, (IMAX-1)
!      x(i+1) = x(i) + dx
!    ENDDO  
!    !
!    !---------------------------------------- 
!    !---------------------------------------- 
!    !
!    ! 2) Initial condition
!    !
!    WRITE(*,'(a)') ' | Assigning initial condition... '
!    !
!    DO i = 1, IMAX
!      IF(x(i).LT.xD) THEN
!        T(i) = TL
!      ELSE
!        T(i) = TR
!      ENDIF  
!    ENDDO  
!    !
!    !---------------------------------------- 
!    !---------------------------------------- 
!    !
!    ! 3) Computation: main loop in time
!    !
!#ifdef IMP_SOLVER
!    WRITE(*,'(a)') ' | '
!    WRITE(*,'(a)') ' | Implicit finite difference solver. '
!#else
!    WRITE(*,'(a)') ' | '
!    WRITE(*,'(a)') ' | Explicit finite difference solver. '
!#endif
!    WRITE(*,'(a)') ' | '
!    WRITE(*,'(a)') ' | START of the COMPUTATION '
!    
!    !
!    DO n = 1, NMAX
!      !
!      ! 3.1) Compute the time step
!      !
!      IF(time.GE.tend) THEN
!        EXIT    
!      ENDIF
!      !
!      dt = d * dx2 / kappa
!      IF((time+dt).GT.tend) THEN
!        dt  = tend - time   
!        tio = tend
!      ENDIF    
!      IF((time+dt).GT.tio) THEN
!        dt = tio - time
!      ENDIF
      !
      ! 3.2) Numerical scheme: FTCS
      !
!#ifdef IMP_SOLVER
!      ! IMPLICIT SOLVER
!      DO i = 1, IMAX
!        av(i)  = -kappa*dt/dx2
!        bv(i)  = 1+2*kappa*dt/dx2
!        cv(i)  = -kappa*dt/dx2
!        rhs(i) = T(i)
!      ENDDO  
!      ! BCs
!      rhs(1)   = rhs(1)    - av(1)   *TL           
!      rhs(IMAX)= rhs(IMAX) - cv(IMAX)*TR     
!      !
!      !CALL Thomas(T1,av,bv,cv,rhs,IMAX)
!      CALL CG(IMAX,T1,rhs,kappa,dt,dx2)
!#else
!      ! EXPLICIT SOLVER
!      DO i = 2, (IMAX-1)
!        T1(i)= T(i) + kappa*dt/(dx2) * ( T(i-1) - 2*T(i) + T(i+1) ) 
!      ENDDO
!      ! BCs
!      T1(1)    = TL   
!      T1(IMAX) = TR  
!#endif      
!      !
!      time = time + dt  ! update time
!      T    = T1         ! overwrite current solution
!      !
!      ! 3.5) Eventually plot the results
!      IF(ABS(time-tio).LT.1e-12) THEN
!        WRITE(*,'(a,f15.7)') ' |   plotting data output at time ', time
!        CALL DataOutput(n,TestName,IMAX,time,x,T)
!        tio = tio + dtio
!      ENDIF  
!      !
!    ENDDO !n
!    !
!    !---------------------------------------- 
!    !---------------------------------------- 
!    !
!    ! Empty memory
!    !
!    DEALLOCATE( x, T, T1 )
!#ifdef IMP_SOLVER
!    DEALLOCATE( av, bv, cv, rhs )
!#endif
!    !
!    WRITE(*,'(a)') ' | '
!    WRITE(*,'(a)') ' |         Finalization was successful. Bye :-)           | '
!    WRITE(*,'(a)') ' | ====================================================== | '
!    !
!END PROGRAM Heat1D
  
 
  

  

  
  