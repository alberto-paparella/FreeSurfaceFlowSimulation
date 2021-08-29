!!======================================================================================================!
!! Progetto II
!! Equazione dei flussi a superficie libera su griglie strutturate
!!======================================================================================================!
!! Università degli studi di Ferrara - Dipartimento di Matematica e Informatica
!! Corso di Algoritmi per il Calcolo Parallelo
!! Prof. Walter Boscheri
!! A.A. 2020/2021
!!======================================================================================================!
!! Studenti.
!! Alberto Paparella - 144261
!! Elena Zerbin      - 145302
!!======================================================================================================!
!! main.f90
!! This is the main file
!!======================================================================================================!
!PROGRAM FS2D
!    !==================================================================================================!
!    USE VarDef_mod  ! In this module we defined the variables we are going to use
!    IMPLICIT NONE
!    !==================================================================================================!
!    INTEGER            :: i, j, n
!    REAL               :: umax, au, av, s0, ct, cs , a , b
!    REAL, ALLOCATABLE  :: amat(:,:), amatj(:,:), bmat(:,:), cmat(:,:), cmatj(:,:)
!    !==================================================================================================!
!    WRITE(*,'(a)') ' | ================================================================== | '
!    WRITE(*,'(a)') ' |  Equazione dei flussi a superficie libera su griglie strutturate.  | '
!    WRITE(*,'(a)') ' | ================================================================== | '
!    WRITE(*,'(a)') ' | '
!    WRITE(*,'(a)') ' | Soluzione con il metodo del gradiente coniugato. '
!    !==================================================================================================!
!    ! SETTINGS OF THE COMPUTATION
!    !==================================================================================================!
!    TestName = 'Gaussian_test'
!    ! Concerning x axys
!    IMAX   = 500    ! Max index for x coordinates
!    xL     = -0.5   ! Left boundary for x coordinates
!    xR     = 0.5    ! Right boundary for y coordinates
!    ! Concerning y axys
!    JMAX   = 500    ! Max index for y coordinates
!    yD     = -0.5   ! Left boundary for y coordinates
!    yU     = 0.5    ! Right boundary for y coordinates
!    ! Initial Conditions
!    s      = 0.1    ! Initial condition
!    ! Concerning time
!    time   = 0.0    ! Starting of the computation time
!    tend   = 0.1    ! Ending of the computation time
!    tio    = 0.0    ! Time iterator for output file
!    dtio   = 5e-2   ! Time step for output file
!    CFL    = 0.9    ! CFL number for time step
!    dt_fix = 2e-2   ! Fixed time step for computation
!    nu     = 0.e-2  ! Kinematic viscosity coefficient
!    !==================================================================================================!
!    CALL Allocate_Var  ! Allocate all variables
!    !==================================================================================================!
!    ! Matrixes for assembling the tridiagonal linear system
!    !==================================================================================================!
!    ALLOCATE( amat  ( IMAX, JMAX ) )    ! Matrix a for i system
!    ALLOCATE( amatj ( IMAX, JMAX ) )    ! Matrix a for j system
!    ALLOCATE( bmat  ( IMAX, JMAX ) )    ! Matrix b for both systems
!    ALLOCATE( cmat  ( IMAX, JMAX ) )    ! Matrix c for i system
!    ALLOCATE( cmatj ( IMAX, JMAX ) )    ! Matrix c for j system
!    ALLOCATE( rhs   ( IMAX, JMAX ) )    ! Matrix rhs 
!    !==================================================================================================!
!    ! 1) Computational domain 
!    !==================================================================================================!
!    WRITE(*,'(a)') ' | Building computational domain... '
!    ! Domain on the x axys
!    dx    = (xR - xL) / REAL(IMAX)
!    dx2   = dx**2
!    x(1)  = xL
!    ! Domain on the y axys
!    dy    = (yU - yD) / REAL(IMAX)
!    dy2   = dy**2
!    y(1)  = yD
!    !==================================================================================================!
!    ! We decided to use two separates dimensions, IMAX and JMAX,
!    ! for x and y coords, to offer more flexibility
!    !==================================================================================================!
!    DO i = 1, IMAX
!        x(i+1) = x(i) + dx
!        xb(i)  = x(i) + dx/2.
!    ENDDO
!    DO j = 1, JMAX
!        y(j+1) = y(j) + dy
!        yb(j)  = y(j) + dy/2.
!    ENDDO
!    !==================================================================================================!
!    ! 2) Initial condition
!    !==================================================================================================!
!    WRITE(*,'(a)') ' | Assigning initial condition... '
!    !==================================================================================================!
!    ! 2.1) Free surface elevation (barycenters)
!    !==================================================================================================!
!    DO i = 1, IMAX
!        DO j = 1, JMAX
!            !===========================================================================================!
!            ! Gaussian profile
!            !===========================================================================================!
!            ! Note: we use xb and yb instead of x and y because x and y represents the vertex
!            ! coords for u and v, while xb and yb represents the barycenter coords for eta
!            !===========================================================================================!
!            eta(i,j) = 1.0 + EXP( (-1.0 / (2 * (s**2) ) ) * ( xb(i)**2 + yb(j)**2 ) )
!        ENDDO
!    ENDDO
!    !==================================================================================================!
!    ! 2.1) Velocity, bottom and total water depth (interfaces)
!    !==================================================================================================!
!    ! Note: u dimensions (IMAX + 1) * JMAX, while v dimensions are IMAX * (JMAX + 1)
!    !==================================================================================================!
!    DO i = 1, IMAX + 1
!        DO j = 1, JMAX
!            u(i, j)  = 0.0  ! Fluid at rest
!            bu(i, j) = 0.0  ! Zero bottom elevation
!        ENDDO
!    ENDDO
!    DO i = 1, IMAX
!        DO j = 1, JMAX + 1
!            v(i, j)  = 0.0  ! Fluid at rest
!            bv(i, j) = 0.0  ! Zero bottom elevation
!        ENDDO
!    ENDDO
!    !==================================================================================================!
!    ! Total water depth
!    !==================================================================================================!
!    DO j = 1, JMAX
!        ! Initialize first and last columns for Hu
!        Hu( 1,      j ) = MAX( 0.0, bu( 1,      j ) + eta( 1,    j ) )
!        Hu( IMAX+1, j ) = MAX( 0.0, bu( IMAX+1, j ) + eta( IMAX, j ) )
!    ENDDO
!    DO i = 1, IMAX
!        ! Initialize first and last rows for Hv
!        Hv( i, 1      ) = MAX( 0.0, bv( i, 1      ) + eta( i, 1    ) )
!        Hv( i, JMAX+1 ) = MAX( 0.0, bv( i, JMAX+1 ) + eta( i, JMAX ) )
!    ENDDO
!    DO i = 2, IMAX
!        DO j = 1, JMAX
!            Hu( i , j ) = MAXVAL( (/ 0.0, bu( i, j ) + eta( i - 1, j ), bu( i, j )+ eta( i, j ) /) )
!        ENDDO
!    ENDDO
!    DO i = 1, IMAX
!        DO j = 2, JMAX
!            Hv( i , j ) = MAXVAL( (/ 0.0, bv( i, j ) + eta( i, j - 1 ), bv( i, j )+ eta( i, j ) /) )
!        ENDDO
!    ENDDO
!    !==================================================================================================!
!    ! Plot initial condition
!    !==================================================================================================!
!    CALL DataOutput(0)
!    tio = tio + dtio
!    !==================================================================================================!
!    ! 3) Computation: main loop in time
!    !==================================================================================================!
!    WRITE(*,'(a)') ' | '
!    WRITE(*,'(a)') ' | START of the COMPUTATION '
!    !==================================================================================================!
!    DO n = 1, NMAX
!        !==============================================================================================!
!        ! 3.1) Compute the time step
!        !==============================================================================================!
!        IF( time.GE.tend ) THEN
!            EXIT    
!        ENDIF
!        !==============================================================================================!
!        umax = MAXVAL( ABS(u) )             ! Semi-implicit time step
!        dt   = MIN( dt_fix, CFL / ( ( umax + 1e-14 ) / dx + 2. * nu / dx2 ) )
!        dt2   = dt**2
!        IF( ( time + dt ).GT.tend ) THEN
!            dt  = tend - time   
!            tio = tend
!        ENDIF    
!        IF( ( time + dt ).GT.tio ) THEN
!            dt = tio - time
!        ENDIF
!        !==============================================================================================!
!        ! 3.2) Compute the operator Fu
!        !==============================================================================================!        
!        ! Fu = u    ! Case without considering Fu
!        ! BC: no-slip wall
!        ! First row and last row initialized to zeros
!        Fu( 1      , : ) = u(1,:)
!        Fu( IMAX+1 , : ) = u(IMAX+1,:)
!        DO i = 2, IMAX
!            DO j = 1, JMAX
!                au = ABS( u(i,j) )
!                ! Explicit upwind
!                ! In our case, nu = 0
!                !Fu(i,j) = ( 1. - dt * ( au/dx + 2.*nu/dx2 ) ) * u(i,j)   &     ! q^N
!                !        + dt * ( nu/dx2 + (au-u(i,j))/(2.*dx) ) * u(i+1,j) &   ! a^+(...)
!                !        + dt * ( nu/dx2 + (au+u(i,j))/(2.*dx) ) * u(i-1,j)     ! a^-(...)
!                Fu(i,j) = ( 1. - dt * ( au/dx ) ) * u(i,j)   &
!                       + dt * ( (au-u(i,j))/(2.*dx) ) * u(i+1,j) &
!                       + dt * ( (au+u(i,j))/(2.*dx) ) * u(i-1,j)
!            ENDDO
!        ENDDO
!        !==============================================================================================!
!        ! 3.3) Compute the operator Fv
!        !==============================================================================================!
!        ! Fv = v    ! Case without considering Fv
!        ! BC: no-slip wall
!        ! First column and last column initialized to zeros
!        Fv( : , 1      ) = v(:,1)
!        Fv( : , JMAX+1 ) = v(:,JMAX+1)
!        DO i = 1, IMAX
!            DO j = 2, JMAX
!                av = ABS( v(i,j) )
!                ! Explicit upwind
!                Fv(i,j) = ( 1. - dt * ( av/dx + 2.*nu/dy2 ) ) * v(i,j)     &
!                        + dt * ( nu/dy2 + (av-v(i,j))/(2.*dy) ) * v(i,j+1) &
!                        + dt * ( nu/dy2 + (av+v(i,j))/(2.*dy) ) * v(i,j-1)
!            ENDDO
!        ENDDO
!        !==============================================================================================!
!        ! 3.4) Solve the free surface equation
!        !==============================================================================================!
!        ! CONJUGATE GRADIENT METHOD
!        !==============================================================================================!
!        DO i = 1, IMAX
!            DO j = 1, JMAX
!                rhs(i,j) = eta(i,j) - dt/dx * ( Hu(i+1,j)*Fu(i+1,j) - Hu(i,j)*Fu(i,j)) &
!                                    - dt/dy * ( Hv(i,j+1)*Fv(i,j+1) - Hv(i,j)*Fv(i,j))
!            ENDDO
!        ENDDO
!        ! Here, we suppose that IMAX and JMAX are the same, at least for the moment
!        CALL CG(IMAX,eta,rhs)
!        !==============================================================================================!
!        ! 3.5) Update the velocity (momentum equation)
!        !==============================================================================================!
!        ct = g*dt/dx    ! Temporary coefficient for x
!        cs = g*dt/dy    ! Temporary coefficient for y
!        !==============================================================================================!
!        u(1,:)      = Fu(1,:) 
!        u(IMAX+1,:) = Fu(IMAX+1,:)
!        DO i = 2, IMAX
!            DO j = 1, JMAX
!                u(i,j) = Fu(i,j) - ct * ( eta(i,j) - eta(i-1,j) )
!            ENDDO
!        ENDDO
!        !==============================================================================================!
!        v(:,1)      = Fv(:,1)
!        v(:,JMAX+1) = Fv(:,JMAX+1)
!        DO i = 1, IMAX
!            DO j = 2, JMAX
!                v(i,j) = Fv(i,j) - cs * ( eta(i,j) - eta(i,j-1) )
!            ENDDO
!        ENDDO
!        !==============================================================================================!
!        ! 3.6) Update total water depth
!        !==============================================================================================!
!        DO j = 1, JMAX
!            ! Initialize first and last columns for Hu
!            Hu( 1,      j ) = MAX( 0.0, bu( 1,      j ) + eta( 1,    j ) )
!            Hu( IMAX+1, j ) = MAX( 0.0, bu( IMAX+1, j ) + eta( IMAX, j ) )
!        ENDDO
!        DO i = 1, IMAX
!            ! Initialize first and last rows for Hv
!            Hv( i, 1      ) = MAX( 0.0, bv( i, 1      ) + eta( i, 1    ) )
!            Hv( i, JMAX+1 ) = MAX( 0.0, bv( i, JMAX+1 ) + eta( i, JMAX ) )
!        ENDDO
!        DO i = 2, IMAX
!            DO j = 1, JMAX
!                Hu( i , j ) = MAXVAL( (/ 0.0, bu( i, j ) + eta( i - 1, j ), bu( i, j )+ eta( i, j ) /) )
!            ENDDO
!        ENDDO
!        DO i = 1, IMAX
!            DO j = 2, JMAX
!                Hv( i , j ) = MAXVAL( (/ 0.0, bv( i, j ) + eta( i, j - 1 ), bv( i, j )+ eta( i, j ) /) )
!            ENDDO
!        ENDDO
!        !==========================================================================================!
!        time = time + dt  ! Update time
!        !==========================================================================================!
!        ! 3.7) Eventually plot the results
!        !CALL DataOutput(n)  ! DEBUG
!        !STOP                ! DEBUG
!        IF(ABS(time-tio).LT.1e-12) THEN
!            WRITE(*,'(a,f15.7)') ' |   plotting data output at time ', time
!            CALL DataOutput(n)
!            tio = tio + dtio
!        ENDIF  
!    !==============================================================================================!
!    ENDDO ! n cycle
!    !==============================================================================================!
!    WRITE(*,'(a)') ' | END of the COMPUTATION '
!    !==============================================================================================!
!    ! Empty memory
!    !==============================================================================================!
!    CALL Deallocate_Var
!    !==============================================================================================!
!    DEALLOCATE( amat,bmat,cmat, amatj,cmatj )
!    !==============================================================================================!
!    WRITE(*,'(a)') ' | '
!    WRITE(*,'(a)') ' |         Finalization was successful. Bye :-)           | '
!    WRITE(*,'(a)') ' | ====================================================== | '
!    !==============================================================================================!    
!    PAUSE
!!======================================================================================================!
!END PROGRAM FS2D
!!======================================================================================================!
!SUBROUTINE matop2D(Ap,p,N)
!    !==================================================================================================!
!    USE VarDef_mod
!    IMPLICIT NONE
!    !==================================================================================================!
!    INTEGER  :: N 
!    REAL     :: Ap(N,N), p(N,N)
!    !==================================================================================================!
!    INTEGER  :: i, j
!    REAL     :: cx, cy
!    !==================================================================================================!
!    cx = g*dt**2/dx2    ! Temporary coefficient
!    cy = g*dt**2/dy2
!    !==================================================================================================!	
!	! 1) Diagonal term
!    !==================================================================================================!	
!	Ap = p
!    !==================================================================================================!	
!	! 2) Fluxes in x-direction
!    !==================================================================================================!	
!	DO j = 1, N
!		DO i = 1, N
!			IF(i.EQ.1) THEN
!				Ap(i,j) = Ap(i,j) - cx * ( Hu(i+1,j)*(p(i+1,j)-p(i,j)) - 0.0 )
!			ELSEIF(i.EQ.N) THEN  
!				Ap(i,j) = Ap(i,j) - cx * ( 0.0 - Hu(i,j)*(p(i,j)-p(i-1,j)) )
!			ELSE  
!				Ap(i,j) = Ap(i,j) - cx * ( Hu(i+1,j)*(p(i+1,j)-p(i,j)) - Hu(i,j)*(p(i,j)-p(i-1,j)) )
!			ENDIF 
!		ENDDO
!    ENDDO
!    !==================================================================================================!	
!    ! 3) Fluxes in y-direction
!    !==================================================================================================!	
!	DO j = 1, N
!		DO i = 1, N
!			IF(j.EQ.1) THEN
!				Ap(i,j) = Ap(i,j) - cy * ( Hv(i,j+1)*(p(i,j+1)-p(i,j)) - 0.0 )
!			ELSEIF(j.EQ.N) THEN  
!				Ap(i,j) = Ap(i,j) - cy * ( 0.0 - Hv(i,j)*(p(i,j)-p(i,j-1)) )
!			ELSE  
!				Ap(i,j) = Ap(i,j) - cy * ( Hv(i,j+1)*(p(i,j+1)-p(i,j)) - Hv(i,j)*(p(i,j)-p(i,j-1)) )
!			ENDIF 
!		ENDDO
!    ENDDO    
!!======================================================================================================!
!    END SUBROUTINE matop2D
    
   
    

!======================================================================================================!
!PARALLELIZATION 
!======================================================================================================!

!#define PARALLEL 
PROGRAM FS2D
    !==================================================================================================!
    USE VarDef_mod  ! In this module we defined the variables we are going to use
    IMPLICIT NONE
#ifdef PARALLEL
	INCLUDE 'mpif.h'
#endif    
    !==================================================================================================!
    INTEGER            :: i, j, n    
    REAL               :: umax, au, av, s0, ct, cs , a , b
    REAL, ALLOCATABLE  :: amat(:,:), amatj(:,:), bmat(:,:), cmat(:,:), cmatj(:,:)
#ifdef PARALLEL   
    INTEGER             :: LCPU, RCPU, MsgLength, nMsg 
    REAL                :: send_messageL, send_messageR, recv_messageL, recv_messageR
    INTEGER             :: send_request(2), recv_request(2) 
    INTEGER             :: send_status_list(MPI_STATUS_SIZE,2),recv_status_list(MPI_STATUS_SIZE,2)
#endif    
    !==================================================================================================!
    TYPE tMPI
       INTEGER :: myrank, nCPU, iErr 
       INTEGER :: AUTO_REAL
       INTEGER :: nElem, istart, iend, jstart, jend 
    END TYPE tMPI 
    TYPE(tMPI) :: MPI
    REAL       :: realtest
    REAL       :: WCT1, WCT2           ! Wall clock times 
    WCT1 = 0. 
    WCT2 = 0.
#ifdef PARALLEL
    CALL MPI_INIT(MPI%iErr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, MPI%myrank, MPI%iErr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, MPI%nCPU,   MPI%iErr)
    !==================================================================================================!
    ! realtest
    ! We differentiate between single-precision and double-precision values
    !==================================================================================================!
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
         
         PRINT *, 'Stopping here.'
         PAUSE
         
         STOP 
    ENDIF
    
    !PRINT *, 'YUPPIE!'
    !PAUSE
    
    !==================================================================================================!
    MPI%nElem  = (IMAX*JMAX)/MPI%nCPU
    MPI%istart = 1 + MPI%myrank*MPI%nElem 
    MPI%iend  = MPI%istart + MPI%nElem - 1 
    MPI%jstart = 1 + MPI%myrank*MPI%nElem 
    MPI%jend  = MPI%jstart + MPI%nElem - 1
    !==================================================================================================!    
#else
    MPI%myrank = 0              ! If the code is not compiled in parallel, then there is only one CPU ...
    MPI%nCPU   = 1              ! ...
    MPI%AUTO_REAL = -1          ! ... and the data type MPI_AUTO_REAL is not needed 
    MPI%istart = 1
    MPI%iend   = IMAX  
    MPI%jstart = 1
    MPI%jend   = JMAX 
#endif
    !==================================================================================================!    
    IF(MPI%myrank.EQ.0) THEN
    WRITE(*,'(a)') ' | ================================================================== | '
    WRITE(*,'(a)') ' |  Equazione dei flussi a superficie libera su griglie strutturate.  | '
    WRITE(*,'(a)') ' | ================================================================== | '
    WRITE(*,'(a)') ' | ' 
    WRITE(*,'(a)') ' | Soluzione con il metodo del gradiente coniugato. '
    WRITE(*,'(a,i10)') '         Total number of CPUs used : ', MPI%nCPU       
    WRITE(*,'(a)') ' | ================================================================== | '
    WRITE(*,'(a)') ' '
    ENDIF   
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
#ifdef PARALLEL
    ALLOCATE( eta(MPI%istart-1:MPI%iend+1, MPI%jstart-1:MPI%jend+1) )
#else
    ALLOCATE( eta( IMAX, JMAX )                           )
    eta = 0.
#endif

    !PRINT *, 'YUPPIE!'
    !PAUSE
   
    !==================================================================================================!
    ! 1) Computational domain 
    !==================================================================================================!
    IF(MPI%myrank.EQ.0) THEN
      WRITE(*,'(a)') '  Building computational domain... '
    ENDIF 
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
     IF(MPI%myrank.EQ.0) THEN
      WRITE(*,'(a)') '  Assigning initial condition... '
     ENDIF
   
    !==================================================================================================!
    ! 2.1) Free surface elevation (barycenters)
    !==================================================================================================!
#ifdef PARALLEL
    DO i= MPI%istart, MPI%iend
        DO j = MPI%jstart, MPI%jend
            !===========================================================================================!
            ! Gaussian profile
            !===========================================================================================!
            ! Note: we use xb and yb instead of x and y because x and y represents the vertex
            ! coords for u and v, while xb and yb represents the barycenter coords for eta
            !===========================================================================================!
            eta(i,j) = 1.0 + EXP( (-1.0 / (2 * (s**2) ) ) * ( xb(i)**2 + yb(j)**2 ) )
        ENDDO
    ENDDO
#else
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
#endif
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
#ifdef PARALLEL
    !==================================================================================================!
    ! Total water depth
    !==================================================================================================!
    DO j = MPI%jstart, MPI%jend
        ! Initialize first and last columns for Hu
        Hu( 1,      j ) = MAX( 0.0, bu( 1,      j ) + eta( 1,    j ) )
        Hu( IMAX+1, j ) = MAX( 0.0, bu( IMAX+1, j ) + eta( IMAX, j ) )
    ENDDO
    DO i= MPI%istart, MPI%iend
        ! Initialize first and last rows for Hv
        Hv( i, 1      ) = MAX( 0.0, bv( i, 1      ) + eta( i, 1    ) )
        Hv( i, JMAX+1 ) = MAX( 0.0, bv( i, JMAX+1 ) + eta( i, JMAX ) )
    ENDDO
    DO i = MPI%istart+1, MPI%iend-1
        DO j = MPI%jstart, MPI%jend
            Hu( i , j ) = MAXVAL( (/ 0.0, bu( i, j ) + eta( i - 1, j ), bu( i, j )+ eta( i, j ) /) )
        ENDDO
    ENDDO
    DO i = MPI%istart, MPI%iend
        DO j = MPI%jstart+1, MPI%jend-1
            Hv( i , j ) = MAXVAL( (/ 0.0, bv( i, j ) + eta( i, j - 1 ), bv( i, j )+ eta( i, j ) /) )
        ENDDO
    ENDDO
#else    
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
#endif    
    !==================================================================================================!
    ! Plot initial condition
    !==================================================================================================!
#ifdef PARALLEL    
    CALL DataOutput(0, MPI%istart, MPI%iend, MPI%jstart, MPI%jend, MPI%myrank)
#else
    CALL DataOutput(0)
#endif    
    tio = tio + dtio    
    !==================================================================================================!
    ! 3) Computation: main loop in time
    !==================================================================================================!
    IF(MPI%myrank.EQ.0) THEN
        WRITE(*,'(a)') ' | '
        WRITE(*,'(a)') ' | START of the COMPUTATION ' 
    ENDIF  
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
    !==================================================================================================!
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
        Fu( 1      , : ) = u(1,:)
        Fu( IMAX+1 , : ) = u(IMAX+1,:)
        !==============================================================================================!
        ! 3.3) Compute the operator Fv
        !==============================================================================================! 
        Fv( : , 1      ) = v(:,1)
        Fv( : , JMAX+1 ) = v(:,JMAX+1)
#ifdef PARALLEL
        DO i = MPI%istart, MPI%iend
          DO j = MPI%jstart+1, MPI%jend  
           Fu(i,j) = ( 1. - dt * ( au/dx ) ) * u(i,j)   &
                       + dt * ( (au-u(i,j))/(2.*dx) ) * u(i+1,j) &
                       + dt * ( (au+u(i,j))/(2.*dx) ) * u(i-1,j)
           ENDDO
        ENDDO  
        DO i = MPI%istart+1, MPI%iend  
          DO j = MPI%jstart, MPI%jend
            av = ABS( v(i,j) )
            ! Explicit upwind
            Fv(i,j) = ( 1. - dt * ( av/dx + 2.*nu/dy2 ) ) * v(i,j)     &
                    + dt * ( nu/dy2 + (av-v(i,j))/(2.*dy) ) * v(i,j+1) &
                    + dt * ( nu/dy2 + (av+v(i,j))/(2.*dy) ) * v(i,j-1)
        ENDDO
      ENDDO
      
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
      send_messageL = eta(MPI%istart, MPI%jstart) 
      CALL MPI_ISEND(send_messageL, MsgLength, MPI%AUTO_REAL, LCPU, 1, MPI_COMM_WORLD, send_request(1), MPI%iErr)
      ! send your rightmost state to your right neighbor CPU => post a send request  
      send_messageR = eta(MPI%iend, MPI%jend)
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
      !==============================================================================================!
        ! 3.4) Solve the free surface equation
        !==============================================================================================!
        ! CONJUGATE GRADIENT METHOD
        !==============================================================================================!  
        DO i = MPI%istart, MPI%iend 
          DO j = MPI%jstart, MPI%jend
                rhs(i,j) = eta(i,j) - dt/dx * ( Hu(i+1,j)*Fu(i+1,j) - Hu(i,j)*Fu(i,j)) &
                                    - dt/dy * ( Hv(i,j+1)*Fv(i,j+1) - Hv(i,j)*Fv(i,j))
            ENDDO
        ENDDO
        CALL CG(IMAX,eta,rhs)
        
        !  
        ! Wait until all communication has finished
        ! 
        nMsg = 2 
        CALL MPI_WAITALL(nMsg,send_request,send_status_list,MPI%ierr)
        CALL MPI_WAITALL(nMsg,recv_request,recv_status_list,MPI%ierr)
        eta(MPI%istart-1, MPI%jstart-1)   = recv_messageL 
        eta(MPI%iend+1, MPI%jend+1)       = recv_messageR
        !==============================================================================================!
        ! 3.5) Update the velocity (momentum equation)
        !==============================================================================================!
        ct = g*dt/dx    ! Temporary coefficient for x
        cs = g*dt/dy    ! Temporary coefficient for y
        !==============================================================================================!
        u(1,:)      = Fu(1,:) 
        u(IMAX+1,:) = Fu(IMAX+1,:)
        v(:,1)      = Fv(:,1)
        v(:,JMAX+1) = Fv(:,JMAX+1)
        
        DO i = MPI%istart +1, MPI%iend
            DO j =MPI%jstart, MPI%jend
                u(i,j) = Fu(i,j) - ct * ( eta(i,j) - eta(i-1,j) )
            ENDDO
        ENDDO
        !==============================================================================================!
    
        DO i = MPI%istart, MPI%iend
            DO j =MPI%jstart+1, MPI%jend
                v(i,j) = Fv(i,j) - cs * ( eta(i,j) - eta(i,j-1) )
            ENDDO
        ENDDO
        
        !==============================================================================================!
        ! 3.6) Update total water depth
        !==============================================================================================!
    
        DO j = MPI%jstart, MPI%jend, (MPI%jend - MPI%jstart)
        ! Initialize first and last columns for Hu
            Hu( 1,      j ) = MAX( 0.0, bu( 1,      j ) + eta( 1,    j ) )
            Hu( IMAX+1, j ) = MAX( 0.0, bu( IMAX+1, j ) + eta( IMAX, j ) )
        ENDDO
        DO i= MPI%istart, MPI%iend , (MPI%iend - MPI%istart)
            ! Initialize first and last rows for Hv
            Hv( i, 1      ) = MAX( 0.0, bv( i, 1      ) + eta( i, 1    ) )
            Hv( i, JMAX+1 ) = MAX( 0.0, bv( i, JMAX+1 ) + eta( i, JMAX ) )
        ENDDO
        DO i = MPI%istart +1, MPI%iend , (MPI%iend - MPI%istart)
            DO j = MPI%jstart, MPI%jend , (MPI%jend - MPI%jstart)
                Hu( i , j ) = MAXVAL( (/ 0.0, bu( i, j ) + eta( i - 1, j ), bu( i, j )+ eta( i, j ) /) )
            ENDDO
        ENDDO
        DO i = MPI%istart, MPI%iend , (MPI%iend - MPI%istart)
            DO j = MPI%jstart+1, MPI%jend , (MPI%jend - MPI%jstart)
                Hv( i , j ) = MAXVAL( (/ 0.0, bv( i, j ) + eta( i, j - 1 ), bv( i, j )+ eta( i, j ) /) )
            ENDDO
        ENDDO
#else      
    
    DO i = 1, IMAX
        DO j = 2, JMAX
            av = ABS( v(i,j) )
            ! Explicit upwind
            Fv(i,j) = ( 1. - dt * ( av/dx + 2.*nu/dy2 ) ) * v(i,j)     &
                    + dt * ( nu/dy2 + (av-v(i,j))/(2.*dy) ) * v(i,j+1) &
                    + dt * ( nu/dy2 + (av+v(i,j))/(2.*dy) ) * v(i,j-1)
        ENDDO
    ENDDO
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
              ! 

        CALL CG(IMAX,eta,rhs)
        
        !==============================================================================================!
        ! 3.5) Update the velocity (momentum equation)
        !==============================================================================================!
        ct = g*dt/dx    ! Temporary coefficient for x
        cs = g*dt/dy    ! Temporary coefficient for y
        !==============================================================================================!
        u(1,:)      = Fu(1,:) 
        u(IMAX+1,:) = Fu(IMAX+1,:)
        v(:,1)      = Fv(:,1)
        v(:,JMAX+1) = Fv(:,JMAX+1)

        DO i = 2, IMAX
            DO j = 1, JMAX
                u(i,j) = Fu(i,j) - ct * ( eta(i,j) - eta(i-1,j) )
            ENDDO
        ENDDO
        !==============================================================================================!
    
        DO i = 1, IMAX
            DO j = 2, JMAX
                v(i,j) = Fv(i,j) - cs * ( eta(i,j) - eta(i,j-1) )
            ENDDO
        ENDDO      
    
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
    
#endif  
    !==========================================================================================!
    time = time + dt  ! Update time
    !==========================================================================================!
    ! 3.7) Eventually plot the results
    !CALL DataOutput(n)  ! DEBUG
    !STOP                ! DEBUG
    IF(ABS(time-tio).LT.1e-12) THEN
      IF(MPI%myrank.EQ.0) THEN
          WRITE(*,'(a,f15.7)') '   => plotting data output at time ', time
           CALL DataOutput(n, MPI%istart, MPI%iend, MPI%jstart, MPI%jend, MPI%myrank)
        ENDIF 
    
#ifdef PARALLEL
    CALL MPI_ALLGATHER(rhs(MPI%istart:MPI%iend, MPI%jstart:MPI%jend),MPI%nElem,MPI%AUTO_REAL,eta,MPI%nElem,MPI%AUTO_REAL,MPI_COMM_WORLD,MPI%iErr)
    IF(MPI%myrank.EQ.0) THEN
    CALL DataOutput(n,MPI%istart,MPI%iend,MPI%jstart,MPI%jend,MPI%myrank)
    ENDIF
     
#else
 CALL DataOutput(n,MPI%istart,MPI%iend,MPI%jstart,MPI%jend,MPI%myrank)
#endif
 tio = tio + dtio
  ENDIF 
    !==============================================================================================!
ENDDO ! n cycle
    !==============================================================================================!
    WRITE(*,'(a)') ' | END of the COMPUTATION '

    !
#ifdef PARALLEL
    WCT2 = MPI_WTime() 
    CALL MPI_FINALIZE(MPI%iErr) 
#else
    CALL CPU_TIME(WCT2) 
#endif  
    !==============================================================================================!
    ! Empty memory
    !==============================================================================================!
    CALL Deallocate_Var
    !==============================================================================================!
    DEALLOCATE( amat,bmat,cmat, amatj,cmatj )
    !==============================================================================================!
    IF(MPI%myrank.EQ.0) THEN
    WRITE(*,'(a)') '| ====================================================== |  '
    WRITE(*,'(a,f15.7)'), '         Total wallclock time needed : ', WCT2-WCT1
    WRITE(*,'(a)') ' | ====================================================== | ' 
    WRITE(*,'(a)') ' |         Finalization was successful. Bye :-)           | '
    WRITE(*,'(a)') ' | ====================================================== | '
    ENDIF  
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
  