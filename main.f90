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
!#define PARALLEL   ! Already defined in project properties
!======================================================================================================!
PROGRAM FS2D
    !==================================================================================================!
    USE VarDef_mod  ! In this module we defined the variables we are going to use
    IMPLICIT NONE
#ifdef PARALLEL
	INCLUDE 'mpif.h'
#endif    
    !==================================================================================================!
    INTEGER            :: i, j, n   ! indices
    REAL               :: umax, au, av, s0, ct, cs , a , b
#ifdef PARALLEL 
    INTEGER             :: LCPU, RCPU, MsgLength, nMsg
    REAL, ALLOCATABLE   :: send_messageL(:), send_messageR(:), recv_messageL(:), recv_messageR(:)
    INTEGER             :: send_request(2), recv_request(2) 
    INTEGER             :: send_status_list(MPI_STATUS_SIZE,2),recv_status_list(MPI_STATUS_SIZE,2)
#endif    
    !==================================================================================================!
    TYPE tMPI
        INTEGER                 :: myrank, nCPU, iErr 
        INTEGER                 :: AUTO_REAL
        INTEGER                 :: nElem, istart, iend, jstart, jend
    END TYPE tMPI 
    TYPE(tMPI) :: MPI
    REAL       :: realtest
    REAL       :: WCT1, WCT2           ! Wall clock times 
    WCT1 = 0. 
    WCT2 = 0.   
    !==================================================================================================!
    ! SETTINGS OF THE COMPUTATION (PART 1)
    !==================================================================================================!
    IMAX   = 200    ! Max index for x coordinates
    JMAX   = 200   ! Max index for y coordinates
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
         STOP 
    ENDIF
    ! Distribute elements among the processors
    ! Every process will have a sub-matrix of dimension IMAX * (JMAX/nCPU)
    MPI%nElem   = IMAX/MPI%nCPU
    MPI%jstart  = 1 + MPI%myrank*MPI%nElem
    MPI%jend    = MPI%jstart + MPI%nElem - 1
#else
    MPI%myrank      = 0             ! If the code is not compiled in parallel, then there is only one CPU ...
    MPI%nCPU        = 1             ! ...
    MPI%AUTO_REAL   = -1            ! ... and the data type MPI_AUTO_REAL is not needed
    MPI%jstart      = 1
    MPI%jend        = JMAX  
#endif
    MPI%istart      = 1
    MPI%iend        = IMAX
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
    ! SETTINGS OF THE COMPUTATION (PART 2)
    !==================================================================================================!
    TestName = 'Gaussian_test'
    ! Concerning x axys
    xL     = -0.5   ! Left boundary for x coordinates
    xR     = 0.5    ! Right boundary for y coordinates
    ! Concerning y axys
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
    ! Maybe rhs must be of the sime dimensions of eta? (for parallelization)
    ALLOCATE( rhs   ( IMAX, JMAX ) )    ! Matrix rhs for assembling the tridiagonal linear system
#ifdef PARALLEL
    ! Distributed matrices
    ALLOCATE( eta ( IMAX,   MPI%jstart-1:MPI%jend+1 ) )    
    ALLOCATE( Hu  ( IMAX+1, MPI%jstart-1:MPI%jend+1 ) )
    ALLOCATE( Hv  ( IMAX,   MPI%jstart-1:MPI%jend+2 ) )
    ALLOCATE( bu  ( IMAX+1, MPI%jstart-1:MPI%jend+1 ) )
    ALLOCATE( bv  ( IMAX,   MPI%jstart-1:MPI%jend+2 ) )    
    ALLOCATE( u   ( IMAX+1, MPI%jstart-1:MPI%jend+1 ) )
    ALLOCATE( Fu  ( IMAX+1, MPI%jstart-1:MPI%jend+1 ) )
    ALLOCATE( v   ( IMAX,   MPI%jstart-1:MPI%jend+2 ) )
    ALLOCATE( Fv  ( IMAX,   MPI%jstart-1:MPI%jend+2 ) )
    ! Buffers to gather from all processors
    ALLOCATE( eta1( IMAX,   JMAX   ) )
    ALLOCATE( u1  ( IMAX+1, JMAX   ) )
    ALLOCATE( v1  ( IMAX,   JMAX+1 ) )
    eta = 0.
    eta1 = 0.
#else
    ALLOCATE( eta( IMAX, JMAX )     )    
    ALLOCATE( Hu(  IMAX+1, JMAX   ) )
    ALLOCATE( Hv(  IMAX,   JMAX+1 ) )
    ALLOCATE( bu( IMAX+1, JMAX   ) )
    ALLOCATE( bv( IMAX,   JMAX+1 ) )
    ALLOCATE( u(   IMAX+1, JMAX   ), Fu( IMAX+1, JMAX   ) )
    ALLOCATE( v(   IMAX,   JMAX+1 ), Fv( IMAX,   JMAX+1 ) )
    eta = 0.   
#endif   
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
    !
    ! We do some useful work here, e.g. update all internal elements of the domain, since we do non-blocking communication ! 
    !
    ! The main finite difference scheme is IDENTICAL to the serial code ! 
    ! This is the simplicify and the beauty of MPI :-) 
    !
    DO i = MPI%istart, MPI%iend
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
#ifdef PARALLEL
    DO i = MPI%istart, MPI%iend + 1
        DO j = MPI%jstart, MPI%jend
            u(i, j)  = 0.0  ! Fluid at rest
            bu(i, j) = 0.0  ! Zero bottom elevation
        ENDDO
    ENDDO
    DO i = MPI%istart, MPI%iend
        DO j = MPI%jstart, MPI%jend + 1
            v(i, j)  = 0.0  ! Fluid at rest
            bv(i, j) = 0.0  ! Zero bottom elevation
        ENDDO
    ENDDO
#else    
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
#endif    
    !==================================================================================================!
    ! Total water depth
    !==================================================================================================!
#ifdef PARALLEL
    DO j = MPI%jstart, MPI%jend
        ! Initialize first and last columns for Hu
        Hu( MPI%istart, j ) = MAX( 0.0, bu( MPI%istart, j ) + eta( MPI%istart, j ) )
        Hu( MPI%iend+1, j ) = MAX( 0.0, bu( MPI%iend+1, j ) + eta( MPI%iend,   j ) )
    ENDDO    
    DO i = MPI%istart, MPI%iend
        ! Initialize first and last rows for Hv
        Hv( i, MPI%jstart ) = MAX( 0.0, bv( i, MPI%jstart ) + eta( i, MPI%jstart ) )
        Hv( i, MPI%jend+1 ) = MAX( 0.0, bv( i, MPI%jend+1 ) + eta( i, MPI%jend   ) )
    ENDDO
    DO i = MPI%istart+1, MPI%iend
        DO j = MPI%jstart, MPI%jend
            Hu( i , j ) = MAXVAL( (/ 0.0, bu( i, j ) + eta( i - 1, j ), bu( i, j )+ eta( i, j ) /) )
        ENDDO
    ENDDO
    DO i = MPI%istart, MPI%iend
        DO j = MPI%jstart+1, MPI%jend
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
    ! Collect all the data from each MPI subdomain
    ! Collect from eta into eta1
    DO i = MPI%istart, MPI%iend
        CALL MPI_ALLGATHER( eta( i, MPI%jstart:MPI%jend ), MPI%nElem, MPI%AUTO_REAL, eta1(i,:), &
                            MPI%nElem, MPI%AUTO_REAL, MPI_COMM_WORLD, MPI%iErr)
        CALL MPI_ALLGATHER( u( i, MPI%jstart:MPI%jend ), MPI%nElem, MPI%AUTO_REAL, u1(i,:), &
                            MPI%nElem, MPI%AUTO_REAL, MPI_COMM_WORLD, MPI%iErr)
        ! v has one more column
        CALL MPI_ALLGATHER( v( i, MPI%jstart:MPI%jend+1 ), MPI%nElem+1, MPI%AUTO_REAL, v1(i,:), &
                            MPI%nElem+1, MPI%AUTO_REAL, MPI_COMM_WORLD, MPI%iErr)
    ENDDO
    ! u has one more line
    CALL MPI_ALLGATHER( u( MPI%iend+1, MPI%jstart:MPI%jend ), MPI%nElem, MPI%AUTO_REAL, u1(MPI%iend+1,:), &
                            MPI%nElem, MPI%AUTO_REAL, MPI_COMM_WORLD, MPI%iErr)
    ! Finally, update eta, u and v
    eta = eta1
    u   = u1
    v   = v1
    !==================================================================================================!
    IF(MPI%myrank.EQ.0) THEN    
        WRITE(*,'(a)') ' | '
        WRITE(*,'(a)') ' | Plotting initial condition ' 
        CALL DataOutput(0, 1, IMAX, 1, JMAX, MPI%myrank)
    ENDIF
#else    
    WRITE(*,'(a)') ' | '
    WRITE(*,'(a)') ' | Plotting initial condition ' 
    CALL DataOutput(0, MPI%istart, MPI%iend, MPI%jstart, MPI%jend, 1)
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
#ifdef PARALLEL
        ! First row and last row initialized to zeros
        Fu( MPI%istart, : ) = u( MPI%istart, : )
        Fu( MPI%iend+1, : ) = u( MPI%iend+1, : )
        DO i = MPI%istart + 1, MPI%iend
            DO j = MPI%jstart, MPI%jend
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
#else
        ! First row and last row initialized to zeros
        Fu( 1      , : ) = u(1,:)
        Fu( IMAX+1 , : ) = u(IMAX+1,:)
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
#endif
        !==============================================================================================!
        ! 3.3) Compute the operator Fv
        !==============================================================================================!
#ifdef PARALLEL
        ! First column and last column initialized to zeros
        Fv( :, MPI%jstart ) = v( :, MPI%jstart )
        Fv( :, MPI%jend+1 ) = v( :, MPI%jend+1 )
        DO i = MPI%istart, MPI%iend
            DO j = MPI%jstart+1, MPI%jend
                av = ABS( v(i,j) )
                ! Explicit upwind
                Fv(i,j) = ( 1. - dt * ( av/dx + 2.*nu/dy2 ) ) * v(i,j)     &
                        + dt * ( nu/dy2 + (av-v(i,j))/(2.*dy) ) * v(i,j+1) &
                        + dt * ( nu/dy2 + (av+v(i,j))/(2.*dy) ) * v(i,j-1)
            ENDDO
        ENDDO
#else
        ! First column and last column initialized to zeros
        Fv( : , 1      ) = v(:,1)
        Fv( : , JMAX+1 ) = v(:,JMAX+1)
        DO i = 1, IMAX
            DO j = 2, JMAX
                av = ABS( v(i,j) )
                ! Explicit upwind
                Fv(i,j) = ( 1. - dt * ( av/dx + 2.*nu/dy2 ) ) * v(i,j)     &
                        + dt * ( nu/dy2 + (av-v(i,j))/(2.*dy) ) * v(i,j+1) &
                        + dt * ( nu/dy2 + (av+v(i,j))/(2.*dy) ) * v(i,j-1)
            ENDDO
        ENDDO
#endif
        !==============================================================================================!
        ! 3.4) Solve the free surface equation
        !==============================================================================================!
        ! CONJUGATE GRADIENT METHOD
        !==============================================================================================!
        ! nb: even rhs is distributed among processors as eta (same dimensions), while Hu and Hv are not
        DO i = MPI%istart, MPI%iend 
            DO j = MPI%jstart, MPI%jend
                rhs(i,j) = eta(i,j) - dt/dx * ( Hu(i+1,j)*Fu(i+1,j) - Hu(i,j)*Fu(i,j)) &
                                    - dt/dy * ( Hv(i,j+1)*Fv(i,j+1) - Hv(i,j)*Fv(i,j))
            ENDDO
        ENDDO
        CALL CG(IMAX, JMAX, eta, rhs , MPI%istart, MPI%iend, MPI%jstart, MPI%jend, MPI%myrank)
        !==============================================================================================!
        ! 3.5) Update the velocity (momentum equation)
        !==============================================================================================!
        ct = g*dt/dx    ! Temporary coefficient for x
        cs = g*dt/dy    ! Temporary coefficient for y
        !==============================================================================================!
#ifdef PARALLEL
        u( MPI%istart, MPI%jstart:MPI%jend ) = Fu( MPI%istart, MPI%jstart:MPI%jend ) 
        u( MPI%iend+1, MPI%jstart:MPI%jend ) = Fu( MPI%iend+1, MPI%jstart:MPI%jend )
        DO i = MPI%istart+1, MPI%iend
            DO j = MPI%jstart, MPI%jend
                u(i,j) = Fu(i,j) - ct * ( eta(i,j) - eta(i-1,j) )
            ENDDO
        ENDDO
        !==============================================================================================!
        v( MPI%istart:MPI%iend, MPI%jstart ) = Fv( MPI%istart:MPI%iend, MPI%jstart )
        v( MPI%istart:MPI%iend, MPI%jend+1 ) = Fv( MPI%istart:MPI%iend, MPI%jend+1 )
        DO i = MPI%istart, MPI%iend
            DO j = MPI%jstart+1, MPI%jend
                v(i,j) = Fv(i,j) - cs * ( eta(i,j) - eta(i,j-1) )
            ENDDO
        ENDDO
#else        
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
#endif
        !==============================================================================================!
        ! 3.6) Update total water depth
        !==============================================================================================!
#ifdef PARALLEL
        DO j = MPI%jstart, MPI%jend
            ! Initialize first and last columns for Hu
            Hu( MPI%istart, j ) = MAX( 0.0, bu( MPI%istart, j ) + eta( MPI%istart, j ) )
            Hu( MPI%iend+1, j ) = MAX( 0.0, bu( MPI%iend+1, j ) + eta( MPI%iend,   j ) )
        ENDDO    
        DO i = MPI%istart, MPI%iend
            ! Initialize first and last rows for Hv
            Hv( i, MPI%jstart ) = MAX( 0.0, bv( i, MPI%jstart ) + eta( i, MPI%jstart ) )
            Hv( i, MPI%jend+1 ) = MAX( 0.0, bv( i, MPI%jend+1 ) + eta( i, MPI%jend   ) )
        ENDDO
        DO i = MPI%istart+1, MPI%iend
            DO j = MPI%jstart, MPI%jend
                Hu( i , j ) = MAXVAL( (/ 0.0, bu( i, j ) + eta( i - 1, j ), bu( i, j )+ eta( i, j ) /) )
            ENDDO
        ENDDO
        DO i = MPI%istart, MPI%iend
            DO j = MPI%jstart+1, MPI%jend
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
        !==========================================================================================!
        time = time + dt  ! Update time
        !==========================================================================================!
        ! 3.7) Eventually plot the results        
#ifdef PARALLEL
        IF(ABS(time-tio).LT.1e-12) THEN
            ! Collect all the data from each MPI subdomain
            ! Collect from eta into eta1
            DO i = MPI%istart, MPI%iend
                CALL MPI_ALLGATHER( eta( i, MPI%jstart:MPI%jend ), MPI%nElem, MPI%AUTO_REAL, eta1(i, MPI%jstart:MPI%jend), &
                                    MPI%nElem, MPI%AUTO_REAL, MPI_COMM_WORLD, MPI%iErr)
                CALL MPI_ALLGATHER( u( i, MPI%jstart:MPI%jend ), MPI%nElem, MPI%AUTO_REAL, u1(i, MPI%jstart:MPI%jend), &
                                    MPI%nElem, MPI%AUTO_REAL, MPI_COMM_WORLD, MPI%iErr)
                ! v has one more column
                CALL MPI_ALLGATHER( v( i, MPI%jstart:MPI%jend+1 ), MPI%nElem+1, MPI%AUTO_REAL, v1(i,MPI%jstart:MPI%jend+1), &
                                    MPI%nElem+1, MPI%AUTO_REAL, MPI_COMM_WORLD, MPI%iErr)
            ENDDO
            ! u has one more line
            CALL MPI_ALLGATHER( u( MPI%iend+1, MPI%jstart:MPI%jend ), MPI%nElem, MPI%AUTO_REAL, u1(MPI%iend+1, MPI%jstart:MPI%jend), &
                                    MPI%nElem, MPI%AUTO_REAL, MPI_COMM_WORLD, MPI%iErr)
            ! Finally, update eta, u and v
            eta = eta1
            u   = u1
            v   = v1
            IF (MPI%myrank.EQ.0) THEN
                ! Only root thread calls dataoutput
                WRITE(*,'(a,f15.7)') ' |   plotting data output at time ', time
                CALL DataOutput(n, 1, IMAX, 1, JMAX, MPI%myrank)
            ENDIF
            tio = tio + dtio
        ENDIF    
#else    
        IF(ABS(time-tio).LT.1e-12) THEN
            WRITE(*,'(a,f15.7)') ' |   plotting data output at time ', time
            CALL DataOutput(n, 1, IMAX, 1, JMAX, MPI%myrank)
            tio = tio + dtio        
        ENDIF    
#endif
    !==============================================================================================!
    ENDDO ! n cycle
    !==============================================================================================!
#ifdef PARALLEL
    IF(MPI%myrank.EQ.0) THEN
        WRITE(*,'(a)') ' | END of the COMPUTATION '
    ENDIF
    WCT2 = MPI_WTime() 
    CALL MPI_FINALIZE(MPI%iErr)
#else
    WRITE(*,'(a)') ' | END of the COMPUTATION '

    CALL CPU_TIME(WCT2) 
#endif  
    !==============================================================================================!
    ! Empty memory
    !==============================================================================================!
    CALL Deallocate_Var
    DEALLOCATE( eta, eta1, rhs )
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
SUBROUTINE matop2D(Ap, p, N, M, istart, iend, jstart, jend)
    !==================================================================================================!
    USE VarDef_mod
    IMPLICIT NONE
    !==================================================================================================!
    INTEGER             :: N, M 
    REAL                :: Ap(N,M), p(N,M)
    !==================================================================================================!
    INTEGER             :: i, j
    INTEGER, INTENT(IN) :: istart, iend, jstart, jend
    REAL                :: cx, cy
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
    DO j = jstart, jend
		DO i = istart, iend
			IF(i.EQ.istart) THEN
				Ap(i,j) = Ap(i,j) - cx * ( Hu(i+1,j)*(p(i+1,j)-p(i,j)) - 0.0 )
			ELSEIF(i.EQ.iend) THEN  
				Ap(i,j) = Ap(i,j) - cx * ( 0.0 - Hu(i,j)*(p(i,j)-p(i-1,j)) )
			ELSE
				Ap(i,j) = Ap(i,j) - cx * ( Hu(i+1,j)*(p(i+1,j)-p(i,j)) - Hu(i,j)*(p(i,j)-p(i-1,j)) )
			ENDIF 
		ENDDO
    ENDDO
    !==================================================================================================!	
    ! 3) Fluxes in y-direction
    !==================================================================================================!	
	DO j = jstart, jend
		DO i = istart, iend
			IF(j.EQ.jstart) THEN
				Ap(i,j) = Ap(i,j) - cy * ( Hv(i,j+1)*(p(i,j+1)-p(i,j)) - 0.0 )
			ELSEIF(j.EQ.jend) THEN  
				Ap(i,j) = Ap(i,j) - cy * ( 0.0 - Hv(i,j)*(p(i,j)-p(i,j-1)) )
			ELSE  
				Ap(i,j) = Ap(i,j) - cy * ( Hv(i,j+1)*(p(i,j+1)-p(i,j)) - Hv(i,j)*(p(i,j)-p(i,j-1)) )
			ENDIF 
		ENDDO
    ENDDO
!======================================================================================================!
END SUBROUTINE matop2D