PROGRAM TEST
    IMPLICIT NONE
    PRINT *, ' Hello World '
    END PROGRAM TEST
    
    
    
# variabili
#tirante idrico # massimo fra 3 valori (funzione?)
#dominio di calcolo omega
# condizione iniziale



# calore una dimensione openmpi

#define PARALLEL

# Soluzione equazione dei flussi a superficie libera su griglie strutturate
PROGRAM FLOW
    !------------------------------------------------------------!
    IMPLICIT NONE
#ifdef PARALLEL     
    INCLUDE 'mpif.h'
#endif     
    !------------------------------------------------------------!
    INTEGER             :: i, n        ! space and time index
    !
    INTEGER, PARAMETER  :: IMAX=4000   ! total number of cells 
    REAL                :: xL, xR      ! left and right domain definition
    REAL                :: xD          ! location of discontinuity
    REAL                :: dx, dx2     ! mesh size and square of mesh size
    REAL, ALLOCATABLE   :: x(:)        ! vertex coords
    !
    REAL                :: CFL         ! Courant-Friedrichs-Lewy number for stability condition (CFL<1) 
    REAL                :: time        ! current time
    REAL                :: dt          ! time step
    REAL                :: tend        ! final time
    INTEGER, PARAMETER  :: NMAX = 1e6  ! maximum number of time steps 
    REAL, ALLOCATABLE   :: T(:), T1(:) ! temperature 
    !
    CHARACTER(LEN=200)  :: TestName    ! name of test problem
    REAL                :: tio         ! output time and time step 
    REAL                :: dtio        ! output time step
    !
    REAL                :: d           ! stability parameter  
    REAL                :: kappa       ! heat conduction coefficient
    REAL                :: TL, TR      ! boundary conditions
#ifdef PARALLEL
    INTEGER             :: LCPU, RCPU, MsgLength, nMsg 
    REAL                :: send_messageL, send_messageR, recv_messageL, recv_messageR
    INTEGER             :: send_request(2), recv_request(2) 
    INTEGER             :: send_status_list(MPI_STATUS_SIZE,2),recv_status_list(MPI_STATUS_SIZE,2)   
#endif     
    !     
    TYPE tMPI
       INTEGER :: myrank, nCPU, iErr 
       INTEGER :: AUTO_REAL
       INTEGER :: nElem,istart,iend 
    END TYPE tMPI 
    TYPE(tMPI) :: MPI
    REAL       :: realtest
    REAL       :: WCT1, WCT2           ! Wall clock times 
    !
    !------------------------------------------------------------!
    !
    WCT1 = 0. 
    WCT2 = 0. 
    !
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
      WRITE(*,'(a)') '      Finite difference code for 1D heat equation       '
      WRITE(*,'(a,i10)') '         Total number of CPUs used : ', MPI%nCPU       
      WRITE(*,'(a)') ' ====================================================== '
      WRITE(*,'(a)') ' '
    ENDIF  
    !
    !======================================== 
    ! SETTINGS OF THE COMPUTATION
    !========================================
    !
    TestName = 'Heat1D_explicit'
    !
    xL     = -1.0
    xR     = 1.0
    xD     = 0.0  ! location of the discontinuity
    !
    time   = 0.0
    tend   = 0.02    
    dtio   = 2e-2
    tio    = dtio
    CFL    = 0.9   ! for stability condition CFL<1
    d      = 0.45  ! for stability condition d<0.5   
    !
    kappa  = 1     ! heat conduction coefficient 
    !
    ! Boundary conditions
    TL = 100
    TR = 50
    !
    !========================================
    !
    ! Allocate variables
    ALLOCATE( x(IMAX)  )
    ALLOCATE( T1(IMAX) )
    ALLOCATE( T(MPI%istart-1:MPI%iend+1) )
    !
    !---------------------------------------- 
    !---------------------------------------- 
    !
    ! 1) Computational domain 
    !
    IF(MPI%myrank.EQ.0) THEN
      WRITE(*,'(a)') '  Building computational domain... '
    ENDIF  
    !    
    dx    = (xR - xL) / REAL(IMAX)
    dx2   = dx**2
    x(1)  = xL
    !
    DO i = 1, IMAX
      x(i) = xL + 0.5*dx + (i-1)*dx
    ENDDO  
    !
    !---------------------------------------- 
    !---------------------------------------- 
    !
    ! 2) Initial condition
    !
    IF(MPI%myrank.EQ.0) THEN
      WRITE(*,'(a)') '  Assigning initial condition... '
    ENDIF  
    !
    DO i = MPI%istart, MPI%iend
      IF(x(i).LT.xD) THEN
        T(i) = TL
      ELSE
        T(i) = TR
      ENDIF  
    ENDDO  
    !    
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
END PROGRAM Heat1D_MPI

SUBROUTINE DataOutput(timestep,TestName,myrank,istart,iend,IMAX,x,T)
    !------------------------------------------------------------!
    IMPLICIT NONE
    !------------------------------------------------------------!
    INTEGER,            INTENT(IN) :: timestep, IMAX, istart, iend, myrank
    CHARACTER(LEN=200), INTENT(IN) :: TestName
    REAL,               INTENT(IN) :: x(IMAX), T(istart:iend)
    !
    INTEGER                        :: i, DataUnit
    CHARACTER(LEN=10)              :: citer, cmyrank
    CHARACTER(LEN=200)             :: IOFileName
    !------------------------------------------------------------!
    !
    WRITE(citer,'(i10.10)') timestep                      ! convert iteration number to string
    WRITE(cmyrank,'(I4.4)') myrank                        ! convert iteration number to string
    IOFileName = TRIM(TestName)//'-'//TRIM(citer)//'-'//TRIM(cmyrank)//'.dat' ! name of output file
    DataUnit   = 100+myrank                               ! unit for output file
    !
    OPEN(UNIT=DataUnit, FILE=TRIM(IOFilename), STATUS='UNKNOWN', ACTION='WRITE')
    !
    ! Header
    WRITE(DataUnit,*) IMAX
    ! Coordinates
    DO i = istart, iend
      WRITE(DataUnit,*) x(i)
    ENDDO  
    ! Temperature
    DO i = istart, iend
      WRITE(DataUnit,*) T(i)
    ENDDO
    !
    CLOSE(DataUnit)
    !
    END SUBROUTINE DataOutput  
  
# Flussi a superficie libera, una dimensione
    
PROGRAM FS1D
    !------------------------------------------------------------!
    USE VarDef_mod
    IMPLICIT NONE
    !------------------------------------------------------------!
    INTEGER            :: i, n
    REAL               :: umax, au, s0, ct
    REAL, ALLOCATABLE  :: av(:), bv(:), cv(:), rhs(:)
    !------------------------------------------------------------!
    !
    WRITE(*,'(a)') ' | ====================================================== | '
    WRITE(*,'(a)') ' |      Semi-implicit code for 1D free surface flows      | '
    WRITE(*,'(a)') ' | ====================================================== | '
    WRITE(*,'(a)') ' | '
    !
#ifdef CGM
    WRITE(*,'(a)') ' | Solver: conjugate gradient method. '
#else
    WRITE(*,'(a)') ' | Solver: Thomas algorithm. '
#endif
    WRITE(*,'(a)') ' | '
    !
    !======================================== 
    ! SETTINGS OF THE COMPUTATION
    !========================================
    !
    TestName = 'Gaussian_test'
    !
    IMAX   = 500    ! Max index for x coordinates
    JMAX   = 500    ! Max index for y coordinates
    xL     = -0.5
    xR     = 0.5
    yB     = -0.5   !
    yU     = 0.5    !
    s      = 0.1    !
    !
    time   = 0.0
    tend   = 0.1
    tio    = 0.0
    dtio   = 5e-2
    CFL    = 0.9
    dt_fix = 1e-2
    !
    nu   = 0.e-2
    !
    !========================================
    !
    CALL Allocate_Var  ! allocate all variables
    !
    ! vectors for assembling the tridiagonal linear system
    ALLOCATE( av(IMAX)  )  
    ALLOCATE( bv(IMAX)  )
    ALLOCATE( cv(IMAX)  )
    ALLOCATE( rhs(IMAX) )
    !
    ! 1) Computational domain 
    !
    WRITE(*,'(a)') ' | Building computational domain... '
    !    
    dx    = (xR - xL) / REAL(IMAX)
    dx2   = dx**2
    x(1)  = xL
    
    !
    dy    = (yU - yB) / REAL(IMAX)  !
    dy    = dy**2                   !
    y(1)  = yB                      !
    !
    
    ! We decided to use two separates dimensions, IMAX and JMAX, for x and y coords
    ! for more flexibility
    DO i = 1, IMAX
      x(i+1) = x(i) + dx
      xb(i)  = x(i) + dx/2.
    ENDDO
    DO j = 1, JMAX
      y(j+1) = y(j) + dy
      yb(j)  = y(j) + dy/2.
    ENDDO
    !
    !---------------------------------------- 
    !---------------------------------------- 
    !
    ! 2) Initial condition
    !
    WRITE(*,'(a)') ' | Assigning initial condition... '
    !
    ! 2.1) Free surface elevation (barycenters)
    !
    DO i = 1, IMAX
      DO j = 1, JMAX
        ! Gaussian profile
        eta(i,j) = 1.0 + EXP( (-1.0 / (2 * (s**2) )) * ( x(i)**2 + y(j)**2 )) ! x or xb ?
      ENDDO
    ENDDO  
    !
    ! 2.1) Velocity, bottom and total water depth (interfaces)
    !
    DO i = 1, IMAX+1
      DO j = 1, JMAX+1
        u(i, j) = 0.0   ! fluid at rest
        h(i, j) = 0.0   ! zero bottom elevation
        v(i, j) = 0.0   !
      ENDDO
    ENDDO
    
    
    
    
    !
    ! Total water depth
    !
    H(1,1)    = MAX( 0.0, h(1,1)+eta(1,1) )
    H(IMAX+1,JMAX+1) = MAX( 0.0, h(IMAX+1,JMAX+1)+eta(IMAX,JMAX) )
    DO i = 2, IMAX
      DO j = 2, JMAX
        H(i,j) = MAXVAL( (/ 0.0, h(i,j)+eta(i-1,j-1), h(i,j)+eta(i,j) /) )
      ENDDO
    ENDDO  
    !
    CALL DataOutput(0)  ! plot initial condition
    tio = tio + dtio
    !
    !---------------------------------------- 
    !---------------------------------------- 
    !
    ! 3) Computation: main loop in time
    !
    WRITE(*,'(a)') ' | '
    WRITE(*,'(a)') ' | START of the COMPUTATION '
    !
    DO n = 1, NMAX
      !
      ! 3.1) Compute the time step
      !
      IF(time.GE.tend) THEN
        EXIT    
      ENDIF
      !
      !umax = MAXVAL( ABS(u) + SQRT(g*H))  ! fully explicit time step
      umax = MAXVAL( ABS(u) )             ! semi-implicit time step
      dt   = MIN( dt_fix, CFL / ( (umax+1e-14)/dx + 2.*nu/dx2 ) )
      IF((time+dt).GT.tend) THEN
        dt  = tend - time   
        tio = tend
      ENDIF    
      IF((time+dt).GT.tio) THEN
        dt = tio - time
      ENDIF
      !
      ! 3.2) Compute the operator Fu
      !
      ! BC: no-slip wall
      Fu(1)      = 0.0  
      Fu(IMAX+1) = 0.0 
      !
      DO i = 2, IMAX
        au = ABS(u(i))
        ! Explicit upwind
        Fu(i) = ( 1. - dt * ( au/dx + 2.*nu/dx2 ) ) * u(i)   &
              + dt * ( nu/dx2 + (au-u(i))/(2.*dx) ) * u(i+1) &
              + dt * ( nu/dx2 + (au+u(i))/(2.*dx) ) * u(i-1) 
      ENDDO  
      !
      ! 3.3) Solve the free surface equation
      !
#ifdef CGM
      !
      ! CONJUGATE GRADIENT METHOD
      !     
      DO i = 1, IMAX
        rhs(i) = eta(i) - dt/dx * ( H(i+1)*Fu(i+1) - H(i)*Fu(i) )
      ENDDO 
      CALL CG(IMAX,eta,rhs)
      !
#else
      !
      ! THOMAS ALGORITHM
      !
      ct = g*dt**2/dx2  ! temporary coefficient
      !
      DO i = 1, IMAX
        IF(i.EQ.1) THEN
          av(i)  = 0.
          bv(i)  = 1. + ct * ( H(i+1) + 0.0 )
          cv(i)  = - ct * H(i+1)
          rhs(i) = eta(i) - dt/dx * ( H(i+1)*Fu(i+1) - H(i)*Fu(i) )
        ELSEIF(i.EQ.IMAX) THEN  
          av(i)  = - ct * H(i)
          bv(i)  = 1. + ct * ( 0.0 + H(i) )
          cv(i)  = 0.
          rhs(i) = eta(i) - dt/dx * ( H(i+1)*Fu(i+1) - H(i)*Fu(i) )
        ELSE  
          av(i)  = - ct * H(i)
          bv(i)  = 1. + ct * ( H(i+1) + H(i) )
          cv(i)  = - ct * H(i+1)
          rhs(i) = eta(i) - dt/dx * ( H(i+1)*Fu(i+1) - H(i)*Fu(i) )
        ENDIF
      ENDDO  
      CALL Thomas(eta,av,bv,cv,rhs,IMAX)
      !
#endif      
      !
      ! 3.4) Update the velocity (momentum equation)
      !
      ct = g*dt/dx  ! temporary coefficient
      !
      u(1)      = Fu(1)
      u(IMAX+1) = Fu(IMAX+1)
      DO i = 2, IMAX
        u(i) = Fu(i) - ct * ( eta(i) - eta(i-1) )  
      ENDDO  
      !
      time = time + dt  ! update time
      !
      ! 3.5) Eventually plot the results
      IF(ABS(time-tio).LT.1e-12) THEN
        WRITE(*,'(a,f15.7)') ' |   plotting data output at time ', time
        CALL DataOutput(n)
        tio = tio + dtio
      ENDIF  
      !
    ENDDO !n
    !
    WRITE(*,'(a)') ' | END of the COMPUTATION '
    !
    !---------------------------------------- 
    !---------------------------------------- 
    !
    ! Empty memory
    !
    CALL Deallocate_Var
    !
    DEALLOCATE( av, bv, cv, rhs )
    !
    WRITE(*,'(a)') ' | '
    WRITE(*,'(a)') ' |         Finalization was successful. Bye :-)           | '
    WRITE(*,'(a)') ' | ====================================================== | '
    !
END PROGRAM FS1D 
  
SUBROUTINE matop1D(Ap,p,N) 
    !------------------------------------------------------------!
    USE VarDef_mod
    IMPLICIT NONE
    !------------------------------------------------------------!
    INTEGER  :: N
    REAL     :: Ap(N), p(N)
    !
    INTEGER  :: i
    REAL     :: ct, av, bv, cv
    !------------------------------------------------------------!
    !
    ct = g*dt**2/dx2  ! temporary coefficient
    !
    DO i = 1, IMAX
      if(i.eq.250) then
        continue
      endif  
      IF(i.EQ.1) THEN
        bv    = 1. + ct * ( H(i+1) + 0.0 )
        cv    = - ct * H(i+1)
        Ap(i) = bv*p(i) + cv*p(i+1)
      ELSEIF(i.EQ.IMAX) THEN  
        av    = - ct * H(i)
        bv    = 1. + ct * ( 0.0 + H(i) )
        Ap(i) = av*p(i-1) + bv*p(i) 
      ELSE  
        av    = - ct * H(i)
        bv    = 1. + ct * ( H(i+1) + H(i) )
        cv    = - ct * H(i+1)
        Ap(i) = av*p(i-1) + bv*p(i) + cv*p(i+1)
      ENDIF  
      !
    ENDDO  
    !
END SUBROUTINE matop1D