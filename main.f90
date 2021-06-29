!------------------------------------------------------------------!
! Progetto II.
! Equazione dei flussi a superficie libera su griglie strutturate.
!------------------------------------------------------------------!
! Università di Ferrara - Dipartimento di Matematica e Informatica.
! Corso di Algoritmi per il Calcolo Parallelo.
! Prof. Walter Boscheri
! A.A. 2020/2021
!------------------------------------------------------------------!
! Studenti.
! Alberto Paparella - 144261
! Elena Zerbin      - 145302
!------------------------------------------------------------------!
    
PROGRAM FS2D
    !------------------------------------------------------------!
    USE VarDef_mod
    IMPLICIT NONE
    !------------------------------------------------------------!
    INTEGER            :: i, j, n
    REAL               :: umax, au, s0, ct
    REAL, ALLOCATABLE  :: avec(:), bvec(:), cvec(:)
    !REAL, ALLOCATABLE  :: rhs(:,:)  ! 2D, not sure, trying a different solution
    !------------------------------------------------------------!
    !
    WRITE(*,'(a)') ' | ========================================================================== | '
    WRITE(*,'(a)') ' |      Equazione dei flussi a superficie libera su griglie strutturate.      | '
    WRITE(*,'(a)') ' | ========================================================================== | '
    WRITE(*,'(a)') ' | '
    !
#ifdef CGM
    WRITE(*,'(a)') ' | Soluzione con il metodo del gradiente coniugato. '
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
    xL     = -0.5
    xR     = 0.5
    !
    JMAX   = 500    ! 2D: Max index for y coordinates
    yD     = -0.5   ! 2D
    yU     = 0.5    ! 2D
    s      = 0.1    ! 2D
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
    ALLOCATE( avec(IMAX)  )  
    ALLOCATE( bvec(IMAX)  )
    ALLOCATE( cvec(IMAX)  )
    ALLOCATE( rhs(IMAX, JMAX) )
    ALLOCATE( chs(IMAX, JMAX) )
    !ALLOCATE( rhs(IMAX, JMAX) ) ! 2D, not sure, trying a different solution
    !
    ! 1) Computational domain 
    !
    WRITE(*,'(a)') ' | Building computational domain... '
    !    
    dx    = (xR - xL) / REAL(IMAX)
    dx2   = dx**2
    x(1)  = xL
    !
    dy    = (yU - yD) / REAL(IMAX)  ! 2D
    dy2   = dy**2                   ! 2D
    y(1)  = yD                      ! 2D
    !    
    ! We decided to use two separates dimensions, IMAX and JMAX,
    ! for x and y coords for more flexibility
    DO i = 1, IMAX
      x(i+1) = x(i) + dx
      xb(i)  = x(i) + dx/2.
    ENDDO
    DO j = 1, JMAX             ! 2D
      y(j+1) = y(j) + dy       ! 2D
      yb(j)  = y(j) + dy/2.    ! 2D
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
        ! Note: it is correct to use xb and yb instead of x and y because x and y represents
        ! the vertex coords for u and v, while xb and yb represents the barycenter coords
        ! for eta.
        eta(i,j) = 1.0 + EXP( (-1.0 / (2 * (s**2) )) * ( xb(i)**2 + yb(j)**2 ))
      ENDDO
    ENDDO
    !---------------------------------------------
    !
    ! 2.1) Velocity, bottom and total water depth (interfaces)
    !
    !---------------------------------------------
    ! Note: u is of dimension (IMAX + 1) * JMAX, while v is of dimension IMAX * (JMAX + 1)
    ! Note: it is convenient to have two matrixes for b now, bu and bv
    ! TODO think about a solution using a single matrix, for eample using a matrix of dimensions IMAX*2+1 * JMAX*2+1
    !      and indices for both axis with step 2 (note the sparse form of the matrix, is it possible a better solution?)
    ! By the way, at the moment I think the 2 vectors solution is faster and also lighter in memory terms
    ! Same valuation values for H 
    DO i = 1, IMAX + 1
      DO j = 1, JMAX
        u(i, j)  = 0.0      ! fluid at rest
        bu(i, j) = 0.0      ! zero bottom elevation
      ENDDO
    ENDDO
    DO i = 1, IMAX
      DO j = 1, JMAX + 1
        v(i, j)  = 0.0      ! fluid at rest
        bv(i, j) = 0.0      ! zero bottom elevation
      ENDDO
    ENDDO
    !
    ! Total water depth
    !
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
    DO j = 2, JMAX
      DO i = 1, IMAX
        Hv( i , j ) = MAXVAL( (/ 0.0, bv( i, j ) + eta( i, j - 1 ), bv( i, j )+ eta( i, j ) /) )
      ENDDO
    ENDDO
    !
    CALL DataOutput(0)  ! plot initial condition
    tio = tio + dtio
    !---------------------------------------------
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
      dt2   = dt**2
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
      Fu(1,1)      = u(1,1)   !2D

      !Fu(IMAX+1) = u(IMAX+1)
      Fu(IMAX+1,JMAX+1) = u(IMAX+1,JMAX+1)  ! 2D
      !
      !is not needed because our formula does not include viscosity,
      !so we have to take the convention equation, with a being our velocity u,
      !and add and subtract the other velocity, depending on whether we take information
      !against redirection, propagation of information
      DO i = 2, IMAX
        DO j = 2, JMAX          ! 2D
            ! au = ABS(u(i))
            !au = ABS(u(i,j))    ! 2D
            ! Explicit upwind
            !Fu(i) = ( 1. - dt * ( au/dx + 2.*nu/dx2 ) ) * u(i)   &
            !      + dt * ( nu/dx2 + (au-u(i))/(2.*dx) ) * u(i+1) &
            !      + dt * ( nu/dx2 + (au+u(i))/(2.*dx) ) * u(i-1)
            Fu(i,j) = (eta(i,j)) - u(i+1,j)*(dt/dx * ( eta(i,j) - eta(i,j) ) )    ! 2D

        ENDDO   ! 2D
      ENDDO  
      ! 3.3) Compute the operator Fv
      !
      ! BC: no-slip wall
      !Fu(1)      = u(1)
      Fv(1, 1)           = v(1,1)  ! 2D
      
      !Fu(IMAX+1) = u(IMAX+1)
      Fv(IMAX+1,JMAX+1) = v(IMAX+1,JMAX+1)  ! 2D
        DO i = 2, IMAX
        DO j = 2, JMAX +1         ! 2D
            ! au = ABS(u(i))
            !au = ABS(u(i,j))    ! 2D
            ! Explicit upwind
            !Fu(i) = ( 1. - dt * ( au/dx + 2.*nu/dx2 ) ) * u(i)   &
            !      + dt * ( nu/dx2 + (au-u(i))/(2.*dx) ) * u(i+1) &
            !      + dt * ( nu/dx2 + (au+u(i))/(2.*dx) ) * u(i-1)
            Fv(i,j) = (eta(i,j)) - v(i,j+1)*(dt/dy * ( eta(i,j) - eta(i,j) )  )   ! 2D to try, if it doesn't work you have to find another formula with the upwind method

        ENDDO   ! 2D
        ENDDO
        !
        ! 3.4) Solve the system for the pressure
        !
        !compute the rhs
         DO j = 1, JMAX    ! 2D
          DO i = 1, IMAX
            !rhs(i) = eta(i) - dt/dx * ( H(i+1)*Fu(i+1) - H(i)*Fu(i) )
            rhs(i,j) = eta(i,j) - dt/dx * ( Hu(i+1,j+1)*Fu(i+1,j+1) - Hu(i,j)*Fu(i,j) )   ! 2D
          ENDDO
          !CALL CG(IMAX,eta,rhs)
          CALL CG(IMAX,eta,rhs) ! 2D
         ENDDO ! 2D
         
         !compute the chs
             DO j = 1, JMAX    ! 2D
          DO i = 1, IMAX
            !rhs(i) = eta(i) - dt/dx * ( H(i+1)*Fu(i+1) - H(i)*Fu(i) )
            chs(i,j) = eta(i,j) - dt/dy * ( Hv(i+1,j+1)*Fv(i+1,j+1) - Hv(i,j)*Fv(i,j) )
          ENDDO
          !CALL CG(IMAX,eta,rhs)
          CALL CG(JMAX,eta,chs) ! 2D
      ENDDO ! 2D
      !
      ! 3.5) Solve the free surface equation
      !
#ifdef CGM
      !
      ! CONJUGATE GRADIENT METHOD
      !
      DO j = 1, JMAX    ! 2D
          DO i = 1, IMAX
            !rhs(i) = eta(i) - dt/dx * ( H(i+1)*Fu(i+1) - H(i)*Fu(i) )
            rhs(i,j) = eta(i,j) - dt/dx * ( Hu(i+1,j+1)*Fu(i+1,j+1) - Hu(i,j)*Fu(i,j) )   ! 2D, TODO use real formula, this is wrong and just for testing
            chs(i,j) = eta(i,j) - dt/dy * ( Hv(i+1,j+1)*Fv(i+1,j+1) - Hv(i,j)*Fv(i,j) )
          ENDDO
          !CALL CG(IMAX,eta,rhs)
          CALL CG(IMAX,eta(:,j),rhs) ! 2D
      ENDDO ! 2D
      !
      
#else
      !
      ! THOMAS ALGORITHM
      !
      !ct = g*dt**2/dx2  ! temporary coefficient
      !
      !DO i = 1, IMAX
      !  IF(i.EQ.1) THEN
      !    avec(i)  = 0.
      !   bvec(i)  = 1. + ct * ( H(i+1) + 0.0 )
      !    cvec(i)  = - ct * H(i+1)
      !    rhs(i) = eta(i) - dt/dx * ( H(i+1)*Fu(i+1) - H(i)*Fu(i) )
      !  ELSEIF(i.EQ.IMAX) THEN  
      !    avec(i)  = - ct * H(i)
      !    bvec(i)  = 1. + ct * ( 0.0 + H(i) )
      !    cvec(i)  = 0.
      !    rhs(i) = eta(i) - dt/dx * ( H(i+1)*Fu(i+1) - H(i)*Fu(i) )
      !  ELSE  
      !    avec(i)  = - ct * H(i)
      !    bvec(i)  = 1. + ct * ( H(i+1) + H(i) )
      !    cvec(i)  = - ct * H(i+1)
      !    rhs(i) = eta(i) - dt/dx * ( H(i+1)*Fu(i+1) - H(i)*Fu(i) )
      !  ENDIF
      !ENDDO  
      !CALL Thomas(eta,avec,bvec,cvec,rhs,IMAX)
      !
#endif      
      !
      ! 3.4) Update the velocity (momentum equation)
      !
      ct = g*dt/dx  ! temporary coefficient
      cs = g*dt/dy ! temporary coefficient
      !
      !u(1)      = Fu(1)
      u(1,1)      = Fu(1,1)                 ! 2D
      !u(IMAX+1) = Fu(IMAX+1)
      u(IMAX+1,JMAX+1) = Fu(IMAX+1,JMAX+1)  ! 2D
      !DO i = 2, IMAX
      !  u(i) = Fu(i) - ct * ( eta(i) - eta(i-1) )  
      !ENDDO
      DO i = 2, IMAX                                            ! 2D
          DO j = 2, JMAX                                        ! 2D
            u(i,j) = Fu(i,j) - ct * ( eta(i,j) - eta(i-1,j-1) ) ! 2D, TODO use real formula, this is wrong and just for testing
          ENDDO                                                 ! 2D
      ENDDO                                                     ! 2D
      !
      time = time + dt  ! update time
      !
      
        !
      !u(1)      = Fu(1)
      v(1,1)      = Fv(1,1)                 ! 2D
      !u(IMAX+1) = Fu(IMAX+1)
      v(IMAX+1,JMAX+1) = Fv(IMAX+1,JMAX+1)  ! 2D
      !DO i = 2, IMAX
      !  u(i) = Fu(i) - ct * ( eta(i) - eta(i-1) )  
      !ENDDO
      DO i = 2, IMAX                                            ! 2D
          DO j = 2, JMAX                                        ! 2D
            v(i,j) = Fu(i,j) - cs * ( eta(i,j) - eta(i-1,j-1) ) ! 2D, TODO use real formula, this is wrong and just for testing
          ENDDO                                                 ! 2D
      ENDDO 
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
    DEALLOCATE( avec, bvec, cvec, rhs )
    !
    WRITE(*,'(a)') ' | '
    WRITE(*,'(a)') ' |         Finalization was successful. Bye :-)           | '
    WRITE(*,'(a)') ' | ====================================================== | '
    !
END PROGRAM FS2D 
  
SUBROUTINE matop1D(Ap,p,N) 
    !------------------------------------------------------------!
    USE VarDef_mod
    IMPLICIT NONE
    !------------------------------------------------------------!
    INTEGER  :: N
    REAL     :: Ap(N), p(N)
    !
    INTEGER  :: i, j
    REAL     :: ct, avec, bvec, cvec
    !------------------------------------------------------------!
    !
    ct = g*dt**2/dx2  ! temporary coefficient
    !
    DO i = 1, IMAX
      DO j = 1, JMAX    ! 2D
          if(i.eq.250) then
            continue
          endif  
          IF(i.EQ.1) THEN
            !bvec    = 1. + ct * ( H(i+1) + 0.0 )
            bvec    = 1. + ct * ( Hu(i+1, j+1) + 0.0 )     ! 2D, TODO use real formula, this is wrong and just for testing
            !cvec    = - ct * H(i+1)
            cvec    = - ct * Hu(i+1, j+1)                  ! 2D, TODO use real formula, this is wrong and just for testing
            Ap(i) = bvec*p(i) + cvec*p(i+1)
          ELSEIF(i.EQ.IMAX) THEN  
            !avec    = - ct * H(i)
            avec    = - ct * Hu(i,j)                       ! 2D, TODO use real formula, this is wrong and just for testing
            !bvec    = 1. + ct * ( 0.0 + H(i) )
            bvec    = 1. + ct * ( 0.0 + Hu(i,j) )          ! 2D, TODO use real formula, this is wrong and just for testing
            Ap(i) = avec*p(i-1) + bvec*p(i) 
          ELSE  
            !avec    = - ct * H(i)
            avec    = - ct * Hu(i,j)                       ! 2D, TODO use real formula, this is wrong and just for testing
            !bvec    = 1. + ct * ( H(i+1) + H(i) )
            bvec    = 1. + ct * ( Hu(i+1,j+1) + Hu(i,j) )   ! 2D, TODO use real formula, this is wrong and just for testing
            !cvec    = - ct * H(i+1)
            cvec    = - ct * Hu(i+1,j+1)                   ! 2D, TODO use real formula, this is wrong and just for testing
            Ap(i) = avec*p(i-1) + bvec*p(i) + cvec*p(i+1)
          ENDIF  
          !
      ENDDO       ! 2D
    ENDDO
    !
END SUBROUTINE matop1D