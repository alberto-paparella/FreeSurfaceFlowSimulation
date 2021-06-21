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
    INTEGER            :: i, n
    REAL               :: umax, au, s0, ct
    REAL, ALLOCATABLE  :: av(:), bv(:), cv(:), rhs(:)
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
#ifdef FS2D
    !JMAX   = 500    ! Max index for y coordinates
    !yB     = -0.5   !
    !yU     = 0.5    !
    !s      = 0.1    !
#endif
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
#ifdef FS2D
    !dy    = (yU - yB) / REAL(IMAX)  !
    !dy    = dy**2                   !
    !y(1)  = yB                      !
#endif
    !
    
    ! We decided to use two separates dimensions, IMAX and JMAX, for x and y coords
    ! for more flexibility
    DO i = 1, IMAX
      x(i+1) = x(i) + dx
      xb(i)  = x(i) + dx/2.
    ENDDO
#ifdef FS2D
    !DO j = 1, JMAX
    !  y(j+1) = y(j) + dy
    !  yb(j)  = y(j) + dy/2.
    !ENDDO
#endif
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
    
#ifdef FS2D
    !DO i = 1, IMAX
    !  DO j = 1, JMAX
    !    ! Gaussian profile
    !    eta(i,j) = 1.0 + EXP( (-1.0 / (2 * (s**2) )) * ( x(i)**2 + y(j)**2 )) ! x or xb ?
    !  ENDDO
    !ENDDO
#else
    DO i = 1, IMAX
      ! Gaussian profile
      eta(i) = 1.0 + EXP(-500.0*((xb(i))**2)/4.)
    ENDDO 
#endif
    !
    ! 2.1) Velocity, bottom and total water depth (interfaces)
    !
#ifdef FS2D
    !DO i = 1, IMAX+1
    !  DO j = 1, JMAX+1
    !    u(i, j) = 0.0   ! fluid at rest
    !    b(i, j) = 0.0   ! zero bottom elevation
    !    v(i, j) = 0.0   !
    !  ENDDO
    !ENDDO
    !
    ! Total water depth
    !
    !H(1,1)    = MAX( 0.0, b(1,1)+eta(1,1) )
    !H(IMAX+1,JMAX+1) = MAX( 0.0, b(IMAX+1,JMAX+1)+eta(IMAX,JMAX) )
    !DO i = 2, IMAX
    !  DO j = 2, JMAX
    !    H(i,j) = MAXVAL( (/ 0.0, b(i,j)+eta(i-1,j-1), b(i,j)+eta(i,j) /) )
    !  ENDDO
    !ENDDO  
    !
    !CALL DataOutput(0)  ! plot initial condition
    !tio = tio + dtio
    !
#else
    DO i = 1, IMAX+1
      u(i) = 0.0   ! fluid at rest
      b(i) = 0.0   ! zero bottom elevation
    ENDDO
    !
    ! Total water depth
    !
    H(1)      = MAX( 0.0, b(1)+eta(1)         )
    H(IMAX+1) = MAX( 0.0, b(IMAX+1)+eta(IMAX) )
    DO i = 2, IMAX
      H(i) = MAXVAL( (/ 0.0, b(i)+eta(i-1), b(i)+eta(i) /) )
    ENDDO  
    !
    CALL DataOutput(0)  ! plot initial condition
    tio = tio + dtio
#endif
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
END PROGRAM FS2D 
  
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