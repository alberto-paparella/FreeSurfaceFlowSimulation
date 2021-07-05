SUBROUTINE Thomas(x,a,b,c,d,N) 
   IMPLICIT NONE 
   ! Argument list 
   INTEGER :: N 
   REAL    :: x(N), a(N), b(N), c(N), d(N) 
   ! Local variables 
   REAL    :: gamma
   INTEGER :: i  
   !
   c(1) = c(1)/b(1)
   d(1) = d(1)/b(1) 
   ! Part I: forward elimination
   DO i = 2, N
      gamma=1.0/(b(i)-c(i-1)*a(i)) 
      c(i) = c(i)*gamma 
      d(i) = (d(i)-a(i)*d(i-1))*gamma 
   ENDDO
   ! Part II: back substitution
   x(N)=d(N); 
   DO i = N-1, 1, -1 
    x(i)=d(i)-c(i)*x(i+1) 
   ENDDO 
END SUBROUTINE Thomas 
  
SUBROUTINE CG(N,x,b)
    IMPLICIT NONE
    !----------------------------------------!
    INTEGER         :: N                        ! size of the linear system
    REAL            :: x(N,N)                  ! solution  
    REAL            :: b(N,N)                  ! right hand side  
    !----------------------------------------!
    INTEGER         :: k, KMAX, iErr
    REAL            :: Ax(N,N), Ap(N,N)
    REAL            :: r(N,N), p(N,N)
    REAL            :: pAp, lambda
    REAL            :: alphak, alpha
    REAL, PARAMETER :: tol = 1e-12           ! tolerance for convergence  
    !----------------------------------------!
    !
    x = b                 ! initial guess   
    CALL matop2D(Ax,x,N)  ! matrix-vector multiplication
    r = b - Ax            ! residual   
    p = r                 ! search direction = max. descent   
    alphak = SUM(r*r) 
    !
    KMAX = N
    !
    DO k = 1, KMAX
      !
      IF(SQRT(alphak).LT.tol) THEN
        WRITE(*,'(a,i3,a,e15.7)') ' |   CG iter: ', k, ' CG res: ', SQRT(alphak)
        RETURN
      ENDIF    
      !   
      CALL matop2D(Ap,p,N)    
      pAp    = SUM(p*Ap)        
      lambda = alphak / pAp
      x      = x + lambda*p
      r      = r - lambda*Ap
      alpha  = SUM(r*r)
      p      = r + alpha/alphak * p
      alphak = alpha
      !
    ENDDO !k    
    !
    IF(k.GE.KMAX) THEN
      PRINT *, ' ERROR. Conjugate gradient did not converge! ', SQRT(alphak)
      STOP     
    ENDIF
    ! 
    END SUBROUTINE CG  
    
  