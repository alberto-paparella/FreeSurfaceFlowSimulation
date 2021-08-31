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
! LinSystem.f90
! In this file we implement the routine for the coniugate gradient method
!======================================================================================================!
SUBROUTINE CG(N,x,b,istart, iend, jstart, jend)
    IMPLICIT NONE
    !==================================================================================================!
    INTEGER             :: N        ! Size of the linear system
    REAL                :: x(N,N)   ! Solution (in matricial form)
    REAL                :: b(N,N)   ! Right hand side (in matricial form)
    !==================================================================================================!
    INTEGER             :: k, KMAX, iErr
    INTEGER, INTENT(IN) :: istart, iend, jstart, jend
    REAL                :: Ax(N,N), Ap(N,N)
    REAL                :: r(N,N), p(N,N)
    REAL                :: pAp, lambda
    REAL                :: alphak, alpha
    REAL, PARAMETER     :: tol = 1e-12      ! Tolerance for convergence  
    !==================================================================================================!
    PRINT *, 'Errore 1'
    x = b                 ! Initial guess   
    PRINT *, 'Errore 2'
!#ifdef PARALLEL    
    CALL matop2D(Ax,x,N,istart,iend,jstart,jend)  ! Matrix-matrix multiplication (it is implemented into the main file)
!#else
    !CALL matop2D(Ax,x,N)  
!#endif
    r = b - Ax            ! Residual   
    p = r                 ! Search direction = max. descent   
    alphak = SUM(r*r) 
    !==================================================================================================!
    KMAX = N*N
    !==================================================================================================!
    DO k = 1, KMAX
        !==============================================================================================!
        IF(SQRT(alphak).LT.tol) THEN
            WRITE(*,'(a,i3,a,e15.7)') ' |   CG iter: ', k, ' CG res: ', SQRT(alphak)
            RETURN
        ENDIF
        !==============================================================================================!
!#ifdef PARALLEL    
        CALL matop2D(Ap,p,N,istart,iend,jstart, jend) 
!#else        
        !CALL matop2D(Ap,p,N) 
!#endif   
        pAp    = SUM(p*Ap)        
        lambda = alphak / pAp
        x      = x + lambda*p
        r      = r - lambda*Ap
        alpha  = SUM(r*r)
        p      = r + alpha/alphak * p
        alphak = alpha
        !==============================================================================================!
    ENDDO   ! k cycle
       IF(k.GE.KMAX) THEN
        PRINT *, ' ERROR. Conjugate gradient did not converge! ', SQRT(alphak)
        PRINT *, k
        !STOP     
    ENDIF
    !==================================================================================================!


!======================================================================================================!
    END SUBROUTINE CG
    