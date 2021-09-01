!======================================================================================================!
! Progetto II
! Equazione dei flussi a superficie libera su griglie strutturate
!======================================================================================================!
! Universit� degli studi di Ferrara - Dipartimento di Matematica e Informatica
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
SUBROUTINE CG(N, M, x, b, istart, iend, jstart, jend, myrank)
    IMPLICIT NONE
    !==================================================================================================!
    INTEGER             :: i
    INTEGER             :: N, M     ! Size of the linear system
    REAL                :: x(N,M)   ! Solution (in matricial form)
    REAL                :: b(N,M)   ! Right hand side (in matricial form)
    !==================================================================================================!
    INTEGER             :: k, KMAX
    INTEGER, INTENT(IN) :: istart, iend, jstart, jend, myrank
    REAL                :: Ax(N,M), Ax1(N,M), Ap(N,M)
    REAL                :: r(N,M), p(N,M)
    REAL                :: pAp, lambda
    REAL                :: alphak, alpha
    REAL, PARAMETER     :: tol = 1e-12      ! Tolerance for convergence  
    !==================================================================================================!
    x = b                 ! Initial guess
    CALL matop2D(Ax, x, N, M, istart, iend, jstart, jend)  ! It is implemented into the main file
    r = b - Ax            ! Residual   
    p = r                 ! Search direction = max. descent   
    alphak = SUM(r*r)
    !==================================================================================================!
    KMAX = N*(jend-jstart)
    !==================================================================================================!
    DO k = 1, KMAX
        !==============================================================================================!
        IF(SQRT(alphak).LT.tol) THEN
            WRITE(*,'(a,i3,a,i3,a,e15.7)') ' |   myrank: ', myrank, ' CG iter: ', k, ' CG res: ', SQRT(alphak)
            RETURN
        ENDIF
        !==============================================================================================!
        CALL matop2D(Ap, p, N, M, istart, iend, jstart, jend)
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
        STOP     
    ENDIF
!======================================================================================================!
END SUBROUTINE CG    