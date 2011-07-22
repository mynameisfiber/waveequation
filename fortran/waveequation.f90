MODULE waveequation
  IMPLICIT NONE

  CONTAINS

  SUBROUTINE step(U, n, m, ng, dx, dt)
    !f2py REAL, DIMENSION(n,m,2), INTENT(in,out) :: U
    !f2py REAL, INTENT(in)                       :: dt, dx
    !f2py INTEGER, INTENT(in)                    :: ng
    !f2py INTEGER, INTENT(hide)                  :: n, m
    REAL, DIMENSION(n,m,2) :: U, Unew
    REAL                   :: dt, dx
    INTEGER                :: ng
    INTEGER                :: n, m

    Unew =         U                +         dt * dU(U   , n, m, ng, dx)
    Unew = 0.75  * U + 0.25  * Unew + 0.25  * dt * dU(Unew, n, m, ng, dx)
    U    = 1./3. * U + 2./3. * Unew + 2./3. * dt * dU(Unew, n, m, ng, dx)
  END SUBROUTINE
  
  FUNCTION dU(U, n, m, ng, dx) 
    !f2py REAL, DIMENSION(n,m,2), INTENT(in) :: U
    !f2py INTEGER, INTENT(hide)              :: n,m
    !f2py INTEGER, INTENT(in)                :: ng
    !f2py REAL, INTENT(in)                   :: dx
    REAL, DIMENSION(n,m,2) :: U, dU
    INTEGER                :: n,m
    INTEGER                :: i,j, ng
    REAL                   :: dx

    !do j=ng+1, m-ng
    !  do i=ng+1, n-ng
    !    dU(i,j,1) = U(i,j,2)
    !    dU(i,j,2) = (U(i-1,j,1) + U(i+1,j,1) - 4*U(i,j,1) + &
    !                 U(i,j-1,1) + U(i,j+1,1)) / (dx*dx)
    !  end do
    !end do
    dU(:,:,1) = U(:,:,2)
    dU(:,:,2) = dU2dx2full(U(:,:,1), n, m, dx, 1) + &
                dU2dx2full(U(:,:,1), n, m, dx, 2) 

    return
  END FUNCTION

  FUNCTION dU2dx2full(U, n, m, dx, dir) RESULT(dU2dx2)
    REAL, DIMENSION(n,m) :: U, dU2dx2
    INTEGER              :: n,m,dir
    REAL                 :: dx

    !dU2dx2 = (CSHIFT(U, 1, dir) + CSHIFT(U, -1, dir) - 2*U) / dx
    
    ! Assume we have 2 ghost zones.
    dU2dx2 = (-1./12. * (CSHIFT(U,2,dir) + CSHIFT(U,-2,dir)) + & 
               4./ 3. * (CSHIFT(U,1,dir) + CSHIFT(U,-1,dir)) - &
                 2.5  * U                                      &
             ) / (dx*dx)

    return
  END FUNCTION

  SUBROUTINE fill(U, n, m)
    !f2py REAL, DIMENSION(n,m,2), INTENT(in,out) :: U
    !f2py INTEGER, INTENT(hide)                  :: n, m
    REAL, DIMENSION(n,m,2) :: U
    INTEGER              :: n,m, i,j
    REAL :: x, y

    do i=0,n
      x = 2./n * i - 1
      do j=0,m
        y = 2./m * j - 1
        U(i,j,1) = 0.
        U(i,j,2) = exp(-(x*x + y*y) / (0.01))
      end do
    end do
  end subroutine


END MODULE
