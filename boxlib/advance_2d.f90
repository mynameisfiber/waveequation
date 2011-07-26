!
      subroutine advance (U, Unew, lo, hi, Ncomp, ng, dx, dt)
        implicit none
        integer lo(2),hi(2),Ncomp,ng, n, m
        double precision    U(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,Ncomp)
        double precision Unew(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,Ncomp)
        double precision dx,dt
        interface
          function dU(U, n, m, ng, dx)
            DOUBLE PRECISION, DIMENSION(n,m,2) :: U, dU
            INTEGER                :: n,m
            INTEGER                :: i,j, ng
            DOUBLE PRECISION       :: dx
          end function
        end interface
           

        n = hi(1) - lo(1) + 2 * ng
        m = hi(2) - lo(2) + 2 * ng

        Unew =         U                +         dt * dU(U   , n, m, ng, dx)
        Unew = 0.75  * U + 0.25  * Unew + 0.25  * dt * dU(Unew, n, m, ng, dx)
        U    = 1./3. * U + 2./3. * Unew + 2./3. * dt * dU(Unew, n, m, ng, dx)
      end subroutine


      FUNCTION dU(U, n, m, ng, dx) 
        implicit none
        DOUBLE PRECISION, DIMENSION(n,m,2) :: U, dU
        INTEGER                :: n,m
        INTEGER                :: i,j, ng
        DOUBLE PRECISION       :: dx
        interface
          function dU2dx2full(U, n, m, dx, dir) result(dU2dx2)
            DOUBLE PRECISION, DIMENSION(n,m) :: U, dU2dx2
            INTEGER              :: n,m,dir
            double precision     :: dx
          end function
        end interface

        dU(:,:,1) = U(:,:,2)
        dU(:,:,2) = dU2dx2full(U(:,:,1), n, m, dx, 1) + &
                    dU2dx2full(U(:,:,1), n, m, dx, 2) 

        return
      END FUNCTION

      FUNCTION dU2dx2full(U, n, m, dx, dir) RESULT(dU2dx2)
        implicit none
        DOUBLE PRECISION, DIMENSION(n,m) :: U, dU2dx2
        INTEGER              :: n,m,dir
        double precision     :: dx

        !dU2dx2 = (CSHIFT(U, 1, dir) + CSHIFT(U, -1, dir) - 2*U) / dx
        
        ! Assume we have 2 ghost zones.
        dU2dx2 = (-1./12. * (CSHIFT(U,2,dir) + CSHIFT(U,-2,dir)) + & 
                   4./ 3. * (CSHIFT(U,1,dir) + CSHIFT(U,-1,dir)) - &
                     2.5  * U                                      &
                 ) / (dx*dx)

        return
      END FUNCTION
