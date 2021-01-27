module matrix_solver
!
!subrutinas para resolver matrices
!
use precision, only: pr=>dp

contains

subroutine tridiag(N,a,d,c,b,x) 
!
!subrutina que resuelve un sistema de ec. lineales correspondiente a una matriz
!tri-diagonal
!
!Ax = b,   donde A es N * N y tiene valores `d' en la diagonal, `c' en la
!                              diagonal superior y `a' en la diagonal inferior
!
implicit none
integer, intent(in)   :: N
real(pr), intent(in)  :: d(N), a(2:N), c(1:N-1), b(N)
real(pr), intent(out) :: x(N)
integer               :: i
real(pr)              :: h(N), p(N), den

!calculo h(i) y p(i)
h(1) = c(1)/d(1) ; p(1) = b(1)/d(1)
do i = 2, N-1
    den = d(i) - a(i)*h(i-1)

    h(i) = c(i)/den
    p(i) = (b(i) - a(i)*p(i-1))/den
enddo
p(N) = (b(N) - a(N)*p(N-1))/(d(N) - a(N)*h(N-1))

!trivial back sustitution
x(N) = p(N)
do i = N-1,1,-1
    x(i) = p(i) - h(i)*x(i+1)
enddo

return
end subroutine tridiag

end module matrix_solver
