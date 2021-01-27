program test
use precision, only     : pr=>dp
use matrix_solver, only : tridiag
implicit none
integer                 :: n
real(pr), allocatable   :: a(:), d(:), c(:), b(:), x(:)

n = 4
allocate(x(n), a(2:n), d(n), c(n-1), b(n)) !los arreglos empiezan en 1 por default

!inicializo valores
d(:) = 1._pr
a(:) = 3._pr
c(:) = 2._pr
b(:) = 4._pr

write(*,*) 'resuelvo matriz la sgte matriz'
write(*,*) nint(d(1)), nint(c(1)), 0         , 0        
write(*,*) nint(a(2)), nint(d(2)), nint(c(2)), 0         
write(*,*) 0         , nint(a(3)), nint(d(3)), nint(c(3))
write(*,*) 0         , 0         , nint(a(4)), nint(d(4))
write(*,*) 'igualada a 4'

call tridiag(n,a,d,c,b,x)

write(*,*) 'según google, el resultado deberia ser'
write(*,*) 'x(1) = -20/19 (=', -20._pr/19._pr,')', x(1)
write(*,*) 'x(2) = 48/19 (=', 48._pr/19._pr,')', x(2)
write(*,*) 'x(3) = 44/19 (=', -20._pr/19._pr,')', x(3)
write(*,*) 'x(4) = -56/19 (=', -56._pr/19._pr,')', x(4)


end program
