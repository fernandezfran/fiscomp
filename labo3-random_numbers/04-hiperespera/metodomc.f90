module metodomc
!
! modulo con métodos de integración de montecarlo
!
use precision, only:     pr => dp
use randomnum, only:   rand => grnd
!                     ^ no-uniforme rand
contains

subroutine montecarlo(ndim,N,V,x,f,I,var)
!método de integración de montecarlo multidimensional
!ndim : dimension
!N    : cantidad de puntos
!V    : volumen de integración
!f(x) : funcion a integrar
!I    : integral de MC (resultado)
!sigma: varianza
implicit none
integer, intent(in)   :: N, ndim
real(pr), intent(in)  :: x(ndim-1), V
real(pr), intent(out) :: I, var
real(pr)              :: rn(ndim-1), sigma, efe, f2
integer               :: k, j
interface
    function f(x)
    use precision, only:   pr => dp
    real(pr), intent(in) :: x(:)
    real(pr)             :: f
    end function f
end interface

I  = 0._pr
f2 = 0._pr
do k = 1, N
    do j = 1, ndim - 1
        rn(j) = 2._pr * rand() - 1._pr
    enddo
    efe = f(rn(:))
    I   = I  + efe
    f2  = f2 + efe*efe
enddo
I     = I/real(N,pr)
sigma = f2/real(N,pr) - I*I 

I   = V*I
var = V*sqrt(sigma/real(N,pr))

end subroutine

end module
