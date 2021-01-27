module metodomc
!
! modulo con métodos de integración de montecarlo
!
use precision, only:     pr => dp
use randomnum, only:   rand => rmzran
use dplpmc, only   : nurand => ldp
!                     ^ no-uniforme rand
contains

subroutine montecarlo(N,a,b,x,f,I)
!método de integración de montecarlo unidimensional
!N: cantidad de puntos
!a y b: extremos de integracion
!f(x): funcion a integrar
!I: integral de MC (resultado)
implicit none
integer, intent(in)   :: N
real(pr), intent(in)  :: a, b, x 
real(pr), intent(out) :: I
real(pr)              :: V, rn
integer               :: k
interface
    function f(x)
    use precision, only:   pr => dp
    real(pr), intent(in) :: x
    real(pr)             :: f
    end function f
end interface

I = 0._pr
V = (b - a)/real(N,pr)
do k = 1, N
    rn = rand()
    I  = I + f(rn) 
enddo
I = I*V

end subroutine


subroutine impsampling(N,a,b,x,f,g,I)
!método de integración de montecarlo unidimensional con important sampling
!N: cantidad de puntos
!f(x): funcion a integrar
!g(x): distribución de probabilidad
!I: integral de MC (resultado)
implicit none
integer, intent(in)   :: N
real(pr), intent(in)  :: a,b,x 
real(pr), intent(out) :: I
real(pr)              :: rn
integer               :: j
interface
    function f(x)
    use precision, only:   pr => dp
    real(pr), intent(in) :: x
    real(pr)             :: f
    end function f
end interface
interface
    function g(x)
    use precision, only:   pr => dp
    real(pr), intent(in) :: x
    real(pr)             :: k, g
    end function g
end interface

I = 0._pr
do j = 1, N
    rn = nurand()
    I  = I + f(rn)/g(rn)
enddo
I = I/real(N,pr)

end subroutine


end module
