program mca
!
! integral de MC con important sampling
!
use precision, only:    pr => dp
use metodomc, only : intmc => impsampling 
implicit none
integer            :: pow, N, Ni, i
real(pr)           :: a, b, x, I_exacta, I_MC, error

!###################            datos de entrada           ########################
pow = 3                                              !exponente de la f(x) = x^pow
a   = 0._pr; b = 1._pr                               !extremos de integración
Ni  = 10                                             !cantidad de evaluaciones de f
!##################################################################################

open(17, file='b-error.dat', status='replace')

do i = 1, 1000
    N = Ni * i
    call intmc(N,a,b,x,f,g,I_MC)
    I_exacta = 1._pr/(real(pow + 1,pr))

    error = abs(I_MC - I_exacta)

    write(17,'(I5,3(E15.6,2x))') N, I_MC, I_exacta, error
enddo

contains


function f(x)
!funcion potencial que lee pow de arriba
real(pr), intent(in) :: x
real(pr)             :: f

f = x**pow

end function f



function g(x)
!funcion potencial que lee pow de arriba
real(pr), intent(in) :: x
real(pr)             :: k, g

k = 2
g = real(k+1,pr)*x**k 

end function g


end program
