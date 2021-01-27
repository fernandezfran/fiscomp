program pendulo_doble
use precision, only : pr=>dp
use odes, only      : integrador=>rk4
real(pr), parameter   :: alfa = 0.05_pr    !acople
integer               :: N, i, i_max
real(pr)              :: h, ti, tf, t, E
real(pr), allocatable :: x(:)

!##############              datos iniciales              #####################

N = 4                                       !dimension del problema

h = 10._pr**(-3)                            !h_opt p/rk4

ti = 0._pr                                  !tiempo inicial
tf = 200._pr                                !       final

!probar energías 5, 20, 100
E = 5._pr                                   !energía inicial

!#############              condicion inicial              ####################

allocate(x(N))

!x(:) = q1, p1, q2, p2
x(:) = (/2._pr,0._pr,0._pr,sqrt(2._pr*E - 4._pr)/)

!##############################################################################


i_max = nint( (tf - ti) / h)
t = ti
do i = 1, i_max

    call integrador(N,h,t,x,f)
    
    t = ti + real(i,pr)*h

enddo


contains

function f(tt,xx)
!
!hamiltoniano de Pullen-Edmonds (dos osciladores de la misma frecuencia w = 1
! y masa m=2, elegir unidades tal que...) acoplados con un factor alfa.
!
! pe(:) = (/q1,p1,q2,p2/)
!
implicit none
real(pr), intent(in) :: tt, xx(:)
real(pr)             :: f(size(xx))

f(1) = xx(2) 
f(2) = xx(1) * (1._pr + 2._pr*alfa*xx(3)*xx(3))
f(3) = xx(4)
f(4) = xx(3) * (1._pr + 2._pr*alfa*xx(1)*xx(1))

end function


end program
