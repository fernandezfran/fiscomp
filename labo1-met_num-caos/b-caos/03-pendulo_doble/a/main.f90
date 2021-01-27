program pendulo_doble
use precision, only : pr=>dp
use odes, only      : integrador=>rk4
integer               :: N, i, i_max
real(pr)              :: h, ti, tf, t, alfa, beta, gama, E
real(pr), allocatable :: x(:)
character(len=32)     :: frmt

!##############              datos iniciales              #####################

N = 4                                       !dimension del problema

alfa = 1._pr/3._pr                          !m2/m1
beta = 0.5_pr                               !L2/L1
gama = 0.5_pr                               !g/L1

h = 10._pr**(-3)                            !h_opt p/rk4

ti = 0._pr                                  !tiempo inicial
tf = 450._pr                                !       final

allocate(x(N))                             

!theta_1, omega_1, theta_2, omega_2 iniciales
x(:) = (/0._pr,0.332_pr,0._pr,0.845_pr/)
!x(:) = (/0._pr,sqrt(1.125_pr),0._pr,0._pr/)

!######       encabezado del archivo de salida con datos iniciales       ######

open(30, file='trayectorias.dat', status='replace')
write(30,*) '#Trayectoria del pendulo doble'
write(30,'(A,x,F6.4)') '#alfa (relación entre masas)', alfa
write(30,'(A,x,F6.4)') '#beta (relación entre largos)', beta
write(30,'(A,x,F6.4)') '#gama (relación entra g respecto a L)', gama
write(30,'(A,x,F6.4)') '#theta_1 inicial', x(1)
write(30,'(A,x,F6.4)') '#omega_1 inicial', x(2)
write(30,'(A,x,F6.4)') '#theta_2 inicial', x(3)
write(30,'(A,x,F6.4)') '#omega_2 inicial', x(4)
write(30,*) '# t [s], thetha_1, omega_1, theta_2, omega_2, Energia'
frmt = '(6(E21.7,2x))'

!##############################################################################


i_max = nint( (tf - ti) / h)
t = ti
call energy(x,E)
do i = 1, i_max

    write(30, frmt) t, x(1), x(2), x(3), x(4), E

    call integrador(N,h,t,x,f)
    call energy(x,E)
    
    t = ti + real(i,pr)*h  

enddo

contains

function f(tt,xx)
!
!funciones a integrar
!
implicit none
real(pr), intent(in) :: tt, xx(:)
real(pr)             :: f(size(xx))
real(pr)             :: sn, cs, den, fc1, fc2

!factores que se repiten
sn  = sin(xx(1) - xx(3))
cs  = cos(xx(1) - xx(3))

den = 1._pr + alfa*sn*sn
fc1 = (1._pr + alfa)*gama*sin(xx(1)) + alfa*beta*xx(4)*xx(4)*sn
fc2 = xx(2)*xx(2)*sn - gama*sin(xx(3))

!funciones
f(1) = xx(2)                                                     !theta_1
f(2) = -(fc1 + alfa*cs*fc2)/den                                  !omega_1
f(3) = xx(4)                                                     !theta_2
f(4) = ((1._pr + alfa)*fc2 + cs*fc1)/(beta*den)                  !omega_2

end function f

subroutine energy(xx,EE)
!
!energia del pendulo doble (L = T - V => cambio el signo de algunos terminos)
!                          (E = H = T + V)
implicit none
real(pr), intent(in)  :: xx(:)
real(pr), intent(out) :: EE
real(pr)              :: T1, T2, T12, V1, V2, TT, VV

T1  = (1._pr + alfa)*xx(2)*xx(2) 
T2  = alfa*beta*beta*xx(4)*xx(4) 
T12 = alfa*beta*xx(2)*xx(4)*cos(xx(1) - xx(3))
V1  = - (1._pr + alfa)*gama*cos(xx(1)) 
V2  = - alfa*beta*gama*cos(xx(3))

TT = (0.5_pr*(T1 + T2) + T12)
VV = V1 + V2

EE = TT + VV

end subroutine energy

end program
