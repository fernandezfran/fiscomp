program pendulo_doble_b
!
!secciones de Poincare para el pendulo doble, fijando la energía, fijandome en
!el plano theta_1 = 0 y omega_1 > 0 para un barrido en p2
!
use precision, only : pr=>dp
use odes, only      : integrador=>rk4
integer               :: N, i, i_max, j, j_max, sgn
real(pr)              :: h, ti, tf, t, alfa, beta, gama
real(pr)              :: E, p1, p2, p2max, delta_p
real(pr), allocatable :: x(:)
character(len=32)     :: frmt

!##############              datos iniciales              #####################

N = 4                                       !dimension del problema
allocate(x(N)) 
!x(:) = (/theta_1, omega_1, theta_2, omega_2/)
x(:) = 0._pr

alfa = 1._pr/3._pr                          !m2/m1
beta = 0.5_pr                               !L2/L1
gama = 0.5_pr                               !g/L1

h = 10._pr**(-3)                            !h_opt p/rk4

ti = 0._pr                                  !tiempo inicial
tf = 5000._pr                               !       final

i_max = nint( (tf - ti) / h)                !cantidad de pasos RK4
j_max = 50                                  !cantidad de pts p2

write(*,*) 'ingrese el valor de energía'
read(*,*) E                                 !valor inicial de la energia


!######       encabezado del archivo de salida con datos iniciales       ######

open(30, file='poincare.dat', status='replace')
write(30,*) '#Trayectoria del pendulo doble'
write(30,'(A,x,F6.4)') '#alfa (relación entre masas)', alfa
write(30,'(A,x,F6.4)') '#beta (relación entre largos)', beta
write(30,'(A,x,F6.4)') '#gama (relación entra g respecto a L)', gama
write(30,*) '# theta2, p2'
frmt = '(2(E21.7,2x))'

!##############################################################################

p2max = g2(x(3))
write(30,*) '#p2max = ', p2max
delta_p = 2._pr*p2max/real(j_max,pr)

do j = 0, j_max

    !vuelvo a fijar las posiciones en cero
    x(:) = 0._pr

    !fijo p2 equispaciado entre -p2max y p2max. Cálculo respectivo de p1
    p2 = -p2max + real(j,pr)*delta_p
    p1 = g1(x(1), x(3), p2)

    !velocidades iniciales
    x(2) = (p1 - cos(x(1) - x(3))*p2/beta)/(1 + alfa*sin(x(1) - x(2))**2)
    x(4) = p2/(alfa*beta*beta) - cos(x(1) - x(3))*x(2)/beta

    !integro temporalmente
    t = ti
    do i = 1, i_max

        call integrador(N,h,t,x,f)
    
        t = ti + real(i,pr)*h

        if ((sgn .ne. sign(1._pr,x(1))) .and. (x(2) .gt. 0._pr)) then
            p2 = alfa*beta*(beta*x(4) + x(2)*cos(x(1) - x(3)))
            write(30,frmt) x(3), p2
        endif
    
        sgn = sign(1._pr,x(1))
   
    enddo
    write(30,*)
    write(30,*)
enddo


contains


!OJO!
!la funciones que vienen a continuación leen parametros que están definidos en el
!main, no se las puede pasar, tal cual están, a un modulo aparte


function g1(theta1,theta2,pp2)
!
!función que devuelve p1
!
implicit none
real(pr), intent(in) :: theta1, theta2, pp2
real(pr)             :: g1

g1 = E + gama*((1._pr + alfa)*cos(theta1) + beta*alfa*cos(theta2))
g1 = g1 * 2._pr*beta*beta*alfa - pp2*pp2
g1 = g1 * (2._pr + alfa - alfa*cos(2._pr*(theta1-theta2)))
g1 = g1 * (0.5_pr/alfa)

g1 =(sqrt(g1) + pp2)/beta

end function g1


function g2(theta2)
!
!función que devuelve p2max
!
implicit none
real(pr), intent(in) :: theta2
real(pr)             :: g2

g2 = E + (1._pr + alfa)*gama + alfa*beta*gama*cos(theta2)
g2 = g2 * 2._pr*alfa*beta*beta*(1._pr + alfa*sin(theta2)**2)
g2 = g2 / (1._pr + alfa)

g2 = sqrt(g2)

end function g2


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

end program
