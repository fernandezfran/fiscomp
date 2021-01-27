program pendulo_doble
use precision, only : pr=>dp
use odes, only      : integrador=>rk4
real(pr), parameter   :: pi = 4._pr * atan(1._pr)
integer               :: N, i, i_max, j, j_max, k, k_max, ctrl
real(pr)              :: h, h1, h2, ti, tf, t, alfa, beta, gama
real(pr)              :: theta1, theta2, theta1_i, theta1_f, theta2_i, theta2_f
real(pr), allocatable :: x(:)
character(len=32)     :: frmt

!##############              datos iniciales              #####################

N = 4                                       !dimension del problema
!theta_1, omega_1, theta_2, omega_2
allocate(x(N))                             

alfa = 1._pr                                !m2/m1
beta = 1._pr                                !L2/L1
gama = 1._pr                                !g/L1
                                            
h = 1.d-3                                   !h_opt p/rk4

ti = 0._pr                                  !tiempo inicial
tf = 10000._pr                              !       final
i_max = nint( (tf - ti) / h)

theta1_i = -3._pr
theta1_f = 3._pr
j_max    = 600                              !cantidad de puntos theta1
h1       = (theta1_f - theta1_i)/real(j_max,pr)

theta2_i = 0._pr
theta2_f = 1._pr
k_max    = 100                              !cantidad de puntos theta2
h2       = (theta2_f - theta2_i)/real(k_max,pr)

!tamaño de la grilla grilla = j_max * k_max

!######       encabezado del archivo de salida con datos iniciales       ######

open(30, file='t_flip01.dat', status='replace')
write(30,*) '#Trayectoria del pendulo doble'
write(30,'(A,x,F6.4)') '#alfa (relación entre masas)', alfa
write(30,'(A,x,F6.4)') '#beta (relación entre largos)', beta
write(30,'(A,x,F6.4)') '#gama (relación entra g respecto a L)', gama
write(30,*) '# theta1, theta2, t_flip'
frmt = '(3(E21.7,2x))'

!##############################################################################

do j = 0, j_max
   
    write(*,*) j, 'que debe llegar a 600' 
    theta1 = theta1_i + h1*real(j,pr)
    
    ctrl = 0

    do k = 0, k_max

        theta2 = theta2_i + h2*real(k,pr)
       
        x(1) = theta1
        x(3) = theta2

        if ((2._pr*cos(x(1)) + cos(x(3))) .gt. 1._pr) cycle
        
        t = ti
        x(2) = 0._pr                          !velocidades iniciales igual a 0
        x(4) = 0._pr

        do i = 1, i_max
        
            call integrador(N,h,t,x,f)
            t = ti + real(i,pr)*h
        
            if ((abs(x(1)).gt.pi) .or. (abs(x(3)).gt.pi)) then
                write(30,frmt) theta1, theta2, t
                ctrl = 1
                exit
            endif

        enddo

    enddo

    if (ctrl .eq. 1) write(30,*)

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

end program
