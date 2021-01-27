program pendulo_doble_c
!
!espectro de potencia para el pendulo doble con las condiciones del inciso (a)
!
use, intrinsic     :: iso_c_binding
use precision, only : pr=>dp
use odes, only      : integrador=>rk4
implicit none
include "/usr/include/fftw3.f03"
real(pr), parameter                    :: pi = 4._pr * atan(1._pr)
real(C_double), allocatable            :: fx(:)
complex(C_double_complex), allocatable :: tfx(:)
type(C_ptr)                            :: plan_rc
integer                                :: N, i, i_max
real(pr)                               :: h, ti, tf, t, alfa, beta, gama, E, factor, rango
real(pr), allocatable                  :: x(:)
character(len=32)                      :: frmt

!##############              datos iniciales              #####################

N = 4                                       !dimension del problema
allocate(x(N))                             

alfa = 1._pr/3._pr                          !m2/m1
beta = 0.5_pr                               !L2/L1
gama = 0.5_pr                               !g/L1

h = 10._pr**(-3)                            !h_opt p/rk4

ti = 0._pr                                  !tiempo inicial
tf = 450._pr                                !       final

i_max = nint( (tf - ti) / h)

!theta_1, omega_1, theta_2, omega_2 iniciales
!x(:) = (/0._pr,0.332_pr,0._pr,0.845_pr/)
x(:) = (/0._pr,sqrt(1.125_pr),0._pr,0._pr/)

!######          definiciones necesarias para realizar la FFTW           ######

rango  = tf - ti
factor = 2._pr*pi/rango                     !frecuencia de Nyquist

allocate( fx(i_max), tfx(i_max/2 + 1) )
plan_rc = fftw_plan_dft_r2c_1d(i_max, fx, tfx, FFTW_MEASURE)

!######       encabezado del archivo de salida con datos iniciales       ######

open(79, file='power_spectrum.dat', status='replace')
write(79,*) '#Espectro de potencia'
write(79,'(A,x,F6.4)') '#alfa (relación entre masas)', alfa
write(79,'(A,x,F6.4)') '#beta (relación entre largos)', beta
write(79,'(A,x,F6.4)') '#gama (relación entra g respecto a L)', gama
write(79,'(A,x,F6.4)') '#theta_1 inicial', x(1)
write(79,'(A,x,F6.4)') '#omega_1 inicial', x(2)
write(79,'(A,x,F6.4)') '#theta_2 inicial', x(3)
write(79,'(A,x,F6.4)') '#omega_2 inicial', x(4)
write(79,*) '# frecuencia, Re(tfx), Im(tfx)'
frmt = '(3(E21.7,2x))'

!##############################################################################


t = ti
do i = 1, i_max

    call integrador(N,h,t,x,f)
    
    t = ti + real(i,pr)*h  
    
    fx(i) = x(3)

enddo

call fftw_execute_dft_r2c(plan_rc, fx, tfx)

do i = 1, i_max/2 + 1
    write(79,frmt) factor*real(i-1,pr), real(tfx(i))/real(i_max,pr), &
                          & aimag(tfx(i))/real(i_max,pr)
enddo
close(79)

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
