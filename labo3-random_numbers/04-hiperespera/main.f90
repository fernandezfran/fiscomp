program hiperesfera
!
! funciona hasta dimensión 20 (con un error de e-2 en este caso)
!  podría mejorarse con un `importance sampling', para que los número aleatorios
!  le peguen a la hiperesfera (a medida que crece la dimensión disminuye su volumen
!                              en relación al hipercubo).
!
use precision, only :    pr=>dp
use metodomc, only  : intmc=>montecarlo
implicit none
real(pr), parameter   :: pi = 4._pr * atan(1._pr)
integer               :: ndim, npts
real(pr)              :: V, vol, V_MC, Var
real(pr), allocatable :: x(:)
real(pr)              :: ti, tf

write(*,*) 'dimensión de la hiperesfera:'
read(*,*) ndim
npts = 2**25
allocate( x(ndim-1) )

!Hipervolumen exacto
if (mod(ndim,2) == 0) then    !si ndim es par
    V = (pi**(ndim/2))/gamma(real(ndim/2 + 1,pr))
else                          !si ndim es impar
    V = (2._pr**ndim)*(pi**((ndim-1)/2))*gamma(real((ndim-1)/2 + 1,pr))/gamma(real(ndim+1,pr))
endif
write(*,*) 'el volumen exacto es'
write(*,*) V
write(*,*)

vol = 2._pr**(ndim - 1)    !volumen de integración
write(*,*) 'volumen relativo de la hiperesfera al hipercubo V[hiperesfera/hipercubo]'
write(*,*) V/vol
write(*,*)

call cpu_time(ti)
call intmc(ndim,npts,vol,x,f,V_MC,Var)
call cpu_time(tf)

write(*,*) 'el volumen cálculado con Monte Carlo es'
write(*,*) V_MC
write(*,*) 'con varianza'
write(*,*) Var
write(*,*)
write(*,*)
write(*,*) 'tiempo de Monte Carlo [s]', tf-ti


contains


function f(x)
implicit none
real(pr), intent(in) :: x(:)
real(pr)             :: f
real(pr)             :: sum2
integer              :: l

f    = 0._pr
sum2 = dot_product(x,x)
if (sum2 <= 1._pr) f = 2._pr*sqrt(1._pr - sum2)

end function f


end program
