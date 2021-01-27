module integrator
!integración para obtener las posiciones y velocidades al tiempo siguiente
use precision, only : pr=> dp
use initialization, only : x, y, z, vx, vy, vz, fx, fy, fz, ix, iy, iz, dt, &
                         & Temp, Ekin, N, pbc
use cell, only  : links
use force, only : forces => forcelj126
implicit none
contains

subroutine velocityverlet
!algoritmo de velocity verlet
implicit none
integer                  :: i
real(pr)                 :: sumv2

sumv2 = 0._pr

do i = 1, N
    !actualizo posiciones
    x(i) = x(i) + vx(i)*dt + 0.5_pr*fx(i)*dt*dt
    y(i) = y(i) + vy(i)*dt + 0.5_pr*fy(i)*dt*dt
    z(i) = z(i) + vz(i)*dt + 0.5_pr*fz(i)*dt*dt

    !hago condiciones periódicas de contorno
    call pbc(x(i), ix(i))
    call pbc(y(i), iy(i))
    call pbc(z(i), iz(i))
    
    !actualizo una parte de las velocidades
    vx(i) = vx(i) + 0.5_pr*fx(i)*dt 
    vy(i) = vy(i) + 0.5_pr*fy(i)*dt 
    vz(i) = vz(i) + 0.5_pr*fz(i)*dt 
enddo

!actualizo fuerzas
call links
call forces

do i = 1, N
    !actualizo lo que falta de las velocidades
    vx(i) = vx(i) + 0.5_pr*fx(i)*dt 
    vy(i) = vy(i) + 0.5_pr*fy(i)*dt 
    vz(i) = vz(i) + 0.5_pr*fz(i)*dt 

    !guardo algunos factores de interes
    sumv2 = sumv2 + vx(i)*vx(i) + vy(i)*vy(i) + vz(i)*vz(i)
enddo

Temp = sumv2/real(3*N,pr)
Ekin = 0.5_pr*sumv2

end subroutine velocityverlet

end module
