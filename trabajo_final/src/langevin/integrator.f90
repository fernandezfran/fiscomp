module integrator
!integración para obtener las posiciones y velocidades al tiempo siguiente
use precision, only : pr=> dp
use initialization, only : x, y, z, vx, vy, vz, fx, fy, fz, ix, iy, iz, dt, &
                         & Temp, T0, zeta, Ekin, N, pbc, pi
use cell, only  : links
use force, only : forces => forcelj126
implicit none
contains

subroutine velocityverlet
!algoritmo de velocity verlet
implicit none
integer                  :: i
real(pr)                 :: sumv2

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

sumv2 = 0._pr
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


subroutine langevinvv
! algoritmo tipo velocity verlet para termostato de langevin
! (podría pensarse como separarlo e introducirlo con flags dentro del vv de arriba) 
use randomnum, only : rand => grnd
implicit none
integer             :: i
real(pr)            :: abc, ac, aa, bb, cc, vxr, vyr, vzr, sumv2

abc = zeta * dt
ac = 2._pr + abc
aa = (2._pr - abc) / ac
bb = sqrt(0.5_pr * T0 * abc) 
cc = 2._pr * dt / ac

do i = 1, N
    vxr = gasdev()
    vyr = gasdev()
    vzr = gasdev()
    
    vx(i) = vx(i) + bb*vxr + 0.5_pr*fx(i)*dt
    vy(i) = vy(i) + bb*vyr + 0.5_pr*fy(i)*dt
    vz(i) = vz(i) + bb*vzr + 0.5_pr*fz(i)*dt

    x(i) = x(i) + cc*vx(i)
    y(i) = y(i) + cc*vy(i)
    z(i) = z(i) + cc*vz(i)

    call pbc(x(i), ix(i))
    call pbc(y(i), iy(i))
    call pbc(z(i), iz(i))
enddo

call links
call forces

sumv2 = 0.0
do i = 1, N
    vxr = gasdev()
    vyr = gasdev()
    vzr = gasdev()
    
    vx(i) = aa*vx(i) + bb*vxr + 0.5_pr*fx(i)*dt
    vy(i) = aa*vy(i) + bb*vyr + 0.5_pr*fy(i)*dt
    vz(i) = aa*vz(i) + bb*vzr + 0.5_pr*fz(i)*dt

    sumv2 = sumv2 + vx(i)*vx(i) + vy(i)*vy(i) + vz(i)*vz(i)
enddo
Temp = sumv2/real(3*N,pr)
Ekin = 0.5_pr*sumv2


    contains
        
        function gasdev()
        ! método de box-muller
        !    se considera x0 = 0 y sigma = 1
        implicit none
        real(pr)               :: gasdev, gset
        real(pr)               :: u1, u2, fs, ang
        integer                :: iset = 0
        save iset, gset
        
        if (iset == 0) then
            u1  = rand(); u2  = rand()
            fs  = sqrt(-2._pr*log(u1))
            ang = 2._pr*pi*u2
        
            gset   = fs*cos(ang)
            gasdev = fs*sin(ang)
        
            iset = 1
        else
            gasdev = gset
            iset = 0
        endif

        endfunction gasdev

end subroutine langevinvv

end module
