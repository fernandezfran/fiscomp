program md_lj
!Dinámica molecular en un sistema de Lennard-Jones en el ensamble
! microcanónico con reescaleo de velocidades
use precision, only  : pr        => dp
use initialization
use cell, only       : maps, links, hochain, cell_list, map
use force, only      : forces    => forcelj126
use integrator, only : integrate => velocityverlet 
implicit none
real(pr), parameter  :: pi43 = (16._pr/3._pr)*atan(1._pr)
real(pr)             :: ti, t, Tsf, tail
real(pr)             :: tstart, tfinish
integer              :: i, teq, tscal, trun, j, s

!### definición de parámetros iniciales
rho    = 0.8_pr
T0     = 1.1_pr
tscal  = 20
dt     = 0.0005_pr
teq    = 5._pr/dt
trun   = 10._pr/dt
rcut   = 2.5_pr
ecut = 4._pr*(1._pr/(rcut**12) - 1._pr/(rcut**6))
tail = 4._pr*pi43*rho*( (2._pr/3._pr)*(1._pr/(rcut**9)) - (1._pr/(rcut**3)) )
Ptail  = tail*rho

open(26, file='time.dat', status='replace')
write(26,'(A)') '# N, t_elapsed [sec]'

!### loop en N
ti = 0._pr
do s = 5, 10
    call cpu_time(tstart)

    N = 4*s**3
    write(*,'(I5)') N
    V = real(N,pr)/rho
    L = V**(1._pr/3._pr)
    Etail  = tail*real(N,pr)

    allocate( x(N), y(N), z(N) )
    allocate( vx(N), vy(N), vz(N) )
    allocate( fx(N), fy(N), fz(N) )

    call maps                    !mapa de celdas vecinas
    call initposfcc              !posiciones iniciales
    call links                   !linked list con las posiciones
    call initvelrand             !velocidades iniciales
    call forces                  !fuerzas iniciales

    !### loop de equilibración
    t = ti
    do i = 1, teq
        call integrate
    
        if (mod(i,tscal) == 0) then !reescaleo las velocidades
            Tsf = sqrt( T0/Temp )
            vx(:) = Tsf*vx(:)
            vy(:) = Tsf*vy(:)
            vz(:) = Tsf*vz(:)
        endif 
    
        t = ti + real(i,pr)*dt
    enddo

    !### loop de medicion
    do i = teq + 1, teq + trun

        call integrate

        t = ti + real(i,pr)*dt

    enddo

    deallocate(x, y, z, vx, vy, vz, fx, fy, fz)
    deallocate( hochain, cell_list, map)

    call cpu_time(tfinish)

    write(26,'(I5,x,E15.6)') N, tfinish-tstart
    write(*,'(I5,x,E15.6)') N, tfinish-tstart
enddo

end program
