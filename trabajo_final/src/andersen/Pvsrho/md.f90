program md_lj
!Dinámica molecular en un sistema de Lennard-Jones en el ensamble
! microcanónico con reescaleo de velocidades
use precision, only  : pr        => dp
use initialization
use force, only      : forces    => forcelj126
use integrator, only : integrate => velocityverlet
use cell, only       : maps, links
use thermostat, only : thermalize=> andersen_temp
use comp_exp, only   : veldist
implicit none
real(pr)             :: t, ti
integer              :: i, teq, trun, j

!### definición de parámetros iniciales
open(15, file='in.lj-data', status='old')
read(15,*) N
read(15,*) rho
read(15,*) T0
read(15,*) dt
read(15,*) rcut
read(15,*) teq
read(15,*) trun
read(15,*) nu
V      = real(N,pr)/rho
L      = V**(1._pr/3._pr)
ecut = 4._pr*(1._pr/(rcut**12) - 1._pr/(rcut**6))
tail = 4._pr*pi43*( (2._pr/3._pr)*(1._pr/(rcut**9)) - (1._pr/(rcut**3)) )
Etail  = tail*real(N,pr)*rho
Ptail  = tail*rho*rho

allocate( x(N), y(N), z(N) )
allocate( ix(N), iy(N), iz(N) ) !imagen de caja en la que está la part
allocate( vx(N), vy(N), vz(N) )
allocate( fx(N), fy(N), fz(N) )

call initposfcc              !posiciones iniciales e imagenes 0
call maps
call links
call initvelrand             !velocidades iniciales
call forces                  !fuerzas iniciales

!open(36, file='trajectory.xyz', status='replace')
!write(36,*) N; write(36,*) !write xyz
!do j = 1, N; write(36,'(A,3(E15.6,2x))') 'Ar', x(j), y(j), z(j); enddo

open(47, file='thermo.dat', status='replace')
write(47,'(A,I4)') '# numero de partículas   : ', N
write(47,'(A,F6.4)') '# densidad               : ', rho
write(47,'(A,F6.4)') '# temperatura inicial    : ', T0
write(47,'(A,E9.3)') '# volumen                : ', V
write(47,'(A,F6.4)') '# rcut                   : ', rcut
write(47,'(A,F6.4)') '# paso temporal          : ', dt
write(47,'(A,I6)') '# pasos de equilibración : ', teq
write(47,'(A,I6)') '# pasos de medición      : ', trun
write(47,'(A)') '# t, Ekin, Epot, Etot, Temp, Pres, Vol, L, rho'

!### loop de equilibración
ti = 0._pr
t = ti
print*, "por termalizar"
do i = 1, teq
    call integrate    !con cálculo de fuerzas
    call thermalize

    !if (mod(i,100) == 0) then !: write thermo & trajectory
    !    !write(36,*) N; write(36,*)
    !    !do j = 1, N; write(36,'(A,3(E15.6,2x))') 'Ar', x(j), y(j), z(j); enddo
    !    write(47,'(9(E15.6,x))') t, Ekin, Epot, Ekin + Epot, Temp, Pres, V, L, real(N,pr)/V
    !endif

    t = ti + real(i,pr)*dt
enddo

!### loop de medición
print*, "por medir"
call veldist(40, 0)
do i = teq + 1, teq + trun
    call integrate
    call thermalize

    if (mod(i,100) == 0) then !: write thermo & trajectory
        call veldist(40, 1)
        !write(36,*) N; write(36,*)
        !do j = 1, N; write(36,'(A,3(E15.6,2x))') 'Ar', x(j), y(j), z(j); enddo
        write(47,'(9(E15.6,x))') t, Ekin, Epot, Ekin + Epot, Temp, Pres, V, L, real(N,pr)/V
    endif

    t = ti + real(i,pr)*dt
enddo

call veldist(40, 2)

end program
