program md_lj
!Dinámica molecular en un sistema de Lennard-Jones en el ensamble
! microcanónico con reescaleo de velocidades
use precision, only  : pr        => dp
use initialization
use cell, only       : maps, links
use force, only      : forces    => forcelj126
use integrator, only : integrate => velocityverlet 
implicit none
real(pr), parameter  :: pi43 = (16._pr/3._pr)*atan(1._pr)
real(pr)             :: ti, t, Tsf, tail
integer              :: i, teq, tscal, trun, j

!### definición de parámetros iniciales
N      = 4*6**3                                   !N debe ser 4*s**3
rho    = 0.8_pr
T0     = 1.1_pr
V      = real(N,pr)/rho
L      = V**(1._pr/3._pr)
dt     = 0.0005_pr
rcut   = 2.5_pr
ecut = 4._pr*(1._pr/(rcut**12) - 1._pr/(rcut**6))
tail = 4._pr*pi43*rho*( (2._pr/3._pr)*(1._pr/(rcut**9)) - (1._pr/(rcut**3)) )
Etail  = tail*real(N,pr)
Ptail  = tail*rho

allocate( x(N), y(N), z(N) )
allocate( vx(N), vy(N), vz(N) )
allocate( fx(N), fy(N), fz(N) )

call maps                    !mapa de celdas vecinas
call initposfcc              !posiciones iniciales
call links                   !linked list con las posiciones
call initvelrand             !velocidades iniciales
call forces                  !fuerzas iniciales
!
!###
!
teq   = 20000
tscal = 20
trun  = 20000
ti = 0._pr
!
open(45, file='thermo.dat', status='replace')
write(45,'(A,I4)') '# numero de partículas   : ', N
write(45,'(A,F6.4)') '# densidad               : ', rho
write(45,'(A,F6.4)') '# temperatura inicial    : ', T0
write(45,'(A,E9.3)') '# volumen                : ', V
write(45,'(A,F6.4)') '# rcut                   : ', rcut
write(45,'(A,F6.4)') '# paso temporal          : ', dt
write(45,'(A,I4)') '# pasos de equilibración : ', teq
write(45,'(A,I4)') '# pasos de medición      : ', trun
write(45,'(A,I4)') '# reescaleo de vel cada  : ', tscal
write(45,'(A)') '# t, Ekin, Epot, Etot, Temp, Pres'

open(36, file='trajectory.xyz', status='replace')
write(36,*) N; write(36,*) !write xyz
do j = 1, N; write(36,'(A,3(E15.6,2x))') 'Ar', x(j), y(j), z(j); enddo

!### loop de equilibración
do i = 1, teq
    call integrate    !con cálculo de fuerzas
    if (mod(i,tscal) == 0) then !reescaleo las velocidades
        Tsf = sqrt( T0/Temp )
        vx(:) = Tsf*vx(:)
        vy(:) = Tsf*vy(:)
        vz(:) = Tsf*vz(:)
    endif 
    if (mod(i,10) == 0) then !: write thermo & trajectory
        write(45,'(6(E15.6,x))') t, Ekin, Epot, Ekin + Epot, Temp, Pres
        write(36,*) N; write(36,*)
        do j = 1, N; write(36,'(A,3(E15.6,2x))') 'Ar', x(j), y(j), z(j); enddo
    endif
    t = ti + real(i,pr)*dt
enddo

!### main loop
do i = teq + 1, teq + trun
    call integrate     !con cálculo de fuerzas
    if (mod(i,10) == 0) then !: write thermo & trajectory
        write(45,'(6(E15.6,x))') t, Ekin, Epot, Ekin + Epot, Temp, Pres
        write(36,*) N; write(36,*)
        do j = 1, N; write(36,'(A,3(E15.6,2x))') 'Ar', x(j), y(j), z(j); enddo
    endif
    t = ti + real(i,pr)*dt
enddo

end program
