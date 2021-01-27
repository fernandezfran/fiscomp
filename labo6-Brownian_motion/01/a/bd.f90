program bd_lj
!Dinámica browniana en un sistema de Lennard-Jones 
use precision, only  : pr        => dp
use randomnum, only  : randseed  => sgrnd
use initialization
use cell, only       : maps, links
use force, only      : forces    => forcelj126
use integrator, only : integrate => ermak
use comp_exp
implicit none
real(pr)             :: t, ti, tail, Epotm, Presm
integer              :: i, teq, trun, j, measurments

!### definición de parámetros iniciales
N      = 512
T0     = 2.5_pr
eta    = 2.87_pr
D0     = T0/pi3/eta
rho    = 0.8_pr
V      = real(N,pr)/rho
L      = V**(1._pr/3._pr)
dt     = 0.0005_pr
rcut   = 2.5_pr
ecut = 4._pr*(1._pr/(rcut**12) - 1._pr/(rcut**6))
tail = 4._pr*pi43*rho*( (2._pr/3._pr)*(1._pr/(rcut**9)) - (1._pr/(rcut**3)) )
Etail  = tail*real(N,pr)
Ptail  = tail*rho

semilla = 9786952 !semilla para integrar
call randseed(semilla)

allocate( x(N), y(N), z(N) )
!allocate( ix(N), iy(N), iz(N) ) !imagen de caja en la que está la part
allocate( fx(N), fy(N), fz(N) )

call maps
call initpossc               !posiciones iniciales e imagenes 0
call links
call forces                  !fuerzas iniciales
!
!definiciones para g(r)
nbin = 1000
allocate( g(nbin), gr(nbin) )
call gdr(0)      !inicializacion de g de r
!
!###
!
teq   = 20000
trun  = 20000
ti = 0._pr
!
open(36, file='trajectory.xyz', status='replace')
write(36,*) N; write(36,*) !write xyz
do j = 1, N; write(36,'(A,3(E15.6,2x))') 'Ar', x(j), y(j), z(j); enddo

open(47, file='thermo.dat', status='replace')
write(47,'(A,I4)') '# numero de partículas   : ', N
write(47,'(A,F6.4)') '# densidad               : ', rho
write(47,'(A,F6.4)') '# temperatura inicial    : ', T0
write(47,'(A,E9.3)') '# volumen                : ', V
write(47,'(A,F6.4)') '# rcut                   : ', rcut
write(47,'(A,F6.4)') '# paso temporal          : ', dt
write(47,'(A,I4)') '# pasos de equilibración : ', teq
write(47,'(A,I4)') '# pasos de medición      : ', trun
write(47,'(A)') '# t, Epot, Pres'

!### loop de equilibración
t = ti
do i = 1, teq
    !calculo nuevas posiciones
    call integrate
    !como cambie posiciones, renuevo vectores de cell_list y hochain
    call links
    !con datos actualizados, calculo fuerzas
    call forces
    
    if (mod(i,100) == 0) then !: write thermo & trajectory
        write(36,*) N; write(36,*)
        do j = 1, N; write(36,'(A,3(E15.6,2x))') 'Ar', x(j), y(j), z(j); enddo
        write(47,'(I6,3(E15.6,x))') i, t, Epot, Pres
    endif

    t = ti + real(i,pr)*dt
enddo

!### loop de medición
measurments = 0
Epotm = 0._pr
Presm = 0._pr
do i = teq + 1, teq + trun
    !calculo nuevas posiciones
    call integrate
    !como cambie posiciones, renuevo vectores de cell_list y hochain
    call links
    !con datos actualizados, calculo fuerzas
    call forces
    
    if (mod(i,100) == 0) then !: write thermo & trajectory & medir g(r)
        call gdr(1)
        write(36,*) N; write(36,*)
        do j = 1, N; write(36,'(A,3(E15.6,2x))') 'Ar', x(j), y(j), z(j); enddo
        write(47,'(I6,3(E15.6,x))') i, t, Epot, Pres
        Epotm = Epotm + Epot
        Presm = Presm + Pres
        measurments = measurments + 1
    endif

    t = ti + real(i,pr)*dt
enddo
call gdr(2)

Epotm = Epotm/real(measurments,pr)
Presm = Presm/real(measurments,pr)
write(47,'(A)') '#Valor medio de la energía potencial y de la presión'
write(47,'(2(E15.6,x))') Epotm, Presm

end program
