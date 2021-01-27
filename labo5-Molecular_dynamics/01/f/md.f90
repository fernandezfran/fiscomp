program md_lj
!Dinámica molecular en un sistema de Lennard-Jones en el ensamble
! microcanónico con reescaleo de velocidades
use precision, only  : pr        => dp
use initialization
use force, only      : forces    => forcelj126
use integrator, only : integrate => velocityverlet 
implicit none
real(pr), parameter  :: pi43 = (16._pr/3._pr)*atan(1._pr)
real(pr)             :: ti, t, Tsf, tail
integer              :: i, teq, tscal, trun, j, k
real(pr)             :: Ekinm, Ekin2, Epotm, Epot2, Etotm, Etot2, sigkin, sigpot, sigtot, Etot

!### definición de parámetros iniciales
N      = 256                                    !N debe ser 4*s**3
rho    = 0.8_pr
T0     = 1.1_pr
V      = real(N,pr)/rho
L      = V**(1._pr/3._pr)

allocate( x(N), y(N), z(N) )
allocate( vx(N), vy(N), vz(N) )
allocate( fx(N), fy(N), fz(N) )

dt    = 0.0005_pr
tscal = 20
teq  = 5._pr/dt 
trun = teq
ti = 0._pr

!### loop en rcut
!
open(45, file='thermo.dat', status='replace')
write(45,'(A)') '#rcut, sigkin, sigpot, sigtot, sigtot/sigkin, sigtot/sigpot'
do k = 1, 9
    rcut = 0.5_pr*real(k+1,pr) !toma valores, cada 0.5, entre 1 y 5
    write(*,'(E15.6)') rcut

    ecut = 4._pr*(1._pr/rcut**12 - 1._pr/rcut**6)
    tail = 4._pr*pi43*rho*( (2._pr/3._pr)*(1._pr/rcut**9) - (1._pr/rcut**3) )

    Etail  = tail*real(N,pr)
    Ptail  = tail*rho

    call initposfcc              !posiciones iniciales
    call initvelrand             !velocidades iniciales
    call forces                  !fuerzas iniciales

    !### loop de equilibración
    do i = 1, teq
        call integrate    !con cálculo de fuerzas
        if (mod(i,tscal) == 0) then !reescaleo las velocidades
            Tsf = sqrt( T0/Temp )
            vx(:) = Tsf*vx(:)
            vy(:) = Tsf*vy(:)
            vz(:) = Tsf*vz(:)
        endif 
        t = ti + real(i,pr)*dt
    enddo

    !### loop de medición
    Ekinm = 0._pr; Epotm = 0._pr; Etotm = 0._pr
    Ekin2 = 0._pr; Epot2 = 0._pr; Etot2 = 0._pr

    do i = teq + 1, teq + trun
        call integrate

        Ekinm = Ekinm + Ekin
        Ekin2 = Ekin2 + Ekin*Ekin

        Epotm = Epotm + Epot
        Epot2 = Epot2 + Epot*Epot

        Etot  = Ekin + Epot
        Etotm = Etotm + Etot
        Etot2 = Etot2 + Etot*Etot 

        t = ti + real(i,pr)*dt
    enddo

    Etotm = Etotm/real(trun,pr) 
    Epotm = Epotm/real(trun,pr) 
    Ekinm = Ekinm/real(trun,pr) 
    Etot2 = Etot2/real(trun,pr) 
    Epot2 = Epot2/real(trun,pr) 
    Ekin2 = Ekin2/real(trun,pr) 

    sigtot = sqrt(abs(Etot2 - Etotm*Etotm))
    sigkin = sqrt(abs(Ekin2 - Ekinm*Ekinm))
    sigpot = sqrt(abs(Epot2 - Epotm*Epotm))

    write(45,'(6(E15.6,x))') rcut, sigkin, sigpot, sigtot, sigtot/sigkin, sigtot/sigpot
enddo


end program
