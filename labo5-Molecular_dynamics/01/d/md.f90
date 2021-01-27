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
integer              :: i, teq, tscal, trun, j
integer                :: ibin, nbin, k
integer, allocatable   :: hist(:)
real(pr)               :: dx

!### definición de parámetros iniciales
N      = 256                                    !N debe ser 4*s**3
rho    = 0.8_pr
T0     = 1.1_pr
V      = real(N,pr)/rho
L      = V**(1._pr/3._pr)
dt     = 0.005_pr
rcut   = 2.5_pr
ecut = 4._pr*(1._pr/(rcut**12) - 1._pr/(rcut**6))
tail = 4._pr*pi43*rho*( (2._pr/3._pr)*(1._pr/(rcut**9)) - (1._pr/(rcut**3)) )
Etail  = tail*real(N,pr)
Ptail  = tail*rho

allocate( x(N), y(N), z(N) )
allocate( vx(N), vy(N), vz(N) )
allocate( fx(N), fy(N), fz(N) )

call initposfcc              !posiciones iniciales
call initvelrand             !velocidades iniciales
call forces                  !fuerzas iniciales
!
!###
!
teq   = 1000
tscal = 20
trun  = 1000
ti = 0._pr
!

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
        write(36,*) N; write(36,*)
        do j = 1, N; write(36,'(A,3(E15.6,2x))') 'Ar', x(j), y(j), z(j); enddo
    endif
    t = ti + real(i,pr)*dt
enddo

!### histograma de velocidades
!
open(64, file='histo.dat', status='replace')
nbin = 100         !100
dx   = 10._pr/nbin !0.1_pr
allocate( hist(nbin) )
! vx
hist(:) = 0
do j = 1, N
   ibin = nint(vx(j)/dx)
   if (ibin >= -nbin/2) hist(ibin + nbin/2) = hist(ibin + nbin/2) + 1
enddo

do k = 1, nbin
    write(64,'(2(E15.6,x))') real(k-nbin/2,pr)*dx, real(hist(k),pr)/(real(N,pr)*dx)
enddo
write(64,*); write(64,*)

!vy
hist(:) = 0
do j = 1, N
   ibin = nint(vy(j)/dx)
   if (ibin >= -nbin/2) hist(ibin + nbin/2) = hist(ibin + nbin/2) + 1
enddo

do k = 1, nbin
    write(64,'(2(E15.6,x))') real(k-nbin/2,pr)*dx, real(hist(k),pr)/(real(N,pr)*dx)
enddo
write(64,*); write(64,*)

!vz
hist(:) = 0
do j = 1, N
   ibin = nint(vz(j)/dx)
   if (ibin >= -nbin/2) hist(ibin + nbin/2) = hist(ibin + nbin/2) + 1
enddo

do k = 1, nbin
    write(64,'(2(E15.6,x))') real(k-nbin/2,pr)*dx, real(hist(k),pr)/(real(N,pr)*dx)
enddo


end program
