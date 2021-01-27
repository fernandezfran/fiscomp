program bd_lj
!Dinámica browniana en un sistema de Lennard-Jones 
use precision, only  : pr        => dp
use randomnum, only  : randseed  => sgrnd, ran2
use initialization
use cell, only       : maps, links, hochain, cell_list, map
use force, only      : forces    => forcelj126
use integrator, only : integrate => ermak
use comp_exp
implicit none
real(pr)             :: t, ti, tail, rhoi, drho, changeseed
integer              :: i, teq, trun, j

!### definición de parámetros iniciales
N      = 512

dt     = 0.0005_pr
ti     = 0._pr
teq    = 20000
trun   = 200000

rcut   = 2.5_pr
ecut   = 4._pr*(1._pr/(rcut**12) - 1._pr/(rcut**6))

T0     = 2.5_pr
eta    = 2.87_pr
D0     = T0/pi3/eta
write(*,*) '#DO =', D0

allocate( x(N), y(N), z(N) )
allocate( ix(N), iy(N), iz(N) ) !imagen de caja en la que está la part
allocate( fx(N), fy(N), fz(N) )
Nmsd = 2000                     !definiciones para MSD
allocate( ntime(Nmsd), r2t(Nmsd) )
allocate( x0(N,Nmsd), y0(N,Nmsd), z0(N,Nmsd) )

semilla = 285469639 
rhoi = 1.05_pr
!rhof = 1.05_pr
drho = -0.05_pr != (rhof - rhoi)/cant. de ptos.
rho = rhoi - drho
write(*,*) '# semilla, rho, int(L/rcut)'
do j = 1, 18
    rho    = rho + drho
    V      = real(N,pr)/rho
    L      = V**(1._pr/3._pr)
    tail   = 4._pr*pi43*rho*( (2._pr/3._pr)*(1._pr/(rcut**9)) - (1._pr/(rcut**3)) )
    Etail  = tail*real(N,pr)
    Ptail  = tail*rho

    write(*,*) semilla, rho, int(L/rcut)
    
    call randseed(semilla) !con ran2 asigno la semilla del MT
    changeseed = ran2(semilla)

    call maps
    call initpossc               !posiciones iniciales e imagenes 0
    call links
    call forces                  !fuerzas iniciales
    
    call diffusion(0,10)
    
    !### loop de equilibración
    t = ti
    do i = 1, teq
        call integrate
        call links
        call forces
    
        t = ti + real(i,pr)*dt
    enddo
    
    !### loop de medición
    do i = teq + 1, teq + trun
        call integrate
        call links
        call forces
     
        if (mod(i,10) == 0) call diffusion(1,10)
     
        t = ti + real(i,pr)*dt
    enddo
       
    call diffusion(2,10)

    deallocate( hochain, cell_list, map )

enddo

end program
