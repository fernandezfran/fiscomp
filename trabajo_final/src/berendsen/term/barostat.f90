module barostat
use precision, only : pr=>dp
use initialization, only : x, y, z, ix, iy, iz, pbc, dt, Pres, P0, tau_P, N, &
                         & L, V, tail, rho, Etail, Ptail
use cell, only           : maps
implicit none
contains

subroutine berendsen_press
implicit none
integer         :: i
real(pr)        :: coup_factor, mu

coup_factor = dt/tau_P
mu = (1._pr - coup_factor*(P0 - Pres))**(1._pr/3._pr)

do i = 1, N
    x(i) = mu*x(i)
    y(i) = mu*y(i)
    z(i) = mu*z(i)
enddo

L = L*mu
V = V*mu**3

do i = 1, N 
    !pbc por si el cambio en los parametros de la caja hizo que alguna
    !particula saliera de la celda de simulacion
    call pbc(x(i), ix(i))
    call pbc(y(i), iy(i))
    call pbc(z(i), iz(i))
enddo

!nueva correción de cola de energía y presión para nuevo tamaño de caja
rho = real(N,pr)/V
Etail = tail*real(N,pr)*rho
Ptail = tail*rho*rho

call maps !cambio tamaño de la caja => cambia linked list

end subroutine berendsen_press 

end module
