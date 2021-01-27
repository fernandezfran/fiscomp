module force
!modulo para determinar las fuerzas, la energía potencial y la presión para 
!algun potencial dado
use precision, only      : pr=>dp
use initialization, only : x, y, z, fx, fy, fz, rcut, ecut, etail, ptail, rho, &
                           & Epot, Pres, Temp, N, L, V
implicit none


contains

subroutine forcelj126
!potencial de Lennard-Jones (12-6)
implicit none
integer                  :: i, j 
real(pr)                 :: rcut2, fr, pres_vir
real(pr)                 :: xi, xj, yi, yj, zi, zj
real(pr)                 :: rx, ry, rz, rij2, r2inv, r6inv

!seteo variables iniciales a 0
Epot = 0._pr
pres_vir = 0._pr
fx(:) = 0._pr;    fy(:) = 0._pr;    fz(:) = 0._pr
rcut2 = rcut*rcut

do i = 1, N-1
    xi = x(i);    yi = y(i);    zi = z(i)
    do j = i+1, N
        xj = x(j);    yj = y(j);    zj = z(j)

        rx = xi - xj
        call minimum_image(rx)
        ry = yi - yj
        call minimum_image(ry)
        rz = zi - zj
        call minimum_image(rz)

        rij2 = rx*rx + ry*ry + rz*rz
 
        if (rij2 <= rcut2) then
            r2inv = 1._pr/rij2
            r6inv = r2inv*r2inv*r2inv
            fr    = 48._pr*r2inv*r6inv*(r6inv - 0.5_pr)

            !ley de newton, acumulo fuerza en i por la particula j y viceversa
            fx(i) = fx(i) + fr*rx
            fx(j) = fx(j) - fr*rx

            fy(i) = fy(i) + fr*ry
            fy(j) = fy(j) - fr*ry

            fz(i) = fz(i) + fr*rz
            fz(j) = fz(j) - fr*rz
            
            !energía potencial del truncado y desplazado (ecut)
            Epot  = Epot + 4._pr*r6inv*(r6inv - 1._pr) - ecut
            pres_vir = pres_vir + rij2*fr
        endif
    enddo
enddo

pres_vir = pres_vir/(3._pr*V)  !(1/dV)*pres_vir
Pres = rho*Temp + pres_vir + Ptail
Epot = Epot + Etail

end subroutine forcelj126


subroutine minimum_image(w)
!donde w= dx, dy ó dz
implicit none
real(pr), intent(inout)    :: w
!w = w - L*anint(w/L)
if (w > 0.5_pr*L) then
    w = w - L
else if (w < -0.5_pr*L) then
    w = w + L
endif
end subroutine minimum_image

end module
