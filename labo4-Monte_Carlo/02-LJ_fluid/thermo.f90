module thermo
!calculo con pontencial de LJ 12-6 con truncado simple
use precision, only      : pr=>dp
use initialization, only : x, y, z, rcut, L, N, rho, T

contains


subroutine deltau(idx,xnew,ynew,znew,dU)
implicit none
integer, intent(in)      :: idx
real(pr), intent(in)     :: xnew, ynew, znew
real(pr), intent(out)    :: dU
real(pr)                 :: rx, ry, rz, rij2, uljnew, uljold, rm2, rm6, rm12
integer                  :: j

uljnew = 0._pr; uljold = 0._pr

!recorro todas las partículas centrandome en la particula idx
!tanto desplazada como anterior
do j = 1, N
    if (j .ne. idx) then
        !### DESPLAZADA
        !calculo la distancia mínima (imagen) de la partícula j
        rx = xnew - x(j)
        call minimum_image(rx)
        ry = ynew - y(j)
        call minimum_image(ry)
        rz = znew - z(j)
        call minimum_image(rz)
        
        rij2 = rx*rx + ry*ry + rz*rz

        !acumulo la energía de LJ:
        if (rij2 <= rcut*rcut) then
            rm2  = 1._pr/rij2
            rm6  = rm2 * rm2 * rm2
            rm12 = rm6 * rm6
            uljnew = uljnew + 4._pr*(rm12 - rm6)
        endif
        
        !### ANTERIOR
        !calculo la distancia mínima (imagen) de la partícula j
        rx = x(idx) - x(j)
        call minimum_image(rx)
        ry = y(idx) - y(j)
        call minimum_image(ry)
        rz = z(idx) - z(j)
        call minimum_image(rz)
        
        rij2 = rx*rx + ry*ry + rz*rz

        !acumulo la energía de LJ:
        if (rij2 <= rcut*rcut) then
            rm2  = 1._pr/rij2
            rm6  = rm2 * rm2 * rm2
            rm12 = rm6 * rm6
            uljold = uljold + 4._pr*(rm12 - rm6)
        endif
    endif
enddo

dU = uljnew - uljold

end subroutine deltau


subroutine calculate(utail,ptail,U,P)
implicit none
real(pr), intent(in)  :: utail, ptail
real(pr), intent(out) :: U, P
real(pr)              :: rx, ry, rz, rij2, ulj, rm2, rm6, rm12, pij
integer               :: i, j

U = 0._pr
P = 0._pr
do i = 1, N-1
    do j = i+1, N
        !calculo la distancia mínima (imagen) de la partícula i a la j y hago PBC
        rx = x(i) - x(j)
        call minimum_image(rx)
        ry = y(i) - y(j)
        call minimum_image(ry)
        rz = z(i) - z(j)
        call minimum_image(rz)

        rij2 = rx*rx + ry*ry + rz*rz

        if (rij2 <= rcut*rcut) then
            rm2  = 1._pr/rij2
            rm6  = rm2 * rm2 * rm2
            rm12 = rm6 * rm6
            ulj = 4._pr*(rm12 - rm6)
            U   = U + ulj

            !pij = -ulj
            !P   = P + pij
        endif

    enddo
enddo

U = U - utail
!P = P + rho*T - ptail

end subroutine calculate


subroutine minimum_image(r)
implicit none
real(pr), intent(inout)        :: r
if (r > 0.5_pr*L) then
    r = r - L 
else if (r < -0.5_pr*L) then
    r = r + L
endif
end subroutine minimum_image

end module
