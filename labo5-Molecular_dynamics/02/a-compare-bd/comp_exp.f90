module comp_exp
use precision, only : pr=>dp
use initialization, only : x, y, z, L, N, rho, pi43, g, gr, nbin 
use force, only : minimum_image

contains


subroutine gdr(switch)
! subrutina para calular la distribución radial de pares, siguiendo el algoritmo
! 7 del Frenkel 
! switch = 0 : inicialización; switch = 1 : sampleo; switch = 2 : resultados
!
! para mayor eficiencia switch = 1 debería estar junto al modulo de fuerza
implicit none
integer, intent(in)    :: switch
integer                :: ig, i, j
real(pr)               :: rx, ry, rz, rij2, rij, xi, xj, yi, yj, zi, zj
real(pr)               :: rr, vb, nid, nidinv
integer, save          :: ini = 0, ngr
real(pr), save         :: dg

if (ini == 0) then !pido memoria para gr
    open(125, file='gder.dat', status='replace')
    write(125,'(A)') '# r, g(r)'
    ini = 1
endif

if (switch == 0) then !inicialización

    ngr = 0
    dg  =  0.5_pr*L/nbin !ancho de bin
    gr(:) = 0._pr
    g(:)  = 0

else if (switch == 1) then !sampleo

    ngr = ngr + 1
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
            rij = sqrt(rij2)

            if (rij <= 0.5*L) then
                !    int(rij/dg) = 0 p/ rij < dg ! => sumo 1 porque mi vector
                !                                     va de 1 hasta N
                ig = int(rij/dg) + 1
                g(ig) = g(ig) + 2                ! contribución de part i y j
            endif
        enddo
    enddo

else !ie switch == 2 !resultado

    do i = 1, nbin
        rr = ( real(i,pr) + 0.5_pr )*dg
        vb = ( real((i+1)**3,pr) - real(i**3,pr) )*dg**3 !volumen entre i+1 e i
        nid = pi43*vb*rho !numero de particulas en vb de un gas ideal
        nidinv = 1._pr/nid
        gr(i) = real(g(i),pr)*nidinv/real(ngr*N,pr)
        write(125,'(2(E15.6,x))') rr, gr(i)
    enddo

endif
return
end subroutine gdr


end module
