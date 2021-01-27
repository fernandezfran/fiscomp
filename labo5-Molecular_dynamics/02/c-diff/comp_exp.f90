module comp_exp
use precision, only : pr=>dp

contains


subroutine gdr(switch)
! subrutina para calular la distribución radial de pares, siguiendo el algoritmo
! 7 del Frenkel 
! switch = 0 : inicialización; switch = 1 : sampleo; switch = 2 : resultados
!
! para mayor eficiencia switch = 1 debería estar junto al modulo de fuerza
use initialization, only : x, y, z, L, N, rho, pi43, g, gr, nbin 
use force, only : minimum_image
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
                ig = nint(rij/dg)
                g(ig) = g(ig) + 2 !contribución de part i y j
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


subroutine diffusion(switch,nsamp)
! subrutina para calcular el desplazamiento cuadrático medio, que puede ser 
! directamente relacionado con el coeficiendo de difusion (autodifusion en LJ).
! Se sigue el algoritmo 8 del Frenkel y las notas (pseudo-codigo) de A. Banchio
!
! switch = 0 : inicialización ; switch = 1 : sampleo ; switch = 2 : resultados
! nsamp  = cada cuando samplea (se utiliza solo en 0 para definir dtime)
use initialization, only : x, y, z, ix, iy, iz, r2t, x0, y0, z0, dt, N, L
implicit none
integer, intent(in)  :: switch, nsamp
integer              :: i, ilast, inext, j
real(pr)             :: time
integer, save        :: init = 0, ntel
real(pr), save       :: dtime

if (init == 0) then
    open(188, file='msdvst.dat', status='replace')
    write(188,'(A)') '# t, msd'
    init = 1
endif

if (switch == 0) then !inicialización

    ntel = 0
    r2t = 0._pr
    x0(:) = x(:);    y0(:) = y(:);    z0(:) = z(:)

else if (switch == 1) then ! sampleo
    
    ntel = ntel + 1
    do i = 1, N
        r2t = (x(i) + real(ix(i),pr)*L - x0(i))**2 + &
            & (y(i) + real(iy(i),pr)*L - y0(i))**2 + &
            & (z(i) + real(iz(i),pr)*L - z0(i))**2
    enddo
    
    r2t = r2t/real(N,pr)
    write(188,'(2(E15.6,x))') ntel*dt, r2t

endif
return
end subroutine diffusion

end module
