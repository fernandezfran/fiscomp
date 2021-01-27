module comp_exp
use precision, only : pr=>dp

contains


subroutine gdr(switch)
! subrutina para calular la distribución radial de pares, siguiendo el algoritmo
! 7 del Frenkel 
! switch = 0 : inicialización; switch = 1 : sampleo; switch = 2 : resultados
!
! para mayor eficiencia switch = 1 debería estar junto al modulo de fuerza
use initialization, only : x, y, z, L, N, rho, pi43 
use force, only : minimum_image
implicit none
integer, intent(in)    :: switch
integer                :: ig, i, j
real(pr)               :: rx, ry, rz, rij2, rij, xi, xj, yi, yj, zi, zj
real(pr)               :: rr, vb, nid, nidinv
integer, save          :: ini = 0, nbin, ngr
real(pr), save         :: dg
real(pr), allocatable, save :: gr(:)
integer, allocatable, save  :: g(:)

if (ini == 0) then !pido memoria para gr
    open(125, file='gder.dat', status='replace')
    write(125,'(A)') '# r, g(r)'
    ini = 1
    nbin = 1000
    allocate( g(nbin), gr(nbin) )
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
use initialization, only : x, y, z, ix, iy, iz, dt, N, L
implicit none
integer, intent(in)  :: switch, nsamp
integer              :: i, ilast, inext, j
real(pr)             :: time
integer, save        :: init = 0, ntel, Nmsd
real(pr), save       :: dtime
integer, allocatable, save  :: ntime(:)
real(pr), allocatable, save :: x0(:,:), y0(:,:), z0(:,:), r2t(:)

if (init == 0) then
    open(188, file='msdvst.dat', status='replace')
    write(188,'(A)') '# t, msd'
    init = 1
    Nmsd = 2000
    allocate( ntime(Nmsd), r2t(Nmsd) )
    allocate( x0(N,Nmsd), y0(N,Nmsd), z0(N,Nmsd) )
endif

if (switch == 0) then !inicialización

    ntel = 0
    dtime = dt*nsamp
    ntime(:) = 0
    r2t(:) = 0._pr
    x0(:,:) = 0._pr;    y0(:,:) = 0._pr;    z0(:,:) = 0._pr

else if (switch == 1) then ! sampleo

    ntel = ntel + 1
    ilast = mod(ntel - 1, Nmsd) + 1

    do i = 1, N
        x0(i,ilast) = x(i) + real(ix(i),pr)*L
        y0(i,ilast) = y(i) + real(iy(i),pr)*L
        z0(i,ilast) = z(i) + real(iz(i),pr)*L
    enddo
    
    if (ntel > Nmsd) then 
        do j = 1, Nmsd
            inext = mod(ntel - j, Nmsd) + 1
            do i = 1, N
                r2t(j) = r2t(j) + (x0(i,ilast) - x0(i,inext))**2 &
                              & + (y0(i,ilast) - y0(i,inext))**2 &
                              & + (z0(i,ilast) - z0(i,inext))**2
            enddo
            ntime(j) = ntime(j) + 1
        enddo
    endif

else !ie switch == 2 !resultado

    do j = 1, Nmsd
        time = real(j-1,pr)*dtime
        r2t(j) = r2t(j)/real(N*ntime(j),pr)
        write(188,'(2(E15.6,x))') time, r2t(j)
    enddo

endif
return
end subroutine diffusion


subroutine veldist(nbin, switch)
use initialization, only : vx, vy, vz, N
implicit none
integer, intent(in)  :: nbin, switch
integer              :: ibin, j
real(pr)             :: histj, histden
integer, save        :: frames
real(pr), save       :: dx
integer, allocatable, save :: hist(:)

if (switch == 0) then

    dx = 10._pr/nbin ! (xman - xmin)/nbin
    allocate( hist(nbin) )
    hist(:) = 0
    frames = 0
    open(163, file='histo.dat', status='replace')

else if (switch == 1) then

    frames = frames + 1
    do j = 1, N
        ibin = nint(vx(j)/dx)
        if ((ibin + nbin/2 > 0) .and. (ibin + nbin/2 <= nbin)) then
            hist(ibin + nbin/2) = hist(ibin + nbin/2) + 1
        endif
        ibin = nint(vy(j)/dx)
        if ((ibin + nbin/2 > 0) .and. (ibin + nbin/2 <= nbin)) then
            hist(ibin + nbin/2) = hist(ibin + nbin/2) + 1
        endif
        ibin = nint(vz(j)/dx)
        if ((ibin + nbin/2 > 0) .and. (ibin + nbin/2 <= nbin)) then
            hist(ibin + nbin/2) = hist(ibin + nbin/2) + 1
        endif
    enddo

else if (switch == 2) then
    histden = real(3*N,pr)*real(frames,pr)*dx
    do j = 1, nbin
        histj = real(hist(j),pr)/histden
        write(163, '(2(E15.6,x))') real(j-nbin/2,pr)*dx, histj
    enddo

endif


end subroutine veldist


end module
