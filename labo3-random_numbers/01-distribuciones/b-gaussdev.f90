program mdbm
!m(étodo) d(e) b(ox) m(üller)
    !PARA IMPLEMENTAR
    !f = (x-x0)/sigma 
    !  => si quiero un x0 y sigma en particular obtengo los nros aleatorios haciendo
    !        x = sigma*f + x0
    !
!distribución no uniforme, a través del método de Box-Müller
!
! quiero la distribución no uniforme gaussiana
!     f(x) = (1/sqrt(2*pi*sigma^2)) * exp(-(x - x0)^2/2*sigma^2)
!
! hay que ir cambiando a mano la forma de generar numeros aleatorios
!  en la funcion gaussdev de abajo
!    solo funciona cuando M/nbin = 100 y nbin >= 1000, no sé porqué...
!
! $ gfortran b-gaussdev.f90 *.o && ./a.out
!
use precision, only : pr=>dp
use randomnum
real(pr), parameter     :: pi = 4._pr*atan(1._pr)
integer                 :: semilla, j, k, M, ibin, nbin
real(pr), allocatable   :: hist(:)
real(pr)                :: f, dx

open(15, file='b-histo.dat', status='replace')

nbin    = 1000                         !cantidad de intervalos del histograma
dx      = 1._pr/real(nbin,pr)
allocate( hist(nbin) )
hist(:) = 0
M       = 100000                       !cantidad de nros. aleatorios por met.

semilla = 13746964
do j = 1,M
    f    = gaussdev(semilla)
    ibin = nint(f*real(M,pr)/real(nbin,pr))
    if (abs(ibin) .le. nbin) hist(ibin+nbin/2) = hist(ibin+nbin/2) + 1
enddo

do k = 1, nbin
    write(15,*) (real(k,pr) - 0.5_pr*nbin)*(real(nbin,pr)/real(M,pr)), hist(k)/real(nbin,pr)
enddo


contains


function gaussdev(fro)
! método de box-muller
!    cambiar a mano el generador de num aleatorios en las variables u1 y u2
!
!    se considera x0 = 0 y sigma = 1
implicit none
integer, intent(inout) :: fro !=semilla en noruego
real(pr)               :: gaussdev, gset
real(pr)               :: u1, u2, fs, ang
integer                :: iset = 0
save iset, gset

if (iset == 0) then
    !u1  = ran0(fro); u2  = ran0(fro)
    !u1  = ran2(fro); u2  = ran2(fro)
    !u1  = rmzran(); u2  = rmzran()
    u1  = grnd(); u2  = grnd()
    fs  = sqrt(-2._pr*log(u1))
    ang = 2._pr*pi*u2

    gset     = fs*cos(ang)
    gaussdev = fs*sin(ang)

    iset = 1

else

    gaussdev = gset
    iset = 0

endif

endfunction


end program
