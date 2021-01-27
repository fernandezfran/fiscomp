program pru
use precision, only :   pr=>dp
use randomnum, only : rand=>grnd
implicit none
integer               :: ibin, nbin, i, j, M
real(pr)              :: p, rn
real(pr), allocatable :: hist(:)

open(9,file='prueba.dat',status='replace')

nbin = 100
allocate( hist(nbin) )
hist(:) = 0._pr
M = 100000

do i = 1, M
    rn = rand()
    p = ldp(rn)
    ibin = nint(p*100)
    hist(ibin) = hist(ibin) + 1 
enddo

do j = 1, nbin
    write(9,*) real(j,pr)-0.5_pr, hist(j)/real(nbin,pr)
enddo

contains

function ldp(u)
!
!distribución de nros aleatorios como ley de potencia, p(x) = (k+1)*x^k
!
implicit none
real(pr), intent(in) :: u
integer              :: k
real(pr)             :: power, ldp

k     = 2
power = 1._pr/(real(k+1,pr))

ldp = u**power

end function ldp


end program
