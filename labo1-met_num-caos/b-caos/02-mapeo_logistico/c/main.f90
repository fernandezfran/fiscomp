program histogram
use precision, only: pr=>dp
implicit none
integer               :: k, i, nbin, ibin 
real(pr)              :: x_0, r, x, f, dx 
integer, allocatable  :: hist(:)

x_0 = 0.6_pr
r = 4._pr
nbin = 100
allocate( hist(nbin) )
dx = 1._pr/real(nbin,pr)
hist(:) = 0

x = x_0
do k = 1, 300 !primeros trecientos puntos que descarto
    call logistic_map(r, x, f)
    x = f
enddo

do k = 301, 10300
    call logistic_map(r, x, f)
    x = f
    ibin = int(x/dx) + 1
    hist(ibin) = hist(ibin) + 1
enddo

open(13, file='histo.dat', status='replace')
do i = 1, nbin
    write(13,*) (real(i-1) + 0.5_pr)*dx, hist(i)/real(nbin,pr)
enddo
close(13)
contains

subroutine logistic_map(mu, xx, fx)
implicit none
real(pr), intent(in)  :: mu, xx
real(pr), intent(out) :: fx

    fx = mu * xx * (1._pr - xx)

    return
end subroutine logistic_map

end program
