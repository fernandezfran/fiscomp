program histo
use precision, only: pr=>dp
implicit none
integer            :: nbin, ibin, i, j
real(pr)           :: x, r, dx
integer, dimension(:), allocatable :: hist

nbin = 80
x = 0.6_pr
r = 4._pr
allocate( hist(nbin) )
dx = 1._pr/real(nbin, pr)
hist(:) = 0

open(22, file='xt.dat', status='replace')
do i = 1,10300
    x = r * x * (1._pr - x)
    write(22,*) i, x
    if (i > 300) then
        ibin = int(x/dx) + 1 !divido y me quedo con la parte entera red. para abajo
        hist(ibin) = hist(ibin) + 1
    endif
enddo
close(22)

open(21, file='hist.dat', status='replace')
do i = 1, nbin !el punto x en medio del bin
    write(21,*) real(i-1,pr)*dx + 0.5_pr*dx, hist(i)
enddo




end program
