program sis_caos
use precision, only: pr=>dp
implicit none
integer            :: i, j, k
real(pr)           :: x_0(4), r(5) 
real(pr)           :: x, f

x_0 = (/0.06_pr, 0.3_pr, 0.6_pr, 0.9_pr/)
r = (/1.5_pr, 3.3_pr, 3.5_pr, 3.55_pr, 4._pr/)

open(23, file='trayectorias.dat', status='replace')

do i = 1, 4
    do j = 1, 5
        write(23,'(A)') '# para r y x_0 iguales a'
        write(23,'(A, 2(E12.3,2x))') '#', x_0(i), r(j)
        x = x_0(i)
        do k = 0, 500
            write(23,'(I3,x,E12.3)') k, x
            call logistic_map(r(j), x, f)
            x = f
        enddo
    enddo
enddo

contains

subroutine logistic_map(mu, xx, fx)
implicit none
real(pr), intent(in)  :: mu, xx
real(pr), intent(out) :: fx

    fx = mu * xx * (1._pr - xx)

    return
end subroutine logistic_map

end program
