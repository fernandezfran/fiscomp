program Lyapunov
!
!programa para obtener el exponente de Lyapunov para r entre 2 y 4, con
!1500 puntos intermedios
!
use precision, only: pr=>dp
implicit none
integer            :: N, M, i, k
real(pr)           :: x, ri, rf, r, exp_ly

open(23, file='lyapunov_vs_r.dat', status='replace')

ri = 2._pr
rf = 4._pr
N  = 1500 !cantidad de r entre ri y rf
M  = 1500

do i = 1, N
    r = ri + real(i,pr)*(rf - ri)/real(N,pr) !para tener N pts entre ri y rf
    x = 0.6_pr
    exp_ly = 0._pr
    do k = 1, M
        call logistic_map(r, x)
        if (k > 300) then
            exp_ly = exp_ly + log(abs(r * (1._pr - 2._pr*x)))
        endif
    enddo
    write(23,*) r, exp_ly/real(M-300,pr)
enddo



contains

subroutine logistic_map(mu, xx)
implicit none
real(pr), intent(in)    :: mu
real(pr), intent(inout) :: xx

    xx = mu * xx * (1._pr - xx)

    return
end subroutine logistic_map

end program
