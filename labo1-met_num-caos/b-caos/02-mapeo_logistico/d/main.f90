program diagrama_orbitas
!
!programa para estudiar el diagrama de orbitas en el plano x, r.
!Los datos de entrada son los r inicial y final, los puntos N que se desea que haya
!entre ellos y los M pasos temporales (de los cuales los primeros 300 se descartan).
!
use precision, only: pr=>dp
implicit none
integer            :: N, M, i, k
real(pr)           :: x, ri, rf, r

open(23, file='x_vs_r.dat', status='replace')

write(*,*) 'factor de crecimiento inicial?'
read(*,*) ri
write(*,*) 'factor de crecimiento final?'
read(*,*) rf
write(*,*) 'cantidad de puntos entre ellos?'
read(*,*) N
write(*,*) 'cantidad de pasos temporales? (notar que los primeros 300 se descartan)'
read(*,*) M

!ri = 3.4_pr
!rf = 4._pr
!N  = 1500 !cantidad de r entre ri y rf

do i = 1, N
    x = 0.6_pr
    r = ri + real(i,pr)*(rf - ri)/real(N,pr) !para tener 1000 pts entre ri y rf
    do k = 1, M
        call logistic_map(r, x)
        if (k > 300) then
            write(23,*) r, x
        endif
    enddo
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
