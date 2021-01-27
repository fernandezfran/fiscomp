program sis_caos
use precision, only: pr=>dp
implicit none
integer            :: i, j, k
real(pr)           :: r, x
character(len=24)   :: aux_r, aux_x
character(len=30)  :: filename

!valores de x_0 y r de la guía
!x_0 = (/0.06_pr, 0.3_pr, 0.6_pr, 0.9_pr/)
!r = (/1.5_pr, 3.3_pr, 3.5_pr, 3.55_pr, 4._pr/)

write(*,*) 'factor de crecimiento?'
read(*,*) r
write(*,*) 'población inicial?'
read(*,*) x

write(aux_r,'(F4.2)') r
write(aux_x,'(F4.2)') x

filename = 'xi_'//trim(adjustl(aux_x))//'-r_'//trim(adjustl(aux_r))//'.dat'
write(*,*) 'ver datos en el archvio de salida ', filename

open(23, file=filename, status='replace')
       
do k = 0, 500 !desde cero para que escriba la cond inicial
    write(23,'(I3,x,E15.6)') k, x
    call logistic_map(r, x)
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
