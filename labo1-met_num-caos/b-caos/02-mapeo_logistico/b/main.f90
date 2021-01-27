program pow_spe_caos
use, intrinsic     :: iso_c_binding
use precision, only: pr=>dp
implicit none
include "/usr/include/fftw3.f03"
real(pr), parameter                    :: pi = 4._pr * atan(1._pr)
real(C_double), allocatable            :: fx(:)
complex(C_double_complex), allocatable :: tfx(:)
integer                                :: N, j, k, i
real(pr)                               :: r(5), x_0, x, rango, factor
type(C_ptr)                            :: plan_rc
character(len=24)                      :: aux_r, aux_x
character(len=30)                      :: filename, formato


! datos de entrada
x_0 = 0.6_pr
r   = (/1.5_pr, 3.3_pr, 3.5_pr, 3.55_pr, 4._pr/)
N   = 10000 !ptos para la t. de f.

rango  = real(N,pr)  
factor = 2._pr * pi / rango !freq de nysquit
formato = '(E13.4,x,2(E21.6,x))'

allocate( fx(N), tfx(N/2 + 1) )
plan_rc = fftw_plan_dft_r2c_1d(N, fx, tfx, FFTW_MEASURE)

do j = 1, 5
    x = x_0
    write(aux_r,'(F4.2)') r(j)
    write(aux_x,'(F4.2)') x   
    filename = 'xi_'//trim(adjustl(aux_x))//'-r_'//trim(adjustl(aux_r))//'.dat'
    
    do k = 0, 299 !puntos descartados para el espectro de potencias
        call logistic_map(r(j), x)
    enddo
    
    fx(:) = 0._pr
    do k = 300, 300 + N 
        call logistic_map(r(j), x)
        fx(k - 299) = x
    enddo

    tfx(:) = 0._pr
    call fftw_execute_dft_r2c(plan_rc, fx, tfx)
    
    open(23, file=filename, status='replace')
    do i = 1, N/2 + 1
        write(23,formato) factor*real(i-1,pr), real(tfx(i))/real(N,pr), &
                          & aimag(tfx(i))/real(N,pr)
    enddo
    close(23)


enddo

deallocate( fx, tfx )

contains

subroutine logistic_map(mu, xx)
implicit none
real(pr), intent(in)    :: mu
real(pr), intent(inout) :: xx

    xx = mu * xx * (1._pr - xx)

    return
end subroutine logistic_map

end program
