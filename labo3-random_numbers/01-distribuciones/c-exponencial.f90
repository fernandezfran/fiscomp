program ldp
!l(ey) d(e) p(otencias)
!distribución no uniforme, a través del método de la función inversa
!
!si quiero la distribución f(x) = (1/lambda)*exp(-x/lambda), para x >= 0,
! entonces, x tiene que tomar valores con la forma -lambda*ln(u)
!
! $ gfortran c-exponencial.f90 *.o && ./a.out
!    con M=5000 ya se mejora considerablemente el grafico.
!
use precision, only : pr=>dp
use randomnum
real(pr), dimension(4)  :: rand
integer, dimension(2)   :: semilla !una semilla por método que la requiere
integer                 :: i, j, k, M, ibin, nbin
real(pr), allocatable   :: hist(:)
real(pr)                :: f

open(15, file='c-histo.dat', status='replace')

nbin    = 15                        !cantidad de intervalos del histograma
allocate( hist(nbin) )
hist(:) = 0
M       = 1000                      !cantidad de nros. aleatorios por met.

semilla = (/ 157234, 259623 /)
!este vector no le gusta mucho... no modifica la seed
rand    = (/ ran0(semilla(1)), ran2(semilla(2)), rmzran(), grnd() /)
do i = 1,4 

    do j = 1,M
        f    = finv(rand(i))
        ibin = int(f)
        if (ibin .le. nbin) hist(ibin) = hist(ibin) + 1
        !se soluciona con esto que no es para nada eficiente
        rand = (/ ran0(semilla(1)), ran2(semilla(2)), rmzran(), grnd() /)
    enddo

    do k = 1, nbin
        write(15,*) (real(k,pr) + 0.5_pr), hist(k)/M
    enddo

    write(15,*)
    write(15,*)
    hist(:) = 0

enddo


contains


function finv(u)
!función inversa para obtener distribución no uniform (1/lambda)*exp(-x/lambda)
implicit none
real(pr), intent(in) :: u
real(pr), parameter  :: lambda = 2._pr
real(pr)             :: finv

finv = - lambda*log(u)

endfunction


end program
