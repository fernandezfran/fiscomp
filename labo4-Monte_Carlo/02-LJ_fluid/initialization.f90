module initialization
use precision, only : pr=>dp
implicit none
real(pr), allocatable :: x(:), y(:), z(:)
real(pr)              :: T, L, V, rho, rcut, utail
integer               :: N

contains


subroutine xyzsc()
implicit none
integer               :: i, j, k, idxp
real(pr)              :: dist

!indice de la partícula
idxp = 1
do i = 1, int(L)
    do j = 1, int(L)
        do k = 1, int(L)
            x(idxp) = real(i-1,pr)
            y(idxp) = real(j-1,pr)
            z(idxp) = real(k-1,pr)
            idxp = idxp + 1
        enddo
    enddo
enddo
!multiplico por la densidad^(1/3) para reescalar las distancias interatómicas
dist = 1._pr/(rho**(1._pr/3._pr))
x(:) = dist*x(:)
y(:) = dist*y(:)
z(:) = dist*z(:)

end subroutine xyzsc


subroutine xyzfcc()
!REVISAR
implicit none
integer                :: i, j, k, s, idxp
real(pr)               :: dist, ss
real(pr), dimension(2) :: shift

shift(1) = 0._pr
shift(2) = 0.5_pr
s = 0 
ss = shift(1)
idxp = 1
do i = 1, int(L)
    do j = 1, int(L)
        do k = 1, int(L)
            x(idxp) = real(i-1,pr) + ss
            y(idxp) = real(j-1,pr) + ss
            z(idxp) = real(k-1,pr) + ss
            idxp = idxp + 1
        enddo
        s = s + 1
        ss = shift(2)
        if (mod(s,2) == 0) ss = shift(1)
    enddo
enddo
!dist = 1._pr/(rho**(1._pr/3._pr))
!x(:) = dist*x(:)
!y(:) = dist*y(:)
!z(:) = dist*z(:)

end subroutine xyzfcc

subroutine xyzrand()
use randomnum, only : rand => rmzran
implicit none
integer             :: i, iend
real(pr)            :: dist

do i = 1, int(L)
    x(i) = rand()*L
    y(i) = rand()*L
    z(i) = rand()*L
enddo
dist = rho**(1._pr/3._pr)
x(:) = dist*x(:)
y(:) = dist*y(:)
z(:) = dist*z(:)

end subroutine xyzrand


subroutine pbc(cordi)
!trae una coordenada al intervalor [0, L)
implicit none
real(pr), intent(inout) :: cordi

if (cordi < 0) then
    cordi = cordi + L
else if (cordi >= L) then
    cordi = cordi - L
endif

end subroutine pbc


end module
