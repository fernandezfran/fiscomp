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

!distancia interatómica
dist = rho**(1._pr/3._pr)
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
x(:) = dist*x(:)
y(:) = dist*y(:)
z(:) = dist*z(:)


end subroutine xyzsc


subroutine xyzrand()
use randomnum, only : rand => rmzran
implicit none
integer             :: i, iend
real(pr)            :: dist

dist = rho**(1._pr/3._pr)
iend = int(L) - 1

do i = 0, iend
    x(i) = rand()*L
    y(i) = rand()*L
    z(i) = rand()*L
enddo
x(:) = dist*x(:)
y(:) = dist*y(:)
z(:) = dist*z(:)


end subroutine xyzrand

end module
