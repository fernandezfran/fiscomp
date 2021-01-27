module initial_cond 
use precision, only   : pr => dp
real(pr)             :: T               !temp
integer, allocatable :: spin(:)         !spines
integer              :: L, N            !largo de la grilla

contains

subroutine order()
!condición inicial todos up
implicit none

spin(:) = 1

end subroutine order


subroutine mess()
!condición inicial aleatoria
use randomnum, only : rand => rmzran
implicit none
integer              :: i

do i = 1, size(spin)
    if (rand() >= 0.5_pr) then
        spin(i) = 1
    else
        spin(i) = -1
    endif
enddo

end subroutine mess


end module initial_cond
