program msd2d
!
! caminata aleatoria en una grilla bidimensional de paso unitario
!
use precision, only :   pr => dp
use randomnum, only : rand => grnd
implicit none
integer               :: M, n, i, j, k, semilla
integer               :: x, y
real(pr)              :: rn
real(pr), allocatable :: msd(:)

semilla = 735612575                !semilla 
M = 1e5                           !realizaciones del experimento
n = 1000                           !n pasos aleatorios
x = 0; y = 0                       !la caminata aleatoria empieza en el origen

allocate( msd(n) )

open(20, file='a-msd.dat', status='replace')
!open(21, file='10-tray.xy', status='replace')

msd(:) = 0._pr
do j = 1, M
    x = 0; y = 0
    do i = 1, n
        
        !write(21,'(2(I6,2x))') x, y

        rn = rand()
        if     (0._pr   <= rn .and. rn < 0.25_pr) then
            x = x + 1._pr
        elseif (0.25_pr <= rn .and. rn < 0.5_pr)  then
            x = x - 1._pr
        elseif (0.5_pr  <= rn .and. rn < 0.75_pr) then
            y = y + 1._pr
        else
            y = y - 1._pr
        endif

        msd(i) = msd(i) + (real(x*x,pr) + real(y*y,pr))

    enddo

    !write(21,*)
    !write(21,*)

enddo

do k = 1, n
    write(20,'(I5,E15.6)') k, msd(k)/real(M,pr)
enddo


end program
