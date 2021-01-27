program cuad2d
!
! caminata aleatoria en una grilla bidimensional de paso unitario
!
use precision, only :   pr => dp
use randomnum, only : rand => ran2
implicit none
integer              :: M, n, i, j, k, semilla
integer              :: x, y
real(pr)             :: rn
integer, allocatable :: cuadrante(:,:)

semilla = 157194652                !semilla 
M = 1e5                            !realizaciones del experimento
n = 100                            !n pasos aleatorios
x = 0; y = 0                       !la caminata aleatoria empieza en el origen

allocate( cuadrante(n,4) )

open(20, file='b-cuadrantes.dat', status='replace')

cuadrante(:,:) = 0
do j = 1, M
    x = 0; y = 0
    do i = 1, n
        rn = rand(semilla)
        !rn = rand()
        if     (0._pr   <= rn .and. rn < 0.25_pr) then
            x = x + 1._pr
        elseif (0.25_pr <= rn .and. rn < 0.5_pr)  then
            x = x - 1._pr
        elseif (0.5_pr  <= rn .and. rn < 0.75_pr) then
            y = y + 1._pr
        else
            y = y - 1._pr
        endif

        if (x >= 0 .and. y > 0) then
            cuadrante(i,1) = cuadrante(i,1) + 1
        elseif (x < 0 .and. y >= 0) then
            cuadrante(i,2) = cuadrante(i,2) + 1
        elseif (x <= 0 .and. y < 0) then
            cuadrante(i,3) = cuadrante(i,3) + 1
        elseif (x > 0 .and. y <= 0) then
            cuadrante(i,4) = cuadrante(i,4) + 1
        endif

    enddo
enddo

!no divido por cuatro en este promedio porque podría no darme un integer
cuadrante(:,1) = cuadrante(:,1) + cuadrante(:,2) + cuadrante(:,3) + cuadrante(:,4)

do k = 1, n
    write(20,'(I5,2x,E15.6)') k, 0.25_pr*real(cuadrante(k,1),pr)/real(M,pr)
enddo

end program
