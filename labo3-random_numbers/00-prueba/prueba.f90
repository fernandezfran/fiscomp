program distribuciones
!
! prueba de como generar los numeros aleatorios
!
use precision, only: pr=>dp
use randomnum
implicit none
integer            :: i, semilla
real(pr)           :: rand


!ran0 --> funciona
semilla = 12345
do i = 1, 100
    rand = ran0(semilla)

    write(*,*) rand
enddo
write(*,*)


!ran2 --> funciona
semilla = 12345
do i = 1, 100
    rand = ran2(semilla)

    write(*,*) rand
enddo
write(*,*)

!mzran --> funciona
do i = 1, 100
    rand = rmzran()

    write(*,*) rand
enddo
write(*,*)


!MT --> funciona
do i = 1, 100
    rand = grnd()

    write(*,*) rand
enddo
write(*,*)


end program
