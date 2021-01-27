program ma
!m(olécula) a(romática)
! particula real en 3d de radio 1 en coordenadas esfericas
!
use precision, only :   pr => dp
use randomnum, only : rand => ran2
implicit none
real(pr), parameter   :: pi = 4._pr * atan(1._pr)
integer               :: semilla, i, j, k, n, M
real(pr)              :: rn1, rn2, dr, theta, phi
real(pr)              :: x, y, z
real(pr), allocatable :: msd(:)

semilla = 684539635                !semilla 
M = 1e2                            !realizaciones del experimento
n = 1000                           !n pasos aleatorios
x = 0; y = 0; z = 0                !la caminata aleatoria empieza en el origen
allocate( msd(n) )
open(18, file='random-walk.xyz', status='replace')
open(19, file='c-msd.dat', status='replace')
!theta entre 0 y pi; phi entre 0 y 2 pi 
dr = .1_pr          !paso del 10 porciento del radio de la particula

do j=1,M
    x = 0; y = 0; z = 0
    do i=1,n
        !rn1 = rand(), rn2 = rand()
        rn1 = rand(semilla); rn2 = rand(semilla)
        !para que estén distribuidos uniformemente
        theta = acos(1._pr - 2._pr*rn1)
        phi   = 2._pr*pi*rn2
    
        x = x + dr*sin(theta)*cos(phi)
        y = y + dr*sin(theta)*sin(phi)
        z = z + dr*cos(theta)

        msd(i) = msd(i) + (real(x*x,pr) + real(y*y,pr)) + real(z*z,pr)

        !trayectoria xyz
        !write(18,*) 1
        !write(18,*) 
        !write(18,'(A3,3(E15.6,2x))') 'H', x, y, z
    enddo
    !write(18,*)
    !write(18,*)
enddo

do k = 1, n
    write(19,'(I5,E15.6)') k, msd(k)/real(M,pr)
enddo

end program
