module initialization
use precision, only : pr=>dp
implicit none
!###DEFINICIÓN DE VARIABLES GLOBALES
!parámetros dinámicos de la simulación
real(pr), allocatable :: x(:), y(:), z(:), vx(:), vy(:), vz(:), fx(:), fy(:), fz(:)
!valores iniciales
real(pr)              :: T0, L, V, rho, dt, rcut, ecut, Etail, Ptail
integer               :: N
!variables a medir
real(pr)              :: Et, Epot, Ekin, Pres, Temp

contains


subroutine initposfcc
!posiciones iniciales en una caja cúbica de lado L, para un material de
!densidad \rho con una estructura cristalina FCC
implicit none
integer                :: i, j, k, idxp, nucell
real(pr)               :: dist

nucell = nint(real(N/4,pr)**(1._pr/3._pr))
idxp = 1
do i = 1, nucell            !number of unit cells
    do j = 1, nucell
        do k = 1, nucell
            x(idxp+0) = real(i-1,pr)
            y(idxp+0) = real(j-1,pr)
            z(idxp+0) = real(k-1,pr)

            x(idxp+1) = real(i-1,pr) + 0.5_pr
            y(idxp+1) = real(j-1,pr) + 0.5_pr
            z(idxp+1) = real(k-1,pr)

            x(idxp+2) = real(i-1,pr) + 0.5_pr
            y(idxp+2) = real(j-1,pr) 
            z(idxp+2) = real(k-1,pr) + 0.5_pr

            x(idxp+3) = real(i-1,pr) 
            y(idxp+3) = real(j-1,pr) + 0.5_pr
            z(idxp+3) = real(k-1,pr) + 0.5_pr

            idxp = idxp + 4
        enddo
    enddo
enddo
dist = (4._pr/rho)**(1._pr/3._pr)
x(:) = dist*x(:)
y(:) = dist*y(:)
z(:) = dist*z(:)

end subroutine initposfcc


subroutine initvelrand
!asigna velocidades aleatorias a las partículas, con una distribucion gaussiana
implicit none
integer             :: i, fro
real(pr)            :: sumvx, sumvy, sumvz, sumv2, sf

fro = 9326953

sumvx = 0._pr;    sumvy = 0._pr;    sumvz = 0._pr
sumv2 = 0._pr

do i = 1, N
    vx(i) = gasdev(fro)
    vy(i) = gasdev(fro)
    vz(i) = gasdev(fro)

    sumvx = sumvx + vx(i)
    sumvy = sumvy + vy(i)
    sumvz = sumvz + vz(i)

    sumv2 = sumv2 + vx(i)*vx(i) + vy(i)*vy(i) + vz(i)*vz(i)
enddo

sumvx = sumvx/real(N,pr)
sumvy = sumvy/real(N,pr)
sumvz = sumvz/real(N,pr)
Temp  = sumv2/real(3*N,pr)

sf = sqrt(T0/Temp)        !factor de escaleo
do i = 1, N
    !extraigo velocidad del centro de masa y escaleo
    vx(i) = sf*(vx(i) - sumvx)
    vy(i) = sf*(vy(i) - sumvy)
    vz(i) = sf*(vz(i) - sumvz)
enddo

Ekin  = 0.5_pr*sumv2

end subroutine initvelrand


function gasdev(semilla)
! método de box-muller
!    se considera x0 = 0 y sigma = 1
use randomnum, only : rand => ran2
implicit none
integer, intent(inout) :: semilla
real(pr)               :: gasdev, gset
real(pr)               :: u1, u2, fs, ang, pi = -4._pr*atan(1._pr)
integer                :: iset = 0
save iset, gset

if (iset == 0) then
    u1  = rand(semilla); u2  = rand(semilla)
    fs  = sqrt(-2._pr*log(u1))
    ang = 2._pr*pi*u2

    gset     = fs*cos(ang)
    gasdev = fs*sin(ang)

    iset = 1
else
    gasdev = gset
    iset = 0
endif

endfunction gasdev


subroutine pbc(cordi)
!trae una coordenada al intervalo [0, L)
implicit none
real(pr), intent(inout) :: cordi

!cordi = cordi - L*anint(cordi/L) 
if (cordi < 0._pr) then
    cordi = cordi + L
else if (cordi >= L) then
    cordi = cordi - L
endif

end subroutine pbc


end module
