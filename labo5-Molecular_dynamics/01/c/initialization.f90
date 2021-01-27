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
!asigna velocidades aleatorias a las partículas
use randomnum, only : ranf => rmzran
implicit none
integer             :: i
real(pr)            :: sumvx, sumvy, sumvz, sumv2, sf

sumvx = 0._pr;    sumvy = 0._pr;    sumvz = 0._pr
sumv2 = 0._pr

do i = 1, N
    vx(i) = ranf() - 0.5_pr
    vy(i) = ranf() - 0.5_pr
    vz(i) = ranf() - 0.5_pr

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
