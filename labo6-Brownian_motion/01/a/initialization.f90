module initialization
use precision, only : pr=>dp
implicit none
!###DEFINICIÓN DE VARIABLES GLOBALES
!parámetros dinámicos de la simulación
real(pr), allocatable :: x(:), y(:), z(:), fx(:), fy(:), fz(:)
!valores iniciales
real(pr)              :: T0, L, V, rho, eta, D0,dt, rcut, ecut, Etail, Ptail
integer               :: N, semilla
!variables a medir
real(pr), parameter   :: pi = 4._pr*atan(1._pr), pi43 = (4._pr*pi)/3._pr, &
                         pi3 = 3._pr*pi
real(pr)              :: Epot, Pres
!variables de g(r)
integer, allocatable  :: g(:)
real(pr), allocatable :: gr(:)
integer               :: nbin
!variables de MSD
real(pr), allocatable :: r2t, x0, y0, z0
integer, allocatable  :: ix, iy, iz, ntime
integer               :: Nmsd

contains

subroutine initpossc
!posiciones iniciales en una caja cúbica de lado L, para un material de
!densidad \rho con una estructura cristalina SC
implicit none
integer                :: i, j, k, idxp, nucell
real(pr)               :: dist

nucell = nint(real(N,pr)**(1._pr/3._pr))
idxp = 1
do i = 1, nucell            !number of unit cells
    do j = 1, nucell
        do k = 1, nucell
            x(idxp) = real(i-1,pr)
            y(idxp) = real(j-1,pr)
            z(idxp) = real(k-1,pr)
            
            idxp = idxp + 1
        enddo
    enddo
enddo
dist = (1._pr/rho)**(1._pr/3._pr)
x(:) = dist*x(:)
y(:) = dist*y(:)
z(:) = dist*z(:)

end subroutine initpossc


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
