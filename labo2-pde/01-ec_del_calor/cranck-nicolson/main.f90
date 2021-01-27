program ec_calor
!
! ecuación del calor 1-d para condiciones de contorno de Dirichlet
! utilizando el método de Cranck-Nicolson
!
use precision, only    : pr=>dp
use matrix_solver
implicit none
real(pr), parameter   :: pi = 4._pr * atan(1._pr)
real(pr)              :: L, T_0, D, dt, dx
integer               :: n, i, j, m, timemes
real(pr)              :: T_inicial, K, C, rho, eta, t
real(pr), allocatable :: x(:), T_antes(:), T_despues(:), T_exacta(:)
real(pr), allocatable :: diag(:), a_i(:), c_i(:), b_i(:)
character(len=32)     :: filename, frmt
real(pr)              :: time1, time2

!########    valores de entrada que podrían pedirse de la terminal    ########
n         = 100                        ! divisiones del eje x

L         = 1._pr                      ! Largo de la barra [metros]
T_0       = 0._pr                      ! Temperatura en los bordes [°C]
T_inicial = 100._pr                    ! Temperatura inicial de la barra [°C]
K         = 237._pr                    ! Conductividad términa [W/(m K)]
C         = 900._pr                    ! Calor especifico [J/(kg K)]
rho       = 2700._pr                   ! densidad (del aluminio) [kg/(m**3)]

D         = K/(C*rho)                  ! Coeficiente de difusion [m**2/s]

dt        = 0.3_pr                     ! delta del tiempo que ingresa desde la
                                       ! terminal [s]

dx        = L/(real(n+1,pr))           ! delta de x
!#############################################################################

! adimensionalización del tiempo y el espacio
dt = (D/(L*L)) * dt                    ! dt'
dx = dx/L                              ! dx'

eta = dt/(dx*dx)


allocate( x(n), T_antes(n), T_despues(n), T_exacta(n) )
allocate( diag(n), a_i(2:n), c_i(n-1), b_i(n) )
!defino los valores de la matriz que van a permanecer constantes
diag(:) = 2._pr*(1._pr + 1._pr/eta) 
a_i(:) = -1._pr
c_i(:) = -1._pr

! discretización del espacio
do i = 1, n
    x(i) = real(i,pr)*dx
enddo

! inicializo la temperatura de la barra (T(0) = T(n+1) = T_0)
T_antes(:) = T_inicial

! inicializo el tiempo
t = 0._pr

filename = 'temp_cn.dat' 
frmt = '(5(E15.6,2x))'
open(55, file=filename, status='replace')

!main loop - cranck-nicolson y exacto
!escribo condiciones iniciales en el archivo
write(55,*) '#tiempo [s], x [m], Temperatura [K], Valor Exacto [K], error absoluto [K]'
write(55,frmt) t, 0._pr, T_0, T_0, T_0-T_0
do m = 1, n
    write(55,frmt) t, x(m), T_antes(m), T_inicial, T_antes(m) - T_inicial
enddo
write(55,frmt) t, L, T_0, T_0, T_0-T_0
write(55,*) !espacio en blanco par GNUplot

do j = 1, 20000

    t = t + dt
    
    b_i(1) = 2._pr*T_0 + 2._pr*(1._pr/eta - 1._pr)*T_antes(1) + T_antes(2)
    do i = 2, n-1
        b_i(i) = T_antes(i-1) + 2._pr*(1._pr/eta - 1._pr)*T_antes(i) + T_antes(i+1)
    enddo
    b_i(n) = T_antes(n-1) + 2._pr*(1._pr/eta - 1._pr)*T_antes(n) + 2._pr*T_0

    call tridiag(n,a_i,diag,c_i,b_i,T_despues)

    if (mod(j,300) .eq. 0) then !escribo cada 300 pasos
        
        call exacta(n,x(:),t,T_exacta(:))
        
        write(55,frmt) L*L*t/D, 0._pr, T_0, T_0, T_0-T_0
        do m = 1, n
            write(55,frmt) L*L*t/D, L*x(m), T_despues(m), T_exacta(m), abs(T_despues(m) - T_exacta(m))
        enddo
        write(55,frmt) L*L*t/D, L, T_0, T_0, T_0-T_0
        write(55,*)

    endif

    T_antes(:) = T_despues(:)

enddo


!                MEDICIÓN DEL TIEMPO DE COMPUTO DE CADA MÉTODO
write(*,*) '¿queres medir el tiempo del explicito y el exacto (0=sí, cualquier otro numero=no)?'
read(*,*) timemes
if (timemes.eq.0) then
    T_antes(:) = T_inicial
    t = 0._pr
    call cpu_time(time1)
    do j = 1, 20000
        t = t + dt
    
        b_i(1) = 2._pr*T_0 + 2._pr*(1._pr/eta - 1._pr)*T_antes(1) + T_antes(2)
        do i = 2, n-1
            b_i(i) = T_antes(i-1) + 2._pr*(1._pr/eta - 1._pr)*T_antes(i) + T_antes(i+1)
        enddo
        b_i(n) = T_antes(n-1) + 2._pr*(1._pr/eta - 1._pr)*T_antes(n) + 2._pr*T_0

        call tridiag(n,a_i,diag,c_i,b_i,T_despues)
       
       T_antes(:) = T_despues(:)
    enddo
    call cpu_time(time2)
    write(*,*) 'metodo c-n   ', time2 - time1

    t = 0._pr
    call cpu_time(time1)
    do j = 1, 20000
        t = t + dt
        call exacta(n,x(:),t,T_exacta(:))
    enddo
    call cpu_time(time2)
    write(*,*) 'metodo exacto', time2 - time1
endif
!


contains


subroutine exacta(nn,xx,tt,T_exa)
!
!la solución exacta lee de arriba la temperatura inicial
!
implicit none
integer, intent(in)   :: nn
real(pr), intent(in)  :: xx(nn), tt
real(pr), intent(out) :: T_exa(nn)
real(pr)              :: sn(nn), ex(nn)
integer               :: mm

T_exa(:) = 0._pr
do mm = 1, 100, 2 !50 términos
    sn(:) = sin(mm*pi*x(:))/(mm*pi)
    ex(:) = exp(-mm*mm*pi*pi*tt)     !adimensional

    T_exa(:) = T_exa(:) + 4._pr*T_inicial*sn(:)*ex(:)
enddo

end subroutine exacta

end program
