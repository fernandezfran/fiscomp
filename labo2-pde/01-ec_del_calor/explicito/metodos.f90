module metodos
!
!modulo con metodos para resolver la ecuación del calor
!
use precision, only : pr=>dp

contains

subroutine explicito(n,T_0,eta,T_antes,T_despues)
!
!este metodo necesita los datos de entrada:
!n   = division del espacio
!T_0 = temperatura en los bordes
!eta = factor (dt/dx/dx)
!T_antes y T_despues
!
implicit none
integer, intent(in)   :: n
real(pr), intent(in)  :: T_0, eta, T_antes(n)
real(pr), intent(out) :: T_despues(n)
integer               :: i

T_despues(1) = T_antes(1) + eta*(T_antes(2) + T_0 - 2._pr*T_antes(1))
do i = 2, n-1
    T_despues(i) = T_antes(i) + eta*(T_antes(i+1) + T_antes(i-1) - 2._pr*T_antes(i))
enddo
T_despues(n) = T_antes(n) + eta*(T_0 + T_antes(n-1) - 2._pr*T_antes(n))


return
end subroutine explicito


end module metodos
