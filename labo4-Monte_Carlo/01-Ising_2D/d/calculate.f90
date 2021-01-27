module calculate
use precision, only   : pr => dp
use initial_cond, only: T, spin, L, N

contains


subroutine energy(H)
implicit none
integer, intent(out) :: H
integer              :: i, j, sp

H = 0

do i = 1, N
    
     j = i - L                     !vecino de arriba
     if (j < 1) j = N - (L - i) 
     sp = spin(j)
 
     !j = i + L                    !vecino de abajo
     !if (j > N) j = j - N
     !sp = sp + spin(j)
     
     !j = i - 1                    !vecino de la izq
     !if (mod(j,L) == 0) j = j + L
     !sp = sp + spin(j)
 
     j = i + 1                    !vecino de la der
     if (mod(i,L) == 0) j = j - L
     sp = sp + spin(j)

     H = H + spin(i)*sp

enddo

!H = -H/2 !estaba recorriendo todos los pares 2 veces, por eso comenté la mitad
H = -H

end subroutine energy


subroutine magnetization(M)
implicit none
integer, intent(out) :: M

M = sum(spin)

end subroutine magnetization


subroutine deltaE(i,dH)
implicit none
integer, intent(in)  :: i
integer, intent(out) :: dH
integer              :: j

j = i - L
if (j < 1) j = N - (L - i)
dH = spin(j)

j = i + L
if (j > N) j = j - N
dH = dH + spin(j)

j = i - 1
if (mod(j,L) == 0) j = j + L
dH = dH + spin(j)

j = i + 1
if (mod(i,L) == 0) j = j - L
dH = dH + spin(j)

!h_bef =  spin(i)*su
!h_aft = -spin(i)*su
!dH    = - (h_aft - h_bef)

dH = 2*spin(i)*dH

end subroutine deltaE


end module
