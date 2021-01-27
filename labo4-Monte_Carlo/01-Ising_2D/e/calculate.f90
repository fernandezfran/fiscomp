module calculate
use precision, only   : ipr => k18, pr=>dp
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


!SUBRUTINA DE AUTOCORRELACION
!#########################comentarios del practico###############################
!funcion autocorrelación para la energía (HH) y el modulo de la magnetización (MM)
!
! C(tau) = <A(t+tau)*A(t)>/<A(t)**2>, para que en tau = 0 sea igual a 1.
! el cálculo de A(t)**2 es trivial, está hecho. No importa a donde empiece a medir
!  luego de que haya hecho la equilibración. t=MCS (conviene todos juntos en vez
!  de separados), quiero Npts (=1000, para empezar) de corr dados desde antes.
!  En el codigo C'(tau) = <A(t)A(t-tau)>
!
!  C(1:n) = 0
!En t=i: 
!  do j = 1, n !pero esto no se hace hasta que C(j) no esté lleno
!      C(j) = C(j) + A(j)*A(i-j)
!################################################################################
!
!Paso la energía y el modulo de la magnetización del paso de MC en el que estoy,
!la catidad de pasos de correlación (tau = dim de Ecorr y Mcorr) y hago que el
!vector sea ciclico utilizando la función mod(,) apropiadamente.
!
!La normalización la hago en el main y ahí imprimo a un archivo de salida
!
subroutine autocorr(E,E2,modM,M2,ncorr,EE,MM,Ecorr,Mcorr)
integer(ipr), intent(in)    :: E, E2, modM, M2
integer, intent(in)         :: ncorr
integer(ipr), intent(inout) :: EE(ncorr), MM(ncorr)
real(pr), intent(inout)     :: Ecorr(ncorr), Mcorr(ncorr)
integer                     :: ilast, inext!, EE(ncorr), MM(ncorr)
integer, save               :: ic

ic = ic + 1
ilast = mod(ic-1,ncorr) + 1
EE(ilast) = E      !la energía en ese punto que entra como argumento
MM(ilast) = modM   !y la magnetización

if (ic <= ncorr) return 

do j = 1, n
    inext    = mod(ic-j,n) + 1
    Ecorr(j) =  Ecorr(j) + real(EE(ilast)*EE(inext),pr)/real(E2,pr)
    Mcorr(j) =  Mcorr(j) + real(MM(ilast)*MM(inext),pr)/real(M2,pr)
enddo

end subroutine autocorr


end module
