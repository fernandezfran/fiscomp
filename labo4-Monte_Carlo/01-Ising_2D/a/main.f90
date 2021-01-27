program ising2d
!
!Ising 2D: mediciones de la energía y la magnetización en función de los pasos
!          de MC para distintos valores de T y L (que ingresan en la terminal)
!          para distintas condiciones iniciales (orden o desorden).
!
!Valores recomendados de la guía: T = 2.0; 3.3; 2.2676
!                                 L = 40
!                               MCS = 10000                       (valor fijo)
!
use precision, only    : pr=>dp
use initial_cond, only : T, spin, L, N, order, mess
use calculate
use randomnum, only    : rand=>grnd
implicit none
real(pr), parameter    :: Tc = 2.269185_pr
integer                :: J, M, H, dH
integer                :: ciflag, i, k, MCS, kk, kkk
character(len=24)      :: aux_temp

!#######################     valores iniciales     ############################
J = 1._pr                                                       !ferromagnetico
!T = 0.8_pr*Tc                                        !temperatura = kB * T / J
write(*,*) 'Temperatura a la que se desea medir [k_B*T/J]'
read(*,*) T
write(aux_temp,'(F4.2)') T

write(*,*) 'Largo de la grilla en x e y (el nro. de part. será L**2)'
read(*,*) L
N = L*L                                                 !cantidad de particulas

allocate( spin(N) )

write(*,*) 'condiciones inicilaes: 0=up, cualquier entero=desorden'
read(*,*) ciflag
if (ciflag == 0) then !T->0
    call order()
    !open(28, file='EM_'//trim(adjustl(aux_temp))//'-order.dat', status='replace')
else                  !T->\infty
    call mess()
    !open(28, file='EM_'//trim(adjustl(aux_temp))//'-mess.dat', status='replace')
endif
!write(28,*) '# MC step, E/N, M/N', T, '(temperatura)'
open(44, file='configs_eq.dat', status='replace')

call energy(H)
write(*,*) 'energía inicial del sistema', H
call magnetization(M)
write(*,*) 'magnetización inicial', M


MCS = 10000                                                  !Monte Carlo Steps
do i = 1, MCS
                                                                    !paso de MC
    do k = 1, N
                                           !cambio de energía si flipeo un spin 
        call deltaE(k,dH)
                                                                    !metropolis
        if (dH <= 0 .or. rand() <= exp(-dH/T)) then
            spin(k) = -spin(k)                                  !actualizo spin
            H = H + dH                        !actualizo la energía del sistema
        endif

    enddo

    if (mod(i,100) == 0) then
        !escribo archivo para hacer película de spines
        kkk = 1
        !spin(kkk) está guardado en filas
        do k = 1, L
            do kk = 1, L
                write(44,*) k, kk, spin(kkk)
                kkk = kkk + 1
            enddo
            write(44,*)
        enddo
        write(44,*)
        write(44,*)
    endif

    call magnetization(M)
    !write(28,'(I6,2(E15.6,2x))') i, real(H,pr)/real(N,pr), real(M,pr)/real(N,pr)

enddo


end program
