program ising2d
!
!Ising 2D: Cálculo del cumulante de Binder en función de la temperatura
!
!Parámetros optimizados con el inciso anterior:
!          MCS de equilibración = 10000
!          MCS de medición      = 1000000 (puntos con los que se hacen los promedios
!                                           <E>, <E^2>, <M>, <M^2> /N)
!
use precision, only    : pr=>dp, ipr=>k20
use initial_cond, only : T, spin, L, N, order, mess
use calculate
use randomnum, only    : rand=>grnd
implicit none
real(pr), parameter    :: Tc = 2.269185_pr
real(pr)               :: T_max, T_min, dT
integer                :: idx, stepT, mceq, MCS, i, k
integer                :: J, M, H, dH, measurments
integer(ipr)           :: Mag2, Mag4, M2
real(pr)               :: Mag2m, Mag4m, U_L

!#######################     valores iniciales     ############################
J = 1._pr                                                       !ferromagnetico

T_max = 1.5_pr*Tc
T_min = 0.5_pr*Tc
dT    = 0.02_pr
stepT = nint((T_max - T_min)/dT)

L = 40
N = L*L                                                 !cantidad de particulas

allocate( spin(N) )

mceq = 10000
MCS  = 1000000                                                !Monte Carlo Steps

! quenching: desorden -> orden
!call mess()
! heating:   orden -> desorden
call order()
call energy(H)
call magnetization(M)
!###############################################################################

!header
open(44, file='tul.dat', status='replace')
write(44,'(A,I5)') '#largo de la grilla en x e y: ', L 
write(44,'(A,I7)') '#cantidad de particulaes: ', N
write(44,'(A,I9)') '#energía inicial del sistema: ', H
write(44,'(A,I9)')'#magnetización inicial: ', M
write(44,'(A,F4.2)')'#temperatura inicial: ', T_min
write(44,'(A,F4.2)') '#temperatura final: ', T_max
write(44,'(A,F4.2)') '#delta de temperatura: ', dT
write(44,'(A,I5)') '#pasos de equilibración: ', mceq
write(44,'(A,I5)') '#pasos de medición: ', MCS
write(44,'(A)') '# Temp, Energia, Mag, |Mag|, Cv, xi'


do idx = 0, stepT                                          !loop en temperatura
    
    T = T_min + real(idx,pr)*dT                                        !heating
    !T = T_max - real(idx,pr)*dT                                     !quenching
    write(*,'(F4.2)') T    
    do i = 1, mceq                                      !pasos de equilibración
        do k = 1, N
            call deltaE(k,dH)
            if (dH <= 0 .or. rand() <= exp(-dH/T)) then
                spin(k) = -spin(k)
                H = H + dH
                M = M + 2*spin(k)
            endif
        enddo
    enddo

    Mag2 = int(M,ipr)*int(M,ipr)
    Mag4 = Mag2*Mag2

    measurments = 1

    do i = mceq + 1, MCS
                                                                    !paso de MC
        do k = 1, N
                                           !cambio de energía si flipeo un spin 
            call deltaE(k,dH)
                                                                    !metropolis
            if (dH <= 0 .or. rand() <= exp(-dH/T)) then
                spin(k) = -spin(k)                              !actualizo spin
                H = H + dH                    !actualizo la energía del sistema
                M = M + 2*spin(k)
            endif

        enddo

        !acumulo Mag2 y Mag4 para el calculo de U_L
        if (mod(i,100) == 0) then
            measurments = measurments + 1
            M2   = int(M,pr)*int(M,ipr)
            Mag2 = Mag2 + M2
            Mag4 = Mag4 + M2*M2
        endif

    enddo

    Mag2m = real(Mag2,pr)/real(measurments,pr)
    Mag4m = real(Mag4,pr)/real(measurments,pr)

    U_L   = 1._pr - (Mag4m/(3._pr*Mag2m*Mag2m))

    write(44,'(2(E15.6,2x))') T/Tc, U_L/real(N,pr)

enddo


end program
