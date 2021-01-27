program ising2d
!
!Ising 2D: Cálculo de la energía, la magnetización, el calor especifico y la
!          suceptibilidad magnética en función de la temperatura en el intervalo
!          [0,3.3]
!
!Parámetros optimizados con el inciso anterior:
!          MCS de equilibración = 10000
!          MCS de medición      = 1000000 (puntos con los que se hacen los promedios
!                                        <E>, <E^2>, <M>, <M^2> /N)
!
use precision, only    : pr=>dp, ipr=>k18
use initial_cond, only : T, spin, L, N, order, mess
use calculate
use randomnum, only    : rand=>grnd
implicit none
real(pr), parameter    :: Tc = 2.269185_pr
real(pr)               :: T_max, T_min, dT
integer                :: idx, stepT, mceq, MCS, i, k
integer                :: J, M, H, dH, measurments
integer(ipr)           :: E, E2, Mag, Mag2, modM
real(pr)               :: Em, E2m, Magm, Mag2m, modMm, CV, xi

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
open(44, file='tmxiecv.dat', status='replace')
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

    E    = int(H,ipr)
    E2   = E*E
    Mag  = int(M,ipr)
    modM = abs(Mag)
    Mag2 = Mag*Mag

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

        !acumulo E, E2, Mag, Mag2, modM
        if (mod(i,100) == 0) then
            measurments = measurments + 1
            E    = E + int(H,ipr)
            E2   = E2 + int(H,ipr)*int(H,ipr)
            Mag  = Mag + int(M,ipr)
            modM = modM + abs(int(M,ipr))
            Mag2 = Mag2 + int(M,ipr)*int(M,ipr)
        endif

    enddo

    Em    = real(E,pr)/real(measurments,pr)
    E2m   = real(E2,pr)/real(measurments,pr)
    Magm  = real(Mag,pr)/real(measurments,pr)
    modMm = real(modM,pr)/real(measurments,pr)
    Mag2m = real(Mag2,pr)/real(measurments,pr)

    CV = (E2m - Em*Em)/(real(N*N,pr)*T*T)
    xi = (Mag2m - modMm*modMm)/(real(N*N,pr)*T)

    write(44,'(6(E15.6,2x))') T/Tc, Em/real(N,pr), Magm/real(N,pr), modMm/real(N,pr), CV, xi

enddo


end program
