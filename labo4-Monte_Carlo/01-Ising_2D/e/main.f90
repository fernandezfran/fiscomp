program ising2d
!
!Ising 2D: Calculo la función autocorrelación para la energía y para el modulo
!            de la magnetización en un tamaño y una temperatura dadas
!
!Parámetros optimizados con el inciso anterior:
!          MCS de equilibración = 10000
!          MCS de medición      = 1000000 (puntos con los que se hacen los promedios
!                                           <E>, <E^2>, <M>, <M^2> /N)
!
use precision, only    : pr=>dp, ipr=>k18
use initial_cond, only : T, spin, L, N, order, mess
use calculate
use randomnum, only    : rand=>grnd
implicit none
real(pr), parameter    :: Tc = 2.269185_pr
real(pr)               :: E2m, M2m
integer                :: mceq, MCS, i, k, ncorr
integer                :: J, M, H, dH, measurments
integer(ipr)           :: E, E2, modM, M2
integer(ipr), allocatable :: EE(:), MM(:)
real(pr), allocatable     :: Ecorr(:), Mcorr(:)

!#######################     valores iniciales     ############################
J = 1._pr                                                       !ferromagnetico
T = 1.5_pr*Tc
L = 10
N = L*L                                                 !cantidad de particulas

allocate( spin(N) )

mceq = 1000
MCS  = 10000                                                 !Monte Carlo Steps

ncorr = 100
allocate( Ecorr(ncorr), Mcorr(ncorr), EE(ncorr), MM(ncorr) )
Ecorr(:) = 0._pr; Mcorr(:) = 0._pr; EE(:) = 0; MM(:) = 0

call mess()
!call order()
call energy(H)
call magnetization(M)
!###############################################################################

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

measurments = 1
E    = int(H,ipr)
E2   = E*E
modM = abs(int(M,ipr))
M2   = modM*modM

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

    measurments = measurments + 1
    E    = E + int(H,ipr)
    E2   = E2 + int(H,ipr)*int(H,ipr)
    modM = modM + abs(int(M,ipr))
    M2   = M2 + int(M,ipr)*int(M,ipr)

    call autocorr(E,E2,modM,M2,ncorr,EE,MM,Ecorr,Mcorr)
enddo

!Ecorr(:) = Ecorr(:)/real(measurments,pr)
!Mcorr(:) = Mcorr(:)/real(measurments,pr)

do i = 1, ncorr
    write(*,'(I5,2(E15.6,2x))') i, Ecorr(i), Mcorr(i)
enddo

end program
