program ljfluid
use precision     , only : pr => dp
use randomnum     , only : ranf => grnd
use initialization, only : x, y, z, T, N, L, V, rho, rcut, initial_cond => xyzfcc, pbc
use thermo        , only : deltau
implicit none
real(pr), parameter      :: pi = (4._pr*atan(1._pr)), pi43 = (4._pr/3._pr)*pi
integer                  :: i, mceq, mcs, j
integer                  :: part, accept, adjust
real(pr)                 :: delta, xnew, ynew, znew, dU, utail, ptail, tasa, cord

T    = 0.9_pr               != T*
rho  = 0.8_pr               != rho* = N/V
N    = 125
V    = real(N,pr)/rho
L    = V**(1._pr/3._pr)
rcut = 2.5_pr
!para estos calculos estoy usando epsilon = sigma = 1, por ahora
utail = 2._pr * pi43 * rho * ((1._pr/3._pr)*rcut**(-9) + rcut**(-3))
ptail = 4._pr * pi43 * rho * rho * ((2._pr/3._pr)*rcut**(-9) + rcut**(-3))

allocate( x(N), y(N), z(N) )
call initial_cond()

!MCS de termalización, ajuste del delta
accept = 0
adjust = 10
mceq   = 2000
delta  = .05_pr
open(30, file='equilibration.xyz', status='replace')
write(30,*) N; write(30,*)
do j = 1, N; write(30,'(A,3(E15.6,2x))') 'Ar', x(j), y(j), z(j); enddo


do i = 1, mceq

    do j = 1, N
        !selecciono particula al azar
        part = int(ranf()*N) + 1
        !desplazo la particula
        xnew = x(part) + delta*(ranf() - 0.5_pr)
        ynew = y(part) + delta*(ranf() - 0.5_pr)
        znew = z(part) + delta*(ranf() - 0.5_pr)
        !calculo el cambio de energía
        call deltau(part,xnew,ynew,znew,dU)
        !metropolis
        if (dU <= 0._pr .or. ranf() <= exp(-dU/T)) then
            !si acepto el movimiento => aplico PBC
            call pbc(xnew)
            x(part) = xnew
            call pbc(ynew)
            y(part) = ynew
            call pbc(znew)
            z(part) = znew
            accept  = accept + 1
        endif
    enddo

    !ajusto cada `adjust' pasos de MC
    if (mod(i,adjust) == 0) then
        tasa = real(accept,pr)/real(adjust*N,pr)
        if (tasa > 0.2_pr) then
            delta = 1.05_pr*delta
        else
            delta = 0.95_pr*delta
        endif
        accept = 0
    endif

    !write .xyz
    if (mod(i,100) == 0) then
        write(30,*) N; write(30,*)
        do j = 1, N; write(30,'(A,3(E15.6,2x))') 'Ar', x(j), y(j), z(j); enddo
    endif

enddo
write(*,*) 'el desplazamiento optimo es,', delta

end program
