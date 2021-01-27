program ising2d
!
!Ising 2D: Histograma de magnetización para temperaturas un poco por debajo
!          y un poco por encima de Tc. 
!
!Valores recomendados de la guía: T1 = 0.88*Tc y T2 = 1.15*Tc
!
use precision, only    : pr=>dp
use initial_cond, only : T, spin, L, N, order, mess
use calculate
use randomnum, only    : rand=>grnd
implicit none
real(pr), parameter    :: Tc = 2.269185_pr
real(pr)               :: a
integer                :: J, M, H, dH
integer                :: ciflag, i, k, MCS, mceq, med, measurments
character(len=24)      :: aux_temp, aux_L
integer                :: ibin, nbin
integer, allocatable   :: hist(:)

!#######################     valores iniciales     ############################
J = 1._pr                                                       !ferromagnetico
write(*,*) 'valor por debajo o por encima de Tc, a=? (T = a*Tc)'
read(*,*) a
T = a*Tc
write(aux_temp,'(F4.2)') T

write(*,*) 'Largo de la grilla en x e y (el nro. de part. será L**2)'
read(*,*) L
write(aux_L,'(I3)') L
N = L*L                                                 !cantidad de particulas

allocate( spin(N) )

write(*,*) 'condiciones inicilaes: 0=up, cualquier entero=desorden'
read(*,*) ciflag
if (ciflag == 0) then !T->0
    call order()
else                  !T->\infty
    call mess()
endif

write(*,*) 'pasos de equilibración'
read(*,*) mceq
write(*,*) 'pasos de Monte Carlo'
read(*,*) MCS
write(*,*) 'cada cuanto se realiza una medición'
read(*,*) med

nbin = 201
allocate( hist(nbin) )
hist(:) = 0
!##############################################################################

open(37, file='histoM_'//trim(adjustl(aux_temp))//'-'//trim(adjustl(aux_L))//'.dat', status='replace')
write(37,*) '# temperatura ', T, 'a la que se realiza el histograma de M'
write(37,*) '# nro de particulas', N
write(37,*) '# cantidad de pasos de eq', mceq
write(37,*) '# cantidad de pasos de MC', MCS
write(37,*) '# cada cuantos pasos se realiza una medición', med
call energy(H)
write(37,*) '# energía inicial del sistema', H
call magnetization(M)
write(37,*) '# magnetización inicial', M

do i = 1, mceq
    do k = 1, N
        call deltaE(k,dH)
        if (dH <= 0 .or. rand() <= exp(-dH/T)) then
            spin(k) = -spin(k)                                  !actualizo spin
            H = H + dH                        !actualizo la energía del sistema
            M = M + 2*spin(k)
        endif
    enddo
enddo

do i = 1, MCS
    do k = 1, N
        call deltaE(k,dH)
        if (dH <= 0 .or. rand() <= exp(-dH/T)) then
            spin(k) = -spin(k)                                  !actualizo spin
            H = H + dH                        !actualizo la energía del sistema
            M = M + 2*spin(k)
        endif
    enddo

    if (mod(i,med) == 0) then
        !M/N puede tomar valores entre -1 y 1 => ibin entre -100 y 100
        ibin = nint(real(100*M,pr)/real(N,pr))
        hist(ibin + 101) = hist(ibin + 101) + 1
    endif
enddo

do k = 1, nbin
    write(37,'(2(E15.6,2x))') real(k-101,pr)/100._pr, real(hist(k),pr)/real(nbin,pr)
enddo


end program
