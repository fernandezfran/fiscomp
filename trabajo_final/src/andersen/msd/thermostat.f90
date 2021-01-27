module thermostat
use precision, only : pr=>dp
use initialization, only : pi, vx, vy, vz, nu, dt, T_damp, T0, Temp, N
implicit none
contains

subroutine berendsen_temp
implicit none
integer             :: i
real(pr)            :: lamda, damp_factor, temp_factor

damp_factor = dt/T_damp
temp_factor = T0/Temp
lamda = sqrt(1._pr + damp_factor*(temp_factor - 1._pr))

do i = 1, N
    vx(i) = vx(i)*lamda
    vy(i) = vy(i)*lamda
    vz(i) = vz(i)*lamda
enddo

end subroutine berendsen_temp


subroutine andersen_temp
use randomnum, only : rmzran, rand => grnd
implicit none   
integer             :: i 
real(pr)            :: sigma

sigma = sqrt(T0)
do i = 1, N
    if (rmzran() .lt. nu*dt) then
        vx(i) = gasdev()*sigma
        vy(i) = gasdev()*sigma
        vz(i) = gasdev()*sigma
    endif
enddo


contains

    function gasdev()
    ! método de box-muller
    !    se considera x0 = 0 y sigma = 1
    implicit none
    real(pr)               :: gasdev, gset
    real(pr)               :: u1, u2, fs, ang
    integer                :: iset = 0
    save iset, gset
    
    if (iset == 0) then
        u1  = rand(); u2  = rand()
        fs  = sqrt(-2._pr*log(u1))
        ang = 2._pr*pi*u2
    
        gset   = fs*cos(ang)
        gasdev = fs*sin(ang)
    
        iset = 1
    else
        gasdev = gset
        iset = 0
    endif

    endfunction gasdev

end subroutine andersen_temp

end module
