module integrator
use precision, only : pr=>dp
use initialization, only : N, x, y, z, ix, iy, iz, fx, fy, fz, dt, pi3, eta, &
                         & D0, pbc, pi, semilla

contains

subroutine ermak
!use randomnum, only : rand => grnd, randseed => sgrnd
use randomnum, only : rand => grnd
implicit none
integer         :: i

do i = 1, N
    x(i) = x(i) + fx(i)*dt/pi3/eta + sqrt(2._pr*D0*dt)*gasdev()
    y(i) = y(i) + fy(i)*dt/pi3/eta + sqrt(2._pr*D0*dt)*gasdev()
    z(i) = z(i) + fz(i)*dt/pi3/eta + sqrt(2._pr*D0*dt)*gasdev()

    call pbc(x(i), ix(i))
    call pbc(y(i), iy(i))
    call pbc(z(i), iz(i))
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

end subroutine ermak

end module
