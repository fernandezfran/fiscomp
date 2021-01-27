module thermostat
use precision, only : pr=>dp
use initialization, only : vx, vy, vz, dt, tau_T, T0, Temp, N
implicit none
contains

subroutine berendsen_temp
implicit none
integer             :: i
real(pr)            :: lamda, coup_factor, temp_factor

coup_factor = dt/tau_T
temp_factor = T0/Temp
lamda = sqrt(1._pr + coup_factor*(temp_factor - 1._pr))

do i = 1, N
    vx(i) = vx(i)*lamda
    vy(i) = vy(i)*lamda
    vz(i) = vz(i)*lamda
enddo

end subroutine berendsen_temp

end module
