module dplpmc
!d(istribucion) (de) p(robabilidad): l(ey) (de) p(otecncias) (para) MC
use precision, only :   pr => dp
use randomnum, only : rand => rmzran

contains

function ldp()
!
!distribución de nros aleatorios como ley de potencia, p(x) = (k+1)*x^k
!
implicit none
integer  :: k
real(pr) :: power, ldp

k     = 2
power = 1._pr/real(k+1,pr)

ldp = rand()**power

end function ldp


end module
