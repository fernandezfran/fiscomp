module odes
use precision, only : pr=>dp
CONTAINS

  subroutine euler(h,x,y)
    implicit none
    real(pr), intent(in)    :: h, x
    real(pr), intent(inout) :: y
    real(pr)                :: fxy

    call func(x,y,fxy)
    y = y + h * fxy

  end subroutine euler


  subroutine rk2(h,x,y)
    implicit none
    real(pr), intent(in)    :: h, x
    real(pr), intent(inout) :: y
    real(pr)                :: k_1, k_2, fxy

    call func(x, y, fxy)
    k_1 = h * fxy

    call func(x + 0.5_pr*h, y + 0.5_pr * k_1 , fxy)
    k_2 = h * fxy

    y = y + k_2

  end subroutine rk2


  subroutine rk4(h,x,y)
    implicit none
    real(pr), intent(in)    :: h, x
    real(pr), intent(inout) :: y
    real(pr)                :: k_1, k_2, k_3, k_4, fxy

    call func(x, y, fxy)
    k_1 = h * fxy

    call func(x + 0.5_pr*h, y + 0.5_pr * k_1 , fxy)
    k_2 = h * fxy

    call func(x + 0.5_pr*h, y + 0.5_pr * k_2 , fxy)
    k_3 = h * fxy

    call func(x + h, y + k_3, fxy)
    k_4 = h * fxy

    y = y + (k_1 + 2._pr * (k_2 + k_3) + k_4)/6._pr

  end subroutine rk4

  subroutine func(x,y,fxy)
    real(pr), intent(in)  :: x, y
    real(pr), intent(out) :: fxy

    fxy = -x * y

    return
  end subroutine func

end module odes
