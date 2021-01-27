module odes
use precision, only : pr=>sp
CONTAINS

  subroutine euler(n,h,t,x,f)
    implicit none
    integer, intent(in)     :: n
    real(pr), intent(in)    :: h, t
    real(pr), intent(inout) :: x(n)
    interface
      function f(t,x)
        use precision, only : pr=>sp
        real(pr), intent(in) :: t, x(:)
        real(pr)             :: f(size(x))
      end function f
    end interface

    x(:) = x(:) + h * f(t,x)

  end subroutine euler


  subroutine rk2(n,h,t,x,f)
    implicit none
    integer, intent(in)     :: n
    real(pr), intent(in)    :: h, t
    real(pr), intent(inout) :: x(n)
    real(pr)                :: k_1(n), k_2(n)
    interface
      function f(t,x)
        use precision, only : pr=>sp
        real(pr), intent(in) :: t, x(:)
        real(pr)             :: f(size(x))
      end function f
    end interface

    k_1(:) = h * f(t,x)
    k_2(:) = h * f(t + 0.5_pr * h, x + 0.5_pr * k_1)

    x(:) = x(:) + k_2(:)

  end subroutine rk2


  subroutine rk4(n,h,t,x,f)
    implicit none
    integer, intent(in)     :: n
    real(pr), intent(in)    :: h, t
    real(pr), intent(inout) :: x(n)
    real(pr)                :: k_1(n), k_2(n), k_3(n), k_4(n)
    interface
      function f(t,x)
        use precision, only : pr=>sp
        real(pr), intent(in) :: t, x(:)
        real(pr)             :: f(size(x))
      end function f
    end interface

    k_1(:) = h * f(t,x)
    k_2(:) = h * f(t + 0.5_pr * h, x + 0.5_pr * k_1)
    k_3(:) = h * f(t + 0.5_pr * h, x + 0.5_pr * k_2)
    k_4(:) = h * f(t + h, x + k_3)

    x(:) = x(:) + (k_1(:) + 2._pr * (k_2(:) + k_3(:)) + k_4(:))/6._pr

  end subroutine rk4

end module odes
