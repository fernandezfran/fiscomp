program odes_p5
  use precision, only : pr=>dp
  use odes
  implicit none
  integer             :: metodo, n, n_max
  real(pr)            :: x0, xf, y0, x, y, h, exa
  character(32)       :: archivo

  !condiciones iniciales (podría leerlo de un archivo o de terminal)
  x0 = 0._pr
  xf = 1._pr

  y0 = 1._pr

  open(13, file = "tar", status = "old")
  read(13,*) metodo
  close(13)

  select case(metodo) !estos h por defecto están bien?
  case(1) !euler
    h = 10._pr**(-4)
    archivo = "resultados-euler.dat"
    print *, "El método de Euler"
  case(2) !rk2
    h = 0.001_pr
    archivo = "resultados-rk2.dat"
    print *, "El método de RK2"
  case(3) !rk4
    h = 0.01_pr
    archivo = "resultados-rk4.dat"
    print *, "El método de RK4"
  end select

  open(14, file = archivo, status = "replace")
  write(14,*) "# n, x, y, exacta, error"

  n = 0
  x = x0
  y = y0
  exa = f(x)

  write(14,'(I5,2x,4(E15.6,2x))') n, x, y, exa, abs(y - exa)/exa
 n_max = nint((xf - x0)/h)
  do n = 1, n_max

    select case(metodo)
    case(1)
      call euler(h,x,y)
    case(2)
      call rk2(h,x,y)
    case(3)
      call rk4(h,x,y)
    end select

    x = x0 + real(n,pr)*h

    exa = f(x)

    write(14,'(I5,2x,4(E15.6,2x))') n, x, y, exa, abs(y - exa)/exa

  enddo

CONTAINS

  function f(xx)
    implicit none
    real(pr), intent(in) :: xx
    real(pr)             :: f

    f = exp(- 0.5_pr * xx * xx)

  end function

end program odes_p5
