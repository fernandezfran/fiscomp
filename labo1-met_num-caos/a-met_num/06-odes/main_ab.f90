program osc_arm
  !
  !Resuelve un oscilador armónico de k = m = 1 (\omega = 1) usando los métodos de
  !Euler, RK2 y RK4.
  !
  !en el modulo precision.f90 están cada uno de los pasos del método que se desee
  !emplear
  !
  use precision, only : pr=>dp
  use odes
  character(32)         :: archivo, formato
  integer               :: metodo, N, i, i_max
  real(pr)              :: ti, tf, t, h, E_t, r, v
  real(pr), allocatable :: x(:)

!##############              datos iniciales              #####################

  N = 2                                       !dimension del problema

  ti = 0._pr                                  !tiempo inicial
  tf = 10._pr                                 !       final

  allocate(x(N))                              !x inicial (después será pisado)
                                              !posición y velocidad inicial,
  x(1) = 1._pr ; x(2) = 1._pr                 !respectivamente

!##############################################################################

  open(13, file = "in.metodo", status = "old")!en este archivo tiene que haber un 
  read(13,*) metodo                           !1, 2 o 3, según el método que se 
  close(13)                                   !desee emplear

  !los h proporsionados aquí son los optimos para cada método

  select case(metodo)
  case(1)                                     !euler
    h = 10._pr**(-6)
    archivo = "resultados-euler.dat"
    print *, "Usando el método de Euler"
  case(2)                                     !rk2
    h = 10._pr**(-5)
    archivo = "resultados-rk2.dat"
    print *, "Usando el método de RK2"
  case(3)                                     !rk4
    h = 10._pr**(-3)
    archivo = "resultados-rk4.dat"
    print *, "Usando el método de RK4"
  end select

  formato = '(3(E15.6,2x),E18.9,2(E15.6,2x))'

  !E_t = 0.5_pr * ( x(1)*x(1)  + x(2)*x(2) )
  !call sol_ex(ti,r,v)

  open(14, file = archivo, status = "replace")
  write(14,*) "# t, x, v, energia, error_r, error_v"
  !write(14,formato) ti, x(1), x(2), E_t, abs(x(1)-r), abs(x(2)-v)

  i_max = nint( (tf - ti) / h)
  t = ti
  do i = 1, i_max

    select case(metodo)
    case(1)
      call euler(N,h,t,x,f)
    case(2)
      call rk2(N,h,t,x,f)
    case(3)
      call rk4(N,h,t,x,f)
    end select

    t = ti + real(i,pr)*h  !avanzo en el tiempo luego de tener las posiciones
                           !y velocidades siguientes

    call sol_ex(t,r,v)
    E_t = 0.5_pr * ( x(1)*x(1)  + x(2)*x(2) )

    write(14,formato) t, x(1), x(2), E_t, abs(x(1)-r), abs(x(2)-v)

  enddo

contains

  function f(tt,xx)
    !
    !función que condiciona el problema, en este caso fuerza del resorte
    !
    implicit none
    real(pr), intent(in) :: tt, xx(:)
    real(pr)             :: f(size(xx))
    real(pr)             :: k

    k = 1._pr                         !cte del resorte

    f(1) = xx(2)
    f(2) = - k * xx(1)

  end function f

  subroutine sol_ex(tt,pos,vel)
    !
    !acá va la solución exacta del problema, si es que se conoce
    !
    implicit none
    real(pr), intent(in)  :: tt
    real(pr), intent(out) :: pos, vel
    real(pr)              :: A, pi, phi !quizas convenga definir estas variables
                                        !como parameter...

    A = sqrt(2._pr)
    pi = acos(-1._pr)
    phi = - 0.25_pr * pi

    pos = A * cos(tt + phi)
    vel = - A * sin(tt + phi)

  end subroutine sol_ex

end program
