program evsh
  !
  !Analisis del error global en función de 1/h para el oscilador armonico para
  !t entre 0 y 10
  !
  use precision, only : pr=>dp
  use odes
  character(32)         :: formato
  integer               :: N, i, j, k, i_max
  real(pr)              :: ti, tf, h, r, v
  real(pr), allocatable :: x_eu(:), x_rk2(:), x_rk4(:), t(:)

!##############              datos iniciales              #####################

  N = 2                                       !dimension del problema

  ti = 0._pr                                  !tiempo inicial
  tf = 10._pr                                 !       final  

!##############################################################################

  formato = '(4(E15.6,2x))'

  allocate(x_eu(N), x_rk2(N), x_rk4(N))  !separo el x de cada método para que no
                                         !se "pisen" las soluciones

  allocate(t(3)) !3 metodos, una entrada del vector para cada uno de ellos 

  open(14, file = "errores.dat", status = "replace")
  write(14,*) "# h, e_euler, e_rk2, e_rk4"

  !este loop devuelve los valores de 1/h desordenados, después hay que hacer un sort
  do k = 0,6
    do j = 1, 9

      h = 0.1_pr * real(j,pr) * 10._pr**(-k)
      
      write(*,*) k, j, h

      i_max = nint( (tf - ti) / h)
    
      t(:) = (/ti,ti,ti/)
      x_eu(:) = (/1._pr,1._pr/)  
      x_rk2(:) = (/1._pr,1._pr/)
      x_rk4(:) = (/1._pr,1._pr/)

      do i = 1, i_max

        call euler(N,h,t(1),x_eu,f)
        call rk2(N,h,t(2),x_rk2,f)
        call rk4(N,h,t(3),x_rk4,f)

        t(:) = ti + real(i,pr)*h

      enddo

      call sol_ex(t(1),r,v)
      write(14,formato) h, abs((x_eu(1) - r)/r), abs((x_rk2(1) - r)/r), abs((x_rk4(1) - r)/r)
    enddo
  enddo

contains

  function f(tt,xx)
    !
    !función que condiciona el problema, en este caso fuerza del resorte
    !
    implicit none
    real(pr), intent(in) :: tt, xx(:)
    real(pr)             :: f(size(xx))

    f(1) = xx(2)
    f(2) = - xx(1)

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
