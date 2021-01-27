program difnum
      !Este programa evalúa la f'(1), donde f(x) = exp(x) con la diferenciación numerica
      !centrada en dos puntos (df) para distintos valores de h que van como 10**(-k)
      !también calcula el error comparando con la derivada exacta
      implicit none
      integer, parameter :: pr = selected_real_kind(13)
      integer            :: k
      real(pr)           :: df, x, h 

      x = 1._pr

      open(27, file='difnum.dat', status='replace')
      write(27, *) "#k, h, f', error"

      do k = 1, 15
          h = 10._pr ** (-k)
          df = (f(x + h) - f(x - h)) / (2._pr * h)
          write(27,'(I2,2x,3(E15.6,2x))') k, h, df, (abs(df - exp(1._pr)))/exp(1._pr)
      enddo


CONTAINS

      function f(xx)
              real(pr), intent(in) :: xx
              real(pr)             :: f

              f = exp(xx)
              
      end function f

end program difnum


