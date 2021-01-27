program g1pb1
  !
  !primer código para calcular la transformada de fourier a una suma de sen y cos
  !y reobtenerla a partir del anti tranformada, para aprender a usar la lib
  !FFTW3 (hay que usar el flag -lfftw3 al compilar)
  !
  use, intrinsic     :: iso_c_binding    !modulo intrinseco que ya está en F (2003)
  use precision, only: pr=>dp
  implicit none
  include "/usr/include/fftw3.f03"
  real(pr), parameter                    :: pi = 4._pr * atan(1._pr)
  real(C_double), allocatable            :: fx(:), atfx(:)
  complex(C_double_complex), allocatable :: tfx(:), cfx(:)
  integer                                :: N, i, j
  type(C_ptr)                            :: plan_rc, plan_cr
  real(pr)                               :: t, rango, factor

  write(*,*) 'número de puntos para la FFT?'
  read(*,*) N

  rango  = 4._pr
  factor = 1._pr/rango

!################################TRANSF DE FOURIER##################################
  allocate( fx(N), tfx(N/2 + 1) ) !importante que esté antes del plan, para que ya
                                  !tenga asignada memoria in and out para el plan...

  plan_rc = fftw_plan_dft_r2c_1d(N, fx, tfx, FFTW_MEASURE)

  allocate( cfx(N/2 + 1), atfx(N) )
  plan_cr = fftw_plan_dft_c2r_1d(N, cfx, atfx, FFTW_MEASURE)

  open(33, file='fx.dat', status='replace')
  do i = 1, N
    t     = rango * real(i-1, pr)/real(N-1, pr)                   !N pts entre 0 y 4
    fx(i) = f(t)
    write(33,*) t,fx(i)
  enddo
  close(33)

  !hacemos la transformada de Fourier fftw_plan_dft_r2c_1d
  write(*,*) "transformando..."
  call fftw_execute_dft_r2c(plan_rc, fx, tfx)

  open(34,file='tfx.dat',status='replace')
  do i = -N/2, N/2
    if ( i < 0 ) then
      write(34,*) factor*real(i, pr), real(conjg(tfx(-i+1)/real(N,pr))), &
                  & aimag(conjg(tfx(-i+1)/real(N,pr)))
    else
      write(34,*) factor*real(i, pr), real(tfx(i+1)/real(N,pr)), &
                  &  aimag(tfx(i+1)/real(N,pr))
    endif
  enddo
  close(34)

  !hacemos la anti-transf. de F fftw_plan_dft_c2r_1d
  cfx(:) = tfx(:)
  write(*,*) "anti-transformando..."
  call fftw_execute_dft_c2r(plan_cr, cfx, atfx)

  open(35, file='atfx.dat', status='replace')
  do i = 1, N
    t = 4._pr * real(i-1,pr)/real(N-1,pr)
    write(35,*) t, atfx(i)/real(N,pr)          !el eje x tiene que ser omega (w)
  enddo
  close(35)

  call fftw_destroy_plan(plan_rc)
  call fftw_destroy_plan(plan_cr)

  contains

  function f(tt)
    !función que se desea transformar
    implicit none
    real(pr), intent(in) :: tt
    real(pr)             :: f

    f = sin(pi * 0.5_pr * tt) + cos(20._pr * pi * tt)

  end function f

end program
