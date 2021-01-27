program vs
  !
  ! COMPARACIÓN DE FFTW vs DFT: graficar cpu_time vs N y ver que va como 
  ! O(n log_2 n) y O(n**2)
  !
  use, intrinsic     :: iso_c_binding    !modulo intrinseco que ya está en F (2003)
  use precision, only: pr=>dp
  use dftmod
  implicit none
  include "/usr/include/fftw3.f03"
  real(pr), parameter                    :: pi = 4._pr * atan(1._pr)
  real(C_double), allocatable            :: fxr(:), fxi(:), tfxr(:), tfxi(:)
  complex(C_double_complex), allocatable :: tfx(:)
  integer                                :: N, i, j, k
  type(C_ptr)                            :: plan_rc
  real(pr)                               :: t, rango, factor, dft_s, dft_f, fftw_s, fftw_f

  
  rango  = 4._pr
  factor = 1._pr/rango

  write(*,*) "N, tiempo fftw, tiempo dft"
  do k = 0, 15
      N = 2**k
 
      allocate( fxr(N), fxi(N), tfx(N/2 + 1), tfxr(N), tfxi(N) )

      fxi(:) = 0._pr
      do i = 1, N
          t = rango * real(i-1, pr)/real(N-1, pr)
          fxr(i) = f(t)
      enddo
      
      call cpu_time(fftw_s)
      plan_rc = fftw_plan_dft_r2c_1d(N, fxr, tfx, FFTW_MEASURE)
      call fftw_execute_dft_r2c(plan_rc, fxr, tfx)
      call fftw_destroy_plan(plan_rc)
      call cpu_time(fftw_f)

      call cpu_time(dft_s)
      call dft(fxr,fxi,tfxr,tfxi,N)
      call cpu_time(dft_f)

      write(*,*) N, fftw_f - fftw_s, dft_f - dft_s

      deallocate( fxr, fxi, tfx, tfxr, tfxi )

  enddo

contains

  function f(tt)
    !función que se desea transformar
    implicit none
    real(pr), intent(in) :: tt
    real(pr)             :: f

    f = sin(pi * 0.5_pr * tt) + cos(20._pr * pi * tt)

  end function f

end program
