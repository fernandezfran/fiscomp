program libr
  use precision, only: pr=>dp
  implicit none
  integer            :: i
  real(pr)           :: x,y
  integer(8)         :: j

  j = fftw_plan_dft_r2c_1d(128, x, y, i) !no la estoy usando bien pero no me importa


end program
