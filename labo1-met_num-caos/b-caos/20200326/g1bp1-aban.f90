     program g1bp1
     use, intrinsic :: iso_c_binding
     use precision, only: pr => dp
     implicit none
     include "/usr/include/fftw3.f03"
     real(pr), parameter                    :: pi=4._pr*atan(1._pr)
     real(C_double), allocatable            :: fx(:),fcr(:)
     complex(C_double_complex), allocatable :: tfx(:),cfx(:)
     integer                                :: N,i,j
     type(C_ptr)                            :: plan_rc, plan_cr
     real(pr)                               :: t,factor
     
     write(*,*) "pi = ", pi
     
     write(*,*) 'numero de puntos para la FFT?'
     read(*,*)  N
     
     allocate( fx(N), tfx(N/2+1), fcr(N), cfx(N/2+1) )
          
     
     plan_rc = fftw_plan_dft_r2c_1d(N, fx, tfx,FFTW_MEASURE)
     plan_cr = fftw_plan_dft_c2r_1d(N, cfx,fcr,FFTW_MEASURE)
     
     open(33,file="fx.d",status='replace')
     do i = 1,N
       t  = 4._pr*real(i-1,pr)/real(N-1,pr)
       fx(i) = 2._pr*sin(pi*0.5_pr*t) + cos(20._pr*pi*t)
       write(33,*) t,fx(i)
     enddo
     close(33)
     
     ! hacemos la transformada de Fourier fftw_plan_dft_r2c_1d
     call fftw_execute_dft_r2c(plan_rc, fx, tfx)
     
     cfx(:) = tfx(:)
     
     open(33,file="tfx.d",status='replace')
     factor = 1._pr/4._pr
     do i = -N/2,N/2
       if ( i < 0) then
         write(33,*) factor*real(i,pr),real( conjg( tfx(-i+1) )/real(N,pr) ),aimag(conjg( tfx(-i+1) )/real(N,pr) )
       else
         write(33,*) factor*real(i,pr),real( tfx(i+1)/real(N,pr) ),aimag( tfx(i+1)/real(N,pr) )
       endif
     enddo
     close(33)
     
     
     write(*,*) 'transforming backward...'
	 call fftw_execute_dft_c2r(plan_cr,cfx,fcr)
     
     
     open(33,file="fx_back.d",status='replace')
     do i = 1,N
       t  = 4._pr*real(i-1,pr)/real(N-1,pr)
       write(33,*) t,fcr(i)/real(N,pr)
     enddo
     close(33)
     
     
     call fftw_destroy_plan( plan_rc )
     call fftw_destroy_plan( plan_cr )
     end program 
     
