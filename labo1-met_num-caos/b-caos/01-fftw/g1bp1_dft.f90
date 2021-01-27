     program g1bp1
     use precision, only: pr => dp
     use dftmod
     implicit none
     real(pr), parameter                       :: pi=4._pr*atan(1._pr)
     real(pr), dimension(:), allocatable       :: fxr,fxi,tfxr,tfxi
     complex(kind((1._pr,1._pr))), allocatable :: tfx(:),cfx(:)
     integer                                   :: N,i,j
     real(pr)                                  :: t,factor
     
     write(*,*) "pi = ", pi
     
     write(*,*) 'numero de puntos para la FFT?'
     read(*,*)  N
     
     allocate( fxr(N), fxi(N),tfxr(N),tfxi(N) )
          
     
     fxi(:) = 0._pr
     open(33,file="fx_dft.d",status='replace')
     do i = 1,N
       t  = 4._pr*real(i-1,pr)/real(N-1,pr)
       fxr(i) = 2._pr*sin(pi*0.5_pr*t) + cos(20._pr*pi*t)
       write(33,*) t,fxr(i),fxi(i)
     enddo
     close(33)
     
     
     ! hacemos la transformada 
     call dft(fxr,fxi,tfxr,tfxi,N)
     
     open(33,file="tfx_dft.d",status='replace')
     factor = 1._pr/4._pr 
     !   notar que    tfx(k) = tfx(N+k)
     do i = 1,N/2
         write(33,*) factor*real(i-1-N/2,pr), tfxr(i+N/2)/real(N,pr), &
                     & tfxi(i+N/2)/real(N,pr)
     enddo
     do i = N/2, N-1
         write(33,*) factor*real(i-N/2,pr), tfxr(i-N/2+1)/real(N,pr), &
                     & tfxi(i-N/2+1)/real(N,pr)
     enddo
     close(33)
     
   
     end program 
     
