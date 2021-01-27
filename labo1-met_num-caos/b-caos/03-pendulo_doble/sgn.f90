program signo
implicit none
integer, parameter :: pr = selected_real_kind(13)
integer            :: sgn
real(pr)           :: a, b

write(*,*) 'de un valor a > 0'
read(*,*) a

write(*,*) 'de un valor b < 0'
read(*,*) b


sgn = sign(1,1)
write(*,*) 'sgn = sign(1,1) = ', sgn
sgn = sign(1._pr,a)
write(*,*) 'sgn = sign(1,a) = ', sgn
sgn = sign(1._pr,b)
write(*,*) 'sgn = sign(1,b) = ', sgn

end program signo
