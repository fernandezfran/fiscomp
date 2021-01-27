program intnum
        implicit none
        integer, parameter :: pr = selected_real_kind(13)
        integer            :: k, n
        real(pr)           :: I_e, I_t, I_s, a, b, h

        I_e = exp(1._pr) - 1._pr
    
        a = 0._pr
        b = 1._pr

        open(49, file="intnum.dat", status="replace") 
        write(49,*) "# n_pts, |I_trap - I_e|, |I_simp - I_e|" 

        do k = 2, 10
            n = 2**k + 1
            h = (b - a) / real(n - 1, pr)
            
            call trap(n,a,h,I_t)
            call simp(n,a,h,I_s)

            write(49,"(I4,x,2(E15.6,x))") n, abs(I_t - I_e), abs(I_s - I_e)
        enddo

CONTAINS

        subroutine trap(nn,aa,hh,It)
                implicit none
                integer, intent(in)   :: nn
                real(pr), intent(in)  :: aa, hh
                real(pr), intent(out) :: It
                integer               :: i
                real(pr)              :: x_i, h, w
                
                It = 0._pr
                do i = 1, nn
                    x_i = aa + real(i - 1, pr) * hh
                    if (i .eq. 1 .or. i .eq. nn) then
                         w = hh/2._pr
                    else
                         w = hh
                    endif
                    It = It + w * f(x_i)
                enddo

                return
        end subroutine trap

        
        subroutine simp(nn,aa,hh,Is)
                implicit none
                integer, intent(in)   :: nn
                real(pr), intent(in)  :: aa, hh
                real(pr), intent(out) :: Is
                integer               :: i
                real(pr)              :: x_i, h, w

                Is = 0._pr
                do i = 1, nn
                    x_i = aa + real(i - 1, pr) * hh
                    if (i .eq. 1 .or. i .eq. nn) then
                        w = hh / 3._pr
                    else if (mod(i,2) .eq. 0) then
                        w = 4._pr * hh / 3._pr
                    else if (mod(i+1,2) .eq. 0) then
                        w = 2._pr * hh / 3._pr
                    endif
                    Is = Is + w * f(x_i)
                enddo

                return
        end subroutine simp


        function f(x)
                implicit none
                real (pr), intent(in)  :: x
                real (pr)              :: f
        
                f = exp(x)
        
        end function f

end program intnum
