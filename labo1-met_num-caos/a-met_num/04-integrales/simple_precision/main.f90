program comp_met
        !programa para comparar los métodos de integración numérica de
        !trapecio, simpson y gauss
        use precision, only: pr=>sp
        use integrales
        implicit none
        integer            :: n
        real(pr)           :: I_e, I_t, I_s, a, b, h

        I_e = 1._pr - exp(-1._pr)

        a = 0._pr
        b = 1._pr

        open(49, file="errores.dat", status="replace")
        write(49,*) "# n_pts, |I_trap - I_e|/I_e, |I_simp - I_e|/I_e"
        
        do n = 3, 9003, 2
            h = (b - a) / real(n - 1, pr)
            call trap(n,a,h,I_t)

            call simp(n,a,h,I_s)
            
            write(49,"(I5,2x,2(E15.6,2x))") n, abs(I_t - I_e)/I_e, &
                                         abs(I_s - I_e)/I_e
        enddo

        close(49)

end program comp_met
