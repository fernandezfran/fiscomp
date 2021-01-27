program comp_met
  !
  !programa para comparar los métodos de integración numérica de trapecio, simpson
  !y gauss a través de sus errores relativos a la solución exacta de una integral
  !conocida.
  !
  !En el modulo precision se encuentra definida la variable dp (doble precisión) y
  !en el modulo integrales están desarrollados los distintos métodos
  !
  use precision, only: pr=>dp
  use integrales
  implicit none
  integer            :: n
  real(pr)           :: I_e, I_t, I_s, I_g, a, b, h

  I_e = 1._pr - exp(-1._pr)

  a = 0._pr
  b = 1._pr

  open(49, file="errores.dat", status="replace")
  write(49,*) "# n_pts, |I_trap - I_e|/I_e, |I_simp - I_e|/I_e, |I_gauss - I_e|/I_e"

  do n = 3, 2003, 2
      h = (b - a) / real(n - 1, pr)

      call trap(n,a,h,I_t)
      call simp(n,a,h,I_s)
      call gausscuad(n,a,b,I_g)

      write(49,"(I5,2x,3(E15.6,2x))") n, abs(I_t - I_e)/I_e, abs(I_s - I_e)/I_e, &
                                      & abs(I_g - I_e)/I_e
  enddo

  close(49)

end program comp_met
