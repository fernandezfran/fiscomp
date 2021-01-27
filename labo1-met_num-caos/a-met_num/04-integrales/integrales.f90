module integrales
  !
  !modulo en el que se encuentran la regla del trapecio, de simpson
  !y la cuadratura de gauss-legendre que está en la pagina de la materia fiscomp
  !y que utiliza gaussmod.f90 (también proporcionado por la materia)
  !
  use precision, only: pr=>dp
  CONTAINS

  subroutine trap(nn,aa,hh,It)
    !
    !como w=h salvo en los extremos, saco estos puntos fuera del loop
    !
    implicit none
    integer, intent(in)   :: nn
    real(pr), intent(in)  :: aa, hh
    real(pr), intent(out) :: It
    integer               :: i
    real(pr)              :: x_i, h, w, fx, fa, fb

    It = 0._pr
    w = hh
    do i = 2, nn - 1
        x_i = aa + real(i - 1, pr) * hh
        call func(x_i,fx)
        It = It + w * fx
    enddo
   
    call func(aa, fa) 
    call func(aa + real(nn - 1, pr) * hh, fb)

    It = It + 0.5_pr * hh * (fa + fb)

    return
  end subroutine trap


  subroutine simp(nn,aa,hh,Is)
    !
    !separo las contribuciones pares de las impartes porque tienen distintso pesos,
    !multiplico por los mismos al final del loop y agrego un punto que falta (ver
    !definición del Landau) y los extremos.
    !
    implicit none
    integer, intent(in)   :: nn
    real(pr), intent(in)  :: aa, hh
    real(pr), intent(out) :: Is
    integer               :: i
    real(pr)              :: x_i, h, w, fx, fa, fb, Is_par, Is_impar

    if (mod(nn,2) .eq. 0) print *,"WARNING: hay problemas con la paridad &
                                  & del n en la regla de Simpson..."

    Is = 0._pr; Is_par = 0._pr; Is_impar = 0._pr
    do i = 2, nn - 3, 2
      x_i = aa + real(i -1, pr) * hh  
      call func(x_i,fx)
      Is_par = Is_par + fx
      call func(x_i + hh, fx)
      Is_impar = Is_impar + fx
    enddo
    call func(aa + real(nn - 2 , pr) * hh ,fx)
    Is_par = Is_par + fx
    
    call func(aa, fa)
    call func(aa + real(nn - 1, pr) * hh, fb)

    Is = hh * (fa + 4._pr * Is_par + 2._pr * Is_impar + fb)

    Is = Is/3._pr

    return
  end subroutine simp


  subroutine gausscuad(n,a,b,gau)
    !calcula la integral de la funcion Func, entre a y b
    !utilizando n puntos de evaluacion, mediante
    !la cuadratura de Gauss-Legendre (gau)

    use gaussmod, only : gauss
    implicit none
    real (pr), INTENT(IN)  :: a,b
    integer, INTENT(IN)            :: n
    real (pr), INTENT(OUT) :: gau
    real (pr)              :: w(0:n),x(0:n),fx
    integer                        :: j

    call gauss(n,0,a,b,x(0:n),w(0:n))  ! calculo nodos y pesos para Gauss

  ! calculo la integral
      gau  = 0._pr

      do j=0,n
        call func(x(j),fx)
        gau = gau + fx*w(j)
      enddo
    return
  end subroutine gausscuad


  subroutine func(x,fx)
    !
    !función que se desea integrar
    !
    implicit none
    real (pr), INTENT(IN)  :: x
    real (pr), INTENT(OUT) :: fx

    fx = exp(-x)
    return
    end subroutine func

end module integrales
