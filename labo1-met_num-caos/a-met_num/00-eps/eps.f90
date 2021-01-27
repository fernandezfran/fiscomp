program eps
      !El epsilon de la maquina (notebook) con doble presicion da
      !2.2204460492503131E-016
      implicit none
      integer, parameter :: pr = selected_real_kind(6)
      real(pr)           :: ep, one
      integer            :: i

      ep = 1._pr
      do i = 1, 150
          ep = ep / 2._pr
          one = 1._pr + ep
          !write(*,*) one
          if (one == 1._pr) then
              write(*,*) "el epsilon de la maquina es", 2._pr * ep
              exit
          endif
      enddo

      write(*,*) "el epsilon (segun la maquina) es:", epsilon(ep)

end program eps
