module cell
use precision, only : pr=>dp
use initialization, only : L, N, rcut 
implicit none !              , number of cells
integer                  :: M, nocells
real(pr)                 :: cell_size, cell_size_inv
integer, allocatable     :: hochain(:), cell_list(:), map(:)
!                        head of chain, linked list

contains


subroutine maps
!enumeración de las celdas en el espacio grillado 3dim y mapa de celdas vecinas
! incluyendo pbc
implicit none
integer                  :: ixc, iyc, izc !indices de las celdas
integer                  :: imap          !numero de celda

! definición de variables y asignación de memoria
M = int(L/rcut)                  !cantidad de celdas de largo mínimo rcut
cell_size = L/real(M,pr)           !largo de cada una de las celdas
cell_size_inv = 1._pr/cell_size
nocells   = M*M*M                  !numero total de celdas
!write(*,*) 'cantidad de celdas en una dirección', M
!write(*,*) 'largo de cada celda', cell_size
!write(*,*) 'cantidad total de celdas', nocells
allocate( hochain(nocells), cell_list(N), map(13*nocells) ) 

do ixc = 1, M
    do iyc = 1, M
        do izc = 1, M
                 !si fuera MC este 13 sería 26   
            imap = 13*(icell(ixc, iyc, izc) - 1)
            
            map(imap + 1 ) = icell(ixc + 1, iyc    , izc    )
            map(imap + 2 ) = icell(ixc + 1, iyc + 1, izc    )
            map(imap + 3 ) = icell(ixc    , iyc + 1, izc    )
            map(imap + 4 ) = icell(ixc - 1, iyc + 1, izc    )
            map(imap + 5 ) = icell(ixc + 1, iyc    , izc - 1)
            map(imap + 6 ) = icell(ixc + 1, iyc + 1, izc - 1)
            map(imap + 7 ) = icell(ixc    , iyc + 1, izc - 1)
            map(imap + 8 ) = icell(ixc - 1, iyc + 1, izc - 1)
            map(imap + 9 ) = icell(ixc + 1, iyc    , izc + 1)
            map(imap + 10) = icell(ixc + 1, iyc + 1, izc + 1)
            map(imap + 11) = icell(ixc    , iyc + 1, izc + 1)
            map(imap + 12) = icell(ixc - 1, iyc + 1, izc + 1)
            map(imap + 13) = icell(ixc    , iyc    , izc + 1)

        enddo
    enddo
enddo

    contains !función que solo va a usar maps

        function icell(a, b, c)
        integer, intent(in)  :: a, b, c
        integer              :: icell
    
            icell = 1 + mod(a - 1 + M, M) + mod(b - 1 + M, M)*M &
                      + mod(c - 1 + M, M)*M*M

        end function icell

end subroutine maps


subroutine links
! llena los arreglos cell_list y hochain
! debe llamarse luego de actualizar las posiciones de las partículas y antes
! de calcular las fuerzas
use initialization, only : x, y, z
implicit none
integer            :: i, cellidx, cellx, celly, cellz, j

hochain(:) = 0
do i = 1, N
    cellx = floor(x(i)*cell_size_inv)
    celly = floor(y(i)*cell_size_inv)
    cellz = floor(z(i)*cell_size_inv)
    cellidx = 1 + cellx + celly*M + cellz*M*M

    cell_list(i) = hochain(cellidx)
    hochain(cellidx) = i
enddo

end subroutine links


end module
