module randomnum
!
!modulo con los distintos generadores de números aleatorios
!
!suponiendo que defino
!real(pr)       :: rand
!integer        :: seed
!entonces a las siguientes funciones (o subrutinas) se las llama así:
!
!línea 45: ran0               :: rand = ran0(seed)
!línea 65: ran2               :: rand = ran2(seed)
!línea 109: mzran             :: rand = rmzran()        !semilla default
!línea 260: Mersenne-Twister  :: rand = grnd()          !          ||
!                                          ^la `d' es de doble presicion
!
use precision, only: pr=>dp

! MZRAN
INTEGER,PARAMETER :: K4B=selected_int_kind(9)
INTEGER,PARAMETER :: DP =KIND(1.0D0)

! MERSENNE - TWISTER
! Default seed
    integer, parameter :: defaultsd = 4357
! Period parameters
    integer, parameter :: N = 624, N1 = N + 1

! the array for the state vector
    integer, save, dimension(0:N-1) :: mt
    integer, save                   :: mti = N1
    
! Overload procedures for saving and getting mt state
    interface mtsave
      module procedure mtsavef
      module procedure mtsaveu
    end interface
    interface mtget
      module procedure mtgetf
      module procedure mtgetu
    end interface


contains

function ran0(idum)
!
!ran0 escrito por mí
!
implicit none
real(pr)               :: ran0
integer, intent(inout) :: idum
integer                :: ia, im, iq, ir, k
real(pr)               :: am
parameter(ia=16807, im=2147483647, am=1._pr/2147483647._pr, iq=127773, ir=2836)

k = idum/iq
idum = ia*(idum - k*iq) - ir*k
if (idum <= 0) idum = idum + im
ran0 = am*real(idum,pr)

end function ran0


!de acá en adelante --> NUMERICAL RECIPIES
      function ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL(pr) ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
      IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, &
      NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      end function ran2




!***************************************************************************C
!***************************************************************************C
!!
!*********************************************************************
!**                Random Number Generator  (period > 2^94)         **
!**  See  G. Marsaglia and A. Zaman, Comp. Phys. 8, 117-121 (1994)  **
!*********************************************************************
        function rmzran()
        implicit none
        integer (k4b) ::  mzran,n,i,j,k,is,js,ks,ns,mzranset
        real (DP)     ::  rmzran
        save i,j,k,n
        data i,j,k,n/521288629,362436069,16163801,1131199299/
        mzran = i-k
        if(mzran.lt.0) mzran = mzran + 2147483579
        i = j
        j = k
        k = mzran
        n = 69069*n + 1013904243
!        n = ishft(3533*ishft(n,-16)+iand(n,65535),16)+3533*iand(n,65535))
!        n = n + 1013904243
        mzran = mzran + n
!       For random reals on  (0,1): UNI() = .5+.2328306E-9*mzran()
!       For random reals on (-1,1): VNI() = .4656613E-9*mzran()
!
        rmzran = real(mzran,kind=DP)*0.2328306E-9_DP + 0.5_DP 
!
        return
        entry mzranset(is,js,ks,ns)
        i = 1+iabs(is)
        j = 1+iabs(js)
        k = 1+iabs(ks)
        n = ns
        mzranset = n
        return
        end function rmzran
!
!*********************************************************************
!*********************************************************************
!
!
!
!*********************************************************************
!**     subrutina de para asegurar que los enteros sean de          **
!**     bits y cumplan los requisitos necesarios para el generador. **
!**     
!**       Extraida del Numerical Recipes                            ** 
!*********************************************************************
!
         SUBROUTINE mzran_init(is,js,ks,n)
           IMPLICIT NONE
           INTEGER(K4B), intent(in), optional :: n,is,js,ks
           INTEGER(K4B),PARAMETER ::hg=huge(1_K4B),hgm=-hg,hgng=hgm-1
           INTEGER(K4B)::new,j,hgt,nl
             hgt=hg
!           The following lines check that kind value K4B is in fact a 
!           32-bit integer with the usual properties that we expect it 
!           to have (under negation and wrap-around addition).If all of 
!           these tests are satisfied,then the routines that use this 
!           module are portable,even though they go eyond Fortran 90 's 
!           integer model.
           if (hg /=2147483647)call nrerror('ran_init:arith assump 1 fails ')
           if (hgng >=0)call nrerror('ran_init:arith assump 2 fails ')
           if (hgt+1 /=hgng)call nrerror('ran_init:arith assump 3 fails ')
           if (not(hg)>=0)call nrerror('ran_init:arith assump 4 fails ')
           if (not(hgng)<0)call nrerror('ran_init:arith assump 5 fails ')
           if (hg+hgng >=0)call nrerror('ran_init:arith assump 6 fails ')
           if (not(-1_k4b)<0)call nrerror('ran_init:arith assump 7 fails ')
           if (not(0_k4b)>=0)call nrerror('ran_init:arith assump 8 fails ')
           if (not(1_k4b)>=0)call nrerror('ran_init:arith assump 9 fails ')
!
!          si se se quiere iniciar la secuencia de mzran con otras semaillas
!          (no usando las definidas el la funcion mzran)
!
           if (present(is) .and. present(js) .and. present(ks) .and. present(n)) then
             nl = mzranset(is,js,ks,n)
           endif
!
           return
         end subroutine mzran_init
!
!
         SUBROUTINE nrerror(char)
           character, intent(in) :: char
           write(*,*) char
           stop
         end subroutine nrerror



!____________________________________________________________________________
! A C-program for MT19937: Real number version
!   genrand() generates one pseudorandom real number (double)
! which is uniformly distributed on [0,1]-interval, for each
! call. sgenrand(seed) set initial values to the working area
! of 624 words. Before genrand(), sgenrand(seed) must be
! called once. (seed is any 32-bit integer except for 0).
! Integer generator is obtained by modifying two lines.
!   Coded by Takuji Nishimura, considering the suggestions by
! Topher Cooper and Marc Rieffel in July-Aug. 1997.
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Library General Public
! License as published by the Free Software Foundation; either
! version 2 of the License, or (at your option) any later
! version.
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
! See the GNU Library General Public License for more details.
! You should have received a copy of the GNU Library General
! Public License along with this library; if not, write to the
! Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
! 02111-1307  USA
!
! Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.
! When you use this, send an email to: matumoto@math.keio.ac.jp
! with an appropriate reference to your work.
!
!***********************************************************************
! Fortran translation by Hiroshi Takano.  Jan. 13, 1999.
!
!   genrand()      -> double precision function grnd()
!   sgenrand(seed) -> subroutine sgrnd(seed)
!                     integer seed
!
! This program uses the following non-standard intrinsics.
!   ishft(i,n): If n>0, shifts bits in i by n positions to left.
!               If n<0, shifts bits in i by n positions to right.
!   iand (i,j): Performs logical AND on corresponding bits of i and j.
!   ior  (i,j): Performs inclusive OR on corresponding bits of i and j.
!   ieor (i,j): Performs exclusive OR on corresponding bits of i and j.
!
!***********************************************************************
! Fortran version rewritten as an F90 module and mt state saving and getting
! subroutines added by Richard Woloshyn. (rwww@triumf.ca). June 30, 1999


!Initialization subroutine
  subroutine sgrnd(seed)
    implicit none
!
!      setting initial seeds to mt[N] using
!      the generator Line 25 of Table 1 in
!      [KNUTH 1981, The Art of Computer Programming
!         Vol. 2 (2nd Ed.), pp102]
!
    integer, intent(in) :: seed

    mt(0) = iand(seed,-1)
    do mti=1,N-1
      mt(mti) = iand(69069 * mt(mti-1),-1)
    enddo
!
    return
  end subroutine sgrnd

!Random number generator
  real(8) function grnd()
    implicit integer(a-z)

! Period parameters
    integer, parameter :: M = 397, MATA  = -1727483681
!                                    constant vector a
    integer, parameter :: LMASK =  2147483647
!                                    least significant r bits
    integer, parameter :: UMASK = -LMASK - 1
!                                    most significant w-r bits
! Tempering parameters
    integer, parameter :: TMASKB= -1658038656, TMASKC= -272236544

    dimension mag01(0:1)
    data mag01/0, MATA/
    save mag01
!                        mag01(x) = x * MATA for x=0,1

    TSHFTU(y)=ishft(y,-11)
    TSHFTS(y)=ishft(y,7)
    TSHFTT(y)=ishft(y,15)
    TSHFTL(y)=ishft(y,-18)

    if(mti.ge.N) then
!                       generate N words at one time
      if(mti.eq.N+1) then
!                            if sgrnd() has not been called,
        call sgrnd( defaultsd )
!                              a default initial seed is used
      endif

      do kk=0,N-M-1
          y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
          mt(kk)=ieor(ieor(mt(kk+M),ishft(y,-1)),mag01(iand(y,1)))
      enddo
      do kk=N-M,N-2
          y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
          mt(kk)=ieor(ieor(mt(kk+(M-N)),ishft(y,-1)),mag01(iand(y,1)))
      enddo
      y=ior(iand(mt(N-1),UMASK),iand(mt(0),LMASK))
      mt(N-1)=ieor(ieor(mt(M-1),ishft(y,-1)),mag01(iand(y,1)))
      mti = 0
    endif

    y=mt(mti)
    mti = mti + 1 
    y=ieor(y,TSHFTU(y))
    y=ieor(y,iand(TSHFTS(y),TMASKB))
    y=ieor(y,iand(TSHFTT(y),TMASKC))
    y=ieor(y,TSHFTL(y))

    if(y .lt. 0) then
      grnd=(dble(y)+2.0d0**32)/(2.0d0**32-1.0d0)
    else
      grnd=dble(y)/(2.0d0**32-1.0d0)
    endif

    return
  end function grnd

!State saving subroutines.
! Usage:  call mtsave( file_name, format_character )
!    or   call mtsave( unit_number, format_character )
! where   format_character = 'u' or 'U' will save in unformatted form, otherwise
!         state information will be written in formatted form.
  subroutine mtsavef( fname, forma )

!NOTE: This subroutine APPENDS to the end of the file "fname".

    character(*), intent(in) :: fname
    character, intent(in)    :: forma

    select case (forma)
      case('u','U')
       open(unit=10,file=trim(fname),status='UNKNOWN',form='UNFORMATTED', &
            position='APPEND')
       write(10)mti
       write(10)mt

      case default
       open(unit=10,file=trim(fname),status='UNKNOWN',form='FORMATTED', &
            position='APPEND')
       write(10,*)mti
       write(10,*)mt

    end select
    close(10)

    return
  end subroutine mtsavef

  subroutine mtsaveu( unum, forma )

    integer, intent(in)    :: unum
    character, intent(in)  :: forma

    select case (forma)
      case('u','U')
       write(unum)mti
       write(unum)mt

      case default
       write(unum,*)mti
       write(unum,*)mt

      end select

    return
  end subroutine mtsaveu

!State getting subroutines.
! Usage:  call mtget( file_name, format_character )
!    or   call mtget( unit_number, format_character )
! where   format_character = 'u' or 'U' will read in unformatted form, otherwise
!         state information will be read in formatted form.
  subroutine mtgetf( fname, forma )

    character(*), intent(in) :: fname
    character, intent(in)    :: forma

    select case (forma)
      case('u','U')
       open(unit=10,file=trim(fname),status='OLD',form='UNFORMATTED')
       read(10)mti
       read(10)mt

      case default
       open(unit=10,file=trim(fname),status='OLD',form='FORMATTED')
       read(10,*)mti
       read(10,*)mt

    end select
    close(10)

    return
  end subroutine mtgetf

  subroutine mtgetu( unum, forma )

    integer, intent(in)    :: unum
    character, intent(in)  :: forma

    select case (forma)
      case('u','U')
       read(unum)mti
       read(unum)mt

      case default
       read(unum,*)mti
       read(unum,*)mt

      end select

    return
  end subroutine mtgetu


end module
