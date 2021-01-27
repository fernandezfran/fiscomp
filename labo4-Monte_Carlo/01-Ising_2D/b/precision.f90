module precision
      integer, parameter  :: sp=selected_real_kind(6),dp=selected_real_kind(13)
      integer, parameter  :: qp=selected_real_kind(21)

      integer, parameter  :: k18 = selected_int_kind(18) ! (10^-18, 10^18)
end module precision
