program main
  use volcgases, only: solve_gases, dp
  implicit none
  call test_volcgases()
contains
  subroutine test_volcgases()
    real(dp) :: T
    real(dp) :: P
    real(dp) :: f_O2
    real(dp) :: mCO2tot
    real(dp) :: mH2Otot
    real(dp) :: P_H2O
    real(dp) :: P_H2
    real(dp) :: P_CO2
    real(dp) :: P_CO
    real(dp) :: P_CH4
    real(dp) :: alphaG 
    real(dp) :: x_CO2
    real(dp) :: x_H2O
    integer :: ierr

    real(dp) :: tn(10)
    integer, parameter :: nt = 100000
    integer :: i

    T = 1473.0_dp
    P = 1.0_dp
    f_O2 = oxygen_fugacity(0.0_dp, T, P)
    mCO2tot = 1000e-6_dp
    mH2Otot = 1000e-6_dp
    
    call cpu_time(tn(1))
    do i = 1,nt
      call solve_gases(T, P, f_O2, mCO2tot, mH2Otot, &
                       P_H2O, P_H2, P_CO2, P_CO, P_CH4, alphaG, x_CO2, x_H2O, &
                       ierr)
      if (ierr /= 0) then
        print*,ierr
        stop 1
      endif
    enddo
    call cpu_time(tn(2))

    print"(a15,es10.2,a)","time = ",(tn(2) - tn(1))/real(nt,dp),' seconds'
    print"(a15,es10.2)","H2O = ",P_H2O/P
    print"(a15,es10.2)","H2 = ",P_H2/P
    print"(a15,es10.2)","CO2 = ",P_CO2/P
    print"(a15,es10.2)","CO = ",P_CO/P
    print"(a15,es10.2)","CH4 = ",P_CH4/P
    print"(a15,es10.2)","alphaG = ",alphaG

    call assert(P_H2O, 0.5252333529267064_dp, 1.0e-10_dp, 'P_H2O')
    call assert(P_H2, 0.01183114571626002_dp, 1.0e-10_dp, 'P_H2')
    call assert(P_CO2, 0.4386535437794375_dp, 1.0e-10_dp, 'P_CO2')
    call assert(P_CO, 0.024281957577105908_dp, 1.0e-10_dp, 'P_CO')
    call assert(P_CH4, 4.903435583386011e-13_dp, 1.0e-10_dp, 'P_CH4')
    call assert(alphaG, 0.00316539763088389_dp, 1.0e-10_dp, 'alphaG')
    call assert(x_CO2, 4.2433333552036306e-07_dp, 1.0e-10_dp, 'x_CO2')
    call assert(x_H2O, 0.0018867859399539038_dp, 1.0e-10_dp, 'x_H2O')
    
  end subroutine

  function oxygen_fugacity(DFMQ, T, P) result(f_O2)
    real(dp), intent(in) :: DFMQ, T, P
    real(dp) :: f_O2
    real(dp), parameter :: A = 25738.0_dp 
    real(dp), parameter :: B = 9.0_dp
    real(dp), parameter :: C = 0.092_dp
    real(dp) :: log10_FMQ
    log10_FMQ = (-A/T+B+C*(P-1.0_dp)/T)
    f_O2 = 10.0_dp**(log10_FMQ + DFMQ)
  end function

  subroutine assert(a, b, tol, msg)
    real(dp), intent(in) :: a, b, tol
    character(*), intent(in) :: msg

    character(:), allocatable :: err_msg
    character(22) :: a_str, b_str

    if (.not.is_close(a, b, tol)) then
      write(a_str,'(es22.15)') a
      write(b_str,'(es22.15)') b
      err_msg = '"'//msg//'" assert failed: a = '//a_str//' and b = '//b_str
      error stop err_msg
    endif
  end subroutine

  !> coppied from fortran stdlib v0.2.0
  elemental function is_close(a, b, tol, abs_tol, equal_nan) result(close)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    real(dp), intent(in) :: a, b
    real(dp), intent(in), optional :: tol, abs_tol
    logical, intent(in), optional :: equal_nan
    logical :: close

    real(dp) :: rel_tol_, abs_tol_
    logical :: equal_nan_

    if (present(tol)) then
      rel_tol_ = tol
    else
      rel_tol_ = 1.0e-5_dp
    endif

    if (present(abs_tol)) then
      abs_tol_ = abs_tol
    else
      abs_tol_ = 0.0_dp
    endif

    if (present(equal_nan)) then
      equal_nan_ = equal_nan
    else
      equal_nan_ = .false.
    endif

    if (ieee_is_nan(a) .or. ieee_is_nan(b)) then
        close = merge(.true., .false., equal_nan_ .and. ieee_is_nan(a) .and. ieee_is_nan(b))
    else
        close = abs(a - b) <= max(abs(rel_tol_*max(abs(a), abs(b))), abs(abs_tol_))
    end if     

  end function
end program