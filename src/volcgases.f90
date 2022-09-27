module volcgases
  use volcgases_const, only: dp
  implicit none
  private

  public :: solve_gases, dp
  public :: solve_gases_c

  type :: SolveGasesData
    real(dp), pointer :: P
    real(dp), pointer :: F1, F2
    real(dp), pointer :: xCO2tot, xH2Otot
    real(dp), pointer :: C1, C2, C3
  end type

  type :: LMData
    integer :: m
    integer :: n
    real(dp), allocatable :: x(:) ! (n) inout
    real(dp), allocatable :: fvec(:) ! (m) out
    real(dp) :: ftol
    real(dp) :: xtol
    real(dp) :: gtol
    integer :: maxfev
    real(dp) :: epsfcn
    real(dp), allocatable :: diag(:) ! (n) inout
    integer :: mode
    real(dp) :: factor
    integer :: nprint
    integer :: info ! out
    integer :: nfev ! out
    integer :: njev
    real(dp), allocatable :: fjac(:,:) ! (ldfjac,n) out
    integer :: ldfjac
    integer, allocatable :: ipvt(:) ! (n) out
    real(dp), allocatable :: qtf(:) ! (n) out
    real(dp), allocatable :: wa1(:) ! (n) ! inout
    real(dp), allocatable :: wa2(:) ! (n) ! inout
    real(dp), allocatable :: wa3(:) ! (n) ! inout
    real(dp), allocatable :: wa4(:) ! (m) ! inout
  end type
  interface LMData
    procedure :: LMData_create
  end interface

contains

  subroutine solve_gases(T, P, f_O2, mCO2tot, mH2Otot, &
                         P_H2O, P_H2, P_CO2, P_CO, P_CH4, alphaG, x_CO2, x_H2O, &
                         ierr)
    use volcgases_const, only: A1, A2, M_H2O, M_CO2, C_CO2, C_H2O, x, a_CO2, a_H2O
    real(dp), intent(in) :: T !! Temperature of degassing (K)
    real(dp), target, intent(in) :: P !! Pressure of degassing in (bars)
    real(dp), intent(in) :: f_O2 !! Oxygen fugacity (bar)
    real(dp), intent(in) :: mCO2tot !! Mass fraction of CO2 in melt
    real(dp), intent(in) :: mH2Otot !! Mass fraction of H2O in melt
    real(dp), intent(out) :: P_H2O !! Partial pressure of H2O in gas (bar)
    real(dp), intent(out) :: P_H2 !! Partial pressure of H2 in gas (bar)
    real(dp), intent(out) :: P_CO2 !! Partial pressure of CO2 in gas (bar)
    real(dp), intent(out) :: P_CO !! Partial pressure of CO2 in gas (bar)
    real(dp), intent(out) :: P_CH4 !! Partial pressure of CO2 in gas (bar)
    real(dp), intent(out) :: alphaG !! (mol gas)/(mol gas and magma)
    real(dp), intent(out) :: x_CO2 !! mol fraction of CO2 in magma after degassing
    real(dp), intent(out) :: x_H2O !! mol fraction of H2O in magma after degassing
    integer, intent(out) :: ierr !! Indicates error.
                                 !!
                                 !! * ierr == 0, means success
                                 !! * ierr == -1, means improper inputs (unphysical negative value)
                                 !! * ierr == -2, means the first, simple, root solve failed

    type(SolveGasesData) :: d
    real(dp), target :: F1, F2
    real(dp), target :: xCO2tot, xH2Otot
    real(dp) :: K1, K2, K3
    real(dp), target :: C1, C2, C3
    real(dp) :: found_high, found_low
    integer :: i

    ierr = 0

    ! check for valid inputs
    if (T <= 0.0_dp) then
      ierr = -1
      return
    endif
    if (P <= 0.0_dp) then
      ierr = -1
      return
    endif
    if (f_O2 <= 0.0_dp) then
      ierr = -1
      return
    endif
    if (mCO2tot <= 0.0_dp) then
      ierr = -1
      return
    endif
    if (mH2Otot <= 0.0_dp) then
      ierr = -1
      return
    endif

    F1 = log(1.0_dp/(M_CO2*x*1.0e6_dp))+C_CO2*P/T+A1
    F2 = log(1.0_dp/(M_H2O*x*100.0_dp))+C_H2O*P/T+A2

    ! calculate mol fraction of CO2 and H2O in the magma
    xCO2tot = (mCO2tot/M_CO2)/x
    xH2Otot = (mH2Otot/M_H2O)/x

    ! equilibrium constants
    ! made with Nasa thermodynamic database (Burcat database)
    K1 = exp(-29755.11319228574_dp/T+6.652127716162998_dp)
    K2 = exp(-33979.12369002451_dp/T+10.418882755464773_dp)
    K3 = exp(-96444.47151911151_dp/T+0.22260815074146403_dp)

    ! constants
    C1 = K1/f_O2**0.5_dp
    C2 = K2/f_O2**0.5_dp
    C3 = K3/f_O2**2.0_dp

    ! pointers to data that needs to be passed around
    d%P => P
    d%F1 => F1
    d%F2 => F2
    d%xCO2tot => xCO2tot
    d%xH2Otot => xH2Otot
    d%C1 => C1
    d%C2 => C2
    d%C3 => C3

    ! Now solve simple system.
    ! This block finds a good starting guess
    block
      real(dp) :: find1, find2, log10P
      real(dp) :: fvec(1)
      integer :: iflag

      log10P = log10(P)
      found_high = 0.0_dp
      do i = 1,1000
        find1 = single_linspace(log10P, -15.0_dp, 1000, i)
        find1 = 10.0_dp**find1
        call fcn_s_(1, 1, [find1], fvec, iflag)
        if (fvec(1) == fvec(1)) then
          found_high = find1
          exit
        endif
      enddo

      found_low = 0.0_dp
      do i = 1,1000
        find2 = single_linspace(-15.0_dp, log10P, 1000, i)
        find2 = 10.0_dp**find2
        call fcn_s_(1, 1, [find2], fvec, iflag)
        if (fvec(1) == fvec(1)) then
          found_low = find2
          exit
        endif
      enddo
    end block

    ! solve the simple system of equations
    block 
      use minpack_module, only: lmdif
      type(LMData) :: ld

      ld = LMData(m=1, n=1, ftol=1.49012e-8_dp, xtol=1.49012e-8_dp, gtol=0.0_dp, maxfev=10000)
      ld%x(1) = found_high
      call lmdif(fcn_s_, ld%m, ld%n, ld%x, ld%fvec, ld%ftol, ld%xtol, ld%gtol, ld%maxfev, ld%epsfcn, ld%diag, &
                 ld%mode, ld%factor, ld%nprint, ld%info, ld%nfev, ld%fjac, ld%ldfjac, ld%ipvt, &
                 ld%qtf, ld%wa1, ld%wa2, ld%wa3, ld%wa4)
      if (ld%info < 1 .or. ld%info > 4) then
        ld%x(1) = found_low
        call lmdif(fcn_s_, ld%m, ld%n, ld%x, ld%fvec, ld%ftol, ld%xtol, ld%gtol, ld%maxfev, ld%epsfcn, ld%diag, &
                   ld%mode, ld%factor, ld%nprint, ld%info, ld%nfev, ld%fjac, ld%ldfjac, ld%ipvt, &
                   ld%qtf, ld%wa1, ld%wa2, ld%wa3, ld%wa4)
        if (ld%info < 1 .or. ld%info > 4) then
          ! Solve for simple system failed. 
          ! We give up
          ierr = -2
          return
        endif
      endif
      P_CO2 = ld%x(1)
    end block

    ! use solution of simple system as a guess for more complicated one
    P_CO = C2*P_CO2
    x_CO2 = exp(F1+a_CO2*log(P_CO2))
    alphaG = P*(x_CO2-xCO2tot)/(-P_CO2+P*x_CO2)
    P_H2O = (P-P_CO2-C2*P_CO2)/(1.0_dp+C1+C3*P_CO2)
    P_H2 = C1*P_H2O
    P_CH4 = C3*P_CO2*P_H2O**2.0_dp
    x_H2O = exp(F2+a_H2O*log(P_H2O))
    ! use different alphaG as inital guess
    alphaG = 0.1_dp

    ! Now do more complicated root solve
    block
      use minpack_module, only: lmder
      type(LMData) :: ld
      real(dp) :: error
      integer :: iflag

      ld = LMData(m=8, n=8, ftol=1.49012e-8_dp, xtol=1.49012e-8_dp, gtol=0.0_dp, maxfev=10000)
      ld%x(:) = [log(x_H2O),log(x_CO2),log(P_H2O),&
                 log(P_CO2),alphaG,log(P_H2),log(P_CH4),log(P_CO)]
      call lmder(fcn_, ld%m, ld%n, ld%x, ld%fvec, ld%fjac, ld%ldfjac, ld%ftol, ld%xtol, ld%gtol, ld%maxfev, &
                 ld%diag, ld%mode, ld%factor, ld%nprint, ld%info, ld%nfev, ld%njev, ld%ipvt, ld%qtf, &
                 ld%wa1, ld%wa2, ld%wa3, ld%wa4)
      if ((ld%info < 1 .and. ld%info > 4) .or. ld%x(5) < 0.0_dp) then
        ld%x(:) = [log(x_H2O),log(x_CO2),log(P_H2O),&
                   log(P_CO2),log(alphaG),log(P_H2),log(P_CH4),log(P_CO)]
        call lmder(fcn1_, ld%m, ld%n, ld%x, ld%fvec, ld%fjac, ld%ldfjac, ld%ftol, ld%xtol, ld%gtol, ld%maxfev, &
                   ld%diag, ld%mode, ld%factor, ld%nprint, ld%info, ld%nfev, ld%njev, ld%ipvt, ld%qtf, &
                   ld%wa1, ld%wa2, ld%wa3, ld%wa4)
        iflag = 1
        call fcn1_(ld%m,ld%n,ld%x,ld%fvec,ld%fjac,ld%ldfjac,iflag)
        error = sum(abs(ld%fvec))
        if ((ld%info < 1 .and. ld%info > 4) .or. error > 1.0e-7_dp) then
          ! Here we are assuming that no outgassing occurs, because we could
          ! not find a physical solution where outgassing did occur.
          x_H2O = xH2Otot
          x_CO2 = xCO2tot
          P_H2O = 0.0_dp
          P_CO2 = 0.0_dp
          alphaG = 0.0_dp
          P_H2 = 0.0_dp
          P_CH4 = 0.0_dp
          P_CO = 0.0_dp
        else
          x_H2O = exp(ld%x(1))
          x_CO2 = exp(ld%x(2))
          P_H2O = exp(ld%x(3))
          P_CO2 = exp(ld%x(4))
          alphaG = exp(ld%x(5))
          P_H2 = exp(ld%x(6))
          P_CH4 = exp(ld%x(7))
          P_CO = exp(ld%x(8))
        endif
      else
        x_H2O = exp(ld%x(1))
        x_CO2 = exp(ld%x(2))
        P_H2O = exp(ld%x(3))
        P_CO2 = exp(ld%x(4))
        alphaG = ld%x(5)
        P_H2 = exp(ld%x(6))
        P_CH4 = exp(ld%x(7))
        P_CO = exp(ld%x(8))
      endif

    end block

  contains

    subroutine fcn_s_(m_, n_, x_, fvec_, iflag_)
      integer, intent(in) :: m_
      integer, intent(in) :: n_
      real(dp), intent(in) :: x_(n_)
      real(dp), intent(out) :: fvec_(m_)
      integer, intent(inout) :: iflag_
      call fcn_s(d, m_, n_, x_, fvec_, iflag_)
    end subroutine

    subroutine fcn_(m_, n_, x_, fvec_, fjac_, ldfjac_, iflag_)
      integer, intent(in) :: m_
      integer, intent(in) :: n_
      integer, intent(in) :: ldfjac_
      integer, intent(inout) :: iflag_
      real(dp), intent(in) :: x_(n_)
      real(dp), intent(inout) :: fvec_(m_)
      real(dp), intent(inout) :: fjac_(ldfjac_,n_)
      call fcn(d, m_, n_, x_, fvec_, fjac_, ldfjac_, iflag_)
    end subroutine

    subroutine fcn1_(m_, n_, x_, fvec_, fjac_, ldfjac_, iflag_)
      integer, intent(in) :: m_
      integer, intent(in) :: n_
      integer, intent(in) :: ldfjac_
      integer, intent(inout) :: iflag_
      real(dp), intent(in) :: x_(n_)
      real(dp), intent(inout) :: fvec_(m_)
      real(dp), intent(inout) :: fjac_(ldfjac_,n_)
      call fcn1(d, m_, n_, x_, fvec_, fjac_, ldfjac_, iflag_)
    end subroutine

  end subroutine

  !> c-interface
  subroutine solve_gases_c(T, P, f_O2, mCO2tot, mH2Otot, &
                           P_H2O, P_H2, P_CO2, P_CO, &
                           P_CH4, alphaG, x_CO2, x_H2O, ierr) bind(c)
    use iso_c_binding, only: c_double, c_int
    real(c_double), intent(in) :: T, P, f_O2, mCO2tot, mH2Otot
    real(c_double), intent(out) :: P_H2O, P_H2, P_CO2, P_CO, P_CH4
    real(c_double), intent(out) :: alphaG, x_CO2, x_H2O
    integer(c_int), intent(out) :: ierr
    call solve_gases(T, P, f_O2, mCO2tot, mH2Otot, &
                     P_H2O, P_H2, P_CO2, P_CO, &
                     P_CH4, alphaG, x_CO2, x_H2O, ierr)
  end subroutine

  subroutine fcn(d, m, n, y, out, fjac, ldfjac, iflag)
    use volcgases_const, only: a_CO2, a_H2O, d_H2O
    type(SolveGasesData), target, intent(in) :: d
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: ldfjac
    integer, intent(inout) :: iflag
    real(dp), intent(in) :: y(n)
    real(dp), intent(inout) :: out(m)
    real(dp), intent(inout) :: fjac(ldfjac,n)

    real(dp) :: ln_x_H2O, ln_x_CO2, ln_H2O, ln_CO2
    real(dp) :: alphaG, ln_H2, ln_CH4, ln_CO
    real(dp) :: lnxH2O, lnxCO2, lnPH2O, lnPCO2
    real(dp) :: lnPH2, lnPCH4, lnPCO
    real(dp), parameter :: e = exp(1.0_dp)

    real(dp), pointer :: P
    real(dp), pointer :: F1, F2
    real(dp), pointer :: xCO2tot, xH2Otot
    real(dp), pointer :: C1, C2, C3

    P => d%P
    F1 => d%F1
    F2 => d%F2
    xCO2tot => d%xCO2tot
    xH2Otot => d%xH2Otot
    C1 => d%C1
    C2 => d%C2
    C3 => d%C3

    if (iflag == 1) then
      ln_x_H2O = y(1)
      ln_x_CO2 = y(2)
      ln_H2O = y(3)
      ln_CO2 = y(4)
      alphaG = y(5)
      ln_H2 = y(6)
      ln_CH4 = y(7)
      ln_CO = y(8)

      out(1) = exp(ln_H2O)+exp(ln_CO2)+exp(ln_H2)+exp(ln_CH4)+exp(ln_CO)-P
      out(2) = -ln_x_CO2+exp(ln_x_H2O)*d_H2O+a_CO2*ln_CO2+F1
      out(3) = -ln_x_H2O+a_H2O*ln_H2O+F2
      out(4) = -xH2Otot*P + (exp(ln_H2O)+exp(ln_H2) &
             +2.0_dp*exp(ln_CH4))*alphaG+(1.0_dp-alphaG)*exp(ln_x_H2O)*P
      out(5) = -xCO2tot*P + (exp(ln_CO2)+exp(ln_CO) &
             +exp(ln_CH4))*alphaG+(1-alphaG)*exp(ln_x_CO2)*P
      out(6) = log(C1)+ln_H2O-ln_H2
      out(7) = log(C2)+ln_CO2-ln_CO
      out(8) = log(C3)+ln_CO2+2.0_dp*ln_H2O-ln_CH4
    else if (iflag == 2) then
      lnxH2O = y(1)
      lnxCO2 = y(2)
      lnPH2O = y(3)
      lnPCO2 = y(4)
      alphaG = y(5)
      lnPH2 = y(6)
      lnPCH4 = y(7)
      lnPCO = y(8)

      fjac(1,:)=[ 0.0_dp,0.0_dp,E**lnPH2O,E**lnPCO2,0.0_dp,E**lnPH2,E**lnPCH4,E**lnPCO ]
      fjac(2,:)=[ d_H2O*E**lnxH2O,-1.0_dp,0.0_dp,a_CO2,0.0_dp,0.0_dp,0.0_dp,0.0_dp]
      fjac(3,:)=[ -1.0_dp,0.0_dp,a_H2O,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp]

      fjac(4,:)=[ (1.0_dp-alphaG)*E**lnxH2O*P,0.0_dp,alphaG*E**lnPH2O,0.0_dp,2.0_dp*E**lnPCH4 &
                  + E**lnPH2 + E**lnPH2O - E**lnxH2O*P, &
                  - alphaG*E**lnPH2,2.0_dp*alphaG*E**lnPCH4,0.0_dp ]

      fjac(5,:)=[ 0.0_dp,(1.0_dp - alphaG)*E**lnxCO2*P,0.0_dp,alphaG*E**lnPCO2,E**lnPCH4 &
                  + E**lnPCO + E**lnPCO2 - E**lnxCO2*P,0.0_dp, &
                  - alphaG*E**lnPCH4,alphaG*E**lnPCO ]
      fjac(6,:)=[0.0_dp,0.0_dp,1.0_dp,0.0_dp,0.0_dp,-1.0_dp,0.0_dp,0.0_dp]
      fjac(7,:)=[0.0_dp,0.0_dp,0.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp,-1.0_dp]
      fjac(8,:)=[0.0_dp,0.0_dp,2.0_dp,1.0_dp,0.0_dp,0.0_dp,-1.0_dp,0.0_dp]

    endif

  end subroutine

  subroutine fcn1(d, m, n, y, out, fjac, ldfjac, iflag)
    use volcgases_const, only: a_CO2, a_H2O, d_H2O
    type(SolveGasesData), target, intent(in) :: d
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: ldfjac
    integer, intent(inout) :: iflag
    real(dp), intent(in) :: y(n)
    real(dp), intent(inout) :: out(m)
    real(dp), intent(inout) :: fjac(ldfjac,n)

    real(dp) :: ln_x_H2O, ln_x_CO2, ln_H2O, ln_CO2
    real(dp) :: lnalphaG, ln_H2, ln_CH4, ln_CO
    real(dp) :: lnxH2O, lnxCO2, lnPH2O, lnPCO2
    real(dp) :: lnPH2, lnPCH4, lnPCO
    real(dp), parameter :: e = exp(1.0_dp)

    real(dp), pointer :: P
    real(dp), pointer :: F1, F2
    real(dp), pointer :: xCO2tot, xH2Otot
    real(dp), pointer :: C1, C2, C3

    P => d%P
    F1 => d%F1
    F2 => d%F2
    xCO2tot => d%xCO2tot
    xH2Otot => d%xH2Otot
    C1 => d%C1
    C2 => d%C2
    C3 => d%C3

    if (iflag == 1) then
      ln_x_H2O = y(1)
      ln_x_CO2 = y(2)
      ln_H2O = y(3)
      ln_CO2 = y(4)
      lnalphaG = y(5)
      ln_H2 = y(6)
      ln_CH4 = y(7)
      ln_CO = y(8)

      out(1) = exp(ln_H2O)+exp(ln_CO2)+exp(ln_H2)+exp(ln_CH4)+exp(ln_CO)-P
      out(2) = -ln_x_CO2+exp(ln_x_H2O)*d_H2O+a_CO2*ln_CO2+F1
      out(1) = -ln_x_H2O+a_H2O*ln_H2O+F2
      out(4) = -xH2Otot*P + (exp(ln_H2O)+exp(ln_H2)+2*exp(ln_CH4))*exp(lnalphaG)+(1-exp(lnalphaG))*exp(ln_x_H2O)*P
      out(5) = -xCO2tot*P + (exp(ln_CO2)+exp(ln_CO)+exp(ln_CH4))*exp(lnalphaG)+(1-exp(lnalphaG))*exp(ln_x_CO2)*P
      out(6) = log(C1)+ln_H2O-ln_H2
      out(7) = log(C2)+ln_CO2-ln_CO
      out(8) = log(C3)+ln_CO2+2*ln_H2O-ln_CH4
    else if (iflag == 2) then
      lnxH2O = y(1)
      lnxCO2 = y(2)
      lnPH2O = y(3)
      lnPCO2 = y(4)
      lnalphaG = y(5)
      lnPH2 = y(6)
      lnPCH4 = y(7)
      lnPCO = y(8)

      fjac(1,:)=[ 0.0_dp,0.0_dp,E**lnPH2O,E**lnPCO2,0.0_dp,E**lnPH2,E**lnPCH4,E**lnPCO ]
      fjac(2,:)=[ d_H2O*E**lnxH2O,-1.0_dp,0.0_dp,a_CO2,0.0_dp,0.0_dp,0.0_dp,0.0_dp]
      fjac(3,:)=[ -1.0_dp,0.0_dp,a_H2O,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp ]

      fjac(4,:)=[ E**lnxH2O*(1.0_dp - E**lnalphaG)*P,0.0_dp,E**(lnalphaG + lnPH2O),0.0_dp, &
                  - E**lnalphaG*(2.0_dp*E**lnPCH4 + E**lnPH2 + E**lnPH2O) - &
                  E**(lnalphaG + lnxH2O)*P,E**(lnalphaG + lnPH2), &
                  - 2.0_dp*E**(lnalphaG + lnPCH4),0.0_dp ]

      fjac(5,:)=[ 0.0_dp,E**lnxCO2*(1.0_dp - E**lnalphaG)*P,0.0_dp,E**(lnalphaG + lnPCO2), &
                  - E**lnalphaG*(E**lnPCH4 + E**lnPCO + E**lnPCO2) &
                  - E**(lnalphaG + lnxCO2)*P,0.0_dp,E**(lnalphaG + lnPCH4), &
                  - E**(lnalphaG + lnPCO) ]
      fjac(6,:)=[0.0_dp,0.0_dp,1.0_dp,0.0_dp,0.0_dp,-1.0_dp,0.0_dp,0.0_dp]
      fjac(7,:)=[0.0_dp,0.0_dp,0.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp,-1.0_dp]
      fjac(8,:)=[0.0_dp,0.0_dp,2.0_dp,1.0_dp,0.0_dp,0.0_dp,-1.0_dp,0.0_dp]

    endif
  end subroutine

  subroutine fcn_s(d, m, n, x, fvec, iflag)
    use volcgases_const, only: a_CO2, a_H2O
    type(SolveGasesData), target, intent(in) :: d
    integer, intent(in) :: m
    integer, intent(in) :: n
    real(dp), intent(in) :: x(n)
    real(dp), intent(out) :: fvec(m)
    integer, intent(inout) :: iflag

    real(dp) :: PCO2
    real(dp), pointer :: P
    real(dp), pointer :: F1, F2
    real(dp), pointer :: xCO2tot, xH2Otot
    real(dp), pointer :: C1, C2, C3
    real(dp), parameter :: e = exp(1.0_dp)

    PCO2 = x(1)

    P => d%P
    F1 => d%F1
    F2 => d%F2
    xCO2tot => d%xCO2tot
    xH2Otot => d%xH2Otot
    C1 => d%C1
    C2 => d%C2
    C3 => d%C3

    fvec(1) = ( -1 * ( e )**( ( F1 + a_CO2 * log( PCO2 ) ) ) * ( P )**( -1 ) &
                * ( P + ( -1 * PCO2 + -1 * C2 * PCO2 ) ) * ( ( 1 + ( C1 + C3 * PCO2 ) &
                ) )**( -1 ) + ( ( P )**( -1 ) * ( P + ( -1 * PCO2 + -1 * C2 * PCO2 ) &
                ) * ( ( 1 + ( C1 + C3 * PCO2 ) ) )**( -1 ) * xCO2tot + ( -1 * ( e &
                )**( ( F2 + a_H2O * log( (P-PCO2-C2*PCO2)/(1+C1+C3*PCO2) ) ) ) * ( -1 * ( P )**( -1 ) * PCO2 + &
                xCO2tot ) + ( ( e )**( ( F1 + a_CO2 * log( PCO2 ) ) ) * xH2Otot &
                + -1 * ( P )**( -1 ) * PCO2 * xH2Otot ) ) ) )
  end subroutine

  function LMData_create(m, n, ftol, xtol, gtol, maxfev) result(ld)
    integer, intent(in) :: m, n
    real(dp), intent(in) :: ftol, xtol, gtol
    integer, intent(in) :: maxfev
    type(LMData) :: ld

    ld%m = m
    ld%n = n
    allocate(ld%x(n))
    allocate(ld%fvec(m))
    ld%ftol = ftol
    ld%xtol = xtol
    ld%gtol = gtol
    ld%maxfev = maxfev
    ld%epsfcn = 0.0_dp
    allocate(ld%diag(n))
    ld%mode = 1
    ld%factor = 100.0_dp
    ld%nprint = -1
    ld%ldfjac = m
    allocate(ld%fjac(ld%ldfjac,n))
    allocate(ld%ipvt(n))
    allocate(ld%qtf(n))
    allocate(ld%wa1(n))
    allocate(ld%wa2(n))
    allocate(ld%wa3(n))
    allocate(ld%wa4(m))

  end function

  pure function single_linspace(from, to, n, i) result(val)
    real(dp), intent(in) :: from, to
    integer, intent(in) :: n, i
    real(dp) :: val

    real(dp) :: range

    range = to - from
    val = from
    if (n == 0 .or. n == 1) return
    val = from + range * (i - 1) / (n - 1)
  end function

end module