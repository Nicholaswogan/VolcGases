
module volc_wrapper
  use iso_c_binding, only: c_double, c_int
  use volc, only: solve_gases
  implicit none
contains
  subroutine solve_gases_wrap(T,PP,f_O2,mCO2tot,mH2Otot, &
                             & P_H2O, P_H2, P_CO2, P_CO, &
                             & P_CH4, alphaG, x_CO2, x_H2O) bind(c)
    real(c_double), intent(in) :: T, PP, f_O2, mCO2tot, mH2Otot
    real(c_double), intent(out) :: P_H2O, P_H2, P_CO2, P_CO, P_CH4
    real(c_double), intent(out) :: alphaG, x_CO2, x_H2O
    call solve_gases(T,PP,f_O2,mCO2tot,mH2Otot, &
                   & P_H2O, P_H2, P_CO2, P_CO, &
                   & P_CH4, alphaG, x_CO2, x_H2O)
  end subroutine
end module