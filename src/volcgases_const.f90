module volcgases_const
  use iso_fortran_env, only: dp => real64
  implicit none

  real(dp), parameter :: A1 = -0.4200250000201988_dp
  real(dp), parameter :: A2 = -2.59560737813789_dp
  real(dp), parameter :: M_H2O = 18.01528_dp
  real(dp), parameter :: M_CO2 = 44.01_dp
  real(dp), parameter :: C_CO2 = 0.14_dp
  real(dp), parameter :: C_H2O = 0.02_dp
  real(dp), parameter :: x = 0.01550152865954013_dp
  real(dp), parameter :: a_H2O = 0.54_dp
  real(dp), parameter :: a_CO2 = 1.0_dp
  real(dp), parameter :: d_H2O = 2.3_dp
  
end module