module volcgases_const
  use iso_fortran_env, only: dp => real64
  implicit none

  !! Constants used in solve_gases. References:
  !!
  !! * [1] https://doi.org/10.3847/PSJ/abb99e
  !! * [2] https://doi.org/10.1016/j.gca.2012.08.035

  real(dp), parameter :: M_H2O = 18.01528_dp !! Molar mass of H2O (g/mol)
  real(dp), parameter :: M_CO2 = 44.01_dp !! Molar mass of CO2 (g/mol)

  !> Solubility constant for CO2 for Mt. Etna composition.
  !> It is equal to the stuff after the second term in Eq. A1 in [1]
  real(dp), parameter :: A1 = -0.4200250000201988_dp 
  !> Solubility constant for H2O for Mt. Etna composition.
  !> It is equal to the stuff after the second term in Eq. A2 in [1]
  real(dp), parameter :: A2 = -2.59560737813789_dp
  !> Inverse of molar mass of magma (mol of magma / g of magma).
  !> See Table 1 in [1]
  real(dp), parameter :: x = 0.01550152865954013_dp

  real(dp), parameter :: C_CO2 = 0.14_dp !! Solubility constant. See Table 4 in [1] 
  real(dp), parameter :: C_H2O = 0.02_dp !! Solubility constant. See Table 6 in [2]
  real(dp), parameter :: a_H2O = 0.54_dp !! Solubility constant. See Table 1 in [1]
  real(dp), parameter :: a_CO2 = 1.0_dp !! Solubility constant. See Table 1 in [1]
  real(dp), parameter :: d_H2O = 2.3_dp !! Solubility constant. See Table 1 in [1]
  
end module