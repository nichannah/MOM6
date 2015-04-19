
module numericalunits

use, intrinsic :: iso_fortran_env, only: error_unit
use MOM_grid, only : ocean_grid_type

implicit none

real, protected :: nu_A  ! electric current
real, protected :: nu_kg ! mass
real, protected :: nu_m  ! length
real, protected :: nu_H  ! thickness
real, protected :: nu_s  ! time
real, protected :: nu_K  ! thermodynamic temperature
real, protected :: nu_cd ! liminous intensity
real, protected :: nu_mole ! amount of a substance

real, protected :: nu_L ! litres
real, protected :: nu_mL ! millilitres
real, protected :: nu_nm ! nanometres
real, protected :: nu_v ! velocity
real, protected :: nu_m2 ! square metres

contains

subroutine nu_init()
    call nu_reset_units()
end subroutine nu_init

subroutine nu_reset_units()

    call random_seed()
    call random_number(nu_A)
    call random_number(nu_kg)
    call random_number(nu_m)
    call random_number(nu_H)
    call random_number(nu_s)
    call random_number(nu_K)
    call random_number(nu_cd)
    call random_number(nu_mole)

    nu_L = 1e-3 * nu_m**3
    nu_mL = 1e-3 * nu_L
    nu_nm = 1e-9 * nu_m
    nu_v = nu_m / nu_s
    nu_m2 = nu_m**2

end subroutine nu_reset_units

subroutine nu_set_grid_units(grid, grid_with_units)
    type(ocean_grid_type), intent(in) :: grid
    type(ocean_grid_type), intent(out) :: grid_with_units

    grid_with_units%dy_Cu = grid%dy_Cu * nu_m
    grid_with_units%IareaT = grid%IareaT * (1 / nu_m2)
    grid_with_units%IdxT = grid%IdxT * (1 / nu_m)

end subroutine nu_set_grid_units

! FIXME: put this somewhere better.
subroutine assert_allclose(actual, desired, msg, abs_tol, rel_tol)

  real, dimension(:), intent(in) :: actual, desired
  real, intent(in), optional :: rel_tol, abs_tol
  character(len=*), intent(in) :: msg
  real :: rtol, atol

  if (present(rel_tol)) then
    rtol = rel_tol
  else
    rtol = 1e-12
  endif

  if (present(abs_tol)) then
    atol = abs_tol
  else
    atol = 0
  endif

  if (any(abs(actual - desired) > (atol + rtol * abs(desired)))) then
     write(error_unit, *) msg
     stop 'assert_allclose triggered, see stderr.'
  endif

end subroutine assert_allclose

end module numericalunits
