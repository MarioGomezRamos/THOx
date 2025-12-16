!------------------------------------------------------------------------------
! Module: precision
! Purpose: Define precision parameters for HPRMAT
!------------------------------------------------------------------------------
module precision
  implicit none

  ! Integer precision
  integer, parameter :: sint = selected_int_kind(9)    ! 4-byte integers
  integer, parameter :: lint = selected_int_kind(19)   ! 8-byte integers

  ! Real precision
  integer, parameter :: sp = selected_real_kind(6)     ! Single precision (4 bytes)
  integer, parameter :: dp = selected_real_kind(15)    ! Double precision (8 bytes)
  integer, parameter :: spreal = sp
  integer, parameter :: dpreal = dp

  ! Mathematical constants
  real(dp), parameter :: pi = 3.14159265358979323846_dp
  complex(dp), parameter :: iu = (0.0_dp, 1.0_dp)

end module precision
