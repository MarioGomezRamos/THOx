program test_nr
  implicit none
  integer :: nrcc, rmaxcc
  real*8 :: hcm
  rmaxcc = 60
  hcm = 0.05
  nrcc = nint(rmaxcc/hcm)+1
  print *, "nrcc = ", nrcc
end program test_nr
