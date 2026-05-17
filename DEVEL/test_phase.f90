program test
  implicit none
  real*8 :: a, b, res
  a = 0.5d0
  b = 1.5d0
  res = (-1)**(a+b)
  print *, "Result is: ", res
end program test
