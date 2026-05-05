program test_err
  implicit none
  include 'mpif.h'
  integer :: ierr, rank

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

  if (rank /= 0) then
    open(unit=6, file='/dev/null', status='unknown')
    open(unit=0, file='/dev/null', status='unknown')
  end if

  write(*,*) 'Hello stdout from rank ', rank
  write(0,*) 'Hello stderr from rank ', rank

  call MPI_FINALIZE(ierr)
end program test_err
