program test_mute
  implicit none
  include 'mpif.h'
  integer :: ierr, rank

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

  if (rank /= 0) then
    open(unit=6, file='/dev/null', status='unknown')
  end if

  write(*,*) 'Hello from rank ', rank
  write(6,*) 'Also from rank ', rank

  call MPI_FINALIZE(ierr)
end program test_mute
