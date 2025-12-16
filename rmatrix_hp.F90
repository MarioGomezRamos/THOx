!------------------------------------------------------------------------------
! HPRMAT: High-Performance R-Matrix Solver
!
! This module provides a drop-in replacement for Pierre Descouvemont's
! R-matrix package (CPC 200, 199-219, 2016) with improved performance.
!
! Key improvements:
!   1. Uses linear equation solving (ZGESV) instead of matrix inversion
!   2. Supports multiple solver backends (CPU/GPU)
!   3. Optional mixed precision for faster computation
!
! Interface compatibility:
!   - rmat_ini(nr, ns, rmax, zrma) - Initialize R-matrix calculation
!   - rmatrix(...) - Compute R-matrix and collision matrix
!   - wf_print(...) - Print wave function
!
! Solver selection (via module variable):
!   solver_type = 1: Dense LAPACK ZGESV (default, reference)
!   solver_type = 2: Mixed Precision (slightly faster)
!   solver_type = 3: Woodbury-Kinetic (CPU optimized)
!   solver_type = 4: GPU cuSOLVER (fastest with GPU)
!
! Author: Jin Lei
! Date: December 2025
!------------------------------------------------------------------------------
module rmat_hp_mod
  implicit none
  real*8 :: rmax0
  real*8, allocatable :: wle(:), xle(:), tc(:,:,:), blo0(:,:), blo1(:,:), blo2(:,:)
  real*8, allocatable :: q1(:,:), q2(:,:)

  ! Solver configuration
  integer, public :: solver_type = 1  ! Default: Dense LAPACK

end module rmat_hp_mod

!------------------------------------------------------------------------------
! rmat_ini: Initialize R-matrix calculation
!
! This subroutine initializes the R-matrix calculation (should be called
! before subroutine rmatrix).
!
! It computes:
!   - The roots and weights of the Legendre quadrature
!   - The matrix elements of the kinetic energy and of the Bloch operator
!
! Inputs:
!   nr   - Number of Lagrange functions per interval
!   ns   - Number of intervals
!   rmax - R-matrix channel radius a
!
! Output:
!   zrma(ns*nr) - Array containing the abscissas of the Lagrange mesh
!------------------------------------------------------------------------------
subroutine rmat_ini_hp(nr, ns, rmax, zrma)
  use rmat_hp_mod
  implicit real*8(a,b,d-h,o-z)
  dimension zrma(ns*nr)

  ! Deallocate if already allocated
  if (allocated(wle)) deallocate(wle)
  if (allocated(xle)) deallocate(xle)
  if (allocated(tc)) deallocate(tc)
  if (allocated(blo0)) deallocate(blo0)
  if (allocated(blo1)) deallocate(blo1)
  if (allocated(blo2)) deallocate(blo2)
  if (allocated(q1)) deallocate(q1)
  if (allocated(q2)) deallocate(q2)

  allocate(tc(nr,nr,ns), blo0(nr,nr), blo1(nr,nr), blo2(nr,nr), &
           wle(nr), xle(nr), q1(nr,ns), q2(nr,ns))

  rmax0 = rmax / ns

  ! Get roots and weights of Legendre quadrature
  call legzo_hp(nr, xle, wle)

  ! Compute matrix elements of kinetic energy and Bloch operator
  do is = 1, ns
    do i1 = 1, nr
      xi = xle(i1)
      xi2 = xi * (1 - xi)

      if (is == 1) then
        xx = 4*nr*(nr+1) + 3 + (1 - 6*xi) / xi2
        tc(i1, i1, 1) = xx / (3 * xi2)
        blo0(i1, i1) = 1 / xi2
      else
        xlb = xi / (1 - xi) * (nr*(nr+1) - 1/(1-xi))
        xla = (1 - xi) / xi * (-nr*(nr+1) + 1/xi)
        tc(i1, i1, is) = (nr*nr + nr + 6 - 2/xi2) / (3*xi2) + xlb - xla
        blo1(i1, i1) = (1 - xi) / xi
        blo2(i1, i1) = xi / (1 - xi)
      end if

      do i2 = 1, i1-1
        xj = xle(i2)
        xj2 = xj * (1 - xj)

        if (is == 1) then
          xx = nr*(nr+1) + 1 + (xi + xj - 2*xi*xj) / (xi-xj)**2 - 1/(1-xi) - 1/(1-xj)
          tc(i1, i2, 1) = xx / sqrt(xi2 * xj2)
          blo0(i1, i2) = 1 / sqrt(xi2 * xj2)
        else
          yy = sqrt(xj2/xi2**3) * (2*xi*xj + 3*xi - xj - 4*xi**2) / (xj-xi)**2
          xlb = sqrt(xi*xj/(1-xi)/(1-xj)) * (nr*(nr+1) - 1/(1-xj))
          xla = sqrt((1-xi)*(1-xj)/xi/xj) * (-nr*(nr+1) + 1/xj)
          tc(i1, i2, is) = yy + xlb - xla
          blo1(i1, i2) = sqrt((1-xi)*(1-xj)/xi/xj)
          blo2(i1, i2) = sqrt(xi*xj/(1-xi)/(1-xj))
        end if

        if (mod(i1+i2, 2) == 1) then
          tc(i1, i2, is) = -tc(i1, i2, is)
          blo0(i1, i2) = -blo0(i1, i2)
          blo1(i1, i2) = -blo1(i1, i2)
          blo2(i1, i2) = -blo2(i1, i2)
        end if

        tc(i2, i1, is) = tc(i1, i2, is)
        blo0(i2, i1) = blo0(i1, i2)
        blo1(i2, i1) = blo1(i1, i2)
        blo2(i2, i1) = blo2(i1, i2)
      end do

      if (is == 1) then
        q2(:, 1) = 1 / sqrt(xle(:) * (1 - xle(:)))
        q1(:, 1) = 0
      else
        q2(:, is) = sqrt(xle(:) / (1 - xle(:)))
        q1(:, is) = -1 / q2(:, is)
      end if

      if (mod(nr, 2) == 1) q2(1:nr, is) = -q2(1:nr, is)
      q1(1:nr:2, is) = -q1(1:nr:2, is)
      q2(1:nr:2, is) = -q2(1:nr:2, is)
    end do
  end do

  tc = tc / rmax0**2
  blo0 = blo0 / rmax0
  blo1 = blo1 / rmax0
  blo2 = blo2 / rmax0
  q1 = q1 / sqrt(rmax0)
  q2 = q2 / sqrt(rmax0)

  nn = 0
  do is = 1, ns
    zrma(nn+1:nn+nr) = ((is-1) + xle(1:nr)) * rmax0
    nn = nn + nr
  end do

end subroutine rmat_ini_hp

!------------------------------------------------------------------------------
! rmatrix: Compute R-matrix and collision matrix
!
! This is the main subroutine that computes the R-matrix and the associated
! collision matrix using high-performance linear equation solvers.
!
! Key improvement: Uses ZGESV (linear solve) instead of matrix inversion,
! which is significantly faster for large matrices.
!
! Extra parameter:
!   isolver - Optional solver type (1=Dense LAPACK, 2=Mixed Precision,
!             3=Woodbury-Kinetic, 4=GPU cuSOLVER)
!------------------------------------------------------------------------------
subroutine rmatrix_hp(nch, lval, qk, eta, rmax, nr, ns, cpot, cu, &
                   ncp1, ndim, nopen, twf, cf, nwf1, nwf2, nc, nvc, ncp2, cpnl, isolver)
  use rmat_hp_mod
  use rmat_solvers, only: solve_rmatrix_woodbury, solve_rmatrix_gpu
  implicit real*8(a,b,d-h,o-z)
  implicit complex*16(c)
  dimension lval(nch), qk(nch), eta(nch), cpot(ncp1, ndim, ndim)
  dimension cu(ndim, ndim), cf(nwf1, nwf2, nc), nvc(nc), cx(3)
  dimension cpnl(ncp2, ndim, ndim)
  integer, intent(in), optional :: isolver
  logical twf, tnl
  integer :: local_solver
  logical, save :: first_call = .true.
  allocatable ch(:,:,:), crma(:,:,:,:), co(:), cop(:)
  allocatable cz(:,:), crma0(:,:,:), crma2(:,:,:), fc(:), dfc(:), gc(:), dgc(:)
  allocatable cfp(:,:), npo(:,:), xc(:,:)

  ! Local variables for high-performance solver
  complex*16, allocatable :: B_vector(:), Rmat_local(:,:)
  integer :: IPIV_local(nch)
  integer :: INFO_local

  ! Set local solver type
  if (present(isolver)) then
    local_solver = isolver
  else
    local_solver = solver_type
  end if

  ! Print solver info (only once)
  if (first_call) then
    write(*,'(a)') ' '
    write(*,'(a)') '=== HPRMAT Solver Selection ==='
    write(*,'(a,i2)') ' Solver type = ', local_solver
    select case (local_solver)
      case (1)
        write(*,'(a)') ' Using: Dense LAPACK ZGESV (reference)'
      case (2)
        write(*,'(a)') ' Using: Mixed Precision (single prec factorization)'
      case (3)
        write(*,'(a)') ' Using: Woodbury-Kinetic (CPU optimized)'
      case (4)
        write(*,'(a)') ' Using: GPU cuSOLVER'
      case default
        write(*,'(a)') ' Using: Default Dense LAPACK ZGESV'
    end select
    write(*,'(a)') '==============================='
    write(*,'(a)') ' '
    first_call = .false.
  end if

  tnl = ncp2 /= 0
  if (tnl .and. ns /= 1) then
    print*, 'ns must be ns=1 for non-local potentials'
    stop
  end if

  lmax = maxval(lval) + 1
  ntot = nch * nr

  allocate(fc(lmax), dfc(lmax), gc(lmax), dgc(lmax))
  allocate(ch(ntot, ntot, ns), crma(nch, nch, 3, ns), co(nch), &
           cop(nch), cz(nch, nch), crma2(nch, nch, ns), &
           crma0(nch, nch, 0:ns), cfp(nch, 0:ns), xc(nch, 2), npo(nch, 2))
  allocate(B_vector(nr), Rmat_local(nch, nch))

  nopen = 0
  nclo = 0
  cu = 0
  cfp = 0

  ! Store Coulomb and Whittaker functions
  do i1 = 1, nch
    l = lval(i1)
    ll = l + 1
    if (qk(i1) > 0) then
      nopen = nopen + 1
      npo(nopen, 1) = i1
      xl = l
      call coulfg(qk(i1)*rmax, eta(i1), xl, xl, fc, gc, dfc, dgc, 1, 0, ifail)
      co(nopen) = dcmplx(gc(ll), fc(ll))
      cop(nopen) = dcmplx(dgc(ll), dfc(ll)) * qk(i1) * rmax
    else
      nclo = nclo + 1
      npo(nclo, 2) = i1
      qk2 = -qk(i1)
      ie = 0
      call whit(eta(i1), rmax, qk2, -qk2**2, l, fc, dfc, ie)
      xc(i1, 1) = rmax * dfc(l+1) / fc(l+1)
      xc(i1, 2) = fc(l+1)
    end if
  end do

  ! Build and solve for each interval
  rma0 = 0
  ms = -nr
  ch = 0

  do is = 1, ns
    rmax1 = (is - 1) * rmax0
    rmax2 = is * rmax0
    ms = ms + nr
    m1 = 0

    do i1 = 1, nch
      ch(m1+1:m1+nr, m1+1:m1+nr, is) = tc(:,:,is)
      fac = lval(i1) * (lval(i1) + 1)

      do ir = 1, nr
        xx = rmax1 + xle(ir) * rmax0
        ch(m1+ir, m1+ir, is) = ch(m1+ir, m1+ir, is) + fac / xx**2
        if (qk(i1) > 0) ch(m1+ir, m1+ir, is) = ch(m1+ir, m1+ir, is) - qk(i1)**2
        if (qk(i1) < 0) ch(m1+ir, m1+ir, is) = ch(m1+ir, m1+ir, is) + qk(i1)**2
      end do

      if (qk(i1) < 0) then
        if (is == 1) then
          ch(m1+1:m1+nr, m1+1:m1+nr, is) = &
            ch(m1+1:m1+nr, m1+1:m1+nr, is) - xc(i1, 1) * blo0 / rmax2
        else
          ch(m1+1:m1+nr, m1+1:m1+nr, is) = &
            ch(m1+1:m1+nr, m1+1:m1+nr, is) - xc(i1, 1) * (blo2/rmax2 - blo1/rmax1)
        end if
      end if

      m2 = 0
      do i2 = 1, nch
        ii = (is - 1) * nr * nr
        do ir = 1, nr
          ch(m1+ir, m2+ir, is) = ch(m1+ir, m2+ir, is) + cpot(ms+ir, i1, i2)
          if (.not. tnl) cycle
          do irp = 1, nr
            ii = ii + 1
            ch(m1+ir, m2+irp, is) = ch(m1+ir, m2+irp, is) + &
              sqrt(wle(ir)*wle(irp)) * rmax0 * cpnl(ii, i1, i2)
          end do
        end do
        m2 = m2 + nr
      end do
      m1 = m1 + nr
    end do

    ! === HIGH-PERFORMANCE SOLVER ===
    ! Use linear equation solving instead of matrix inversion
    B_vector(:) = q2(:, is)

    ! DEBUG - check ch matrix has imaginary part
    if (is == 1) then
      write(*,*) '=== ch matrix check ==='
      write(*,'(a,2e15.6)') ' ch(1,1,1) = ', real(ch(1,1,1)), aimag(ch(1,1,1))
      write(*,'(a,2e15.6)') ' ch(nr,nr,1)=', real(ch(nr,nr,1)), aimag(ch(nr,nr,1))
    endif

    select case (local_solver)
    case (1)
      ! Dense LAPACK - solve using ZGESV
      call solve_rmatrix_hp(ch(:,:,is), B_vector, nch, nr, 1.0d0, Rmat_local)
    case (2)
      ! Mixed precision
      call solve_rmatrix_hp_mixed(ch(:,:,is), B_vector, nch, nr, 1.0d0, Rmat_local)
    case (3)
      ! Woodbury-Kinetic - exploits matrix structure for O(n^2) instead of O(n^3)
      call solve_rmatrix_woodbury(ch(:,:,is), B_vector, nch, nr, 1.0d0, Rmat_local)
    case (4)
      ! GPU cuSOLVER
      call solve_rmatrix_gpu(ch(:,:,is), B_vector, nch, nr, 1.0d0, Rmat_local)
    case default
      ! Default to dense LAPACK
      call solve_rmatrix_hp(ch(:,:,is), B_vector, nch, nr, 1.0d0, Rmat_local)
    end select

    crma(:,:,1,is) = Rmat_local

    ! For propagation, also need R01, R10
    if (is > 1 .or. ns > 1) then
      ! Compute R matrices with q1 and q2 boundaries
      ! solve_rmatrix_hp_4 outputs: R00=<q1|C^-1|q1>, R01=<q1|C^-1|q2>, R11=<q2|C^-1|q2>
      ! But Pierre's convention is: crma(:,:,1)=R11, crma(:,:,2)=R01, crma(:,:,3)=R00
      call solve_rmatrix_hp_4(ch(:,:,is), q1(:,is), q2(:,is), nch, nr, 1.0d0, &
                              crma(:,:,3,is), crma(:,:,2,is), crma(:,:,1,is))
    end if

    ! Propagate R-matrix
    if (is == 1) then
      crma0(:,:,1) = crma(:,:,1,1) / rmax2
    else
      cz = crma(:,:,3,is) + crma0(:,:,is-1) * rmax1
      call cminv_sym_hp(cz, nch, nch)
      crma0(:,:,is) = matmul(transpose(crma(:,:,2,is)), cz)
      crma0(:,:,is) = matmul(crma0(:,:,is), crma(:,:,2,is))
      crma0(:,:,is) = crma(:,:,1,is) - crma0(:,:,is)
      crma0(:,:,is) = crma0(:,:,is) / rmax2
    end if
  end do

  ! Calculation of matrices Z_I and Z_O
  do m1 = 1, nopen
    i1 = npo(m1, 1)
    do m2 = 1, nopen
      i2 = npo(m2, 1)
      crma2(m1, m2, 1:ns) = crma0(i1, i2, 1:ns) * sqrt(qk(i1)/qk(i2))
    end do
  end do

  do m1 = 1, nopen
    cz(1:nopen, m1) = -crma2(1:nopen, m1, ns) * cop(m1)
    cu(1:nopen, m1) = -crma2(1:nopen, m1, ns) * conjg(cop(m1))
    cz(m1, m1) = cz(m1, m1) + co(m1)
    cu(m1, m1) = cu(m1, m1) + conjg(co(m1))
  end do

  ! Solve for collision matrix using linear solve instead of inversion
  call ZGESV(nopen, nopen, cz, nch, IPIV_local, cu(1:nopen, 1:nopen), nch, INFO_local)
  if (INFO_local /= 0) stop 'Error in ZGESV for collision matrix'

  ! Wave function calculation
  if (.not. twf) then
    deallocate(fc, dfc, gc, dgc, ch, crma, co, cop, cz, crma2, crma0, cfp, xc, npo)
    deallocate(B_vector, Rmat_local)
    return
  end if

  ! Loop over entrance channels
  do ic = 1, nc
    nvex = nvc(ic)
    i2 = npo(nvex, 1)

    do m1 = 1, nopen
      i1 = npo(m1, 1)
      cfp(i1, ns) = -cu(m1, nvex) * cop(m1) * sqrt(qk(i2)/qk(i1))
      if (i1 == i2) cfp(i1, ns) = cfp(i1, ns) + conjg(cop(m1))
    end do
    cfp(:, ns) = cfp(:, ns) / rmax

    do is = ns-1, 1, -1
      rin = is * rmax0
      cz = crma(:,:,3,is+1) + rin * crma0(:,:,is)
      call cminv_nsym_hp(cz, nch, nch)
      cz = matmul(cz, crma(:,:,2,is+1))
      cfp(:, is) = matmul(cz, cfp(:, is+1))
    end do

    do is = 1, ns
      m1 = 0
      do i1 = 1, nch
        do ir1 = 1, nr
          m1 = m1 + 1
          cy = 0
          m2 = 0
          do i2 = 1, nch
            do ir2 = 1, nr
              m2 = m2 + 1
              cy = cy + ch(m1, m2, is) * (q2(ir2, is)*cfp(i2, is) - q1(ir2, is)*cfp(i2, is-1))
            end do
          end do
          mm1 = (is - 1)*nr + ir1
          cf(mm1, i1, ic) = cy
        end do
      end do
    end do

    ! Bound state amplitudes
    do m1 = 1, nclo
      i1 = npo(m1, 2)
      cy = sum(cf((ns-1)*nr+1:ns*nr, i1, ic) * q2(:, ns))
      cu(m1+nopen, ic) = cy / xc(i1, 2)
    end do

    ! Wave function normalization
    do i1 = 1, nch
      mm = 0
      do is = 1, ns
        cf(mm+1:mm+nr, i1, ic) = cf(mm+1:mm+nr, i1, ic) / sqrt(rmax0 * wle(1:nr))
        mm = mm + nr
      end do
    end do
  end do

  deallocate(fc, dfc, gc, dgc, ch, crma, co, cop, cz, crma2, crma0, cfp, xc, npo)
  deallocate(B_vector, Rmat_local)

end subroutine rmatrix_hp

!------------------------------------------------------------------------------
! Internal solver routines
!------------------------------------------------------------------------------
subroutine solve_rmatrix_hp(cmat, B_vector, nch, nlag, normfac, Rmat)
  implicit none
  integer, intent(in) :: nch, nlag
  real*8, intent(in) :: normfac
  complex*16, intent(in) :: cmat(nch*nlag, nch*nlag), B_vector(nlag)
  complex*16, intent(out) :: Rmat(nch, nch)

  complex*16, allocatable :: A_work(:,:), X_vector(:,:)
  integer, allocatable :: IPIV(:)
  integer :: ich, ichp, ir, ntotal, INFO

  ntotal = nch * nlag
  allocate(A_work(ntotal, ntotal), X_vector(ntotal, nch), IPIV(ntotal))

  A_work = cmat
  X_vector = (0.d0, 0.d0)
  do ich = 1, nch
    do ir = 1, nlag
      X_vector((ich-1)*nlag + ir, ich) = B_vector(ir)
    end do
  end do

  call ZGESV(ntotal, nch, A_work, ntotal, IPIV, X_vector, ntotal, INFO)
  if (INFO /= 0) stop 'Error in solve_rmatrix_hp ZGESV'

  Rmat = (0.d0, 0.d0)
  do ichp = 1, nch
    do ich = 1, nch
      do ir = 1, nlag
        Rmat(ich, ichp) = Rmat(ich, ichp) + B_vector(ir) * X_vector(ir + (ich-1)*nlag, ichp)
      end do
    end do
  end do
  Rmat = Rmat * normfac

  deallocate(A_work, X_vector, IPIV)
end subroutine solve_rmatrix_hp

subroutine solve_rmatrix_hp_mixed(cmat, B_vector, nch, nlag, normfac, Rmat)
  implicit none
  integer, intent(in) :: nch, nlag
  real*8, intent(in) :: normfac
  complex*16, intent(in) :: cmat(nch*nlag, nch*nlag), B_vector(nlag)
  complex*16, intent(out) :: Rmat(nch, nch)

  integer, parameter :: sp = kind(1.0e0)
  complex(sp), allocatable :: cmat_sp(:,:), X_sp(:,:)
  complex*16, allocatable :: X_dp(:,:)
  integer, allocatable :: IPIV(:)
  integer :: ich, ichp, ir, ntotal, INFO, i, j

  ntotal = nch * nlag
  allocate(cmat_sp(ntotal, ntotal), X_sp(ntotal, nch), X_dp(ntotal, nch), IPIV(ntotal))

  ! Convert to single precision
  do j = 1, ntotal
    do i = 1, ntotal
      cmat_sp(i, j) = cmplx(real(cmat(i, j)), aimag(cmat(i, j)), sp)
    end do
  end do

  X_sp = (0.0, 0.0)
  do ich = 1, nch
    do ir = 1, nlag
      X_sp((ich-1)*nlag + ir, ich) = cmplx(real(B_vector(ir)), aimag(B_vector(ir)), sp)
    end do
  end do

  call CGETRF(ntotal, ntotal, cmat_sp, ntotal, IPIV, INFO)
  if (INFO /= 0) then
    ! Fallback to double precision
    deallocate(cmat_sp, X_sp, X_dp, IPIV)
    call solve_rmatrix_hp(cmat, B_vector, nch, nlag, normfac, Rmat)
    return
  end if

  call CGETRS('N', ntotal, nch, cmat_sp, ntotal, IPIV, X_sp, ntotal, INFO)

  ! Convert back to double
  do ich = 1, nch
    do i = 1, ntotal
      X_dp(i, ich) = dcmplx(real(X_sp(i, ich)), aimag(X_sp(i, ich)))
    end do
  end do

  Rmat = (0.d0, 0.d0)
  do ichp = 1, nch
    do ich = 1, nch
      do ir = 1, nlag
        Rmat(ich, ichp) = Rmat(ich, ichp) + B_vector(ir) * X_dp(ir + (ich-1)*nlag, ichp)
      end do
    end do
  end do
  Rmat = Rmat * normfac

  deallocate(cmat_sp, X_sp, X_dp, IPIV)
end subroutine solve_rmatrix_hp_mixed

subroutine solve_rmatrix_hp_4(cmat, B0, B1, nch, nlag, normfac, R00, R01, R11)
  implicit none
  integer, intent(in) :: nch, nlag
  real*8, intent(in) :: normfac
  complex*16, intent(in) :: cmat(nch*nlag, nch*nlag)
  real*8, intent(in) :: B0(nlag), B1(nlag)
  complex*16, intent(out) :: R00(nch, nch), R01(nch, nch), R11(nch, nch)

  complex*16, allocatable :: A_work(:,:), X_vector(:,:)
  integer, allocatable :: IPIV(:)
  integer :: ich, ichp, ir, ntotal, INFO

  ntotal = nch * nlag
  allocate(A_work(ntotal, ntotal), X_vector(ntotal, 2*nch), IPIV(ntotal))

  A_work = cmat
  X_vector = (0.d0, 0.d0)
  do ich = 1, nch
    do ir = 1, nlag
      X_vector((ich-1)*nlag + ir, ich) = B0(ir)
      X_vector((ich-1)*nlag + ir, ich + nch) = B1(ir)
    end do
  end do

  call ZGESV(ntotal, 2*nch, A_work, ntotal, IPIV, X_vector, ntotal, INFO)
  if (INFO /= 0) stop 'Error in solve_rmatrix_hp_4'

  R00 = (0.d0, 0.d0)
  R01 = (0.d0, 0.d0)
  R11 = (0.d0, 0.d0)

  do ichp = 1, nch
    do ich = 1, nch
      do ir = 1, nlag
        R00(ich, ichp) = R00(ich, ichp) + B0(ir) * X_vector(ir+(ich-1)*nlag, ichp) * normfac
        R01(ich, ichp) = R01(ich, ichp) + B0(ir) * X_vector(ir+(ich-1)*nlag, ichp+nch) * normfac
        R11(ich, ichp) = R11(ich, ichp) + B1(ir) * X_vector(ir+(ich-1)*nlag, ichp+nch) * normfac
      end do
    end do
  end do

  deallocate(A_work, X_vector, IPIV)
end subroutine solve_rmatrix_hp_4

!------------------------------------------------------------------------------
! Matrix inversion using linear solve (more efficient than explicit inversion)
!------------------------------------------------------------------------------
subroutine cminv_sym_hp(c, n, mdim)
  implicit complex*16(c)
  dimension c(mdim, mdim)
  integer :: IPIV(n), INFO

  complex*16, allocatable :: identity(:,:)
  allocate(identity(n, n))

  identity = (0.d0, 0.d0)
  do i = 1, n
    identity(i, i) = (1.d0, 0.d0)
  end do

  call ZGESV(n, n, c, mdim, IPIV, identity, n, INFO)
  c(1:n, 1:n) = identity

  deallocate(identity)
end subroutine cminv_sym_hp

subroutine cminv_nsym_hp(c, n, mdim)
  implicit complex*16(c)
  dimension c(mdim, mdim)
  integer :: IPIV(n), INFO

  complex*16, allocatable :: identity(:,:)
  allocate(identity(n, n))

  identity = (0.d0, 0.d0)
  do i = 1, n
    identity(i, i) = (1.d0, 0.d0)
  end do

  call ZGESV(n, n, c, mdim, IPIV, identity, n, INFO)
  c(1:n, 1:n) = identity

  deallocate(identity)
end subroutine cminv_nsym_hp

!------------------------------------------------------------------------------
! wf_print: Print wave function on uniform mesh
!------------------------------------------------------------------------------
subroutine wf_print(nch, lval, qk, eta, rmax, nr, ns, cu, &
                    ndim, nopen, cf, nwf1, nwf2, zrma, iv, nom, npoin, h, cwftab)
  implicit none
  integer, intent(in) :: nch, nr, ns, ndim, nopen, nwf1, nwf2, nom, iv, npoin
  integer, intent(in) :: lval(nch)
  real*8, intent(in) :: qk(nch), eta(nch), rmax, zrma(ns*nr), h
  complex*16, intent(in) :: cu(ndim, ndim), cf(nwf1, nwf2, nom)
  complex*16, intent(out) :: cwftab(npoin)

  integer :: nsr, nop, nclo, iw, i, ll, ifail, jw
  real*8 :: r, wfr, wfi, xl, yp1, yp2
  real*8 :: xfr(ns*nr), xfi(ns*nr), y2r(ns*nr), y2i(ns*nr)
  real*8 :: fc(500), dfc(500), gc(500), dgc(500)
  complex*16 :: co, cfx

  yp1 = 1.0d30
  yp2 = 1.0d30

  nsr = ns * nr
  nop = 0
  nclo = nopen

  do iw = 1, iv
    if (qk(iw) > 0) nop = nop + 1
    if (qk(iw) < 0) nclo = nclo + 1
  end do

  xfr(1:nsr) = real(cf(1:nsr, iv, nom))
  xfi(1:nsr) = aimag(cf(1:nsr, iv, nom))

  call spline(zrma, xfr, nsr, yp1, yp2, y2r)
  call spline(zrma, xfi, nsr, yp1, yp2, y2i)

  do i = 1, npoin
    r = i * h
    if (r <= rmax) then
      call splint(zrma, xfr, y2r, nsr, r, wfr)
      call splint(zrma, xfi, y2i, nsr, r, wfi)
      cwftab(i) = dcmplx(wfr, wfi)
    else
      ll = lval(iv)
      if (qk(iv) > 0) then
        xl = lval(iv)
        call coulfg(qk(iv)*r, eta(iv), xl, xl, fc, gc, dfc, dgc, 1, 0, ifail)
        co = dcmplx(gc(ll+1), fc(ll+1)) * sqrt(qk(nom)/qk(iv))
        cfx = -cu(nop, nom) * co
        if (iv == nom) cfx = cfx + conjg(co)
      else
        jw = 0
        call whit(eta(iv), r, -qk(iv), -qk(iv)**2, ll, fc, dfc, jw)
        cfx = cu(nclo, nom) * fc(ll+1)
      end if
      cwftab(i) = cfx
    end if
  end do

end subroutine wf_print

!------------------------------------------------------------------------------
! Auxiliary routines (from Pierre's original code)
!------------------------------------------------------------------------------

! Legendre-Gauss quadrature
SUBROUTINE LEGZO_HP(N, X, W)
  IMPLICIT REAL*8 (A-H,O-Z)
  DIMENSION X(N), W(N)
  data pi/3.1415926535898d0/, one/1/

  N0 = (N+1) / 2
  DO NR = 1, N0
    Z = COS(pi * (NR - 0.25D0) / (N + 0.5d0))
10  Z0 = Z
    P = 1
    DO I = 1, NR-1
      P = P * (Z - X(I))
    END DO
    F0 = 1
    IF (NR == N0 .AND. N /= 2*INT(N/2)) Z = 0
    F1 = Z
    DO K = 2, N
      PF = (2 - one/K) * Z * F1 - (1 - one/K) * F0
      PD = K * (F1 - Z*PF) / (1 - Z*Z)
      F0 = F1
      F1 = PF
    END DO
    IF (Z == 0) GO TO 40
    FD = PF / P
    Q = 0
    DO I = 1, NR-1
      WP = 1
      DO J = 1, NR-1
        IF (J /= I) WP = WP * (Z - X(J))
      END DO
      Q = Q + WP
    END DO
    GD = (PD - Q*FD) / P
    Z = Z - FD/GD
    IF (ABS(Z - Z0) > ABS(Z) * 1.0D-15) GO TO 10
40  X(NR) = Z
    X(N+1-NR) = -Z
    W(NR) = 2 / ((1 - Z*Z) * PD * PD)
    W(N+1-NR) = W(NR)
  END DO
  x(1:n) = (1 + x(n:1:-1)) / 2
  w(1:n) = w(n:1:-1) / 2
END SUBROUTINE LEGZO_HP

! Spline interpolation routines
SUBROUTINE spline(x, y, n, yp1, ypn, y2)
  implicit real*8(a-h,o-z)
  dimension x(n), y(n), y2(n), u(n)

  if (yp1 > .99d30) then
    y2(1) = 0
    u(1) = 0
  else
    y2(1) = -0.5d0
    u(1) = (3/(x(2)-x(1))) * ((y(2)-y(1))/(x(2)-x(1)) - yp1)
  endif

  do i = 2, n-1
    sig = (x(i) - x(i-1)) / (x(i+1) - x(i-1))
    p = sig * y2(i-1) + 2
    y2(i) = (sig - 1.) / p
    u(i) = (6.*((y(i+1)-y(i))/(x(i+1)-x(i)) - (y(i)-y(i-1))/ &
           (x(i)-x(i-1))) / (x(i+1)-x(i-1)) - sig*u(i-1)) / p
  end do

  if (ypn > .99d30) then
    qn = 0
    un = 0
  else
    qn = 0.5d0
    un = (3/(x(n)-x(n-1))) * (ypn - (y(n)-y(n-1))/(x(n)-x(n-1)))
  endif

  y2(n) = (un - qn*u(n-1)) / (qn*y2(n-1) + 1)
  do k = n-1, 1, -1
    y2(k) = y2(k)*y2(k+1) + u(k)
  end do
END SUBROUTINE spline

SUBROUTINE splint(xa, ya, y2a, n, x, y)
  implicit real*8(a-h,o-z)
  dimension xa(n), y2a(n), ya(n)

  klo = 1
  khi = n
1 if (khi - klo > 1) then
    k = (khi + klo) / 2
    if (xa(k) > x) then
      khi = k
    else
      klo = k
    endif
    goto 1
  endif

  h = xa(khi) - xa(klo)
  if (h == 0) pause 'bad xa input in splint'
  a = (xa(khi) - x) / h
  b = (x - xa(klo)) / h
  y = a*ya(klo) + b*ya(khi) + ((a**3-a)*y2a(klo) + (b**3-b)*y2a(khi)) * h**2/6
END SUBROUTINE splint
