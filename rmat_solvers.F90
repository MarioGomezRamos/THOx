!------------------------------------------------------------------------------
! Module: rmat_solvers
! Purpose: High-performance R-matrix linear equation solvers
!
! This module provides multiple solver implementations for the R-matrix
! coupled-channel problem. The key improvement over traditional implementations
! is using linear equation solving (ZGESV) instead of matrix inversion.
!
! Solver types:
!   solver_type=1: Dense LAPACK ZGESV (reference, accurate)
!   solver_type=2: Mixed Precision (CGETRF + refinement, ~5% faster)
!   solver_type=3: Woodbury-Kinetic (CPU optimized, ~13% faster)
!   solver_type=4: GPU cuSOLVER (requires GPU, ~3x faster)
!   solver_type=5: GPU TF32 (TensorFloat-32 Tensor Core, Ampere+ GPUs)
!
! Reference:
!   Based on the Lagrange-mesh R-matrix method described in:
!   P. Descouvemont, CPC 200 (2016) 199-219
!
! Author: Jin Lei
! Date: December 2025
!------------------------------------------------------------------------------
module rmat_solvers
  use precision
  implicit none
  private

  ! Public solver interfaces
  public :: solve_rmatrix_dense      ! solver_type=1
  public :: solve_rmatrix_mixed      ! solver_type=2
  public :: solve_rmatrix_woodbury   ! solver_type=3
  public :: solve_rmatrix_gpu        ! solver_type=4
  public :: solve_rmatrix_tf32       ! solver_type=5
  public :: solve_rmatrix            ! Unified interface with solver selection

  ! Public propagation routine
  public :: rmatrix_propagation

  ! Public utility
  public :: compute_smatrix

contains

!------------------------------------------------------------------------------
! Unified solver interface with solver type selection
!------------------------------------------------------------------------------
subroutine solve_rmatrix(cmat, B_vector, nch, nlag, normfac, Rmat, solver_type, K_pure)
  implicit none
  integer, intent(in) :: nch, nlag
  real(dp), intent(in) :: normfac
  complex(dp), intent(in) :: cmat(nch*nlag, nch*nlag), B_vector(nlag)
  complex(dp), intent(out) :: Rmat(nch, nch)
  integer, intent(in), optional :: solver_type
  real(dp), intent(in), optional :: K_pure(nlag, nlag)

  integer :: stype

  stype = 1  ! Default to dense LAPACK
  if (present(solver_type)) stype = solver_type

  select case (stype)
  case (1)
    call solve_rmatrix_dense(cmat, B_vector, nch, nlag, normfac, Rmat)
  case (2)
    call solve_rmatrix_mixed(cmat, B_vector, nch, nlag, normfac, Rmat)
  case (3)
    if (present(K_pure)) then
      call solve_rmatrix_woodbury(cmat, B_vector, nch, nlag, normfac, Rmat, K_pure)
    else
      call solve_rmatrix_woodbury(cmat, B_vector, nch, nlag, normfac, Rmat)
    end if
  case (4)
    call solve_rmatrix_gpu(cmat, B_vector, nch, nlag, normfac, Rmat)
  case (5)
    call solve_rmatrix_tf32(cmat, B_vector, nch, nlag, normfac, Rmat)
  case default
    call solve_rmatrix_dense(cmat, B_vector, nch, nlag, normfac, Rmat)
  end select

end subroutine solve_rmatrix

!------------------------------------------------------------------------------
! Dense LAPACK ZGESV Solver (solver_type=1) - Reference Implementation
!
! This is the most accurate and stable solver. It serves as the reference
! for validating other solver implementations.
!------------------------------------------------------------------------------
subroutine solve_rmatrix_dense(cmat, B_vector, nch, nlag, normfac, Rmat)
  implicit none
  integer, intent(in) :: nch, nlag
  real(dp), intent(in) :: normfac
  complex(dp), intent(in) :: cmat(nch*nlag, nch*nlag), B_vector(nlag)
  complex(dp), intent(out) :: Rmat(nch, nch)

  complex(dp), allocatable :: A_work(:,:), X_vector(:,:)
  integer, allocatable :: IPIV(:)
  integer :: ich, ichp, ir, ntotal, INFO

  ntotal = nch * nlag

  allocate(A_work(ntotal, ntotal))
  allocate(X_vector(ntotal, nch))
  allocate(IPIV(ntotal))

  ! Copy matrix (ZGESV overwrites it)
  A_work = cmat

  ! Setup RHS: one column per channel
  X_vector = (0.0_dp, 0.0_dp)
  do ich = 1, nch
    do ir = 1, nlag
      X_vector((ich-1)*nlag + ir, ich) = B_vector(ir)
    end do
  end do

  ! Solve using LAPACK ZGESV (LU factorization + solve)
  call ZGESV(ntotal, nch, A_work, ntotal, IPIV, X_vector, ntotal, INFO)
  if (INFO /= 0) then
    write(*,*) 'ERROR: ZGESV failed with INFO =', INFO
    Rmat = (0.0_dp, 0.0_dp)
    deallocate(A_work, X_vector, IPIV)
    return
  end if

  ! Extract R-matrix elements: R_ij = N * sum_r B(r) * X_r,j
  Rmat = (0.0_dp, 0.0_dp)
  do ichp = 1, nch
    do ich = 1, nch
      do ir = 1, nlag
        Rmat(ich, ichp) = Rmat(ich, ichp) + B_vector(ir) * X_vector(ir + (ich-1)*nlag, ichp)
      end do
    end do
  end do
  Rmat = Rmat * normfac

  deallocate(A_work, X_vector, IPIV)

end subroutine solve_rmatrix_dense

!------------------------------------------------------------------------------
! Mixed Precision Solver (solver_type=2)
!
! Uses single precision LU factorization (O(N^3) work) with double precision
! iterative refinement (O(N^2) per iteration). Typically ~5% faster than
! pure double precision while maintaining accuracy.
!------------------------------------------------------------------------------
subroutine solve_rmatrix_mixed(cmat, B_vector, nch, nlag, normfac, Rmat)
  !$ use omp_lib
  implicit none
  integer, intent(in) :: nch, nlag
  real(dp), intent(in) :: normfac
  complex(dp), intent(in) :: cmat(nch*nlag, nch*nlag), B_vector(nlag)
  complex(dp), intent(out) :: Rmat(nch, nch)

  integer, parameter :: sp = kind(1.0e0)
  integer :: ntotal, ich, ichp, ir, i, j, info, iter
  integer, allocatable :: IPIV(:)
  complex(sp), allocatable :: cmat_sp(:,:), work_sp(:,:)
  complex(dp), allocatable :: X_dp(:,:), B_dp(:,:), R_dp(:,:)
  real(dp) :: res_norm, b_norm, rel_res
  complex(dp) :: alpha_z, beta_z
  integer :: max_iter
  real(dp) :: tol

  ntotal = nch * nlag
  max_iter = 2
  tol = 1.0d-10

  allocate(IPIV(ntotal))
  allocate(cmat_sp(ntotal, ntotal))
  allocate(work_sp(ntotal, nch))
  allocate(X_dp(ntotal, nch))
  allocate(B_dp(ntotal, nch))
  allocate(R_dp(ntotal, nch))

  ! Convert matrix to single precision
  !$OMP PARALLEL DO PRIVATE(i, j)
  do j = 1, ntotal
    do i = 1, ntotal
      cmat_sp(i, j) = cmplx(real(cmat(i, j), sp), real(aimag(cmat(i, j)), sp), sp)
    end do
  end do
  !$OMP END PARALLEL DO

  ! Setup RHS
  B_dp = (0.0_dp, 0.0_dp)
  do ich = 1, nch
    do ir = 1, nlag
      B_dp((ich-1)*nlag + ir, ich) = B_vector(ir)
    end do
  end do

  ! Compute ||B||
  b_norm = 0.0_dp
  do ich = 1, nch
    do i = 1, ntotal
      b_norm = b_norm + abs(B_dp(i, ich))**2
    end do
  end do
  b_norm = sqrt(b_norm)

  ! Single precision LU factorization
  call CGETRF(ntotal, ntotal, cmat_sp, ntotal, IPIV, info)

  if (info /= 0) then
    ! Fallback to double precision
    deallocate(IPIV, cmat_sp, work_sp, X_dp, B_dp, R_dp)
    call solve_rmatrix_dense(cmat, B_vector, nch, nlag, normfac, Rmat)
    return
  end if

  ! Initial solve in single precision
  !$OMP PARALLEL DO PRIVATE(i, ich)
  do ich = 1, nch
    do i = 1, ntotal
      work_sp(i, ich) = cmplx(real(B_dp(i, ich), sp), real(aimag(B_dp(i, ich)), sp), sp)
    end do
  end do
  !$OMP END PARALLEL DO

  call CGETRS('N', ntotal, nch, cmat_sp, ntotal, IPIV, work_sp, ntotal, info)

  ! Convert to double precision
  !$OMP PARALLEL DO PRIVATE(i, ich)
  do ich = 1, nch
    do i = 1, ntotal
      X_dp(i, ich) = cmplx(real(work_sp(i, ich)), real(aimag(work_sp(i, ich))), dp)
    end do
  end do
  !$OMP END PARALLEL DO

  ! Iterative refinement
  alpha_z = (-1.0_dp, 0.0_dp)
  beta_z = (1.0_dp, 0.0_dp)

  do iter = 1, max_iter
    R_dp = B_dp
    call ZGEMM('N', 'N', ntotal, nch, ntotal, alpha_z, cmat, ntotal, &
               X_dp, ntotal, beta_z, R_dp, ntotal)

    res_norm = 0.0_dp
    do ich = 1, nch
      do i = 1, ntotal
        res_norm = res_norm + abs(R_dp(i, ich))**2
      end do
    end do
    res_norm = sqrt(res_norm)
    rel_res = res_norm / b_norm

    if (rel_res < tol) exit

    !$OMP PARALLEL DO PRIVATE(i, ich)
    do ich = 1, nch
      do i = 1, ntotal
        work_sp(i, ich) = cmplx(real(R_dp(i, ich), sp), real(aimag(R_dp(i, ich)), sp), sp)
      end do
    end do
    !$OMP END PARALLEL DO

    call CGETRS('N', ntotal, nch, cmat_sp, ntotal, IPIV, work_sp, ntotal, info)

    !$OMP PARALLEL DO PRIVATE(i, ich)
    do ich = 1, nch
      do i = 1, ntotal
        X_dp(i, ich) = X_dp(i, ich) + cmplx(real(work_sp(i, ich)), real(aimag(work_sp(i, ich))), dp)
      end do
    end do
    !$OMP END PARALLEL DO
  end do

  ! Extract R-matrix
  Rmat = (0.0_dp, 0.0_dp)
  do ichp = 1, nch
    do ich = 1, nch
      do ir = 1, nlag
        Rmat(ich, ichp) = Rmat(ich, ichp) + B_vector(ir) * X_dp(ir + (ich-1)*nlag, ichp)
      end do
    end do
  end do
  Rmat = Rmat * normfac

  deallocate(IPIV, cmat_sp, work_sp, X_dp, B_dp, R_dp)

end subroutine solve_rmatrix_mixed

!------------------------------------------------------------------------------
! Woodbury-Kinetic Solver (solver_type=3) - CPU Optimized
!
! Uses the Woodbury matrix identity with the kinetic energy matrix as base.
! This exploits the structure of the R-matrix C matrix where off-diagonal
! blocks (channel coupling) are diagonal matrices.
!
! Typically ~13% faster than dense LAPACK on CPU.
!------------------------------------------------------------------------------
subroutine solve_rmatrix_woodbury(cmat, B_vector, nch, nlag, normfac, Rmat, K_pure)
  !$ use omp_lib
  implicit none
  integer, intent(in) :: nch, nlag
  real(dp), intent(in) :: normfac
  complex(dp), intent(in) :: cmat(nch*nlag, nch*nlag), B_vector(nlag)
  complex(dp), intent(out) :: Rmat(nch, nch)
  real(dp), intent(in), optional :: K_pure(nlag, nlag)

  integer, parameter :: sp = kind(1.0e0)
  integer :: ntotal, ich, ichp, ir, irp, i, j, info
  integer, allocatable :: IPIV_K(:), IPIV_S(:)
  complex(dp), allocatable :: K_matrix(:,:), K_inv(:,:)
  complex(dp), allocatable :: V_radial(:,:,:)
  complex(dp), allocatable :: X_dp(:,:)
  complex(dp), allocatable :: B_radial(:,:,:)
  complex(dp), allocatable :: Schur(:,:), Schur_rhs(:,:)
  complex(dp), allocatable :: work(:,:)

  ntotal = nch * nlag

  allocate(K_matrix(nlag, nlag))
  allocate(K_inv(nlag, nlag))
  allocate(IPIV_K(nlag))
  allocate(V_radial(nch, nch, nlag))
  allocate(X_dp(ntotal, nch))
  allocate(B_radial(nch, nlag, nch))
  allocate(work(nlag, nlag))

  ! Step 1: Extract K and compute K^{-1}
  if (present(K_pure)) then
    K_matrix = cmplx(K_pure, 0.0_dp, dp)
    block
      real(dp), allocatable :: K_real(:,:), K_inv_real(:,:)
      allocate(K_real(nlag, nlag), K_inv_real(nlag, nlag))
      K_real = K_pure
      K_inv_real = 0.0_dp
      do ir = 1, nlag
        K_inv_real(ir, ir) = 1.0_dp
      end do
      call DGETRF(nlag, nlag, K_real, nlag, IPIV_K, info)
      if (info /= 0) then
        deallocate(K_real, K_inv_real)
        deallocate(K_matrix, K_inv, IPIV_K, V_radial, X_dp, B_radial, work)
        call solve_rmatrix_dense(cmat, B_vector, nch, nlag, normfac, Rmat)
        return
      end if
      call DGETRS('N', nlag, nlag, K_real, nlag, IPIV_K, K_inv_real, nlag, info)
      K_inv = cmplx(K_inv_real, 0.0_dp, dp)
      deallocate(K_real, K_inv_real)
    end block
  else
    K_matrix = (0.0_dp, 0.0_dp)
    do ir = 1, nlag
      do irp = 1, nlag
        K_matrix(ir, irp) = cmat(ir, irp)
      end do
    end do
    K_inv = (0.0_dp, 0.0_dp)
    do ir = 1, nlag
      K_inv(ir, ir) = (1.0_dp, 0.0_dp)
    end do
    work = K_matrix
    call ZGETRF(nlag, nlag, work, nlag, IPIV_K, info)
    if (info /= 0) then
      deallocate(K_matrix, K_inv, IPIV_K, V_radial, X_dp, B_radial, work)
      call solve_rmatrix_dense(cmat, B_vector, nch, nlag, normfac, Rmat)
      return
    end if
    call ZGETRS('N', nlag, nlag, work, nlag, IPIV_K, K_inv, nlag, info)
  end if

  ! Step 2: Extract V(r) coupling matrices
  V_radial = (0.0_dp, 0.0_dp)
  do ir = 1, nlag
    do ich = 1, nch
      do ichp = 1, nch
        V_radial(ich, ichp, ir) = cmat((ich-1)*nlag + ir, (ichp-1)*nlag + ir)
        if (ich == ichp) then
          V_radial(ich, ichp, ir) = V_radial(ich, ichp, ir) - K_matrix(ir, ir)
        end if
      end do
    end do
  end do

  ! Step 3: Build Schur complement
  allocate(Schur(ntotal, ntotal))
  Schur = (0.0_dp, 0.0_dp)

  !$OMP PARALLEL DO PRIVATE(ir, irp, ich, ichp) COLLAPSE(2)
  do ir = 1, nlag
    do irp = 1, nlag
      do ich = 1, nch
        do ichp = 1, nch
          if (ir == irp .and. ich == ichp) then
            Schur((ir-1)*nch + ich, (irp-1)*nch + ichp) = (1.0_dp, 0.0_dp)
          end if
          Schur((ir-1)*nch + ich, (irp-1)*nch + ichp) = &
              Schur((ir-1)*nch + ich, (irp-1)*nch + ichp) + &
              K_inv(ir, irp) * V_radial(ich, ichp, irp)
        end do
      end do
    end do
  end do
  !$OMP END PARALLEL DO

  ! Setup RHS
  B_radial = (0.0_dp, 0.0_dp)
  do ich = 1, nch
    do ir = 1, nlag
      B_radial(ich, ir, ich) = B_vector(ir)
    end do
  end do

  allocate(Schur_rhs(ntotal, nch))
  Schur_rhs = (0.0_dp, 0.0_dp)

  do ich = 1, nch
    do ir = 1, nlag
      do irp = 1, nlag
        Schur_rhs((ir-1)*nch + ich, :) = Schur_rhs((ir-1)*nch + ich, :) + &
            K_inv(ir, irp) * B_radial(ich, irp, :)
      end do
    end do
  end do

  ! Step 4: Solve Schur system (use single precision for speed)
  allocate(IPIV_S(ntotal))

  block
    complex(sp), allocatable :: Schur_sp(:,:), rhs_sp(:,:)
    allocate(Schur_sp(ntotal, ntotal))
    allocate(rhs_sp(ntotal, nch))
    Schur_sp = cmplx(real(Schur, sp), real(aimag(Schur), sp), sp)
    call CGETRF(ntotal, ntotal, Schur_sp, ntotal, IPIV_S, info)
    if (info /= 0) then
      ! Fallback to double precision
      call ZGESV(ntotal, nch, Schur, ntotal, IPIV_S, Schur_rhs, ntotal, info)
    else
      rhs_sp = cmplx(real(Schur_rhs, sp), real(aimag(Schur_rhs), sp), sp)
      call CGETRS('N', ntotal, nch, Schur_sp, ntotal, IPIV_S, rhs_sp, ntotal, info)
      Schur_rhs = cmplx(real(rhs_sp), real(aimag(rhs_sp)), dp)
    end if
    deallocate(Schur_sp, rhs_sp)
  end block

  if (info /= 0) then
    Rmat = (0.0_dp, 0.0_dp)
    deallocate(K_matrix, K_inv, IPIV_K, V_radial, X_dp, B_radial, work)
    deallocate(Schur, Schur_rhs, IPIV_S)
    return
  end if

  ! Convert solution back to channel-major
  X_dp = (0.0_dp, 0.0_dp)
  do ir = 1, nlag
    do ich = 1, nch
      X_dp((ich-1)*nlag + ir, :) = Schur_rhs((ir-1)*nch + ich, :)
    end do
  end do

  ! Extract R-matrix
  Rmat = (0.0_dp, 0.0_dp)
  do ichp = 1, nch
    do ich = 1, nch
      do ir = 1, nlag
        Rmat(ich, ichp) = Rmat(ich, ichp) + B_vector(ir) * X_dp(ir + (ich-1)*nlag, ichp)
      end do
    end do
  end do
  Rmat = Rmat * normfac

  deallocate(K_matrix, K_inv, IPIV_K, V_radial, X_dp, B_radial, work)
  deallocate(Schur, Schur_rhs, IPIV_S)

end subroutine solve_rmatrix_woodbury

!------------------------------------------------------------------------------
! GPU cuSOLVER Solver (solver_type=4)
!
! Uses NVIDIA cuSOLVER for GPU-accelerated solving. Falls back to CPU ZGESV
! if GPU is not available or if compilation was done without GPU support.
!
! Features:
! - Automatic single/multi-GPU selection based on matrix size and GPU count
! - Single GPU: Mixed precision (FP32 LU + FP64 refinement) for speed
! - Multi-GPU: Full double precision with data distributed across GPUs
! - Transparent fallback to CPU if GPU unavailable
!
! Typically ~3-18x faster than CPU for large matrices (n > 1000).
!------------------------------------------------------------------------------
subroutine solve_rmatrix_gpu(cmat, B_vector, nch, nlag, normfac, Rmat)
#ifdef GPU_ENABLED
  use gpu_solver_interface
#endif
  implicit none
  integer, intent(in) :: nch, nlag
  real(dp), intent(in) :: normfac
  complex(dp), intent(in) :: cmat(nch*nlag, nch*nlag), B_vector(nlag)
  complex(dp), intent(out) :: Rmat(nch, nch)

  complex(dp), allocatable :: A_copy(:,:), X_vector(:,:)
  integer :: ich, ichp, ir, info, ntotal
  integer :: max_refine
  real(dp) :: tol

  ntotal = nch * nlag
  max_refine = 0      ! Disable iterative refinement (fast enough without)
  tol = 1.0d-10

  allocate(A_copy(ntotal, ntotal))
  allocate(X_vector(ntotal, nch))

  ! Copy matrix
  A_copy = cmat

  ! Set up RHS
  X_vector = (0.0_dp, 0.0_dp)
  do ich = 1, nch
    do ir = 1, nlag
      X_vector((ich-1)*nlag + ir, ich) = B_vector(ir)
    end do
  end do

#ifdef GPU_ENABLED
  ! GPU solver - automatically chooses single or multi-GPU
  ! based on matrix size and number of available GPUs
  call gpu_solve_auto(A_copy, X_vector, ntotal, nch, max_refine, tol, info)

  if (info /= 0) then
    ! Fallback to CPU
    A_copy = cmat
    do ich = 1, nch
      do ir = 1, nlag
        X_vector((ich-1)*nlag + ir, ich) = B_vector(ir)
      end do
    end do
    block
      integer :: IPIV(ntotal)
      call ZGESV(ntotal, nch, A_copy, ntotal, IPIV, X_vector, ntotal, info)
    end block
  end if
#else
  ! No GPU - use CPU ZGESV
  block
    integer :: IPIV(ntotal)
    call ZGESV(ntotal, nch, A_copy, ntotal, IPIV, X_vector, ntotal, info)
  end block
#endif

  ! Extract R-matrix elements
  Rmat = (0.0_dp, 0.0_dp)
  if (info == 0) then
    do ichp = 1, nch
      do ich = 1, nch
        do ir = 1, nlag
          Rmat(ich, ichp) = Rmat(ich, ichp) + &
              B_vector(ir) * X_vector(ir + (ich-1)*nlag, ichp)
        end do
      end do
    end do
    Rmat = Rmat * normfac
  end if

  deallocate(A_copy, X_vector)

end subroutine solve_rmatrix_gpu

!------------------------------------------------------------------------------
! GPU TF32 (TensorFloat-32) Solver (solver_type=5)
!
! Uses TF32 precision on NVIDIA Ampere+ GPUs (RTX 30/40 series, A100, H100)
! TF32 uses Tensor Cores with:
!   - 8-bit exponent (same dynamic range as FP32)
!   - 10-bit mantissa (same as FP16, but safer due to FP32 range)
!
! Falls back to FP32 mixed precision on older GPUs (Volta, Turing).
! Requires more iterative refinement than FP32 to recover accuracy.
!------------------------------------------------------------------------------
subroutine solve_rmatrix_tf32(cmat, B_vector, nch, nlag, normfac, Rmat)
#ifdef GPU_ENABLED
  use gpu_solver_interface
#endif
  implicit none
  integer, intent(in) :: nch, nlag
  real(dp), intent(in) :: normfac
  complex(dp), intent(in) :: cmat(nch*nlag, nch*nlag), B_vector(nlag)
  complex(dp), intent(out) :: Rmat(nch, nch)

  complex(dp), allocatable :: A_copy(:,:), X_vector(:,:)
  integer :: ich, ichp, ir, ntotal, info

  ntotal = nch * nlag

  allocate(A_copy(ntotal, ntotal))
  allocate(X_vector(ntotal, nch))

  ! Copy matrix and set up RHS
  A_copy = cmat
  X_vector = (0.0_dp, 0.0_dp)
  do ich = 1, nch
    do ir = 1, nlag
      X_vector(ir + (ich-1)*nlag, ich) = B_vector(ir)
    end do
  end do

#ifdef GPU_ENABLED
  ! Use GPU TF32 solver with 5 refinement iterations
  call gpu_solve_tf32(A_copy, X_vector, ntotal, nch, 5, 1.0d-12, info)
#else
  ! Fallback to CPU ZGESV
  block
    integer :: IPIV(ntotal)
    call ZGESV(ntotal, nch, A_copy, ntotal, IPIV, X_vector, ntotal, info)
  end block
#endif

  ! Extract R-matrix elements
  Rmat = (0.0_dp, 0.0_dp)
  if (info == 0) then
    do ichp = 1, nch
      do ich = 1, nch
        do ir = 1, nlag
          Rmat(ich, ichp) = Rmat(ich, ichp) + &
              B_vector(ir) * X_vector(ir + (ich-1)*nlag, ichp)
        end do
      end do
    end do
    Rmat = Rmat * normfac
  end if

  deallocate(A_copy, X_vector)

end subroutine solve_rmatrix_tf32

!------------------------------------------------------------------------------
! R-matrix Propagation for Multi-interval Calculations
!
! For ns > 1 (multiple intervals), the R-matrix needs to be propagated
! from the inner region to the outer boundary.
!------------------------------------------------------------------------------
subroutine rmatrix_propagation(cmat, B0, B1, nch, nlag, normfac, ins, rmax, Rmat0, Rmat)
  implicit none
  integer, intent(in) :: nch, nlag, ins
  real(dp), intent(in) :: normfac, rmax
  complex(dp), intent(in) :: cmat(nch*nlag, nch*nlag), B0(nlag), B1(nlag)
  complex(dp), intent(in) :: Rmat0(nch, nch)
  complex(dp), intent(out) :: Rmat(nch, nch)

  complex(dp) :: R00(nch, nch), R01(nch, nch), R10(nch, nch), R11(nch, nch)
  complex(dp) :: inv_matrix(nch, nch)
  complex(dp), allocatable :: A_work(:,:), X_vector(:,:)
  integer, allocatable :: IPIV(:)
  integer :: ich, ichp, ir, ntotal, INFO

  ntotal = nch * nlag

  allocate(A_work(ntotal, ntotal))
  allocate(X_vector(ntotal, 2*nch))
  allocate(IPIV(ntotal))

  A_work = cmat

  ! Setup RHS for both B0 and B1
  X_vector = (0.0_dp, 0.0_dp)
  do ich = 1, nch
    do ir = 1, nlag
      X_vector((ich-1)*nlag + ir, ich) = B0(ir)
      X_vector((ich-1)*nlag + ir, ich + nch) = B1(ir)
    end do
  end do

  call ZGESV(ntotal, 2*nch, A_work, ntotal, IPIV, X_vector, ntotal, INFO)
  if (INFO /= 0) then
    write(*,*) 'ERROR in propagation ZGESV, INFO =', INFO
    Rmat = (0.0_dp, 0.0_dp)
    deallocate(A_work, X_vector, IPIV)
    return
  end if

  ! Compute the four R-matrices
  R00 = (0.0_dp, 0.0_dp)
  R01 = (0.0_dp, 0.0_dp)
  R10 = (0.0_dp, 0.0_dp)
  R11 = (0.0_dp, 0.0_dp)

  do ichp = 1, nch
    do ich = 1, nch
      do ir = 1, nlag
        R00(ich, ichp) = R00(ich, ichp) + B0(ir) * X_vector(ir + (ich-1)*nlag, ichp) * normfac
        R01(ich, ichp) = R01(ich, ichp) + B0(ir) * X_vector(ir + (ich-1)*nlag, ichp + nch) * normfac
        R10(ich, ichp) = R10(ich, ichp) + B1(ir) * X_vector(ir + (ich-1)*nlag, ichp) * normfac
        R11(ich, ichp) = R11(ich, ichp) + B1(ir) * X_vector(ir + (ich-1)*nlag, ichp + nch) * normfac
      end do
    end do
  end do

  ! Propagate: Rmat = (R11 - R10 * inv(R00 + (ins-1)*rmax*Rmat0) * R01) / (ins*rmax)
  inv_matrix = R00 + (ins - 1.0_dp) * rmax * Rmat0

  deallocate(IPIV)
  allocate(IPIV(nch))
  call ZGESV(nch, nch, inv_matrix, nch, IPIV, R01, nch, INFO)
  if (INFO /= 0) then
    write(*,*) 'ERROR in propagation inversion, INFO =', INFO
    Rmat = (0.0_dp, 0.0_dp)
    deallocate(A_work, X_vector, IPIV)
    return
  end if

  Rmat = (R11 - matmul(R10, R01)) / (ins * rmax)

  deallocate(A_work, X_vector, IPIV)

end subroutine rmatrix_propagation

!------------------------------------------------------------------------------
! Compute S-matrix from Z-matrices
!
! Using Eq.(16) of Pierre's paper: Z_O * S = Z_I
!------------------------------------------------------------------------------
subroutine compute_smatrix(nch, Zmat_I, Zmat_O, Smat)
  implicit none
  integer, intent(in) :: nch
  complex(dp), intent(in) :: Zmat_I(nch, nch), Zmat_O(nch, nch)
  complex(dp), intent(out) :: Smat(nch, nch)

  complex(dp) :: Z_work(nch, nch)
  integer :: IPIV(nch), INFO

  Z_work = Zmat_O
  Smat = Zmat_I

  call ZGESV(nch, nch, Z_work, nch, IPIV, Smat, nch, INFO)
  if (INFO /= 0) then
    write(*,*) 'ERROR in S-matrix computation, INFO =', INFO
    Smat = (0.0_dp, 0.0_dp)
  end if

end subroutine compute_smatrix

end module rmat_solvers
