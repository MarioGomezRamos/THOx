/*
 * cusolver_interface.cu - Optimized GPU solver for HPRMAT
 *
 * Uses FP32 (single precision) for LU factorization on GPU
 * with optional FP64 iterative refinement for accuracy.
 *
 * Key optimizations:
 * - FP64→FP32 conversion on GPU using CUDA kernels
 * - Persistent memory allocation (reused across calls)
 * - Minimized CPU-GPU data transfers
 * - Multi-GPU support via cusolverMg for large matrices
 *
 * Multi-GPU notes:
 * - Automatically detects number of GPUs
 * - Falls back to single-GPU if only 1 GPU available
 * - Uses 1D block-cyclic distribution across GPUs
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <cusolverMg.h>
#include <cublas_v2.h>
#include <cuComplex.h>

// Error checking macros
#define CUDA_CHECK(call) do { \
    cudaError_t err = call; \
    if (err != cudaSuccess) { \
        fprintf(stderr, "CUDA error at %s:%d: %s\n", __FILE__, __LINE__, \
                cudaGetErrorString(err)); \
        return -1; \
    } \
} while(0)

#define CUSOLVER_CHECK(call) do { \
    cusolverStatus_t status = call; \
    if (status != CUSOLVER_STATUS_SUCCESS) { \
        fprintf(stderr, "cuSOLVER error at %s:%d: %d\n", __FILE__, __LINE__, status); \
        return -1; \
    } \
} while(0)

#define CUBLAS_CHECK(call) do { \
    cublasStatus_t status = call; \
    if (status != CUBLAS_STATUS_SUCCESS) { \
        fprintf(stderr, "cuBLAS error at %s:%d: %d\n", __FILE__, __LINE__, status); \
        return -1; \
    } \
} while(0)

// Global handles - single GPU
static cusolverDnHandle_t cusolver_handle = NULL;
static cublasHandle_t cublas_handle = NULL;
static int gpu_initialized = 0;
static int current_device = 0;

// Global handles - multi-GPU
static cusolverMgHandle_t cusolvermg_handle = NULL;
static int num_gpus = 0;
static int *device_ids = NULL;
static int multi_gpu_initialized = 0;

// Multi-GPU grid descriptor
static cudaLibMgMatrixDesc_t mg_desc_A = NULL;
static cudaLibMgMatrixDesc_t mg_desc_B = NULL;
static cudaLibMgGrid_t mg_grid = NULL;

// Multi-GPU memory arrays (one per GPU)
static void **d_A_mg = NULL;
static void **d_B_mg = NULL;
static int **d_ipiv_mg = NULL;  // Array of pivot arrays, one per GPU
static void **d_work_mg = NULL; // Array of workspaces, one per GPU
static int64_t mg_lwork = 0;

// Block size for multi-GPU distribution (must be power of 2, typically 256-1024)
#define MG_BLOCK_SIZE 256

// Persistent GPU memory
static cuComplex *d_A_sp = NULL;       // Single precision matrix
static cuComplex *d_B_sp = NULL;       // Single precision RHS
static cuDoubleComplex *d_A_dp = NULL; // Double precision matrix (for refinement)
static cuDoubleComplex *d_B_dp = NULL; // Double precision RHS
static cuDoubleComplex *d_R = NULL;    // Residual
static cuDoubleComplex *d_X = NULL;    // Solution
static int *d_ipiv = NULL;
static int *d_info = NULL;
static void *d_work = NULL;
static size_t work_size = 0;
static int alloc_n = 0;
static int alloc_nrhs = 0;

// ============================================
// CUDA Kernels for FP64 <-> FP32 conversion
// ============================================

// Convert complex double to complex float (on GPU)
__global__ void convert_z2c_kernel(const cuDoubleComplex* __restrict__ src,
                                    cuComplex* __restrict__ dst,
                                    int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        dst[idx] = make_cuFloatComplex((float)src[idx].x, (float)src[idx].y);
    }
}

// Convert complex float to complex double (on GPU)
__global__ void convert_c2z_kernel(const cuComplex* __restrict__ src,
                                    cuDoubleComplex* __restrict__ dst,
                                    int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        dst[idx] = make_cuDoubleComplex((double)src[idx].x, (double)src[idx].y);
    }
}

// Compute residual R = B - A*X on GPU
__global__ void compute_residual_kernel(const cuDoubleComplex* __restrict__ A,
                                         const cuDoubleComplex* __restrict__ X,
                                         const cuDoubleComplex* __restrict__ B,
                                         cuDoubleComplex* __restrict__ R,
                                         int n, int nrhs) {
    int row = blockIdx.x * blockDim.x + threadIdx.x;
    int col = blockIdx.y;

    if (row < n && col < nrhs) {
        cuDoubleComplex sum = make_cuDoubleComplex(0.0, 0.0);
        for (int k = 0; k < n; k++) {
            cuDoubleComplex a = A[row + k * n];
            cuDoubleComplex x = X[k + col * n];
            sum.x += a.x * x.x - a.y * x.y;
            sum.y += a.x * x.y + a.y * x.x;
        }
        R[row + col * n].x = B[row + col * n].x - sum.x;
        R[row + col * n].y = B[row + col * n].y - sum.y;
    }
}

// Add correction: X += dX
__global__ void add_correction_kernel(cuDoubleComplex* __restrict__ X,
                                       const cuComplex* __restrict__ dX_sp,
                                       int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        X[idx].x += (double)dX_sp[idx].x;
        X[idx].y += (double)dX_sp[idx].y;
    }
}

// ============================================
// Initialization and cleanup
// ============================================

extern "C" {

int gpu_solver_init_(int *device_id) {
    if (gpu_initialized) {
        return 0;
    }

    current_device = *device_id;
    CUDA_CHECK(cudaSetDevice(current_device));

    cudaDeviceProp prop;
    CUDA_CHECK(cudaGetDeviceProperties(&prop, current_device));
    printf("GPU Solver initialized on: %s (%.1f GB)\n",
           prop.name, prop.totalGlobalMem / 1e9);

    CUSOLVER_CHECK(cusolverDnCreate(&cusolver_handle));
    CUBLAS_CHECK(cublasCreate(&cublas_handle));

    gpu_initialized = 1;
    return 0;
}

int gpu_solver_finalize_() {
    if (!gpu_initialized) {
        return 0;
    }

    if (d_A_sp) cudaFree(d_A_sp);
    if (d_B_sp) cudaFree(d_B_sp);
    if (d_A_dp) cudaFree(d_A_dp);
    if (d_B_dp) cudaFree(d_B_dp);
    if (d_R) cudaFree(d_R);
    if (d_X) cudaFree(d_X);
    if (d_ipiv) cudaFree(d_ipiv);
    if (d_info) cudaFree(d_info);
    if (d_work) cudaFree(d_work);

    d_A_sp = d_B_sp = NULL;
    d_A_dp = d_B_dp = d_R = d_X = NULL;
    d_ipiv = d_info = NULL;
    d_work = NULL;
    alloc_n = alloc_nrhs = 0;

    if (cusolver_handle) cusolverDnDestroy(cusolver_handle);
    if (cublas_handle) cublasDestroy(cublas_handle);

    cusolver_handle = NULL;
    cublas_handle = NULL;
    gpu_initialized = 0;

    return 0;
}

static int ensure_gpu_memory(int n, int nrhs) {
    if (n <= alloc_n && nrhs <= alloc_nrhs) {
        return 0;
    }

    // Free old allocations
    if (d_A_sp) cudaFree(d_A_sp);
    if (d_B_sp) cudaFree(d_B_sp);
    if (d_A_dp) cudaFree(d_A_dp);
    if (d_B_dp) cudaFree(d_B_dp);
    if (d_R) cudaFree(d_R);
    if (d_X) cudaFree(d_X);
    if (d_ipiv) cudaFree(d_ipiv);
    if (d_info) cudaFree(d_info);
    if (d_work) cudaFree(d_work);

    // Allocate new memory
    CUDA_CHECK(cudaMalloc(&d_A_sp, n * n * sizeof(cuComplex)));
    CUDA_CHECK(cudaMalloc(&d_B_sp, n * nrhs * sizeof(cuComplex)));
    CUDA_CHECK(cudaMalloc(&d_A_dp, n * n * sizeof(cuDoubleComplex)));
    CUDA_CHECK(cudaMalloc(&d_B_dp, n * nrhs * sizeof(cuDoubleComplex)));
    CUDA_CHECK(cudaMalloc(&d_R, n * nrhs * sizeof(cuDoubleComplex)));
    CUDA_CHECK(cudaMalloc(&d_X, n * nrhs * sizeof(cuDoubleComplex)));
    CUDA_CHECK(cudaMalloc(&d_ipiv, n * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_info, sizeof(int)));

    // Get workspace size for CGETRF
    int lwork;
    CUSOLVER_CHECK(cusolverDnCgetrf_bufferSize(cusolver_handle, n, n, d_A_sp, n, &lwork));

    work_size = lwork * sizeof(cuComplex);
    CUDA_CHECK(cudaMalloc(&d_work, work_size));

    alloc_n = n;
    alloc_nrhs = nrhs;

    return 0;
}

// ============================================
// Optimized Mixed Precision Solver
// ============================================

int gpu_solve_mixed_(
    double *h_A,
    double *h_B,
    int *n_ptr,
    int *nrhs_ptr,
    int *max_refine_ptr,
    double *tol_ptr,
    int *info_ptr
) {
    int n = *n_ptr;
    int nrhs = *nrhs_ptr;
    int max_refine = *max_refine_ptr;
    double tol = *tol_ptr;
    int h_info;

    if (!gpu_initialized) {
        int device = 0;
        if (gpu_solver_init_(&device) != 0) {
            *info_ptr = -1;
            return -1;
        }
    }

    if (ensure_gpu_memory(n, nrhs) != 0) {
        *info_ptr = -2;
        return -1;
    }

    cuDoubleComplex *h_A_dp = (cuDoubleComplex*)h_A;
    cuDoubleComplex *h_B_dp = (cuDoubleComplex*)h_B;

    // Copy FP64 data to GPU
    CUDA_CHECK(cudaMemcpy(d_A_dp, h_A_dp, n * n * sizeof(cuDoubleComplex),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_B_dp, h_B_dp, n * nrhs * sizeof(cuDoubleComplex),
                          cudaMemcpyHostToDevice));

    // Convert FP64 -> FP32 on GPU (fast!)
    int blockSize = 256;
    int numBlocks_A = (n * n + blockSize - 1) / blockSize;
    int numBlocks_B = (n * nrhs + blockSize - 1) / blockSize;

    convert_z2c_kernel<<<numBlocks_A, blockSize>>>(d_A_dp, d_A_sp, n * n);
    convert_z2c_kernel<<<numBlocks_B, blockSize>>>(d_B_dp, d_B_sp, n * nrhs);
    CUDA_CHECK(cudaGetLastError());

    // LU factorization in FP32
    CUSOLVER_CHECK(cusolverDnCgetrf(cusolver_handle, n, n, d_A_sp, n,
                                     (cuComplex*)d_work, d_ipiv, d_info));
    CUDA_CHECK(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        fprintf(stderr, "CGETRF failed with info = %d\n", h_info);
        *info_ptr = h_info;
        return -1;
    }

    // Solve in FP32
    CUSOLVER_CHECK(cusolverDnCgetrs(cusolver_handle, CUBLAS_OP_N, n, nrhs,
                                     d_A_sp, n, d_ipiv, d_B_sp, n, d_info));

    // Convert solution FP32 -> FP64 on GPU
    convert_c2z_kernel<<<numBlocks_B, blockSize>>>(d_B_sp, d_X, n * nrhs);
    CUDA_CHECK(cudaGetLastError());

    // Iterative refinement on GPU (if requested)
    if (max_refine > 0) {
        cuDoubleComplex alpha = make_cuDoubleComplex(-1.0, 0.0);
        cuDoubleComplex beta = make_cuDoubleComplex(1.0, 0.0);

        for (int iter = 0; iter < max_refine; iter++) {
            // Compute residual R = B - A*X on GPU using cuBLAS
            CUDA_CHECK(cudaMemcpy(d_R, d_B_dp, n * nrhs * sizeof(cuDoubleComplex),
                                  cudaMemcpyDeviceToDevice));
            CUBLAS_CHECK(cublasZgemm(cublas_handle, CUBLAS_OP_N, CUBLAS_OP_N,
                                      n, nrhs, n,
                                      &alpha, d_A_dp, n, d_X, n,
                                      &beta, d_R, n));

            // Check convergence (compute max norm on GPU)
            int max_idx;
            CUBLAS_CHECK(cublasIzamax(cublas_handle, n * nrhs, d_R, 1, &max_idx));
            cuDoubleComplex max_val;
            CUDA_CHECK(cudaMemcpy(&max_val, d_R + max_idx - 1, sizeof(cuDoubleComplex),
                                  cudaMemcpyDeviceToHost));
            double max_res = sqrt(max_val.x * max_val.x + max_val.y * max_val.y);

            if (max_res < tol) {
                break;
            }

            // Convert residual to FP32 and solve for correction
            convert_z2c_kernel<<<numBlocks_B, blockSize>>>(d_R, d_B_sp, n * nrhs);
            CUSOLVER_CHECK(cusolverDnCgetrs(cusolver_handle, CUBLAS_OP_N, n, nrhs,
                                             d_A_sp, n, d_ipiv, d_B_sp, n, d_info));

            // Add correction: X += dX (on GPU)
            add_correction_kernel<<<numBlocks_B, blockSize>>>(d_X, d_B_sp, n * nrhs);
            CUDA_CHECK(cudaGetLastError());
        }
    }

    // Copy solution back to host
    CUDA_CHECK(cudaMemcpy(h_B_dp, d_X, n * nrhs * sizeof(cuDoubleComplex),
                          cudaMemcpyDeviceToHost));

    *info_ptr = 0;
    return 0;
}

// ============================================
// TF32 (TensorFloat-32) Mixed Precision Solver
// Uses Tensor Core with TF32 precision for faster LU
// Requires Ampere or newer GPU (sm_80+)
// ============================================

int gpu_solve_tf32_(
    double *h_A,
    double *h_B,
    int *n_ptr,
    int *nrhs_ptr,
    int *max_refine_ptr,
    double *tol_ptr,
    int *info_ptr
) {
    int n = *n_ptr;
    int nrhs = *nrhs_ptr;
    int max_refine = *max_refine_ptr;
    double tol = *tol_ptr;
    int h_info;

    if (!gpu_initialized) {
        int device = 0;
        if (gpu_solver_init_(&device) != 0) {
            *info_ptr = -1;
            return -1;
        }
    }

    // Check if GPU supports TF32 (Ampere sm_80+)
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, current_device);
    int sm_version = prop.major * 10 + prop.minor;

    if (sm_version < 80) {
        // Fall back to regular FP32 mixed precision on older GPUs
        printf("TF32 not supported on sm_%d, using FP32 mixed precision\n", sm_version);
        return gpu_solve_mixed_(h_A, h_B, n_ptr, nrhs_ptr, max_refine_ptr, tol_ptr, info_ptr);
    }

    if (ensure_gpu_memory(n, nrhs) != 0) {
        *info_ptr = -2;
        return -1;
    }

    cuDoubleComplex *h_A_dp = (cuDoubleComplex*)h_A;
    cuDoubleComplex *h_B_dp = (cuDoubleComplex*)h_B;

    // Enable TF32 Tensor Core math mode for cuBLAS
    // This affects GEMM operations used internally by cuSOLVER
    cublasSetMathMode(cublas_handle, CUBLAS_TF32_TENSOR_OP_MATH);

    // Copy FP64 data to GPU
    CUDA_CHECK(cudaMemcpy(d_A_dp, h_A_dp, n * n * sizeof(cuDoubleComplex),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_B_dp, h_B_dp, n * nrhs * sizeof(cuDoubleComplex),
                          cudaMemcpyHostToDevice));

    // Convert FP64 -> FP32 on GPU
    int blockSize = 256;
    int numBlocks_A = (n * n + blockSize - 1) / blockSize;
    int numBlocks_B = (n * nrhs + blockSize - 1) / blockSize;

    convert_z2c_kernel<<<numBlocks_A, blockSize>>>(d_A_dp, d_A_sp, n * n);
    convert_z2c_kernel<<<numBlocks_B, blockSize>>>(d_B_dp, d_B_sp, n * nrhs);
    CUDA_CHECK(cudaGetLastError());

    // LU factorization in FP32 with TF32 Tensor Core acceleration
    // cuSOLVER uses cuBLAS GEMM internally, which now uses TF32
    CUSOLVER_CHECK(cusolverDnCgetrf(cusolver_handle, n, n, d_A_sp, n,
                                     (cuComplex*)d_work, d_ipiv, d_info));
    CUDA_CHECK(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        fprintf(stderr, "TF32 CGETRF failed with info = %d\n", h_info);
        cublasSetMathMode(cublas_handle, CUBLAS_DEFAULT_MATH);
        *info_ptr = h_info;
        return -1;
    }

    // Solve in FP32/TF32
    CUSOLVER_CHECK(cusolverDnCgetrs(cusolver_handle, CUBLAS_OP_N, n, nrhs,
                                     d_A_sp, n, d_ipiv, d_B_sp, n, d_info));

    // Convert solution FP32 -> FP64 on GPU
    convert_c2z_kernel<<<numBlocks_B, blockSize>>>(d_B_sp, d_X, n * nrhs);
    CUDA_CHECK(cudaGetLastError());

    // Iterative refinement is more important for TF32 due to lower precision
    // TF32 has ~10-bit mantissa vs FP32's 23-bit, so we need more refinement
    int actual_refine = (max_refine > 0) ? max_refine : 5;  // Default to 5 iterations for TF32

    cuDoubleComplex alpha = make_cuDoubleComplex(-1.0, 0.0);
    cuDoubleComplex beta = make_cuDoubleComplex(1.0, 0.0);

    for (int iter = 0; iter < actual_refine; iter++) {
        // Compute residual R = B - A*X on GPU using cuBLAS (in FP64 for accuracy)
        CUDA_CHECK(cudaMemcpy(d_R, d_B_dp, n * nrhs * sizeof(cuDoubleComplex),
                              cudaMemcpyDeviceToDevice));

        // Use default math mode for residual computation (FP64)
        cublasSetMathMode(cublas_handle, CUBLAS_DEFAULT_MATH);
        CUBLAS_CHECK(cublasZgemm(cublas_handle, CUBLAS_OP_N, CUBLAS_OP_N,
                                  n, nrhs, n,
                                  &alpha, d_A_dp, n, d_X, n,
                                  &beta, d_R, n));

        // Check convergence
        int max_idx;
        CUBLAS_CHECK(cublasIzamax(cublas_handle, n * nrhs, d_R, 1, &max_idx));
        cuDoubleComplex max_val;
        CUDA_CHECK(cudaMemcpy(&max_val, d_R + max_idx - 1, sizeof(cuDoubleComplex),
                              cudaMemcpyDeviceToHost));
        double max_res = sqrt(max_val.x * max_val.x + max_val.y * max_val.y);

        if (max_res < tol) {
            break;
        }

        // Convert residual to FP32 and solve for correction
        convert_z2c_kernel<<<numBlocks_B, blockSize>>>(d_R, d_B_sp, n * nrhs);

        // Back to TF32 mode for the solve
        cublasSetMathMode(cublas_handle, CUBLAS_TF32_TENSOR_OP_MATH);
        CUSOLVER_CHECK(cusolverDnCgetrs(cusolver_handle, CUBLAS_OP_N, n, nrhs,
                                         d_A_sp, n, d_ipiv, d_B_sp, n, d_info));

        // Add correction: X += dX (on GPU)
        add_correction_kernel<<<numBlocks_B, blockSize>>>(d_X, d_B_sp, n * nrhs);
        CUDA_CHECK(cudaGetLastError());
    }

    // Restore default math mode
    cublasSetMathMode(cublas_handle, CUBLAS_DEFAULT_MATH);

    // Copy solution back to host
    CUDA_CHECK(cudaMemcpy(h_B_dp, d_X, n * nrhs * sizeof(cuDoubleComplex),
                          cudaMemcpyDeviceToHost));

    *info_ptr = 0;
    return 0;
}

// Double precision solver (for comparison)
int gpu_solve_double_(
    double *h_A,
    double *h_B,
    int *n_ptr,
    int *nrhs_ptr,
    int *info_ptr
) {
    int n = *n_ptr;
    int nrhs = *nrhs_ptr;
    int h_info;

    if (!gpu_initialized) {
        int device = 0;
        if (gpu_solver_init_(&device) != 0) {
            *info_ptr = -1;
            return -1;
        }
    }

    if (ensure_gpu_memory(n, nrhs) != 0) {
        *info_ptr = -2;
        return -1;
    }

    cuDoubleComplex *h_A_dp = (cuDoubleComplex*)h_A;
    cuDoubleComplex *h_B_dp = (cuDoubleComplex*)h_B;

    CUDA_CHECK(cudaMemcpy(d_A_dp, h_A_dp, n * n * sizeof(cuDoubleComplex),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_B_dp, h_B_dp, n * nrhs * sizeof(cuDoubleComplex),
                          cudaMemcpyHostToDevice));

    // Get workspace for ZGETRF
    int lwork;
    CUSOLVER_CHECK(cusolverDnZgetrf_bufferSize(cusolver_handle, n, n, d_A_dp, n, &lwork));

    if (lwork * sizeof(cuDoubleComplex) > work_size) {
        if (d_work) cudaFree(d_work);
        work_size = lwork * sizeof(cuDoubleComplex);
        CUDA_CHECK(cudaMalloc(&d_work, work_size));
    }

    CUSOLVER_CHECK(cusolverDnZgetrf(cusolver_handle, n, n, d_A_dp, n,
                                     (cuDoubleComplex*)d_work, d_ipiv, d_info));
    CUDA_CHECK(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        *info_ptr = h_info;
        return -1;
    }

    CUSOLVER_CHECK(cusolverDnZgetrs(cusolver_handle, CUBLAS_OP_N, n, nrhs,
                                     d_A_dp, n, d_ipiv, d_B_dp, n, d_info));

    CUDA_CHECK(cudaMemcpy(h_B_dp, d_B_dp, n * nrhs * sizeof(cuDoubleComplex),
                          cudaMemcpyDeviceToHost));

    *info_ptr = 0;
    return 0;
}

int gpu_is_available_() {
    int device_count = 0;
    cudaError_t err = cudaGetDeviceCount(&device_count);
    if (err != cudaSuccess || device_count == 0) {
        return 0;
    }
    return 1;
}

// Get number of available GPUs
int gpu_get_device_count_() {
    int device_count = 0;
    cudaError_t err = cudaGetDeviceCount(&device_count);
    if (err != cudaSuccess) {
        return 0;
    }
    return device_count;
}

void gpu_get_info_(char *info_str, int *len) {
    if (!gpu_is_available_()) {
        snprintf(info_str, *len, "No GPU available");
        return;
    }

    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    snprintf(info_str, *len, "%s (%.1f GB)", prop.name, prop.totalGlobalMem / 1e9);
}

// ============================================
// Multi-GPU Support using cusolverMg
// ============================================

// Initialize multi-GPU solver
int gpu_multi_init_(int *ngpu_ptr) {
    if (multi_gpu_initialized) {
        *ngpu_ptr = num_gpus;
        return 0;
    }

    // Get device count
    cudaError_t err = cudaGetDeviceCount(&num_gpus);
    if (err != cudaSuccess || num_gpus == 0) {
        fprintf(stderr, "No GPUs available for multi-GPU solver\n");
        *ngpu_ptr = 0;
        return -1;
    }

    // Limit to available GPUs
    if (*ngpu_ptr > 0 && *ngpu_ptr < num_gpus) {
        num_gpus = *ngpu_ptr;
    }

    printf("Multi-GPU solver: Using %d GPUs\n", num_gpus);

    // Allocate device ID array
    device_ids = (int*)malloc(num_gpus * sizeof(int));
    for (int i = 0; i < num_gpus; i++) {
        device_ids[i] = i;
    }

    // Print GPU info and enable peer access
    for (int i = 0; i < num_gpus; i++) {
        cudaSetDevice(device_ids[i]);
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, i);
        printf("  GPU %d: %s (%.1f GB)\n", i, prop.name, prop.totalGlobalMem / 1e9);

        // Enable peer access between all GPU pairs
        for (int j = 0; j < num_gpus; j++) {
            if (i != j) {
                int canAccess = 0;
                cudaDeviceCanAccessPeer(&canAccess, device_ids[i], device_ids[j]);
                if (canAccess) {
                    cudaDeviceEnablePeerAccess(device_ids[j], 0);
                }
            }
        }
    }
    cudaSetDevice(device_ids[0]);  // Reset to first device

    // Create cusolverMg handle
    cusolverStatus_t status = cusolverMgCreate(&cusolvermg_handle);
    if (status != CUSOLVER_STATUS_SUCCESS) {
        fprintf(stderr, "cusolverMgCreate failed: %d\n", status);
        free(device_ids);
        device_ids = NULL;
        *ngpu_ptr = 0;
        return -1;
    }

    // Set devices for the handle
    status = cusolverMgDeviceSelect(cusolvermg_handle, num_gpus, device_ids);
    if (status != CUSOLVER_STATUS_SUCCESS) {
        fprintf(stderr, "cusolverMgDeviceSelect failed: %d\n", status);
        cusolverMgDestroy(cusolvermg_handle);
        cusolvermg_handle = NULL;
        free(device_ids);
        device_ids = NULL;
        *ngpu_ptr = 0;
        return -1;
    }

    // Allocate per-GPU memory pointer arrays
    d_A_mg = (void**)malloc(num_gpus * sizeof(void*));
    d_B_mg = (void**)malloc(num_gpus * sizeof(void*));
    for (int i = 0; i < num_gpus; i++) {
        d_A_mg[i] = NULL;
        d_B_mg[i] = NULL;
    }

    multi_gpu_initialized = 1;
    *ngpu_ptr = num_gpus;
    return 0;
}

// Finalize multi-GPU solver
int gpu_multi_finalize_() {
    if (!multi_gpu_initialized) {
        return 0;
    }

    // Free per-GPU memory
    if (d_A_mg) {
        for (int i = 0; i < num_gpus; i++) {
            if (d_A_mg[i]) {
                cudaSetDevice(device_ids[i]);
                cudaFree(d_A_mg[i]);
            }
        }
        free(d_A_mg);
        d_A_mg = NULL;
    }

    if (d_B_mg) {
        for (int i = 0; i < num_gpus; i++) {
            if (d_B_mg[i]) {
                cudaSetDevice(device_ids[i]);
                cudaFree(d_B_mg[i]);
            }
        }
        free(d_B_mg);
        d_B_mg = NULL;
    }

    if (d_ipiv_mg) {
        for (int i = 0; i < num_gpus; i++) {
            if (d_ipiv_mg[i]) {
                cudaSetDevice(device_ids[i]);
                cudaFree(d_ipiv_mg[i]);
            }
        }
        free(d_ipiv_mg);
        d_ipiv_mg = NULL;
    }

    if (d_work_mg) {
        for (int i = 0; i < num_gpus; i++) {
            if (d_work_mg[i]) {
                cudaSetDevice(device_ids[i]);
                cudaFree(d_work_mg[i]);
            }
        }
        free(d_work_mg);
        d_work_mg = NULL;
    }
    mg_lwork = 0;

    // Destroy descriptors
    if (mg_desc_A) {
        cusolverMgDestroyMatrixDesc(mg_desc_A);
        mg_desc_A = NULL;
    }
    if (mg_desc_B) {
        cusolverMgDestroyMatrixDesc(mg_desc_B);
        mg_desc_B = NULL;
    }
    if (mg_grid) {
        cusolverMgDestroyGrid(mg_grid);
        mg_grid = NULL;
    }

    // Destroy handle
    if (cusolvermg_handle) {
        cusolverMgDestroy(cusolvermg_handle);
        cusolvermg_handle = NULL;
    }

    if (device_ids) {
        free(device_ids);
        device_ids = NULL;
    }

    multi_gpu_initialized = 0;
    num_gpus = 0;
    return 0;
}

// Multi-GPU solver using cusolverMg
// Uses double precision for maximum accuracy
int gpu_solve_multi_(
    double *h_A,
    double *h_B,
    int *n_ptr,
    int *nrhs_ptr,
    int *info_ptr
) {
    int n = *n_ptr;
    int nrhs = *nrhs_ptr;
    cusolverStatus_t status;

    // Initialize if not done
    if (!multi_gpu_initialized) {
        int ngpu = 0;  // Use all available
        if (gpu_multi_init_(&ngpu) != 0 || ngpu < 2) {
            // Fall back to single-GPU
            fprintf(stderr, "Multi-GPU init failed, falling back to single GPU\n");
            int max_refine = 0;
            double tol = 1e-10;
            return gpu_solve_mixed_(h_A, h_B, n_ptr, nrhs_ptr,
                                    &max_refine, &tol, info_ptr);
        }
    }

    cuDoubleComplex *h_A_z = (cuDoubleComplex*)h_A;
    cuDoubleComplex *h_B_z = (cuDoubleComplex*)h_B;

    // Block size for column distribution
    // For getrf, rows are not distributed (T_row = n), only columns are distributed
    int T_col = MG_BLOCK_SIZE;

    // Create grid (1 x num_gpus column distribution)
    if (mg_grid == NULL) {
        // Use 1 row, num_gpus columns - columns distributed across GPUs
        status = cusolverMgCreateDeviceGrid(&mg_grid, 1, num_gpus, device_ids,
                                             CUDALIBMG_GRID_MAPPING_COL_MAJOR);
        if (status != CUSOLVER_STATUS_SUCCESS) {
            fprintf(stderr, "cusolverMgCreateDeviceGrid failed: %d\n", status);
            *info_ptr = -1;
            return -1;
        }
    }

    // Create matrix descriptors
    // Both row and col block sizes should be the same (T_col) - grid determines distribution
    if (mg_desc_A) {
        cusolverMgDestroyMatrixDesc(mg_desc_A);
        mg_desc_A = NULL;
    }
    if (mg_desc_B) {
        cusolverMgDestroyMatrixDesc(mg_desc_B);
        mg_desc_B = NULL;
    }

    // A is n x n with T_col x T_col blocks
    status = cusolverMgCreateMatrixDesc(&mg_desc_A, n, n, T_col, T_col,
                                         CUDA_C_64F, mg_grid);
    if (status != CUSOLVER_STATUS_SUCCESS) {
        fprintf(stderr, "cusolverMgCreateMatrixDesc for A failed: %d\n", status);
        *info_ptr = -3;
        return -1;
    }

    // B is n x nrhs
    int T_B_col = (nrhs > T_col) ? T_col : nrhs;
    status = cusolverMgCreateMatrixDesc(&mg_desc_B, n, nrhs, T_col, T_B_col,
                                         CUDA_C_64F, mg_grid);
    if (status != CUSOLVER_STATUS_SUCCESS) {
        fprintf(stderr, "cusolverMgCreateMatrixDesc for B failed: %d\n", status);
        *info_ptr = -3;
        return -1;
    }

    // Calculate local storage size for each GPU using block-cyclic column distribution
    // Following NVIDIA's cusolverMg example formula:
    // local_n = ceil(ceil(n / T_col) / num_gpus) * T_col
    int64_t num_col_blocks_A = (n + T_col - 1) / T_col;
    int64_t num_col_blocks_B = (nrhs + T_B_col - 1) / T_B_col;

    // Each GPU gets approximately the same number of column blocks
    int64_t blocks_per_gpu_A = (num_col_blocks_A + num_gpus - 1) / num_gpus;
    int64_t blocks_per_gpu_B = (num_col_blocks_B + num_gpus - 1) / num_gpus;

    int64_t local_ncols_A_val = blocks_per_gpu_A * T_col;
    int64_t local_ncols_B_val = blocks_per_gpu_B * T_B_col;

    int64_t *local_ncols_A = (int64_t*)malloc(num_gpus * sizeof(int64_t));
    int64_t *local_ncols_B = (int64_t*)malloc(num_gpus * sizeof(int64_t));

    for (int i = 0; i < num_gpus; i++) {
        local_ncols_A[i] = local_ncols_A_val;
        local_ncols_B[i] = local_ncols_B_val;
    }

    // Allocate memory on each GPU
    if (d_A_mg == NULL) {
        d_A_mg = (void**)malloc(num_gpus * sizeof(void*));
        for (int i = 0; i < num_gpus; i++) d_A_mg[i] = NULL;
    }
    if (d_B_mg == NULL) {
        d_B_mg = (void**)malloc(num_gpus * sizeof(void*));
        for (int i = 0; i < num_gpus; i++) d_B_mg[i] = NULL;
    }
    if (d_ipiv_mg == NULL) {
        d_ipiv_mg = (int**)malloc(num_gpus * sizeof(int*));
        for (int i = 0; i < num_gpus; i++) d_ipiv_mg[i] = NULL;
    }

    for (int i = 0; i < num_gpus; i++) {
        cudaSetDevice(device_ids[i]);

        // Free old allocations
        if (d_A_mg[i]) cudaFree(d_A_mg[i]);
        if (d_B_mg[i]) cudaFree(d_B_mg[i]);
        if (d_ipiv_mg[i]) cudaFree(d_ipiv_mg[i]);

        // Allocate new - local leading dimension is n (full rows)
        size_t size_A = (size_t)n * local_ncols_A[i] * sizeof(cuDoubleComplex);
        size_t size_B = (size_t)n * local_ncols_B[i] * sizeof(cuDoubleComplex);

        if (size_A > 0) {
            cudaError_t err = cudaMalloc(&d_A_mg[i], size_A);
            if (err != cudaSuccess) {
                fprintf(stderr, "GPU %d: Failed to allocate A (%.2f GB)\n", i, size_A / 1e9);
                *info_ptr = -2;
                free(local_ncols_A); free(local_ncols_B);
                return -1;
            }
        }
        if (size_B > 0) {
            cudaError_t err = cudaMalloc(&d_B_mg[i], size_B);
            if (err != cudaSuccess) {
                fprintf(stderr, "GPU %d: Failed to allocate B\n", i);
                *info_ptr = -2;
                free(local_ncols_A); free(local_ncols_B);
                return -1;
            }
        }
        // Each GPU needs its own pivot array
        // Size is n (since row block = n, no row distribution)
        cudaMalloc(&d_ipiv_mg[i], n * sizeof(int));
    }

    // Copy data to GPUs using block-cyclic column distribution
    // GPU i gets column blocks: i, i+P, i+2P, ... (block-cyclic)
    for (int i = 0; i < num_gpus; i++) {
        cudaSetDevice(device_ids[i]);
        int64_t local_col = 0;

        // Copy A blocks
        for (int64_t b = i; b < num_col_blocks_A; b += num_gpus) {
            int64_t global_col_start = b * T_col;
            int64_t block_width = T_col;
            if (global_col_start + block_width > n) {
                block_width = n - global_col_start;
            }
            // Copy column block (full n rows for each column)
            for (int64_t c = 0; c < block_width; c++) {
                cudaMemcpy((cuDoubleComplex*)d_A_mg[i] + local_col * n,
                          h_A_z + (global_col_start + c) * n,
                          n * sizeof(cuDoubleComplex),
                          cudaMemcpyHostToDevice);
                local_col++;
            }
        }

        // Copy B blocks
        local_col = 0;
        for (int64_t b = i; b < num_col_blocks_B; b += num_gpus) {
            int64_t global_col_start = b * T_B_col;
            int64_t block_width = T_B_col;
            if (global_col_start + block_width > nrhs) {
                block_width = nrhs - global_col_start;
            }
            for (int64_t c = 0; c < block_width; c++) {
                cudaMemcpy((cuDoubleComplex*)d_B_mg[i] + local_col * n,
                          h_B_z + (global_col_start + c) * n,
                          n * sizeof(cuDoubleComplex),
                          cudaMemcpyHostToDevice);
                local_col++;
            }
        }
    }

    // Synchronize all GPUs after data transfer
    for (int i = 0; i < num_gpus; i++) {
        cudaSetDevice(device_ids[i]);
        cudaDeviceSynchronize();
    }
    cudaSetDevice(device_ids[0]);

    printf("Multi-GPU: n=%d, T_col=%d, blocks_per_gpu=%ld, local_cols=%ld\n",
           n, T_col, blocks_per_gpu_A, local_ncols_A_val);

    // Get workspace size for LU factorization
    int64_t lwork_getrf = 0;
    status = cusolverMgGetrf_bufferSize(cusolvermg_handle, n, n,
                                         d_A_mg, 1, 1, mg_desc_A,
                                         d_ipiv_mg, CUDA_C_64F, &lwork_getrf);
    if (status != CUSOLVER_STATUS_SUCCESS) {
        fprintf(stderr, "cusolverMgGetrf_bufferSize failed: %d\n", status);
        *info_ptr = -4;
        free(local_ncols_A); free(local_ncols_B);
        return -1;
    }
    printf("Multi-GPU: workspace size = %ld bytes per GPU\n", lwork_getrf);

    // Allocate workspace on each GPU if needed
    if (d_work_mg == NULL || lwork_getrf > mg_lwork) {
        if (d_work_mg) {
            for (int i = 0; i < num_gpus; i++) {
                if (d_work_mg[i]) {
                    cudaSetDevice(device_ids[i]);
                    cudaFree(d_work_mg[i]);
                }
            }
            free(d_work_mg);
        }
        d_work_mg = (void**)malloc(num_gpus * sizeof(void*));
        for (int i = 0; i < num_gpus; i++) {
            cudaSetDevice(device_ids[i]);
            cudaMalloc(&d_work_mg[i], lwork_getrf);
        }
        mg_lwork = lwork_getrf;
    }

    // LU factorization
    int h_info_mg = 0;
    status = cusolverMgGetrf(cusolvermg_handle, n, n,
                              d_A_mg, 1, 1, mg_desc_A,
                              d_ipiv_mg, CUDA_C_64F,
                              d_work_mg, lwork_getrf, &h_info_mg);
    if (status != CUSOLVER_STATUS_SUCCESS || h_info_mg != 0) {
        fprintf(stderr, "cusolverMgGetrf failed: status=%d, info=%d\n", status, h_info_mg);
        *info_ptr = (h_info_mg != 0) ? h_info_mg : -5;
        free(local_ncols_A); free(local_ncols_B);
        return -1;
    }

    // Get workspace for solve
    int64_t lwork_getrs = 0;
    status = cusolverMgGetrs_bufferSize(cusolvermg_handle, CUBLAS_OP_N, n, nrhs,
                                         d_A_mg, 1, 1, mg_desc_A,
                                         d_ipiv_mg, d_B_mg, 1, 1, mg_desc_B,
                                         CUDA_C_64F, &lwork_getrs);
    if (status != CUSOLVER_STATUS_SUCCESS) {
        fprintf(stderr, "cusolverMgGetrs_bufferSize failed: %d\n", status);
        *info_ptr = -6;
        free(local_ncols_A); free(local_ncols_B);
        return -1;
    }

    // Reallocate workspace if needed
    if (lwork_getrs > mg_lwork) {
        for (int i = 0; i < num_gpus; i++) {
            cudaSetDevice(device_ids[i]);
            cudaFree(d_work_mg[i]);
            cudaMalloc(&d_work_mg[i], lwork_getrs);
        }
        mg_lwork = lwork_getrs;
    }

    // Solve
    status = cusolverMgGetrs(cusolvermg_handle, CUBLAS_OP_N, n, nrhs,
                              d_A_mg, 1, 1, mg_desc_A,
                              d_ipiv_mg, d_B_mg, 1, 1, mg_desc_B,
                              CUDA_C_64F, d_work_mg, lwork_getrs, &h_info_mg);
    if (status != CUSOLVER_STATUS_SUCCESS || h_info_mg != 0) {
        fprintf(stderr, "cusolverMgGetrs failed: status=%d, info=%d\n", status, h_info_mg);
        *info_ptr = (h_info_mg != 0) ? h_info_mg : -7;
        free(local_ncols_A); free(local_ncols_B);
        return -1;
    }

    // Copy solution back from GPUs (reverse block-cyclic distribution)
    for (int i = 0; i < num_gpus; i++) {
        cudaSetDevice(device_ids[i]);
        int64_t local_col = 0;

        for (int64_t b = i; b < num_col_blocks_B; b += num_gpus) {
            int64_t global_col_start = b * T_B_col;
            int64_t block_width = T_B_col;
            if (global_col_start + block_width > nrhs) {
                block_width = nrhs - global_col_start;
            }
            for (int64_t c = 0; c < block_width; c++) {
                cudaMemcpy(h_B_z + (global_col_start + c) * n,
                          (cuDoubleComplex*)d_B_mg[i] + local_col * n,
                          n * sizeof(cuDoubleComplex),
                          cudaMemcpyDeviceToHost);
                local_col++;
            }
        }
    }

    free(local_ncols_A);
    free(local_ncols_B);

    *info_ptr = 0;
    return 0;
}

// Unified solver that automatically chooses single or multi-GPU
int gpu_solve_auto_(
    double *h_A,
    double *h_B,
    int *n_ptr,
    int *nrhs_ptr,
    int *max_refine_ptr,
    double *tol_ptr,
    int *info_ptr
) {
    int n = *n_ptr;
    int nrhs = *nrhs_ptr;

    // Check GPU memory availability
    // Complex double matrix: n*n*16 bytes for A, n*nrhs*16 bytes for B
    // Need ~3x for workspace (A in single precision + workspace + solution)
    size_t required_mem = (size_t)n * n * 16 * 3 + (size_t)n * nrhs * 16 * 2;

    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    size_t available_mem = prop.totalGlobalMem;

    if (required_mem > available_mem * 0.9) {
        fprintf(stderr, "Warning: Matrix too large for single GPU (%.2f GB required, %.2f GB available)\n",
                required_mem / 1e9, available_mem / 1e9);
        fprintf(stderr, "Consider using a smaller matrix or a GPU with more memory.\n");
        fprintf(stderr, "Maximum recommended matrix size for this GPU: %d x %d\n",
                (int)sqrt(available_mem * 0.9 / 48.0), (int)sqrt(available_mem * 0.9 / 48.0));
    }

    // Multi-GPU support is currently disabled due to cusolverMg API complexity
    // Single GPU can handle up to ~30k x 30k matrices on RTX 3090 (24GB)
    // TODO: Re-enable multi-GPU once cusolverMg is properly configured
#if 0 && CUDART_VERSION >= 12000
    int device_count = gpu_get_device_count_();
    if (device_count >= 2 && n > 20000) {
        printf("Using multi-GPU solver (%d GPUs) for n=%d\n", device_count, n);
        int ret = gpu_solve_multi_(h_A, h_B, n_ptr, nrhs_ptr, info_ptr);
        if (ret == 0) {
            return 0;  // Success
        }
        printf("Multi-GPU failed, falling back to single-GPU\n");
    }
#endif

    // Use single-GPU mixed precision solver
    return gpu_solve_mixed_(h_A, h_B, n_ptr, nrhs_ptr, max_refine_ptr, tol_ptr, info_ptr);
}

} // extern "C"
