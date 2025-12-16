/*
 * cusolver_stub.c - Stub functions for CPU-only builds (no CUDA)
 */

int gpu_solver_init_(int *device_id) { return -1; }
int gpu_solver_finalize_(void) { return 0; }
int gpu_is_available_(void) { return 0; }
int gpu_get_device_count_(void) { return 0; }
void gpu_get_info_(char *info, int *len) {
    const char *msg = "No GPU available";
    int i;
    for (i = 0; msg[i] && i < *len-1; i++) info[i] = msg[i];
    info[i] = '\0';
}
int gpu_solve_double_(void *A, void *B, int *n, int *nrhs) { return -1; }
int gpu_solve_mixed_(void *A, void *B, int *n, int *nrhs) { return -1; }
int gpu_solve_tf32_(void *A, void *B, int *n, int *nrhs) { return -1; }
int gpu_solve_auto_(void *A, void *B, int *n, int *nrhs) { return -1; }
int gpu_multi_init_(int *num_gpus) { return -1; }
int gpu_multi_finalize_(void) { return 0; }
int gpu_solve_multi_(void *A, void *B, int *n, int *nrhs) { return -1; }
