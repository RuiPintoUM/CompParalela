#include "fluid_solver.h"
#include <cmath>

#define IX(i, j, k) ((i) + (M + 2) * (j) + (M + 2) * (N + 2) * (k))
#define SWAP(x0, x)                                                            \
  {                                                                            \
    float *tmp = x0;                                                           \
    x0 = x;                                                                    \
    x = tmp;                                                                   \
  }
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define LINEARSOLVERTIMES 20

#define TILE_SIZE 6

// Add sources (density or velocity)
void add_source(int M, int N, int O, float *x, float *s, float dt) {
  int size = (M + 2) * (N + 2) * (O + 2);
  for (int i = 0; i < size; i++) {
    x[i] += dt * s[i];
  }
}

// Set boundary conditions
void set_bnd(int M, int N, int O, int b, float *x) {
  int i, j;

  // Set boundary on faces
  for (i = 1; i <= M; i++) {
    for (j = 1; j <= N; j++) {
      x[IX(i, j, 0)] = b == 3 ? -x[IX(i, j, 1)] : x[IX(i, j, 1)];
      x[IX(i, j, O + 1)] = b == 3 ? -x[IX(i, j, O)] : x[IX(i, j, O)];
    }
  }
  for (i = 1; i <= N; i++) {
    for (j = 1; j <= O; j++) {
      x[IX(0, i, j)] = b == 1 ? -x[IX(1, i, j)] : x[IX(1, i, j)];
      x[IX(M + 1, i, j)] = b == 1 ? -x[IX(M, i, j)] : x[IX(M, i, j)];
    }
  }
  for (i = 1; i <= M; i++) {
    for (j = 1; j <= O; j++) {
      x[IX(i, 0, j)] = b == 2 ? -x[IX(i, 1, j)] : x[IX(i, 1, j)];
      x[IX(i, N + 1, j)] = b == 2 ? -x[IX(i, N, j)] : x[IX(i, N, j)];
    }
  }

  // Set corners
  x[IX(0, 0, 0)] = 0.33f * (x[IX(1, 0, 0)] + x[IX(0, 1, 0)] + x[IX(0, 0, 1)]);
  x[IX(M + 1, 0, 0)] =
      0.33f * (x[IX(M, 0, 0)] + x[IX(M + 1, 1, 0)] + x[IX(M + 1, 0, 1)]);
  x[IX(0, N + 1, 0)] =
      0.33f * (x[IX(1, N + 1, 0)] + x[IX(0, N, 0)] + x[IX(0, N + 1, 1)]);
  x[IX(M + 1, N + 1, 0)] = 0.33f * (x[IX(M, N + 1, 0)] + x[IX(M + 1, N, 0)] +
                                    x[IX(M + 1, N + 1, 1)]);
}


void lin_solve(int M, int N, int O, int b, float *x, float *x0, float a, float c) {
  float new_c = 1 / c;
  // Loop over linear solve steps
  for (int l = 0; l < LINEARSOLVERTIMES; l++) {
    // Iterate over tiles in the k, j, i directions
    for (int bk = 1; bk <= O; bk += TILE_SIZE) {  
      int bk_min = (bk + TILE_SIZE < O + 1) ? bk + TILE_SIZE : O + 1;
      
      for (int bj = 1; bj <= N; bj += TILE_SIZE) {
        int bj_min = (bj + TILE_SIZE < N + 1) ? bj + TILE_SIZE : N + 1;
        
        for (int bi = 1; bi <= M; bi += TILE_SIZE) {
          int bi_min = (bi + TILE_SIZE < M + 1) ? bi + TILE_SIZE : M + 1;
          
          // Iterate over elements within each tile
          for (int k = bk; k < bk_min; k++) {
            for (int j = bj; j < bj_min; j++) {
              for (int i = bi; i < bi_min; i++) {
                x[IX(i, j, k)] = (x0[IX(i, j, k)] +
                                  a * (x[IX(i - 1, j, k)] + x[IX(i + 1, j, k)] +
                                      x[IX(i, j - 1, k)] + x[IX(i, j + 1, k)] +
                                      x[IX(i, j, k - 1)] + x[IX(i, j, k + 1)])) * new_c;
              }
            }
          }
        }
      }
    }
    // Apply boundary conditions
    set_bnd(M, N, O, b, x);
  }
}


// Diffusion step (uses implicit method)
void diffuse(int M, int N, int O, int b, float *x, float *x0, float diff,
             float dt) {
  int max = (M > N ? (M > O ? M : O) : (N > O ? N : O)); 
  float a = dt * diff * max * max;
  float c = 1 + 6 * a; 
  lin_solve(M, N, O, b, x, x0, a, c);
}

// Advection step (uses velocity field to move quantities)
void advect(int M, int N, int O, int b, float *d, float *d0, float *u, float *v,
            float *w, float dt) {
  float dtX = dt * M, dtY = dt * N, dtZ = dt * O;

  for (int k = 1; k <= O; k++) {
    for (int j = 1; j <= N; j++) {
      for (int i = 1; i <= M; i++) {
        float x = i - dtX * u[IX(i, j, k)];
        float y = j - dtY * v[IX(i, j, k)];
        float z = k - dtZ * w[IX(i, j, k)];

        // Clamp to grid boundaries
        if (x < 0.5f)
          x = 0.5f;
        if (x > M + 0.5f)
          x = M + 0.5f;
        if (y < 0.5f)
          y = 0.5f;
        if (y > N + 0.5f)
          y = N + 0.5f;
        if (z < 0.5f)
          z = 0.5f;
        if (z > O + 0.5f)
          z = O + 0.5f;

        int i0 = (int)x, i1 = i0 + 1;
        int j0 = (int)y, j1 = j0 + 1;
        int k0 = (int)z, k1 = k0 + 1;

        float s1 = x - i0, s0 = 1 - s1;
        float t1 = y - j0, t0 = 1 - t1;
        float u1 = z - k0, u0 = 1 - u1;

        d[IX(i, j, k)] =
            s0 * (t0 * (u0 * d0[IX(i0, j0, k0)] + u1 * d0[IX(i0, j0, k1)]) +
                  t1 * (u0 * d0[IX(i0, j1, k0)] + u1 * d0[IX(i0, j1, k1)])) +
            s1 * (t0 * (u0 * d0[IX(i1, j0, k0)] + u1 * d0[IX(i1, j0, k1)]) +
                  t1 * (u0 * d0[IX(i1, j1, k0)] + u1 * d0[IX(i1, j1, k1)]));
      }
    }
  }
  set_bnd(M, N, O, b, d);
}


void project(int M, int N, int O, float *u, float *v, float *w, float *p, float *div) {
  float halfInvMaxDim = 0.5f / (M > N ? (M > O ? M : O) : (N > O ? N : O));

  // Tiled loop for better cache locality
  for (int k0 = 1; k0 <= O; k0 += TILE_SIZE) {
    int k0_min = (k0 + TILE_SIZE < O + 1) ? k0 + TILE_SIZE : O + 1;
    
    for (int j0 = 1; j0 <= N; j0 += TILE_SIZE) {
      int j0_min = (j0 + TILE_SIZE < N + 1) ? j0 + TILE_SIZE : N + 1;
      
      for (int i0 = 1; i0 <= M; i0 += TILE_SIZE) {
        int i0_min = (i0 + TILE_SIZE < M + 1) ? i0 + TILE_SIZE : M + 1;

        // Process each tile/block
        for (int k = k0; k < k0_min; k++) {
          for (int j = j0; j < j0_min; j++) {
            for (int i = i0; i < i0_min; i++) {
              div[IX(i, j, k)] =
                  -halfInvMaxDim * (u[IX(i + 1, j, k)] - u[IX(i - 1, j, k)] +
                                    v[IX(i, j + 1, k)] - v[IX(i, j - 1, k)] +
                                    w[IX(i, j, k + 1)] - w[IX(i, j, k - 1)]);
              p[IX(i, j, k)] = 0.0f; // Clear p array as before
            }
          }
        }
      }
    }
  }

  // Boundary handling (set_bnd)
  set_bnd(M, N, O, 0, div);
  set_bnd(M, N, O, 0, p);

  // Solve pressure equation
  lin_solve(M, N, O, 0, p, div, 1, 6);

  // Adjust velocities
  for (int k0 = 1; k0 <= O; k0 += TILE_SIZE) {
    int k0_min = (k0 + TILE_SIZE < O + 1) ? k0 + TILE_SIZE : O + 1;
    
    for (int j0 = 1; j0 <= N; j0 += TILE_SIZE) {
      int j0_min = (j0 + TILE_SIZE < N + 1) ? j0 + TILE_SIZE : N + 1;
      
      for (int i0 = 1; i0 <= M; i0 += TILE_SIZE) {
        int i0_min = (i0 + TILE_SIZE < M + 1) ? i0 + TILE_SIZE : M + 1;

        for (int k = k0; k < k0_min; k++) {
          for (int j = j0; j < j0_min; j++) {
            for (int i = i0; i < i0_min; i++) {
              u[IX(i, j, k)] -= 0.5f * (p[IX(i + 1, j, k)] - p[IX(i - 1, j, k)]);
              v[IX(i, j, k)] -= 0.5f * (p[IX(i, j + 1, k)] - p[IX(i, j - 1, k)]);
              w[IX(i, j, k)] -= 0.5f * (p[IX(i, j, k + 1)] - p[IX(i, j, k - 1)]);
            }
          }
        }
      }
    }
  }

  // Boundary handling for u, v, w (velocity fields)
  set_bnd(M, N, O, 1, u);
  set_bnd(M, N, O, 2, v);
  set_bnd(M, N, O, 3, w);
}

// Step function for density
void dens_step(int M, int N, int O, float *x, float *x0, float *u, float *v,
               float *w, float diff, float dt) {
  add_source(M, N, O, x, x0, dt);
  SWAP(x0, x);
  diffuse(M, N, O, 0, x, x0, diff, dt);
  SWAP(x0, x);
  advect(M, N, O, 0, x, x0, u, v, w, dt);
}

// Step function for velocity
void vel_step(int M, int N, int O, float *u, float *v, float *w, float *u0,
              float *v0, float *w0, float visc, float dt) {
  add_source(M, N, O, u, u0, dt);
  add_source(M, N, O, v, v0, dt);
  add_source(M, N, O, w, w0, dt);
  SWAP(u0, u);
  diffuse(M, N, O, 1, u, u0, visc, dt);
  SWAP(v0, v);
  diffuse(M, N, O, 2, v, v0, visc, dt);
  SWAP(w0, w);
  diffuse(M, N, O, 3, w, w0, visc, dt);
  project(M, N, O, u, v, w, u0, v0);
  SWAP(u0, u);
  SWAP(v0, v);
  SWAP(w0, w);
  advect(M, N, O, 1, u, u0, u0, v0, w0, dt);
  advect(M, N, O, 2, v, v0, u0, v0, w0, dt);
  advect(M, N, O, 3, w, w0, u0, v0, w0, dt);
  project(M, N, O, u, v, w, u0, v0);
}