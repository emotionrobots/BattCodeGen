/*!
 *============================================================================================================
 * 
 *  @file	integrator.h
 *
 *  @brief	ODE integrator header file
 *
 *============================================================================================================
 */
#ifndef __INTEGRATOR_H__
#define __INTEGRATOR_H__
#include <stddef.h>
#include "global.h"


typedef struct {
   size_t nx; // number of differential states x
   size_t nz; // number of algebraic states z (can be 0 for ODE)
   size_t np; // number of inputs u (can be 0)
   void *user; // user data passed to callbacks
   void (*f)(double *t,  ModelCtx *ctx);
   void (*g)(double *t,  ModelCtx *ctx);
   // Optional Jacobian callback to fill block matrix J of size (nx+nz)x(nx+nz):
   // J = [ I - h*fx -h*fz ]
   // [ gx gz ] evaluated at (t_{n+1}, x, z)
   // If NULL, numeric finite differences are used.
   void (*J)(double *t, double h, ModelCtx *ctx, double *J);
} DaeSystem;


// Integrator settings (tolerances & limits for Newton iterations)
typedef struct {
   double newton_tol; // L2 norm tolerance on residual [F1;F2]
   int newton_max_iter; // max Newton iterations per step
   double fd_eps; // finite-difference epsilon for numeric Jacobians
   int verbose; // print Newton diagnostics if >0
} DaeSettings;


// One implicit backward-Euler step. Returns 0 on success, non-zero on failure.
int dae_step_backward_euler(const DaeSystem *sys, const DaeSettings *cfg,
                            double *t, double h, ModelCtx *ctx);


#endif // __INTEGRATOR_H__
