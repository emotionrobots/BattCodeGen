/*!
 *===========================================================================================================
 *
 *  @file	 integrator.c 
 *
 *  @brief	ODE integrator implementation 
 *
 *===========================================================================================================
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "integrator.h"
#include "linalg.h"


/*!
 *----------------------------------------------------------------------------------------------------------
 *
 *  @fn		void build_residual(const DaeSystem *sys, double *t, double h, const double *x_old, ModelCtx *ctx, 
		    double *R) 
 *
 *  @brief	Helper to build residual vector R = [F1;F2]
 *
 *----------------------------------------------------------------------------------------------------------
 */
static 
void build_residual(const DaeSystem *sys, double *t, double h, const double *x_old, ModelCtx *ctx, 
		    double *R) 
{
   size_t nx = sys->nx;
   size_t nz = sys->nz;

   // F1 = x - x_old - h*f(t,x,z,u)
   sys->f(t, ctx);
   for (size_t i = 0; i < nx; ++i) 
      R[i] = ctx->x[i] - x_old[i] - h*ctx->f[i];

   // F2 = g(t,x,z,u) 
   if (nz) 
   {
      sys->g(t, ctx);
      for (size_t j = 0; j < nz; ++j) 
	 R[nx + j] = ctx->g[j];
   }
}


/*!
 *----------------------------------------------------------------------------------------------------------
 *
 *  @fn		double l2_norm(const double *v, size_t n)
 *
 *  @brief	Compute L2 norm of vector 'v'
 *
 *----------------------------------------------------------------------------------------------------------
 */
static 
double l2_norm(const double *v, size_t n)
{
   double s=0.0; 

   for(size_t i=0;i<n;++i) 
      s += v[i]*v[i];

   return sqrt(s);
}




/*!
 *----------------------------------------------------------------------------------------------------------
 *
 *  @fn		int dae_step_backward_euler(const DaeSystem *sys, const DaeSettings *cfg, 
 *                                          double *t, double h, ModelCtx *ctx) 
 *
 *  @brief	Perform a DAE Backward Euler step
 *
 *----------------------------------------------------------------------------------------------------------
 */
int dae_step_backward_euler(const DaeSystem *sys, const DaeSettings *cfg, double *t, double h, ModelCtx *ctx) 
{
   int iter; int rc=0;
   const size_t nx=sys->nx, nz=sys->nz, n=nx+nz;
   double x_old[nx];
   memcpy(x_old, ctx->x, nx*sizeof(double));
   double R[n];
   double J[n*n];
   int piv[n];
   double wk[(nx + (nz? nz:1) + 2*n)];


   for(iter=0; iter<cfg->newton_max_iter; ++iter)
   {
      // Residual at current (x,z)
      build_residual(sys, t, h, x_old, ctx, R);
      double nr = l2_norm(R, n);

      if (nr < cfg->newton_tol) 
      { 
	 rc=0; 
	 break; 
      }

      // Jacobian of residual wrt [x;z]
      sys->J(t, h, ctx, J);

      // Solve J * delta = -R
      for(size_t i=0;i<n;++i) 
	 R[i] = -R[i];

      if (lu_decompose(J, (int)n, piv)!=0 || lu_solve(J,(int)n,piv,R)!=0) 
      { 
	 rc=2; 
	 break; 
      }

      // Apply update
      for(size_t i=0;i<nx;++i) 
	 ctx->x[i] += R[i];

      for(size_t j=0;j<nz;++j) 
	 ctx->z[j] += R[nx+j];

      if (cfg->verbose) 
      {
         fprintf(stderr, " Newton iter %d, ||R||=%e\n", iter+1, nr);
      }

      if (l2_norm(R, n) < cfg->newton_tol*0.1) 
      { /* small step */ 
   //      rc = 0;
//	 break;
      }
   }

   if (iter==cfg->newton_max_iter) 
      rc=1;

   return rc;
}



