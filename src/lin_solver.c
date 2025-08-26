/*!
 *====================================================================================================
 *
 *  @file 	lin_solver.c 
 *
 *  @brief	Robust sparse linear solver in C:
 *  			- CSR matrix
 *  			- Jacobi-preconditioned BiCGSTAB (stabilized)
 *
 *====================================================================================================
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>


#ifndef MIN
#define MIN(a,b) ((a)<(b)?(a):(b))
#endif
#ifndef MAX
#define MAX(a,b) ((a)>(b)?(a):(b))
#endif


/* ===================== CSR Sparse Matrix ===================== */
typedef struct {
    int n;          // dimension (n x n)
    int nnz;        // number of nonzeros
    int *rowptr;    // length n+1
    int *colind;    // length nnz
    double *val;    // length nnz
} CSR;


/* ===================== Jacobi Preconditioner ===================== */
typedef struct {
    int n;
    double *dinv; // inverse of diagonal entries (safe, sets 1 if diag=0)
} Jacobi;



/* ===================== BiCGSTAB Solver (Left-Preconditioned) ===================== */
typedef struct {
    int maxit;
    double tol;     // relative residual tolerance (||r||/||b||)
    int iters;      // output
    double relres;  // output
} BiCGSTAB_Options;



/*!
 *----------------------------------------------------------------------------------------------------
 *
 *  @fn		void csr_free(CSR *A)
 *
 *  @brief  	Free CSR
 *
 *----------------------------------------------------------------------------------------------------
 */
static 
void csr_free(CSR *A)
{
    if(!A) return;

    free(A->rowptr); 
    free(A->colind); 
    free(A->val);
    A->rowptr = NULL; 
    A->colind = NULL; 
    A->val = NULL; 
    A->n = 0; 
    A->nnz = 0;
}


/*!
 *----------------------------------------------------------------------------------------------------
 *
 *  @fn		void csr_matvec(const CSR *A, const double *x, double *y)
 *
 *  @brief	CSR matrix to vector
 *
 *----------------------------------------------------------------------------------------------------
 */
static 
void csr_matvec(const CSR *A, const double *x, double *y)
{
    int n = A->n;

    for (int i=0; i<n; ++i)
    {
        double s = 0.0;
        for(int k=A->rowptr[i]; k<A->rowptr[i+1]; ++k)
	{
            s += A->val[k] * x[A->colind[k]];
        }
        y[i] = s;
    }
}


/*!
 *----------------------------------------------------------------------------------------------------
 *
 *  @fn		void vec_fill(double *a, int n, double v)
 *
 *  @brief 	Simple vector fill value
 *
 *----------------------------------------------------------------------------------------------------
 */
static 
void vec_fill(double *a, int n, double v)
{ 
   for (int i=0; i<n; ++i) 
      a[i]=v; 
}


/*!
 *----------------------------------------------------------------------------------------------------
 *
 *  @fn		void vec_copy(double *dst, const double *src, int n)
 *
 *  @brief	Simple vector copy
 *
 *----------------------------------------------------------------------------------------------------
 */
static 
void vec_copy(double *dst, const double *src, int n)
{ 
   memcpy(dst, src, n*sizeof(double)); 
}


/*!
 *----------------------------------------------------------------------------------------------------
 *
 *  @fn		void vec_axpby(double *y, double a, const double *x, double b, int n)
 *  
 *  @brief  	Compute  y = a*x + b*y
 *
 *----------------------------------------------------------------------------------------------------
 */
static 
void vec_axpby(double *y, double a, const double *x, double b, int n)
{
    for (int i=0;i<n;++i) 
       y[i] = a*x[i] + b*y[i];
}

/*!
 *----------------------------------------------------------------------------------------------------
 *
 *  @fn		void vec_add(double *z, const double *x, const double *y, int n)
 *  
 *  @brief	Adding two vectors
 *
 *----------------------------------------------------------------------------------------------------
 */
static 
void vec_add(double *z, const double *x, const double *y, int n)
{
   for (int i=0; i<n; ++i) 
      z[i] = x[i] + y[i];
}


/*!
 *----------------------------------------------------------------------------------------------------
 *
 *  @fn		void vec_sub(double *z, const double *x, const double *y, int n)
 *
 *  @brief	Subtract two vectors
 *
 *----------------------------------------------------------------------------------------------------
 */
static 
void vec_sub(double *z, const double *x, const double *y, int n)
{
   for (int i=0; i<n; ++i) 
      z[i] = x[i] - y[i];
}


/*!
 *----------------------------------------------------------------------------------------------------
 *
 *  @fn		double vec_dot(const double *x, const double *y, int n)
 *
 *  @brief	Perform dot product of two vectors
 *
 *----------------------------------------------------------------------------------------------------
 */
static 
double vec_dot(const double *x, const double *y, int n)
{
   double s = 0.0;

   for (int i=0; i<n; ++i) 
      s += x[i]*y[i];
   return s;
}


/*!
 *----------------------------------------------------------------------------------------------------
 *
 *  @fn		double vec_nrm2(const double *x, int n)
 *  
 *  @brief	Compute L2 norm of a vector
 *
 *----------------------------------------------------------------------------------------------------
 */
static 
double vec_nrm2(const double *x, int n)
{
   return sqrt(vec_dot(x,x,n));
}



/*!
 *----------------------------------------------------------------------------------------------------
 *
 *  @fn		void jacobi_build(const CSR *A, Jacobi *M)
 *  
 *  @brief	Build Jacobian M from CSR A
 *
 *----------------------------------------------------------------------------------------------------
 */
static 
void jacobi_build(const CSR *A, Jacobi *M)
{
    M->n = A->n;
    M->dinv = (double*)malloc(sizeof(double)*A->n);
    for(int i=0;i<A->n;++i)
    {
        double di = 0.0;
        for(int k=A->rowptr[i]; k<A->rowptr[i+1]; ++k)
	{
            if (A->colind[k] == i)
	    { 
	       di = A->val[k]; 
	       break; 
	    }
        }

        if (fabs(di) < 1e-30) 
           M->dinv[i] = 1.0; 
	else 
           M->dinv[i] = 1.0/di;
    }
}


/*!
 *----------------------------------------------------------------------------------------------------
 *
 *  @fn		void jacobi_apply(const Jacobi *M, const double *r, double *z)
 *  
 *  @brief	Apply Jacobian: z = M^{-1} r
 *
 *----------------------------------------------------------------------------------------------------
 */
static 
void jacobi_apply(const Jacobi *M, const double *r, double *z)
{
    for (int i=0; i<M->n; ++i) 
       z[i] = M->dinv[i] * r[i];
}


/*!
 *----------------------------------------------------------------------------------------------------
 *
 *  @fn		void jacobi_free(Jacobi *M)
 *
 *  @brief	Free jacobian
 *
 *----------------------------------------------------------------------------------------------------
 */
static 
void jacobi_free(Jacobi *M)
{
    free(M->dinv); 
    M->dinv = NULL; 
    M->n = 0;
}



/*!
 *----------------------------------------------------------------------------------------------------
 *
 *  @fn		int safe_isfinite(double x)
 *
 *  @brief	Check if x is infinite
 *
 *----------------------------------------------------------------------------------------------------
 */
static 
int safe_isfinite(double x)
{ 
   return isfinite(x); 
}


/*!
 *----------------------------------------------------------------------------------------------------
 *
 *  @fn		int bicgstab_solve(const CSR *A, const Jacobi *M, const double *b, double *x, 
 *                                 BiCGSTAB_Options *opt)
 *
 *  @brief	Linear solve private
 *
 *----------------------------------------------------------------------------------------------------
 */
static 
int bicgstab_solve(const CSR *A, const Jacobi *M, const double *b, double *x, BiCGSTAB_Options *opt)
{
    const int n = A->n;
    const int maxit = opt->maxit;
    const double tol = opt->tol;
    const double eps = 1e-30;

    double *r  = (double*)malloc(sizeof(double)*n);
    double *rhat = (double*)malloc(sizeof(double)*n); // shadow r̂
    double *p  = (double*)malloc(sizeof(double)*n);
    double *v  = (double*)malloc(sizeof(double)*n);
    double *s  = (double*)malloc(sizeof(double)*n);
    double *t  = (double*)malloc(sizeof(double)*n);
    double *ph = (double*)malloc(sizeof(double)*n);
    double *sh = (double*)malloc(sizeof(double)*n);
    double *Ax = (double*)malloc(sizeof(double)*n);

    // r = b - A*x
    csr_matvec(A, x, Ax);
    for (int i=0;i<n;++i) r[i] = b[i] - Ax[i];
    vec_copy(rhat, r, n);
    vec_fill(p, n, 0.0);
    vec_fill(v, n, 0.0);

    const double bnorm = fmax(vec_nrm2(b,n), 1e-30);
    double rho_old = 1.0, alpha = 1.0, omega = 1.0, rho = 1.0;

    int iter;
    for(iter=1; iter<=maxit; ++iter)
    {
        // Periodic residual replacement (polishing)
        if (iter % 50 == 0)
	{
            csr_matvec(A, x, Ax);
            for (int i=0;i<n;++i) 
	       r[i] = b[i] - Ax[i];
        }

        rho = vec_dot(rhat, r, n);
        if (fabs(rho) < eps)
	{
            // pick a new r̂ and restart
            vec_copy(rhat, r, n);
            rho = vec_dot(rhat, r, n);
            if (fabs(rho) < eps) 
	       break; // give up
        }

        double beta;
        if (fabs(omega) < 1e-20)
	{
            // avoid blow-up in beta: restart direction
            vec_copy(p, r, n);
            vec_fill(v, n, 0.0);
            beta = 0.0;
        } 
	else 
	{
            beta = (rho/rho_old) * (alpha/omega);
        }

        // p = r + beta*(p - omega*v)
        for (int i=0;i<n;++i) 
	   p[i] = r[i] + beta*(p[i] - omega*v[i]);

        // ph = M^{-1} p
        jacobi_apply(M, p, ph);

        // v = A*ph
        csr_matvec(A, ph, v);

        double rhv = vec_dot(rhat, v, n);
        if (fabs(rhv) < eps)
	{
            // restart with fresh r̂
            vec_copy(rhat, r, n);
            rhv = vec_dot(rhat, v, n);
            if (fabs(rhv) < eps) 
	       break;
        }
        alpha = rho / rhv;

        // s = r - alpha*v
        for (int i=0;i<n;++i) 
           s[i] = r[i] - alpha*v[i];

        // Early convergence on s
        double snorm = vec_nrm2(s,n);
        if (snorm/bnorm < tol)
	{
            for (int i=0;i<n;++i) 
	       x[i] += alpha*ph[i];
            break;
        }

        // sh = M^{-1} s
        jacobi_apply(M, s, sh);

        // t = A*sh
        csr_matvec(A, sh, t);

        double tt = vec_dot(t, t, n);
        double ts = vec_dot(t, s, n);

        if (tt < eps)
	{
            // direction degenerate; force a restart: finalize x with alpha*ph, recompute r, continue
            for (int i=0;i<n;++i) 
	       x[i] += alpha*ph[i];

            csr_matvec(A, x, Ax);
            for (int i=0;i<n;++i) 
	       r[i] = b[i] - Ax[i];

            omega = 1.0; // neutral to keep beta finite
            rho_old = rho;
            continue;
        }

        omega = ts / tt;

        // x = x + alpha*ph + omega*sh
        for (int i=0;i<n;++i) 
	   x[i] += alpha*ph[i] + omega*sh[i];

        // r = s - omega*t
        for (int i=0;i<n;++i) 
           r[i] = s[i] - omega*t[i];

        // Convergence & health checks
        double rnorm = vec_nrm2(r,n);
        double relres = rnorm / bnorm;
        if (relres < tol) 
	   break;

        if (!safe_isfinite(relres) || !safe_isfinite(omega) || !safe_isfinite(alpha)) 
	   break;

        if (fabs(omega) < eps) 
	{
            // restart direction
            vec_copy(p, r, n);
            vec_fill(v, n, 0.0);
            omega = 1.0; // neutralize for next beta
        }

        rho_old = rho;
    }

    // Final residual (for reporting)
    csr_matvec(A, x, Ax);
    for (int i=0;i<n;++i) 
       r[i] = b[i] - Ax[i];

    opt->relres = vec_nrm2(r,n) / fmax(vec_nrm2(b,n),1e-30);
    opt->iters = iter;

    free(r); 
    free(rhat); 
    free(p); 
    free(v); 
    free(s); 
    free(t);
    free(ph); 
    free(sh); 
    free(Ax);

    return (opt->relres <= opt->tol) ? 0 : 1;
}



/*!
 *----------------------------------------------------------------------------------------------------
 *
 *  @fn		int lin_solve(const CSR *A, const double *b, double *x, const int n)
 *
 *  @brief	Public API of linear solver  Solving  A*x=b, given A (n x n) and b (n), solve x.
 *
 *  @return	0 if success; negative otherwise
 *
 *----------------------------------------------------------------------------------------------------
 */
int lin_solve(const CSR *A, const double *b, double *x, const int n)
{
   Jacobi M;
   int rc = -1;

   BiCGSTAB_Options opt;
   opt.maxit=2000;
   opt.tol=1e-8;
   opt.iters=0;
   opt.relres=0.0;

   if (A==NULL || b==NULL || x==NULL)
      goto err_ret;

   // Precondition
   jacobi_build(A, &M);

   // Solve
   if (bicgstab_solve(A, &M, b, x, &opt)==0)
      rc = 0; 
    
err_ret:
   jacobi_free(&M);
   return rc;
}


#if 0
//====================================================================================================


/*!
 *----------------------------------------------------------------------------------------------------
 *
 *  @fn		CSR build_random_dd_csr(int n, int k_off, unsigned int seed, double **diag_out)
 *
 *  @brief	Random Sparse Matrix Generator - build an n x n sparse matrix with 'k_off' 
 *              random off-diagonals per row plus a strictly diagonally dominant diagonal. 
 *              Reproducible with fixed seed.
 *
 *----------------------------------------------------------------------------------------------------
 */
static 
CSR build_random_dd_csr(int n, int k_off, unsigned int seed, double **diag_out)
{
    srand(seed);

    const int nnz = n * (k_off + 1);
    CSR A;
    A.n = n; A.nnz = nnz;
    A.rowptr = (int*)malloc(sizeof(int)*(n+1));
    A.colind = (int*)malloc(sizeof(int)*nnz);
    A.val    = (double*)malloc(sizeof(double)*nnz);

    double *diag = NULL;
    if (diag_out){ diag = (double*)malloc(sizeof(double)*n); *diag_out = diag; }

    int cursor = 0; A.rowptr[0] = 0;
    unsigned char *used = (unsigned char*)calloc(n,1);

    const double off_scale = 0.1; // tamer off-diagonals ~(-0.05, 0.05)
    for (int i=0;i<n;++i){
        memset(used, 0, n);
        used[i] = 1;
        double sumabs = 0.0;
        int added = 0;
        while (added < k_off){
            int j = rand() % n;
            if (!used[j]){
                used[j] = 1;
                double v = off_scale * ( ((double)rand()/(double)RAND_MAX) - 0.5 ); // ~(-0.05,0.05)
                if (fabs(v) < 1e-8) v = (v >= 0 ? 1e-8 : -1e-8);
                A.colind[cursor] = j;
                A.val[cursor] = v;
                cursor++;
                sumabs += fabs(v);
                added++;
            }
        }
        double d = sumabs + 1.0; // strictly diagonally dominant
        if (diag) diag[i] = d;
        A.colind[cursor] = i;
        A.val[cursor] = d;
        cursor++;
        A.rowptr[i+1] = cursor;
    }
    free(used);
    return A;
}


//====================================================================================================

/*!
 *----------------------------------------------------------------------------------------------------
 *  Main()
 *----------------------------------------------------------------------------------------------------
 */
void main()
{
   int n = 500;
   const int k_off = 100;
   double *diag = NULL;
   CSR A = build_random_dd_csr(n, k_off, 24u, &diag);

   double *x_true = (double *)malloc(sizeof(double)*n);
   for (int i=0;i<n;++i) 
      x_true[i] = sin(0.01*i) + 0.3*cos(0.007*i);

   double *b = (double *)malloc(sizeof(double)*n);
   csr_matvec(&A, x_true, b);

   double *x = (double *)malloc(sizeof(double)*n);
   vec_fill(x, n, 0.0);


   clock_t t0 = clock();
   int rc = lin_solve(&A, b, x, n);
   if (rc < 0)
      printf("Error!!!!!\n");

   clock_t t1 = clock();
   double secs = (double)(t1 - t0) / CLOCKS_PER_SEC;

   // Validate against x_true (relative error)
   double *diff = (double*)malloc(sizeof(double)*n);
   vec_sub(diff, x, x_true, n);
   double rel_err = vec_nrm2(diff, n) / fmax(vec_nrm2(x_true, n), 1e-30);

   // Report
   printf("=== BiCGSTAB (Jacobi) on %dx%d sparse system ===\n", n, n);
   printf("nnz = %d (avg %.1f per row)\n", A.nnz, (double)A.nnz / n);
   printf("relative solution error = %lf\n", rel_err);
   printf("elapsed time = %.3f s\n", secs);
   printf("UNIT TEST: %s\n", (rc==0 && rel_err < 1e-6) ? "PASS" : "FAIL");

   // Cleanup
   free(diff); free(x); free(b); free(x_true); free(diag);
   csr_free(&A); 
}

#endif

