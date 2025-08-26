/*!
 *=================================================================================================
 *
 *  @fn		main.c
 *  
 *  @brief	Run DFN at 10 ms steps with GITT 1A pulses; stop at 3.0 V.
 *
 *   Build with Makefile below. Assumes code generated into batt_model.c
 *
 *=================================================================================================
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "batt_model.h"  


//------------------------------------------------------------------------------------------------
//  Global variables
//------------------------------------------------------------------------------------------------

// Finite-difference step for numerical Jacobians 
static const casadi_real EPSJ = 1e-8;



/*!
 *-------------------------------------------------------------------------------------------------
 *
 *  @fn		casadi_real gitt_current(double t) 
 *
 *  @brief	GITT current profile
 *
 *-------------------------------------------------------------------------------------------------
 */
static 
casadi_real gitt_current(casadi_real t) 
{
    // Simple GITT: 10 min pulse @ +1 A, then 50 min rest @ 0 A; repeat.
    const casadi_real period = 3600.0;     // 60 min
    // const casadi_real tpulse = 1200.0;      // 10 min
    const casadi_real tpulse = 600.0;      // 10 min
    casadi_real m = fmod(t, period);
    // return (m < tpulse) ? 1.0 : 0.0;
    return (m < tpulse) ? 10.0 : 0.0;
}


#if 0
/*!
 *-------------------------------------------------------------------------------------------------
 *
 *  @fn		void dump_signature(void) 
 *
 *  @brief	Print out name of arg[] and res[] content names
 *
 *  @note	Usually just print out 'i0', 'i1',..., and 'o0','o1',...
 *
 *-------------------------------------------------------------------------------------------------
 */
static
void dump_signature(void) 
{
  printf("rhs: n_in=%lld n_out=%lld\n", rhs__n_in(), rhs__n_out());
  for (int i=0;i<rhs__n_in();++i) printf("  in[%d] = %s\n", i, rhs__name_in(i));
  for (int i=0;i<rhs__n_out();++i) printf("  out[%d] = %s\n", i, rhs__name_out(i));

  printf("alg: n_in=%lld n_out=%lld\n", alg__n_in(), alg__n_out());
  for (int i=0;i<alg__n_in();++i) printf("  in[%d] = %s\n", i, alg__name_in(i));
  for (int i=0;i<alg__n_out();++i) printf("  out[%d] = %s\n", i, alg__name_out(i));

  printf("x0: n_in=%lld n_out=%lld\n", x0_n_in(), x0_n_out());
  for (int i=0;i<x0_n_in();++i) printf("  in[%d] = %s\n", i, x0_name_in(i));
  for (int i=0;i<x0_n_out();++i) printf("  out[%d] = %s\n", i, x0_name_out(i));

  printf("z0: n_in=%lld n_out=%lld\n", z0_n_in(), z0_n_out());
  for (int i=0;i<z0_n_in();++i) printf("  in[%d] = %s\n", i, z0_name_in(i));
  for (int i=0;i<z0_n_out();++i) printf("  out[%d] = %s\n", i, z0_name_out(i));

  printf("variables: n_in=%lld n_out=%lld\n", variables_n_in(), variables_n_out());
  for (int i=0;i<variables_n_in();++i) printf("  in[%d] = %s\n", i, variables_name_in(i));
  for (int i=0;i<variables_n_out();++i) printf("  out[%d] = %s\n", i, variables_name_out(i));
}
#endif

/*!
 *-------------------------------------------------------------------------------------------------
 *
 *  @fn		casadi_int len_from_sp(const casadi_int* sp) 
 *
 *  @brief	Helper function get dense vector length from sparsity (rows*cols)
 *
 *-------------------------------------------------------------------------------------------------
 */
static 
casadi_int len_from_sp(const casadi_int* sp) 
{
    return sp[0] * sp[1]; // for dense vectors/matrices
}


/*!
 *-------------------------------------------------------------------------------------------------
 *  
 *  @fn		void model_setup(casadi_real *x, casadi_real *z, casadi_real *p, 
 *                               casadi_real *dx, casadi_real *g, casadi_real *vars,
 *                               casadi_int *nx, casadi_int *nz, casadi_int *np, casadi_int *nvar) 
 *
 *  @brief	Setup model to run
 *
 *  @param	t	time casadi_real pointer (allocated)
 *  @param	x	x vector pointer (unallocated)
 *  @param	z	z vector pointer (unallocated)	
 *  @param	p	p vector pointer (unallocated)	
 *  @param	dx	dx vector pointer (unallocated)
 *  @param	g	parameter vector pointer (unallocated)	
 *  @param	vars	variables vector pointer (unallocated)	
 *  @param      nx	pointer nx casadi_int (allocated)
 *  @param      nz	pointer nz casadi_int (allocated)
 *  @param      np	pointer np casadi_int (allocated)
 *  @param      nvar	pointer nvar casadi_int (allocated)
 *
 *-------------------------------------------------------------------------------------------------
 */
static
void model_setup(casadi_real *t, casadi_real **x, casadi_real **z, casadi_real **p, 
		 casadi_real **dx, casadi_real **g, casadi_real **vars,
		 casadi_int *nx, casadi_int *nz, casadi_int *np, casadi_int *nvar)
{
    // --- Query sizes from one representative function ---
    *t = 0.0;
    *nx  = len_from_sp(rhs__sparsity_out(0));            // len(dx/dt) == #diff. states
    *nz  = len_from_sp(alg__sparsity_out(0));            // #alg residuals == #alg states
    *np  = len_from_sp(x0_sparsity_in(0));               // x0 takes only p
    *nvar= len_from_sp(variables_sparsity_out(0));       // length of variables() stack


    // --- Allocate model buffers ---
    *x   = (casadi_real*)calloc(*nx,   sizeof(casadi_real));         // rhs
    *z   = (casadi_real*)calloc(*nz,   sizeof(casadi_real));         // residual state
    *p   = (casadi_real*)calloc(*np,   sizeof(casadi_real));         // parameters
    *dx  = (casadi_real*)calloc(*nx,   sizeof(casadi_real));         // dx rhs
    *g   = (casadi_real*)calloc(*nz,   sizeof(casadi_real));         // algebraic residual
    *vars= (casadi_real*)calloc(*nvar, sizeof(casadi_real));         // variables() output


    // --- Initial conditions for x0() ---
    {
        casadi_int sz_arg=0, sz_res=0, sz_iw=0, sz_w=0;
        x0_work(&sz_arg, &sz_res, &sz_iw, &sz_w);           	  
        casadi_int *iw = (casadi_int*)calloc(sz_iw, sizeof(casadi_int));
        casadi_real *w  = (casadi_real*)calloc(sz_w, sizeof(casadi_real));

	printf("x0_arg_sz=%lld, x0_res_sz=%lld\n", sz_arg, sz_res);
	printf("x0_n_in()=%lld, x0_n_out()=%lld\n", x0_n_in(), x0_n_out());

        const casadi_real* x0_arg[x0_n_in()]; 
        casadi_real*       x0_res[x0_n_out()];
	x0_arg[0] = *p;
	x0_res[0] = *x;

        x0(x0_arg, x0_res, iw, w, 0);

        free(iw); free(w);
    }


    // --- Initial conditions for z0() ---
    {
        casadi_int sz_arg=0, sz_res=0, sz_iw=0, sz_w=0;
        z0_work(&sz_arg, &sz_res, &sz_iw, &sz_w);           	  
        casadi_int *iw = (casadi_int*)calloc(sz_iw, sizeof(casadi_int));
        casadi_real *w  = (casadi_real*)calloc(sz_w, sizeof(casadi_real));

	printf("z0_arg_sz=%lld, z0_res_sz=%lld\n", sz_arg, sz_res);
	printf("z0_n_in()=%lld, z0_n_out()=%lld\n", z0_n_in(), z0_n_out());
	printf("rhs_n_in()=%lld, rhs_n_out()=%lld\n", rhs__n_in(), rhs__n_out());
	printf("alg_n_in()=%lld, alg_n_out()=%lld\n", alg__n_in(), alg__n_out());
	printf("var_n_in()=%lld, var_n_out()=%lld\n", variables_n_in(), variables_n_out());
	printf("sizeof(casadi_real *) = %ld, sizeof(casadi_real) = %ld\n", 
			sizeof(casadi_real *), sizeof(casadi_real));

        const casadi_real* z0_arg[z0_n_in()];
        casadi_real*       z0_res[z0_n_out()];
        z0_arg[0] = *p;
        z0_res[0] = *z;

        z0(z0_arg, z0_res, iw, w, 0);

        free(iw); free(w);
    }
}



/*!
 *-------------------------------------------------------------------------------------------------

 *  @fn		int model_eval_rhs_alg(casadi_real *t, 
 *                                     casadi_real *x, 
 *                                     casadi_real *z, 
 *                                     casadi_real *p, 
 *                                     casadi_real *dx, 
 *                                     casadi_real *g)
 *
 *  @brief	Evaluate the rhs and alg once
 *
 *
 *  @param	t	time casadi_real pointer (allocated)
 *  @param	x	x vector pointer (allocated)
 *  @param	z	z vector pointer (allocated)	
 *  @param	p	p vector pointer (allocated)	
 *
 *  @param	dx	dx vector pointer (allocated)
 *  @param	g	parameter vector pointer (allocated)	
 *
 *
 *  @return	0 if successful; negative otherwise
 *
 *-------------------------------------------------------------------------------------------------
 */
int model_eval_rhs_alg(casadi_real *t, casadi_real *x, casadi_real *z, casadi_real *p, 
                       casadi_real *dx, casadi_real *g)
{
    if (t==NULL || x==NULL || z==NULL || p==NULL || dx==NULL || g==NULL)
       return -1;

    // --- One evaluation of rhs_ ---
    {
        casadi_int sz_arg=0, sz_res=0, sz_iw=0, sz_w=0;
        rhs__work(&sz_arg, &sz_res, &sz_iw, &sz_w);           	  
        casadi_int *iw = (casadi_int*)calloc(sz_iw, sizeof(casadi_int));
        casadi_real *w  = (casadi_real*)calloc(sz_w, sizeof(casadi_real));

        const casadi_real* rhs_arg[rhs__n_in()];   // [4]
        casadi_real*       rhs_res[rhs__n_out()];  // [1]
        rhs_arg[0] = t;
        rhs_arg[1] = x;
        rhs_arg[2] = z;
        rhs_arg[3] = p;
        rhs_res[0] = dx;

        rhs_(rhs_arg, rhs_res, iw, w, 0);

        free(iw); free(w);
    }

    // --- One evaluation of alg_ ---
    {
        casadi_int sz_arg=0, sz_res=0, sz_iw=0, sz_w=0;
        alg__work(&sz_arg, &sz_res, &sz_iw, &sz_w);           	  
        casadi_int *iw = (casadi_int*)calloc(sz_iw, sizeof(casadi_int));
        casadi_real *w  = (casadi_real*)calloc(sz_w, sizeof(casadi_real));

        const casadi_real* alg_arg[alg__n_in()];
        casadi_real*       alg_res[alg__n_out()];
        alg_arg[0] = t;
        alg_arg[1] = x;
        alg_arg[2] = z;
        alg_arg[3] = p;
        alg_res[0] = g;

        alg_(alg_arg, alg_res, iw, w, 0);

        free(iw); free(w);
    }

    return 0;
}

/*!
 *-------------------------------------------------------------------------------------------------
 *                             
 *  @fn		int model_eval_vars(casadi_real *t, casadi_real *x, casadi_real *z, casadi_real *p
 *                                  casadi_real *vars) 
 *
 *  @brief	Evaluate model variables once
 *
 *  @param	t	time pointer     (allocated)
 *  @param	x	x vector pointer (allocated)
 *  @param	z	z vector pointer (allocated)
 *  @param	p 	p vector pointer (allocated)
 *
 *  @param	vars	variables vector pointer (allocated)
 * 
 *
 *  @return	0 if successful; negative otherwise
 *
 *-------------------------------------------------------------------------------------------------
 */ 
int model_eval_vars(casadi_real *t, casadi_real *x, casadi_real *z, casadi_real *p,
                     casadi_real *vars) 
{
   if (t==NULL || x==NULL || z==NULL || p==NULL || vars==NULL)
      return -1;

   casadi_int sz_arg=0, sz_res=0, sz_iw=0, sz_w=0;
   variables_work(&sz_arg, &sz_res, &sz_iw, &sz_w);           	  
   casadi_int *iw = (casadi_int*)calloc(sz_iw, sizeof(casadi_int));
   casadi_real *w  = (casadi_real*)calloc(sz_w, sizeof(casadi_real));

   const casadi_real* vars_arg[variables_n_in()];  // [4]
   casadi_real*       vars_res[variables_n_out()]; // [1]
   vars_arg[0] = t;
   vars_arg[1] = x;
   vars_arg[2] = z;
   vars_arg[3] = p;
   vars_res[0] = vars;

   variables(vars_arg, vars_res, iw, w, 0);

   free(iw); free(w);

   return 0;
}


/*!
 *-------------------------------------------------------------------------------------------------
 *
 *  @fn		void model_takedown(casadi_real *x, casadi_real *z, casadi_real *p, 
 *      	                    casadi_real *dx, casadi_real *g, casadi_real *vars)
 *
 *  @brief	Take down the model by freeing all model memories
 *
 *-------------------------------------------------------------------------------------------------
 */
void model_takedown(casadi_real *x, casadi_real *z, casadi_real *p, 
		    casadi_real *dx, casadi_real *g, casadi_real *vars)
{
   if (vars) free(vars); 
   if (g) free(g); 
   if (dx) free(dx); 
   if (p) free(p); 
   if (z) free(z); 
   if (x) free(x);
}



/*!
 *-------------------------------------------------------------------------------------------------
 *
 *  @fn		solve_linear 
 *
 *  @brief	Linear solve (Gaussian elimination, no pivoting – ok for small systems)
 *
 *  Solve A*x = b
 *
 *  Input: 	A, b
 *  Output	x (via b)
 *
 *-------------------------------------------------------------------------------------------------
 */
static 
int solve_linear(int n, casadi_real *A, casadi_real *b) 
{
    for (int k = 0; k < n; ++k) 
    {
        casadi_real piv = A[k*n + k];

        printf("---piv=%lf  n=%d   k=%d\n", piv, n, k);
        if (fabs(piv) < 1e-14) 
	{
           printf("piv=%lf  n=%d   k=%d\n", piv, n, k);
           return -1; // singular
	}

        casadi_real inv = 1.0 / piv;

        for (int j = k; j < n; ++j) 
           A[k*n + j] *= inv;

        b[k] *= inv;

        for (int i = 0; i < n; ++i) 
	{
	   if (i != k) 
	   {
               casadi_real f = A[i*n + k];
               for (int j = k; j < n; ++j) 
		  A[i*n + j] -= f * A[k*n + j];
               b[i] -= f * b[k];
           }
     	}
    }
    return 0;
}



/*!
 *-------------------------------------------------------------------------------------------------
 *
 *  @fn		int step_implicit_euler(casadi_real *t, 
 *                                      casadi_real *x, 
 *                                      casadi_real *z, 
 *                                      casadi_real *p, 
 *                                      casadi_int nx, 
 *                                      casadi_int nz, 
 *                                      casadi_real dt)
 *
 *  @brief	Newton solve for implicit Euler step:
 *                      F([x+,z+]) = [ x+ - x - dt * f(x+,z+) ; g(x+,z+) ] = 0
 *
 *
 *  @param	t	pointer to time (allocated)
 *  @param	x	pointer to x vector (allocated)
 *  @param	z	pointer to z vector (allocated)
 *  @param	p	pointer to p vector (allocated) I and Tamb
 *  @param	nx	x vector size 
 *  @param	nz	z vector size 
 *  @param	dt	time increment 
 *
 *
 *  @return	0 if success; negative otherwise
 *
 *-------------------------------------------------------------------------------------------------
 */
static 
int step_implicit_euler(casadi_real *t, casadi_real *x, casadi_real *z, casadi_real *p, 
		        casadi_int nx, casadi_int nz, casadi_real dt)
{
    casadi_real tnext = *t + dt;
    const casadi_int n = nx + nz;
    casadi_real *F  = (casadi_real*)calloc(n, sizeof(casadi_real));
    casadi_real *J  = (casadi_real*)calloc(n*n, sizeof(casadi_real));
    casadi_real *dY = (casadi_real*)calloc(n, sizeof(casadi_real));
    casadi_real *xk = (casadi_real*)malloc(nx*sizeof(casadi_real));
    casadi_real *zk = (casadi_real*)malloc(nz*sizeof(casadi_real));
    casadi_real *fx = (casadi_real*)malloc(nx*sizeof(casadi_real));
    casadi_real *gx = (casadi_real*)malloc(nz*sizeof(casadi_real));
    casadi_real *fx2= (casadi_real*)malloc(nx*sizeof(casadi_real));
    casadi_real *gx2= (casadi_real*)malloc(nz*sizeof(casadi_real));

    if (F==NULL || J==NULL || dY==NULL || xk==NULL || zk==NULL || 
             fx==NULL || gx==NULL || fx2==NULL || gx2==NULL)
    {
       perror("step_implicit_euler malloc() failed.");
       goto err_ret;
    }


    if (t==NULL || x==NULL || z==NULL || p==NULL)
    {
       perror("step_implicit_euler malloc() failed.");
       goto err_ret;
    }

    memcpy(xk, x, nx*sizeof(casadi_real));
    memcpy(zk, z, nz*sizeof(casadi_real));


    int iter, maxit=20;

    for (iter=0; iter<maxit; ++iter) 
    {
        // F(y) with y=[xk, zk]
        model_eval_rhs_alg(&tnext, xk, zk, p, fx, gx);
        for (int i=0;i<nx;++i) 
           F[i] = xk[i] - x[i] - dt*fx[i];

        for (int j=0;j<nz;++j) 
	   F[nx+j]  = gx[j];
    

        // Check convergence
        casadi_real nrm=0;
        for (int i=0;i<n;++i) 
           nrm = fmax(nrm, fabs(F[i]));

        if (nrm < 1e-8) 
	{ 
	   break; 
	}

        // Numerical Jacobian (dense): columns for [x; z]
        // Perturb x
        for (int c=0;c<nx;++c) 
	{
            casadi_real tmp = xk[c];
            xk[c] = tmp + EPSJ;
            model_eval_rhs_alg(&tnext, xk, zk, p, fx2, gx2);
            xk[c] = tmp;

            for (int r=0;r<nx;++r) 
	       J[r*n + c] = ((xk[r] - x[r] - dt*fx2[r]) - F[r]) / EPSJ;

            for (int r=0;r<nz;++r) 
	       J[(nx+r)*n + c] = (gx2[r] - gx[r]) / EPSJ;
        }

        // Perturb z
        for (int c=0;c<nz;++c) 
	{
            casadi_real tmp = zk[c];
            zk[c] = tmp + EPSJ;
            model_eval_rhs_alg(&tnext, xk, zk, p, fx2, gx2);
            zk[c] = tmp;

            for (int r=0;r<nx;++r) 
	       J[r*n + (nx+c)] = ((xk[r] - x[r] - dt*fx2[r]) - F[r]) / EPSJ;

            for (int r=0;r<nz;++r) 
	       J[(nx+r)*n + (nx+c)] = (gx2[r] - gx[r]) / EPSJ;
        }


        // Solve J * dY = -F
        for (int i=0;i<n;++i) dY[i] = -F[i];
        if (solve_linear(n, J, dY)) 
	{ 
	   iter = -1; 
	   break; 
	}

        // Update
        for (int i=0;i<nx;++i) 
	   xk[i] += dY[i];

        for (int j=0;j<nz;++j) 
	   zk[j] += dY[nx+j];
    }

    int ok = (iter>=0 && iter<maxit);

    if (ok) 
    { 
       memcpy(x, xk, nx*sizeof(casadi_real)); 
       memcpy(z, zk, nz*sizeof(casadi_real)); 
    }
    else
       printf("Not ok.....iter=%d, maxit=%d\n", iter, maxit);


    free(gx2);
    free(fx2);
    free(gx);
    free(fx);
    free(zk);
    free(xk);
    free(dY);
    free(J); 
    free(F);  

    return 0;

err_ret:
    return -1;
}



/*!
 *-------------------------------------------------------------------------------------------------
 *  Main program
 *-------------------------------------------------------------------------------------------------
 */
int main(void) 
{
    //const casadi_real dt = 0.01;       // 10 ms
    const casadi_real dt = 1.0;       // 1 s
    casadi_real t = 0.0;
    casadi_real Tamb = 298.15;         // allow changing per step
    casadi_real I = gitt_current(t);
    casadi_real *x = NULL;
    casadi_real *z = NULL;
    casadi_real *p = NULL;
    casadi_real *dx = NULL;
    casadi_real *g = NULL;
    casadi_real *vars = NULL;
    casadi_int nx = 0; 
    casadi_int nz = 0; 
    casadi_int np = 0; 
    casadi_int nvar = 0; 

    // Model setup
    model_setup(&t, &x, &z, &p, &dx, &g, &vars, &nx, &nz, &np, &nvar);
    printf("nx=%lld, nz=%lld, np=%lld, nvar=%lld\n", nx, nz, np, nvar);

    // Print CSV header
    printf("time_s,current_A,curr_meas_A,SOC,SOH,voltage_V,anode_V,cathode_V,QdisAh,QmaxAh,LAM_neg_pct,LAM_pos_pct,LLI_pct\n");

    // Coulomb counting 
    casadi_real nominal_Ah = 0.0;



    //----------------------------------------------------------------------------
    //  We’ll read Nominal cell capacity [A.h] from first variables eval below
    //  Run until voltage < 3.0 V
    //----------------------------------------------------------------------------
    for (int step = 0; step < 100000000; ++step) 
    {
	//------------------------------------------------------------------------
        //  Apply current based on GITT protocol and update parameter vector 
	//------------------------------------------------------------------------
        I = gitt_current(t);

        p[0] = I; 
        p[1] = Tamb; 

	//------------------------------------------------------------------------
        // Integration (implicit Euler)
	//------------------------------------------------------------------------
        if (step > 0) // at t=0 we print initial values first
        {    
            if (step_implicit_euler(&t, x, z, p, nx, nz, dt)) 
	    {
                fprintf(stderr, "Newton failed at t=%.6f s\n", t);
                goto err_ret;
            }
        }


	//------------------------------------------------------------------------
        // Evaluate variables
	//
        // vars[] meaning:
	//
        //    0: Voltage [V]
        //    1: Measured Current [A]
        //    2: Discharge capacity [A.h]
        //    3: Nominal cell capacity [A.h]
        //    4: X-averaged cell temperature [K]
        //    5: Anode potential [V]
        //    6: Cathode potential [V]
        //    7: Loss of lithium inventory [%]
        //    8: Loss of active material in negative electrode [%]
        //    9: Loss of active material in positive electrode [%]
	//
	//------------------------------------------------------------------------
        if (model_eval_vars(&t, x, z, p, vars) != 0)
           goto err_ret;
	
        casadi_real V      = vars[0];
        casadi_real Imeas  = vars[1];
        casadi_real QdisAh = vars[2];
        nominal_Ah         = (step==0 ? vars[3] : nominal_Ah);
        casadi_real Va     = vars[5];
        casadi_real Vc     = vars[6];
        casadi_real LLIpct = vars[7];
        casadi_real LAMneg = vars[8];
        casadi_real LAMpos = vars[9];


	//------------------------------------------------------------------------
        // Compute SOC:   SOC = 1 - Q_dis / Nominal
	//------------------------------------------------------------------------
        casadi_real SOC = (nominal_Ah > 0) ?  fmax(0.0, 1.0 - (casadi_real)(QdisAh / nominal_Ah)) : NAN;


	//------------------------------------------------------------------------
        // SOH (capacity based) ≈ 1 - max(LAMneg, LAMpos)/100  (simple, illustrative)
	//------------------------------------------------------------------------
        casadi_real SOH = 1.0 - (casadi_real)(fmax(LAMneg, LAMpos) / 100.0);

	//------------------------------------------------------------------------
        //  Output data
	//------------------------------------------------------------------------
        printf("%.3f,%.6f,%.6f,%.6f, %.6f,%.6f,%.6f,%.6f, %.6f, %.6f,%.6f,%.6f,%.6f\n",
               t, I, Imeas, SOC, SOH, V, Va, Vc, QdisAh, nominal_Ah, LAMneg, LAMpos, LLIpct);


	//------------------------------------------------------------------------
	// Check stop condition 
	//------------------------------------------------------------------------
        if (V < 3.0) break;  

	
	//------------------------------------------------------------------------
	//  Increment simulation time
	//------------------------------------------------------------------------
        t += dt;
    }


    // Free state parameters and variables
    model_takedown(x, z, p, dx, g, vars);

    return 0;

err_ret:
    return -1;
}

