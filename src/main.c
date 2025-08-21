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

        if (fabs(piv) < 1e-14) 
           return -1; // singular
		     
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
 *  @fn		void eval_rhs_alg(const casadi_real *x, const casadi_real *z,
 *                                casadi_real t, casadi_real I, casadi_real Tamb,
 *                                casadi_real *f_out, casadi_real *g_out)
 *
 *  @brief	Evaluate rhs f(x,z,t,u) and alg g(x,z,t,u) using generated code
 *
 *-------------------------------------------------------------------------------------------------
 */
static 
int eval_rhs_alg(const casadi_real *x, const casadi_real *z, casadi_real t, 
                 casadi_real I, casadi_real Tamb,
                 casadi_real *f_out, casadi_real *g_out)
{
    const casadi_real **arg;
    casadi_real **res; 
    casadi_real *w; 
    casadi_int *iw;
    casadi_int sz_arg=0, sz_res=0, sz_iw=0, sz_w=0;
    casadi_int mem;


    //----------------------------------------------------------------------------
    //  Input vector --  order must match input_parameter_order in Python
    //----------------------------------------------------------------------------
    casadi_real u[2] = { I, Tamb };


    //----------------------------------------------------------------------------
    // Allocate rhs_ working memories and call rhs_()
    //----------------------------------------------------------------------------
    rhs__work(&sz_arg,&sz_res,&sz_iw,&sz_w);
    arg = sz_arg ? (const casadi_real **)malloc(sz_arg * sizeof(casadi_real*)) : NULL;
    res = sz_res ? (casadi_real **)malloc(sz_res * sizeof(casadi_real*)) : NULL;
    iw = sz_iw ? (casadi_int *)malloc(sz_iw * sizeof(casadi_int)) : NULL;
    w = sz_w ? (casadi_real *)malloc(sz_w * sizeof(casadi_real)) : NULL;
    mem = rhs__alloc_mem();

    if (sz_arg < 5 || sz_res < 1)
       goto err_exit;
   
    arg[0] = x;
    arg[1] = z;
    arg[2] = &t;
    arg[3] = u;   
    arg[4] = NULL; // (if codegen includes parameters vector, leave NULL here)
    res[0] = f_out;
    rhs_(arg,res,iw,w,mem);
    
    rhs__free_mem(mem);
    free(w); 
    free(iw); 
    free(res); 
    free(arg);
   

    //----------------------------------------------------------------------------
    // Allocate alg_ working memories and all alg_()
    //----------------------------------------------------------------------------
    alg__work(&sz_arg,&sz_res,&sz_iw,&sz_w);
    arg = sz_arg ? (const casadi_real **)malloc(sz_arg * sizeof(casadi_real *)) : NULL;
    res = sz_res ? (casadi_real **)malloc(sz_res * sizeof(casadi_real *)) : NULL;
    iw = sz_iw ? (casadi_int *) malloc(sz_iw * sizeof(casadi_int)) : NULL;
    w = sz_w ? (casadi_real *) malloc(sz_w * sizeof(casadi_real)) : NULL;
    mem = alg__alloc_mem();

    if (sz_arg < 5 || sz_res < 1)
       goto err_exit;
   
    arg[0] = x;
    arg[1] = z;
    arg[2] = &t;
    arg[3] = u;   
    arg[4] = NULL; // (if codegen includes parameters vector, leave NULL here)
    res[0] = g_out;
    alg_(arg,res,iw,w,mem);

    alg__free_mem(mem);
    free(w); 
    free(iw); 
    free(res); 
    free(arg);

    return 0;

err_exit:
    return -1;
}


/*!
 *-------------------------------------------------------------------------------------------------
 *
 *  @fn		int step_implicit_euler(casadi_real *x, casadi_real *z, 
 *                                      casadi_int nx, casadi_int nz,
 *                                      casadi_real tnext, casadi_real dt, 
 *				        casadi_real I, casadi_real Tamb)
 *
 *  @brief	Newton solve for implicit Euler step:
 *                      F([x+,z+]) = [ x+ - x - dt * f(x+,z+) ; g(x+,z+) ] = 0
 *
 *-------------------------------------------------------------------------------------------------
 */
static 
int step_implicit_euler(casadi_real *x, casadi_real *z, 
		        casadi_int nx, casadi_int nz,
                        casadi_real tnext, casadi_real dt, 
			casadi_real I, casadi_real Tamb)
{
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


#if 1
    memcpy(xk, x, nx*sizeof(casadi_real));
    memcpy(zk, z, nz*sizeof(casadi_real));

    int iter, maxit=20;
    for (iter=0; iter<maxit; ++iter) {
        // F(y) with y=[xk, zk]
        eval_rhs_alg(xk, zk, tnext, I, Tamb, fx, gx);
        for (int i=0;i<nx;++i) F[i]     = xk[i] - x[i] - dt*fx[i];
        for (int j=0;j<nz;++j) F[nx+j]  = gx[j];
    

        // Check convergence
        casadi_real nrm=0;
        for (int i=0;i<n;++i) nrm = fmax(nrm, fabs(F[i]));
        if (nrm < 1e-8) { break; }

        // Numerical Jacobian (dense): columns for [x; z]
        // Perturb x
        for (int c=0;c<nx;++c) {
            casadi_real tmp = xk[c];
            xk[c] = tmp + EPSJ;
            eval_rhs_alg(xk, zk, tnext, I, Tamb, fx2, gx2);
            xk[c] = tmp;

            for (int r=0;r<nx;++r) J[r*n + c]      = ((xk[r] - x[r] - dt*fx2[r]) - F[r]) / EPSJ;
            for (int r=0;r<nz;++r) J[(nx+r)*n + c] = (gx2[r] - gx[r]) / EPSJ;
        }
        // Perturb z
        for (int c=0;c<nz;++c) {
            casadi_real tmp = zk[c];
            zk[c] = tmp + EPSJ;
            eval_rhs_alg(xk, zk, tnext, I, Tamb, fx2, gx2);
            zk[c] = tmp;

            for (int r=0;r<nx;++r) J[r*n + (nx+c)]      = ((xk[r] - x[r] - dt*fx2[r]) - F[r]) / EPSJ;
            for (int r=0;r<nz;++r) J[(nx+r)*n + (nx+c)] = (gx2[r] - gx[r]) / EPSJ;
        }


        // Solve J * dY = -F
        for (int i=0;i<n;++i) dY[i] = -F[i];
        if (solve_linear(n, J, dY)) { iter = -1; break; }

        // Update
        for (int i=0;i<nx;++i) xk[i] += dY[i];
        for (int j=0;j<nz;++j) zk[j] += dY[nx+j];
    }

    int ok = (iter>=0 && iter<maxit);
    if (ok) { memcpy(x, xk, nx*sizeof(casadi_real)); memcpy(z, zk, nz*sizeof(casadi_real)); }
#endif

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
 *  
 *  @fn		void eval_vars(const casadi_real *x, const casadi_real *z, casadi_real t, 
 *		               casadi_real I, casadi_real Tamb,
 *                             casadi_real *vars_out ) 
 *
 *  @brief	Evaluation variables
 *
 *  @param	'x', 'z' and 't' are inputs: 'f' and 'g' vectors, and time 't'
 *  @param	'I' and 'Tamb' are inputs
 *  @param	'vars_out' is output: pointer to pre-allocated 'variables' array
 *
 *-------------------------------------------------------------------------------------------------
 */
static
int eval_vars(const casadi_real *x, const casadi_real *z, casadi_real t, 
              casadi_real I, casadi_real Tamb,
              casadi_real *vars_out)
{
   casadi_int sz_arg=0, sz_res=0, sz_iw=0, sz_w=0;
   casadi_real u[2] = { I, Tamb };
   casadi_int mem = 0;

   variables_work(&sz_arg, &sz_res, &sz_iw, &sz_w);

   const casadi_real **arg = sz_arg ? (const casadi_real **)malloc(sz_arg * sizeof(casadi_real *)) : NULL;
   casadi_real **res = sz_res ? (casadi_real **) malloc(sz_res * sizeof(casadi_real *)) : NULL;
   casadi_int *iw = sz_iw ? malloc(sz_iw * sizeof(casadi_int)) : NULL;
   casadi_real *w = sz_w ? malloc(sz_w * sizeof(casadi_real)) : NULL;

   if (sz_arg < 5 || sz_res < 1)
   {
      perror("sz_arg < 5 or sz_res < 1 error");
      goto err_ret;
   }

   arg[0] = x; 
   arg[1] = z; 
   arg[2] = &t; 
   arg[3] = u; 
   arg[4] = NULL;
   res[0] = vars_out;


   mem = variables_alloc_mem();
   if (mem < 0)
   {
      perror("variables_alloc_mem() failed");
      goto err_ret;
   }

   // Call variables();
   int rc = variables(arg, res, iw, w, mem);
   if (rc)
   {
      perror("variables() failed");
      goto err_ret;
   }

   // Free memories
   variables_free_mem(mem);
   if (iw) free(iw);
   if (w) free(w);
   if (res) free(res);
   if (arg) free(arg);
   
   return 0;

err_ret:
   return -1;
}


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
    const casadi_real tpulse = 600.0;      // 10 min
    casadi_real m = fmod(t, period);
    return (m < tpulse) ? 1.0 : 0.0;
}



/*!
 *-------------------------------------------------------------------------------------------------
 *  Main program
 *-------------------------------------------------------------------------------------------------
 */
int main(void) 
{
    int rc = 0;
    const casadi_real dt = 0.01;       // 10 ms
    casadi_real t = 0.0;
    casadi_real Tamb = 298.15;         // allow changing per step
    const casadi_real **arg = NULL; 
    casadi_real **res = NULL;
    casadi_real *w = NULL;
    casadi_int *iw = NULL;
    casadi_int sz_arg=0, sz_res=0, sz_iw=0, sz_w=0, mem=0;
    casadi_real I = gitt_current(t);


    //----------------------------------------------------------------------------------
    // Initialize x = x0
    //----------------------------------------------------------------------------------
    // casadi_int nx = x0_n_out();   // length of x state
    casadi_int nx = rhs__sparsity_in(1)[0];   // length of x state
    printf("nx = %lld\n", nx);

    casadi_real *x = (casadi_real*)malloc(nx * sizeof(casadi_real));
    if (!x) 
    {
       perror("x = malloc() error.");
       goto err_ret;
    }

    x0_work(&sz_arg, &sz_res, &sz_iw, &sz_w);
    arg = sz_arg ? (const casadi_real **)malloc(sz_arg * sizeof (casadi_real *)) : NULL;
    res = sz_res ? (casadi_real **)malloc(sz_res * sizeof (casadi_real *)) : NULL;
    iw = sz_iw ? (casadi_int *) malloc(sz_iw * sizeof(casadi_int)) : NULL;
    w = sz_w ? (casadi_real *) malloc(sz_w * sizeof(casadi_real)) : NULL;
    printf("sz_arg=%lld, sz_res=%lld, sz_iw=%lld, sz_w=%lld\n", sz_arg, sz_res, sz_iw, sz_w);

    if (sz_res < 1) 
    {
       perror("sz_res < 1 error.");
       goto err_ret;
    }

    // Set pointer to returned values in 'x'
    res[0] = x;                    
		  
    mem = x0_alloc_mem();
    if (mem < 0) 
    {
       perror("x0_alloc_mem() failed");
       goto err_ret;
    }


    // Call x0()
    rc = x0(arg, res, iw, w, mem);
    if (rc)  
    {
       perror("x0() failed");
       goto err_ret;
    }

    // Print x0 results
    for (int i=0; i < nx; i++) 
       printf("x0[%d] = %.9g\n", i, x[i]);


    // Free memories
    if (mem > 0) x0_free_mem(mem);
    if (w) free(w);
    if (iw) free(iw);
    if (res) free(res);
    if (arg) free(arg);


    //----------------------------------------------------------------------------------
    //  Initiaize z0
    //----------------------------------------------------------------------------------
    // casadi_int nz = z0_n_out();   // length of z state
    casadi_int nz = rhs__sparsity_in(2)[0];   // length of z state
    casadi_real *z = (casadi_real*)malloc(nz * sizeof(casadi_real));
    if (!z) 
    {
       perror("z = malloc() error.");
       goto err_ret;
    }

    z0_work(&sz_arg, &sz_res, &sz_iw, &sz_w);
    printf("sz_arg=%lld, sz_res=%lld, sz_iw=%lld, sz_w=%lld\n", sz_arg, sz_res, sz_iw, sz_w);
    arg = sz_arg ? (const casadi_real **)malloc(sz_arg * sizeof(casadi_real *)) : NULL;
    res = sz_res ? (casadi_real **)malloc(sz_res * sizeof(casadi_real *)) : NULL;
    iw = sz_iw ? (casadi_int *) malloc(sz_iw * sizeof(casadi_int)) : NULL;
    w = sz_w ? (casadi_real *) malloc(sz_w * sizeof(casadi_real)) : NULL;

    if (sz_res < 1)
    {
       perror("z: sz_res < 1 error.");
       goto err_ret;
    }

    // Set pointer to receive returned value in 'z'
    res[0] = z;

    mem = z0_alloc_mem();
    if (mem < 0) 
    {
       perror("z0_alloc_mem() failed");
       goto err_ret;
    }

    // Call z0 
    rc = z0(arg, res, iw, w, mem);
    if (rc)
    {
       perror("z0() failed.");
       goto err_ret;
    }

    // Print z0 results
    for (int i=0; i < nz; i++)
       printf("z0]%d] = %.9g\n", i, z[i]);

    // Free memories
    z0_free_mem(mem);
    if (w) free(w);
    if (iw) free(iw);
    if (res) free(res);
    if (arg) free(arg);


    //----------------------------------------------------------------------------------
    // Allocate variable storage 
    //----------------------------------------------------------------------------------
    //casadi_int varN = variables_n_out();
    casadi_int varN = variables_sparsity_out(0)[0];
    casadi_real *vars = (casadi_real*)malloc(varN*sizeof(casadi_real));

    // CSV header
    printf("time_s,current_A,SOC,SOH,voltage_V,anode_V,cathode_V,LAM_neg_pct,LAM_pos_pct,LLI_pct\n");

    // Coulomb counting 
    casadi_real nominal_Ah = 0.0;


#if 1
    //----------------------------------------------------------------------------
    //  We’ll read Nominal cell capacity [A.h] from first variables eval below
    //  Run until voltage < 3.0 V
    //----------------------------------------------------------------------------
    for (int step = 0; step < 100000000; ++step) 
    {
	//------------------------------------------------------------------------
        //  Apply current based on GITT protocol 
	//------------------------------------------------------------------------
        I = gitt_current(t);


	//------------------------------------------------------------------------
        // Integration (implicit Euler)
	//------------------------------------------------------------------------
        if (step > 0) // at t=0 we print initial values first
        {    
            if (step_implicit_euler(x, z, nx, nz, t + dt, dt, I, Tamb)) 
	    {
                fprintf(stderr, "Newton failed at t=%.6f s\n", t);
                break;
            }
        }


	//------------------------------------------------------------------------
        // Evaluate variables
	//------------------------------------------------------------------------
        eval_vars(x, z, t, I, Tamb, vars);


	//------------------------------------------------------------------------
        // The export_vars[] order from Python:
        // 0: Voltage [V]
        // 1: Current [A]
        // 2: Discharge capacity [A.h]
        // 3: Nominal cell capacity [A.h]
        // 4: X-averaged cell temperature [K]
        // 5: Anode potential [V]
        // 6: Cathode potential [V]
        // 7: Loss of lithium inventory [%]
        // 8: Loss of active material in negative electrode [%]
        // 9: Loss of active material in positive electrode [%]
	//------------------------------------------------------------------------
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
        printf("%.3f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
               t, I, SOC, SOH, V, Va, Vc, LAMneg, LAMpos, LLIpct);


	//------------------------------------------------------------------------
	// Check stop condition 
	//------------------------------------------------------------------------
        if (V < 3.0) break;  

	
	//------------------------------------------------------------------------
	//  Increment simulation time
	//------------------------------------------------------------------------
        t += dt;
    }
#endif

    // Free state parameters and variables
    if (x) free(x); 
    if (z) free(z); 
    if (vars) free(vars);

    return 0;

err_ret:
    return -1;
}

