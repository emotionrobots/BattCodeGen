// cgen/runner.c : run DFN at 10 ms steps with GITT 1A pulses; stop at 3.0 V.
// Build with Makefile below. Assumes code generated into batt_model.c

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "batt_model.h"         // generated header (because with_header=true)

#ifndef CASADI_REAL
#define CASADI_REAL double
#endif

#ifndef CASADI_INT
#define CASADI_INT  long long int
#endif


// Finite-difference step for numerical Jacobians
static const CASADI_REAL EPSJ = 1e-8;

// Helper: dense linear solve (Gaussian elimination, no pivoting – ok for small systems)
static int solve_linear(int n, CASADI_REAL *A, CASADI_REAL *b) {
    for (int k = 0; k < n; ++k) {
        CASADI_REAL piv = A[k*n + k];
        if (fabs(piv) < 1e-14) return 1; // singular
        CASADI_REAL inv = 1.0 / piv;
        for (int j = k; j < n; ++j) A[k*n + j] *= inv;
        b[k] *= inv;
        for (int i = 0; i < n; ++i) if (i != k) {
            CASADI_REAL f = A[i*n + k];
            for (int j = k; j < n; ++j) A[i*n + j] -= f * A[k*n + j];
            b[i] -= f * b[k];
        }
    }
    return 0;
}

// Evaluate rhs f(x,z,t,u) and alg g(x,z,t,u) using generated code
static void eval_rhs_alg(
    const CASADI_REAL *x, const CASADI_REAL *z,
    CASADI_REAL t, CASADI_REAL I, CASADI_REAL Tamb,
    CASADI_REAL *f_out, CASADI_REAL *g_out
){
    const CASADI_REAL *arg[5];
    CASADI_REAL *res_f[1]; CASADI_REAL *res_g[1];
    // CASADI_REAL iw_dummy[1]; CASADI_REAL w_dummy[1];

    CASADI_REAL u[2] = { I, Tamb };

    arg[0] = x;
    arg[1] = z;
    arg[2] = &t;
    arg[3] = u;   // order must match input_parameter_order in Python
    arg[4] = NULL; // (if codegen includes parameters vector, leave NULL here)
    res_f[0] = f_out;

    rhs_(arg, res_f, NULL, NULL, 0);

    res_g[0] = g_out;
    alg_(arg, res_g, NULL, NULL, 0);
}

// Newton solve for implicit Euler step:
// F([x+,z+]) = [ x+ - x - dt * f(x+,z+) ; g(x+,z+) ] = 0
static int step_implicit_euler(
    CASADI_REAL *x, CASADI_REAL *z, int nx, int nz,
    CASADI_REAL tnext, CASADI_REAL dt, CASADI_REAL I, CASADI_REAL Tamb
){
    const int n = nx + nz;
    CASADI_REAL *F  = (CASADI_REAL*)calloc(n, sizeof(CASADI_REAL));
    CASADI_REAL *J  = (CASADI_REAL*)calloc(n*n, sizeof(CASADI_REAL));
    CASADI_REAL *dY = (CASADI_REAL*)calloc(n, sizeof(CASADI_REAL));
    CASADI_REAL *xk = (CASADI_REAL*)malloc(nx*sizeof(CASADI_REAL));
    CASADI_REAL *zk = (CASADI_REAL*)malloc(nz*sizeof(CASADI_REAL));
    CASADI_REAL *fx = (CASADI_REAL*)malloc(nx*sizeof(CASADI_REAL));
    CASADI_REAL *gx = (CASADI_REAL*)malloc(nz*sizeof(CASADI_REAL));
    CASADI_REAL *fx2= (CASADI_REAL*)malloc(nx*sizeof(CASADI_REAL));
    CASADI_REAL *gx2= (CASADI_REAL*)malloc(nz*sizeof(CASADI_REAL));

    memcpy(xk, x, nx*sizeof(CASADI_REAL));
    memcpy(zk, z, nz*sizeof(CASADI_REAL));

    int iter, maxit=20;
    for (iter=0; iter<maxit; ++iter) {
        // F(y) with y=[xk, zk]
        eval_rhs_alg(xk, zk, tnext, I, Tamb, fx, gx);
        for (int i=0;i<nx;++i) F[i]     = xk[i] - x[i] - dt*fx[i];
        for (int j=0;j<nz;++j) F[nx+j]  = gx[j];

        // Check convergence
        CASADI_REAL nrm=0;
        for (int i=0;i<n;++i) nrm = fmax(nrm, fabs(F[i]));
        if (nrm < 1e-8) { break; }

        // Numerical Jacobian (dense): columns for [x; z]
        // Perturb x
        for (int c=0;c<nx;++c) {
            CASADI_REAL tmp = xk[c];
            xk[c] = tmp + EPSJ;
            eval_rhs_alg(xk, zk, tnext, I, Tamb, fx2, gx2);
            xk[c] = tmp;

            for (int r=0;r<nx;++r) J[r*n + c]      = ((xk[r] - x[r] - dt*fx2[r]) - F[r]) / EPSJ;
            for (int r=0;r<nz;++r) J[(nx+r)*n + c] = (gx2[r] - gx[r]) / EPSJ;
        }
        // Perturb z
        for (int c=0;c<nz;++c) {
            CASADI_REAL tmp = zk[c];
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
    if (ok) { memcpy(x, xk, nx*sizeof(CASADI_REAL)); memcpy(z, zk, nz*sizeof(CASADI_REAL)); }

    free(F); free(J); free(dY); free(xk); free(zk); free(fx); free(gx); free(fx2); free(gx2);
    return ok ? 0 : 1;
}

// Evaluate exported variables[] (Voltage, Current, capacity, temps, LLI/LAM…)
static void eval_vars(
    const CASADI_REAL *x, const CASADI_REAL *z,
    CASADI_REAL t, CASADI_REAL I, CASADI_REAL Tamb,
    CASADI_REAL *vars_out /* length = variables_n_out( ) */
){
    const CASADI_REAL *arg[5];
    CASADI_REAL *res[1];
    CASADI_REAL u[2] = { I, Tamb };

    arg[0] = x; arg[1] = z; arg[2] = &t; arg[3] = u; arg[4] = NULL;
    res[0] = vars_out;
    variables(arg, res, NULL, NULL, 0);
}

static CASADI_REAL gitt_current(double t) {
    // Simple GITT: 10 min pulse @ +1 A, then 50 min rest @ 0 A; repeat.
    const double period = 3600.0;     // 60 min
    const double tpulse = 600.0;      // 10 min
    double m = fmod(t, period);
    return (m < tpulse) ? 1.0 : 0.0;
}


/*!
 *--------------------------------------------------------------------------------------------------------
 *  Main program
 *--------------------------------------------------------------------------------------------------------
 */
int main(void) 
{
    int rc = 0;
    const CASADI_REAL dt = 0.01;       // 10 ms
    CASADI_REAL t = 0.0;
    CASADI_REAL Tamb = 298.15;         // allow changing per step
    CASADI_REAL I = gitt_current(t);
    CASADI_REAL u[2] = { I, Tamb };
    const CASADI_REAL **arg = NULL; 
    CASADI_REAL *res[1] = {0};
    CASADI_REAL *w = NULL;
    CASADI_INT *iw = NULL;
    CASADI_INT sz_arg=0, sz_res=0, sz_iw=0, sz_w=0, mem=0;


    //----------------------------------------------------------------------------------
    // Query problem sizes from generated API:
    CASADI_INT nx = x0_n_out();   // length of x state
    CASADI_REAL *x = (CASADI_REAL*)malloc(nx * sizeof(CASADI_REAL));
    if (!x) 
    {
       perror("malloc x");
       return 1;
    }
    printf("nx = %lld\n", nx);

    // Initial conditions from generated code
    x0_work(&sz_arg, &sz_res, &sz_iw, &sz_w);
    iw = sz_iw ? (CASADI_INT *) malloc(sz_iw * sizeof(CASADI_INT)) : NULL;
    w = sz_w ? (CASADI_REAL *) malloc(sz_w * sizeof(CASADI_REAL)) : NULL;

    res[0] = x;
    arg = NULL;
    mem = x0_alloc_mem();
    if (mem < 0) 
    {
       perror("x0_alloc_mem() failed");
       return 1;
    }

    // X0 initialization 
    rc = x0(arg, res, iw, w, mem);
    if (rc) 
    {
       perror("x0()");
       x0_free_mem(mem);
       free(iw);
       free(w);
       free(x);
       return 1;
    }
    printf("x0[0] = %.9g\n", (CASADI_REAL)x[0]);
    if (nx > 1) printf("x0[1] = %.9g\n", (CASADI_REAL)x[1]);


    printf("Step 1\n");

    //----------------------------------------------------------------------------------
    // Query problem sizes from generated API:
    CASADI_INT nz = z0_n_out();   // length of z state
    CASADI_REAL *z = (CASADI_REAL*)malloc(nz * sizeof(CASADI_REAL));
    if (!z) 
    {
       perror("malloc z");
       return 1;
    }
    printf("nz = %lld\n", nz);

    // Initial conditions from generated code
    z0_work(&sz_arg, &sz_res, &sz_iw, &sz_w);
    iw = sz_iw ? (CASADI_INT *) malloc(sz_iw * sizeof(CASADI_INT)) : NULL;
    w = sz_w ? (CASADI_REAL *) malloc(sz_w * sizeof(CASADI_REAL)) : NULL;

    res[0] = z;
    arg = NULL;
    mem = z0_alloc_mem();
    if (mem < 0) 
    {
       perror("x0_alloc_mem() failed");
       return 1;
    }

    // Z0 initialization 
    rc = z0(arg, res, iw, w, 0);
    if (rc)
    {
       perror("z0()");
       z0_free_mem(mem);
       free(iw);
       free(w);
       free(z);
       return 1;
    }
    printf("z0[0] = %.9g\n", (CASADI_REAL)z[0]);
    if (nz > 1) printf("z0[1] = %.9g\n", (CASADI_REAL)z[1]);

    printf("Step 3\n");

#if 0
    // Time loop
    int varN = variables_n_out();
    CASADI_REAL *vars = (CASADI_REAL*)malloc(varN*sizeof(CASADI_REAL));

    // CSV header
    printf("time_s,current_A,SOC,SOUH,voltage_V,anode_V,cathode_V,LAM_neg_pct,LAM_pos_pct,LLI_pct\n");

    // Helper accumulators
    CASADI_REAL nominal_Ah = 0.0;

    // We’ll read Nominal cell capacity [A.h] from first variables eval below
    // Run until voltage < 3.0 V
    for (int step = 0; step < 100000000; ++step) {
        I = gitt_current(t);

        // Integrate to t+dt (implicit Euler)
        if (step > 0) { // at t=0 we print initial values first
            if (step_implicit_euler(x, z, nx, nz, t + dt, dt, I, Tamb)) {
                fprintf(stderr, "Newton failed at t=%.6f s\n", t);
                break;
            }
        }

        // Evaluate variables
        eval_vars(x, z, t, I, Tamb, vars);
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

        CASADI_REAL V      = vars[0];
        CASADI_REAL Imeas  = vars[1];
        CASADI_REAL QdisAh = vars[2];
        nominal_Ah         = (step==0 ? vars[3] : nominal_Ah);
        CASADI_REAL Va     = vars[5];
        CASADI_REAL Vc     = vars[6];
        CASADI_REAL LLIpct = vars[7];
        CASADI_REAL LAMneg = vars[8];
        CASADI_REAL LAMpos = vars[9];

        // SOC and a simple SOH proxy:
        // SOC = 1 - Q_dis / Nominal
        double SOC = (nominal_Ah > 0) ? fmax(0.0, 1.0 - (double)(QdisAh / nominal_Ah)) : NAN;
        // SOH (capacity based) ≈ 1 - max(LAMneg, LAMpos)/100  (simple, illustrative)
        double SOH = 1.0 - (double)(fmax(LAMneg, LAMpos) / 100.0);

        printf("%.3f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
               (double)t, (double)I, SOC, SOH, (double)V,
               (double)Va, (double)Vc, (double)LAMneg, (double)LAMpos, (double)LLIpct);

        if (V < 3.0) break;  // stop condition per prompt

        t += dt;
    }
#endif

    free(x); 
    free(z); 

    // free(vars);

    return 0;
}

