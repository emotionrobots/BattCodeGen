/*!
 *=============================================================================================================
 *
 *  @file 	global.h	
 *
 *  @brief	Global defines for Simulate DFN model of NMC battery
 *
 *=============================================================================================================
 */
#ifndef __GLOBAL_H__
#define __GLOBAL_H__
#include "batt_model.h"

// Variable array index
enum {
   IDX_VOLTAGE,		// 0
   IDX_CURRENT,		// 1
   IDX_Q_DIS,		// 2
   IDX_Q_NOM,		// 3
   IDX_AVG_CELL_TEMP,	// 4
   IDX_N_V,   		// 5 anode V (-)
   IDX_P_V,   		// 6 cathode V (+)
   IDX_LLI,		// 7
   IDX_N_LAM,		// 8
   IDX_P_LAM		// 9
};


// Model context holding pre-allocated workspaces shared across calls
typedef struct {
   int nx, nz, np, nvars;
   double *x, *z, *p, *vars;
   double *f, *g;  
   casadi_int iw_sz, w_sz;
   casadi_int *iw; 
   double *w; 
   const char **var_names;
} ModelCtx;

#endif // __GLOBAL_H__


