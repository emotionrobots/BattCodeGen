/*!
 *=====================================================================================================
 *
 *  @file	linalg.h 
 *
 *  @brief	Linear algebra header
 *
 *=====================================================================================================
 */
#ifndef __LINALG_H__
#define __LINALG_H__

#ifdef __cplusplus
extern "C" {
#endif

int lu_decompose(double *A, int n, int *piv); // in-place LU with partial pivoting
int lu_solve(const double *LU, int n, const int *piv, double *b);

#ifdef __cplusplus
}
#endif

#endif // __LINALG_H__
