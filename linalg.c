/*!
 *=======================================================================================================
 *
 *  @file	linalg.c 
 *
 *  @brief	linear algebraic solver
 *
 *=======================================================================================================
 */
#include <math.h>
#include <stdlib.h>
#include "linalg.h"


/*!
 *-------------------------------------------------------------------------------------------------------
 *
 *  @fn		int lu_decompose(double *A, int n, int *piv)
 *  
 *  @brief	LU decomposition
 *
 *-------------------------------------------------------------------------------------------------------
 */
int lu_decompose(double *A, int n, int *piv)
{
   for (int i=0;i<n;i++) 
      piv[i]=i;

   for (int k=0;k<n;k++)
   {
      // pivot
      int piv_row=k; 
      double amax=fabs(A[k*n+k]);
      for (int i=k+1;i<n;i++)
      { 
	 double v=fabs(A[i*n+k]); 
	 if (v>amax)
	 {
	    amax=v; 
	    piv_row=i;
         } 
      }

      if (amax<1e-18) 
      {
	 return -1; // singular
      }

      if (piv_row!=k)
      {
         for (int j=0;j<n;j++)
	 { 
            double tmp=A[k*n+j]; 
	    A[k*n+j]=A[piv_row*n+j]; 
	    A[piv_row*n+j]=tmp; 
	 }

         int t=piv[k]; 
	 piv[k]=piv[piv_row]; 
	 piv[piv_row]=t;
      }

      // factorize
      for (int i=k+1;i<n;i++)
      {
         A[i*n+k] /= A[k*n+k];
         double lik = A[i*n+k];
         for (int j=k+1;j<n;j++) 
            A[i*n+j] -= lik*A[k*n+j];
      }
   }

   return 0;
}


/*!
 *-------------------------------------------------------------------------------------------------------
 *
 *  @fn		int lu_solve(const double *LU, int n, const int *piv, double *b)
 *
 *  @brief	LU solver
 *
 *-------------------------------------------------------------------------------------------------------
 */
int lu_solve(const double *LU, int n, const int *piv, double *b)
{
   double *y = (double*)malloc(n*sizeof(double));

   // Apply permutation to b
   for (int i=0;i<n;i++) 
      y[i]=b[piv[i]];

   // Forward solve Ly = Pb
   for (int i=0;i<n;i++)
   {
      for (int j=0;j<i;j++) 
	 y[i] -= LU[i*n+j]*y[j];
   }

   // Backward solve Ux = y
   for (int i=n-1;i>=0;i--)
   {
      for (int j=i+1;j<n;j++) 
	 y[i] -= LU[i*n+j]*y[j];

      y[i] /= LU[i*n+i];
   }

   // Copy to b
   for (int i=0;i<n;i++) 
      b[i]=y[i];

   return 0;
}
