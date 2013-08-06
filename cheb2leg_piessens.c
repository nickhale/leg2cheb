/*******************************************************************************
 * This file is part of LEG2CHEB.
 * Copyright (c) 2013 Nick Hale and Alex Townsend, University of Oxford
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 ******************************************************************************/

/*******************************************************************************
* C = cheb2leg(X)
*  Converts the chebyshev coefficients X to the Legendre coefficients C.
*  X and C are ordered such that
*    p(x) = X[0]T_0(x) + x[1]T_1(x) + ... x[N-1]T_{N-1}(x)
*         = C[0]P_0(x) + C[1]P_1(x) + ... C[N-1]P_{N-1}(x)
******************************************************************************/

#include <mex.h>
#include <math.h>
#include <stdio.h> 

void cheb2leg(int N, double *x, double *c);
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

void cheb2leg(int N, double *x, double *c)
{   
    long int i, j;
    
    /* Note the rows of L are scaled by 1/(i+.5), which is aplied at the end */
    double Lii, Lij;
    
    /* First coefficient: */
    c[0] = x[0];
    for ( j = 2 ; j < N ; j += 2 ) {
        /* Update off-daigonal via recurrence: */
        Lij = -1. / (j*j - 1.);
        /* Matrix multiply: */
        c[0] += Lij * x[j];
    }
    
    /* Remaining coefficients: */
    Lii = 1.;
    for ( i = 1 ; i < N ; i++ ) {
        /* Update diagonal via recurrence: */
        Lii = (double)i / (i+.5)*Lii;
        /* Initialise off-diagonal: */
        Lij = Lii;
        /* Initialise coefficient: */
        c[i] = 0;
        /* Loop over columns: */
        for ( j = i ; j < N ; j += 2 ) {
            /* Matrix multiply: */
            c[i] += Lij * x[j];
            /* Update off-daigonal via recurrence: */
            Lij = (1 - 2.*(2.*j*(j+2.)+i*(i+1.))/((j-i+2.)*(i+j+3.))/(double)j) * Lij;
        }
        /* Scale properly: */
        c[i] *= i + 0.5;
    }
    
    return;
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *x, *c;
  mwSize N;
  
  /* Create a pointer to the input matrix x and determine its size: */
  N = mxGetNumberOfElements(prhs[0]);
  x = mxGetPr(prhs[0]);
  
  /* Set the output pointer to the output matrix: */
  plhs[0] = mxCreateDoubleMatrix(N, 1, mxREAL);
  
  /* Create a C++ pointer to a copy of the output matrix */
  c = mxGetPr(plhs[0]);
  
  /* Call the C++ subroutine */
  cheb2leg( N, x, c );
 
  return;
}

