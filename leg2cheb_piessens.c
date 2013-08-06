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
*    p(x) = X[0]P_0(x) + x[1]P_1(x) + ... x[N-1]P_{N-1}(x)
*         = C[0]T_0(x) + C[1]T_1(x) + ... C[N-1]T_{N-1}(x)
******************************************************************************/

#include <mex.h>
#include <math.h>
#include <stdio.h> 

void leg2cheb(int N, double *x, double *c);
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

void leg2cheb(int N, double *x, double *c)
{   
    long int i, j;
    double Mii, Mij;
    
    for ( i = 0 ; i < N ; i++ ) {
        if ( i <= 1 ) {
            Mii = 1.;
        }
        /* Initialise off-diagonal: */
        Mij = Mii;
        c[i] = 0;
        /* Loop over columns: */
        for ( j = i ; j < N ; j += 2 ) {
            /* Matrix multiply: */
            c[i] += Mij*x[j];
            /* Update off-daigonal via recurrence: */
            Mij = (1.-((2.*j+3.)/(j-i+2.))/(j+i+2.))*Mij;
        }
        /* Update diagonal via recurrence: */
        Mii = (i+.5)/(i+1.)*Mii;
    }
   
    return;
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *x, *c;
  mwSize N;
  
  /*  Create a pointer to the input matrix x  */
  N = mxGetNumberOfElements(prhs[0]);
  x = mxGetPr(prhs[0]);
  
  /*  Set the output pointer to the output matrix */
  plhs[0] = mxCreateDoubleMatrix(N, 1, mxREAL);
  
  /*  Create a C++ pointer to a copy of the output matrix */
  c = mxGetPr(plhs[0]);
  
  /*  Call the C++ subroutine */
  leg2cheb( N, x, c );
 
  return;
}

