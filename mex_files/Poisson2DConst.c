//
//  Poisson2DConst.c
//  
//
//  Created by Giacomo on 7/8/15.
//
//

#include "mex.h"
#include "math.h"
#include <stdio.h>

void Poisson2DConst(double *x, double *y, double *x0, double *y0, double *G, double *Gx, double *Gy, int field, int elem)
{
  int i, j, dolog;
  double P;
  double Px, Py;
  double r, r2, dx, dy, DX, DY, GPx, GPy;
  double dl, xMid, yMid, distX, distY, dist;
  
  //define pi
  double pi = 3.141592653589793;
  
  //Gauss points and Gauss weigths
  double GP[] = {-0.932469514203152, -0.661209386466265, -0.238619186083197, 0.238619186083197, 0.661209386466265, 0.932469514203152};
  double GW[] = {0.171324492379170, 0.360761573048139, 0.467913934572691, 0.467913934572691, 0.360761573048139, 0.171324492379170};
  
  for (j=0; j<elem; j++) {  
      //elements variables
      dx = x[j+1]-x[j];
      dy = y[j+1]-y[j];
      dl = sqrt(dx*dx+dy*dy);
      xMid = (x[j+1]+x[j])/2;
      yMid = (y[j+1]+y[j])/2;
      
      for (i=0; i<field; i++) {
          
          //decide if doing singular tretment
          distX = xMid-x0[i];
          distY = yMid-y0[i];
          dist = sqrt(distX*distX+distY*distY);
          dolog = (dist>10e-8);
          //dolog = 1;
          
          //gauss integration
          for (int l=0; l<6; l++) {
                GPx = GP[l]/2*dx+xMid;
                GPy = GP[l]/2*dy+yMid;
                
                //compute green function
                DX = GPx-x0[i];
                DY = GPy-y0[i];
                r = sqrt(DX*DX+DY*DY);
                r2 = r*r;

                //single layer
                P = 0.5/pi*log(r)*dolog;
                
                //printf(" %f \n",r);

                //double layer
                Px = -0.5/pi*DX/r2;
                Py = -0.5/pi*DY/r2;
              
                //integration SINGLE LAYER
                G[j*field+i] += GW[l]*dl*P/2;
                
                //integration DOUBLE LAYER
                Gx[j*field+i] += GW[l]*dl*Px/2*dolog;
                Gy[j*field+i] += GW[l]*dl*Py/2*dolog;
                
          }
          
          //singular treatment
          G[j*field+i] += (-0.5/pi*(dl-dl*log(dl/2)))*(1-dolog);
                    
      }
  }
}

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
/* variable declarations here */
double *inMatrixX;       /* 1xN input matrix */
double *inMatrixY;       /* 1xN input matrix */
mwSize nrows;           /* size of matrix */
double *inMatrixX0;       /* 1xM input matrix */
double *inMatrixY0;       /* 1xM input matrix */
mwSize ncols;           /* size of matrix */
double *G;      /* output matrix 1x(NxM) */
double *Gx;      /* output matrix 1x(NxM) */
double *Gy;      /* output matrix 1x(NxM) */

/* code here */
/* create a pointer to the real data in the input matrix  */
inMatrixX = mxGetPr(prhs[0]);
inMatrixY = mxGetPr(prhs[1]);
inMatrixX0 = mxGetPr(prhs[2]);
inMatrixY0 = mxGetPr(prhs[3]);
    
/* get dimensions of the input matrix */
nrows = mxGetN(prhs[0]);
ncols = mxGetN(prhs[2]);
    
/* create the output matrix */
plhs[0] = mxCreateDoubleMatrix(1,nrows*ncols-ncols,mxREAL);
plhs[1] = mxCreateDoubleMatrix(1,nrows*ncols-ncols,mxREAL);
plhs[2] = mxCreateDoubleMatrix(1,nrows*ncols-ncols,mxREAL);

/* get a pointer to the real data in the output matrix */
G = mxGetPr(plhs[0]);
Gx = mxGetPr(plhs[1]);
Gy = mxGetPr(plhs[2]);

/* call the computational routine */
Poisson2DConst(inMatrixX, inMatrixY, inMatrixX0, inMatrixY0, G, Gx, Gy, ncols, nrows-1);
}