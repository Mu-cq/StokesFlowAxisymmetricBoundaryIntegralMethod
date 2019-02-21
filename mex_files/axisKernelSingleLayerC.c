//
//  axisKernelSingleLayerC.c
//  
//
//  Created by Giacomo on 13/1/17.
//
//

#include "mex.h"
#include "math.h"
#include "stdio.h"

void axisKernelSingleLayerC(double *x, double *y, double *x0, double *y0, double *SXX, double *SXY, double *SYX, double *SYY, int field, int elem)
{
  int i, j, axis;
  double xHere, yHere;
  double F, E, G, B, D, C, k, k2;   //elliptic integrals
  double I10, I11, I30, I31, I32;
  double r, r0, d, d2, d3;
  double DX, DY, DX2, r2;
  
  //FILE *pFile1;
  //FILE *pFile2;
  //pFile1 = fopen ("myfile1.txt","w");
  //pFile2 = fopen ("myfile2.txt","w");
  
  //define pi and tol
  double pi = M_PI;
  double tol = 0.0000000000001;
  
  for (j=0; j<elem; j++) {  
      
      //x and y here
      xHere = x[j];
      yHere = y[j];
      
      for (i=0; i<field; i++) {
          
          //check if I'm on the axis
          axis = (y0[i] < 10e-8);
       
          //compute green function physical variables
          DX = xHere-x0[i];
          DY = yHere-y0[i];
          d2 = DX*DX+DY*DY;
          d = sqrt(d2);
          r0 = y0[i];
          r = yHere;
          k2 = (4*r0*r)/(DX*DX+pow((r+r0),2));
          k = sqrt(k2);
          r2 = r*r;
                
          //dummy variables
          DX2 = DX*DX;
                
          //compute elliptic integral of first and second kind
          F = 0.5*pi;
          E = 1.0;
          G = 1.0;
          B = k;
          D = 10*tol;

          while (fabs(D) > tol) {
                  C = sqrt(1.0-B*B);
                  B = (1.0-C)/(1.0+C);
                  D = F*B;
                  F = F+D;
                  G = 0.50*G*B;
                  E = E+G;
          }

          E = F*(1.0-0.50*k2*E);
                
          //often used variables for single layer
          I10 = 2*k/sqrt(r0*r)*F;
          I11 = k/pow(r0*r,1.5)*((DX2+r0*r0+r2)*F-(DX2+r0*r0+r2+2*r*r0)*E);
          I30 = 2*k/pow(r0*r,0.5)*E/d2;
          I31 = k/pow(r0*r,1.5)*(-F+(r0*r0+r2+DX2)/d2*E);
          I32 = k/pow(r0*r,2.5)*(-F*(DX2+r0*r0+r2) + (DX2*DX2+2.0*DX2*r2+2.0*DX2*r0*r0+r2*r2+r0*r0*r0*r0)*E/d2);
                
          if (axis) {
                    
                     d3  = d2*d;
                     I10 = 2*pi/d;
                     I11 = 0.0;
                     I30 = 2*pi/d3;
                     I31 = 0.0;
                     I32 = pi/d3;

          }
                
          //single layer
          SXX[j*field+i] = r*(I10+DX*DX*I30);
          SXY[j*field+i] = r*DX*(r*I30-r0*I31);
          SYX[j*field+i] = r*DX*(r*I31-r0*I30);
          SYY[j*field+i] = r*(I11+(r*r+r0*r0)*I31-r*r0*(I30+I32));
          
          //build output matrix
          //GXX[j*field+i] += GXXa[j*field+i];
          //GXX[(j+1)*field+i] += GXXb[j*field+i];
          
      }
  }
  
}

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
/* variable declarations here */
mwSize nrows;           /* size of matrix */
double *inMatrixX0;       /* 1xM input matrix */
double *inMatrixY0;       /* 1xM input matrix */
double *inMatrixX;       /* 1xM input matrix */
double *inMatrixY;       /* 1xM input matrix */
mwSize ncols;           /* size of matrix */
double *SXX;      /* output matrix 1x(NxM) */
double *SXY;      /* output matrix 1x(NxM) */
double *SYX;      /* output matrix 1x(NxM) */
double *SYY;      /* output matrix 1x(NxM) */


/* get the value of the scalar input  */
//multiplier = mxGetScalar(prhs[0]);

/* code here */
/* create a pointer to the real data in the input matrix  */
inMatrixX = mxGetPr(prhs[0]);
inMatrixY = mxGetPr(prhs[1]);
inMatrixX0 = mxGetPr(prhs[2]);
inMatrixY0 = mxGetPr(prhs[3]);
    
/* get dimensions of the input matrix */
nrows = mxGetN(prhs[2]);
ncols = mxGetN(prhs[0]);
    
/* create the output matrix (their dimension depending on the inputs) */
plhs[0] = mxCreateDoubleMatrix(1,nrows*ncols,mxREAL);
plhs[1] = mxCreateDoubleMatrix(1,nrows*ncols,mxREAL);
plhs[2] = mxCreateDoubleMatrix(1,nrows*ncols,mxREAL);
plhs[3] = mxCreateDoubleMatrix(1,nrows*ncols,mxREAL);

/* get a pointer to the real data in the output matrix */
SXX = mxGetPr(plhs[0]);
SXY = mxGetPr(plhs[1]);
SYX = mxGetPr(plhs[2]);
SYY = mxGetPr(plhs[3]);

/* call the computational routine */
axisKernelSingleLayerC(inMatrixX, inMatrixY, inMatrixX0, inMatrixY0, SXX, SXY, SYX, SYY, ncols, nrows);
}
