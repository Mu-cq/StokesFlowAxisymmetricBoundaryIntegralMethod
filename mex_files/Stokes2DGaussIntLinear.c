//
//  Stokes2DGaussIntLinear.c
//  
//
//  Created by Giacomo on 6/17/15.
//
//

#include "mex.h"
#include "math.h"
#include "stdio.h"

void Stokes2DGaussIntLinear(double *x, double *y, double *x0, double *y0, double *GXXa, double *GXYa, double *GYYa, double *TXXXa, double *TXXYa, double *TXYYa, double *TYYYa, double *GXXb, double *GXYb, double *GYYb, double *TXXXb, double *TXXYb, double *TXYYb, double *TYYYb, int field, int elem)
{
  int i, j;
  double SXX, SXY, SYY, QXXX, QXXY, QXYY, QYYY; //green's functions
  double dx, dy, dl, xMid, yMid;    //elements variables
  double distXA, distXB, distYA, distYB, distA, distB; //for sing treatment
  int dologLeft, dologRight, dolog;
  double DX, DY, r, phiA, phiB, GPx, GPy;   //green function variables
  
  //Gauss points and Gauss weigths
  double GP[] = {-0.932469514203152, -0.661209386466265, -0.238619186083197, 0.238619186083197, 0.661209386466265, 0.932469514203152};
  double GW[] = {0.171324492379170, 0.360761573048139, 0.467913934572691, 0.467913934572691, 0.360761573048139, 0.171324492379170};
  
  for (j=0; j<elem; j++) {  
      //elemnts variables      
      dx = x[j+1]-x[j];
      dy = y[j+1]-y[j];
      dl = sqrt(dx*dx+dy*dy);
      xMid = (x[j+1]+x[j])/2;
      yMid = (y[j+1]+y[j])/2;
      
      for (i=0; i<field; i++) {
          
          //decide if doing singular treatment
          distXA = x[j]-x0[i];
          distYA = y[j]-y0[i];
          distA = sqrt(distXA*distXA+distYA*distYA);
          dologLeft = (distA>10e-8);
          
          //decide if doing singular treatment
          distXB = x[j+1]-x0[i];
          distYB = y[j+1]-y0[i];
          distB = sqrt(distXB*distXB+distYB*distYB);
          dologRight = (distB>10e-8);
          
          dolog = (dologLeft+dologRight > 1.9);
          
          //gauss integration
          for (int l=0; l<6; l++) {
                //gauss point
                GPx = GP[l]/2*dx+xMid;
                GPy = GP[l]/2*dy+yMid;
                
                //compute green function
                DX = GPx-x0[i];
                DY = GPy-y0[i];
                r = sqrt(DX*DX+DY*DY);
                
                //hat function
                phiA = 1-(GP[l]+1)/2;
                phiB = (GP[l]+1)/2;

                //single layer
                SXX = -log(r)*dolog+(DX*DX)/pow(r,2);
                SXY = (DX*DY)/pow(r,2);
                SYY = -log(r)*dolog+(DY*DY)/pow(r,2);

                //double layer
                QXXX = -4*(DX*DX*DX)/pow(r,4)*dolog;
                QXXY = -4*(DX*DX*DY)/pow(r,4)*dolog;
                QXYY = -4*(DX*DY*DY)/pow(r,4)*dolog;
                QYYY = -4*(DY*DY*DY)/pow(r,4)*dolog;
              
                //integration phiA SINGLE LAYER
                GXXa[j*field+i] += phiA*GW[l]*dl*SXX/2;
                GXYa[j*field+i] += phiA*GW[l]*dl*SXY/2;
                GYYa[j*field+i] += phiA*GW[l]*dl*SYY/2;
                
                //integration phiA DOUBLE LAYER
                TXXXa[j*field+i] += phiA*GW[l]*dl*QXXX/2;
                TXXYa[j*field+i] += phiA*GW[l]*dl*QXXY/2;
                TXYYa[j*field+i] += phiA*GW[l]*dl*QXYY/2;
                TYYYa[j*field+i] += phiA*GW[l]*dl*QYYY/2;
                
                //integration phiB SINGLE LAYER
                GXXb[j*field+i] += phiB*GW[l]*dl*SXX/2;
                GXYb[j*field+i] += phiB*GW[l]*dl*SXY/2;
                GYYb[j*field+i] += phiB*GW[l]*dl*SYY/2;
                
                //integration phiB DOUBLE LAYER
                TXXXb[j*field+i] += phiB*GW[l]*dl*QXXX/2;
                TXXYb[j*field+i] += phiB*GW[l]*dl*QXXY/2;
                TXYYb[j*field+i] += phiB*GW[l]*dl*QXYY/2;
                TYYYb[j*field+i] += phiB*GW[l]*dl*QYYY/2;
          }
          
          //singular treatment: add part form analytical integration (as in research notes)          
          GXXa[j*field+i] += (0.75*dl-0.5*dl*log(dl))*(1-dologLeft);
          GYYa[j*field+i] += (0.75*dl-0.5*dl*log(dl))*(1-dologLeft);
          
          GXXb[j*field+i] += (0.25*dl-0.5*dl*log(dl))*(1-dologLeft);
          GYYb[j*field+i] += (0.25*dl-0.5*dl*log(dl))*(1-dologLeft);
          
          GXXa[j*field+i] += (0.25*dl-0.5*dl*log(dl))*(1-dologRight);
          GYYa[j*field+i] += (0.25*dl-0.5*dl*log(dl))*(1-dologRight);
          
          GXXb[j*field+i] += (0.75*dl-0.5*dl*log(dl))*(1-dologRight);
          GYYb[j*field+i] += (0.75*dl-0.5*dl*log(dl))*(1-dologRight);
          
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
double *GXXa;      /* output matrix 1x(NxM) */
double *GXYa;      /* output matrix 1x(NxM) */
double *GYYa;      /* output matrix 1x(NxM) */
double *TXXXa;      /* output matrix 1x(NxM) */
double *TXXYa;      /* output matrix 1x(NxM) */
double *TXYYa;      /* output matrix 1x(NxM) */
double *TYYYa;      /* output matrix 1x(NxM) */
double *GXXb;      /* output matrix 1x(NxM) */
double *GXYb;      /* output matrix 1x(NxM) */
double *GYYb;      /* output matrix 1x(NxM) */
double *TXXXb;      /* output matrix 1x(NxM) */
double *TXXYb;      /* output matrix 1x(NxM) */
double *TXYYb;      /* output matrix 1x(NxM) */
double *TYYYb;      /* output matrix 1x(NxM) */

/* get the value of the scalar input  */
//multiplier = mxGetScalar(prhs[0]);

/* code here */
/* create a pointer to the real data in the input matrix  */
inMatrixX = mxGetPr(prhs[0]);
inMatrixY = mxGetPr(prhs[1]);
inMatrixX0 = mxGetPr(prhs[2]);
inMatrixY0 = mxGetPr(prhs[3]);
    
/* get dimensions of the input matrix */
nrows = mxGetN(prhs[0]);
ncols = mxGetN(prhs[2]);
    
/* create the output matrix (their dimension depending on the inputs) */
plhs[0] = mxCreateDoubleMatrix(1,nrows*ncols-ncols,mxREAL);
plhs[1] = mxCreateDoubleMatrix(1,nrows*ncols-ncols,mxREAL);
plhs[2] = mxCreateDoubleMatrix(1,nrows*ncols-ncols,mxREAL);
plhs[3] = mxCreateDoubleMatrix(1,nrows*ncols-ncols,mxREAL);
plhs[4] = mxCreateDoubleMatrix(1,nrows*ncols-ncols,mxREAL);
plhs[5] = mxCreateDoubleMatrix(1,nrows*ncols-ncols,mxREAL);
plhs[6] = mxCreateDoubleMatrix(1,nrows*ncols-ncols,mxREAL);
plhs[7] = mxCreateDoubleMatrix(1,nrows*ncols-ncols,mxREAL);
plhs[8] = mxCreateDoubleMatrix(1,nrows*ncols-ncols,mxREAL);
plhs[9] = mxCreateDoubleMatrix(1,nrows*ncols-ncols,mxREAL);
plhs[10] = mxCreateDoubleMatrix(1,nrows*ncols-ncols,mxREAL);
plhs[11] = mxCreateDoubleMatrix(1,nrows*ncols-ncols,mxREAL);
plhs[12] = mxCreateDoubleMatrix(1,nrows*ncols-ncols,mxREAL);
plhs[13] = mxCreateDoubleMatrix(1,nrows*ncols-ncols,mxREAL);

/* get a pointer to the real data in the output matrix */
GXXa = mxGetPr(plhs[0]);
GXYa = mxGetPr(plhs[1]);
GYYa = mxGetPr(plhs[2]);
TXXXa = mxGetPr(plhs[3]);
TXXYa = mxGetPr(plhs[4]);
TXYYa = mxGetPr(plhs[5]);
TYYYa = mxGetPr(plhs[6]);
GXXb = mxGetPr(plhs[7]);
GXYb = mxGetPr(plhs[8]);
GYYb = mxGetPr(plhs[9]);
TXXXb = mxGetPr(plhs[10]);
TXXYb = mxGetPr(plhs[11]);
TXYYb = mxGetPr(plhs[12]);
TYYYb = mxGetPr(plhs[13]);

/* call the computational routine */
Stokes2DGaussIntLinear(inMatrixX, inMatrixY, inMatrixX0, inMatrixY0, GXXa, GXYa, GYYa, TXXXa, TXXYa, TXYYa, TYYYa, GXXb, GXYb, GYYb, TXXXb, TXXYb, TXYYb, TYYYb, ncols, nrows-1);
}
