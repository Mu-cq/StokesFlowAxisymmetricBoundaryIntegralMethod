//
//  GT_2DStokes_cfunction.c
//  
//
//  Created by Giacomo on 6/16/15.
//
//

#include "mex.h"
#include "math.h"

void GT_2DStokes_cfunction(double *x, double *y, double *x0, double *y0,double *SXX, double *SXY, double *SYX, double *SYY, double *QXXX, double *QXXY, double *QXYX, double *QXYY, double *QYXX, double *QYXY, double *QYYX, double *QYYY, int m, int n)
{
  int i;
  int j;
  
  for (i=0; i<m; i++) {      
      for (j=0; j<n; j++) {
          double dx = x[i]-x0[j];
          double dy = y[i]-y0[j];
          double r = sqrt(dx*dx+dy*dy);
          
          //single layer
          SXX[i*n+j] = -log(r)+(dx*dx)/pow(r,2);
          SXY[i*n+j] = (dx*dy)/pow(r,2);
          SYX[i*n+j] = SXY[i*n+j];
          //SYX[i*n+j] = (dy*dx)/pow(r,2);
          SYY[i*n+j] = -log(r)+(dy*dy)/pow(r,2);
          
          //double layer
          QXXX[i*n+j] = -4*(dx*dx*dx)/pow(r,4);
          QXXY[i*n+j] = -4*(dx*dx*dy)/pow(r,4);
          QXYX[i*n+j] = QXXY[i*n+j];
          //QXYX[i*n+j] = -4*(dx*dy*dx)/pow(r,4);
          QXYY[i*n+j] = -4*(dx*dy*dy)/pow(r,4);
          QYXX[i*n+j] = QXXY[i*n+j];
          QYXY[i*n+j] = QXYY[i*n+j];
          QYYX[i*n+j] = QXYY[i*n+j];
          //QYXX[i*n+j] = -4*(dx*dx*dy)/pow(r,4);
          //QYXY[i*n+j] = -4*(dy*dx*dy)/pow(r,4);
          //QYYX[i*n+j] = -4*(dy*dy*dx)/pow(r,4);
          QYYY[i*n+j] = -4*(dy*dy*dy)/pow(r,4);
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
double *SXX;      /* output matrix 1x(NxM) */
double *SXY;      /* output matrix 1x(NxM) */
double *SYX;      /* output matrix 1x(NxM) */
double *SYY;      /* output matrix 1x(NxM) */
double *QXXX;      /* output matrix 1x(NxM) */
double *QXXY;      /* output matrix 1x(NxM) */
double *QXYX;      /* output matrix 1x(NxM) */
double *QXYY;      /* output matrix 1x(NxM) */
double *QYXX;      /* output matrix 1x(NxM) */
double *QYXY;      /* output matrix 1x(NxM) */
double *QYYX;      /* output matrix 1x(NxM) */
double *QYYY;      /* output matrix 1x(NxM) */

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
plhs[0] = mxCreateDoubleMatrix(1,nrows*ncols,mxREAL);
plhs[1] = mxCreateDoubleMatrix(1,nrows*ncols,mxREAL);
plhs[2] = mxCreateDoubleMatrix(1,nrows*ncols,mxREAL);
plhs[3] = mxCreateDoubleMatrix(1,nrows*ncols,mxREAL);
plhs[4] = mxCreateDoubleMatrix(1,nrows*ncols,mxREAL);
plhs[5] = mxCreateDoubleMatrix(1,nrows*ncols,mxREAL);
plhs[6] = mxCreateDoubleMatrix(1,nrows*ncols,mxREAL);
plhs[7] = mxCreateDoubleMatrix(1,nrows*ncols,mxREAL);
plhs[8] = mxCreateDoubleMatrix(1,nrows*ncols,mxREAL);
plhs[9] = mxCreateDoubleMatrix(1,nrows*ncols,mxREAL);
plhs[10] = mxCreateDoubleMatrix(1,nrows*ncols,mxREAL);
plhs[11] = mxCreateDoubleMatrix(1,nrows*ncols,mxREAL);

/* get a pointer to the real data in the output matrix */
SXX = mxGetPr(plhs[0]);
SXY = mxGetPr(plhs[1]);
SYX = mxGetPr(plhs[2]);
SYY = mxGetPr(plhs[3]);
QXXX = mxGetPr(plhs[4]);
QXXY = mxGetPr(plhs[5]);
QXYX = mxGetPr(plhs[6]);
QXYY = mxGetPr(plhs[7]);
QYXX = mxGetPr(plhs[8]);
QYXY = mxGetPr(plhs[9]);
QYYX = mxGetPr(plhs[10]);
QYYY = mxGetPr(plhs[11]);

/* call the computational routine */
GT_2DStokes_cfunction(inMatrixX,inMatrixY,inMatrixX0,inMatrixY0,SXX,SXY,SYX,SYY,QXXX,QXXY,QXYX,QXYY,QYXX,QYXY,QYYX,QYYY,nrows,ncols);
}