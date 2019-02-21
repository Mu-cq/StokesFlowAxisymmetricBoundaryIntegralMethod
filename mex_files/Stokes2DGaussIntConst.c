//
//  Stokes2DGaussIntConst.c
//  
//
//  Created by Giacomo on 6/16/15.
//
//

#include "mex.h"
#include "math.h"

void Stokes2DGaussIntConst(double *x, double *y, double *x0, double *y0, double *GXX, double *GXY, double *GYY, double *TXXX, double *TXXY, double *TXYY, double *TYYY, int field, int elem)
{
  int i;
  int j;
  double SXX;  double SXY;    double SYY;
  double QXXX; double QXXY;   double QXYY;   double QYYY;
  
  //Gauss points and Gauss weigths
  double GP[] = {-0.932469514203152, -0.661209386466265, -0.238619186083197, 0.238619186083197, 0.661209386466265, 0.932469514203152};
  double GW[] = {0.171324492379170, 0.360761573048139, 0.467913934572691, 0.467913934572691, 0.360761573048139, 0.171324492379170};
  
  for (j=0; j<elem; j++) {  
      //elements variables
      double dx = x[j+1]-x[j];
      double dy = y[j+1]-y[j];
      double dl = sqrt(dx*dx+dy*dy);
      double xMid = (x[j+1]+x[j])/2;
      double yMid = (y[j+1]+y[j])/2;
      double distX; double distY;   double dist;
      int dolog;
      
      for (i=0; i<field; i++) {
          
          //decide if doing singular tretment
          distX = xMid-x0[i];
          distY = yMid-y0[i];
          dist = sqrt(distX*distX+distY*distY);
          int dolog = (dist>10e-8);
          
          //gauss integration
          for (int l=0; l<6; l++) {
                double GPx = GP[l]/2*dx+xMid;
                double GPy = GP[l]/2*dy+yMid;
                
                //compute green function
                double DX = GPx-x0[i];
                double DY = GPy-y0[i];
                double r = sqrt(DX*DX+DY*DY);

                //single layer
                SXX = -log(r)*dolog+(DX*DX)/pow(r,2);
                SXY = (DX*DY)/pow(r,2);
                SYY = -log(r)*dolog+(DY*DY)/pow(r,2);

                //double layer
                QXXX = -4*(DX*DX*DX)/pow(r,4);
                QXXY = -4*(DX*DX*DY)/pow(r,4);
                QXYY = -4*(DX*DY*DY)/pow(r,4);
                QYYY = -4*(DY*DY*DY)/pow(r,4);
                //GT_2DStokes_cfunction(GPx, GPy, x0[i], y0[i], SXX, SXY, SYY, QXXX, QXXY, QXYY, QYYY);
              
                //integration SINGLE LAYER
                GXX[j*field+i] += GW[l]*dl*SXX/2;
                GXY[j*field+i] += GW[l]*dl*SXY/2;
                GYY[j*field+i] += GW[l]*dl*SYY/2;
                
                //integration DOUBLE LAYER
                TXXX[j*field+i] += GW[l]*dl*QXXX/2*dolog;
                TXXY[j*field+i] += GW[l]*dl*QXXY/2*dolog;
                TXYY[j*field+i] += GW[l]*dl*QXYY/2*dolog;
                TYYY[j*field+i] += GW[l]*dl*QYYY/2*dolog;
          }
          
          //singular treatment
          GXX[j*field+i] += (-dl*log(dl/2)+dl)*(1-dolog);
          GYY[j*field+i] += (-dl*log(dl/2)+dl)*(1-dolog);
          
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
double *GXX;      /* output matrix 1x(NxM) */
double *GXY;      /* output matrix 1x(NxM) */
double *GYY;      /* output matrix 1x(NxM) */
double *TXXX;      /* output matrix 1x(NxM) */
double *TXXY;      /* output matrix 1x(NxM) */
double *TXYY;      /* output matrix 1x(NxM) */
double *TYYY;      /* output matrix 1x(NxM) */
//int elem;   int field;

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
plhs[3] = mxCreateDoubleMatrix(1,nrows*ncols-ncols,mxREAL);
plhs[4] = mxCreateDoubleMatrix(1,nrows*ncols-ncols,mxREAL);
plhs[5] = mxCreateDoubleMatrix(1,nrows*ncols-ncols,mxREAL);
plhs[6] = mxCreateDoubleMatrix(1,nrows*ncols-ncols,mxREAL);

/* get a pointer to the real data in the output matrix */
GXX = mxGetPr(plhs[0]);
GXY = mxGetPr(plhs[1]);
GYY = mxGetPr(plhs[2]);
TXXX = mxGetPr(plhs[3]);
TXXY = mxGetPr(plhs[4]);
TXYY = mxGetPr(plhs[5]);
TYYY = mxGetPr(plhs[6]);

/* call the computational routine */
Stokes2DGaussIntConst(inMatrixX, inMatrixY, inMatrixX0, inMatrixY0, GXX, GXY, GYY, TXXX, TXXY, TXYY, TYYY, ncols, nrows-1);
}