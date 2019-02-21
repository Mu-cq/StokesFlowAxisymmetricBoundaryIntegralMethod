//
//  Stokes2DAxisIntConstParallel.c
//  
//
//  Created by Giacomo on 10/26/15.
//
//

#include <stdio.h>
#include "mex.h"
#include "math.h"
#include <omp.h>

void Stokes2DAxisIntConstParallel(double *x, double *y, double *x0, double *y0, double *GXX, double *GXY, double *GYX, double *GYY, double *TXXX, double *TXXY, double *TXYY, double *TYXX, double *TYXY, double *TYYY, double *PYX, double *PYY, int field, int elem, int nnn)
{
  int dolog, axis;
  double SXX, SXY, SYX, SYY;
  double QXXX, QXXY, QXYY, QYXX, QYXY, QYYY;
  double pyx, pyy;
  double F, E, G, B, D, C, k, k2;   //elliptic integrals
  double I10, I11, I30, I31, I32, I50, I51, I52, I53;
  double r, r0, r02, r03, d, d2, GPx, GPy, d3, d5;
  double k2p, k4, k5, k6, yy5;
  double FCTR, RL30, RL50, RL52, RL54, RL56;
  double DX, DY, DX2, DX3, r2, r3;
  double dx, dy, dl, xMid, yMid, distX, distY, dist;
  double tempGXX, tempGXY, tempGYX, tempGYY, tempTXXX, tempTXXY, tempTXYY, tempTYXX, tempTYXY, tempTYYY, tempPYX, tempPYY; 
  
  //define pi and tol
  double pi = 3.141592653589793;
  double tol = 0.0000000000001;
  
  //Gauss points and Gauss weigths
  double GP[] = {-0.932469514203152, -0.661209386466265, -0.238619186083197, 0.238619186083197, 0.661209386466265, 0.932469514203152};
  double GW[] = {0.171324492379170, 0.360761573048139, 0.467913934572691, 0.467913934572691, 0.360761573048139, 0.171324492379170};
  
  omp_set_num_threads(nnn);
  
  #pragma omp parallel for private(dx, dy, dl, xMid, yMid, axis, distX, distY, dist, dolog, GPx, GPy, DX, DY, d2, d, r0, r, k2, k, r2, r02, r03, r3, k2p, k4, k6, k5, yy5, DX2, DX3, F, E, G, B, D, C, I10, I11, I30, I31, I32, RL30, RL50, RL52, RL54, RL56, FCTR, I50, I51, I52, I53, SXX, SXY, SYX, SYY, QXXX, QXXY, QXYY, QYXX, QYXY, QYYY, pyx, pyy) reduction(+:tempGXX,tempGXY,tempGYX,tempGYY,tempTXXX,tempTXXY,tempTXYY,tempTYXX,tempTYXY,tempTYYY,tempPYX,tempPYY)
  for (int j=0; j<elem; j++) {  
      //elements variables
      dx = x[j+1]-x[j];
      dy = y[j+1]-y[j];
      dl = sqrt(dx*dx+dy*dy);
      xMid = (x[j+1]+x[j])/2;
      yMid = (y[j+1]+y[j])/2;
      
      for (int i=0; i<field; i++) {
          
          //check if I'm on the axis
          axis = (y0[i] < 10e-8);
          
          //decide if doing singular tretment
          distX = xMid-x0[i];
          distY = yMid-y0[i];
          dist = sqrt(distX*distX+distY*distY);
          dolog = (dist > 10e-8);
          
          tempGXX = 0;
          tempGXY = 0;
          tempGYX = 0;
          tempGYY = 0;
          tempTXXX = 0;
          tempTXXY = 0;
          tempTXYY = 0;
          tempTYXX = 0;
          tempTYXY = 0;
          tempTYYY = 0;
          tempPYX = 0;
          tempPYY = 0;
          
          //gauss integration
          for (int l=0; l<6; l++) {
                GPx = GP[l]/2*dx+xMid;
                GPy = GP[l]/2*dy+yMid;
                
                //compute green function physical variables
                DX = GPx-x0[i];
                DY = GPy-y0[i];
                d2 = DX*DX+DY*DY;
                d = sqrt(d2);
                r0 = y0[i];
                r = GPy;
                k2 = (4*r0*r)/(DX*DX+pow((r+r0),2));
                k = sqrt(k2);
                r2 = r*r;
                r02 = r0*r0;
                r03 = r0*r0*r0;
                r3 = r*r*r;
                
                //dummy variables
                k2p = 1-k2;
                k4 = k2*k2;
                k6 = k2*k2*k2;
                k5 = k4*k;
                yy5 = pow(r0*r,2.5);
                DX2 = DX*DX;
                DX3 = DX*DX*DX;
                
                //compute elliptic integral of first and second kind
                F = 0.5*pi;
                E = 1.0;
                G = 1.0;
                B = k;
                D = 10*tol;

                while (sqrt(D*D) > tol) {
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
                
                //dummy variables for double layer
                RL30 = E/k2p;
                RL50 =  (2.0 * (2.0 - k2)*RL30 - F) / (3.0*k2p);
                RL52 =  (RL50 - RL30)/k2;
                RL54 =  (F-RL50+2.0*k2*RL52)/k4;
                RL56 = -(E-RL50+3.0*k2*RL52-3.0*k4*RL54)/ k6;
                
                //often used variables for double layer
                FCTR = k5/(8.0*yy5);
                I50 = FCTR * RL50;
                I51 = FCTR * (2.0 * RL52-RL50);
                I52 = FCTR * (4.0 * (RL54-RL52)+RL50);
                I53 = FCTR * (8.0 * RL56 - 12.0 *RL54 +6.0 * RL52 - RL50);
                
                if (axis) {
                    
                     d3  = d2*d;
                     d5  = d3*d2;
                     I10 = 2*pi/d;
                     I11 = 0.0;
                     I30 = 2*pi/d3;
                     I31 = 0.0;
                     I32 = pi/d3;

                     I50 = 2*pi/d5;
                     I51 = 0.0;
                     I52 = pi/d5;
                     I53 = 0.0;
                }
                
                //single layer
                SXX = r*(I10+DX*DX*I30) + 2*(1-dolog)*log(d)*(1-axis);
                SXY = r*DX*(r*I30-r0*I31);
                SYX = r*DX*(r*I31-r0*I30);
                SYY = r*(I11+(r*r+r0*r0)*I31-r*r0*(I30+I32)) + 2*(1-dolog)*log(d)*(1-axis);

                //double layer
                QXXX = -6*r*DX3*I50;
                QXXY = -6*r*DX2*(r*I50-r0*I51);
                QXYY = -6*r*DX*(r0*r0*I52+r*r*I50-2*r*r0*I51);
                QYXX = -6*r*DX2*(r*I51-r0*I50);
                QYXY = -6*r*DX*((r*r+r0*r0)*I51-r*r0*(I50+I52));
                QYYY = -6*r*(r3*I51-r2*r0*(I50+2*I52)+r*r0*r0*(I53+2*I51)-r0*r0*r0*I52);
                
                //desingularization double layer
                pyx = - 6*r*DX * (r2*I52 - 2*r*r0*I51 + r02*I50);
                pyy = - 6*r*(r02*r * (2.0 * I52+I50) - r03 * I51 - r2*r0 * (2.0*I51+I53) + r3 * I52);
              
                //integration SINGLE LAYER
                tempGXX += GW[l]*dl*SXX/2;
                tempGXY += GW[l]*dl*SXY/2;
                tempGYX += GW[l]*dl*SYX/2;
                tempGYY += GW[l]*dl*SYY/2;
                
                //integration DOUBLE LAYER
                tempTXXX += GW[l]*dl*QXXX/2;
                tempTXXY += GW[l]*dl*QXXY/2;
                tempTXYY += GW[l]*dl*QXYY/2;
                tempTYXX += GW[l]*dl*QYXX/2;
                tempTYXY += GW[l]*dl*QYXY/2;
                tempTYYY += GW[l]*dl*QYYY/2;
                
                //integration for desigularization
                tempPYX += GW[l]*dl*pyx/2;
                tempPYY += GW[l]*dl*pyy/2;
                
          }
          
          //integration SINGLE LAYER
          GXX[j*field+i] = tempGXX;
          GXY[j*field+i] = tempGXY;
          GYX[j*field+i] = tempGYX;
          GYY[j*field+i] = tempGYY;
                
          //integration DOUBLE LAYER
          TXXX[j*field+i] = tempTXXX;
          TXXY[j*field+i] = tempTXXY;
          TXYY[j*field+i] = tempTXYY;
          TYXX[j*field+i] = tempTYXX;
          TYXY[j*field+i] = tempTYXY;
          TYYY[j*field+i] = tempTYYY;
                
          //integration for desigularization
          PYX[j*field+i] = tempPYX;
          PYY[j*field+i] = tempPYY;
          
          //singular treatment
          GXX[j*field+i] = tempGXX + 2*(-dl*log(dl/2)+dl)*(1-dolog)*(1-axis);
          GYY[j*field+i] = tempGYY + 2*(-dl*log(dl/2)+dl)*(1-dolog)*(1-axis);
          
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
double *GYX;      /* output matrix 1x(NxM) */
double *GYY;      /* output matrix 1x(NxM) */
double *TXXX;      /* output matrix 1x(NxM) */
double *TXXY;      /* output matrix 1x(NxM) */
double *TXYY;      /* output matrix 1x(NxM) */
double *TYXX;      /* output matrix 1x(NxM) */
double *TYXY;      /* output matrix 1x(NxM) */
double *TYYY;      /* output matrix 1x(NxM) */
double *PYX;      /* output matrix 1x(NxM) */
double *PYY;      /* output matrix 1x(NxM) */
//int elem;   int field;

/* code here */
/* create a pointer to the real data in the input matrix  */
inMatrixX = mxGetPr(prhs[0]);
inMatrixY = mxGetPr(prhs[1]);
inMatrixX0 = mxGetPr(prhs[2]);
inMatrixY0 = mxGetPr(prhs[3]);
int n = mxGetScalar(prhs[4]);
    
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
plhs[7] = mxCreateDoubleMatrix(1,nrows*ncols-ncols,mxREAL);
plhs[8] = mxCreateDoubleMatrix(1,nrows*ncols-ncols,mxREAL);
plhs[9] = mxCreateDoubleMatrix(1,nrows*ncols-ncols,mxREAL);
plhs[10] = mxCreateDoubleMatrix(1,nrows*ncols-ncols,mxREAL);
plhs[11] = mxCreateDoubleMatrix(1,nrows*ncols-ncols,mxREAL);

/* get a pointer to the real data in the output matrix */
GXX = mxGetPr(plhs[0]);
GXY = mxGetPr(plhs[1]);
GYX = mxGetPr(plhs[2]);
GYY = mxGetPr(plhs[3]);
TXXX = mxGetPr(plhs[4]);
TXXY = mxGetPr(plhs[5]);
TXYY = mxGetPr(plhs[6]);
TYXX = mxGetPr(plhs[7]);
TYXY = mxGetPr(plhs[8]);
TYYY = mxGetPr(plhs[9]);
PYX = mxGetPr(plhs[10]);
PYY = mxGetPr(plhs[11]);

/* call the computational routine */
Stokes2DAxisIntConstParallel(inMatrixX, inMatrixY, inMatrixX0, inMatrixY0, GXX, GXY, GYX, GYY, TXXX, TXXY, TXYY, TYXX, TYXY, TYYY, PYX, PYY, ncols, nrows-1, n);
}