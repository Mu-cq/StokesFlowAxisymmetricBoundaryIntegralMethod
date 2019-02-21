//
//  Stokes2DAxisIntLinear.c
//  
//
//  Created by Giacomo on 6/25/15.
//
//

#include "mex.h"
#include "math.h"
#include "stdio.h"

void Stokes2DGaussIntLinear(double *x, double *y, double *x0, double *y0, double *GXXa, double *GXYa, double *GYXa, double *GYYa, double *TXXXa, double *TXXYa, double *TXYYa, double *TYXXa, double *TYXYa, double *TYYYa, double *GXXb, double *GXYb, double *GYXb, double *GYYb, double *TXXXb, double *TXXYb, double *TXYYb, double *TYXXb, double *TYXYb, double *TYYYb, int field, int elem)
{
  int i, j, dolog1, dolog2, dolog, axis;
  double SXX, SXY, SYX, SYY;
  double QXXX, QXXY, QXYY, QYXX, QYXY, QYYY;
  double F, E, G, B, D, C, k, k2;   //elliptic integrals
  double I10, I11, I30, I31, I32, I50, I51, I52, I53;
  double r, r0, d, d2, d3, d5, GPx, GPy, phiA, phiB;
  double k2p, k4, k5, k6, yy5;
  double FCTR, RL30, RL50, RL52, RL54, RL56;
  double DX, DY, DX2, DX3, r2, r3;
  double dx, dy, dl, xMid, yMid, distXA, distYA, distA, distXB, distYB, distB;
  
  //define pi and tol
  double pi = 3.141592653589793;
  double tol = 0.0000000000001;
  
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
          
          //check if I'm on the axis
          axis = (y0[i] < 10e-8);
          
          //decide if doing singular treatment
          distXA = x[j]-x0[i];
          distYA = y[j]-y0[i];
          distA = sqrt(distXA*distXA+distYA*distYA);
          dolog1 = (distA>10e-8);
          
          //decide if doing singular treatment
          distXB = x[j+1]-x0[i];
          distYB = y[j+1]-y0[i];
          distB = sqrt(distXB*distXB+distYB*distYB);
          dolog2 = (distB>10e-8);
          
          dolog = (dolog1+dolog2 > 1.9);
          
          //gauss integration
          for (int l=0; l<6; l++) {
                //gauss point
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
                    
                
                //hat function
                phiA = 1-(GP[l]+1)/2;
                phiB = (GP[l]+1)/2;
                
                //single layer
                SXX = r*(I10+DX*DX*I30) + 2*log(d)*(1-dolog)*(1-axis);
                SXY = r*DX*(r*I30-r0*I31);
                SYX = r*DX*(r*I31-r0*I30);
                SYY = r*(I11+(r*r+r0*r0)*I31-r*r0*(I30+I32)) + 2*log(d)*(1-dolog)*(1-axis);
                
                //double layer
                QXXX = -6*r*DX3*I50;
                QXXY = -6*r*DX2*(r*I50-r0*I51);
                QXYY = -6*r*DX*(r0*r0*I52+r*r*I50-2*r*r0*I51);
                QYXX = -6*r*DX2*(r*I51-r0*I50);
                QYXY = -6*r*DX*((r*r+r0*r0)*I51-r*r0*(I50+I52));
                QYYY = -6*r*(r3*I51-r2*r0*(I50+2*I52)+r*r0*r0*(I53+2*I51)-r0*r0*r0*I52);
              
                //integration phiA SINGLE LAYER
                GXXa[j*field+i] += phiA*GW[l]*dl*SXX/2;
                GXYa[j*field+i] += phiA*GW[l]*dl*SXY/2;
                GYXa[j*field+i] += phiA*GW[l]*dl*SYX/2;
                GYYa[j*field+i] += phiA*GW[l]*dl*SYY/2;
                
                //integration phiA DOUBLE LAYER
                TXXXa[j*field+i] += phiA*GW[l]*dl*QXXX/2;
                TXXYa[j*field+i] += phiA*GW[l]*dl*QXXY/2;
                TXYYa[j*field+i] += phiA*GW[l]*dl*QXYY/2;
                TYXXa[j*field+i] += phiA*GW[l]*dl*QYXX/2;
                TYXYa[j*field+i] += phiA*GW[l]*dl*QYXY/2;
                TYYYa[j*field+i] += phiA*GW[l]*dl*QYYY/2;
                
                //integration phiB SINGLE LAYER
                GXXb[j*field+i] += phiB*GW[l]*dl*SXX/2;
                GXYb[j*field+i] += phiB*GW[l]*dl*SXY/2;
                GYXb[j*field+i] += phiB*GW[l]*dl*SYX/2;
                GYYb[j*field+i] += phiB*GW[l]*dl*SYY/2;
                
                //integration phiB DOUBLE LAYER
                TXXXb[j*field+i] += phiB*GW[l]*dl*QXXX/2;
                TXXYb[j*field+i] += phiB*GW[l]*dl*QXXY/2;
                TXYYb[j*field+i] += phiB*GW[l]*dl*QXYY/2;
                TYXXb[j*field+i] += phiB*GW[l]*dl*QYXX/2;
                TYXYb[j*field+i] += phiB*GW[l]*dl*QYXY/2;
                TYYYb[j*field+i] += phiB*GW[l]*dl*QYYY/2;
          }
                    
          //singular treatment: add part form analytical integration (as in research notes) IF I'M ON THE AXIS I DON'T NEED SINGULAR TREATMENT
          GXXa[j*field+i] += (1.5*dl-dl*log(dl))*(1-dolog1)*(1-axis);
          GYYa[j*field+i] += (1.5*dl-dl*log(dl))*(1-dolog1)*(1-axis);
          
          GXXb[j*field+i] += (0.5*dl-dl*log(dl))*(1-dolog1)*(1-axis);
          GYYb[j*field+i] += (0.5*dl-dl*log(dl))*(1-dolog1)*(1-axis);
          
          GXXa[j*field+i] += (0.5*dl-dl*log(dl))*(1-dolog2)*(1-axis);
          GYYa[j*field+i] += (0.5*dl-dl*log(dl))*(1-dolog2)*(1-axis);
          
          GXXb[j*field+i] += (1.5*dl-dl*log(dl))*(1-dolog2)*(1-axis);
          GYYb[j*field+i] += (1.5*dl-dl*log(dl))*(1-dolog2)*(1-axis);
          
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
double *GYXa;      /* output matrix 1x(NxM) */
double *GYYa;      /* output matrix 1x(NxM) */
double *TXXXa;      /* output matrix 1x(NxM) */
double *TXXYa;      /* output matrix 1x(NxM) */
double *TXYYa;      /* output matrix 1x(NxM) */
double *TYXXa;      /* output matrix 1x(NxM) */
double *TYXYa;      /* output matrix 1x(NxM) */
double *TYYYa;      /* output matrix 1x(NxM) */
double *GXXb;      /* output matrix 1x(NxM) */
double *GXYb;      /* output matrix 1x(NxM) */
double *GYXb;      /* output matrix 1x(NxM) */
double *GYYb;      /* output matrix 1x(NxM) */
double *TXXXb;      /* output matrix 1x(NxM) */
double *TXXYb;      /* output matrix 1x(NxM) */
double *TXYYb;      /* output matrix 1x(NxM) */
double *TYXXb;      /* output matrix 1x(NxM) */
double *TYXYb;      /* output matrix 1x(NxM) */
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
plhs[14] = mxCreateDoubleMatrix(1,nrows*ncols-ncols,mxREAL);
plhs[15] = mxCreateDoubleMatrix(1,nrows*ncols-ncols,mxREAL);
plhs[16] = mxCreateDoubleMatrix(1,nrows*ncols-ncols,mxREAL);
plhs[17] = mxCreateDoubleMatrix(1,nrows*ncols-ncols,mxREAL);
plhs[18] = mxCreateDoubleMatrix(1,nrows*ncols-ncols,mxREAL);
plhs[19] = mxCreateDoubleMatrix(1,nrows*ncols-ncols,mxREAL);

/* get a pointer to the real data in the output matrix */
GXXa = mxGetPr(plhs[0]);
GXYa = mxGetPr(plhs[1]);
GYXa = mxGetPr(plhs[2]);
GYYa = mxGetPr(plhs[3]);
TXXXa = mxGetPr(plhs[4]);
TXXYa = mxGetPr(plhs[5]);
TXYYa = mxGetPr(plhs[6]);
TYXXa = mxGetPr(plhs[7]);
TYXYa = mxGetPr(plhs[8]);
TYYYa = mxGetPr(plhs[9]);
GXXb = mxGetPr(plhs[10]);
GXYb = mxGetPr(plhs[11]);
GYXb = mxGetPr(plhs[12]);
GYYb = mxGetPr(plhs[13]);
TXXXb = mxGetPr(plhs[14]);
TXXYb = mxGetPr(plhs[15]);
TXYYb = mxGetPr(plhs[16]);
TYXXb = mxGetPr(plhs[17]);
TYXYb = mxGetPr(plhs[18]);
TYYYb = mxGetPr(plhs[19]);

/* call the computational routine */
Stokes2DGaussIntLinear(inMatrixX, inMatrixY, inMatrixX0, inMatrixY0, GXXa, GXYa, GYXa, GYYa, TXXXa, TXXYa, TXYYa, TYXXa, TYXYa, TYYYa, GXXb, GXYb, GYXb, GYYb, TXXXb, TXXYb, TXYYb, TYXXb, TYXYb, TYYYb, ncols, nrows-1);
}
